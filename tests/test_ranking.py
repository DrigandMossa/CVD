"""Ranked, prebuilt vector options: scoring, ordering, modifications and the text built from them."""
import unittest

from cvd.design import build_design
from cvd.model import Config
from cvd.parts import STRATEGIES, WPRE_SHORT_LEN
from cvd.ranking import (SETTING_CHOICES, basis_for_option, compare_line, config_for_option,
                         explain_option, explain_option_parts, fit_text, modification_parts,
                         options_report, options_table, parse_explain, parse_setting,
                         profile_insert, rank_vector_options, ranking_summary)
from cvd.ranking.tables import CRITERIA


def orf(codons):
    return "ATG" + "GCC" * (codons - 1)


def cassette(design):
    return sum(f.length for f in design.features if f.in_cassette)


def rank(strategy="rep", codons=394, setting="", optimize=True, **kw):
    cfg = Config(strategy=strategy, gene="X", sequence=orf(codons), **kw)
    return rank_vector_options(cfg, setting, optimize)


class Setting(unittest.TestCase):
    def test_parse_setting_accepts_the_choices_and_blank(self):
        self.assertEqual(parse_setting("In vivo"), "in vivo")
        self.assertEqual(parse_setting("ex vivo"), "ex vivo")
        self.assertEqual(parse_setting(SETTING_CHOICES[0]), "")
        self.assertEqual(parse_setting(""), "")

    def test_parse_setting_rejects_anything_else(self):
        with self.assertRaises(ValueError):
            parse_setting("in the lab")

    def test_setting_decides_the_winning_vector(self):
        self.assertEqual(rank(setting="in vivo").options[0].vector, "aav")
        self.assertEqual(rank(setting="ex vivo").options[0].vector, "lenti")

    def test_car_is_assumed_ex_vivo(self):
        res = rank("car")
        self.assertEqual(res.setting, "ex vivo")
        self.assertEqual(res.options[0].vector, "lenti")


class Invariants(unittest.TestCase):
    def test_every_ranked_set_is_well_formed(self):
        for strategy in STRATEGIES:
            for codons in (100, 394, 1000, 2000, 3700):
                for setting in ("", "in vivo", "ex vivo"):
                    with self.subTest(strategy=strategy, codons=codons, setting=setting):
                        res = rank(strategy, codons, setting, target="ACGTACGTACGTACGTACGT")
                        self.check(res)

    def check(self, res):
        opts = res.options
        self.assertEqual([o.rank for o in opts], list(range(1, len(opts) + 1)))
        self.assertEqual(len({o.key for o in opts}), len(opts))
        fits = [o.feasible for o in opts]
        self.assertEqual(fits, sorted(fits, reverse=True))          # options that fit come first
        feasible = [o for o in opts if o.feasible]
        self.assertEqual([o.total for o in feasible], sorted((o.total for o in feasible), reverse=True))
        for o in opts:
            self.assertEqual(sum(s.maximum for s in o.scores), 100)
            self.assertEqual(o.total, sum(s.points for s in o.scores))
            self.assertEqual(o.applicability, o.total if o.feasible else 0)
            self.assertEqual(o.cassette, cassette(o.design))
            self.assertEqual([s.key for s in o.scores], [c[0] for c in CRITERIA])
        late = [o for o in opts if not o.feasible]
        self.assertEqual([o.cassette / o.limit for o in late], sorted(o.cassette / o.limit for o in late))

    def test_ranking_is_deterministic(self):
        a, b = rank(), rank()
        self.assertEqual([o.key for o in a.options], [o.key for o in b.options])
        self.assertEqual([o.total for o in a.options], [o.total for o in b.options])

    def test_a_modification_never_makes_a_fitting_option_stop_fitting(self):
        for strategy in STRATEGIES:
            for codons in (394, 1000, 2000):
                plain = {o.key: o for o in rank(strategy, codons, "in vivo", False).options}
                for o in rank(strategy, codons, "in vivo", True).options:
                    if plain[o.key].feasible:
                        self.assertTrue(o.feasible, (strategy, codons, o.key))


class KnownAnswers(unittest.TestCase):
    def test_a_3_kb_insert_does_not_fit_aav_with_the_strong_layout_but_the_compact_one_does(self):
        by = {o.key: o for o in rank(codons=1000, setting="in vivo", optimize=False).options}
        self.assertFalse(by["aav|strong"].feasible)
        self.assertTrue(by["aav|compact"].feasible)
        self.assertEqual(by["aav|strong"].applicability, 0)
        self.assertEqual(by["aav|strong"].rank, len(by))
        tag, line = compare_line(by["aav|compact"], by["aav|strong"])
        self.assertIn("fits while", line)

    def test_a_6_kb_insert_leaves_only_lentivirus(self):
        res = rank(codons=2000, setting="in vivo", optimize=False)
        self.assertEqual(res.options[0].vector, "lenti")
        self.assertTrue(all(not o.feasible for o in res.options if o.vector == "aav"))

    def test_nothing_fits_an_11_kb_insert_and_the_summary_says_so(self):
        res = rank(codons=3700)
        self.assertFalse(any(o.feasible for o in res.options))
        self.assertEqual(res.options[0].vector, "lenti")
        self.assertIn("No option fits", ranking_summary(res))

    def test_a_suicide_gene_prefers_the_tumour_selective_promoter(self):
        res = rank("sui", 394, "in vivo")
        self.assertEqual(res.options[0].arrangement, "selective")

    def test_a_guide_has_the_single_u6_layout_per_vector(self):
        res = rank("guide", target="ACGTACGTACGTACGTACGT")
        self.assertEqual({o.arrangement for o in res.options}, {"standard"})
        self.assertEqual(len(res.options), 2)

    def test_car_adds_a_t_cell_promoter_layout(self):
        self.assertIn("tcell", {o.arrangement for o in rank("car").options})


class InsertProfile(unittest.TestCase):
    def test_no_sequence_means_a_typical_length_and_an_info_note(self):
        p = profile_insert(Config(strategy="rep", gene="TP53"))
        self.assertFalse(p.known)
        self.assertEqual(p.length, 1182)
        self.assertTrue(any(level == "info" for level, _ in p.issues))

    def test_lentiviral_motifs_are_found_and_reported_by_position(self):
        seq = "ATG" + "GCC" * 20 + "AATAAA" + "GCC" * 20 + "CAGGTAAGT" + "GCC" * 20
        p = profile_insert(Config(strategy="rep", gene="X", sequence=seq))
        self.assertEqual(p.polya, [64])
        self.assertEqual(len(p.donors), 1)

    def test_lentivirus_pays_for_cryptic_sites_and_aav_does_less(self):
        clean = rank(setting="ex vivo")
        dirty = rank_vector_options(
            Config(strategy="rep", gene="X", sequence="ATG" + "GCC" * 20 + "AATAAA" + "GCC" * 370), "ex vivo")
        for vec in ("lenti", "aav"):
            a = next(o for o in clean.options if o.vector == vec).score("insert").points
            b = next(o for o in dirty.options if o.vector == vec).score("insert").points
            self.assertLess(b, a)
        lenti_drop = (next(o for o in clean.options if o.vector == "lenti").score("insert").points
                      - next(o for o in dirty.options if o.vector == "lenti").score("insert").points)
        aav_drop = (next(o for o in clean.options if o.vector == "aav").score("insert").points
                    - next(o for o in dirty.options if o.vector == "aav").score("insert").points)
        self.assertGreater(lenti_drop, aav_drop)


class Modifications(unittest.TestCase):
    def test_every_modification_is_measured_on_the_rebuilt_design(self):
        for strategy in STRATEGIES:
            for codons in (100, 700, 1000, 1100, 3700):
                for extra in ({}, dict(marker="EGFP", link="IRES", tag="HA", polya="bGH",
                                       resistance="AmpR", kozak=False)):
                    res = rank(strategy, codons, "in vivo", target="ACGTACGTACGTACGTACGT", **extra)
                    for o in res.options:
                        with self.subTest(strategy=strategy, codons=codons, key=o.key):
                            self.assertEqual(sum(m.bp for m in o.applied_mods), o.cassette - o.base_cassette)
                            self.assertTrue(all(m.bp == 0 for m in o.mods if not m.applied))
                            self.assertEqual(cassette(build_design(o.config)), o.cassette)

    def test_in_vivo_therapy_loses_its_marker_tag_and_gains_the_leaner_parts(self):
        res = rank(codons=700, setting="in vivo", marker="EGFP", link="IRES", tag="HA", kozak=False,
                   polya="bGH")
        o = {x.key: x for x in res.options}["aav|strong"]
        keys = [m.key for m in o.applied_mods]
        self.assertEqual(keys[:2], ["kozak", "polya"])
        self.assertIn("tag", keys)
        self.assertIn("marker", keys)
        self.assertNotIn("link", keys)           # the marker goes, so it is not re-linked first
        self.assertEqual(keys[-1], "resistance")
        c = o.config
        self.assertEqual((c.marker, c.tag, c.polya, c.resistance), ("none", "none", "SV40", "KanR"))

    def test_ex_vivo_keeps_a_marker_but_links_it_with_p2a(self):
        res = rank(codons=700, setting="ex vivo", marker="EGFP", link="IRES")
        o = {x.key: x for x in res.options}["lenti|strong"]
        self.assertEqual((o.config.marker, o.config.link), ("EGFP", "P2A"))
        kept = [m for m in o.mods if m.key == "marker"]
        self.assertEqual(len(kept), 1)
        self.assertFalse(kept[0].applied)

    def test_research_strategies_keep_their_marker_and_ampicillin(self):
        for o in rank("model", 700, "in vivo").options:
            self.assertIn(o.config.marker, (None, "Puro"))
            self.assertEqual(o.config.resistance, "AmpR")

    def test_aav_near_capacity_switches_to_the_short_wpre(self):
        o = {x.key: x for x in rank(codons=1000, setting="in vivo").options}["aav|strong"]
        self.assertTrue(o.config.wpre_short)
        self.assertTrue(o.feasible)
        self.assertEqual(o.base_cassette - o.cassette, 589 - WPRE_SHORT_LEN)
        short = next(m for m in o.applied_mods if m.key == "wpre-short")
        self.assertTrue(short.tradeoff)

    def test_wpre_is_dropped_only_when_that_is_what_makes_it_fit(self):
        o = {x.key: x for x in rank(codons=1100, setting="in vivo").options}["aav|strong"]
        self.assertTrue(o.feasible)
        self.assertIs(o.config.wpre, False)
        self.assertIn("wpre-drop", [m.key for m in o.applied_mods])

    def test_lentivirus_is_never_given_the_short_wpre(self):
        for o in rank(codons=1000, setting="ex vivo").options:
            if o.vector == "lenti":
                self.assertFalse(o.config.wpre_short)

    def test_cryptic_sites_are_listed_as_not_applied_for_lentivirus(self):
        seq = "ATG" + "GCC" * 20 + "AATAAA" + "GCC" * 370
        res = rank_vector_options(Config(strategy="rep", gene="X", sequence=seq), "ex vivo")
        lenti = next(o for o in res.options if o.vector == "lenti")
        recode = [m for m in lenti.mods if m.key == "recode"]
        self.assertEqual(len(recode), 1)
        self.assertFalse(recode[0].applied)

    def test_switching_modifications_off_gives_the_plain_layouts(self):
        res = rank(codons=700, setting="in vivo", optimize=False)
        self.assertFalse(res.optimized)
        for o in res.options:
            self.assertEqual(o.mods, [])
            self.assertEqual(o.base_cassette, o.cassette)
        self.assertIn("switched off", modification_parts(res, 0)[1][1])

    def test_a_guide_keeps_its_selection_marker(self):
        res = rank("guide", 394, "ex vivo", target="ACGTACGTACGTACGTACGT")
        for o in res.options:
            self.assertNotIn("marker", [m.key for m in o.applied_mods])


class WpreShortInTheDesign(unittest.TestCase):
    def test_short_wpre_is_247_bp_and_named(self):
        d = build_design(Config(vector="aav", wpre=True, wpre_short=True))
        f = next(f for f in d.features if f.id == "wpre")
        self.assertEqual((f.name, f.length), ("WPRE3 (short)", WPRE_SHORT_LEN))

    def test_the_flag_does_nothing_without_wpre(self):
        d = build_design(Config(vector="aav", wpre=False, wpre_short=True))
        self.assertNotIn("wpre", [f.id for f in d.features])

    def test_features_still_tile_the_plasmid(self):
        d = build_design(Config(vector="aav", wpre=True, wpre_short=True))
        pos = 1
        for f in d.features:
            self.assertEqual(f.start, pos)
            pos = f.end + 1
        self.assertEqual(pos - 1, d.total)


class Text(unittest.TestCase):
    def setUp(self):
        self.res = rank(codons=1000, setting="in vivo")

    def test_every_option_has_an_explanation_that_names_it(self):
        for i, o in enumerate(self.res.options):
            parts = explain_option_parts(self.res, i)
            self.assertEqual(parts[0][0], "title")
            self.assertIn(o.title, parts[0][1])
            self.assertTrue({tag for tag, _ in parts} <= {"title", "head", "text", "good", "bad", "mute"})

    def test_each_option_is_compared_with_every_other_option(self):
        text = explain_option(self.res, 0)
        self.assertEqual(sum(1 for line in text.splitlines() if line.startswith(("Above", "Below"))),
                         len(self.res.options) - 1)

    def test_modification_text_lists_applied_and_skipped_changes(self):
        o = {x.key: x for x in self.res.options}["aav|strong"]
        text = " ".join(t for _, t in modification_parts(self.res, self.res.options.index(o)))
        self.assertIn("Use the short WPRE3", text)
        self.assertIn("Considered but not applied", text)

    def test_the_report_covers_every_option_and_how_scores_work(self):
        report = options_report(self.res)
        for o in self.res.options:
            self.assertIn(o.title, report)
        self.assertIn("How the score works", report)
        self.assertIn(options_table(self.res).splitlines()[0], report)

    def test_the_text_is_ascii_so_every_console_can_print_it(self):
        self.assertTrue(options_report(self.res).isascii())

    def test_fit_text_tells_the_three_cases_apart(self):
        by = {o.key: o for o in rank(codons=2000, setting="in vivo", optimize=False).options}
        self.assertTrue(fit_text(by["lenti|compact"]).startswith("fits"))
        self.assertTrue(fit_text(by["aav|strong"]).startswith("too large"))

    def test_parse_explain(self):
        self.assertEqual(parse_explain("all", 3), [0, 1, 2])
        self.assertEqual(parse_explain("2", 3), [1])
        for bad in ("0", "4", "x"):
            with self.assertRaises(ValueError):
                parse_explain(bad, 3)


class ApplyingAnOption(unittest.TestCase):
    def test_config_for_option_rebuilds_the_prebuilt_design_exactly(self):
        cfg = Config(strategy="rep", gene="X", sequence=orf(700), marker="EGFP", tag="HA")
        res = rank_vector_options(cfg, "in vivo")
        for o in res.options:
            built = build_design(config_for_option(cfg, res, o))
            self.assertEqual(built.name, o.design.name)
            self.assertEqual(cassette(built), o.cassette)
            self.assertEqual(built.total, o.design.total)

    def test_the_basis_records_the_option_and_its_modifications(self):
        cfg = Config(strategy="rep", gene="X", sequence=orf(700), marker="EGFP")
        res = rank_vector_options(cfg, "in vivo")
        o = res.options[0]
        text = basis_for_option(res, o)
        self.assertIn("Vector option #1", text)
        self.assertIn("Modifications:", text)
        self.assertTrue(o.design.basis.startswith("Vector option #1"))
        self.assertIn("Design basis", __import__("cvd.export.summary", fromlist=["x"]).to_summary(o.design))

    def test_an_earlier_basis_is_kept(self):
        cfg = Config(strategy="rep", gene="X", sequence=orf(100), basis="Aneuploidy analysis: 14 of 39")
        res = rank_vector_options(cfg)
        self.assertTrue(config_for_option(cfg, res, res.options[0]).basis.startswith("Aneuploidy analysis"))

    def test_candidates_ignore_the_users_element_overrides(self):
        cfg = Config(strategy="rep", gene="X", sequence=orf(100), overrides={"prom": "ACGT"})
        res = rank_vector_options(cfg)
        self.assertTrue(all(not o.config.overrides for o in res.options))


if __name__ == "__main__":
    unittest.main()
