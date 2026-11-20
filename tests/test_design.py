import unittest

from cvd.design import build_design
from cvd.model import Config
from cvd.parts import SGRNA_SCAFFOLD, SHRNA_LOOP
from cvd.seqops import revcomp

ORF = "ATG" + "GCG" * 10          # 33 bp, in frame, no internal stop
TARGET = "GTTGGAGCTGATGGCGTAGG"   # 20 nt, no TTTT


def ids(design):
    return [f.id for f in design.features]


def levels(design, level):
    return [text for lv, text in design.messages if lv == level]


class Backbone(unittest.TestCase):
    def test_lentivirus_has_ltrs_and_packaging_elements(self):
        got = ids(build_design(Config(vector="lenti", gene="TP53")))
        self.assertEqual(got[:4], ["ltr5", "psi", "rre", "cppt"])
        self.assertIn("ltr3", got)

    def test_aav_has_itrs_and_no_ltrs(self):
        got = ids(build_design(Config(vector="aav", gene="TP53")))
        self.assertEqual(got[0], "itr5")
        self.assertIn("itr3", got)
        self.assertNotIn("ltr5", got)

    def test_plasmid_has_no_viral_elements_and_no_wpre(self):
        got = ids(build_design(Config(vector="plasmid", gene="TP53")))
        for absent in ("ltr5", "ltr3", "itr5", "itr3", "f1", "wpre"):
            self.assertNotIn(absent, got)

    def test_every_vector_carries_bacterial_selection_and_origin(self):
        for vector in ("lenti", "aav", "plasmid"):
            got = ids(build_design(Config(vector=vector, gene="TP53")))
            self.assertEqual(got[-2:], ["res", "ori"], vector)


class Layout(unittest.TestCase):
    def test_features_tile_the_sequence_without_gaps_or_overlaps(self):
        design = build_design(Config(gene="TP53", sequence=ORF, marker="EGFP"))
        position = 1
        for feature in design.features:
            self.assertEqual(feature.start, position, feature.id)
            self.assertEqual(feature.end, feature.start + feature.length - 1, feature.id)
            position = feature.end + 1
        self.assertEqual(position - 1, design.total)

    def test_total_equals_sum_of_feature_lengths(self):
        design = build_design(Config(gene="TP53", sequence=ORF))
        self.assertEqual(design.total, sum(f.length for f in design.features))

    def test_assembled_sequence_is_the_concatenation_of_features(self):
        design = build_design(Config(gene="TP53", sequence=ORF))
        self.assertEqual(design.seq, "".join(f.seq for f in design.features))

    def test_missing_sequences_become_n_placeholders(self):
        design = build_design(Config(gene="TP53"))
        goi = next(f for f in design.features if f.id == "goi")
        self.assertEqual(goi.seq, "N" * goi.length)
        self.assertIn("placeholder", " ".join(levels(design, "warn")))

    def test_cassette_flag_covers_only_elements_between_the_itrs(self):
        design = build_design(Config(vector="aav", gene="TP53"))
        inside = [f.id for f in design.features if f.in_cassette]
        self.assertNotIn("itr5", inside)
        self.assertNotIn("itr3", inside)
        self.assertIn("goi", inside)
        self.assertIn("prom", inside)


class Cassettes(unittest.TestCase):
    def test_shrna_hairpin_is_sense_loop_antisense(self):
        design = build_design(Config(strategy="guide", gene="KRAS", target=TARGET))
        guide = next(f for f in design.features if f.id == "guide")
        self.assertEqual(guide.seq, TARGET + SHRNA_LOOP + revcomp(TARGET))
        self.assertEqual(guide.length, len(guide.seq))

    def test_sgrna_is_protospacer_plus_scaffold(self):
        design = build_design(Config(strategy="guide", gene="KRAS", target=TARGET,
                                     guide_type="sgRNA"))
        guide = next(f for f in design.features if f.id == "guide")
        self.assertEqual(guide.seq, TARGET + SGRNA_SCAFFOLD)

    def test_guide_strategy_uses_u6_and_a_pol_iii_terminator(self):
        got = ids(build_design(Config(strategy="guide", gene="KRAS", target=TARGET)))
        self.assertIn("u6", got)
        self.assertIn("term", got)
        self.assertNotIn("prom", got)
        self.assertNotIn("kozak", got)

    def test_guide_marker_gets_its_own_pgk_promoter(self):
        got = ids(build_design(Config(strategy="guide", gene="KRAS", target=TARGET,
                                      marker="Puro")))
        self.assertEqual(got[got.index("term") + 1:got.index("term") + 3], ["pgk", "mark2"])

    def test_car_chains_the_domains_in_frame_order(self):
        got = ids(build_design(Config(strategy="car", gene="CD19")))
        self.assertEqual([i for i in got if i in ("lead", "goi", "hinge", "costim", "cd3z")],
                         ["lead", "goi", "hinge", "costim", "cd3z"])

    def test_car_without_costimulatory_domain_omits_it(self):
        got = ids(build_design(Config(strategy="car", gene="CD19", costim="none")))
        self.assertNotIn("costim", got)

    def test_p2a_marker_puts_the_single_stop_after_the_reporter(self):
        got = ids(build_design(Config(gene="TP53", marker="EGFP", link="P2A")))
        self.assertEqual(got.index("stop"), got.index("mark") + 1)
        self.assertEqual(got.count("stop"), 1)

    def test_ires_marker_keeps_the_stop_before_the_linker(self):
        got = ids(build_design(Config(gene="TP53", marker="EGFP", link="IRES")))
        self.assertLess(got.index("stop"), got.index("ires"))

    def test_kozak_can_be_switched_off(self):
        self.assertNotIn("kozak", ids(build_design(Config(gene="TP53", kozak=False))))


class Overrides(unittest.TestCase):
    def test_override_replaces_length_and_marks_the_feature_custom(self):
        design = build_design(Config(gene="TP53", overrides={"prom": "acgt acgt 12"}))
        prom = next(f for f in design.features if f.id == "prom")
        self.assertEqual(prom.seq, "ACGTACGT")
        self.assertEqual(prom.length, 8)
        self.assertTrue(prom.custom)

    def test_built_in_element_cannot_be_overridden(self):
        with self.assertRaises(ValueError) as caught:
            build_design(Config(gene="TP53", overrides={"goi": "ACGT"}))
        self.assertIn("goi", str(caught.exception))

    def test_unknown_element_id_is_rejected_and_lists_valid_ids(self):
        with self.assertRaises(ValueError) as caught:
            build_design(Config(gene="TP53", overrides={"nope": "ACGT"}))
        self.assertIn("prom", str(caught.exception))

    def test_override_without_nucleotides_is_rejected(self):
        with self.assertRaises(ValueError):
            build_design(Config(gene="TP53", overrides={"prom": "XYZ"}))


class Capacity(unittest.TestCase):
    def test_oversized_aav_cassette_is_an_error(self):
        design = build_design(Config(vector="aav", gene="APC"))
        self.assertTrue(any("packaging limit" in t for t in levels(design, "error")))

    def test_aav_cassette_near_the_limit_is_a_warning(self):
        design = build_design(Config(vector="aav", gene="EGFR", promoter="EFS"))
        self.assertTrue(any("close to" in t for t in levels(design, "warn")))

    def test_small_aav_cassette_is_only_a_note(self):
        design = build_design(Config(vector="aav", gene="CDKN2A"))
        self.assertFalse([t for t in levels(design, "error") if "limit" in t])
        self.assertTrue(any("limit about 4.7 kb" in t for t in levels(design, "info")))

    def test_large_lentiviral_insert_warns_about_titer(self):
        design = build_design(Config(vector="lenti", gene="APC"))
        self.assertTrue(any("titer" in t for t in levels(design, "warn")))


class Naming(unittest.TestCase):
    def test_generated_name_encodes_vector_promoter_gene_tag_and_marker(self):
        design = build_design(Config(gene="TP53", tag="HA", marker="EGFP", link="IRES"))
        self.assertEqual(design.name, "pLV-EF1a>TP53-HA-IRES-EGFP")

    def test_guide_names_use_the_u6_promoter_and_a_sh_or_sg_prefix(self):
        self.assertTrue(build_design(Config(strategy="guide", gene="KRAS", target=TARGET,
                                            guide_type="sgRNA")).name.startswith("pLV-U6>sgKRAS"))

    def test_variant_label_is_appended_to_the_gene(self):
        self.assertIn("BRAF_V600E", build_design(Config(strategy="model", gene="BRAF",
                                                        variant="V600E")).name)

    def test_explicit_name_wins(self):
        self.assertEqual(build_design(Config(gene="TP53", name="  pMyVector ")).name,
                         "pMyVector")


class Validation(unittest.TestCase):
    def test_unknown_vector_is_rejected(self):
        with self.assertRaises(ValueError):
            build_design(Config(vector="retro", gene="TP53"))

    def test_unknown_strategy_is_rejected(self):
        with self.assertRaises(ValueError):
            build_design(Config(strategy="nope", gene="TP53"))

    def test_messages_are_ordered_error_then_warning_then_note(self):
        design = build_design(Config(vector="aav", gene="APC", sequence="ATGTAACGT"))
        rank = [("error", "warn", "info").index(lv) for lv, _ in design.messages]
        self.assertEqual(rank, sorted(rank))


if __name__ == "__main__":
    unittest.main()
