import unittest

from cvd.cnv import (analyze_aneuploidy, basis_text, config_from_recommendation,
                     example_cnv_text, parse_segments, recommend)
from cvd.model import Config
from cvd.parts import STRATEGIES


def analyse_example():
    segments, vtype, notes, sample = parse_segments(example_cnv_text())
    return analyze_aneuploidy(segments, value_type=vtype, sample=sample, notes=notes)


class Suggestions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = analyse_example()
        cls.recs = recommend(cls.res)

    def pairs(self):
        return [(r.strategy, r.gene) for r in self.recs]

    def test_a_flag_payload_is_always_offered(self):
        self.assertEqual(self.pairs()[0], ("flag", "ULBP2"))
        self.assertIn(("flag", "tEGFR (surface tag)"), self.pairs())

    def test_mhc_loss_is_cited_in_the_nkg2d_rationale(self):
        reason = self.recs[0].reason
        self.assertIn("6p", reason)
        self.assertIn("MHC-I", reason)

    def test_3q_gain_pairs_calreticulin_with_cd47_knockdown(self):
        self.assertIn(("flag", "CALR"), self.pairs())
        self.assertIn(("guide", "CD47"), self.pairs())

    def test_high_score_adds_cytokine_support(self):
        self.assertIn(("cyt", "IL-12 (single chain)"), self.pairs())

    def test_lost_arms_suggest_replacing_their_tumour_suppressors(self):
        self.assertIn(("rep", "TP53"), self.pairs())
        self.assertIn(("rep", "PTEN"), self.pairs())

    def test_gained_arms_suggest_knocking_down_their_oncogenes(self):
        self.assertIn(("guide", "MYC"), self.pairs())

    def test_at_most_two_suppressors_and_two_oncogenes_are_offered(self):
        self.assertEqual(len([r for r in self.recs if r.strategy == "rep"]), 2)
        self.assertEqual(len([r for r in self.recs if r.strategy == "guide"
                              and r.gene not in ("CD47", "CD274")]), 2)

    def test_no_duplicate_strategy_gene_pairs(self):
        self.assertEqual(len(self.pairs()), len(set(self.pairs())))

    def test_every_suggestion_names_a_real_strategy_and_carries_a_reason(self):
        for rec in self.recs:
            self.assertIn(rec.strategy, STRATEGIES, rec.gene)
            self.assertEqual(rec.promoter, STRATEGIES[rec.strategy]["promoter"])
            self.assertTrue(rec.title and rec.reason, rec.gene)

    def test_oversized_suppressors_are_never_suggested(self):
        for gene in ("APC", "NF1", "BRCA1"):
            self.assertNotIn(("rep", gene), self.pairs())

    def test_pdl1_knockdown_needs_9p_gain_which_the_example_lacks(self):
        self.assertNotIn(("guide", "CD274"), self.pairs())

    def test_a_quiet_genome_still_offers_the_mhc_independent_flags_only(self):
        quiet = analyze_aneuploidy(parse_segments("1\t1\t120000000\t0.0\n")[0],
                                   value_type="log2", center=False)
        self.assertEqual([(r.strategy, r.gene) for r in recommend(quiet)],
                         [("flag", "ULBP2"), ("flag", "tEGFR (surface tag)")])


class Basis(unittest.TestCase):
    def test_basis_text_records_the_score_and_the_chosen_design(self):
        res = analyse_example()
        rec = recommend(res)[0]
        text = basis_text(res, rec)
        self.assertIn("score 14", text)
        self.assertIn("7 gains, 7 losses", text)
        self.assertIn(rec.title, text)
        self.assertIn(rec.reason, text)

    def test_sample_name_is_included_when_known(self):
        res = analyse_example()
        res.sample = "TCGA-XX"
        self.assertIn("TCGA-XX", basis_text(res, recommend(res)[0]))


class ApplyToConfig(unittest.TestCase):
    def setUp(self):
        self.rec = recommend(analyse_example())[0]

    def test_strategy_gene_promoter_and_basis_are_applied(self):
        cfg = config_from_recommendation(Config(), self.rec, "because")
        self.assertEqual((cfg.strategy, cfg.gene, cfg.basis), (self.rec.strategy,
                                                               self.rec.gene, "because"))
        self.assertEqual(cfg.promoter, self.rec.promoter)

    def test_an_explicit_promoter_can_be_kept(self):
        cfg = config_from_recommendation(Config(promoter="CMV"), self.rec, "because",
                                         keep_promoter=True)
        self.assertEqual(cfg.promoter, "CMV")


if __name__ == "__main__":
    unittest.main()
