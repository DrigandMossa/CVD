import unittest

from cvd.cnv import analyze_aneuploidy, aneuploidy_level, example_cnv_text, parse_segments
from cvd.model import Segment
from cvd.parts import ACROCENTRIC_P, BUILDS, N_SCORED_ARMS

CHR1_CEN = BUILDS["hg38"]["1"][1]


def example(build="hg38", **kwargs):
    segments, vtype, notes, sample = parse_segments(example_cnv_text(build))
    options = dict(build=build, value_type=vtype, center=True, sample=sample, notes=notes)
    options.update(kwargs)
    return analyze_aneuploidy(segments, **options)


def call_of(result, arm):
    return next(a.call for a in result.arms if a.arm == arm)


class SyntheticExample(unittest.TestCase):
    """The shipped example encodes a known answer, so it pins the whole pipeline."""

    @classmethod
    def setUpClass(cls):
        cls.res = example()

    def test_score_is_the_number_of_called_arms(self):
        self.assertEqual((self.res.gains, self.res.losses), (7, 7))
        self.assertEqual(self.res.score, 14)

    def test_gained_and_lost_arms_match_the_encoded_values(self):
        self.assertEqual([a.arm for a in self.res.arms if a.call == "gain"],
                         ["1q", "3q", "5p", "7p", "7q", "8q", "20q"])
        self.assertEqual([a.arm for a in self.res.arms if a.call == "loss"],
                         ["6p", "9p", "10q", "13q", "15q", "17p", "18q"])

    def test_partly_altered_arm_is_not_called(self):
        self.assertEqual(call_of(self.res, "4q"), "partial")

    def test_both_arms_altered_reports_a_whole_chromosome_event(self):
        self.assertEqual(self.res.whole_chromosomes, ["7 gain"])

    def test_acrocentric_p_arms_are_never_scored(self):
        for arm in ACROCENTRIC_P:
            self.assertEqual(call_of(self.res, arm), "not scored", arm)

    def test_every_arm_is_reported_exactly_once(self):
        arms = [a.arm for a in self.res.arms]
        self.assertEqual(len(arms), 44)
        self.assertEqual(len(set(arms)), 44)

    def test_fraction_of_genome_altered_is_a_percentage(self):
        self.assertAlmostEqual(self.res.fga, 40.3, places=0)

    def test_immune_findings_are_raised_for_altered_annotated_arms(self):
        headlines = [f[0] for f in self.res.findings]
        self.assertIn("6p loss: HLA class I region (6p21.3)", headlines)
        self.assertIn("3q gain: CD47 (3q13.12)", headlines)
        self.assertFalse(any(h.startswith("9p gain") for h in headlines))

    def test_log2_input_cannot_claim_whole_genome_doubling(self):
        self.assertFalse(self.res.wgd_suspected)
        self.assertTrue(any("doubling" in n for n in self.res.notes))

    def test_the_other_build_gives_the_same_calls(self):
        other = example("hg19")
        self.assertEqual(other.score, self.res.score)


class ArmCallRule(unittest.TestCase):
    """An arm is called at 80% of its covered length, neutral at 20% or less."""

    def arm_with(self, altered_fraction):
        cut = int(altered_fraction * CHR1_CEN)
        segments = [Segment("1", 0, cut, 0.5), Segment("1", cut, CHR1_CEN, 0.0)]
        return call_of(analyze_aneuploidy(segments, value_type="log2", center=False), "1p")

    def test_exactly_eighty_percent_is_called(self):
        self.assertEqual(self.arm_with(0.80), "gain")

    def test_just_under_eighty_percent_is_partial(self):
        self.assertEqual(self.arm_with(0.79), "partial")

    def test_twenty_percent_or_less_is_neutral(self):
        self.assertEqual(self.arm_with(0.20), "neutral")

    def test_just_over_twenty_percent_is_partial(self):
        self.assertEqual(self.arm_with(0.21), "partial")

    def test_thin_coverage_is_reported_as_no_data(self):
        segments = [Segment("1", 0, int(0.4 * CHR1_CEN), 0.5)]
        self.assertEqual(call_of(analyze_aneuploidy(segments, value_type="log2",
                                                    center=False), "1p"), "no data")


class Baselines(unittest.TestCase):
    def test_copy_number_input_centres_on_the_genome_median(self):
        segments = [Segment("1", 0, CHR1_CEN, 4.0),
                    Segment("2", 0, BUILDS["hg38"]["2"][0], 2.0)]
        result = analyze_aneuploidy(segments, value_type="cn")
        self.assertIn("median", result.baseline_desc)
        self.assertEqual(call_of(result, "1p"), "gain")

    def test_explicit_ploidy_overrides_the_median(self):
        segments = [Segment("1", 0, CHR1_CEN, 4.0), Segment("2", 0, 90000000, 4.0)]
        result = analyze_aneuploidy(segments, value_type="cn", ploidy=2.0)
        self.assertIn("given", result.baseline_desc)
        self.assertEqual(call_of(result, "1p"), "gain")

    def test_high_baseline_copy_number_suspects_genome_doubling(self):
        segments = [Segment(c, 0, 40000000, 4.0) for c in ("1", "2", "3")]
        self.assertTrue(analyze_aneuploidy(segments, value_type="cn").wgd_suspected)

    def test_centring_log2_ratios_can_change_the_calls(self):
        segments = [Segment("1", 0, CHR1_CEN, 0.5), Segment("2", 0, 90000000, 0.5)]
        self.assertEqual(call_of(analyze_aneuploidy(segments, value_type="log2",
                                                    center=False), "1p"), "gain")
        self.assertEqual(call_of(analyze_aneuploidy(segments, value_type="log2",
                                                    center=True), "1p"), "neutral")

    def test_sparse_input_notes_that_the_score_is_a_lower_bound(self):
        result = analyze_aneuploidy([Segment("1", 0, CHR1_CEN, 0.5)], value_type="log2",
                                    center=False)
        self.assertTrue(any("lower bound" in n for n in result.notes))
        self.assertTrue(any(str(N_SCORED_ARMS) in n for n in result.notes))


class Validation(unittest.TestCase):
    def test_thresholds_must_straddle_zero(self):
        for gain, loss in ((0.2, 0.1), (-0.2, -0.3), (0.0, -0.2)):
            with self.assertRaises(ValueError):
                analyze_aneuploidy([Segment("1", 0, 1000, 0.0)], gain=gain, loss=loss)

    def test_unknown_build_is_rejected(self):
        with self.assertRaises(ValueError):
            analyze_aneuploidy([Segment("1", 0, 1000, 0.0)], build="hg18")

    def test_ploidy_must_be_positive(self):
        for ploidy in (0, -1.0):
            with self.assertRaises(ValueError):
                analyze_aneuploidy([Segment("1", 0, CHR1_CEN, 2.0)], value_type="cn",
                                   ploidy=ploidy)

    def test_segments_outside_the_build_are_an_error(self):
        with self.assertRaises(ValueError) as caught:
            analyze_aneuploidy([Segment("X", 0, 1000, 0.5)])
        self.assertIn("build", str(caught.exception))


class Levels(unittest.TestCase):
    def test_bins(self):
        self.assertEqual(aneuploidy_level(0), "low")
        self.assertEqual(aneuploidy_level(2), "low")
        self.assertEqual(aneuploidy_level(3), "intermediate")
        self.assertEqual(aneuploidy_level(9), "intermediate")
        self.assertEqual(aneuploidy_level(10), "high")


if __name__ == "__main__":
    unittest.main()
