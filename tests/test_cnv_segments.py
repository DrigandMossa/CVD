import unittest

from cvd.cnv import norm_chrom, parse_segments, weighted_median

HEADER = "chromosome\tloc.start\tloc.end\tseg.mean\nchr1\t1\t1000\t0.5\nchr2\t1\t1000\t-0.4\n"


class ChromosomeNames(unittest.TestCase):
    def test_strips_chr_prefix_case_insensitively(self):
        self.assertEqual(norm_chrom("chr7"), "7")
        self.assertEqual(norm_chrom("CHR7"), "7")
        self.assertEqual(norm_chrom(" 7 "), "7")

    def test_numeric_sex_chromosome_aliases(self):
        self.assertEqual(norm_chrom("23"), "X")
        self.assertEqual(norm_chrom("24"), "Y")
        self.assertEqual(norm_chrom("chrX"), "X")

    def test_mitochondrial_and_unplaced_contigs_are_dropped(self):
        for name in ("chrM", "MT", "chr1_KI270706v1_random", "25"):
            self.assertIsNone(norm_chrom(name), name)


class WeightedMedian(unittest.TestCase):
    def test_weights_dominate_the_value_order(self):
        self.assertEqual(weighted_median([(0.0, 1), (1.0, 100)]), 1.0)

    def test_equal_weights_give_the_middle_value(self):
        self.assertEqual(weighted_median([(1.0, 10), (2.0, 10), (3.0, 10)]), 2.0)

    def test_zero_weights_are_ignored_and_empty_is_zero(self):
        self.assertEqual(weighted_median([(5.0, 0)]), 0.0)
        self.assertEqual(weighted_median([]), 0.0)


class Headers(unittest.TestCase):
    def test_recognises_dnacopy_style_column_aliases(self):
        segments, vtype, _notes, sample = parse_segments(HEADER)
        self.assertEqual(vtype, "log2")
        self.assertEqual(sample, "")
        self.assertEqual([(s.chrom, s.start, s.end, s.value) for s in segments],
                         [("1", 1, 1000, 0.5), ("2", 1, 1000, -0.4)])

    def test_comma_separated_input(self):
        segments, _vt, _n, _s = parse_segments("chrom,start,end,log2\n1,1,500,0.3\n")
        self.assertEqual(len(segments), 1)

    def test_whitespace_separated_input(self):
        segments, _vt, _n, _s = parse_segments("chrom start end log2\n1 1 500 0.3\n")
        self.assertEqual(len(segments), 1)

    def test_comment_lines_are_skipped(self):
        segments, _vt, _n, _s = parse_segments("# a note\n" + HEADER)
        self.assertEqual(len(segments), 2)

    def test_copy_number_column_is_preferred_and_reported(self):
        _s, vtype, _n, _sample = parse_segments("chrom\tstart\tend\tcopy_number\n1\t1\t900\t3\n")
        self.assertEqual(vtype, "cn")

    def test_missing_position_columns_is_an_error_naming_the_header(self):
        with self.assertRaises(ValueError) as caught:
            parse_segments("name\tvalue\nfoo\t1\n")
        self.assertIn("name", str(caught.exception))

    def test_missing_value_column_is_an_error(self):
        with self.assertRaises(ValueError):
            parse_segments("chrom\tstart\tend\n1\t1\t900\n")

    def test_requested_value_type_that_is_absent_is_an_error(self):
        with self.assertRaises(ValueError):
            parse_segments("chrom\tstart\tend\tlog2\n1\t1\t900\t0.2\n", value_type="cn")


class Headerless(unittest.TestCase):
    def test_four_columns_are_chrom_start_end_value(self):
        segments, _vt, _n, _s = parse_segments("1\t1\t1000\t0.5\n")
        self.assertEqual((segments[0].chrom, segments[0].value), ("1", 0.5))

    def test_six_columns_are_read_as_dnacopy_output(self):
        text = "TCGA-A\t1\t1\t1000\t42\t0.5\nTCGA-A\t2\t1\t1000\t37\t-0.5\n"
        segments, vtype, _n, sample = parse_segments(text)
        self.assertEqual(sample, "TCGA-A")
        self.assertEqual(vtype, "log2")
        self.assertEqual(len(segments), 2)

    def test_other_column_counts_are_rejected_with_guidance(self):
        with self.assertRaises(ValueError) as caught:
            parse_segments("1\t1\t1000\n")
        self.assertIn("header", str(caught.exception))


class Samples(unittest.TestCase):
    TWO = ("sample\tchrom\tstart\tend\tlog2\nA\t1\t1\t900\t0.5\nB\t1\t1\t900\t-0.5\n")

    def test_several_samples_without_a_choice_is_an_error_listing_them(self):
        with self.assertRaises(ValueError) as caught:
            parse_segments(self.TWO)
        self.assertIn("A, B", str(caught.exception))

    def test_named_sample_is_selected(self):
        segments, _vt, _n, sample = parse_segments(self.TWO, sample="B")
        self.assertEqual(sample, "B")
        self.assertEqual(segments[0].value, -0.5)

    def test_unknown_sample_is_an_error(self):
        with self.assertRaises(ValueError):
            parse_segments(self.TWO, sample="C")


class Robustness(unittest.TestCase):
    def test_value_type_is_inferred_from_negative_values_and_noted(self):
        _s, vtype, notes, _sample = parse_segments("1\t1\t900\t-0.3\n")
        self.assertEqual(vtype, "log2")
        self.assertTrue(any("log2" in n for n in notes))

    def test_all_positive_values_are_read_as_copy_number(self):
        _s, vtype, _n, _sample = parse_segments("1\t1\t900\t3\n")
        self.assertEqual(vtype, "cn")

    def test_unreadable_rows_are_skipped_and_counted(self):
        _s, _vt, notes, _sample = parse_segments(HEADER + "chr3\tbad\t1000\t0.1\n")
        self.assertTrue(any("skipped" in n for n in notes))

    def test_reversed_coordinates_are_normalised(self):
        segments, _vt, _n, _s = parse_segments("1\t1000\t1\t0.5\n")
        self.assertEqual((segments[0].start, segments[0].end), (1, 1000))

    def test_zero_length_and_non_autosomal_segments_are_dropped(self):
        segments, _vt, _n, _s = parse_segments("1\t500\t500\t0.5\nchrM\t1\t900\t0.5\n"
                                               "2\t1\t900\t0.5\n")
        self.assertEqual([s.chrom for s in segments], ["2"])

    def test_empty_input_is_an_error(self):
        with self.assertRaises(ValueError):
            parse_segments("\n# only a comment\n")

    def test_no_usable_rows_is_an_error(self):
        with self.assertRaises(ValueError):
            parse_segments("chrom\tstart\tend\tlog2\nbad\trow\there\tx\n")


if __name__ == "__main__":
    unittest.main()
