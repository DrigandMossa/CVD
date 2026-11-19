import unittest

from cvd.seqops import clean_sequence, count_overlapping, gc_percent, pick, revcomp


class CleanSequence(unittest.TestCase):
    def test_strips_fasta_header_whitespace_and_digits(self):
        seq, removed = clean_sequence(">chr1 TP53 ORF\nacgt 123\nACGT\n")
        self.assertEqual(seq, "ACGTACGT")
        self.assertEqual(removed, 0)

    def test_strips_multiple_fasta_headers(self):
        seq, _ = clean_sequence(">one\nACGT\n>two\nTTTT\n")
        self.assertEqual(seq, "ACGTTTTT")

    def test_uracil_becomes_thymine(self):
        self.assertEqual(clean_sequence("AUGC")[0], "ATGC")

    def test_counts_removed_non_nucleotides(self):
        seq, removed = clean_sequence("ACGTXYZN")
        self.assertEqual(seq, "ACGTN")
        self.assertEqual(removed, 3)

    def test_empty_and_none(self):
        self.assertEqual(clean_sequence(""), ("", 0))
        self.assertEqual(clean_sequence(None), ("", 0))


class Helpers(unittest.TestCase):
    def test_revcomp(self):
        self.assertEqual(revcomp("ATGCN"), "NGCAT")

    def test_revcomp_is_an_involution(self):
        self.assertEqual(revcomp(revcomp("ACGTTGCA")), "ACGTTGCA")

    def test_gc_percent_rounds(self):
        self.assertEqual(gc_percent("GGCC"), 100)
        self.assertEqual(gc_percent("ATAT"), 0)
        self.assertEqual(gc_percent("GCAT"), 50)
        self.assertEqual(gc_percent("GCA"), 67)

    def test_gc_percent_of_empty_is_zero(self):
        self.assertEqual(gc_percent(""), 0)

    def test_count_overlapping_counts_overlaps(self):
        self.assertEqual(count_overlapping("AAAA", "AA"), 3)
        self.assertEqual(count_overlapping("GAATTCGAATTC", "GAATTC"), 2)
        self.assertEqual(count_overlapping("ACGT", "TTTT"), 0)


class Pick(unittest.TestCase):
    def test_case_insensitive_match_returns_canonical_option(self):
        self.assertEqual(pick("ef1A", ["CMV", "EF1a"], "promoter"), "EF1a")

    def test_works_on_dict_keys(self):
        self.assertEqual(pick("AMPR", {"AmpR": 1, "KanR": 2}, "selection"), "AmpR")

    def test_unknown_value_lists_the_options(self):
        with self.assertRaises(ValueError) as caught:
            pick("T7", ["CMV", "EF1a"], "promoter")
        message = str(caught.exception)
        self.assertIn("promoter", message)
        self.assertIn("CMV, EF1a", message)


if __name__ == "__main__":
    unittest.main()
