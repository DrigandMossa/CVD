import unittest

from cvd.design import check_insert


def run(seq, label="Insert", vector="lenti", check_atg=True):
    msgs = []
    check_insert(msgs, seq, label, vector, check_atg)
    return msgs


def texts(msgs, level):
    return [t for lv, t in msgs if lv == level]


class ReadingFrame(unittest.TestCase):
    def test_length_not_a_multiple_of_three_warns(self):
        self.assertTrue(any("multiple of 3" in t for t in texts(run("ATGGC"), "warn")))

    def test_in_frame_stop_is_an_error_reporting_the_codon(self):
        errors = texts(run("ATGTAAGCGGCG"), "error")
        self.assertTrue(errors)
        self.assertIn("codon 2", errors[0])

    def test_clean_orf_has_no_errors(self):
        self.assertEqual(texts(run("ATG" + "GCG" * 5), "error"), [])

    def test_final_stop_codon_is_not_counted_as_internal(self):
        self.assertEqual(texts(run("ATGGCGTAA"), "error"), [])

    def test_missing_start_codon_warns_when_checked(self):
        self.assertTrue(any("ATG" in t for t in texts(run("GCGGCGGCG"), "warn")))

    def test_missing_start_codon_is_ignored_for_scfv(self):
        self.assertFalse(any("ATG" in t for t in texts(run("GCGGCG", check_atg=False), "warn")))

    def test_trailing_stop_codon_warns_about_tags_and_p2a(self):
        self.assertTrue(any("ends in a stop codon" in t for t in texts(run("ATGGCGTAA"), "warn")))


class Composition(unittest.TestCase):
    def test_note_reports_length_with_thousands_separator_and_gc(self):
        notes = texts(run("ATG" + "GC" * 999), "info")
        self.assertTrue(any("2,001 bp" in t for t in notes))

    def test_cryptic_polya_warns_only_for_lentivirus(self):
        seq = "ATG" + "AATAAA" + "GCG"
        self.assertTrue(any("polyadenylation" in t for t in texts(run(seq, vector="lenti"), "warn")))
        self.assertFalse(any("polyadenylation" in t for t in texts(run(seq, vector="aav"), "warn")))

    def test_restriction_sites_are_split_into_free_and_cutting_enzymes(self):
        note = [t for t in texts(run("ATGGAATTCGCG"), "info") if "Restriction sites" in t][0]
        self.assertIn("EcoRI x1", note)
        self.assertIn("BamHI", note.split("Cut by")[0])

    def test_insert_free_of_every_listed_enzyme_says_so(self):
        note = [t for t in texts(run("ATG" + "GCG" * 5), "info") if "Restriction sites" in t][0]
        self.assertNotIn("Cut by", note)

    def test_label_appears_in_every_message(self):
        for _lv, text in run("ATGGC", label="scFv"):
            self.assertIn("scFv", text)


if __name__ == "__main__":
    unittest.main()
