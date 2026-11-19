import os
import re
import tempfile
import unittest
import xml.etree.ElementTree as ET

from cvd.cnv import (analyze_aneuploidy, aneuploidy_report, arms_tsv, example_cnv_text,
                     parse_segments, recommend, write_aneuploidy_files)
from cvd.design import build_design
from cvd.export import to_fasta, to_genbank, to_summary, to_svg, write_outputs
from cvd.model import Config
from cvd.parts import KINDS, N_SCORED_ARMS

ORF = "ATG" + "GCG" * 10


def design(**kwargs):
    options = dict(gene="TP53", sequence=ORF, marker="EGFP")
    options.update(kwargs)
    return build_design(Config(**options))


class GenBank(unittest.TestCase):
    def setUp(self):
        self.d = design()
        self.text = to_genbank(self.d)

    def test_locus_line_reports_length_and_circular_dna(self):
        first = self.text.splitlines()[0]
        self.assertTrue(first.startswith("LOCUS"))
        self.assertIn("%d bp" % self.d.total, first)
        self.assertIn("circular", first)

    def test_locus_name_is_truncated_and_free_of_separators(self):
        long_name = to_genbank(design(name="a-very-long-vector-name-indeed"))
        locus = long_name.splitlines()[0].split()[1]
        self.assertLessEqual(len(locus), 16)
        self.assertNotIn(" ", locus)

    def test_every_feature_appears_with_its_label_and_span(self):
        for feature in self.d.features:
            self.assertIn('/label="%s"' % feature.name, self.text)
            self.assertIn("%d..%d" % (feature.start, feature.end), self.text)

    def test_feature_keys_come_from_the_kind_table(self):
        for feature in self.d.features:
            self.assertIn(KINDS[feature.kind][1], self.text)

    def test_origin_block_holds_the_whole_sequence_in_sixty_base_rows(self):
        body = self.text.split("ORIGIN\n")[1]
        rows = [r for r in body.splitlines() if r and not r.startswith("//")]
        self.assertEqual(len(rows), -(-self.d.total // 60))
        joined = "".join("".join(r.split()[1:]) for r in rows)
        self.assertEqual(joined, self.d.seq.lower())

    def test_record_is_terminated(self):
        self.assertTrue(self.text.rstrip().endswith("//"))

    def test_non_ascii_names_are_transliterated(self):
        self.assertNotIn("é", to_genbank(design(name="café vector")))


class Fasta(unittest.TestCase):
    def test_header_carries_the_name_and_length(self):
        d = design()
        header = to_fasta(d).splitlines()[0]
        self.assertTrue(header.startswith(">"))
        self.assertIn(str(d.total), header)

    def test_header_has_no_whitespace_in_the_identifier(self):
        self.assertNotIn(" ", to_fasta(design(name="my vector")).splitlines()[0].split()[0])

    def test_sequence_wraps_at_seventy_columns_and_round_trips(self):
        d = design()
        rows = to_fasta(d).splitlines()[1:]
        self.assertTrue(all(len(r) <= 70 for r in rows))
        self.assertEqual("".join(rows), d.seq)


class Svg(unittest.TestCase):
    def test_output_is_well_formed_xml(self):
        root = ET.fromstring(to_svg(design()))
        self.assertTrue(root.tag.endswith("svg"))

    def test_one_arc_and_one_legend_row_per_feature(self):
        d = design()
        root = ET.fromstring(to_svg(d))
        paths = [e for e in root.iter() if e.tag.endswith("path")]
        self.assertEqual(len(paths), len(d.features))

    def test_feature_names_with_xml_characters_are_escaped(self):
        root = ET.fromstring(to_svg(design(gene="A&B <test>")))
        titles = [e.text for e in root.iter() if e.tag.endswith("title")]
        self.assertTrue(any("A&B <test>" in (t or "") for t in titles))

    def test_canvas_grows_with_the_legend(self):
        small = ET.fromstring(to_svg(design(vector="plasmid", marker="none")))
        large = ET.fromstring(to_svg(design(strategy="car", gene="CD19", marker="EGFP")))
        self.assertLess(int(small.get("height")), int(large.get("height")))

    def test_arc_colours_come_from_the_kind_table(self):
        svg = to_svg(design())
        for feature in design().features:
            self.assertIn(KINDS[feature.kind][0], svg)


class Summary(unittest.TestCase):
    def test_table_lists_every_element_id_position_and_length(self):
        d = design()
        text = to_summary(d)
        for feature in d.features:
            pattern = r"^%s\s+%d\.\.%d\s+%d\s" % (re.escape(feature.id), feature.start,
                                                    feature.end, feature.length)
            self.assertRegex(text, re.compile(pattern, re.M))

    def test_design_checks_are_included_with_their_levels(self):
        text = to_summary(design())
        for level, message in design().messages:
            self.assertIn("[%s] %s" % (level, message), text)

    def test_user_supplied_elements_are_flagged(self):
        self.assertIn("(your sequence)", to_summary(design()))

    def test_placeholder_only_design_is_not_flagged_as_user_supplied(self):
        self.assertNotIn("(your sequence)", to_summary(design(sequence="")))

    def test_design_basis_is_shown_when_present(self):
        self.assertIn("picked from the CNV run", to_summary(design(basis="picked from the CNV run")))

    def test_summary_mentions_how_to_replace_placeholders(self):
        self.assertIn("--override", to_summary(design()))


class WriteOutputs(unittest.TestCase):
    def test_writes_four_files_named_from_a_sanitised_stem(self):
        d = design(name="pLV/My Vector*1")
        with tempfile.TemporaryDirectory() as folder:
            paths = write_outputs(d, folder)
            self.assertEqual([os.path.splitext(p)[1] for p in paths],
                             [".gb", ".fasta", ".svg", ".txt"])
            for path in paths:
                self.assertTrue(os.path.isfile(path))
                stem = os.path.basename(path)
                self.assertNotRegex(stem, r"[/*\s]")

    def test_creates_the_folder_if_missing(self):
        with tempfile.TemporaryDirectory() as folder:
            nested = os.path.join(folder, "a", "b")
            write_outputs(design(), nested)
            self.assertTrue(os.path.isdir(nested))

    def test_files_use_unix_line_endings(self):
        with tempfile.TemporaryDirectory() as folder:
            path = write_outputs(design(), folder)[0]
            with open(path, "rb") as fh:
                self.assertNotIn(b"\r\n", fh.read())


class AneuploidyFiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        segments, vtype, notes, sample = parse_segments(example_cnv_text())
        cls.res = analyze_aneuploidy(segments, value_type=vtype, sample=sample, notes=notes)
        cls.recs = recommend(cls.res)

    def test_report_states_the_score_out_of_the_scorable_arms(self):
        text = aneuploidy_report(self.res, self.recs)
        self.assertIn("Aneuploidy score: 14 of %d" % N_SCORED_ARMS, text)

    def test_report_lists_called_arms_and_the_whole_chromosome_event(self):
        text = aneuploidy_report(self.res, self.recs)
        self.assertIn("Gained arms: 1q, 3q, 5p, 7p, 7q, 8q, 20q", text)
        self.assertIn("Whole-chromosome events: 7 gain", text)

    def test_report_numbers_the_suggestions_for_design_from_cnv(self):
        text = aneuploidy_report(self.res, self.recs)
        for i, rec in enumerate(self.recs, 1):
            self.assertIn("%d. %s: %s" % (i, rec.title, rec.gene), text)

    def test_report_always_carries_the_limits_section(self):
        self.assertIn("Limits", aneuploidy_report(self.res, self.recs))

    def test_arms_tsv_has_a_header_and_one_row_per_arm(self):
        rows = arms_tsv(self.res).strip().split("\n")
        self.assertEqual(rows[0].split("\t")[0], "arm")
        self.assertEqual(len(rows), len(self.res.arms) + 1)
        self.assertTrue(all(len(r.split("\t")) == 6 for r in rows))

    def test_writes_report_and_table_with_a_sanitised_stem(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = write_aneuploidy_files(self.res, self.recs, folder, "TCGA/XX 01")
            self.assertEqual([os.path.basename(p) for p in paths],
                             ["TCGA_XX_01_aneuploidy.txt", "TCGA_XX_01_arms.tsv"])

    def test_blank_stem_falls_back_to_sample(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = write_aneuploidy_files(self.res, self.recs, folder, "***")
            self.assertTrue(os.path.basename(paths[0]).startswith("sample_"))


class ExampleFile(unittest.TestCase):
    def test_example_is_parseable_and_labelled_synthetic(self):
        text = example_cnv_text()
        self.assertIn("not patient data", text.splitlines()[0])
        segments, vtype, _notes, _sample = parse_segments(text)
        self.assertEqual(vtype, "log2")
        self.assertTrue(segments)

    def test_example_is_available_for_both_builds(self):
        self.assertNotEqual(example_cnv_text("hg38"), example_cnv_text("hg19"))


if __name__ == "__main__":
    unittest.main()
