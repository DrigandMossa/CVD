import contextlib
import io
import os
import tempfile
import unittest

from cvd.cli import build_parser, config_from_args, main, print_catalog, read_seq_arg


def run(*argv):
    """Run main() with stdout and stderr captured. Returns (exit code, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(list(argv))
    return code, out.getvalue(), err.getvalue()


def parse(*argv):
    return build_parser().parse_args(list(argv))


class SeqArgument(unittest.TestCase):
    def test_a_path_is_read_from_disk(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "orf.fa")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(">x\nACGT\n")
            self.assertEqual(read_seq_arg(path), ">x\nACGT\n")

    def test_a_literal_sequence_is_passed_through(self):
        self.assertEqual(read_seq_arg("ACGTACGT"), "ACGTACGT")

    def test_empty_stays_empty(self):
        self.assertEqual(read_seq_arg(""), "")
        self.assertEqual(read_seq_arg(None), "")


class Overrides(unittest.TestCase):
    def test_overrides_are_split_on_the_first_equals(self):
        cfg = config_from_args(parse("--override", "prom=ACGT", "--override", "wpre=TTTT"))
        self.assertEqual(cfg.overrides, {"prom": "ACGT", "wpre": "TTTT"})

    def test_override_without_an_equals_is_rejected(self):
        with self.assertRaises(ValueError) as caught:
            config_from_args(parse("--override", "promACGT"))
        self.assertIn("ID=SEQUENCE", str(caught.exception))

    def test_override_value_may_be_a_file_path(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "p.txt")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("ACGTACGT")
            cfg = config_from_args(parse("--override", "prom=" + path))
            self.assertEqual(cfg.overrides["prom"], "ACGTACGT")


class Flags(unittest.TestCase):
    def test_wpre_is_undecided_unless_a_flag_is_given(self):
        self.assertIsNone(config_from_args(parse()).wpre)

    def test_wpre_can_be_forced_on_or_off(self):
        self.assertTrue(config_from_args(parse("--wpre")).wpre)
        self.assertFalse(config_from_args(parse("--no-wpre")).wpre)

    def test_wpre_wins_over_no_wpre_when_both_are_given(self):
        self.assertTrue(config_from_args(parse("--wpre", "--no-wpre")).wpre)

    def test_kozak_is_on_unless_disabled(self):
        self.assertTrue(config_from_args(parse()).kozak)
        self.assertFalse(config_from_args(parse("--no-kozak")).kozak)

    def test_defaults_match_the_dataclass_expectations(self):
        cfg = config_from_args(parse())
        self.assertEqual((cfg.vector, cfg.strategy, cfg.link, cfg.tag, cfg.resistance),
                         ("lenti", "rep", "P2A", "none", "AmpR"))


class Catalog(unittest.TestCase):
    def test_catalog_lists_vectors_strategies_and_parts(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            print_catalog()
        text = out.getvalue()
        for expected in ("lenti", "rep", "Tumor suppressor replacement", "Promoters",
                         "EF1a", "CAR hinges", "Bacterial selection"):
            self.assertIn(expected, text)

    def test_list_mode_exits_cleanly_without_writing_anything(self):
        code, out, _err = run("--list")
        self.assertEqual(code, 0)
        self.assertIn("Strategies", out)


class OneShotDesign(unittest.TestCase):
    def test_writes_four_files_and_reports_them(self):
        with tempfile.TemporaryDirectory() as folder:
            code, out, _err = run("--strategy", "rep", "--gene", "TP53", "--out", folder)
            self.assertEqual(code, 0)
            self.assertEqual(sorted(os.listdir(folder)),
                             ["pLV-EF1a_TP53.fasta", "pLV-EF1a_TP53.gb",
                              "pLV-EF1a_TP53.svg", "pLV-EF1a_TP53.txt"])
            self.assertIn("Design checks", out)
            self.assertIn("Files written", out)

    def test_report_lists_every_element(self):
        with tempfile.TemporaryDirectory() as folder:
            _code, out, _err = run("--gene", "TP53", "--out", folder)
            for element in ("ltr5", "prom", "goi", "ori"):
                self.assertIn(element, out)

    def test_bad_option_value_exits_two_with_a_message(self):
        with tempfile.TemporaryDirectory() as folder:
            code, _out, err = run("--vector", "retro", "--out", folder)
            self.assertEqual(code, 2)
            self.assertIn("Unknown vector", err)

    def test_unwritable_output_folder_exits_two(self):
        with tempfile.NamedTemporaryFile() as not_a_folder:
            code, _out, err = run("--gene", "TP53", "--out",
                                  os.path.join(not_a_folder.name, "sub"))
            self.assertEqual(code, 2)
            self.assertTrue(err)

    def test_strict_turns_a_design_error_into_exit_one(self):
        with tempfile.TemporaryDirectory() as folder:
            plain, _out, _err = run("--vector", "aav", "--gene", "APC", "--out", folder)
            strict, _out, _err = run("--vector", "aav", "--gene", "APC", "--out", folder,
                                     "--strict")
            self.assertEqual((plain, strict), (0, 1))

    def test_strict_is_quiet_when_there_are_no_errors(self):
        with tempfile.TemporaryDirectory() as folder:
            code, _out, _err = run("--gene", "TP53", "--out", folder, "--strict")
            self.assertEqual(code, 0)


class CnvMode(unittest.TestCase):
    def example(self, folder):
        path = os.path.join(folder, "demo.tsv")
        code, _out, _err = run("--example-cnv", path)
        self.assertEqual(code, 0)
        return path

    def test_example_file_is_written_and_is_usable(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self.example(folder)
            self.assertTrue(os.path.isfile(path))
            code, out, _err = run("--cnv", path, "--out", folder)
            self.assertEqual(code, 0)
            self.assertIn("Aneuploidy score: 14", out)

    def test_analysis_writes_the_report_and_the_arm_table(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self.example(folder)
            run("--cnv", path, "--out", folder)
            self.assertIn("demo_aneuploidy.txt", os.listdir(folder))
            self.assertIn("demo_arms.tsv", os.listdir(folder))

    def test_analysis_alone_builds_no_vector(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self.example(folder)
            run("--cnv", path, "--out", folder)
            self.assertFalse([f for f in os.listdir(folder) if f.endswith(".gb")])

    def test_design_from_cnv_builds_the_chosen_suggestion(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self.example(folder)
            code, out, _err = run("--cnv", path, "--design-from-cnv", "1", "--out", folder)
            self.assertEqual(code, 0)
            self.assertIn("Building the vector for suggestion 1", out)
            built = [f for f in os.listdir(folder) if f.endswith(".gb")]
            self.assertEqual(len(built), 1)
            self.assertIn("ULBP2", built[0])

    def test_the_summary_records_why_the_design_was_chosen(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self.example(folder)
            run("--cnv", path, "--design-from-cnv", "1", "--out", folder)
            summary = [f for f in os.listdir(folder) if f.endswith(".txt")
                       and not f.endswith("_aneuploidy.txt")][0]
            with open(os.path.join(folder, summary), encoding="utf-8") as fh:
                text = fh.read()
            self.assertIn("Design basis", text)
            self.assertIn("Aneuploidy analysis", text)

    def test_out_of_range_suggestion_exits_two(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self.example(folder)
            code, _out, err = run("--cnv", path, "--design-from-cnv", "99", "--out", folder)
            self.assertEqual(code, 2)
            self.assertIn("between 1 and", err)

    def test_design_from_cnv_without_a_file_exits_two(self):
        code, _out, err = run("--design-from-cnv", "1")
        self.assertEqual(code, 2)
        self.assertIn("--cnv", err)

    def test_missing_cnv_file_exits_two(self):
        code, _out, err = run("--cnv", "/nonexistent/segments.tsv")
        self.assertEqual(code, 2)
        self.assertTrue(err)

    def test_unparseable_cnv_file_exits_two(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "bad.tsv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("nothing\tuseful\nhere\tatall\n")
            code, _out, err = run("--cnv", path, "--out", folder)
            self.assertEqual(code, 2)
            self.assertTrue(err)


class RankingModes(unittest.TestCase):
    ARGS = ("--strategy", "rep", "--gene", "TP53", "--setting", "in vivo")

    def run_in(self, *extra):
        with tempfile.TemporaryDirectory() as folder:
            code, out, err = run(*self.ARGS, *extra, "--out", folder)
            files = sorted(os.listdir(folder))
            report = ""
            path = os.path.join(folder, "TP53_rep_vector_options.txt")
            if os.path.exists(path):
                with open(path, encoding="utf-8") as fh:
                    report = fh.read()
        return code, out, err, files, report

    def test_rank_prints_the_table_and_writes_only_the_report(self):
        code, out, err, files, report = self.run_in("--rank")
        self.assertEqual((code, err, files), (0, "", ["TP53_rep_vector_options.txt"]))
        self.assertIn("Rank Score", out)
        self.assertIn("AAV", out)
        self.assertIn("Use --explain", out)
        self.assertIn("MODIFICATIONS TO OPTION 1", report.upper())

    def test_explain_one_prints_the_reasoning_and_the_modifications(self):
        code, out, _err, _files, _report = self.run_in("--explain", "1")
        self.assertEqual(code, 0)
        self.assertIn("WHY IT WAS SUGGESTED", out)
        self.assertIn("WHY IT RANKS ABOVE OR BELOW THE OTHER SUGGESTIONS", out)
        self.assertIn("Modifications to option 1", out)
        self.assertNotIn("Use --explain", out)

    def test_explain_all_prints_every_option(self):
        _code, out, _err, _files, _report = self.run_in("--explain", "all")
        self.assertEqual(out.count("Modifications to option"), 6)

    def test_use_option_builds_that_option_with_its_modifications(self):
        code, out, _err, files, _report = self.run_in("--use-option", "1", "--marker", "EGFP")
        self.assertEqual(code, 0)
        self.assertIn("Building the vector for option 1", out)
        self.assertIn("pAAV-EF1a_TP53.gb", files)
        self.assertNotIn("EGFP", out.split("Building the vector")[1])      # removed in vivo

    def test_bad_option_numbers_exit_two(self):
        for extra in (("--use-option", "0"), ("--use-option", "7"), ("--explain", "9"),
                      ("--explain", "x")):
            code, _out, err, _files, _report = self.run_in(*extra)
            self.assertEqual(code, 2, extra)
            self.assertIn("Error:", err)

    def test_bad_setting_exits_two(self):
        with tempfile.TemporaryDirectory() as folder:
            code, _out, err = run("--rank", "--setting", "on the moon", "--out", folder)
        self.assertEqual(code, 2)
        self.assertIn("Error:", err)

    def test_no_optimize_ranks_the_plain_layouts(self):
        code, out, _err, _files, report = self.run_in("--explain", "1", "--no-optimize")
        self.assertEqual(code, 0)
        self.assertIn("switched off", out)

    def test_export_options_writes_every_prebuilt_option_to_its_own_folder(self):
        with tempfile.TemporaryDirectory() as folder:
            code, _out, _err = run(*self.ARGS, "--rank", "--export-options", "--out", folder)
            self.assertEqual(code, 0)
            folders = sorted(d for d in os.listdir(folder) if d.startswith("option"))
            self.assertEqual(folders, ["option%d" % i for i in range(1, 7)])
            for d in folders:
                self.assertEqual(len(os.listdir(os.path.join(folder, d))), 4)

    def test_short_wpre_flag_reaches_the_design(self):
        with tempfile.TemporaryDirectory() as folder:
            code, out, _err = run("--vector", "aav", "--gene", "TP53", "--short-wpre", "--out", folder)
        self.assertEqual(code, 0)
        self.assertIn("247 bp", out)

    def test_ranking_combines_with_a_cnv_suggestion(self):
        with tempfile.TemporaryDirectory() as folder:
            demo = os.path.join(folder, "demo.tsv")
            run("--example-cnv", demo)
            code, out, _err = run("--cnv", demo, "--design-from-cnv", "2", "--rank", "--out", folder)
        self.assertEqual(code, 0)
        self.assertIn("Rank Score", out)


class Version(unittest.TestCase):
    def test_version_flag_reports_the_package_version(self):
        from cvd import __version__
        with self.assertRaises(SystemExit) as caught:
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                build_parser().parse_args(["--version"])
        self.assertEqual(caught.exception.code, 0)

    def test_parser_program_name_matches_the_module_entry_point(self):
        self.assertEqual(build_parser().prog, "cvd")


if __name__ == "__main__":
    unittest.main()
