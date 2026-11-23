"""The worked example in the README: a real TCGA tumor from open data."""
import contextlib
import io
import os
import tempfile
import unittest

from cvd.cli import main

EXAMPLE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "examples", "tcga_chol_W5-AA2X-01.seg")


def run(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(list(argv))
    return code, out.getvalue(), err.getvalue()


@unittest.skipUnless(os.path.exists(EXAMPLE), "example file not present")
class ReadmeExample(unittest.TestCase):
    def test_the_aneuploidy_numbers_in_the_readme_hold(self):
        with tempfile.TemporaryDirectory() as folder:
            code, out, _err = run("--cnv", EXAMPLE, "--build", "hg19", "--out", folder)
        self.assertEqual(code, 0)
        self.assertIn("Aneuploidy score: 27 of 39 scorable arms (16 gains, 11 losses)", out)
        self.assertIn("Fraction of the genome gained or lost: 76%", out)
        self.assertIn("1. NKG2D ligand flag: ULBP2", out)
        self.assertIn("3q gain: CD47", out)

    def test_the_ranking_in_the_readme_holds(self):
        with tempfile.TemporaryDirectory() as folder:
            code, out, _err = run("--cnv", EXAMPLE, "--build", "hg19", "--design-from-cnv", "1",
                                  "--setting", "in vivo", "--rank", "--out", folder)
        self.assertEqual(code, 0)
        rows = [line.split() for line in out.splitlines() if line[:1] in "123456" and "  " in line]
        scores = [int(r[1]) for r in rows[:6]]
        self.assertEqual(scores, [97, 92, 92, 80, 75, 75])
        self.assertIn("tumor-selective promoter (Survivin)", out.split("Rank Score")[1].splitlines()[1])

    def test_building_option_one_records_the_basis(self):
        with tempfile.TemporaryDirectory() as folder:
            code, _out, _err = run("--cnv", EXAMPLE, "--build", "hg19", "--design-from-cnv", "1",
                                   "--setting", "in vivo", "--rank", "--use-option", "1", "--out", folder)
            self.assertEqual(code, 0)
            with open(os.path.join(folder, "pAAV-Survivin_ULBP2.txt"), encoding="utf-8") as fh:
                summary = fh.read()
        self.assertIn("Design basis", summary)
        self.assertIn("TCGA-W5-AA2X-01", summary)
        self.assertIn("Vector option #1 of 6", summary)


if __name__ == "__main__":
    unittest.main()
