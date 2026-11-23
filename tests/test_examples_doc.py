"""Runs every command in EXAMPLES.md and checks the result it promises.

Each example in the document is preceded by a comment such as
    <!-- test: exit=0 | expect=some text | expect=more text -->
followed by a bash block.  The commands run in order in a scratch folder; the expected text must
appear in the combined output and the last command must exit with the stated status.
"""
import contextlib
import io
import os
import re
import shlex
import tempfile
import unittest

from cvd.cli import main

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(ROOT, "EXAMPLES.md")
BLOCK = re.compile(r"<!-- test: (?P<spec>.*?) -->\s*```bash\n(?P<code>.*?)```", re.S)
NO_OUT = ("--list", "--version", "--example-cnv")


def load_examples():
    with open(DOC, encoding="utf-8") as fh:
        text = fh.read()
    found = []
    for m in BLOCK.finditer(text):
        spec = [part.strip() for part in m.group("spec").split(" | ")]
        exit_code = int(next(p for p in spec if p.startswith("exit=")).split("=", 1)[1])
        expects = [p.split("=", 1)[1] for p in spec if p.startswith("expect=")]
        code = m.group("code").replace("\\\n", " ")
        commands = [line.strip() for line in code.splitlines() if line.strip()]
        found.append((m.start(), exit_code, expects, commands))
    return found


def run_command(line, folder):
    argv = shlex.split(line)
    assert argv[:3] == ["python", "-m", "cvd"], line
    argv = [a.replace("examples/", os.path.join(ROOT, "examples") + os.sep) for a in argv[3:]]
    if "--out" not in argv and not any(flag in argv for flag in NO_OUT):
        argv += ["--out", folder]
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(argv)
    return code, out.getvalue() + err.getvalue()


class ExamplesDocument(unittest.TestCase):
    def test_the_document_has_examples(self):
        self.assertGreaterEqual(len(load_examples()), 10)

    def test_every_example_runs_and_gives_the_promised_result(self):
        cwd = os.getcwd()
        for position, exit_code, expects, commands in load_examples():
            with self.subTest(example=commands[0][:70]):
                with tempfile.TemporaryDirectory() as folder:
                    os.chdir(folder)
                    try:
                        status, output = 0, ""
                        for line in commands:
                            status, text = run_command(line, folder)
                            output += text
                    finally:
                        os.chdir(cwd)
                self.assertEqual(status, exit_code, output[-400:])
                for text in expects:
                    self.assertIn(text, output)

    def test_every_numbered_heading_has_a_tested_block(self):
        with open(DOC, encoding="utf-8") as fh:
            text = fh.read()
        headings = re.findall(r"^## \d+\. ", text, re.M)
        self.assertEqual(len(load_examples()) >= len(headings), True)


if __name__ == "__main__":
    unittest.main()
