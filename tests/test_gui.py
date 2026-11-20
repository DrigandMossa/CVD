"""GUI-layer tests.

Everything that needs a real tkinter is skipped when tkinter is missing, so the suite
still passes on a bare Linux container; CI runners have tkinter and do run these.
"""
import contextlib
import importlib.util
import io
import unittest

import cvd.gui

HAS_TK = importlib.util.find_spec("tkinter") is not None


class Fallback(unittest.TestCase):
    def test_launch_gui_reports_missing_tkinter_instead_of_raising(self):
        original = cvd.gui.tk
        cvd.gui.tk = None
        try:
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                code = cvd.gui.launch_gui()
        finally:
            cvd.gui.tk = original
        self.assertEqual(code, 2)
        self.assertIn("tkinter", err.getvalue())
        self.assertIn("python3-tk", err.getvalue())

    def test_importing_the_gui_package_never_needs_tkinter(self):
        self.assertTrue(hasattr(cvd.gui, "launch_gui"))


@unittest.skipUnless(HAS_TK, "tkinter is not installed")
class ElementNotes(unittest.TestCase):
    """The 'Sequence' column tells the user where each element's bases come from."""

    def setUp(self):
        from cvd.design import build_design
        from cvd.gui.app import element_note
        from cvd.model import Config
        self.note = element_note
        self.features = {f.id: f for f in build_design(
            Config(gene="TP53", marker="EGFP")).features}

    def test_insert_reflects_whether_the_user_supplied_a_sequence(self):
        goi = self.features["goi"]
        self.assertEqual(self.note(goi, True, False), "your sequence")
        self.assertEqual(self.note(goi, False, False), "placeholder (N)")

    def test_built_in_elements_say_so(self):
        self.assertEqual(self.note(self.features["p2a"], False, False), "built in")
        self.assertEqual(self.note(self.features["kozak"], False, False), "built in")

    def test_placeholder_elements_invite_a_paste(self):
        self.assertIn("double-click", self.note(self.features["prom"], False, False))

    def test_an_overridden_element_is_marked_as_the_users_own(self):
        from cvd.design import build_design
        from cvd.model import Config
        design = build_design(Config(gene="TP53", overrides={"prom": "ACGT"}))
        prom = next(f for f in design.features if f.id == "prom")
        self.assertEqual(self.note(prom, False, False), "your sequence")

    def test_guide_element_follows_the_target_box(self):
        from cvd.design import build_design
        from cvd.model import Config
        design = build_design(Config(strategy="guide", gene="KRAS", target="ACGTACGTACGTACGTACGT"))
        guide = next(f for f in design.features if f.id == "guide")
        self.assertEqual(self.note(guide, False, True), "your sequence")
        self.assertEqual(self.note(guide, False, False), "placeholder (N)")


@unittest.skipUnless(HAS_TK, "tkinter is not installed")
class WindowClasses(unittest.TestCase):
    def test_the_window_modules_import_and_expose_their_classes(self):
        from cvd.gui.aneuploidy import AneuploidyWindow
        from cvd.gui.app import VectorDesignerApp
        self.assertTrue(callable(VectorDesignerApp))
        self.assertTrue(callable(AneuploidyWindow))

    def test_message_levels_cover_every_level_the_designer_emits(self):
        from cvd.design import build_design
        from cvd.gui.app import VectorDesignerApp
        from cvd.model import Config
        design = build_design(Config(vector="aav", gene="APC", sequence="ATGTAACGT"))
        for level, _text in design.messages:
            self.assertIn(level, VectorDesignerApp.LEVELS)

    def test_the_biosafety_notice_is_present_in_the_window_text(self):
        from cvd.gui.app import VectorDesignerApp
        self.assertIn("biosafety", VectorDesignerApp.BIOSAFETY)


if __name__ == "__main__":
    unittest.main()
