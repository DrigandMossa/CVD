"""The Vector options tab, driven through stand-in widgets.

These tests replace tkinter with small fakes, so they run on machines without a display or
without tkinter at all.  They check the logic that connects the ranking to the window (what is
listed, what a click shows, what "use this option" changes), not how the widgets look.
"""
import importlib
import sys
import tempfile
import os
import unittest
from unittest import mock


class Var:
    def __init__(self, value=None):
        self.value, self.callbacks = value, []

    def get(self):
        return self.value

    def set(self, value):
        self.value = value
        for callback in self.callbacks:
            callback("n", "", "write")

    def trace_add(self, mode, callback):
        self.callbacks.append(callback)


class FakeText:
    def __init__(self, *args, **kwargs):
        self.content, self.parts, self.state, self._m = "", [], None, mock.MagicMock()

    def __getattr__(self, name):
        return getattr(self._m, name)

    def get(self, *args):
        return self.content

    def configure(self, **kwargs):
        if "state" in kwargs:
            self.state = kwargs["state"]

    def insert(self, index, text, *tags):
        if tags:
            assert self.state == "normal", "inserting into a read-only Text"
            self.parts.append((tags[0][0] if isinstance(tags[0], tuple) else tags[0], text))
            self.content += text
        else:
            self.content = self.content + text if index == "end" else text

    def delete(self, *args):
        self.content, self.parts = "", []


class FakeTree:
    def __init__(self, *args, **kwargs):
        self.rows, self.order, self.sel = {}, [], ()

    def __getattr__(self, name):
        return mock.MagicMock()

    def get_children(self):
        return tuple(self.order)

    def delete(self, *ids):
        self.rows, self.order, self.sel = {}, [], ()

    def insert(self, parent, where, iid=None, text="", values=(), tags=(), **kwargs):
        self.rows[iid] = (text, values, tags)
        self.order.append(iid)

    def exists(self, iid):
        return iid in self.rows

    def selection(self):
        return self.sel

    def selection_set(self, iid):
        self.sel = (iid,)

    def focus(self, iid):
        pass


class OptionsTab(unittest.TestCase):
    def setUp(self):
        fake_tk = mock.MagicMock()
        fake_tk.StringVar = fake_tk.BooleanVar = Var
        fake_tk.TclError = type("TclError", (Exception,), {})
        fake_tk.Text = FakeText
        fake_ttk = mock.MagicMock()
        fake_ttk.Treeview.side_effect = FakeTree
        fake_ttk.Label.side_effect = lambda *a, **k: mock.MagicMock()
        fake_tk.ttk, fake_tk.font = fake_ttk, mock.MagicMock()
        modules = {"tkinter": fake_tk, "tkinter.font": fake_tk.font, "tkinter.ttk": fake_ttk,
                   "tkinter.filedialog": fake_tk.filedialog, "tkinter.messagebox": fake_tk.messagebox}
        for name in ("cvd.gui.app", "cvd.gui.aneuploidy"):
            sys.modules.pop(name, None)
        patcher = mock.patch.dict(sys.modules, modules)
        patcher.start()
        self.addCleanup(self.restore, patcher)
        self.fake_tk = fake_tk
        self.module = importlib.import_module("cvd.gui.app")
        self.module.messagebox = mock.MagicMock()
        self.module.filedialog = mock.MagicMock()
        self.root = mock.MagicMock()
        self.root.after.side_effect = lambda ms, fn: 1
        self.app = self.module.VectorDesignerApp(self.root)
        self.app.refresh()
        self.tree = self.app.opt_tree

    @staticmethod
    def restore(patcher):
        patcher.stop()
        for name in ("cvd.gui.app", "cvd.gui.aneuploidy"):
            sys.modules.pop(name, None)

    def select(self, key):
        self.tree.sel = (key,)
        self.app.on_option_select()

    def test_every_option_is_listed_in_rank_order_with_the_first_selected(self):
        res = self.app.ranking
        self.assertEqual(self.tree.order, [o.key for o in res.options])
        self.assertEqual([self.tree.rows[k][0][:2] for k in self.tree.order[:2]], ["1.", "2."])
        self.assertEqual(self.tree.sel, (res.options[0].key,))

    def test_clicking_an_option_shows_its_own_reasoning_and_modifications(self):
        seen = set()
        for key in self.tree.order:
            self.select(key)
            option = self.app.current_option()
            self.assertEqual(option.key, key)
            self.assertIn(option.title, self.app.why_box.content)
            self.assertIn("Modifications to option %d" % option.rank, self.app.mod_box.content)
            seen.add(self.app.why_box.content)
        self.assertEqual(len(seen), len(self.tree.order))
        self.assertEqual(self.app.why_box.state, "disabled")
        self.assertEqual(self.app.mod_box.state, "disabled")

    def test_the_modifications_tab_title_carries_the_count(self):
        self.select(self.tree.order[0])
        text = self.app.detail.tab.call_args.kwargs["text"]
        count = len(self.app.current_option().applied_mods)
        self.assertEqual(text, "Modifications (%d)" % count)

    def test_the_setting_reranks_and_keeps_a_valid_selection(self):
        self.app.var["setting"].set("In vivo")
        self.app.refresh()
        self.assertEqual(self.app.ranking.options[0].vector, "aav")
        self.app.var["setting"].set("Ex vivo")
        self.app.refresh()
        self.assertEqual(self.app.ranking.options[0].vector, "lenti")
        self.assertTrue(self.tree.exists(self.tree.sel[0]))

    def test_the_selection_survives_a_refresh(self):
        self.tree.sel = (self.tree.order[2],)
        self.app.refresh()
        self.assertEqual(self.tree.sel, (self.tree.order[2],))

    def test_one_small_map_per_option_and_a_click_selects_and_expands_it(self):
        g = self.app.gallery
        clicks = [c for c in g.tag_bind.call_args_list if c.args[1] == "<Button-1>"]
        clicks = clicks[-len(self.app.ranking.options):]
        self.assertEqual(len(clicks), len(self.app.ranking.options))
        self.fake_tk.Toplevel.reset_mock()
        for i, click in enumerate(clicks):
            click.args[2](None)
            self.assertEqual(self.tree.sel, (self.app.ranking.options[i].key,))
        self.assertEqual(self.fake_tk.Toplevel.call_count, len(clicks))

    def test_every_option_can_be_expanded(self):
        for key in list(self.tree.order):
            self.app.expand_option(key)

    def prepare_in_vivo_with_marker_and_tag(self):
        app = self.app
        app.var["setting"].set("In vivo")
        app.var["marker"].set("EGFP")
        app.var["link"].set("IRES")
        app.var["tag"].set("HA")
        app.insert_box.content = "ATG" + "GCC" * 700
        app.refresh()

    def test_using_a_modified_option_carries_every_modification_into_the_designer(self):
        self.prepare_in_vivo_with_marker_and_tag()
        option = next(o for o in self.app.ranking.options if o.key == "aav|strong")
        applied = {m.key for m in option.applied_mods}
        self.assertTrue({"tag", "marker", "wpre-short"} <= applied)
        self.select(option.key)
        self.app.use_option()
        v = self.app.var
        self.assertEqual((v["vector"].get(), v["marker"].get(), v["tag"].get(), v["resistance"].get()),
                         ("aav", "none", "none", "KanR"))
        self.assertTrue(self.app.wpre_short.get())
        design = self.module.build_design(self.app.make_config())
        self.assertEqual(sum(f.length for f in design.features if f.in_cassette), option.cassette)
        self.assertIn("Modifications:", self.app.make_config().basis)
        self.app.tabs.select.assert_called_with(self.app.design_tab)

    def test_a_manual_change_afterwards_forgets_the_option_basis(self):
        self.select(self.tree.order[0])
        self.app.use_option()
        self.assertTrue(self.app.option_basis)
        self.app.on_vector()
        self.assertEqual(self.app.option_basis, "")

    def test_an_option_that_does_not_fit_asks_first(self):
        self.app.insert_box.content = "ATG" + "GCC" * 3700
        self.app.refresh()
        bad = next(o for o in self.app.ranking.options if not o.feasible)
        self.select(bad.key)
        self.module.messagebox.askyesno.return_value = False
        self.app.use_option()
        self.assertEqual(self.app.option_basis, "")
        self.module.messagebox.askyesno.return_value = True
        self.app.use_option()
        self.assertTrue(self.app.option_basis)

    def test_switching_modifications_off_shows_the_plain_layouts(self):
        self.app.optimize.set(False)
        self.app.refresh()
        self.assertTrue(all(not o.mods for o in self.app.ranking.options))
        self.assertIn("switched off", self.app.mod_box.content)

    def test_copy_and_save(self):
        self.select(self.tree.order[0])
        self.app.copy_option()
        self.root.clipboard_append.assert_called()
        self.assertIn("Modifications to option 1", self.root.clipboard_append.call_args.args[0])
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "o.txt")
            self.module.filedialog.asksaveasfilename.return_value = path
            self.app.save_options()
            with open(path, encoding="utf-8") as fh:
                self.assertIn("SCORE BREAKDOWN", fh.read())

    def test_nothing_selected_asks_for_a_selection(self):
        self.tree.sel = ()
        self.app.use_option()
        self.module.messagebox.showinfo.assert_called()

    def test_every_strategy_and_setting_ranks_and_renders(self):
        from cvd.parts import STRATEGIES
        for key in STRATEGIES:
            for setting in ("Not specified", "In vivo", "Ex vivo"):
                self.app.var["strategy"].set(STRATEGIES[key]["name"])
                self.app.on_strategy_chosen()
                self.app.var["setting"].set(setting)
                self.app.insert_box.content = ""
                self.app.refresh()
                self.assertTrue(self.tree.sel)
                for option_key in self.tree.order:
                    self.select(option_key)


if __name__ == "__main__":
    unittest.main()
