"""Main window: options on the left, live map, element table and checks on the right."""

import dataclasses
import sys
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, ttk
from typing import Dict, Optional

from ..design.builder import build_design
from ..export.fasta import to_fasta
from ..export.files import write_outputs
from ..export.genbank import to_genbank
from ..model import Config, Design, Feature, Recommendation
from ..parts.elements import (BUILT_IN_IDS, COSTIM, HINGES, KINDS, LINKS, MARKERS, POLYA,
                              PROMOTERS, RESISTANCE, TAGS)
from ..parts.strategies import STRATEGIES
from ..seqops import clean_sequence
from ..version import __version__
from .aneuploidy import AneuploidyWindow


def element_note(f: Feature, has_insert: bool, has_target: bool) -> str:
    """Short description of where an element's sequence comes from."""
    if f.custom:
        return "your sequence"
    if f.id == "goi":
        return "your sequence" if has_insert else "placeholder (N)"
    if f.id == "guide":
        return "your sequence" if has_target else "placeholder (N)"
    if f.id in BUILT_IN_IDS:
        return "built in"
    return "placeholder (N), double-click to paste"


class VectorDesignerApp:
    """Main window: options on the left; live map, element table and checks on the right."""

    MAP_SIZE = 360
    LEVELS = {"error": ("Error", "#b3261e"),
              "warn": ("Warning", "#a55a00"),
              "info": ("Note", "#587083")}
    BIOSAFETY = ("Research design tool. Producing virus requires institutional biosafety approval, and clinical "
                 "gene therapy products need GMP manufacturing, preclinical safety data and regulatory authorization.")

    def __init__(self, root) -> None:
        self.root = root
        root.title("Cancer Gene Therapy Vector Designer")
        root.geometry("1220x780")
        root.minsize(1000, 640)
        if sys.platform.startswith("linux"):
            try:
                ttk.Style().theme_use("clam")
            except tk.TclError:
                pass
        self.bold = tkfont.nametofont("TkDefaultFont").copy()
        self.bold.configure(weight="bold")
        self.mono = tkfont.nametofont("TkFixedFont")
        self.overrides: Dict[str, str] = {}
        self.design: Optional[Design] = None
        self.pending = None
        self.arcs: Dict[str, int] = {}
        self.swatches: Dict[str, object] = {}
        self.summary_line = ""
        self.has_insert = False
        self.has_target = False
        self.basis = ""
        self.status = tk.StringVar(value="")
        self._make_vars()
        self._make_menu()
        self._make_left()
        self._make_right()
        ttk.Label(root, textvariable=self.status, anchor="w", padding=(12, 4)).grid(
            row=1, column=0, columnspan=2, sticky="ew")
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)
        self.on_strategy()

    # ---- state ---------------------------------------------------------
    def _make_vars(self) -> None:
        first = STRATEGIES["rep"]
        names = dict(vector="lenti", strategy=first["name"], gene="", variant="", promoter=first["promoter"],
                     tag="none", marker=first["marker"], link="P2A", polya="SV40", resistance="AmpR",
                     hinge="CD8a hinge/TM", costim="4-1BB", guide_type="shRNA", target="", name="")
        self.var = {k: tk.StringVar(value=v) for k, v in names.items()}
        self.wpre = tk.BooleanVar(value=True)
        self.kozak = tk.BooleanVar(value=True)
        for var in list(self.var.values()) + [self.wpre, self.kozak]:
            var.trace_add("write", self.schedule)

    def strategy_key(self) -> str:
        by_name = {S["name"]: k for k, S in STRATEGIES.items()}
        return by_name.get(self.var["strategy"].get(), "rep")

    def make_config(self) -> Config:
        v = self.var
        return Config(vector=v["vector"].get(), strategy=self.strategy_key(), gene=v["gene"].get(),
                      variant=v["variant"].get(), sequence=self.insert_box.get("1.0", "end"),
                      promoter=v["promoter"].get(), marker=v["marker"].get(), link=v["link"].get(),
                      tag=v["tag"].get(), polya=v["polya"].get(), resistance=v["resistance"].get(),
                      wpre=self.wpre.get(), kozak=self.kozak.get(), hinge=v["hinge"].get(),
                      costim=v["costim"].get(), guide_type=v["guide_type"].get(),
                      target=v["target"].get(), name=v["name"].get(), overrides={}, basis=self.basis)

    # ---- widgets -------------------------------------------------------
    def _make_menu(self) -> None:
        bar = tk.Menu(self.root)
        file_menu = tk.Menu(bar, tearoff=0)
        file_menu.add_command(label="Export all files...", command=self.export)
        file_menu.add_command(label="Copy GenBank", command=lambda: self.copy_text("genbank"))
        file_menu.add_command(label="Copy FASTA", command=lambda: self.copy_text("fasta"))
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self.root.destroy)
        bar.add_cascade(label="File", menu=file_menu)
        analysis_menu = tk.Menu(bar, tearoff=0)
        analysis_menu.add_command(label="Aneuploidy analysis...", command=self.open_aneuploidy)
        bar.add_cascade(label="Analysis", menu=analysis_menu)
        help_menu = tk.Menu(bar, tearoff=0)
        help_menu.add_command(label="About", command=self.about)
        bar.add_cascade(label="Help", menu=help_menu)
        self.root.config(menu=bar)

    def _wheel(self, canvas, area) -> None:
        """Mouse-wheel scrolling for the option panel while the pointer is over it."""
        def scroll(event):
            if isinstance(event.widget, tk.Text):
                return
            if getattr(event, "num", 0) == 4:
                canvas.yview_scroll(-1, "units")
            elif getattr(event, "num", 0) == 5:
                canvas.yview_scroll(1, "units")
            elif getattr(event, "delta", 0):
                canvas.yview_scroll(-1 if event.delta > 0 else 1, "units")

        def enter(_event):
            for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                canvas.bind_all(seq, scroll)

        def leave(_event):
            for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                canvas.unbind_all(seq)

        area.bind("<Enter>", enter)
        area.bind("<Leave>", leave)

    def _section(self, parent, row: int, title: str):
        box = ttk.LabelFrame(parent, text=title, padding=(8, 4))
        box.grid(row=row, column=0, sticky="ew", pady=(0, 8))
        box.columnconfigure(0, weight=1)
        box.columnconfigure(1, weight=1)
        return box

    def _field(self, parent, row: int, col: int, text: str, var, values=None, readonly: bool = True, span: int = 1):
        """A label above a drop-down (or a text entry when values is None)."""
        frame = ttk.Frame(parent)
        frame.grid(row=row, column=col, columnspan=span, sticky="ew", padx=3, pady=2)
        frame.columnconfigure(0, weight=1)
        ttk.Label(frame, text=text).grid(row=0, column=0, sticky="w")
        if values is None:
            widget = ttk.Entry(frame, textvariable=var)
        else:
            widget = ttk.Combobox(frame, textvariable=var, values=list(values),
                                  state="readonly" if readonly else "normal")
        widget.grid(row=1, column=0, sticky="ew")
        return frame, widget

    def _make_left(self) -> None:
        outer = ttk.Frame(self.root)
        outer.grid(row=0, column=0, sticky="ns")
        canvas = tk.Canvas(outer, width=410, highlightthickness=0)
        bar = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=bar.set)
        bar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        form = ttk.Frame(canvas, padding=(12, 10))
        window = canvas.create_window((0, 0), window=form, anchor="nw")
        form.bind("<Configure>", lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window, width=e.width))
        self._wheel(canvas, outer)
        form.columnconfigure(0, weight=1)

        box = self._section(form, 0, "Vector type")
        for i, (value, label) in enumerate((("lenti", "Lentivirus"), ("aav", "AAV"), ("plasmid", "Plasmid"))):
            ttk.Radiobutton(box, text=label, value=value, variable=self.var["vector"],
                            command=self.on_vector).grid(row=0, column=i, padx=4, sticky="w")
            box.columnconfigure(i, weight=1)

        box = self._section(form, 1, "Therapeutic strategy")
        _, combo = self._field(box, 0, 0, "Strategy", self.var["strategy"], [s["name"] for s in STRATEGIES.values()], span=2)
        combo.bind("<<ComboboxSelected>>", lambda _e: self.on_strategy_chosen())
        self.tip = ttk.Label(box, wraplength=360, foreground="#587083", justify="left")
        self.tip.grid(row=1, column=0, columnspan=2, sticky="w", padx=3, pady=(0, 4))
        _, self.gene_box = self._field(box, 2, 0, "Gene (choose or type a name)", self.var["gene"], [], readonly=False)
        self._field(box, 2, 1, "Variant label (optional)", self.var["variant"])

        self.insert_frame = self._section(form, 2, "Insert sequence")
        self.insert_label = ttk.Label(self.insert_frame, wraplength=360, justify="left")
        self.insert_label.grid(row=0, column=0, columnspan=2, sticky="w", padx=3)
        self.insert_box = tk.Text(self.insert_frame, height=7, width=44, wrap="char", undo=True, font=self.mono)
        self.insert_box.grid(row=1, column=0, columnspan=2, sticky="ew", padx=3, pady=3)
        for seq in ("<KeyRelease>", "<<Paste>>", "<<Cut>>"):
            self.insert_box.bind(seq, self.schedule)
        ttk.Button(self.insert_frame, text="Load file...", command=self.load_insert).grid(row=2, column=0, sticky="w", padx=3)
        ttk.Button(self.insert_frame, text="Clear", command=self.clear_insert).grid(row=2, column=1, sticky="e", padx=3)

        self.guide_frame = self._section(form, 3, "Guide RNA")
        for i, label in enumerate(("shRNA", "sgRNA")):
            ttk.Radiobutton(self.guide_frame, text=label, value=label,
                            variable=self.var["guide_type"]).grid(row=0, column=i, padx=4, sticky="w")
        ttk.Label(self.guide_frame, wraplength=360, justify="left",
                  text="Target: shRNA sense strand (19 to 25 nt) or sgRNA protospacer (20 nt, no PAM)."
                  ).grid(row=1, column=0, columnspan=2, sticky="w", padx=3)
        self.target_entry = ttk.Entry(self.guide_frame, textvariable=self.var["target"], font=self.mono)
        self.target_entry.grid(row=2, column=0, columnspan=2, sticky="ew", padx=3, pady=3)

        self.car_frame = self._section(form, 4, "CAR domains")
        self._field(self.car_frame, 0, 0, "Hinge / transmembrane", self.var["hinge"], list(HINGES))
        self._field(self.car_frame, 0, 1, "Costimulatory domain", self.var["costim"], list(COSTIM))

        box = self._section(form, 5, "Expression")
        self.prom_row, _ = self._field(box, 0, 0, "Promoter", self.var["promoter"], list(PROMOTERS))
        self.tag_row, _ = self._field(box, 0, 1, "C-terminal tag", self.var["tag"], list(TAGS))
        self._field(box, 1, 0, "Reporter / selection", self.var["marker"], list(MARKERS))
        self.link_row, _ = self._field(box, 1, 1, "Linked by", self.var["link"], list(LINKS))
        self._field(box, 2, 0, "Polyadenylation", self.var["polya"], list(POLYA))
        self._field(box, 2, 1, "Bacterial selection", self.var["resistance"], list(RESISTANCE))
        ttk.Checkbutton(box, text="Include WPRE", variable=self.wpre).grid(row=3, column=0, sticky="w", padx=3, pady=2)
        self.kozak_check = ttk.Checkbutton(box, text="Add Kozak (GCCACC)", variable=self.kozak)
        self.kozak_check.grid(row=3, column=1, sticky="w", padx=3, pady=2)

        box = self._section(form, 6, "Vector name")
        self._field(box, 0, 0, "Leave empty for a generated name", self.var["name"], span=2)

        ttk.Label(form, text=self.BIOSAFETY, wraplength=370, foreground="#587083", justify="left"
                  ).grid(row=7, column=0, sticky="w", pady=(0, 6))

    def _make_right(self) -> None:
        right = ttk.Frame(self.root, padding=(0, 10, 12, 0))
        right.grid(row=0, column=1, sticky="nsew")
        right.columnconfigure(1, weight=1)
        right.rowconfigure(0, weight=3)
        right.rowconfigure(1, weight=2)
        self.canvas = tk.Canvas(right, width=self.MAP_SIZE, height=self.MAP_SIZE, background="white",
                                highlightthickness=1, highlightbackground="#d3dde5")
        self.canvas.grid(row=0, column=0, sticky="n", padx=(0, 12))

        table = ttk.Frame(right)
        table.grid(row=0, column=1, sticky="nsew")
        table.columnconfigure(0, weight=1)
        table.rowconfigure(0, weight=1)
        self.tree = ttk.Treeview(table, columns=("pos", "bp", "note"), height=14, selectmode="browse")
        for col, text, width, anchor in (("#0", "Element", 190, "w"), ("pos", "Position", 105, "w"),
                                         ("bp", "bp", 55, "e"), ("note", "Sequence", 210, "w")):
            self.tree.heading(col, text=text)
            self.tree.column(col, width=width, anchor=anchor, stretch=(col == "note"))
        vsb = ttk.Scrollbar(table, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        self.tree.bind("<Double-1>", self.on_row_double_click)
        self.tree.bind("<<TreeviewSelect>>", self.on_select)

        notes = ttk.Frame(right)
        notes.grid(row=1, column=0, columnspan=2, sticky="nsew", pady=(10, 0))
        notes.columnconfigure(0, weight=1)
        notes.rowconfigure(0, weight=1)
        self.msg_box = tk.Text(notes, height=9, wrap="word", state="disabled", relief="flat", padx=6, pady=4)
        msb = ttk.Scrollbar(notes, orient="vertical", command=self.msg_box.yview)
        self.msg_box.configure(yscrollcommand=msb.set)
        self.msg_box.grid(row=0, column=0, sticky="nsew")
        msb.grid(row=0, column=1, sticky="ns")
        for level, (_label, color) in self.LEVELS.items():
            self.msg_box.tag_configure(level, foreground=color)
            self.msg_box.tag_configure(level + "_label", foreground=color, font=self.bold)

        buttons = ttk.Frame(right)
        buttons.grid(row=2, column=0, columnspan=2, sticky="ew", pady=8)
        ttk.Button(buttons, text="Export files...", command=self.export).pack(side="left")
        ttk.Button(buttons, text="Copy GenBank", command=lambda: self.copy_text("genbank")).pack(side="left", padx=6)
        ttk.Button(buttons, text="Copy FASTA", command=lambda: self.copy_text("fasta")).pack(side="left")
        ttk.Button(buttons, text="Aneuploidy analysis...", command=self.open_aneuploidy).pack(side="left", padx=6)
        ttk.Label(buttons, text="Double-click a table row to paste a verified sequence for that element.",
                  foreground="#587083").pack(side="right")

    # ---- events --------------------------------------------------------
    def on_vector(self) -> None:
        vec = self.var["vector"].get()
        self.wpre.set(vec != "plasmid")
        self.var["polya"].set("bGH" if vec == "plasmid" else "SV40")

    def on_strategy(self) -> None:
        S = STRATEGIES[self.strategy_key()]
        self.gene_box["values"] = list(S["genes"])
        self.var["gene"].set(S["genes"][0])
        self.var["promoter"].set(S["promoter"])
        self.var["marker"].set(S["marker"])
        self.tip.configure(text=S["tip"])
        self.schedule()

    def schedule(self, *_args) -> None:
        """Re-draw shortly after the last change, so typing and pasting stay responsive."""
        if self.pending is not None:
            self.root.after_cancel(self.pending)
        self.pending = self.root.after(200, self.refresh)

    @staticmethod
    def _show(widget, visible: bool) -> None:
        widget.grid() if visible else widget.grid_remove()

    def update_layout(self) -> None:
        key = self.strategy_key()
        guide, car = key == "guide", key == "car"
        self._show(self.insert_frame, not guide)
        self._show(self.guide_frame, guide)
        self._show(self.car_frame, car)
        self._show(self.prom_row, not guide)
        self._show(self.tag_row, not guide and not car)
        self._show(self.link_row, not guide and self.var["marker"].get() != "none")
        self._show(self.kozak_check, not guide)
        self.insert_label.configure(
            text="scFv sequence (VH-linker-VL, no ATG or stop). Paste DNA or FASTA." if car
            else "ORF sequence (starts with ATG, no stop codon). Paste DNA or FASTA from NCBI or CCDS. "
                 "Leave empty to use a typical-length placeholder.")

    def refresh(self) -> None:
        self.pending = None
        self.update_layout()
        cfg = self.make_config()
        try:
            base = build_design(dataclasses.replace(cfg, overrides={}))
            ids = {f.id for f in base.features} - BUILT_IN_IDS
            cfg.overrides = {k: v for k, v in self.overrides.items() if k in ids}
            self.design = build_design(cfg)
        except ValueError as exc:
            self.set_status("Error: %s" % exc)
            return
        d = self.design
        self.summary_line = "%s  |  %s bp  |  %d elements" % (d.name, format(d.total, ","), len(d.features))
        self.has_insert = bool(clean_sequence(cfg.sequence)[0])
        self.has_target = bool(clean_sequence(cfg.target)[0])
        self.set_status(self.summary_line)
        self.draw_map()
        self.fill_tree()
        self.fill_messages()

    def set_status(self, text: str) -> None:
        self.status.set(text)

    def draw_map(self) -> None:
        c, d, size = self.canvas, self.design, self.MAP_SIZE
        c.delete("all")
        self.arcs = {}
        box = (size / 2 - 140, size / 2 - 140, size / 2 + 140, size / 2 + 140)
        c.create_oval(*box, outline="#d3dde5", width=30)
        for f in d.features:
            start = 360.0 * (f.start - 1) / d.total
            extent = max(360.0 * f.length / d.total, 1.5)
            item = c.create_arc(*box, start=90 - start, extent=-extent, style="arc",
                                outline=KINDS[f.kind][0], width=26)
            self.arcs[f.id] = item
            c.tag_bind(item, "<Enter>", lambda _e, f=f: self.set_status(
                "%s   %d..%d   %d bp" % (f.name, f.start, f.end, f.length)))
            c.tag_bind(item, "<Leave>", lambda _e: self.set_status(self.summary_line))
        c.create_line(size / 2, 14, size / 2, 30, fill="#14202b", width=2)
        name = d.name if len(d.name) <= 40 else d.name[:39] + "..."
        c.create_text(size / 2, size / 2 - 10, text=name, width=200, justify="center", font=self.bold, fill="#14202b")
        c.create_text(size / 2, size / 2 + 28, text="%s bp" % format(d.total, ","), fill="#587083")

    def swatch(self, kind: str):
        if kind not in self.swatches:
            img = tk.PhotoImage(width=12, height=12)
            img.put(KINDS[kind][0], to=(0, 0, 12, 12))
            self.swatches[kind] = img
        return self.swatches[kind]

    def fill_tree(self) -> None:
        previous = self.tree.selection()
        self.tree.delete(*self.tree.get_children())
        for f in self.design.features:
            self.tree.insert("", "end", iid=f.id, text=f.name, image=self.swatch(f.kind),
                             values=("%d..%d" % (f.start, f.end), f.length,
                                     element_note(f, self.has_insert, self.has_target)))
        if previous and self.tree.exists(previous[0]):
            self.tree.selection_set(previous[0])

    def fill_messages(self) -> None:
        box = self.msg_box
        box.configure(state="normal")
        box.delete("1.0", "end")
        for level, text in self.design.messages:
            box.insert("end", self.LEVELS[level][0] + ": ", (level + "_label",))
            box.insert("end", text + "\n", (level,))
        box.configure(state="disabled")

    def on_select(self, _event=None) -> None:
        for item in self.arcs.values():
            self.canvas.itemconfigure(item, width=26)
        selected = self.tree.selection()
        if selected and selected[0] in self.arcs:
            self.canvas.itemconfigure(self.arcs[selected[0]], width=40)

    def on_row_double_click(self, event) -> None:
        fid = self.tree.identify_row(event.y)
        if not fid:
            return
        if fid in ("goi", "guide"):
            (self.target_entry if fid == "guide" else self.insert_box).focus_set()
            self.set_status("Paste the sequence for this element in the box on the left.")
        elif fid in BUILT_IN_IDS:
            self.set_status("This element has a built-in sequence.")
        else:
            self.edit_element(fid)

    # ---- actions -------------------------------------------------------
    def edit_element(self, fid: str) -> None:
        feature = next(f for f in self.design.features if f.id == fid)
        top = tk.Toplevel(self.root)
        top.title("Sequence for " + feature.name)
        top.transient(self.root)
        top.grab_set()
        ttk.Label(top, wraplength=440, justify="left", padding=8,
                  text="Paste a verified sequence for %s (plain DNA or FASTA). Leave it empty to go back to "
                       "the N placeholder." % feature.name).pack(anchor="w")
        box = tk.Text(top, width=64, height=10, wrap="char", font=self.mono)
        box.pack(padx=8, fill="both", expand=True)
        box.insert("1.0", self.overrides.get(fid, ""))
        box.focus_set()

        def load():
            path = filedialog.askopenfilename(title="Open a sequence file", filetypes=[
                ("Sequence files", "*.fasta *.fa *.txt"), ("All files", "*.*")])
            if path:
                with open(path, "r", encoding="utf-8", errors="replace") as fh:
                    box.delete("1.0", "end")
                    box.insert("1.0", fh.read())

        def accept():
            seq, _bad = clean_sequence(box.get("1.0", "end"))
            if seq:
                self.overrides[fid] = seq
            else:
                self.overrides.pop(fid, None)
            top.destroy()
            self.schedule()

        row = ttk.Frame(top, padding=8)
        row.pack(fill="x")
        ttk.Button(row, text="Load file...", command=load).pack(side="left")
        ttk.Button(row, text="Cancel", command=top.destroy).pack(side="right")
        ttk.Button(row, text="Use sequence", command=accept).pack(side="right", padx=6)

    def load_insert(self) -> None:
        path = filedialog.askopenfilename(title="Open a sequence file", filetypes=[
            ("Sequence files", "*.fasta *.fa *.txt *.gb"), ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                text = fh.read()
        except OSError as exc:
            messagebox.showerror("Could not open file", str(exc))
            return
        self.insert_box.delete("1.0", "end")
        self.insert_box.insert("1.0", text)
        self.schedule()

    def clear_insert(self) -> None:
        self.insert_box.delete("1.0", "end")
        self.schedule()

    def export(self) -> None:
        if self.design is None:
            return
        folder = filedialog.askdirectory(title="Choose a folder for the exported files")
        if not folder:
            return
        try:
            paths = write_outputs(self.design, folder)
        except OSError as exc:
            messagebox.showerror("Could not write files", str(exc))
            return
        messagebox.showinfo("Files written", "\n".join(paths))
        self.set_status("Exported %d files to %s" % (len(paths), folder))

    def copy_text(self, kind: str) -> None:
        if self.design is None:
            return
        text = to_genbank(self.design) if kind == "genbank" else to_fasta(self.design)
        self.root.clipboard_clear()
        self.root.clipboard_append(text)
        self.set_status("Copied %s to the clipboard." % ("GenBank" if kind == "genbank" else "FASTA"))

    def on_strategy_chosen(self) -> None:
        """The user picked a strategy by hand, so any earlier aneuploidy-based reason no longer applies."""
        self.basis = ""
        self.on_strategy()

    def open_aneuploidy(self) -> None:
        AneuploidyWindow(self)

    def apply_recommendation(self, rec: Recommendation, basis: str) -> None:
        """Set the designer to a suggestion from the aneuploidy analysis."""
        self.var["strategy"].set(STRATEGIES[rec.strategy]["name"])
        self.on_strategy()
        self.var["gene"].set(rec.gene)
        if rec.promoter:
            self.var["promoter"].set(rec.promoter)
        self.basis = basis
        self.schedule()
        self.set_status("Designer set from the aneuploidy analysis: %s (%s)." % (rec.title, rec.gene))
        self.root.lift()

    def about(self) -> None:
        messagebox.showinfo("About", "Cancer Gene Therapy Vector Designer %s\n\n%s\n\nLibrary lengths are typical "
                            "values. Elements shown as N need verified sequences before cloning."
                            % (__version__, self.BIOSAFETY))
