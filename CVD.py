#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cancer Gene Therapy Vector Designer
===================================

Design viral vectors (lentivirus, AAV or plain plasmid) around cancer-relevant
inserts, check the insert sequence, and export a GenBank file, a FASTA file, a
circular SVG map and a text summary.  Standard library only (Python 3.8+).

USAGE
    python vector_designer.py                     open the graphical window
    python vector_designer.py --wizard            answer text prompts instead
    python vector_designer.py --list              show strategies, genes and parts
    python vector_designer.py [options]           one-shot design (see --help)
    python vector_designer.py --cnv segments.tsv  aneuploidy analysis and design suggestions
    python vector_designer.py --cnv segments.tsv --design-from-cnv 1
                                                  ...and build the vector for suggestion 1
    python vector_designer.py --example-cnv demo.tsv    write a synthetic input file to try

    The window needs tkinter, which ships with the python.org installers for
    Windows and macOS.  On Debian/Ubuntu install it with: sudo apt install python3-tk
    Without tkinter the program falls back to the text prompts.

EXAMPLES
    # p53 replacement in AAV, ORF read from a FASTA file
    python vector_designer.py --vector aav --strategy rep --gene TP53 \\
        --seq tp53_orf.fasta --promoter EF1a --out designs

    # CAR T construct against CD19 in a lentivirus, scFv pasted from a file
    python vector_designer.py --vector lenti --strategy car --gene CD19 \\
        --seq cd19_scfv.txt --hinge "CD8a hinge/TM" --costim 4-1BB

    # shRNA knockdown of KRAS with a puromycin marker
    python vector_designer.py --strategy guide --gene KRAS --guide-type shRNA \\
        --target GTTGGAGCTGATGGCGTAGG --marker Puro

ANEUPLOIDY AND IMMUNE FLAGGING
    --cnv reads a copy-number segment file (chrom, start, end and a copy_number or
    log2 column; DNAcopy and most CNV caller outputs work), calls gained and lost
    chromosome arms (an arm is called at 80% or more of its length, as in Taylor et
    al., Cancer Cell 2018), computes an aneuploidy score, and suggests designs: an
    immune "flag" payload (NKG2D ligand, surface tag for cetuximab, calreticulin),
    checkpoint or don't-eat-me knockdown, cytokine support, and tumor suppressor
    replacement or oncogene knockdown for the altered arms.  The suggestions come
    from known gene locations and published associations; they do not predict
    response.  The vector itself cannot sense aneuploidy: selectivity comes from
    delivery and the tumor-selective promoter, and aneuploid cells can also occur
    in some normal tissues.

IMPORTANT
    Library genes and most vector elements use typical lengths and are written
    as runs of N until you supply verified sequences (--seq for the insert,
    --override ID=SEQUENCE for other elements, for example from an Addgene
    reference plasmid).  This is a research design tool.  Producing virus needs
    institutional biosafety approval, and a clinical gene therapy product needs
    GMP manufacturing, preclinical safety data and regulatory authorization.

BUILD A STAND-ALONE EXECUTABLE (optional)
    pip install pyinstaller
    pyinstaller --onefile --windowed --name vector_designer vector_designer.py
    The program appears in the dist/ folder.  --windowed hides the console, which
    suits the graphical program; leave it out to keep the command-line output.
"""
from __future__ import annotations

import argparse
import dataclasses
import datetime
import math
import os
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from xml.sax.saxutils import escape as xml_escape

try:  # the command-line modes keep working when tkinter is not installed
    import tkinter as tk
    import tkinter.font as tkfont
    from tkinter import filedialog, messagebox, ttk
except ImportError:
    tk = None

__version__ = "1.0.0"

# --------------------------------------------------------------------------
# Parts library
# --------------------------------------------------------------------------
VECTORS = ["lenti", "aav", "plasmid"]

PROMOTERS = {"CMV": 508, "EF1a": 1178, "EFS": 212, "CAG": 1733, "PGK": 516,
             "UbC": 1220, "hSyn": 469, "SV40": 330, "MND": 255,
             "hTERT": 450, "Survivin": 270}
TUMOR_SELECTIVE = {"hTERT", "Survivin"}
MARKERS = {"none": 0, "EGFP": 720, "mCherry": 711, "Puro": 600, "Neo": 795,
           "Blast": 399, "Hygro": 1026}
POLYA = {"bGH": 225, "SV40": 122, "none": 0}
RESISTANCE = {"AmpR": 861, "KanR": 816}
TAGS = {"none": "",
        "FLAG3": "GACTACAAAGACCATGACGGTGATTATAAAGATCATGACATCGACTACAAGGATGACGATGACAAG",
        "HA": "TACCCATACGATGTTCCAGATTACGCT",
        "Myc": "GAACAAAAACTCATCTCAGAAGAGGATCTG",
        "His6": "CATCACCATCACCATCAC"}
LINKS = {"P2A": 0, "IRES": 580, "separate PGK promoter": 516}
HINGES = {"CD8a hinge/TM": 204, "CD28 hinge/TM": 201}
COSTIM = {"4-1BB": 126, "CD28": 123, "none": 0}
P2A_SEQ = "GGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAGACGTGGAGGAGAACCCTGGACCT"
SGRNA_SCAFFOLD = ("GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC")
SHRNA_LOOP = "TTCAAGAGA"

# typical ORF length (bp, without stop codon) and a short note
GENES: Dict[str, Tuple[int, str]] = {
    "Custom": (1000, ""),
    "TP53": (1182, "TP53 is the most frequently mutated tumor suppressor. Use the wild-type sequence for replacement."),
    "PTEN": (1212, ""), "RB1": (2787, ""), "CDKN2A": (471, "Encodes p16INK4a."),
    "VHL": (639, ""), "SMAD4": (1659, ""), "STK11": (1302, ""),
    "CDKN1A": (495, "Encodes p21."),
    "BRCA1": (5592, "Large ORF: check the capacity message for your vector type."),
    "APC": (8532, "Very large ORF: too big for AAV and difficult for lentivirus."),
    "NF1": (8520, "Very large ORF: too big for AAV and difficult for lentivirus."),
    "HSV-TK": (1131, "Pair with ganciclovir to kill transduced cells."),
    "iCasp9": (1170, "The dimerizer AP1903 (rimiducid) triggers apoptosis."),
    "CD::UPRT": (1100, "Converts 5-fluorocytosine (5-FC) into toxic 5-FU."),
    "IL-12 (single chain)": (1600, "p40-linker-p35 fusion."),
    "IL-15": (486, ""), "IL-2": (462, ""), "GM-CSF": (435, ""),
    "CD19": (750, ""), "BCMA": (750, ""), "HER2": (750, ""), "GD2": (750, ""),
    "Mesothelin": (750, ""), "EGFRvIII": (750, ""), "CD22": (750, ""),
    "KRAS": (570, ""), "MYC": (1311, ""), "BRAF": (2301, ""), "EGFR": (3633, ""),
    "PIK3CA": (3207, ""), "CTNNB1": (2346, ""), "BCL2": (720, ""),
    # immune-flag payloads (typical ORF lengths without the stop codon)
    "ULBP2": (738, "NKG2D ligand. Marks cells for NK and CD8 T-cell attack without needing MHC-I. Proteases can shed "
                   "NKG2D ligands from the surface and the soluble form can blunt NKG2D, so check surface retention."),
    "ULBP1": (732, "NKG2D ligand, GPI-anchored. Marks cells for NK and CD8 T-cell attack without needing MHC-I."),
    "MICA": (1149, "NKG2D ligand. Proteases (ADAM10, ADAM17) can shed it from tumor cells, which blunts the effect."),
    "MICB": (1149, "NKG2D ligand. Can be shed from tumor cells in the same way as MICA."),
    "CALR": (1251, "Calreticulin. Surface calreticulin is an 'eat me' signal for macrophages. Often paired with CD47 knockdown."),
    "tEGFR (surface tag)": (1071, "Truncated EGFR without its signaling domain. Cetuximab binds it and can kill cells "
                                  "through ADCC and CDC. Length is approximate; paste your sequence."),
    "Bispecific T-cell engager (secreted)": (1600, "Secreted anti-CD3 by tumor-antigen engager. Paste your scFv sequences; "
                                                   "length is approximate."),
    "CD47": (969, "Don't-eat-me signal for macrophages. Knock it down with an shRNA."),
    "CD274": (870, "PD-L1. Knock it down with an shRNA to release T cells."),
}

STRATEGIES = {
    "rep": dict(name="Tumor suppressor replacement",
                genes=["TP53", "PTEN", "RB1", "CDKN2A", "VHL", "SMAD4", "STK11",
                       "CDKN1A", "BRCA1", "APC", "NF1", "Custom"],
                promoter="EF1a", marker="none",
                tip="Restores a lost tumor suppressor. Use the wild-type ORF. AAV suits in vivo delivery; lentivirus suits ex vivo work."),
    "sui": dict(name="Suicide gene and prodrug",
                genes=["HSV-TK", "iCasp9", "CD::UPRT", "Custom"],
                promoter="Survivin", marker="none",
                tip="Tumor-selective promoters (hTERT, Survivin) limit expression to cancer cells. A prodrug or dimerizer then kills the transduced cells."),
    "cyt": dict(name="Immunostimulatory cytokine",
                genes=["IL-12 (single chain)", "IL-15", "IL-2", "GM-CSF", "Custom"],
                promoter="EF1a", marker="none",
                tip="Local cytokine expression recruits and activates immune cells. Watch for systemic toxicity."),
    "flag": dict(name="Immune flagging of aneuploid tumor cells",
                 genes=["ULBP2", "ULBP1", "MICA", "MICB", "CALR", "tEGFR (surface tag)",
                        "Bispecific T-cell engager (secreted)", "Custom"],
                 promoter="Survivin", marker="none",
                 tip="Expresses a protein that marks tumor cells for immune attack: an NKG2D ligand (NK and CD8 T cells), "
                     "calreticulin (macrophage 'eat me' signal), a surface tag for cetuximab, or a secreted T-cell engager. "
                     "Run the aneuploidy analysis to choose the payload. The vector cannot sense aneuploidy itself; "
                     "selectivity comes from delivery and the tumor-selective promoter."),
    "car": dict(name="CAR T-cell therapy",
                genes=["CD19", "BCMA", "HER2", "GD2", "Mesothelin", "EGFRvIII", "CD22", "Custom"],
                promoter="MND", marker="none",
                tip="Lentivirus is the usual choice for ex vivo T-cell engineering. Supply your scFv sequence and pick the hinge and costimulatory domains."),
    "guide": dict(name="Oncogene knockdown or editing",
                  genes=["KRAS", "MYC", "BCL2", "BRAF", "EGFR", "CTNNB1", "PIK3CA", "CD274", "CD47", "Custom"],
                  promoter="EF1a", marker="Puro",
                  tip="U6-driven vector to knock down (shRNA) or edit (sgRNA) an oncogene, a mutant allele or an immune "
                      "checkpoint or don't-eat-me gene (CD274, CD47). Supply your target sequence."),
    "model": dict(name="Cancer modeling (research use)",
                  genes=["KRAS", "MYC", "BRAF", "EGFR", "PIK3CA", "CTNNB1", "BCL2", "Custom"],
                  promoter="EF1a", marker="Puro",
                  tip="Oncogene expression for transformation assays and animal models. Use the variant label for hotspot mutations."),
}

RESTRICTION_SITES = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
                     "XhoI": "CTCGAG", "NotI": "GCGGCCGC", "NheI": "GCTAGC",
                     "AgeI": "ACCGGT", "XbaI": "TCTAGA", "KpnI": "GGTACC",
                     "MluI": "ACGCGT", "SalI": "GTCGAC", "PacI": "TTAATTAA"}

# element kind -> (SVG color, GenBank feature key)
KINDS = {
    "ltr": ("#6b7f90", "LTR"), "itr": ("#6b7f90", "repeat_region"),
    "pack": ("#9b7bd8", "misc_feature"), "promoter": ("#e8a317", "promoter"),
    "cds": ("#1fa971", "CDS"), "tag": ("#18a5b8", "misc_feature"),
    "marker": ("#e0483f", "CDS"), "reg": ("#7a63d1", "misc_feature"),
    "polya": ("#8497a5", "polyA_signal"), "ori": ("#a7b5c0", "rep_origin"),
    "bact": ("#d6479a", "CDS"), "link": ("#3b82d6", "misc_feature"),
    "stop": ("#c0392b", "misc_feature"), "guide": ("#e56b2f", "misc_feature"),
}
# these elements get their sequence from the insert options or are built in
BUILT_IN_IDS = {"goi", "guide", "kozak", "stop", "tag", "p2a", "term"}


# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------
@dataclass
class Config:
    vector: str = "lenti"
    strategy: str = "rep"
    gene: str = ""
    variant: str = ""
    sequence: str = ""            # insert ORF / scFv (text, FASTA allowed)
    promoter: Optional[str] = None
    marker: Optional[str] = None
    link: str = "P2A"
    tag: str = "none"
    polya: Optional[str] = None
    resistance: str = "AmpR"
    wpre: Optional[bool] = None
    kozak: bool = True
    hinge: str = "CD8a hinge/TM"
    costim: str = "4-1BB"
    guide_type: str = "shRNA"
    target: str = ""
    name: str = ""
    overrides: Dict[str, str] = field(default_factory=dict)
    basis: str = ""               # why this design was chosen; shown in the summary


@dataclass
class Feature:
    id: str
    name: str
    kind: str
    length: int
    seq: Optional[str] = None
    in_cassette: bool = False
    custom: bool = False
    start: int = 0
    end: int = 0


@dataclass
class Design:
    name: str
    vector: str
    strategy: str
    features: List[Feature]
    seq: str
    messages: List[Tuple[str, str]]   # (level, text); level: error / warn / info
    basis: str = ""

    @property
    def total(self) -> int:
        return len(self.seq)


# --------------------------------------------------------------------------
# Sequence helpers
# --------------------------------------------------------------------------
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def clean_sequence(text: str) -> Tuple[str, int]:
    """Strip FASTA headers, whitespace and digits; U -> T. Returns (ACGTN, removed)."""
    text = re.sub(r"^>.*$", "", text or "", flags=re.M)
    text = re.sub(r"[\s\d]", "", text).upper().replace("U", "T")
    kept = re.sub(r"[^ACGTN]", "", text)
    return kept, len(text) - len(kept)


def revcomp(seq: str) -> str:
    return "".join(_COMP[c] for c in reversed(seq))


def gc_percent(seq: str) -> int:
    return round(100 * sum(c in "GC" for c in seq) / len(seq)) if seq else 0


def count_overlapping(seq: str, motif: str) -> int:
    n, i = 0, seq.find(motif)
    while i >= 0:
        n += 1
        i = seq.find(motif, i + 1)
    return n


def pick(value: str, options, what: str) -> str:
    """Case-insensitive choice from a dict/list, with a helpful error."""
    for opt in options:
        if str(opt).lower() == str(value).lower():
            return opt
    raise ValueError("Unknown %s %r. Choose from: %s" % (what, value, ", ".join(map(str, options))))


# --------------------------------------------------------------------------
# Design logic
# --------------------------------------------------------------------------
def check_insert(msgs: list, seq: str, label: str, vector: str, check_atg: bool) -> None:
    n = len(seq)
    if n % 3:
        msgs.append(("warn", "%s length (%d bp) is not a multiple of 3." % (label, n)))
    else:
        stops, first = 0, 0
        for i in range(0, n - 3, 3):
            if seq[i:i + 3] in ("TAA", "TAG", "TGA"):
                stops += 1
                first = first or i // 3 + 1
        if stops:
            msgs.append(("error", "%d in-frame stop codon(s) inside %s; the first is at codon %d. "
                         "Check the reading frame, or confirm the premature stop is intended." % (stops, label, first)))
    if check_atg and not seq.startswith("ATG"):
        msgs.append(("warn", "%s does not start with ATG." % label))
    if n % 3 == 0 and seq[-3:] in ("TAA", "TAG", "TGA"):
        msgs.append(("warn", "%s ends in a stop codon. Remove it if you use a tag or a P2A reporter." % label))
    msgs.append(("info", "%s: %s bp, %d%% GC." % (label, format(n, ","), gc_percent(seq))))
    if vector == "lenti" and re.search("AATAAA|ATTAAA", seq):
        msgs.append(("warn", "%s contains a polyadenylation signal (AATAAA or ATTAAA). In a lentiviral genome "
                     "this can truncate the vector RNA and lower titer." % label))
    free, cut = [], []
    for enzyme, site in RESTRICTION_SITES.items():
        c = count_overlapping(seq, site)
        (cut.append("%s x%d" % (enzyme, c)) if c else free.append(enzyme))
    text = "Restriction sites in %s: none for %s." % (label, ", ".join(free) or "any listed enzyme")
    if cut:
        text += " Cut by %s." % ", ".join(cut)
    msgs.append(("info", text))


def build_design(cfg: Config) -> Design:
    vector = pick(cfg.vector, VECTORS, "vector")
    strat_key = pick(cfg.strategy, STRATEGIES, "strategy")
    S = STRATEGIES[strat_key]
    is_guide, is_car = strat_key == "guide", strat_key == "car"

    prom = pick(cfg.promoter or S["promoter"], PROMOTERS, "promoter")
    marker = pick(cfg.marker or S["marker"], MARKERS, "marker")
    link = "separate PGK promoter" if is_guide else pick(cfg.link, LINKS, "link")
    tag = "none" if (is_guide or is_car) else pick(cfg.tag, TAGS, "tag")
    polya = pick(cfg.polya or ("bGH" if vector == "plasmid" else "SV40"), POLYA, "polyadenylation signal")
    res = pick(cfg.resistance, RESISTANCE, "bacterial selection")
    wpre = cfg.wpre if cfg.wpre is not None else vector != "plasmid"
    hinge = pick(cfg.hinge, HINGES, "hinge")
    costim = pick(cfg.costim, COSTIM, "costimulatory domain")
    sgrna = pick(cfg.guide_type, ["shRNA", "sgRNA"], "guide type") == "sgRNA"

    gene_in = (cfg.gene or "").strip()
    lib_name = next((g for g in GENES if g.lower() == gene_in.lower()
                     or g.lower().startswith(gene_in.lower() + " (")), None)
    gene_info = GENES[lib_name] if lib_name and lib_name != "Custom" else GENES["Custom"]
    gene = lib_name if lib_name and lib_name != "Custom" else (gene_in or ("Target" if is_guide else "Transgene"))
    if cfg.variant and not is_guide:
        gene += " " + cfg.variant.strip()

    insert, bad = (("", 0) if is_guide else clean_sequence(cfg.sequence))
    target = clean_sequence(cfg.target)[0].replace("N", "") if is_guide else ""
    msgs: List[Tuple[str, str]] = []

    feats: List[Feature] = []
    cassette = False

    def add(fid, name, kind, length, seq=None):
        feats.append(Feature(fid, name, kind, length, seq, cassette))

    if vector == "lenti":
        add("ltr5", "5' LTR", "ltr", 634)
        add("psi", "Psi packaging signal", "pack", 124)
        add("rre", "RRE", "pack", 234)
        add("cppt", "cPPT/CTS", "pack", 118)
    if vector == "aav":
        add("itr5", "5' ITR", "itr", 145)
    cassette = True

    p2a = marker != "none" and link == "P2A" and not is_guide
    if is_guide:
        add("u6", "U6 promoter", "promoter", 241)
    else:
        add("prom", prom + " promoter", "promoter", PROMOTERS[prom])
        if cfg.kozak:
            add("kozak", "Kozak", "reg", 6, "GCCACC")

    if is_guide:
        n = target or "N" * (20 if sgrna else 21)
        hairpin = n + SGRNA_SCAFFOLD if sgrna else n + SHRNA_LOOP + revcomp(n)
        add("guide", gene + (" sgRNA" if sgrna else " shRNA"), "guide", len(hairpin), hairpin)
        add("term", "Pol III terminator", "polya", 6, "TTTTTT")
    elif is_car:
        add("lead", "CD8a signal peptide", "tag", 63)
        add("goi", "anti-%s scFv" % gene, "cds", len(insert) or gene_info[0], insert or None)
        add("hinge", hinge, "tag", HINGES[hinge])
        if costim != "none":
            add("costim", costim + " costimulatory domain", "tag", COSTIM[costim])
        add("cd3z", "CD3zeta signaling domain", "tag", 336)
    else:
        add("goi", gene, "cds", len(insert) or gene_info[0], insert or None)

    if tag != "none":
        add("tag", tag + " tag", "tag", len(TAGS[tag]), TAGS[tag])
    if not p2a and not is_guide:
        add("stop", "Stop (TAA)", "stop", 3, "TAA")
    if marker != "none" and not is_guide:
        if link == "P2A":
            add("p2a", "P2A", "link", 66, P2A_SEQ)
        if link == "IRES":
            add("ires", "IRES (EMCV)", "link", 580)
        add("mark", marker, "marker", MARKERS[marker])
        if p2a:
            add("stop", "Stop (TAA)", "stop", 3, "TAA")
    if wpre and not is_guide:
        add("wpre", "WPRE", "reg", 589)
    if polya != "none" and not is_guide:
        add("pa", polya + " polyA", "polya", POLYA[polya])
    if marker != "none" and link.startswith("separate"):
        add("pgk", "PGK promoter", "promoter", 516)
        add("mark2", marker + " (PGK)", "marker", MARKERS[marker])
    cassette = False
    if vector == "lenti":
        add("ltr3", "3' dU3 LTR", "ltr", 234)
    if vector == "aav":
        add("itr3", "3' ITR", "itr", 145)
    if vector != "plasmid":
        add("f1", "f1 ori", "ori", 429)
    add("res", res, "bact", RESISTANCE[res])
    add("ori", "pUC ori", "ori", 589)

    # user overrides (verified sequences for placeholder elements)
    ids = {f.id for f in feats}
    for fid, raw in cfg.overrides.items():
        if fid in BUILT_IN_IDS:
            raise ValueError("Element %r takes its sequence from --seq, --target or the built-in library; "
                             "it cannot be overridden." % fid)
        if fid not in ids:
            raise ValueError("No element with id %r in this design. Ids in this design: %s"
                             % (fid, ", ".join(sorted(ids - BUILT_IN_IDS))))
        s, _ = clean_sequence(raw)
        if not s:
            raise ValueError("Override for %r contains no A, C, G, T or N." % fid)
        f = next(f for f in feats if f.id == fid)
        f.seq, f.length, f.custom = s, len(s), True

    pos, parts = 1, []
    for f in feats:
        if not f.seq:
            f.seq = "N" * f.length
        f.start, f.end = pos, pos + f.length - 1
        pos += f.length
        parts.append(f.seq)
    seq = "".join(parts)

    # ---- design checks -------------------------------------------------
    if not is_guide and not insert:
        msgs.append(("warn", "No sequence supplied, so a %s bp placeholder (typical length) is used for %s."
                     % (format(gene_info[0], ","), gene)))
    if bad:
        msgs.append(("warn", "%d non-nucleotide character(s) were removed from the insert sequence." % bad))
    if not is_guide and insert:
        check_insert(msgs, insert, "scFv" if is_car else "Insert", vector, check_atg=not is_car)
    if not is_guide and not is_car and gene_info[1]:
        msgs.append(("info", gene_info[1]))
    if is_guide:
        if not target:
            msgs.append(("warn", "No target sequence supplied, so an N-filled placeholder is used."))
        else:
            if sgrna and len(target) != 20:
                msgs.append(("warn", "sgRNA protospacers are usually 20 nt (yours: %d)." % len(target)))
            if not sgrna and not 19 <= len(target) <= 25:
                msgs.append(("warn", "shRNA sense strands are usually 19 to 25 nt (yours: %d)." % len(target)))
            if "TTTT" in target:
                msgs.append(("error", "The target contains TTTT, which ends Pol III transcription early. Choose another site."))
            if not sgrna and not 30 <= gc_percent(target) <= 60:
                msgs.append(("warn", "Target GC content is %d%%. Most shRNA designs aim for 30 to 55%%." % gc_percent(target)))
            if target[0] != "G":
                msgs.append(("info", "U6 transcribes best from a 5' G. Consider adding one to the 5' end."))
        if sgrna:
            msgs.append(("info", "This vector carries the sgRNA only. Pair it with Cas9 (SaCas9 fits AAV)."))
    if strat_key == "model":
        msgs.append(("info", "Research use only: oncogene expression for transformation assays and animal models. "
                             "This is not a therapeutic design."))
    if is_car:
        msgs.append(("info", "CAR domains use typical lengths and are written as N until you supply verified sequences. "
                             "Keep the scFv, hinge, costimulatory and CD3zeta domains in frame."))
    if strat_key == "flag":
        msgs.append(("info", "This vector cannot sense aneuploidy. Selectivity comes from delivery and the "
                             "tumor-selective promoter, not from the payload."))
        msgs.append(("warn", "Aneuploid cells can also occur in some normal tissues. Test expression and immune killing "
                             "in normal cells before trusting the selectivity."))
    if p2a:
        msgs.append(("info", "P2A selected: keep the insert (and tag) in frame with P2A and the reporter."))
    if prom in TUMOR_SELECTIVE and not is_guide:
        msgs.append(("info", "%s is a tumor-selective promoter. Expression is limited to cells where it is active "
                             "and is usually weaker than EF1a or CMV, so confirm activity in your cell model." % prom))
    if vector != "plasmid":
        cas = sum(f.length for f in feats if f.in_cassette)
        n = format(cas, ",")
        if vector == "aav":
            if cas > 4700:
                msgs.append(("error", "Cassette between the ITRs is %s bp, above the roughly 4.7 kb AAV packaging limit." % n))
            elif cas > 4400:
                msgs.append(("warn", "Cassette between the ITRs is %s bp, close to the roughly 4.7 kb AAV limit." % n))
            else:
                msgs.append(("info", "Cassette between the ITRs: %s bp (limit about 4.7 kb)." % n))
        else:
            if cas > 8000:
                msgs.append(("warn", "Insert region between the LTRs is %s bp. Large inserts lower lentiviral titer." % n))
            else:
                msgs.append(("info", "Insert region between the LTRs: %s bp." % n))
        if vector == "lenti":
            msgs.append(("info", "Lentiviral vectors integrate. This design uses a self-inactivating (dU3) 3' LTR."))

    msgs.sort(key=lambda m: ("error", "warn", "info").index(m[0]))

    prefix = {"lenti": "pLV", "aav": "pAAV", "plasmid": "pEx"}[vector]
    tagpart = "-" + tag if tag != "none" else ""
    markpart = ""
    if marker != "none":
        markpart = "-%s-%s" % ({"P2A": "P2A", "IRES": "IRES"}.get(link, "PGK"), marker)
    lead = ("sg" if sgrna else "sh") if is_guide else ("CAR-" if is_car else "")
    short = re.sub(r"\s+", "_", re.sub(r"\s*\([^)]*\)", "", gene).strip())
    auto = "%s-%s>%s%s%s%s" % (prefix, "U6" if is_guide else prom, lead, short, tagpart, markpart)
    return Design(cfg.name.strip() or auto, vector, strat_key, feats, seq, msgs, cfg.basis)


# --------------------------------------------------------------------------
# Aneuploidy analysis
# --------------------------------------------------------------------------
# (chromosome length, approximate centromere position) in bp. Centromere positions are
# approximate (within a few Mb). Arm calls need 80% of an arm, so this does not change them.
BUILDS = {
    "hg38": {"1": (248956422, 123400000), "2": (242193529, 93900000), "3": (198295559, 90900000),
             "4": (190214555, 50000000), "5": (181538259, 48800000), "6": (170805979, 59800000),
             "7": (159345973, 60100000), "8": (145138636, 45200000), "9": (138394717, 43000000),
             "10": (133797422, 39800000), "11": (135086622, 53400000), "12": (133275309, 35500000),
             "13": (114364328, 17700000), "14": (107043718, 17200000), "15": (101991189, 19000000),
             "16": (90338345, 36800000), "17": (83257441, 25100000), "18": (80373285, 18500000),
             "19": (58617616, 26200000), "20": (64444167, 28100000), "21": (46709983, 12000000),
             "22": (50818468, 15000000)},
    "hg19": {"1": (249250621, 125000000), "2": (243199373, 93300000), "3": (198022430, 91000000),
             "4": (191154276, 50400000), "5": (180915260, 48400000), "6": (171115067, 61000000),
             "7": (159138663, 59900000), "8": (146364022, 45600000), "9": (141213431, 49000000),
             "10": (135534747, 40200000), "11": (135006516, 53700000), "12": (133851895, 35800000),
             "13": (115169878, 17900000), "14": (107349540, 17600000), "15": (102531392, 19000000),
             "16": (90354753, 36600000), "17": (81195210, 24000000), "18": (78077248, 17200000),
             "19": (59128983, 26500000), "20": (63025520, 27500000), "21": (48129895, 13200000),
             "22": (51304566, 14700000)},
}
ACROCENTRIC_P = {"13p", "14p", "15p", "21p", "22p"}          # not scored (no usable p-arm sequence)
AUTOSOMES = [str(i) for i in range(1, 23)]
ARM_ORDER = [c + s for c in AUTOSOMES for s in "pq"]
N_SCORED_ARMS = len([a for a in ARM_ORDER if a not in ACROCENTRIC_P])   # 39

# well-known cancer genes by chromosome arm: tumor suppressors (tsg) and oncogenes (onc)
ARM_GENES = {"TP53": ("17p", "tsg"), "PTEN": ("10q", "tsg"), "RB1": ("13q", "tsg"), "CDKN2A": ("9p", "tsg"),
             "VHL": ("3p", "tsg"), "SMAD4": ("18q", "tsg"), "STK11": ("19p", "tsg"), "CDKN1A": ("6p", "tsg"),
             "BRCA1": ("17q", "tsg"), "APC": ("5q", "tsg"), "NF1": ("17q", "tsg"),
             "KRAS": ("12p", "onc"), "MYC": ("8q", "onc"), "BRAF": ("7q", "onc"), "EGFR": ("7p", "onc"),
             "PIK3CA": ("3q", "onc"), "CTNNB1": ("3p", "onc"), "BCL2": ("18q", "onc")}

# arm, direction of change, what lives there, why it matters for immune attack
IMMUNE_ARM_NOTES = [
    ("6p", "loss", "HLA class I region (6p21.3)",
     "Loss here can lower MHC-I and weaken CD8 T-cell recognition. A flag that does not need MHC-I is attractive."),
    ("15q", "loss", "B2M (15q21.1)",
     "B2M is needed for MHC-I on the cell surface. Its loss is a known route of resistance to checkpoint blockade."),
    ("6q", "loss", "ULBP genes (6q25)",
     "The ULBP genes encode NKG2D ligands, so loss lowers the cells' own ligands. A vector-supplied ligand can offset that."),
    ("9p", "loss", "type I interferon cluster with CDKN2A (9p21.3)",
     "Deletions that take the interferon cluster have been linked to immune evasion and checkpoint-blockade resistance."),
    ("9p", "gain", "CD274 (PD-L1) and JAK2 (9p24.1)",
     "Extra copies can raise PD-L1, a checkpoint signal that restrains T cells."),
    ("3q", "gain", "CD47 (3q13.12)",
     "CD47 is the don't-eat-me signal for macrophages. Extra copies may raise it."),
]


@dataclass
class Segment:
    chrom: str
    start: int
    end: int
    value: float


@dataclass
class ArmResult:
    arm: str
    coverage: float      # fraction of the arm covered by input segments
    mean_lr: float       # mean log2 ratio relative to the baseline
    frac_gain: float
    frac_loss: float
    call: str            # gain, loss, neutral, partial, no data or not scored


@dataclass
class AneuploidyResult:
    sample: str
    build: str
    value_type: str
    baseline_desc: str
    arms: List[ArmResult]
    whole_chromosomes: List[str]
    gains: int
    losses: int
    fga: float           # percent of the covered autosomal genome that is gained or lost
    wgd_suspected: bool
    findings: List[Tuple[str, str]]
    notes: List[str]

    @property
    def score(self) -> int:
        return self.gains + self.losses


@dataclass
class Recommendation:
    strategy: str
    gene: str
    promoter: Optional[str]
    title: str
    reason: str


_CHROM_NAMES = ("chrom", "chromosome", "chr", "seqnames", "seqname", "contig")
_START_NAMES = ("start", "loc.start", "begin", "startpos", "start_pos", "chromstart", "segment_start")
_END_NAMES = ("end", "loc.end", "stop", "endpos", "end_pos", "chromend", "segment_end")
_CN_NAMES = ("cn", "copy_number", "copynumber", "total_cn", "tcn", "total_copy_number", "copy.number", "tcn.em", "totalcn")
_LOG2_NAMES = ("log2", "log2ratio", "log2_ratio", "seg.mean", "segmean", "seg_mean", "logr",
               "log2_copy_ratio", "segment_mean", "mean", "value")
_SAMPLE_NAMES = ("sample", "sample_id", "sampleid", "id", "tumor_sample_barcode", "case", "barcode")


def _num(text: str) -> Optional[float]:
    try:
        return float(text)
    except ValueError:
        return None


def norm_chrom(name: str) -> Optional[str]:
    """'chr7' -> '7', '23' -> 'X'. Returns None for mitochondrial and unplaced contigs."""
    s = re.sub(r"^chr", "", name.strip(), flags=re.I).upper()
    if s in ("X", "23"):
        return "X"
    if s in ("Y", "24"):
        return "Y"
    return str(int(s)) if s.isdigit() and 1 <= int(s) <= 22 else None


def weighted_median(pairs: List[Tuple[float, float]]) -> float:
    pairs = sorted(p for p in pairs if p[1] > 0)
    total, acc = sum(w for _, w in pairs), 0.0
    for value, weight in pairs:
        acc += weight
        if acc >= total / 2:
            return value
    return 0.0


def parse_segments(text: str, value_type: str = "auto", sample: Optional[str] = None
                   ) -> Tuple[List[Segment], str, List[str], str]:
    """Read a copy-number segment file (TSV, CSV or whitespace separated).

    Needs chromosome, start, end and one value column (copy number or log2 ratio).
    Returns (segments, value type 'cn' or 'log2', notes, sample name).
    """
    value_type = pick(value_type, ["auto", "cn", "log2"], "value type")
    rows = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not rows:
        raise ValueError("The copy-number file is empty.")
    delim = "\t" if "\t" in rows[0] else ("," if "," in rows[0] else None)

    def split(line: str) -> List[str]:
        return [c.strip().strip('"') for c in (line.split(delim) if delim else line.split())]

    first = split(rows[0])
    notes: List[str] = []
    i_sample: Optional[int] = None
    if not any(_num(c) is not None for c in first):          # header row present
        names = [c.lower() for c in first]

        def find(options):
            return next((names.index(o) for o in options if o in names), None)

        i_chr, i_start, i_end = find(_CHROM_NAMES), find(_START_NAMES), find(_END_NAMES)
        if None in (i_chr, i_start, i_end):
            raise ValueError("Could not find chromosome, start and end columns. Header found: " + ", ".join(first))
        i_sample, i_cn, i_log = find(_SAMPLE_NAMES), find(_CN_NAMES), find(_LOG2_NAMES)
        vt = value_type
        if vt == "cn":
            i_val = i_cn
        elif vt == "log2":
            i_val = i_log
        else:
            vt, i_val = ("cn", i_cn) if i_cn is not None else ("log2", i_log)
        if i_val is None:
            raise ValueError("Could not find a %s column. Header found: %s"
                             % ("copy-number" if vt == "cn" else "log2 ratio" if vt == "log2" else "copy-number or log2 ratio",
                                ", ".join(first)))
        body = rows[1:]
    else:                                                     # no header: guess from the column count
        vt = value_type
        if len(first) == 4:
            i_chr, i_start, i_end, i_val = 0, 1, 2, 3
        elif len(first) == 6:                                 # DNAcopy: ID, chrom, start, end, n.marks, seg.mean
            i_sample, i_chr, i_start, i_end, i_val = 0, 1, 2, 3, 5
        else:
            raise ValueError("No header row and %d columns. Add a header (chrom, start, end, log2 or copy_number) "
                             "or use 4 columns: chrom, start, end, value." % len(first))
        body = rows

    parsed, skipped = [], 0
    for line in body:
        c = split(line)
        try:
            parsed.append((c[i_sample] if i_sample is not None else "", c[i_chr],
                           int(float(c[i_start])), int(float(c[i_end])), float(c[i_val])))
        except (ValueError, IndexError):
            skipped += 1
    if not parsed:
        raise ValueError("No usable rows were found. Check the delimiter and the column order.")
    if skipped:
        notes.append("%d row(s) could not be read and were skipped." % skipped)

    chosen = ""
    if i_sample is not None:
        names_in_file = list(dict.fromkeys(p[0] for p in parsed))
        if len(names_in_file) > 1 and not sample:
            shown = ", ".join(names_in_file[:5]) + (" ..." if len(names_in_file) > 5 else "")
            raise ValueError("The file holds %d samples (%s). Choose one with --sample." % (len(names_in_file), shown))
        chosen = sample or names_in_file[0]
        if chosen not in names_in_file:
            raise ValueError("Sample %r is not in the file." % chosen)
        parsed = [p for p in parsed if p[0] == chosen]
    if vt == "auto":
        vt = "log2" if any(p[4] < 0 for p in parsed) else "cn"
        notes.append("The value column was not named, so it was read as %s." % ("log2 ratios" if vt == "log2" else "copy numbers"))

    segments = []
    for _sample, chrom, start, end, value in parsed:
        ch = norm_chrom(chrom)
        if ch is not None and end != start:
            segments.append(Segment(ch, min(start, end), max(start, end), value))
    return segments, vt, notes, chosen


def aneuploidy_level(score: int) -> str:
    """Rough orientation bins only; compare with the distribution for your tumor type."""
    return "low" if score <= 2 else ("intermediate" if score <= 9 else "high")


def analyze_aneuploidy(segments: List[Segment], build: str = "hg38", value_type: str = "log2",
                       ploidy: Optional[float] = None, center: bool = True, gain: float = 0.2,
                       loss: float = -0.2, sample: str = "", notes: Optional[List[str]] = None) -> AneuploidyResult:
    """Call arm-level gains and losses and compute an aneuploidy score.

    An arm is called gained or lost when at least 80% of its covered length changes in one
    direction, and neutral when at most 20% does (the rule used by Taylor et al., Cancer Cell 2018).
    The score is the number of gained or lost arms among the 39 scorable autosomal arms.
    """
    build = pick(build, BUILDS, "genome build")
    value_type = pick(value_type, ["cn", "log2"], "value type")
    if not (gain > 0 > loss):
        raise ValueError("The gain threshold must be above 0 and the loss threshold below 0.")
    dims = BUILDS[build]
    notes = list(notes or [])
    by_chrom: Dict[str, List[Segment]] = {}
    for s in segments:
        if s.chrom in dims:
            a, b = max(0, s.start), min(dims[s.chrom][0], s.end)
            if b > a:
                by_chrom.setdefault(s.chrom, []).append(Segment(s.chrom, a, b, s.value))
    clipped = [s for segs in by_chrom.values() for s in segs]
    if not clipped:
        raise ValueError("No autosomal segments overlap the %s chromosomes. Check the genome build and chromosome names." % build)

    weights = [(s.value, s.end - s.start) for s in clipped]
    if value_type == "cn":
        base = float(ploidy) if ploidy else float(max(1, round(weighted_median(weights))))
        if base <= 0:
            raise ValueError("Ploidy must be above 0.")
        desc = "baseline copy number %g (%s)" % (base, "given" if ploidy else "genome-wide median")

        def to_lr(v: float) -> float:
            return math.log2(max(v, 0.1) / base)
    else:
        centre = weighted_median(weights) if center else 0.0
        desc = "log2 ratios centred on %+.3f" % centre if center else "log2 ratios, not centred"

        def to_lr(v: float) -> float:
            return v - centre
        base = 2.0
        notes.append("Log2 ratios cannot show absolute ploidy, so whole-genome doubling cannot be called."
                     + (" The --ploidy value is only used with copy-number input." if ploidy else ""))

    arms: List[ArmResult] = []
    for arm in ARM_ORDER:
        chrom, side = arm[:-1], arm[-1]
        length, cen = dims[chrom]
        lo, hi = (0, cen) if side == "p" else (cen, length)
        cov = g = l = 0
        weighted = 0.0
        for s in by_chrom.get(chrom, []):
            ov = min(s.end, hi) - max(s.start, lo)
            if ov <= 0:
                continue
            r = to_lr(s.value)
            cov += ov
            weighted += r * ov
            if r >= gain:
                g += ov
            elif r <= loss:
                l += ov
        cov_eff = min(cov, hi - lo)
        fg, fl = (min(g / cov, 1.0), min(l / cov, 1.0)) if cov else (0.0, 0.0)
        coverage = cov_eff / (hi - lo)
        if arm in ACROCENTRIC_P:
            call = "not scored"
        elif coverage < 0.5:
            call = "no data"
        elif fg >= 0.8:
            call = "gain"
        elif fl >= 0.8:
            call = "loss"
        elif max(fg, fl) <= 0.2:
            call = "neutral"
        else:
            call = "partial"
        arms.append(ArmResult(arm, coverage, weighted / cov if cov else 0.0, fg, fl, call))

    by_arm = {a.arm: a for a in arms}
    whole = []
    for chrom in AUTOSOMES:
        if chrom + "p" in ACROCENTRIC_P:      # the p arm cannot be assessed, so the q-arm call stands alone
            continue
        sides = [by_arm[chrom + s].call for s in "pq"]
        for direction in ("gain", "loss"):
            if all(c == direction for c in sides):
                whole.append("%s %s" % (chrom, direction))
    gains = sum(a.call == "gain" for a in arms)
    losses = sum(a.call == "loss" for a in arms)
    covered = sum(s.end - s.start for s in clipped)
    altered = sum(s.end - s.start for s in clipped if to_lr(s.value) >= gain or to_lr(s.value) <= loss)
    calls = {a.arm: a.call for a in arms}
    findings = [("%s %s: %s" % (arm, direction, label), note)
                for arm, direction, label, note in IMMUNE_ARM_NOTES if calls.get(arm) == direction]
    with_data = sum(a.call not in ("no data", "not scored") for a in arms)
    if with_data < 20:
        notes.append("Only %d of the %d scorable arms had enough data, so the score is a lower bound."
                     % (with_data, N_SCORED_ARMS))
    if gains + losses < 3:
        notes.append("Few arm-level events were found. Check the baseline or ploidy before relying on the "
                     "aneuploidy-directed choices below.")
    return AneuploidyResult(sample, build, value_type, desc, arms, whole, gains, losses,
                            100.0 * altered / covered, value_type == "cn" and base >= 3, findings, notes)


def recommend(res: AneuploidyResult) -> List[Recommendation]:
    """Rule-based design suggestions from the arm calls. They follow known gene locations and
    published associations; they are not predictions of response."""
    calls = {a.arm: a.call for a in res.arms}
    recs: List[Recommendation] = []

    def add(strategy: str, gene: str, title: str, reason: str) -> None:
        if not any(r.strategy == strategy and r.gene == gene for r in recs):
            recs.append(Recommendation(strategy, gene, STRATEGIES[strategy]["promoter"], title, reason))

    why = []
    if calls.get("6p") == "loss" or calls.get("15q") == "loss":
        why.append("6p or 15q is lost (HLA class I region, B2M), which points to weaker MHC-I presentation, so a flag "
                   "that does not need MHC-I is attractive.")
    if calls.get("6q") == "loss":
        why.append("6q is lost and carries the ULBP genes, so the cells probably have fewer NKG2D ligands of their own.")
    add("flag", "ULBP2", "NKG2D ligand flag",
        " ".join(why + ["NKG2D ligands mark stressed or transformed cells for attack by NK cells and CD8 T cells. Shed "
                        "soluble ligand can blunt the effect, so check that it stays on the surface in your model."]))
    add("flag", "tEGFR (surface tag)", "Surface tag for cetuximab",
        "Adds an inert surface epitope that cetuximab can target (ADCC and CDC), and that a matching CAR could target. "
        "It does not rely on MHC-I or on NKG2D.")
    if calls.get("3q") == "gain":
        add("flag", "CALR", "Calreticulin eat-me signal",
            "3q is gained and carries CD47, the don't-eat-me signal. Surface calreticulin is the opposing eat-me signal "
            "for macrophages.")
        add("guide", "CD47", "CD47 knockdown", "3q is gained, so CD47 dosage may be raised. Knocking it down exposes "
                                              "the cells to macrophages. Pair with an eat-me flag.")
    if res.score >= 10 or res.wgd_suspected:
        add("cyt", "IL-12 (single chain)", "Cytokine to warm up the tumor",
            "Tumors with many arm-level events tend to show markers of immune evasion and respond less well to "
            "immunotherapy (Davoli et al., Science 2017). Local IL-12 can help recruit and activate immune cells. "
            "Watch for systemic toxicity.")
    if calls.get("9p") == "gain":
        add("guide", "CD274", "PD-L1 knockdown", "9p is gained and carries CD274 (PD-L1) and JAK2. Knocking PD-L1 down "
                                                "can release T cells.")
    n_tsg = n_onc = 0
    for gene, (arm, role) in ARM_GENES.items():
        if role == "tsg" and calls.get(arm) == "loss" and n_tsg < 2 and GENES[gene][0] <= 3000:
            n_tsg += 1
            add("rep", gene, "Replace %s" % gene,
                "%s is lost and carries the tumor suppressor %s. Arm-level loss does not prove the gene is inactive, so "
                "confirm with mutation and expression data." % (arm, gene))
        elif role == "onc" and calls.get(arm) == "gain" and n_onc < 2:
            n_onc += 1
            add("guide", gene, "Knock down %s" % gene,
                "%s is gained and carries the oncogene %s. Confirm that %s is a dependency in your model before "
                "committing." % (arm, gene, gene))
    return recs


def basis_text(res: AneuploidyResult, rec: Recommendation) -> str:
    where = "sample %s, " % res.sample if res.sample else ""
    return ("Aneuploidy analysis (%s%s): score %d (%d gains, %d losses), %.0f%% of the genome altered.\n"
            "Chosen design: %s (%s). %s" % (where, res.build, res.score, res.gains, res.losses, res.fga,
                                          rec.title, rec.gene, rec.reason))


def config_from_recommendation(cfg: Config, rec: Recommendation, basis: str, keep_promoter: bool = False) -> Config:
    cfg.strategy, cfg.gene, cfg.basis = rec.strategy, rec.gene, basis
    if rec.promoter and not keep_promoter:
        cfg.promoter = rec.promoter
    return cfg


def aneuploidy_report(res: AneuploidyResult, recs: List[Recommendation]) -> str:
    def arms_with(call: str) -> str:
        return ", ".join(a.arm for a in res.arms if a.call == call) or "none"

    title = "Aneuploidy analysis" + (": " + res.sample if res.sample else "")
    lines = [title, "=" * len(title),
             "Genome build %s. Input: %s; %s." % (res.build, "copy number" if res.value_type == "cn" else "log2 ratios",
                                                 res.baseline_desc),
             "",
             "Aneuploidy score: %d of %d scorable arms (%d gains, %d losses). Level: %s (orientation only)."
             % (res.score, N_SCORED_ARMS, res.gains, res.losses, aneuploidy_level(res.score)),
             "Fraction of the genome gained or lost: %.0f%%." % res.fga,
             "Gained arms: " + arms_with("gain"),
             "Lost arms: " + arms_with("loss")]
    if res.whole_chromosomes:
        lines.append("Whole-chromosome events: " + ", ".join(res.whole_chromosomes))
    lines.append("Partly altered, not called (20 to 80% of the arm): " + arms_with("partial"))
    lines.append("Not enough data (under 50% of the arm covered): " + arms_with("no data"))
    if res.wgd_suspected:
        lines.append("The baseline copy number is 3 or more, which suggests whole-genome doubling.")
    lines += ["Note: " + n for n in res.notes]
    if res.findings:
        lines += ["", "Immune-relevant findings"] + ["  - %s. %s" % f for f in res.findings]
    lines += ["", "Recommended designs (use --design-from-cnv N, or the button in the window)"]
    for i, r in enumerate(recs, 1):
        lines += ["  %d. %s: %s [%s]" % (i, r.title, r.gene, STRATEGIES[r.strategy]["name"]), "     " + r.reason]
    lines += ["", "Limits",
              "  - Arm-level calls from segments are a screen. Confirm with purity- and ploidy-corrected, allele-specific data.",
              "  - The vector cannot sense aneuploidy. Selectivity comes from delivery and the promoter.",
              "  - Aneuploid cells can occur in some normal tissues. Test expression and killing in normal cells.",
              "  - Suggestions follow known gene locations and published associations. They do not predict response."]
    return "\n".join(lines) + "\n"


def arms_tsv(res: AneuploidyResult) -> str:
    rows = ["arm\tcoverage\tmean_log2_ratio\tfraction_gain\tfraction_loss\tcall"]
    rows += ["%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s" % (a.arm, a.coverage, a.mean_lr, a.frac_gain, a.frac_loss, a.call)
             for a in res.arms]
    return "\n".join(rows) + "\n"


def write_aneuploidy_files(res: AneuploidyResult, recs: List[Recommendation], outdir: str, stem: str) -> List[str]:
    os.makedirs(outdir, exist_ok=True)
    stem = re.sub(r"[^\w.-]+", "_", stem).strip("_") or "sample"
    paths = []
    for suffix, text in (("_aneuploidy.txt", aneuploidy_report(res, recs)), ("_arms.tsv", arms_tsv(res))):
        path = os.path.join(outdir, stem + suffix)
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        paths.append(path)
    return paths


# arm-level log2 ratios for the synthetic example (all other arms are neutral)
EXAMPLE_ARM_VALUES = {"1q": 0.40, "3q": 0.50, "5p": 0.45, "6p": -0.45, "7p": 0.50, "7q": 0.50, "8q": 0.60,
                      "9p": -0.70, "10q": -0.50, "13q": -0.45, "15q": -0.50, "17p": -0.80, "18q": -0.60, "20q": 0.40}


def example_cnv_text(build: str = "hg38") -> str:
    """A synthetic segment file for trying the analysis. It is not patient data."""
    lines = ["# Synthetic example for testing only. This is not patient data.", "chrom\tstart\tend\tlog2"]
    for chrom in AUTOSOMES:
        length, cen = BUILDS[build][chrom]
        for side, lo, hi in (("p", 1, cen), ("q", cen, length)):
            arm = chrom + side
            if arm in ACROCENTRIC_P:
                continue
            if arm == "4q":                       # only 60% of the arm is altered, so it stays uncalled
                cut = cen + int(0.6 * (length - cen))
                lines += ["%s\t%d\t%d\t-0.50" % (chrom, lo, cut), "%s\t%d\t%d\t0.02" % (chrom, cut, hi)]
            else:
                lines.append("%s\t%d\t%d\t%.2f" % (chrom, lo, hi, EXAMPLE_ARM_VALUES.get(arm, 0.02)))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Exporters
# --------------------------------------------------------------------------
def _ascii(text: str) -> str:
    return text.encode("ascii", "replace").decode().replace('"', "")


def to_genbank(d: Design) -> str:
    today = datetime.date.today()
    mon = "JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC".split()[today.month - 1]
    locus = re.sub(r"[^\w.-]", "_", d.name)[:16]
    out = ["LOCUS       %-16s %7d bp    DNA     circular SYN %02d-%s-%d" % (locus, d.total, today.day, mon, today.year),
           "DEFINITION  %s, designed with Cancer Gene Therapy Vector Designer." % _ascii(d.name),
           "FEATURES             Location/Qualifiers"]
    for f in d.features:
        out.append("     %-16s%d..%d" % (KINDS[f.kind][1], f.start, f.end))
        out.append('                     /label="%s"' % _ascii(f.name))
    out.append("ORIGIN")
    s = d.seq.lower()
    for i in range(0, len(s), 60):
        line = s[i:i + 60]
        out.append("%9d %s" % (i + 1, " ".join(line[j:j + 10] for j in range(0, len(line), 10))))
    out.append("//")
    return "\n".join(out) + "\n"


def to_fasta(d: Design) -> str:
    header = ">%s %d bp circular" % (re.sub(r"\s+", "_", d.name), d.total)
    return header + "\n" + "\n".join(d.seq[i:i + 70] for i in range(0, d.total, 70)) + "\n"


def to_svg(d: Design) -> str:
    T, cx, cy, r = d.total, 180.0, 180.0, 130.0
    height = max(360, 48 + 18 * len(d.features))

    def pt(a):
        return cx + r * math.sin(a), cy - r * math.cos(a)

    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 700 %d" width="700" height="%d" '
           'font-family="Helvetica, Arial, sans-serif">' % (height, height),
           '<rect width="700" height="%d" fill="#ffffff"/>' % height,
           '<circle cx="180" cy="180" r="130" fill="none" stroke="#d3dde5" stroke-width="30" opacity=".6"/>']
    for f in d.features:
        a0 = 2 * math.pi * (f.start - 1) / T
        a1 = max(2 * math.pi * f.end / T, a0 + 0.025)
        (x0, y0), (x1, y1) = pt(a0), pt(a1)
        svg.append('<path d="M%.2f %.2fA%d %d 0 %d 1 %.2f %.2f" fill="none" stroke="%s" stroke-width="26">'
                   '<title>%s: %d..%d (%d bp)</title></path>'
                   % (x0, y0, r, r, 1 if a1 - a0 > math.pi else 0, x1, y1, KINDS[f.kind][0],
                      xml_escape(f.name), f.start, f.end, f.length))
    title = d.name if len(d.name) <= 26 else d.name[:25] + "..."
    svg.append('<line x1="180" y1="30" x2="180" y2="12" stroke="#14202b" stroke-width="1.5"/>'
               '<text x="186" y="18" font-size="10" fill="#587083">1</text>')
    svg.append('<text x="180" y="172" text-anchor="middle" font-size="14" font-weight="600" fill="#14202b">%s</text>'
               % xml_escape(title))
    svg.append('<text x="180" y="194" text-anchor="middle" font-size="13" fill="#587083">%s bp</text>' % format(d.total, ","))
    for i, f in enumerate(d.features):
        y = 30 + 18 * i
        svg.append('<rect x="372" y="%d" width="11" height="11" rx="2" fill="%s"/>' % (y - 10, KINDS[f.kind][0]))
        svg.append('<text x="390" y="%d" font-size="12" fill="#14202b">%s <tspan fill="#587083">%d..%d</tspan></text>'
                   % (y, xml_escape(f.name), f.start, f.end))
    svg.append("</svg>")
    return "\n".join(svg) + "\n"


def to_summary(d: Design) -> str:
    lines = [d.name, "=" * len(d.name),
             "Vector: %s    Strategy: %s    Size: %s bp" % (d.vector, STRATEGIES[d.strategy]["name"], format(d.total, ",")),
             "", "%-8s %-14s %7s  %s" % ("id", "position", "bp", "element")]
    for f in d.features:
        note = "  (your sequence)" if f.custom or f.id in ("goi", "guide") and "N" not in f.seq else ""
        lines.append("%-8s %-14s %7d  %s%s" % (f.id, "%d..%d" % (f.start, f.end), f.length, f.name, note))
    lines += ["", "Design checks"] + ["  [%s] %s" % (lv, tx) for lv, tx in d.messages]
    if d.basis:
        lines += ["", "Design basis", d.basis]
    lines += ["", "Elements written as N have a typical length only. Supply verified sequences with",
              "--override ID=SEQUENCE (ids are in the first column) or --seq for the insert."]
    return "\n".join(lines) + "\n"


def write_outputs(d: Design, outdir: str) -> List[str]:
    os.makedirs(outdir, exist_ok=True)
    stem = re.sub(r"[^\w.-]+", "_", d.name).strip("_") or "vector"
    paths = []
    for ext, text in ((".gb", to_genbank(d)), (".fasta", to_fasta(d)), (".svg", to_svg(d)), (".txt", to_summary(d))):
        p = os.path.join(outdir, stem + ext)
        with open(p, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        paths.append(p)
    return paths


# --------------------------------------------------------------------------
# Command line and interactive mode
# --------------------------------------------------------------------------
def read_seq_arg(value: Optional[str]) -> str:
    """A path to a text/FASTA file, or the sequence itself."""
    if not value:
        return ""
    if os.path.isfile(value):
        with open(value, "r", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    return value


def print_catalog() -> None:
    print("Vectors:      " + ", ".join(VECTORS))
    print("\nStrategies and their genes:")
    for k, v in STRATEGIES.items():
        print("  %-6s %s\n         genes: %s\n         %s" % (k, v["name"], ", ".join(v["genes"]), v["tip"]))
    for label, opts in (("Promoters", PROMOTERS), ("Markers", MARKERS), ("Linkers", LINKS), ("Tags", TAGS),
                        ("PolyA", POLYA), ("Bacterial selection", RESISTANCE),
                        ("CAR hinges", HINGES), ("CAR costimulatory", COSTIM)):
        print("%s: %s" % (label, ", ".join(opts)))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="vector_designer", description="Cancer gene therapy viral vector designer.",
        epilog="Run with no arguments to open the graphical window. Use --list to see the available parts.")
    p.add_argument("--gui", action="store_true", help="open the graphical window (the default with no arguments)")
    p.add_argument("--wizard", action="store_true", help="answer text prompts instead of using the window")
    p.add_argument("--list", action="store_true", help="list strategies, genes and parts, then exit")
    p.add_argument("--vector", default="lenti", help="lenti, aav or plasmid (default: lenti)")
    p.add_argument("--strategy", default="rep", help="rep, sui, cyt, car, guide or model (default: rep)")
    p.add_argument("--gene", default="", help="gene or antigen name (library entry or any custom name)")
    p.add_argument("--variant", default="", help="variant label, for example R175H")
    p.add_argument("--seq", default="", help="insert ORF (or CAR scFv) as a file path or a sequence")
    p.add_argument("--promoter", help="promoter (default depends on strategy)")
    p.add_argument("--marker", help="reporter or selection marker (default depends on strategy)")
    p.add_argument("--link", default="P2A", help="how the marker is linked: P2A, IRES or 'separate PGK promoter'")
    p.add_argument("--tag", default="none", help="C-terminal tag")
    p.add_argument("--polya", help="polyadenylation signal: bGH, SV40 or none")
    p.add_argument("--resistance", default="AmpR", help="bacterial selection: AmpR or KanR")
    p.add_argument("--no-wpre", action="store_true", help="leave out WPRE")
    p.add_argument("--wpre", action="store_true", help="force WPRE on (default: on except for plasmid)")
    p.add_argument("--no-kozak", action="store_true", help="leave out the Kozak sequence")
    p.add_argument("--hinge", default="CD8a hinge/TM", help="CAR hinge/transmembrane domain")
    p.add_argument("--costim", default="4-1BB", help="CAR costimulatory domain: 4-1BB, CD28 or none")
    p.add_argument("--guide-type", default="shRNA", help="shRNA or sgRNA (guide strategy)")
    p.add_argument("--target", default="", help="shRNA sense strand or sgRNA protospacer (no PAM), or a file path")
    p.add_argument("--override", action="append", default=[], metavar="ID=SEQ",
                   help="replace a placeholder element with a verified sequence or file (repeatable)")
    p.add_argument("--cnv", metavar="FILE", help="copy-number segment file to analyze for aneuploidy")
    p.add_argument("--build", default="hg38", help="genome build of the segment file: hg38 or hg19 (default: hg38)")
    p.add_argument("--cnv-type", default="auto", help="value column: auto, cn (copy number) or log2")
    p.add_argument("--sample", help="sample name, when the segment file holds more than one")
    p.add_argument("--ploidy", type=float, help="baseline copy number for copy-number input (default: genome median)")
    p.add_argument("--gain-threshold", type=float, default=0.2, help="log2 ratio for a gain (default: 0.2)")
    p.add_argument("--loss-threshold", type=float, default=-0.2, help="log2 ratio for a loss (default: -0.2)")
    p.add_argument("--no-center", action="store_true", help="do not centre log2 ratios on their median")
    p.add_argument("--design-from-cnv", type=int, nargs="?", const=1, metavar="N",
                   help="after the analysis, build the vector for suggestion N (default 1)")
    p.add_argument("--example-cnv", metavar="FILE", help="write a synthetic example segment file and exit")
    p.add_argument("--name", default="", help="vector name (default: generated)")
    p.add_argument("--out", default=".", help="output folder (default: current folder)")
    p.add_argument("--strict", action="store_true", help="exit with status 1 if any design check reports an error")
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    return p


def config_from_args(a: argparse.Namespace) -> Config:
    overrides = {}
    for item in a.override:
        if "=" not in item:
            raise ValueError("--override needs ID=SEQUENCE, got %r" % item)
        key, val = item.split("=", 1)
        overrides[key.strip()] = read_seq_arg(val.strip())
    wpre = True if a.wpre else (False if a.no_wpre else None)
    return Config(vector=a.vector, strategy=a.strategy, gene=a.gene, variant=a.variant,
                  sequence=read_seq_arg(a.seq), promoter=a.promoter, marker=a.marker, link=a.link,
                  tag=a.tag, polya=a.polya, resistance=a.resistance, wpre=wpre, kozak=not a.no_kozak,
                  hinge=a.hinge, costim=a.costim, guide_type=a.guide_type,
                  target=read_seq_arg(a.target), name=a.name, overrides=overrides)


def _choose(label: str, options: List[str], default: str) -> str:
    print("\n%s" % label)
    for i, o in enumerate(options, 1):
        print("  %2d) %s%s" % (i, o, "   [default]" if o == default else ""))
    while True:
        raw = input("Choose 1-%d or type a name (Enter for default): " % len(options)).strip()
        if not raw:
            return default
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return options[int(raw) - 1]
        for o in options:
            if o.lower() == raw.lower():
                return o
        print("Not an option, try again.")


def wizard() -> Tuple[Config, str]:
    print("Cancer Gene Therapy Vector Designer %s" % __version__)
    print("Answer the prompts, or press Enter to accept a default.")
    cfg = Config()
    cfg.vector = _choose("Vector type", VECTORS, "lenti")
    keys = list(STRATEGIES)
    labels = ["%s: %s" % (k, STRATEGIES[k]["name"]) for k in keys]
    chosen = _choose("Therapeutic strategy", labels, labels[0])
    cfg.strategy = keys[labels.index(chosen)]
    S = STRATEGIES[cfg.strategy]
    print("\n" + S["tip"])
    cfg.gene = _choose("Gene", S["genes"], S["genes"][0])
    if cfg.gene == "Custom":
        cfg.gene = input("Gene name: ").strip() or "Custom"
    if cfg.strategy == "guide":
        cfg.guide_type = _choose("Guide type", ["shRNA", "sgRNA"], "shRNA")
        cfg.target = read_seq_arg(input("\nTarget sequence (or file path; Enter for a placeholder): ").strip())
    else:
        if cfg.strategy != "guide":
            cfg.variant = input("\nVariant label, for example R175H (Enter to skip): ").strip()
        cfg.sequence = read_seq_arg(input("Insert sequence or FASTA file path (Enter for a placeholder): ").strip())
        if cfg.strategy == "car":
            cfg.hinge = _choose("CAR hinge/transmembrane", list(HINGES), "CD8a hinge/TM")
            cfg.costim = _choose("Costimulatory domain", list(COSTIM), "4-1BB")
        else:
            cfg.tag = _choose("C-terminal tag", list(TAGS), "none")
    cfg.promoter = _choose("Promoter", list(PROMOTERS), S["promoter"]) if cfg.strategy != "guide" else None
    cfg.marker = _choose("Reporter / selection marker", list(MARKERS), S["marker"])
    if cfg.marker != "none" and cfg.strategy != "guide":
        cfg.link = _choose("Marker linked by", list(LINKS), "P2A")
    cfg.name = input("\nVector name (Enter for a generated name): ").strip()
    out = input("Output folder (Enter for the current folder): ").strip() or "."
    return cfg, out


# --------------------------------------------------------------------------
# Graphical user interface (tkinter)
# --------------------------------------------------------------------------
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


class AneuploidyWindow:
    """Load a copy-number segment file, review the arm calls, and send a suggestion to the designer."""

    def __init__(self, app: "VectorDesignerApp") -> None:
        self.app = app
        self.result: Optional[AneuploidyResult] = None
        self.recs: List[Recommendation] = []
        top = self.top = tk.Toplevel(app.root)
        top.title("Aneuploidy analysis")
        top.geometry("1020x780")
        top.minsize(880, 620)
        defaults = dict(path="", build="hg38", vtype="auto", ploidy="", sample="", gain="0.2", loss="-0.2")
        self.v = {k: tk.StringVar(value=val) for k, val in defaults.items()}
        self.center = tk.BooleanVar(value=True)
        top.columnconfigure(0, weight=1)
        top.rowconfigure(3, weight=2)
        top.rowconfigure(5, weight=1)

        row = ttk.Frame(top, padding=(12, 10, 12, 0))
        row.grid(row=0, column=0, sticky="ew")
        row.columnconfigure(1, weight=1)
        ttk.Label(row, text="Copy-number segment file").grid(row=0, column=0, sticky="w")
        ttk.Entry(row, textvariable=self.v["path"]).grid(row=0, column=1, sticky="ew", padx=6)
        ttk.Button(row, text="Browse...", command=self.browse).grid(row=0, column=2)
        ttk.Button(row, text="Save example file...", command=self.save_example).grid(row=0, column=3, padx=(6, 0))

        opts = ttk.Frame(top, padding=(12, 8, 12, 0))
        opts.grid(row=1, column=0, sticky="ew")
        col = 0
        for label, key, width, values in (("Genome build", "build", 7, list(BUILDS)),
                                          ("Values", "vtype", 6, ["auto", "cn", "log2"]),
                                          ("Ploidy", "ploidy", 6, None), ("Sample", "sample", 12, None),
                                          ("Gain at", "gain", 6, None), ("Loss at", "loss", 6, None)):
            ttk.Label(opts, text=label).grid(row=0, column=col, padx=(0, 4))
            if values:
                widget = ttk.Combobox(opts, textvariable=self.v[key], values=values, width=width, state="readonly")
            else:
                widget = ttk.Entry(opts, textvariable=self.v[key], width=width)
            widget.grid(row=0, column=col + 1, padx=(0, 12))
            col += 2
        ttk.Checkbutton(opts, text="Centre log2 on median", variable=self.center).grid(row=0, column=col, padx=(0, 12))
        ttk.Button(opts, text="Analyze", command=self.run).grid(row=0, column=col + 1)

        self.summary = ttk.Label(top, wraplength=980, justify="left", padding=(12, 8))
        self.summary.grid(row=2, column=0, sticky="ew")

        arms = ttk.Frame(top, padding=(12, 0))
        arms.grid(row=3, column=0, sticky="nsew")
        arms.columnconfigure(0, weight=1)
        arms.rowconfigure(0, weight=1)
        self.arm_tree = ttk.Treeview(arms, columns=("cov", "mean", "gain", "loss", "call"), height=11, selectmode="browse")
        for name, text, width in (("#0", "Arm", 80), ("cov", "Covered", 80), ("mean", "Mean log2", 90),
                                  ("gain", "Gained", 80), ("loss", "Lost", 80), ("call", "Call", 110)):
            self.arm_tree.heading(name, text=text)
            self.arm_tree.column(name, width=width, anchor="w")
        self.arm_tree.tag_configure("gain", foreground="#0b5cad")
        self.arm_tree.tag_configure("loss", foreground="#b3261e")
        asb = ttk.Scrollbar(arms, orient="vertical", command=self.arm_tree.yview)
        self.arm_tree.configure(yscrollcommand=asb.set)
        self.arm_tree.grid(row=0, column=0, sticky="nsew")
        asb.grid(row=0, column=1, sticky="ns")

        ttk.Label(top, text="Suggested designs", font=self.app.bold, padding=(12, 8, 12, 2)).grid(row=4, column=0, sticky="w")
        recs = ttk.Frame(top, padding=(12, 0))
        recs.grid(row=5, column=0, sticky="nsew")
        recs.columnconfigure(0, weight=1)
        recs.rowconfigure(0, weight=1)
        self.rec_tree = ttk.Treeview(recs, columns=("kind",), height=7, selectmode="browse")
        self.rec_tree.heading("#0", text="Design")
        self.rec_tree.heading("kind", text="Strategy")
        self.rec_tree.column("#0", width=480, anchor="w")
        self.rec_tree.column("kind", width=300, anchor="w")
        self.rec_tree.bind("<<TreeviewSelect>>", self.on_select)
        rsb = ttk.Scrollbar(recs, orient="vertical", command=self.rec_tree.yview)
        self.rec_tree.configure(yscrollcommand=rsb.set)
        self.rec_tree.grid(row=0, column=0, sticky="nsew")
        rsb.grid(row=0, column=1, sticky="ns")

        self.reason = ttk.Label(top, wraplength=980, justify="left", padding=(12, 6), foreground="#587083")
        self.reason.grid(row=6, column=0, sticky="ew")
        bar = ttk.Frame(top, padding=(12, 0, 12, 10))
        bar.grid(row=7, column=0, sticky="ew")
        ttk.Button(bar, text="Use selected in vector designer", command=self.use_selected).pack(side="left")
        ttk.Button(bar, text="Save report...", command=self.save_report).pack(side="left", padx=6)
        ttk.Label(bar, foreground="#587083", text="The vector cannot sense aneuploidy. Selectivity comes from "
                  "delivery and the promoter.").pack(side="right")

    def browse(self) -> None:
        path = filedialog.askopenfilename(title="Open a copy-number segment file", filetypes=[
            ("Segment files", "*.seg *.tsv *.txt *.csv"), ("All files", "*.*")])
        if path:
            self.v["path"].set(path)

    def save_example(self) -> None:
        path = filedialog.asksaveasfilename(title="Save the synthetic example file", defaultextension=".tsv",
                                            initialfile="example_segments.tsv")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(example_cnv_text(self.v["build"].get() if self.v["build"].get() in BUILDS else "hg38"))
        except OSError as exc:
            messagebox.showerror("Could not write the file", str(exc))
            return
        self.v["path"].set(path)

    def run(self) -> None:
        path = self.v["path"].get().strip()
        if not path:
            messagebox.showinfo("Choose a file", "Choose a copy-number segment file first, or save the example file "
                                "to try the analysis.")
            return
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                text = fh.read()
            ploidy = float(self.v["ploidy"].get()) if self.v["ploidy"].get().strip() else None
            gain, loss = float(self.v["gain"].get()), float(self.v["loss"].get())
            sample = self.v["sample"].get().strip() or None
            segments, vtype, notes, sample_name = parse_segments(text, self.v["vtype"].get(), sample)
            result = analyze_aneuploidy(segments, self.v["build"].get(), vtype, ploidy, self.center.get(),
                                        gain, loss, sample_name, notes)
        except (ValueError, OSError) as exc:
            messagebox.showerror("Could not analyze the file", str(exc))
            return
        self.result, self.recs = result, recommend(result)
        self.show()

    def show(self) -> None:
        res = self.result
        whole = "; whole chromosomes: " + ", ".join(res.whole_chromosomes) if res.whole_chromosomes else ""
        self.summary.configure(text="Aneuploidy score %d of %d (%d gains, %d losses; %s, for orientation only). "
                                    "%.0f%% of the genome is gained or lost%s. %s" % (
                                        res.score, N_SCORED_ARMS, res.gains, res.losses, aneuploidy_level(res.score),
                                        res.fga, whole, " ".join(res.notes)))
        self.arm_tree.delete(*self.arm_tree.get_children())
        for a in res.arms:
            self.arm_tree.insert("", "end", text=a.arm, tags=(a.call.replace(" ", "_"),),
                                 values=("%.0f%%" % (a.coverage * 100), "%+.2f" % a.mean_lr,
                                         "%.0f%%" % (a.frac_gain * 100), "%.0f%%" % (a.frac_loss * 100), a.call))
        self.rec_tree.delete(*self.rec_tree.get_children())
        for i, r in enumerate(self.recs):
            self.rec_tree.insert("", "end", iid=str(i), text="%d. %s: %s" % (i + 1, r.title, r.gene),
                                 values=(STRATEGIES[r.strategy]["name"],))
        self.reason.configure(text="Select a suggestion to read why it was made.")

    def on_select(self, _event=None) -> None:
        selected = self.rec_tree.selection()
        if selected:
            self.reason.configure(text=self.recs[int(selected[0])].reason)

    def use_selected(self) -> None:
        selected = self.rec_tree.selection()
        if self.result is None or not selected:
            messagebox.showinfo("Choose a suggestion", "Run the analysis and select a suggestion first.")
            return
        rec = self.recs[int(selected[0])]
        self.app.apply_recommendation(rec, basis_text(self.result, rec))
        self.top.destroy()

    def save_report(self) -> None:
        if self.result is None:
            return
        path = filedialog.asksaveasfilename(title="Save the report", defaultextension=".txt",
                                            initialfile="aneuploidy_report.txt")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(aneuploidy_report(self.result, self.recs))
        except OSError as exc:
            messagebox.showerror("Could not write the file", str(exc))


def launch_gui() -> int:
    """Open the window. Returns 0 after it closes, or 2 if no window could be opened."""
    if tk is None:
        print("Error: tkinter is not installed. On Debian/Ubuntu run: sudo apt install python3-tk",
              file=sys.stderr)
        return 2
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        print("Error: could not open a window (%s)." % exc, file=sys.stderr)
        return 2
    VectorDesignerApp(root)
    root.mainloop()
    return 0


def report(d: Design, paths: List[str]) -> None:
    print("\n%s  |  %s  |  %s bp" % (d.name, d.vector, format(d.total, ",")))
    print("-" * 60)
    for f in d.features:
        print("%-10s %5d..%-6d %6d bp  %s" % (f.id, f.start, f.end, f.length, f.name))
    print("\nDesign checks")
    for level, text in d.messages:
        print("  [%s] %s" % (level, text))
    print("\nFiles written:")
    for p in paths:
        print("  " + p)


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    no_arguments = argv is None and len(sys.argv) == 1
    args = parser.parse_args(argv)
    if args.list:
        print_catalog()
        return 0
    if args.example_cnv:
        try:
            with open(args.example_cnv, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(example_cnv_text(args.build if args.build in BUILDS else "hg38"))
        except OSError as exc:
            print("Error: %s" % exc, file=sys.stderr)
            return 2
        print("Wrote a synthetic example segment file to %s. Try: python vector_designer.py --cnv %s"
              % (args.example_cnv, args.example_cnv))
        return 0
    cnv_rec: Optional[Recommendation] = None
    cnv_basis = ""
    if args.cnv:
        try:
            with open(args.cnv, "r", encoding="utf-8", errors="replace") as fh:
                segments, vtype, notes, sample_name = parse_segments(fh.read(), args.cnv_type, args.sample)
            result = analyze_aneuploidy(segments, args.build, vtype, args.ploidy, not args.no_center,
                                        args.gain_threshold, args.loss_threshold, sample_name, notes)
            recs = recommend(result)
            if args.design_from_cnv is not None and not 1 <= args.design_from_cnv <= len(recs):
                raise ValueError("--design-from-cnv must be between 1 and %d." % len(recs))
            stem = sample_name or os.path.splitext(os.path.basename(args.cnv))[0]
            files = write_aneuploidy_files(result, recs, args.out, stem)
        except (ValueError, OSError) as exc:
            print("Error: %s" % exc, file=sys.stderr)
            return 2
        print(aneuploidy_report(result, recs))
        print("Files written:\n  " + "\n  ".join(files))
        if args.design_from_cnv is None:
            return 0
        cnv_rec = recs[args.design_from_cnv - 1]
        cnv_basis = basis_text(result, cnv_rec)
        print("\nBuilding the vector for suggestion %d: %s (%s)\n" % (args.design_from_cnv, cnv_rec.title, cnv_rec.gene))
    elif args.design_from_cnv is not None:
        print("Error: --design-from-cnv needs --cnv FILE.", file=sys.stderr)
        return 2
    wizard_mode = args.wizard
    if not wizard_mode and (args.gui or no_arguments):
        status = launch_gui()
        if status == 0 or args.gui:
            return status
        print("Falling back to text prompts.\n", file=sys.stderr)
        wizard_mode = True
    try:
        if wizard_mode:
            cfg, out = wizard()
        else:
            cfg, out = config_from_args(args), args.out
            if cnv_rec is not None:
                cfg = config_from_recommendation(cfg, cnv_rec, cnv_basis, keep_promoter=bool(args.promoter))
        design = build_design(cfg)
        paths = write_outputs(design, out)
    except (ValueError, OSError) as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 2
    except (KeyboardInterrupt, EOFError):
        print("\nCancelled.", file=sys.stderr)
        return 130
    report(design, paths)
    if args.strict and any(level == "error" for level, _ in design.messages):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
