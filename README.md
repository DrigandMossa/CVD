# Cancer Gene Therapy Vector Designer (`cvd`)

A research tool for designing viral vectors around cancer-relevant inserts. You choose a
delivery vector and a therapeutic strategy; the tool assembles the plasmid element by
element from a parts library, checks your insert sequence, and exports a GenBank file, a
FASTA file, a circular SVG map and a text summary.

It also reads copy-number segment files, calls arm-level gains and losses, computes an
aneuploidy score, and suggests designs based on what it finds — including immune-flagging
payloads that do not depend on MHC-I.

For any insert it can also generate **ranked, prebuilt vector options** (lentivirus and AAV,
several promoter layouts). Each option is a complete design, modified where the part library
allows, and every score, comparison and modification is explained.

Python 3.8+, standard library only, no runtime dependencies. Run it straight from a
checkout with `python -m cvd`, or `pip install .` to get a `cvd` command on your PATH.

> **Research use.** Library elements use typical lengths and are written as runs of `N`
> until you supply verified sequences. Producing virus requires institutional biosafety
> approval; a clinical gene therapy product requires GMP manufacturing, preclinical safety
> data and regulatory authorization.

---

## Contents

- [Quick start](#quick-start)
- [What it builds](#what-it-builds)
- [Running it](#running-it)
  - [Graphical window](#graphical-window)
  - [Text wizard](#text-wizard)
  - [One-shot command line](#one-shot-command-line)
- [Ranked vector options](#ranked-vector-options)
- [Aneuploidy analysis](#aneuploidy-analysis)
- [Supplying real sequences](#supplying-real-sequences)
- [Outputs](#outputs)
- [Parts library](#parts-library)
- [Design checks](#design-checks)
- [Option reference](#option-reference)
- [Exit codes](#exit-codes)
- [Project layout](#project-layout)
- [Development](#development)
- [Limitations](#limitations)

---

## Quick start

```bash
# straight from a checkout -- nothing to install
python -m cvd --list

# or install the console script
pip install .
cvd --list

# open the graphical window
python -m cvd

# see every strategy, gene and part
python -m cvd --list

# p53 replacement in AAV, ORF read from a FASTA file
python -m cvd --vector aav --strategy rep --gene TP53 \
    --seq tp53_orf.fasta --promoter EF1a --out designs

# rank the vector options for an insert and see why each ranks where it does
python -m cvd --strategy rep --gene TP53 --setting "in vivo" --rank --explain 1

# try the aneuploidy analysis on synthetic data
python -m cvd --example-cnv demo.tsv
python -m cvd --cnv demo.tsv
```

The window needs `tkinter`, which ships with the python.org installers for Windows and
macOS. On Debian/Ubuntu: `sudo apt install python3-tk`. Without `tkinter` the tool falls
back to the text prompts.

---

## What it builds

A design is a linear list of annotated elements laid out around a circular plasmid. The
backbone depends on the vector; the cassette depends on the strategy.

| Vector | Backbone elements | Cassette capacity |
|---|---|---|
| `lenti` | 5' LTR, Psi, RRE, cPPT/CTS, 3' dU3 LTR (self-inactivating), f1 ori | warns above 8 kb — large inserts lower titer |
| `aav` | 5' ITR, 3' ITR, f1 ori | warns at 4.4 kb, errors above 4.7 kb |
| `plasmid` | no viral elements | unlimited |

All three get a bacterial resistance marker and a pUC origin.

| Strategy | Name | Default promoter | Typical use |
|---|---|---|---|
| `rep` | Tumor suppressor replacement | EF1a | restore a lost suppressor with the wild-type ORF |
| `sui` | Suicide gene and prodrug | Survivin | HSV-TK + ganciclovir, iCasp9 + rimiducid, CD::UPRT + 5-FC |
| `cyt` | Immunostimulatory cytokine | EF1a | local IL-12, IL-15, IL-2, GM-CSF |
| `flag` | Immune flagging of aneuploid cells | Survivin | NKG2D ligands, calreticulin, cetuximab surface tag, T-cell engager |
| `car` | CAR T-cell therapy | MND | signal peptide – scFv – hinge/TM – costim – CD3ζ |
| `guide` | Oncogene knockdown or editing | U6 (fixed) | shRNA hairpin or sgRNA + scaffold |
| `model` | Cancer modeling (research only) | EF1a | oncogene expression for transformation assays |

---

## Running it

### Graphical window

```bash
python -m cvd            # default with no arguments
python -m cvd --gui      # force it; fails rather than falling back
```

Options on the left, live output on the right: a circular map, a table of every element
with its coordinates, and the design checks. Everything redraws about 200 ms after your
last keystroke.

- Hover an arc on the map to see the element's name, position and length.
- Select a table row to highlight the matching arc.
- Double-click a placeholder row to paste a verified sequence for that element.
- The **Vector options** tab ranks prebuilt options for the current insert; see
  [Ranked vector options](#ranked-vector-options).
- **Analysis → Aneuploidy analysis…** opens the CNV window; a suggestion there can be
  pushed straight into the designer, and the reason travels with it into the summary file.
- **File → Export all files…**, or **Copy GenBank** / **Copy FASTA** for the clipboard.

### Text wizard

```bash
python -m cvd --wizard
```

Prompts for each choice in turn. Press Enter to accept the default shown in brackets; type
either the number or the name.

### One-shot command line

```bash
# CAR T construct against CD19 in a lentivirus
python -m cvd --vector lenti --strategy car --gene CD19 \
    --seq cd19_scfv.txt --hinge "CD8a hinge/TM" --costim 4-1BB

# shRNA knockdown of KRAS with a puromycin marker
python -m cvd --strategy guide --gene KRAS --guide-type shRNA \
    --target GTTGGAGCTGATGGCGTAGG --marker Puro

# fail the build if any check reports an error
python -m cvd --strategy rep --gene TP53 --seq tp53.fasta --strict
```

---

## Ranked vector options

```bash
python -m cvd --strategy rep --gene TP53 --seq tp53.fasta --setting "in vivo" --rank
python -m cvd ... --rank --explain 1          # why option 1, and how it compares
python -m cvd ... --rank --explain all
python -m cvd ... --rank --use-option 2       # build option 2 (GenBank, FASTA, map, summary)
python -m cvd ... --rank --export-options     # every prebuilt option in its own folder
python -m cvd ... --rank --no-optimize        # rank the plain layouts, without modifications
```

**Candidates.** Lentivirus and AAV, each with a strong constitutive promoter (EF1a + WPRE),
a compact promoter without WPRE (EFS) and a tumor-selective promoter (Survivin + WPRE).
CAR constructs add a T-cell promoter (MND); guides use their single U6 cassette.

**Score.** Five criteria out of 100: packaging capacity (30), therapeutic setting (25),
expression pattern (20), insert versus vector genome (15) and vector safety profile (10).
An option whose cassette exceeds its vector's capacity scores 0 and ranks after every
option that fits, closest fit first. The delivery setting comes from `--setting`
(`in vivo` or `ex vivo`); CAR and modeling constructs are assumed ex vivo. Scores rank the
generated options against each other using fixed rules — they are not predictions of titer
or efficacy.

**Explanations.** `--explain N` (and the window's click-through) shows why the option was
suggested, what holds it back, the per-criterion breakdown, and one line for every other
option saying why this one ranks above or below it. Those lines are generated from the
actual score differences, so they cannot drift from the ranking.

**Modifications.** Every option is a complete prebuilt vector, and the tool modifies it
where the part library allows. Each change is rebuilt and measured on the design, and each
is listed with its reason, its effect on the cassette size and any cost; changes that were
considered and not made are listed too.

| Modification | When |
|---|---|
| Kozak sequence added | it was switched off |
| Marker linked by P2A instead of IRES or a second promoter | a marker is kept |
| SV40 late polyA instead of bGH | AAV (about 100 bp smaller, equal or higher expression) |
| C-terminal tag removed | therapeutic strategies |
| Reporter or selection marker removed | therapeutic strategies in vivo, or any option that does not fit |
| Short WPRE3 (247 bp) instead of WPRE (589 bp) | AAV above 75% of capacity |
| WPRE removed | only when that is what makes the option fit |
| Kanamycin instead of ampicillin selection | therapeutic strategies |

The AAV figures come from one published comparison in neurons (Choi et al. 2014, *Mol
Brain* 7:17): the short WPRE kept about 83% of the full element's expression, removing WPRE
entirely gave about 29%, and a 49 bp synthetic polyA roughly halved expression, so it is
listed as considered and not used. Other tissues may differ. Cryptic polyA and splice sites
in a lentiviral insert are reported but not fixed, because that means editing your sequence.

In the window, the **Vector options** tab shows a small map for each option; clicking a map
opens the full prebuilt vector with every element. The breakdown has two tabs: *Why this
ranking* and *Modifications*. **Use this option in the designer** carries the modifications
into the designer and records them in the design basis of the summary file.

---

## Aneuploidy analysis

```bash
python -m cvd --cnv segments.tsv
python -m cvd --cnv segments.tsv --design-from-cnv 1 --out designs
```

**Input.** TSV, CSV or whitespace-separated. Needs chromosome, start, end and one value
column — either copy number or a log2 ratio. Column names are matched against the usual
aliases (`chrom`/`chromosome`/`seqnames`, `loc.start`, `seg.mean`, `log2`, `copy_number`,
`tcn`, …), so DNAcopy output and most CNV caller formats work unchanged. A headerless file
is accepted with 4 columns (`chrom start end value`) or 6 columns (DNAcopy's
`ID chrom start end n.marks seg.mean`). Comment lines starting with `#` are skipped. If the
file holds several samples, pick one with `--sample`.

**Calling.** Log2 ratios are centred on their length-weighted genome median by default
(`--no-center` to disable); copy numbers are converted against a baseline that defaults to
the genome-wide median and can be set with `--ploidy`. An arm is called **gain** or **loss**
when at least 80% of its covered length moves in one direction, and **neutral** when at
most 20% does — the rule from Taylor et al., *Cancer Cell* 2018. Between those it is
reported as **partial** and not counted. Arms covered under 50% are **no data**; the five
acrocentric p arms (13p, 14p, 15p, 21p, 22p) are **not scored**, leaving 39 scorable
autosomal arms.

**Score.** The number of gained plus lost arms out of 39. Bins are for orientation only —
low (≤2), intermediate (3–9), high (≥10). Compare against the distribution for your tumor
type rather than treating the bin as meaningful on its own.

**Suggestions.** Rule-based, from known gene locations and published associations:

- 6p or 15q loss (HLA class I region, B2M) → a flag that does not need MHC-I
- 6q loss (ULBP cluster) → vector-supplied NKG2D ligand to replace the cells' own
- 3q gain (CD47) → calreticulin eat-me flag, and CD47 knockdown
- 9p gain (CD274, JAK2) → PD-L1 knockdown
- score ≥10 or suspected whole-genome doubling → local IL-12
- arm loss over a tumor suppressor ≤3 kb → replacement; arm gain over an oncogene → knockdown

These are not predictions of response.

---

## Supplying real sequences

Every element that is not built in starts as a run of `N` at a typical length. Three ways
to replace them:

```bash
# the insert ORF or CAR scFv — a file path or the sequence itself
--seq tp53_orf.fasta

# the shRNA sense strand or sgRNA protospacer (no PAM)
--target GTTGGAGCTGATGGCGTAGG

# any other element, by the id shown in the first column of the summary (repeatable)
--override wpre=wpre_verified.txt --override prom=ACGT...
```

FASTA headers, whitespace, digits and non-nucleotide characters are stripped; `U` becomes
`T`. The ids `goi`, `guide`, `kozak`, `stop`, `tag`, `p2a` and `term` take their sequence
from `--seq`, `--target` or the built-in library and cannot be overridden. In the GUI,
double-click a table row instead.

Addgene reference plasmids are a reasonable source for backbone elements; NCBI or CCDS for
ORFs.

---

## Outputs

Written to `--out` (default: current folder), named from the design name with unsafe
characters replaced:

| File | Contents |
|---|---|
| `<name>.gb` | GenBank with one labelled feature per element, circular topology |
| `<name>.fasta` | the full assembled sequence, 70 bp per line |
| `<name>.svg` | circular map with a legend listing every element and its coordinates |
| `<name>.txt` | element table, design checks, and the design basis if one came from a CNV run |

An aneuploidy run also writes `<sample>_aneuploidy.txt` (the full report) and
`<sample>_arms.tsv` (per-arm coverage, mean log2 ratio, gain/loss fractions and call).

Names are generated when you do not pass `--name`, following
`p{LV,AAV,Ex}-<promoter>><lead><gene><tag><marker>` — for example
`pLV-MND>CAR-CD19-P2A-EGFP` or `pAAV-EF1a>TP53`.

---

## Parts library

| Category | Options |
|---|---|
| Promoters | CMV, EF1a, EFS, CAG, PGK, UbC, hSyn, SV40, MND, **hTERT**, **Survivin** |
| Markers | none, EGFP, mCherry, Puro, Neo, Blast, Hygro |
| Marker linkage | P2A, IRES (EMCV), separate PGK promoter |
| C-terminal tags | none, FLAG3, HA, Myc, His6 |
| PolyA | bGH, SV40, none |
| Bacterial selection | AmpR, KanR |
| CAR hinge/TM | CD8a hinge/TM, CD28 hinge/TM |
| CAR costimulatory | 4-1BB, CD28, none |
| Restriction sites checked | EcoRI, BamHI, HindIII, XhoI, NotI, NheI, AgeI, XbaI, KpnI, MluI, SalI, PacI |

**hTERT** and **Survivin** are tumor-selective — expression is limited to cells where they
are active, and is usually weaker than EF1a or CMV.

Genes are grouped per strategy but any name is accepted; `Custom` or an unrecognized name
gives a 1000 bp placeholder. Several library entries carry a note that appears in the
checks (APC and NF1 are flagged as too large for AAV, for instance).

---

## Design checks

Messages are sorted error → warning → note.

**Insert sequence** — length not a multiple of 3; in-frame stop codons, with the first
codon position; missing start ATG; a trailing stop codon that will break a tag or P2A
fusion; length and GC content; cryptic polyadenylation signals (`AATAAA`, `ATTAAA`) in
lentiviral inserts; which of the listed enzymes cut and which do not.

**Guide sequences** — protospacer length (20 nt for sgRNA, 19–25 nt for shRNA); `TTTT`,
which terminates Pol III early (reported as an error); GC outside 30–55%; a reminder that
U6 transcribes best from a 5' G, and that an sgRNA vector needs Cas9 supplied separately.

**Capacity** — AAV cassette measured between the ITRs against the ~4.7 kb packaging limit;
lentiviral insert region measured between the LTRs.

**Strategy notes** — frame warnings for P2A and CAR domain fusions, tumor-selective
promoter caveats, the reminder that a `flag` vector cannot itself sense aneuploidy, and a
research-use-only note for `model`.

---

## Option reference

```
--gui                     open the graphical window (default with no arguments)
--wizard                  answer text prompts instead
--list                    list strategies, genes and parts, then exit

--vector VECTOR           lenti, aav or plasmid (default: lenti)
--strategy STRATEGY       rep, sui, cyt, flag, car, guide or model (default: rep)
--gene GENE               gene or antigen name (library entry or any custom name)
--variant LABEL           variant label, for example R175H
--seq SEQ                 insert ORF or CAR scFv — file path or sequence
--name NAME               vector name (default: generated)
--out DIR                 output folder (default: current folder)

--promoter P              promoter (default depends on strategy)
--marker M                reporter or selection marker (default depends on strategy)
--link L                  P2A, IRES or "separate PGK promoter"
--tag T                   C-terminal tag
--polya P                 bGH, SV40 or none
--resistance R            AmpR or KanR
--wpre / --no-wpre        force WPRE on or off (default: on except for plasmid)
--no-kozak                leave out the Kozak sequence
--short-wpre              use the 247 bp short WPRE (saves space in AAV)

--hinge H                 CAR hinge/transmembrane domain
--costim C                CAR costimulatory domain
--guide-type T            shRNA or sgRNA
--target SEQ              shRNA sense strand or sgRNA protospacer, or a file path
--override ID=SEQ         replace a placeholder element (repeatable)

--cnv FILE                copy-number segment file to analyze
--build B                 hg38 or hg19 (default: hg38)
--cnv-type T              auto, cn or log2
--sample NAME             sample name, when the file holds more than one
--ploidy N                baseline copy number for copy-number input
--gain-threshold N        log2 ratio for a gain (default: 0.2)
--loss-threshold N        log2 ratio for a loss (default: -0.2)
--no-center               do not centre log2 ratios on their median
--design-from-cnv [N]     build the vector for suggestion N (default 1)
--example-cnv FILE        write a synthetic example segment file and exit

--rank                    rank vector options for the insert; write an options report
--explain N|all           why option N (or all) was suggested and how it compares; implies --rank
--use-option N            build ranked option N; implies --rank
--setting S               delivery setting for ranking: "in vivo" or "ex vivo"
--no-optimize             rank the plain layouts, without modifications
--export-options          with --rank, write every prebuilt option to its own folder

--strict                  exit 1 if any design check reports an error
--version
```

---

## Exit codes

| Code | Meaning |
|---|---|
| 0 | success |
| 1 | `--strict` and at least one check reported an error |
| 2 | bad input (including a bad `--setting`, `--explain` or `--use-option`), unwritable output, or no window could be opened under `--gui` |
| 130 | cancelled (Ctrl-C or EOF in the wizard) |

---

## Project layout

```
cvd/
  __main__.py        python -m cvd
  version.py         single source of the version string
  seqops.py          sequence cleaning, revcomp, GC, case-insensitive option matching
  model.py           Config, Feature, Design, Segment, ArmResult, AneuploidyResult, ...
  parts/             pure data -- no imports from the rest of the package
    elements.py        promoters, markers, tags, linkers, polyA, element kinds
    genes.py           gene library: typical ORF length and notes
    strategies.py      the seven strategies, their genes and defaults
    enzymes.py         restriction sites screened in inserts
    genome.py          chromosome/centromere coordinates, cancer genes by arm
  design/
    builder.py         build_design: element layout, overrides, naming
    checks.py          insert-level checks (frame, stops, polyA, enzymes)
  cnv/
    segments.py        segment-file parsing and column aliasing
    analysis.py        arm calling and the aneuploidy score
    advice.py          rule-based design suggestions
    report.py          text report, arm TSV, synthetic example file
  export/
    genbank.py  fasta.py  svgmap.py  summary.py   pure Design -> str
    files.py           writes every format to a folder
  ranking/
    tables.py          capacities, criteria weights, layouts, rule text (pure data)
    model.py           InsertProfile, CriterionScore, Modification, VectorOption, RankingResult
    scoring.py         insert profile and the per-criterion points
    optimize.py        the modifications, each rebuilt and measured
    engine.py          rank_vector_options: candidates, modify, score, order
    explain.py         tables, explanations, reports, option-to-config
  cli/
    args.py            argument parsing and the parts catalogue
    wizard.py          interactive text prompts
    rank.py            --rank, --explain, --use-option, --export-options
    run.py             mode dispatch and main()
  gui/
    __init__.py        launch_gui, the only place that knows tkinter may be missing
    app.py             main window
    aneuploidy.py      CNV window
tests/               unittest suite, standard library only
```

Dependencies point one way: `parts` imports nothing from the package, `seqops` and `model`
import only `parts`, `design`/`cnv`/`export` sit on those, `ranking` sits on `design`, and
`cli` and `gui` sit on top.
`gui` is the only package that imports tkinter, and only `gui/__init__.py` tolerates its
absence -- so every command-line mode works on a machine with no Tk at all.

This README is maintained alongside those modules -- update the option reference and the
parts tables when they change.

### Known internal debt

Carried over deliberately rather than changed in the same pass as the split:

- `build_design()` still branches on `is_guide` / `is_car` rather than dispatching to a
  per-strategy cassette builder, so adding a strategy means editing branch logic.
- Planning the feature list, assigning coordinates and filling in placeholder sequence all
  happen in one pass, which is why the GUI builds a design twice per refresh (once just to
  learn which element ids exist).
- `KINDS` carries both an SVG colour and a GenBank feature key, so two exporters' concerns
  share one table in the parts library.
- Design-check messages are `(level, text)` tuples ordered by index into a literal tuple,
  with the level-to-label mapping repeated in the GUI, the terminal report and the summary.

### Building a stand-alone executable

```bash
pip install -r requirements-dev.txt
pyinstaller cvd.spec
```

The result lands in `dist/cvd` (`dist/cvd.exe` on Windows). The spec builds a console
executable on purpose: a windowed build opens the designer fine but silently discards the
output of `--list`, `--cnv` and one-shot runs, which makes the command-line modes useless.
The cost is a console window behind the designer.

`.github/workflows/release.yml` builds the same spec for Windows, Linux and macOS on every
`v*` tag (or on demand), smoke-tests each binary, uploads them as workflow artifacts and
attaches them to the release.

## Development

```bash
python -m unittest discover          # the whole suite, no dependencies
python -m unittest discover -v       # with test names
python -m unittest tests.test_design # one module
pytest -q                            # same tests, nicer output (pip install pytest)
```

The suite uses `unittest` from the standard library so it runs against a bare checkout.
Tests that need a real tkinter skip themselves when it is missing. The Vector options tab is
also tested through stand-in widgets (`tests/test_gui_options.py`), so that logic is covered
on machines with no display and no Tk.

`.github/workflows/ci.yml` runs the suite on Python 3.8, 3.11 and 3.13 across Linux,
Windows and macOS, then smoke-tests every command-line mode end to end (including that
`--strict` still exits non-zero on an oversized AAV cassette) and checks the installed
`cvd` console script. flake8 runs for real errors (undefined names, syntax) and reports
style separately without failing the build.

When adding a part to the library, add it to the table in `cvd/parts/` and to the
corresponding table in this README; the tests assert on behaviour, not on catalogue
contents, so a new promoter will not break them.

---

## Limitations

- **Lengths are typical, not real.** Any element shown as `N` needs a verified sequence
  before cloning. The tool checks frames and sizes, not biology.
- **Arm-level CNV calls are a screen.** Confirm with purity- and ploidy-corrected,
  allele-specific data. Arm loss does not prove a tumor suppressor is inactive, and arm
  gain does not prove an oncogene is a dependency.
- **Log2 input cannot show absolute ploidy**, so whole-genome doubling is not called from
  it; `--ploidy` applies only to copy-number input.
- **The vector cannot sense aneuploidy.** Selectivity comes from delivery and the
  tumor-selective promoter, not from the payload. Aneuploid cells occur in some normal
  tissues — test expression and immune killing in normal cells before relying on it.
- **Shed ligands blunt flagging.** MICA, MICB and ULBP2 can be cleaved from the cell
  surface by ADAM10/ADAM17, and soluble ligand can desensitize NKG2D. Check surface
  retention in your model.
- **Ranking scores are rule-based.** They compare the options generated here against each
  other; they do not predict titer or efficacy. Only lentivirus and AAV are ranked —
  adenovirus, dual-AAV and self-complementary AAV are not modeled, and element overrides
  are ignored when ranking.
- **Modification figures are from one study.** The short-WPRE and polyA numbers come from
  AAV cassettes in neurons; check them in your tissue.
- **No codon optimization, no primer design, no cloning strategy.** Export the GenBank
  file and take it to a dedicated tool for those.
