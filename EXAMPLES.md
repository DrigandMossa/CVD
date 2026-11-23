# Tested examples

Every command in this file is run by the test suite (`tests/test_examples_doc.py`), so the
commands and the results they promise stay correct. Outputs are abridged. All designs are
computer plans: elements shown as `N` are only the right length, and nothing here has been tested
in cells, animals or people.

Run any example from the folder that holds `cvd/`. Files go to the `--out` folder (the current
folder by default).

- [1. A tumor suppressor in AAV](#1-a-tumor-suppressor-in-aav)
- [2. A CAR T-cell vector](#2-a-car-t-cell-vector)
- [3. A knockdown vector for KRAS](#3-a-knockdown-vector-for-kras)
- [4. An insert that is too big](#4-an-insert-that-is-too-big)
- [5. Ranking vector options for a suicide gene](#5-ranking-vector-options-for-a-suicide-gene)
- [6. When no AAV option fits](#6-when-no-aav-option-fits)
- [7. Seeing the modifications](#7-seeing-the-modifications)
- [8. Building a ranked option](#8-building-a-ranked-option)
- [9. Aneuploidy on the built-in practice data](#9-aneuploidy-on-the-built-in-practice-data)
- [10. A real tumor from open data](#10-a-real-tumor-from-open-data)

---

## 1. A tumor suppressor in AAV

Put back a working copy of TP53 with a strong promoter.

<!-- test: exit=0 | expect=pAAV-EF1a>TP53 | expect=5,249 bp | expect=Design checks -->
```bash
python -m cvd --vector aav --strategy rep --gene TP53 --promoter EF1a
```

```text
pAAV-EF1a>TP53  |  aav  |  5,249 bp
itr5           1..145       145 bp  5' ITR
prom         146..1323     1178 bp  EF1a promoter
...
```

It writes a GenBank file, a FASTA file, a circular map (SVG) and a text summary.

## 2. A CAR T-cell vector

A chimeric antigen receptor against CD19, with an EGFP tracking marker, in a lentivirus.

<!-- test: exit=0 | expect=pLV-MND>CAR-CD19-P2A-EGFP -->
```bash
python -m cvd --vector lenti --strategy car --gene CD19 --marker EGFP
```

The promoter defaults to MND, which works well in T cells. The scFv is a placeholder until you
supply the real sequence with `--seq`.

## 3. A knockdown vector for KRAS

An shRNA under the U6 promoter with a puromycin marker.

<!-- test: exit=0 | expect=pLV-U6>shKRAS-PGK-Puro -->
```bash
python -m cvd --strategy guide --gene KRAS --guide-type shRNA --target GTTGGAGCTGATGGCGTAGG --marker Puro
```

## 4. An insert that is too big

APC is too large for AAV. `--strict` turns the design error into exit code 1, so scripts and CI
can catch it.

<!-- test: exit=1 | expect=above the roughly 4.7 kb AAV packaging limit -->
```bash
python -m cvd --vector aav --gene APC --strict
```

```text
[error] Cassette between the ITRs is 10,430 bp, above the roughly 4.7 kb AAV packaging limit.
```

## 5. Ranking vector options for a suicide gene

For a suicide gene, the program wants expression limited to tumor cells, so the tumor-selective
promoter wins in AAV.

<!-- test: exit=0 | expect=1    100    AAV         tumor-selective promoter (Survivin) -->
```bash
python -m cvd --strategy sui --gene HSV-TK --setting "in vivo" --rank
```

```text
Rank Score  Vector      Layout                                         Cassette  Capacity
1    100    AAV         tumor-selective promoter (Survivin)            2,121     fits (45%)
2    92     AAV         strong constitutive promoter (EF1a)            3,029     fits (64%)
3    92     AAV         compact promoter, no WPRE (EFS)                1,474     fits (31%)
...
```

Add `--explain 1` to see why option 1 is on top and why it ranks above or below every other
option.

## 6. When no AAV option fits

BRCA1 has a very long coding sequence. Even for in vivo delivery, only the lentivirus options fit,
and the best one is the compact promoter without WPRE.

<!-- test: exit=0 | expect=1    83     Lentivirus  compact promoter, no WPRE (EFS) | expect=too large -->
```bash
python -m cvd --strategy rep --gene BRCA1 --vector aav --setting "in vivo" --rank
```

AAV options that do not fit score 0 and are listed after the options that do.

## 7. Seeing the modifications

Ask for an AAV therapy but start with an EGFP marker linked by IRES and an HA tag. For in vivo
therapy the program removes the tag and the marker, and lists every change with its reason.

<!-- test: exit=0 | expect=Remove the C-terminal tag | expect=Remove the reporter or selection marker | expect=Modifications to option 1 -->
```bash
python -m cvd --strategy rep --gene TP53 --marker EGFP --link IRES --tag HA --setting "in vivo" --explain 1
```

Compare with the plain layouts by adding `--no-optimize`.

<!-- test: exit=0 | expect=Modifications are switched off -->
```bash
python -m cvd --strategy rep --gene TP53 --marker EGFP --link IRES --tag HA --setting "in vivo" --explain 1 --no-optimize
```

## 8. Building a ranked option

`--use-option N` builds the chosen option, with its modifications, as GenBank, FASTA, map and
summary files.

<!-- test: exit=0 | expect=Building the vector for option 1 | expect=pAAV-EF1a>TP53 -->
```bash
python -m cvd --strategy rep --gene TP53 --setting "in vivo" --use-option 1
```

## 9. Aneuploidy on the built-in practice data

The program can write a made-up copy-number file for you to practice on. It is synthetic, not
from a patient.

<!-- test: exit=0 | expect=Aneuploidy score: 14 of 39 scorable arms (7 gains, 7 losses) | expect=Building the vector for suggestion 1 -->
```bash
python -m cvd --example-cnv demo.tsv
python -m cvd --cnv demo.tsv --design-from-cnv 1
```

## 10. A real tumor from open data

A real TCGA cholangiocarcinoma from the cBioPortal PanCancer Atlas, with its copy-number segments in
`examples/tcga_chol_W5-AA2X-01.seg`. The README walks through it in full.

<!-- test: exit=0 | expect=Aneuploidy score: 27 of 39 scorable arms (16 gains, 11 losses) | expect=pAAV-Survivin>ULBP2 -->
```bash
python -m cvd --cnv examples/tcga_chol_W5-AA2X-01.seg --build hg19 --design-from-cnv 1 --setting "in vivo" --rank --use-option 1
```
