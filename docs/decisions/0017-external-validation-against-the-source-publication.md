# ADR 0017 — External validation against the source publication's own files

**Date:** 2026-09-02
**Status:** Accepted
**Task:** P3-T2b
**Implements:** CLAUDE.md hard constraint 5, ADR 0008's reporting rule
**Does not supersede:** ADR 0015 (panel), ADR 0016 (pre-registration), ADR 0007
(detection threshold). Nothing here changes a panel or a threshold.

## Context

### Every check Phase 3 had was internal

P3-T2 shipped with 25 cross-checks and the commit message called two of them
"independent implementations". They are independent *implementations* and they
are not independent *evidence*. `checkpoint_detection.py`, `export_tsv.py` and
`score_signatures.py` all read the same `obs['negprobe']` column, written once
at `export_tsv.py:65`, and all apply the same `q3 > negprobe × 2.0`. Three code
paths, one input, one hypothesis.

What that catches: indexing errors, AOI subsetting errors, gene-symbol
mismatches, transcription slips. What it cannot catch, and would pass silently
in the presence of: the detection rule being wrong, the multiple being mis-set,
`negprobe` being the wrong statistic, or the matrix not being what
`docs/data-provenance.md` says it is.

The second cross-check compares against ADR 0008's reconnaissance table. That
table is this project's own Gate 0 computation. "Published" there means
published in an ADR.

The project's only genuinely external check was ADR 0007's reproduction of the
paper's sequencing-saturation claim, which validates DCC parsing rather than
values, and ADR 0007 says so in terms.

### The paper deposits Source Data

PMID 36216799 deposits its per-figure values. The `Figure 1b-e` sheet of
Supplementary Data 5 is the complete 120-AOI matrix. That makes the whole input
checkable against the authors' file rather than against ourselves — the one
thing 25 more internal assertions could never buy.

## Decision

**P3-T2b fetches three supplementary files under ADR 0005's rules and runs four
checks, ordered so each licenses the next.**

### 1. The matrix, exactly

All 120 × 18,694 = **2,243,280 values are identical** to the deposited Source
Data. AOI label sets match 120/120. The gene sets differ by exactly one row,
`NegProbe-WTX`, which this project holds in `obs`; it is compared rather than
dropped, and it equals `obs['negprobe']` exactly.

This is the first external evidence that this project's inputs are the published
inputs. It is a hard failure, and it has to be: everything else in the file is
an argument about what those numbers mean.

### 2. The 119-vs-120 gap is closed

`docs/limitations.md` recorded that the paper analyses 119 ROIs where GEO
deposits 120, and that **which** AOI the authors excluded was unknown.

The fibrosis sheet carries a single `TBME15` where the full matrix carries
`TBME15a` and `TBME15b`. The rule asserts that column against every candidate:
agreement **1.0000 against TBME15a, 0.0003 against TBME15b**.

**The answer is TBME15b, dropped without comment.** No AOI failed QC upstream;
the paper's 19 TBME "cases" are 20 AOIs with patient 15's second AOI omitted.

Agreement here is within tolerance rather than exact, and that is a finding in
itself: **the two deposited sheets round differently.** `Figure 1b-e` reproduces
GEO bit-for-bit, while `Fig.4-6 & S4,5` carries an extra significant figure for
large values — 20375.33 where GEO has 20375.0. 64 of 18,694 values differ, none
of them near background. Requiring exactness would fail on the publisher's
rounding rather than on anything about the data, so the alias test is
discriminating (≥ 99% for the right candidate, ≤ 50% for every other) rather
than bit-exact. The exact comparison is check 1 and stays exact.

### 3. Count quantisation, measured against the raw counts

The matrix is Q3-scaled counts with no zeros (Q1), so an AOI's values sit on a
lattice whose spacing is one raw count times that AOI's Q3 factor.

**The quantum is the spacing between adjacent distinct values, not the smallest
value.** The smallest value is 1 count only in the shallowest AOIs and is **1–35
counts, median 9**, across the cohort. The first implementation of this check
used the minimum as the quantum, inflating every count by that factor; it failed
visibly — the lattice check dropped to 0.15 and the AOIs it accepted were the
handful where the minimum happened to be 1.

The recovered scale is then **checked against P0-T2's `.dcc` archive**, which
carries integer per-probe counts: `median(value / quantum)` must equal the AOI's
own median raw count. It agrees for **120 of 120 AOIs, median ratio 1.0005**.
No PKC is needed and none exists (ADR 0006) — the comparison is
distribution-to-distribution, so probe identity never has to be resolved.

That makes the count scale a measurement. What it measures:

| compartment | background | CD274 | CD276 | CTLA4 | HAVCR2 | IDO1 | LAG3 | PDCD1 | TIGIT | VSIR |
|---|---|---|---|---|---|---|---|---|---|---|
| `L` | 29 | 95 | 145 | 34 | 54 | 47 | 38 | 64 | 37 | 57 |
| `LB` | 50 | 106 | 255 | 59 | 69 | 58 | 62 | 97 | 57 | 92 |
| `TIME-L` | 29 | 95 | 120 | 81 | 76 | 61 | 61 | 60 | 60 | 93 |
| `TIME-B` | 30 | 78 | 132 | 44 | 73 | 41 | 48 | 53 | 37 | 85 |
| `TBME` | 34 | 83 | 104 | 53 | 73 | 43 | 44 | 57 | 45 | 94 |

(median implied raw counts per AOI; `mLN` and `BC` in the table.)

And the resolution that goes with it: an AOI holds a median of **645 distinct
values across 18,694 genes**, and a panel gene's exact value is shared by a
median of 100–360 other genes, up to **2,607**.

**This is a statement about resolution, where ADR 0007's floor is a statement
about background.** They agree on the substance — `CTLA4` at 34 counts against a
29-count background in `L` is the same fact as 2/30 detected — and the count
version is the one that survives an argument about where the threshold should
sit.

### 4. The published contrast, reproduced exactly

All **99 published log2FC values reproduce to 1.6e-5** and all 99 p-values to
2.0e-4 relative. That establishes the estimator: **mean-of-log2 with Student's
t**, positive = up in F(h). Log2-of-mean matches only 4 of 99, so this was not
guessable and is why the check is a hard failure.

Having reproduced the pipeline, the rule reports what it excludes. Of the 15
checkpoint genes examined — this project's nine plus the six the paper names —
**none reaches the paper's own |log2FC| > 1.5 threshold, and none appears in the
99-gene list.** Largest is VTCN1 at 0.888. The smallest |log2FC| in the published
list is 1.519, so the threshold was applied exactly as stated.

The paper's checkpoint claim — PDCD1LG2, BTLA, VTCN1 and IDO1 elevated in
non-fibrotic TBME, feeding a combination-therapy recommendation — is
reproducible and nominally significant (p = 0.0028–0.0111) at effect sizes of
0.39–0.89 log2, on genes carrying 31–74 counts against a 34-count background,
in an n = 8 vs n = 6 subgroup contrast. It is **not** from the thresholded DEG
analysis whose threshold the Methods state.

## What this does NOT license

- **No panel change.** The six genes the paper names are transcribed into
  `external_validation.py` as the paper's claim and reported beside ours. They
  are not added to `config/checkpoints.yaml`. Panel membership is pre-registered
  (ADR 0015) and a stop-and-ask decision; adding genes after seeing which ones
  move is precisely the failure P3-T2 exists to prevent.
- **No threshold change.** ADR 0016 §3 and ADR 0007 are untouched. The
  |log2FC| > 1.5 threshold is applied only to the paper's own list.
- **No claim about biology.** Check 4 describes the *published* analysis on the
  *published* groups. It is not a result about brain metastasis and must never
  be reported as one. ADR 0008 point 4 is unchanged: a gene below the floor in a
  compartment is "not assessable in <compartment>", never "lower in
  <compartment>".

## Consequences

- **The Option C secondary table gains its strongest column.** The decision
  still open at the top of `docs/NEXT_STEPS.md` is whether to relax ADR 0016 §3.
  Nothing here relaxes it, and the case for keeping it is now concrete: the
  published analysis is what a checkpoint claim looks like when detection is not
  audited first, and this project can now show that rather than assert it.
- **`resources/` gains a third sanctioned writer.** `raw/` (P0-T2, GEO),
  `reference/` (P2-T5, SpatialDecon), now `supplementary/` (P3-T2b). Same
  two-rule split, same `protected()` outputs, same reuse of `acquire_geo.py` and
  `record_provenance.py` unmodified — editing either would fire the code
  rerun-trigger on `p0t2_fetch_geo` and abort the DAG (ADR 0005).
- **`config/external_validation.yaml` is a fourth config file, deliberately not
  in `config.yaml`,** whose SHA-256 is in the `.h5ad` `uns` and whose edit costs
  a Phase 1 + 2 rebuild (ADR 0013). Its path is a literal in
  `05_checkpoints.smk` for the same reason. Phase 3 has now added zero
  `config.yaml` edits beyond the one paid at `f1f55b8`.
- **No new Python dependency.** `.xlsx` is parsed with `zipfile` +
  `ElementTree`. openpyxl would have to go into `py-analysis.yaml`, which is
  `p0t2_fetch_geo`'s env, and that fires the software-env trigger onto protected
  outputs. The env file also warns that anything it lacks falls through to the
  user's base install, which would make this rule pass here and fail everywhere
  else. The reader uses `float()`, which is correctly rounded, where ADR 0013's
  postscript records that pandas' C parser is not — necessary because check 1
  asserts exact equality.
- **Two reader bugs are now guarded by the checks that caught them.** Excel
  omits blank rows from the sheet XML, so a sequentially-appending reader is
  shifted by one wherever a sheet has a spacer row; rows are placed by their
  declared index. And the quantum error above. Both were caught by assertions
  rather than by inspection, which is the argument for writing the assertion
  first.
- **Gate 3's record moves to ADR 0018.**

## Alternatives rejected

**Do the analysis in a notebook and cite it.** Hard constraint 5: a number that
did not come from a rule is not a result. It also would not have caught the
quantum error, because a notebook has no `.dcc` cross-check unless someone
thinks to write one.

**Fetch only the Source Data.** The fibrosis scores are what make the group
assignment the authors' own rather than this project's reconstruction, and the
DEG list is what turns "these genes are not in the list" from an assertion into
a check.

**Assert the alias by exact equality.** Fails on the publisher's rounding.
Discrimination is what the test actually needs, and the observed margin —
1.0000 vs 0.0003 — is three orders wider than any tolerance question.

**Add the paper's six checkpoint genes to the panel.** Rejected on ADR 0015
grounds and it is not close. They are reported, never adopted.
