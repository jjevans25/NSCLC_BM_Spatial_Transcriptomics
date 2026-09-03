# Next steps — Phase 3 (checkpoint landscape by compartment)

**Last session:** 2026-09-03
**Branch:** `P3` — **not pushed, not merged.** Branched from `main`
at `9721279`. **Working tree clean, `snakemake -n` says nothing to do, `--lint`
clean, report builds.** Tomorrow starts from green.
**Gate 0:** PASSED. **Gate 1:** PASSED. **Gate 2:** PASSED (ADR 0014).
**Gate 3:** open — but **both its clauses are already secured** (see Gate 3 below).
**Phase 3:** **P3-T1, P3-T2, P3-T2b, P3-T3, P3-T4, P3-T5 and P3-T6 COMPLETE.**
**P3-T7 is all that remains**, and it is writing.
**Next free ADR number: 0020** (0017 = P3-T2b, 0018 = the Option C decision,
0019 = the P3-T6 slider; **Gate 3's record is now 0020**, moved because ADR 0019
had to exist before the notebook did).

**ONE THING P3-T6 DID NOT VERIFY, and it is the acceptance test ADR 0010 names:
nobody has opened the export in a browser.** Everything mechanical passes — the
figure decoded out of the built `index.html` is pixel-correct, `public/` serves
over HTTP, no `UnhashableStub` renders — but export-time execution is not
Pyodide execution, and ADR 0010's whole point is that the difference is invisible
until a browser runs it. Do this first, it takes two minutes:

```bash
python -m http.server --directory results/reports/checkpoint_explorer
```

Check the four model tabs render as *tables* rather than stub text, and that the
dot plot appears. If something is broken it will be API skew in an unpinned
`marimo` under Pyodide, not the data path.

Phase 3's goal (§6): Aim **A4** — which compartment carries each checkpoint
gene, and whether that shifts brain vs. lung. ADR 0008 demoted it to
**exploratory** at Gate 0 and nothing since has changed that.

---

## START HERE — Phase 3 is COMPLETE; the open item is the merge

`P3-T7` landed: `docs/analysis-notes.md` gained ~700 words on where PD-L1 sits
and what it would mean, `docs/limitations.md` §10 gained the gene-level
confirmation, and **ADR 0020 is the Gate 3 record — GATE 3 PASSES**.

What remains is bookkeeping, not analysis:

1. **Open the P3-T6 export in a browser** (see below). Still the one unverified
   acceptance criterion in the phase.
2. **Squash-merge `P3` into `main` with the gate verdict in the message**
   (ADR 0004). Not yet done — it is a decision, not a chore.
3. Then Phase 4 (A5, crosstalk), which is **exploratory and anticipated null**:
   the ligand side is largely undetected and Phase 3 confirmed the same pattern
   gene by gene (ADR 0014, ADR 0020 §Consequences).

**Gate 3's verdict, in one line:** `CD274` (PD-L1) is enriched in the lung immune
compartment, **+0.612 SD [+0.361, +0.863], q = 0.0001**, surviving both
pre-registered sensitivities — and separately, **40 of 63 gene × compartment
cells are not assessable, the panel clearing the floor in `TIME-L` alone**. Both
Gate 3 clauses, satisfied by different tasks.

**One acceptance criterion could not be checked and was not rounded up.**
P3-T6's "interactive offline" needs a human with a browser, so by §0 it is not
an acceptance criterion as written. The export references zero external origins
and is structurally self-contained, which *is* checkable and does pass. Recorded
as a defect in the plan, not in the work.

---

## RESOLVED — the modelling restriction is Option C (ADR 0018)

The decision that blocked P3-T3 is made. **The primary restriction stands
exactly as ADR 0016 §3 pre-registered it** — a gene enters a primary fit only if
it clears the floor in *both* groups — **and a declared exploratory secondary
carries the genes it excludes**, with no FDR and a per-row background-gradient
direction. ADR 0018 records it, including the two rejected options.

**Then P3-T3 measured something that changes how that decision reads.**

The case for Option B was that where the background gradient is large it runs
*against* the observed direction, so the artefact could only shrink those
effects, never manufacture them. **That was wrong**, and it was wrong in a way
nobody could see until the models ran. It conflated two quantities that move in
opposite directions under the same gradient:

| | rule | higher background → |
|---|---|---|
| **detection** | `q3 > 2 × negprobe` | bar rises, gene detected less → looks **depleted** |
| **expression** | `log2(q3 + 1)` | background adds to signal → looks **enriched** |

The models are on expression. Background is higher in the immune compartments
(`L` 4.96 → `TIME-L` 5.36) and every carrier estimate is positive, so the
gradient is **permissive** — it could have contributed to the very effects
Option B would have promoted on the grounds that it could not.

`carrier_background_sensitivity` quantifies it. Adjusting for `negprobe_log2`
shifts the estimate by a mean of **−0.155** for the genes the primary excludes
and **+0.010** for the genes it admits:

| gene × stratum | primary | bg-adjusted | shift | |
|---|---|---|---|---|
| IDO1, lung | +0.524 | +0.277 | **−0.247** | p 0.059 → 0.329 |
| PDCD1, brain | +0.245 | +0.017 | **−0.227** | |
| TIGIT, lung | +1.036 | +0.840 | −0.196 | |
| HAVCR2, lung | +0.963 | +0.780 | −0.182 | |
| PDCD1, lung | +0.408 | +0.233 | −0.176 | p 0.0008 → 0.029 |
| CD276, brain *(primary)* | +0.159 | +0.245 | +0.086 | |

**The genes the restriction excludes are exactly the genes the adjustment
moves.** That is the pre-registration working. Do not reopen it, and do not
repeat the detection-vs-expression conflation — it is now guarded by a block in
`fit_checkpoint_models.R`'s report and by ADR 0018.

---

---

## P3-T2b — external validation (ADR 0017)

**Added after P3-T2, because P3-T2's 25 cross-checks are not what the commit
message called them.** Three implementations, yes — but all three read the same
`obs['negprobe']` column and test the same hypothesis. They catch indexing and
transcription errors; they cannot catch the rule being wrong. Cross-check 2
compares against ADR 0008's reconnaissance, which is this project's own Gate 0
computation. Before P3-T2b the only external check in the project was ADR 0007's
saturation reproduction, which validates DCC parsing, not values.

The paper deposits Source Data. `p3t2b_external_validation` uses it:

1. **All 2,243,280 values are identical** to the deposited matrix; 120/120 AOI
   labels match; the published `NegProbe-WTX` row equals `obs['negprobe']`
   exactly. This project's inputs are the published inputs — now evidence, not
   assumption.
2. **The 119-vs-120 gap is closed: the paper dropped `TBME15b`.** Agreement
   1.0000 against `TBME15a`, 0.0003 against `TBME15b`. `docs/limitations.md` and
   `docs/data-provenance.md` are updated; the open thread is gone.
3. **Count quantisation, validated against the `.dcc` raw counts** for 120/120
   AOIs (median ratio 1.0005). The panel in raw counts against background:

   | compartment | bg | CD274 | CD276 | CTLA4 | HAVCR2 | IDO1 | LAG3 | PDCD1 | TIGIT | VSIR |
   |---|---|---|---|---|---|---|---|---|---|---|
   | `L` | 29 | 95 | 145 | 34 | 54 | 47 | 38 | 64 | 37 | 57 |
   | `LB` | 50 | 106 | 255 | 59 | 69 | 58 | 62 | 97 | 57 | 92 |
   | `TIME-L` | 29 | 95 | 120 | 81 | 76 | 61 | 61 | 60 | 60 | 93 |
   | `TIME-B` | 30 | 78 | 132 | 44 | 73 | 41 | 48 | 53 | 37 | 85 |
   | `TBME` | 34 | 83 | 104 | 53 | 73 | 43 | 44 | 57 | 45 | 94 |

   An AOI holds a median of **645 distinct values across 18,694 genes**, and a
   panel gene's exact value is shared by a median of 100–360 other genes, up to
   2,607. **This is the same fact as the detection floor, stated as resolution
   rather than as background** — and it is the version that survives an argument
   about where the threshold should sit, which is exactly the argument at the
   top of this file.

4. **The published F(h) vs F(-) contrast reproduces exactly** — all 99 log2FC to
   1.6e-5, all 99 p-values to 2.0e-4. Estimator: mean-of-log2 with Student's t
   (log2-of-mean matches only 4 of 99). Having reproduced it: **none of the 15
   checkpoint genes examined reaches the paper's own |log2FC| > 1.5 threshold**
   (largest, VTCN1, 0.888) and none appears in the 99-gene list. The paper's
   checkpoint claim — PDCD1LG2, BTLA, VTCN1, IDO1 up in non-fibrotic TBME,
   feeding a combination-therapy recommendation — is nominally significant
   (p = 0.0028–0.0111) at 0.39–0.89 log2 on genes carrying 31–74 counts against
   a 34-count background, n = 8 vs 6. It is not from the thresholded DEG
   analysis whose threshold the Methods state.

**No panel change, no threshold change, no biological claim.** The six genes the
paper names are transcribed into the script as the paper's claim and reported
beside ours; they are NOT in `config/checkpoints.yaml` (ADR 0015 is
pre-registered and stop-and-ask). ADR 0016 §3 and ADR 0007 are untouched.

New: `config/external_validation.yaml`, `workflow/schemas/external_validation.schema.yaml`,
`workflow/scripts/external_validation.py`, three rules in `05_checkpoints.smk`,
`resources/supplementary/` (third sanctioned writer, ADR 0005 rules), outputs
`external_validation_matrix.tsv`, `aoi_quantisation.tsv`,
`checkpoint_count_quantisation.tsv`, `published_deg_reproduction.tsv`,
`published_fibrosis_contrast.tsv`, `external_validation_summary.json`.

**Two traps it walked into, both caught by assertions rather than by reading:**

- **Excel omits blank rows from the sheet XML.** A reader that appends
  sequentially is shifted by one wherever a sheet has a spacer row between title
  and header. Rows are now placed by their declared `r` index, and `header_row`
  in the config is a true sheet position.
- **The count quantum is the SPACING between adjacent distinct values, not the
  smallest value.** The smallest value is 1 count only in the shallowest AOIs
  and is 1–35 counts (median 9) across the cohort. The first version used the
  minimum and inflated every count by that factor. **This is why the `.dcc`
  cross-check exists** — the lattice can look clean at the wrong spacing, and
  only the raw counts settle it. Do not weaken that check.

**`openpyxl` must not be added to `py-analysis.yaml`.** That is
`p0t2_fetch_geo`'s env and its outputs are `protected()`; a dependency change
fires the software-env trigger and aborts the DAG. `.xlsx` is parsed with
`zipfile` + `ElementTree`, and the reader uses `float()` (correctly rounded)
where pandas' C parser is not — check 1 asserts exact equality.

---

## What Phase 3 has completed

### P3-T1 — the panel (commit `f1f55b8`)

`config/checkpoints.yaml` holds PROJECT_PLAN §6's nine genes: `CD274`, `PDCD1`,
`CTLA4`, `LAG3`, `HAVCR2`, `TIGIT`, `IDO1`, `VSIR`, `CD276`. Each carries
`alias`, a controlled `trial_stage` enum, one line of clinical `rationale` and a
`source`; `source_note` states the gap wherever the citation does not establish
exactly what the other fields claim (LAG3's approval is in melanoma, IDO1's
phase 3 is in melanoma, CD274's cited trial uses an anti-PD-1 agent). All nine
resolve against the 18,694 measured genes.

**Citation verification is not a formality — it caught more than Phase 2's did.**
Three of the nine PMIDs first proposed pointed at entirely unrelated papers: a
digital-health regulatory guide (proposed for TIGIT), a review of autoimmune
encephalitides (VSIR), and a paper on the hydroaminomethylation of α-olefins
(CD276). A fourth error surfaced one level down — the CD276 replacement's first
author is Malapelle, not the name first written. Every citation is now verified
against PubMed by title, journal, year **and first author**. ADR 0011 said a
schema can enforce that a source is present, never that it is correct; ADR 0015
records the second, worse instance. **Verify any future addition the same way.**

Also landed: `workflow/schemas/checkpoints.schema.yaml` (source pattern admits
`PMID:` and `NCT`, deliberately not `MSigDB:`), the `CHECKPOINTS` loader in
`common.smk`, `p3t1_resolve_checkpoints`, `workflow/scripts/resolve_checkpoints.py`,
ADR 0015 and ADR 0016.

**The `config.yaml` edit is PAID FOR.** It happened once, in `f1f55b8`, and it
forced the expected Phase 1 + 2 rebuild. **Do not touch `config/config.yaml`
again for the rest of Phase 3** — everything P3-T3 … T7 needs is already in the
`checkpoints:` block. As a free regression check, that rebuild returned **36 of
37 tables byte-identical**; the only diff was `h5ad_summary.json`'s `git_sha`
and `config_sha256`, which is the provenance working correctly.

Guards were negative-tested: 13 deliberate violations, all rejected, real files
still valid — each of the six formula fields with `(1|patient_id)` removed, a
detection fraction > 1, `fdr_method: bonferroni`, a stray key, the block deleted
entirely, a gene with no source, a free-text source, an `MSigDB:` source, a
missing rationale, `trial_stage: phase_4`, and an unknown per-gene key.

### P3-T2 — the detection audit (commit `3d1c4b4`)

`workflow/scripts/checkpoint_detection.py` + `p3t2_checkpoint_detection`. It
reuses the project's **one** detection rule — `q3 > 2.0 ×` that AOI's
`NegProbe-WTX` (ADR 0007) — lifted from `score_signatures.py`, not rewritten.

**Keyed on `aoi_code`, never on `compartment`.** `compartment` is degenerate
across sites — `L`, `LB` and `mLN` are all labelled `tumour` — so grouping by it
would silently pool lung and brain tumour AOIs and destroy the contrast the
phase exists to measure. Anything written from here on must key the same way.

**All 25 cross-checks agree exactly**, two of them against independent
implementations:

- 9 per-gene totals vs `var['detected_in_n_aoi']`, computed by P0-T7's
  `export_tsv.py` with the same rule but separate code.
- 16 `TIME-L`/`TIME-B` counts vs ADR 0008's Gate 0 reconnaissance table, **plus**
  its published negprobe ratios, which match to the published 2 dp.

Outputs: `checkpoint_detection.tsv` (63 rows, 9 genes × 7 compartments),
`checkpoint_expression.tsv` (720 rows, the four modelled compartments only —
restricting the file is what makes accidental pooling impossible rather than
merely discouraged), `checkpoint_detection_summary.json`.

### The result

Detected AOIs per gene per compartment, at 2× background:

| gene | alias | `L` | `LB` | `mLN` | `TBME` | `TIME-L` | `TIME-B` | `BC` |
|---|---|---|---|---|---|---|---|---|
| CD274 | PD-L1 | 19/30 | 17/27 | 11/13 | 7/20 | **13/15** | **5/8** | 0/7 |
| CD276 | B7-H3 | 29/30 | 26/27 | 13/13 | 13/20 | **15/15** | **8/8** | 1/7 |
| CTLA4 | CTLA-4 | 2/30 | 2/27 | 4/13 | 3/20 | **10/15** | 2/8 | 0/7 |
| HAVCR2 | TIM-3 | 2/30 | 2/27 | 5/13 | 5/20 | **10/15** | **7/8** | 0/7 |
| IDO1 | IDO1 | 9/30 | 4/27 | 5/13 | 0/20 | **8/15** | 2/8 | 0/7 |
| LAG3 | LAG-3 | 2/30 | 2/27 | 2/13 | 2/20 | 6/15 | 3/8 | 0/7 |
| PDCD1 | PD-1 | 12/30 | 14/27 | 7/13 | 6/20 | **8/15** | 3/8 | 0/7 |
| TIGIT | TIGIT | 0/30 | 2/27 | 3/13 | 1/20 | **8/15** | 1/8 | 0/7 |
| VSIR | VISTA | 8/30 | 11/27 | 6/13 | 12/20 | **14/15** | **7/8** | 6/7 |

Panel assessability at the pre-registered 0.5 / 0.5 floors — **the panel clears
in one compartment out of seven**:

| compartment | assessable | not assessable there |
|---|---|---|
| `TIME-L` | **8/9** | LAG3 |
| `TIME-B` (n = 8) | 4/9 | CTLA4, IDO1, LAG3, PDCD1, TIGIT |
| `LB` | 3/9 | CTLA4, HAVCR2, IDO1, LAG3, TIGIT, VSIR |
| `mLN` | 3/9 | the same six |
| `L` | 2/9 | those six plus PDCD1 |
| `TBME` | 2/9 | seven, including CD274 |
| `BC` | 1/9 | eight |

**The lung immune compartment is the only place this panel can be read.** That
is the deliverable Gate 3 anticipated, not a failure.

QC-flagged AOIs were retained and scored (flag-don't-drop) with a
`detection_rate_qc_clean` column alongside. Excluding `TBME`'s four flagged AOIs
moves its per-gene detection by at most 0.15 and **flips no assessability
verdict** — so `TBME`'s deficit is the compartment, not those AOIs. That
question is now settled and does not need re-opening.

### Two records added while doing it

- **ADR 0016 gained §5:** the P3-T5 dot plot colours by `median_negprobe_ratio`,
  not §6's literal "mean expression". A mean-log2 scale renders *undetected*
  genes **brighter in brain than in lung**, which carries the ADR 0008 artefact
  straight into the deliverable figure. The background-normalised ratio is
  immune by construction and is the statistic ADR 0008 itself reported. Both
  mean-log2 columns are still in the table.
- **ADR 0013 gained a postscript:** `pandas.read_csv`'s default C parser is
  **not correctly rounded**, the same property §3 recorded for R's `R_strtod`.
  Reading P3-T2's own `%.17g` output back reported 32 of 720 values unequal —
  the file was exact, the reader was not. `float_precision="round_trip"`
  recovers it. **Any Python rule needing exact values back from a project TSV
  must pass it.** No existing result is affected (Phase 2's tolerance is 1e-4
  against a ~1e-16 discrepancy), and the Phase 2 scripts were deliberately left
  unmodified rather than fire their rerun triggers for no numerical gain.

---

### P3-T3 — the carrier models (ADR 0018)

`fit_checkpoint_models.R` + `p3t3_carrier_models`, and
`crosscheck_checkpoint_mixedlm.py` + `p3t3_checkpoint_crosscheck`.
`expression ~ compartment + (1|patient_id)` within lung (`L` vs `TIME-L`) and
within brain (`LB` vs `TIME-B`), never pooled — `compartment` is degenerate
across sites, so subsetting by site first is what makes it a clean two-level
factor. Reference level is the tumour compartment and the script asserts it.

**The script is parameterised and P3-T4 drives it unchanged** — pass
`contrast_var="site"`, `reference_level="lung"`, `test_level="brain"`,
`strata_var="compartment"`. Do not write a second copy; a second copy is a
second place for the reference level to drift.

**Primary — 2 of 9 genes per site, exactly as pre-registered.** The only Phase 3
table with a q-value.

| stratum | gene | estimate [95% CI] | p | q | detection |
|---|---|---|---|---|---|
| lung | **CD274** | **+0.612 [+0.361, +0.863]** | 0.0001 | **0.0001** | 19/30 → 13/15 |
| lung | CD276 | +0.147 [−0.139, +0.432] | 0.290 | 0.290 | 29/30 → 15/15 |
| brain | CD274 | +0.222 [−0.384, +0.829] | 0.461 | 0.461 | 17/27 → 5/8 |
| brain | CD276 | +0.159 [−0.190, +0.508] | 0.325 | 0.461 | 26/27 → 8/8 |

**CD274 (PD-L1) is enriched in the lung immune compartment, q = 0.0001**, and it
survives both sensitivities (batch +0.599, background +0.635). It is the one
result Gate 3 asked for: an interpretable compartment-resolved pattern with a CI
worth showing. The brain estimate is the same sign at a third the size with
`TIME-B` n = 8 and a 1.1–1.3 SD power floor — uninformative, not negative.

**Exploratory — 8 of 9 in lung, 5 of 9 in brain, no FDR** (ADR 0018). Large
apparent immune enrichment for CTLA4 +1.425, VSIR +1.336, TIGIT +1.036,
HAVCR2 +0.963 in lung. **Every one is a gene detected in 0–9 of 30 tumour AOIs**,
and the background adjustment moves them (mean −0.155). Not reportable as a
finding; the detection counts are.

The overlap assertion passed: the 12 rows fitted in both tables agree to 0.0
(max |difference|). 21 singular fits, surfaced — not a reason to drop
`(1|patient_id)` (ADR 0009 §1).

Cross-check: lme4 vs statsmodels on CD276 agree to **7.3e-9** (brain) and
**1.1e-6** (lung). **The SE is reported but NOT gated**, matching
`p2t3_model_crosscheck`, and P3-T3 measured why: the two libraries use different
fixed-effect SE estimators, agreeing only when the variance components are
sharply determined — they are not, with 17/30 lung and 19/27 brain patients
contributing one AOI. statsmodels reaches the same variance components under
lbfgs, powell and bfgs, so it is not an optimiser artefact. **The SE, CI and
p-value of record are lmerTest's**; recomputing the CI in statsmodels gives a
wider one (lung SE 0.160 vs 0.134).

Outputs: `checkpoint_carrier_models.tsv`, `checkpoint_carrier_exploratory.tsv`,
`checkpoint_carrier_summary.json`, `checkpoint_carrier_crosscheck.json`.


### P3-T4 — the shift models and the paired check

`p3t4_shift_models` drives `fit_checkpoint_models.R` **unchanged** with
`contrast_var="site"`, plus `p3t4_checkpoint_crosscheck` and
`p3t4_checkpoint_paired_check` (`checkpoint_paired_check.py`, a new file — an
edit to `paired_check.py` would fire the rerun-trigger on `p2t4_paired_check`).

**Primary — 4 genes in immune, 2 in tumour. Everything is null.**

| stratum | gene | brain − lung [95% CI] | q | detection |
|---|---|---|---|---|
| immune | CD274 | −0.260 [−0.595, +0.075] | 0.462 | 13/15 → 5/8 |
| immune | VSIR | +0.150 [−0.141, +0.440] | 0.512 | 14/15 → 7/8 |
| immune | HAVCR2 | +0.180 [−0.241, +0.601] | 0.512 | 10/15 → 7/8 |
| immune | CD276 | +0.072 [−0.217, +0.362] | 0.595 | 15/15 → 8/8 |
| tumour | CD276 | −0.043 [−0.247, +0.161] | 0.875 | 29/30 → 26/27 |
| tumour | CD274 | +0.023 [−0.278, +0.324] | 0.875 | 19/30 → 17/27 |

**`TIME-B` n = 8, power floor 1.1–1.3 SD. Every one of these is far below it, so
they are uninformative, not negative.** CD274 leans brain-lower (−0.26, and
−0.247 under the background adjustment), the same direction as Phase 2's
antigen-presentation reduction — but at a quarter the power floor it is a
direction, not a result.

**Every row's gradient is `negligible`**, as ADR 0018 predicted: the site
background deltas are +0.03 (tumour) and −0.06 (immune) against 0.31–0.40 for
the compartment contrasts. **The background artefact is a P3-T3 problem, not a
P3-T4 problem** — measured, not assumed.

**Paired check — direction only, no test statistic in EITHER compartment.**
The immune set is the 5 patients §2.3 names (P5, P12, P15, P19, P35); the tumour
set is 23. The tumour set would support a paired test and **deliberately does not
get one** — a test introduced when it becomes available, on a contrast nobody
pre-registered a test for, is the move ADR 0018 rejected. The script asserts no
p-value, q-value or CI column can reach the table.

Direction agrees in **9 of 11** comparable cells (immune 6/8, tumour 3/3); 7 had
no unpaired fit. The two disagreements are immune CD276 (paired −0.161 vs
unpaired +0.072) and immune PDCD1 (+0.032 vs −0.055) — both effects within noise
of zero at n = 5, which is what a disagreement at this size means.

**The cross-check confirmed why the SE is not gated.** Same gene, same
tolerance: the SE gap is **17.1% in immune** (23 AOIs, 5 paired patients) and
**0.4% in tumour** (57 AOIs, 23 paired), while the estimates agree to 4.7e-6 and
2.4e-6. Where the random intercept is well identified the two libraries
converge; where it is not, they do not — the prediction made at P3-T3, tested on
a contrast it was not derived from.

Outputs: `checkpoint_shift_models.tsv`, `checkpoint_shift_exploratory.tsv`,
`checkpoint_shift_summary.json`, `checkpoint_shift_crosscheck.json`,
`checkpoint_paired_deltas.tsv`, `checkpoint_paired_concordance.tsv`,
`checkpoint_paired_summary.json`.


### P3-T5 — the deliverable dot plot

`workflow/scripts/checkpoint_dotplot.py` + `p3t5_checkpoint_dotplot` +
`workflow/report/checkpoint_dotplot.rst`. Genes × (compartment × site), all
seven compartments, all 120 AOIs. Dot area = detection rate; colour =
`median_negprobe_ratio` (ADR 0016 §5); ring = detected in 0 AOIs; hatch = below
the pre-registered floor. **No config change, no schema change, no ADR** — every
number it uses was already pre-registered.

**The colour midpoint is read from `qc.detection_background_multiple`, not
hardcoded**, so it cannot drift from the rule it depicts. `TwoSlopeNorm` centred
there on `PRGn` — deliberately not P2-T7's `RdBu_r`, which means a *signed
effect*; two deliverables must not share a colour language for different
quantities.

**40 of 63 cells hatched, 23 clear, 9 rings.** Six assertions pass, including
re-deriving `assessable` and `modelled` at plot time and checking P3-T2's
`figure_colour_basis` contract string. The log reprints the 9×7 grid so the PNG
is verifiable without opening it.

**Two properties logged as observations, never asserted:** `ratio > 2.0` and
`assessable` agree on all 63 cells (they are the same rule seen twice, and can
legitimately diverge only where `n_aoi` is even and the rate sits exactly at the
floor — no such cell exists); and 0 cells change assessability under
`detection_rate_qc_clean`.

**Three things went wrong in the build and are worth not repeating:**

- **The title clipped at both canvas edges on the first render** — the exact
  failure P2-T7's docstring records, reproduced despite the plan warning about
  it. Fixed by splitting into seven short lines. **Open the PNG; do not assume.**
- **The colourbar, the right-hand `k/7` axis and the key all collided.** The
  twinx attaches to the main axes' right spine, so anything immediately beside
  it lands on those labels. Layout is now `main | key | colourbar` with
  `wspace=0.34`, and the threshold/background annotations live in the tick
  labels rather than as free text beside the bar.
- **Greyscale sends both ends of a diverging map to dark grey**, so 0.84 and
  5.11 look alike without colour. Documented in the script and caption rather
  than designed away: every claim is carried redundantly by dot area, the hatch,
  the bold edge and the printed counts, so a greyscale reader loses the ratio
  and nothing else.

Verified: `--lint` clean, dry run proposes exactly one job, report builds (so the
`.rst` is valid), and the PNG is byte-identical on rebuild.


### P3-T6 — the marimo explorer (ADR 0019)

`notebooks/apps/checkpoint_explorer.py` + `p3t6_checkpoint_explorer` +
`workflow/report/checkpoint_explorer.rst`. A detection-floor slider, a
compartment multiselect and a site filter over the P3-T5 dot plot, with the
P3-T3/T4 model tables in their own section below it. Exported to WASM HTML,
`results/reports/checkpoint_explorer`, in the report under **Interactive**.

**ADR 0019 was written before the notebook, because a slider over a
pre-registered threshold is either a serious problem or the most useful thing in
the phase depending on what it is for, and nothing in the repo said which.** The
floor of record stays 0.5; the app defaults there, reproduces the static figure
there, and shows a banner the instant it leaves. `qc.detection_background_multiple`
is deliberately NOT exposed — `detection_rate` is already computed at it, so a
control over it could only relabel an axis, and recomputing detection in a
notebook is forbidden twice over.

**Model estimates are adjacent to the dot plot, never on it** (ADR 0018). The
plot is detection, the tables are expression, and the two move opposite ways
under the same background gradient. Four tabs: primary (the only Phase 3 output
with a q-value), the two sensitivities, the genes the restriction **excluded** —
the exclusion is the finding — and the exploratory table under a standing
no-claim banner.

**Defaults were verified by diff, not by squinting.** A scratchpad harness ran
`app.run()` and compared the app's derived sets against
`results/logs/p3t5_checkpoint_dotplot.log`: column order, row order, the `k/7`
sort key, the 40-cell hatch set, the 9-cell ring set, the full 9×7 grid, the
40/23/9 counts, and both self-checks. **8 of 8 matched.** The same harness drove
the app to floors 0.5/0.6/0.75/0.9/1.0 and through every site filter and an
empty selection, which is where the crashes would have been.

**One real bug, and it was only visible in the export:** the WASM sandbox renders
under a **dark matplotlib theme**, so the figure came back with black row bands
and white tick labels — decoded straight out of the built `index.html`. The
notebook now resets `plt.rcParams` to matplotlib's defaults and sets explicit
facecolors before building the figure, so ambient rcParams cannot decide what a
deliverable looks like. **`results/reports/landscape_explorer` has the same black
background and was deliberately left alone** — it is P1-T5's notebook and editing
it fires that rule's rerun trigger. Worth fixing when something else touches
Phase 1.

**Still unverified: nobody has opened it in a browser.** See the top of this file.


## What remains

### P3-T7 — clinical interpretation (~2 h)

~600 words appended to `docs/analysis-notes.md`; P2-T6's section is the format.

Four things it is *required* to say, each because a reader would otherwise draw
the wrong conclusion from silence:

- **The glial compartment was excluded from the models, and why** (ADR 0016 §2)
  — so nobody concludes `TBME` was overlooked.
- **The Option C split** (ADR 0018 §Consequences): that the primary answers for
  two genes per contrast, that a separate exploratory table exists, why the
  excluded genes were excluded, and **that the exclusion is itself the finding**.
- **`TIME-B` n = 8 inline on every brain claim** (hard constraint 8), with the
  1.1–1.3 SD power floor, so the P3-T4 nulls read as uninformative rather than
  negative.
- **Q2's caveat**: a `TIME` AOI is the PanCK-negative segment of an ROI sited in
  a CD45-rich region, not a CD45-sorted population. The compartment label is not
  a cell-type label.

Then `docs/limitations.md` §10, ADR 0020 as the Gate 3 record, `/gate 3`, and the
squash-merge with the gate verdict in the message (ADR 0004).

---

## What to reuse

| need | reuse |
|---|---|
| detection per gene × group | **`checkpoint_detection.tsv` already has it.** Do not recompute |
| mixed models | `workflow/scripts/fit_checkpoint_models.R` (Phase 3's own; `fit_signature_models.R` was its template) |
| implementation cross-check | `workflow/scripts/crosscheck_checkpoint_mixedlm.py` — **estimate only**; SE, df and p are different estimators and are reported, not gated |
| paired direction check | `workflow/scripts/checkpoint_paired_check.py` — direction only, no test statistic in either compartment |
| rendering not-assessable | `workflow/scripts/checkpoint_dotplot.py` (dots) or `contexture_heatmap.py` (cells) |
| app-tier WASM export | `p1t5_landscape_explorer` in `03_landscape.smk` — read ADR 0010 first |
| exact TSV reads in Python | `pd.read_csv(..., float_precision="round_trip")` — ADR 0013 postscript |
| reading an `.xlsx` with no new dependency | `read_sheet` / `sheet_frame` in `workflow/scripts/external_validation.py` |
| raw counts per AOI | `dcc_median_counts` in the same file — the `.dcc` archive, no PKC needed |
| an external comparator at all | `results/tables/external_validation_summary.json` (ADR 0017) |
| model fits, either contrast | `workflow/scripts/fit_checkpoint_models.R` — parameterised, drives both P3-T3 and P3-T4 |
| another app-tier notebook | `notebooks/apps/checkpoint_explorer.py` — ADR 0010's data path plus the rcParams reset the WASM sandbox needs |
| a control over a pre-registered number | ADR 0019's four parts: default at the registered value, no write path, self-labelling off-default, expose nothing that cannot actually be recomputed |
| the numbers for P3-T5 and P3-T7 | `checkpoint_carrier_models.tsv` (primary) and `checkpoint_shift_models.tsv`. The `_exploratory` pair carries no q-value and no claim may rest on it |

---

## Stop and ask

1. **The modelling restriction** — **DECIDED, Option C (ADR 0018).** The
   primary stays as ADR 0016 §3 pre-registered it; a declared exploratory
   secondary carries the rest, with no FDR. Reopening it is a new stop-and-ask,
   and P3-T3 measured that the argument for reopening it was inverted.
2. **The panel membership** (ADR 0015). Adding or removing a gene remains a
   scientific decision, *especially* now that the audit shows which ones are
   undetected.
3. **The detection floor** (ADR 0016). Pre-registered. Changing it after seeing
   the audit is the failure P3-T2 was designed to prevent. **P3-T6's slider is
   not an exception and does not weaken this** — it is a sensitivity display
   with no write path (ADR 0019), and the floor of record is still 0.5.
4. **Any change to the random-effects structure.** `(1|patient_id)` is mandatory
   on non-independence grounds (ADR 0009) and the config schema enforces it with
   a regex. Singular fits are *not* a reason to drop it — half of Phase 2's fits
   were singular and it stayed.
5. **Dropping an AOI.** Flag, don't drop; exclusions go in
   `results/tables/qc_excluded.tsv` with a reason.

## Gate 3

> **GATE 3** — At least one checkpoint gene shows an interpretable
> compartment-resolved pattern with a CI you're willing to show a stranger. If
> everything is not-assessable, that's still a legitimate finding.

**BOTH CLAUSES ARE SECURED. P3-T7 records the verdict; it does not have to go
looking for one.**

**First clause — CD274 (PD-L1) in the lung immune compartment:**
**+0.612 SD [+0.361, +0.863], q = 0.0001**, from
`checkpoint_carrier_models.tsv` (P3-T3, primary). Survives both pre-registered
sensitivities — batch +0.599, background +0.635 — and the statsmodels refit
agrees with lme4 to 1.1e-6. That is a compartment-resolved pattern with a CI
worth showing a stranger.

**Second clause — the audit.** 40 of 63 gene × compartment cells are not
assessable; the panel clears the floor in `TIME-L` alone. That is a legitimate
finding *because* the detection table licenses it, and the table is now
cross-validated three ways: against `export_tsv.py`'s independent
implementation, against ADR 0008's Gate 0 reconnaissance, and — since P3-T2b —
against the source publication's own deposited Source Data, where **all
2,243,280 values match exactly**.

**What P3-T7 must not do.** The lung-vs-brain shift (P3-T4) is null in every
assessable gene, and with `TIME-B` n = 8 against a 1.1–1.3 SD power floor those
nulls are **uninformative, not negative**. Writing them up as evidence of no
difference would fail the gate more surely than having no result at all.

---

## This session's commits, newest first

| commit | task |
|---|---|
| *(this session)* | **P3-T6** — the marimo explorer; the slider is a sensitivity display (ADR 0019) |
| `9754b81` | **P3-T5** — the deliverable dot plot; detection, deliberately not effects |
| `dab467f` | **P3-T4** — shift models, paired check, a confirmed and uninformative null |
| `322e27a` | **P3-T3** — carrier models; Option C, and the gradient argument was inverted |
| `5b657cd` | **P3-T2b** — external validation against the source publication's own files |

Four ADR-level facts these established, in rough order of how much they change
what the project can say:

1. **All 2,243,280 values of `layers['q3']` are identical to the paper's
   deposited Source Data** (ADR 0017). Every check before this was internal.
2. **The 119-vs-120 gap is closed: the paper dropped `TBME15b`** — recorded for
   two phases as "not recoverable from GEO", which was true, and recoverable
   from Source Data, which GEO does not carry.
3. **Detection and expression move OPPOSITE ways under the same background
   gradient**, so the argument for relaxing the pre-registered restriction was
   inverted (ADR 0018). Measured, not reasoned: adjusting for `negprobe_log2`
   shifts the excluded genes by −0.155 and the admitted ones by +0.010.
4. **The count scale is measured**, validated against the `.dcc` raw counts for
   120/120 AOIs. `TIME-B` background is ~30 counts and CTLA4/TIGIT/IDO1 carry
   37–44 — the detection floor restated as resolution.

---

## State you'll have forgotten

- **`--conda-prefix "$HOME/nsclc-envs"` is MANDATORY on this machine**, alongside
  `--use-conda`. That symlink points *at* `.snakemake/conda`, so nothing moves —
  but this working directory contains a space, and **three** separate scripts
  interpolate the conda prefix unquoted. Recreate with
  `ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"`. **Do not** relocate the
  envs — that fires the software-env trigger on `p0t2_fetch_geo` and its
  `protected()` outputs abort the DAG. ADR 0013.
- **`config/config.yaml` must not be edited again this phase.** Its SHA-256 goes
  into `uns`, so any edit rebuilds the `.h5ad` and re-runs Phases 1 and 2.
  Phase 3 paid that once in `f1f55b8`. `config/checkpoints.yaml` is *not* an
  input of that rule and is free to edit.
- **A clean Phase 3 dry run is `p3tN… + all` only.** If `snakemake -n` proposes
  anything from Phase 0/1/2, `config.yaml` was touched by accident — stop rather
  than spend 20 minutes.
- **`rule p2t0_repair_r_env` must run before any R rule.** Declare
  `results/interim/env_repair/<env>.ok` as an input of every R rule. Both R envs
  are built and verified: `r-stats` (lme4 2.0.6 / Matrix 1.7.5) and `r-geomx`.
  `r-stats.yaml`'s header already says "Phases 2-3" — Phase 3 needs no new env.
- **A green `--conda-create-envs-only` is not evidence an env works.** Load the
  libraries.
- **Neither R's nor pandas' default float parser is correctly rounded** (ADR 0013
  §3 and its P3-T2 postscript). Exactness is a property of the *reader*. Assert
  structure across the boundary; check numbers end-to-end.
- **`py-analysis.conda-lock.yml` is stale** and `conda-lock` is not installed.
  Still the P6-T1 blocker. Phase 3 has added no new Python dependency.
- **P6-T1 must be run on a path containing a space**, or it will not exercise
  `p2t0_repair_r_env` and will report a false pass.
- **`Markdowns/PROJECT_PLAN.md` is gitignored** — ADRs are the durable record.

## Useful commands

```bash
conda activate nsclc_bm_spatial
ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"   # once per machine
snakemake -n --use-conda --conda-prefix "$HOME/nsclc-envs"   # must stay clean
snakemake --lint                                             # must stay clean
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" --cores 4
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" \
  --report results/reports/workflow-report.html
uvx marimo check notebooks/apps/checkpoint_explorer.py
python -m http.server --directory results/reports/checkpoint_explorer
/gate 3
```

---

## Phase 2, for reference

Gate 2 (ADR 0014): the published direction is reproduced. Brain minus lung,
within `TIME`, `TIME-B` n = 8:

| signature | coverage | estimate (z-score) | q | |
|---|---|---|---|---|
| antigen presentation | 13/13 both | −0.477 SD [−1.030, +0.075] | 0.098 | direction recovered, underpowered |
| cytotoxicity | 2/4 both | −0.946 SD [−1.590, −0.302] | **0.018** | recovered, significant |
| myeloid M2 | 6/9 both | +0.091 SD [−0.333, +0.514] | 0.660 | null |

−0.48 SD sits below the 1.1–1.3 SD power floor, so antigen presentation's
non-significance is uninformative, not negative. Paired direction agreed 12/12;
lme4 vs statsmodels 5e-9.

Three of six signatures failed their coverage floor — `exhaustion` 1/6 and `tls`
1/5 in brain, `myeloid_m1` 2/10 at *both* sites. `exhaustion` and `tls` returned
large, nominally significant "reductions in brain" that are **not reportable**.
**P3-T2 has now confirmed the same pattern gene by gene for the checkpoint
panel**, which overlaps `exhaustion` in five of nine genes.
