# ADR 0009 — `(1|patient)` is mandatory, but not for the reason the plan assumed

**Date:** 2026-08-26
**Status:** Accepted
**Task:** P1-T6 (records P1-T3 and P1-T4)
**Supersedes:** the ADR sketched as `0002-mixed-models.md` in PROJECT_PLAN §6
**Gate:** this is the Gate 1 record

## Context

PROJECT_PLAN §6 anticipated this ADR as a formality: measure the variance split,
find patient dominant, and `(1|patient)` follows obviously. The plan says as much
elsewhere — §6's P1-T4 expects that patients "cluster together across
compartments (*they usually do*)".

**Both expectations are wrong in this dataset, and the number the ADR was
supposed to lean on is not there.** P1-T3 and P1-T4 measured it two independent
ways and agree.

The number `0002-mixed-models.md` is also taken — ADR 0002 is the `.ipynb`
PreToolUse hook — so this is 0009.

## What was measured

**P1-T3** (`variance_partition.tsv`, `variance_partition_summary.json`): one
crossed random-effects model per gene, REML, over
patient / site / compartment / `dsp_run` + residual, on all 18,694 genes
(18,691 converged), plus PVCA over the 14 leading PCs (60.4% of variance).

| Component | Per-gene median | PVCA | Permutation null |
|---|---|---|---|
| compartment | 0.211 | **0.412** | 0.003 |
| patient | 0.124 | 0.218 | 0.026 |
| `dsp_run` | 0.078 | 0.096 | 0.002 |
| site | 0.006 | 0.023 | 0.003 |
| residual | 0.495 | 0.252 | 0.889 |

**P1-T4** (`cluster_agreement.tsv`, `sample_clustering_summary.json`): AOI–AOI
Spearman correlation over the same 2000 HVGs, average linkage on 1 − ρ. Peak
adjusted Rand index against the dendrogram cut — compartment **0.741** (k = 7),
site 0.034, patient **0.003**. At k = 7, 34% of the 157 AOI pairs sharing a
patient fall in the same cluster, against a **37%** baseline over all pairs:
patient co-membership is at or slightly *below* chance.

**Gate 1's sentence: 41% of variance is compartment-attributable.**

## Decision

### 1. `(1|patient)` is mandatory in every downstream model

It is **not** justified by variance share. Patient is 12% per-gene and 22% by
PVCA — real (5–8× its permutation null) but nowhere near dominant, and it
produces no cluster structure at all.

It is justified by **non-independence**, on two legs:

- **Design (P0-T3).** AOIs within a patient are not independent by construction:
  P12 and P24 each contribute two `TIME-L`, P15 two `TBME`, and five patients
  (5, 12, 15, 19, 35) contribute paired immune AOIs at both sites. 33 of 42
  subjects contribute more than one AOI. Treating 120 AOIs as 120 independent
  observations overstates the effective n regardless of how small the patient
  variance component is.
- **Direct measurement (P1-T4).** Restricted to the **102 same-patient AOI pairs
  whose compartments differ** — so no part of the contrast can be compartment
  agreement wearing a patient label — mean ρ is **0.425 against 0.318** for
  different-patient pairs, **+0.107**. This is pairwise evidence of correlation
  within a subject that owes nothing to a variance decomposition, and it is
  measured *across* compartments, which is exactly the structure the random
  intercept absorbs.

**These two findings are not in tension and must not be reported as if one
refutes the other.** Patient ARI ≈ 0 while within-patient correlation is
elevated: compartment owns the top-level partition, and patient is a real
second-order effect that never becomes cluster structure. A reader who saw only
T4's ARI would conclude the random intercept could be dropped. It cannot.

### 2. A naive η² is not reportable for any factor in this design

`patient_id` has **42 levels across 120 AOIs**, whose η² null expectation is
already **0.345**. The observed median naive η² for patient is **0.447** —
barely clear of its own floor. Reported naively that reads as "patient explains
45% of variance"; the penalised estimate is 0.124 against a permutation null of
0.026, so the crude number overstates it roughly 3.5×.

Every variance share this project reports is accompanied by a null baseline.
Both are computed and stored: the analytic η² floor `(k−1)/(n−1)`, and an
empirical permutation null (labels shuffled independently, model refit).

### 3. `dsp_run` stays in the variance model; batch is *not* a headline confound

P1-T2 found PC1 shifting consistently positive in run B within every compartment
spanning both runs (L +19.8, LB +9.7, TIME-L +4.7, TIME-B +7.6), so omitting
batch would push that variance into `compartment` — the number the gate turns on.

Quantified (sensitivity S2, HVG subset): dropping the term moves **+0.014** into
compartment and **+0.023** into patient. Real, directionally as predicted, and
about one and a half percentage points — not the distortion the PC1 shifts
implied. Most of what it holds falls to *patient*, which is what a term
near-nested in patient should do.

Two caveats travel with this term permanently:

- It has **two levels** and is **near-nested in patient** — only 1 of 42 patients
  spans both runs — so its own component is poorly estimated.
- It is estimable only in `L`, `LB`, `TIME-L` and `TIME-B`. `mLN`, `TBME` and
  `BC` sit wholly within one run (Q3), so **the partition is conditional on
  that** and any statement about those three compartments cannot separate batch
  from biology.

### 4. `site` and `compartment` are reported together or not at all

Site carries 0.006 per-gene and 0.023 by PVCA — at or barely above its null
everywhere. The two are structurally confounded (`lymph_node` is tumour-only;
`glial_stroma` and `normal_control` are brain-only), so a model carrying both
splits their shared variance by convention rather than evidence. Sensitivity S1
collapses them into the single 7-level `aoi_code`: 0.190, against 0.205 for
compartment alone. **Practically all the site signal is inseparable from
compartment.** No result may attribute an effect to site without stating that.

### 5. Method choices, recorded because each could have gone the other way

- **Variance components, not η²** — see §2. Crossed random intercepts via an
  explicitly constructed `VCSpec`, **not** `vc_formula=`: the formula route
  returns `vcomp` in alphabetical name order regardless of the dict order passed
  in, which silently mislabels every component. This was caught only because a
  prototype reported "patient 0.52 dominant" when the value belonged to
  compartment.
- **`lbfgs` with a `powell` fallback**, method recorded per gene. Plain `lbfgs`
  converged on only 64% of genes and failed on 6 of 11 PVCA PCs. Where both
  converge they agree to ~3e-4 at the median with a log-likelihood gap of ~1e-7,
  so the failures are flat likelihood ridges, not a broken optimiser. **The null
  baseline is itself optimiser-sensitive**: under `lbfgs` alone the `dsp_run`
  permutation null came out at 0.23 rather than 0.00, which would have read as
  "batch is mostly a null artefact".
- **PVCA retains PCs to 0.60 cumulative variance** — Bushel's default, kept
  rather than tuned. 14 PCs here.
- **All 18,694 genes**, HVGs flagged rather than pre-filtered: the highest-
  variance genes are immunoglobulins, and letting feature selection choose the
  gene set would let it choose the answer.
- **`BC` controls retained** in the primary fit, matching the design table.
  Without them (S3) compartment falls 0.205 → 0.184 and the ordering is
  unchanged.
- **T4 uses Spearman on the 2000 HVGs, average linkage on 1 − ρ.** Over all
  18,694 genes every AOI pair correlates at ~0.99 and no structure is visible.
  Spearman because the top-variance immunoglobulins would otherwise drive a
  Pearson correlation. Average linkage (UPGMA) because Ward assumes a Euclidean
  geometry that a correlation distance does not have.

## Consequences

- Every Phase 2+ model carries `(1|patient)`. Changing that random-effects
  structure is a "stop and ask" decision (`CLAUDE.md`).
- Compartment is the primary fixed effect of interest; the lung-vs-brain
  contrast is estimated **within** compartment, never pooled across.
- `dsp_run` is carried as a covariate where estimable, and results for `mLN`,
  `TBME` and `BC` state that batch is inseparable from biology there.
- No p-value is reported without n, effect size and a confidence interval
  (hard constraint 7), and no variance share without its null.
- **`TIME-B` n = 8** remains the binding constraint. The lung-vs-brain immune
  contrast detects ~1.1–1.3 SD at 80% power (P0-T8); a null in Phase 2 is
  uninformative, not negative.
- PROJECT_PLAN §6's expectation that patients cluster together is recorded here
  as **measured false**, so a later reader does not reinstate it from the plan.
