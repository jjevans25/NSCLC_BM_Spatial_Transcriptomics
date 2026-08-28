# ADR 0012 — `dsp_run` is a reported sensitivity in Phase 2, not the primary model

**Date:** 2026-08-28
**Status:** Accepted
**Task:** P2-T3
**Builds on:** ADR 0009 (the model contract), ADR 0008 (detection and brain background)

## Context

`docs/NEXT_STEPS.md` flagged this as the one stop-and-ask Phase 2 could not
avoid, and the evidence genuinely pulls both ways:

- **For adjusting.** P1-T3 measured `dsp_run` at ~10% of variance (0.078 per-gene
  median, 0.096 by PVCA) against a permutation null of 0.002. It is real.
- **Against adjusting.** `TIME-B` splits **6 / 2** across the two DSP runs
  (`TIME-L` splits 9 / 6). At n = 2 in one cell the batch adjustment is barely
  identified, and the site estimate it produces leans on those two AOIs.

PROJECT_PLAN §6 specifies `score ~ site + (1|patient)` with no batch term.

## Decision

**Fit the plan's model as primary and the batch-adjusted model as a reported
sensitivity, per signature and per scoring method, with both fits' singularity
surfaced.** This is the shape P1-T3's sensitivity S2 used, and it is chosen for
the same reason: it puts a number on batch absorption instead of asserting that
batch does or does not matter.

```
primary      score ~ site + (1|patient_id)
sensitivity  score ~ site + dsp_run + (1|patient_id)
```

`(1|patient_id)` is mandatory in both, on non-independence grounds (ADR 0009 §1).
The config schema enforces it: both formula fields carry a regex requiring
`(1|patient_id)`, so a config edit that drops the random intercept fails at
Snakefile load rather than producing a silently wrong model.

`site` is releveled to `lung`, so the reported coefficient is **brain minus
lung**. Gate 2 expects a *negative* antigen-presentation estimate, and a
reversed reference level would inverse the gate's verdict, so the script asserts
the level order rather than trusting factor construction.

## What it measured

**Batch absorption is small.** Across all 24 fits, |sensitivity − primary| has a
median of **0.0076** and a maximum of **0.0465**. The direction, significance
and ordering of every signature are unchanged between the two models. This is
consistent with P1-T3's S2, which found dropping the term moved +0.014 into
compartment — real, and about one and a half percentage points.

**Half the fits are singular: 12 of 24.** `cytotoxicity`, `myeloid_m2` and `tls`
go singular under *both* models and both scoring methods — their patient
variance component collapses to zero, and the Satterthwaite df come back at 21
(= 23 − 2), i.e. the random intercept absorbs nothing. This is expected with 23
AOIs over 18 patients and little within-patient replication, and it is
**reported, never suppressed**: a singular fit is why the term is justified on
design rather than on estimated variance, which is precisely ADR 0009's
argument.

Singularity is *not* a reason to drop `(1|patient_id)`. Dropping it would treat
23 AOIs as 23 independent observations when 5 patients contribute two, and would
overstate the effective n. Changing the random-effects structure remains a
stop-and-ask.

## Consequences

- Every Phase 2 result table carries both models. The **primary** is the result
  of record; the sensitivity is reported beside it.
- BH-FDR is applied across signatures **within** each method × model family, not
  across models. The sensitivity refit is the same hypothesis under a different
  adjustment, not an additional one.
- Batch remains inseparable from biology in `mLN`, `TBME` and `BC` (Q3), which
  is unchanged by anything here — those compartments are not in Phase 2.
- The estimate is cross-checked against `statsmodels.MixedLM` for the Gate 2
  signature and agrees to ~5e-9 (ssGSEA) and ~2e-8 (z-score). Degrees of freedom
  and p-values are deliberately **not** compared: lmerTest uses Satterthwaite,
  statsmodels a Wald z with no finite-sample correction, so they differ by
  construction and the p-value of record is lmerTest's.
