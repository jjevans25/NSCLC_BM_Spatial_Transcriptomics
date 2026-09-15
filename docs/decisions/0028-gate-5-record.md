# ADR 0028 — Gate 5 record: the timebox was respected, and nothing survived the family

**Date:** 2026-09-15
**Status:** Accepted
**Task:** P5-T5 (Phase 5 complete)
**Gate:** GATE 5 — *"Timebox respected. The project is complete and coherent
without this phase; treat it as a bonus, not an obligation."*
**Verdict:** **PASSED**
**Records:** the outcome of ADR 0024 (pre-registration), ADR 0025 (the censored
patient), ADR 0026 (the corrected family), ADR 0027 (the external cohort)

## The gate

PROJECT_PLAN §6's Gate 5 does not turn on finding an effect. It turns on the
phase being bounded, and on P5-T1's acceptance criterion — *"a joined table **or**
a documented abandonment"* — being met honestly either way.

**It passed.** Phase 5 ran P5-T1 through P5-T5 in one session, produced every
output through a Snakemake rule, and stopped where it was designed to stop.

## The result

**0 of 14 pre-registered tests survive Benjamini–Hochberg at FDR 0.05.**
Smallest q = **0.2563**, on TCGA `antigen_presentation` (HR 0.70 [0.52, 0.94],
raw p = 0.018, n = 502 with 182 events).

| | |
|---|---|
| cells considered | **18** — 12 GeoMx signature × arm, 6 TCGA |
| cells tested | **14** — 8 GeoMx assessable + 6 TCGA |
| cells not assessable | **4** — `exhaustion` brain, `tls` brain, `myeloid_m1` both sites |
| survive FDR 0.05 | **0** |

**The one suggestive signal does not survive, and must not be reported as
though it did.** TCGA `antigen_presentation` is the strongest result in the
phase and it lands at q = 0.26. Higher antigen presentation associates with
better overall survival in the direction one would expect, at a cohort size that
can actually detect such a thing — and that is a *hypothesis*, not a finding.
P5-T5 exists precisely so that the raw p = 0.018 cannot be quoted alone.

## Three things this phase establishes

### 1. The GeoMx null is a power limit; the TCGA null is closer to evidence

These are **not the same kind of null**, and collapsing them would be the
phase's easiest error.

- **This study**: 12 lung and 7 brain patients enter a fit, splits 6/6 and 3/4,
  hazard-ratio intervals spanning eight- to thirty-fold. A null here is
  **uninformative, not negative** (ADR 0024 §7, fixed before the first curve).
- **TCGA**: 502 patients, 182 events, intervals roughly 0.5–1.2 wide. A null
  here is a **much stronger** statement — though still not proof of absence.

**The asymmetry is the deliverable**, not any single q-value. PROJECT_PLAN §6
anticipated exactly this: *"a negative here mostly tells you about your n, not
about the biology."*

### 2. The comparison is bounded, and the bound is structural

Per ADR 0027 §5, and it applies to every cross-cohort sentence:

- **Not the same measurement.** GeoMx `TIME-L` is the **PanCK-negative segment**
  of an ROI sited in a CD45-rich region — not a CD45-sorted population (Q2).
  TCGA is **whole bulk tumour**. A disagreement between the cohorts is not
  necessarily a disagreement about biology; an agreement is not necessarily a
  replication.
- **No brain comparator exists at all.** TCGA LUAD is primary lung. The brain
  arm — `TIME-B` **n = 8**, the binding constraint on the entire project —
  received **no external check**. P6-T6 must say so.

One concrete consequence: the two cohorts' `antigen_presentation` point
estimates point **opposite ways** (GeoMx HR 2.19 [0.65, 7.35]; TCGA HR 0.70
[0.52, 0.94]). **That is not a contradiction.** The GeoMx interval spans 1 and is
consistent with both directions, so it is one uninformative estimate beside one
informative one, on two different measurements. Writing it up as "the cohorts
disagree" would be wrong twice over.

### 3. The family was counted, and the first count was wrong

ADR 0024 §8 declared "6 signatures × 2 cohorts = 12". That forgot ADR 0024 §2's
own rule that this study's two arms are never pooled. ADR 0026 corrected it; the
true family is **14**, and P5-T5 **derives it by counting rows in the upstream
tables** rather than quoting any recorded number.

The direction of the error is what matters: an undercounted denominator makes
every q **too small**. Had the declared 12 been used, the smallest q would have
been computed against the wrong family — and the failure would have been
invisible, because a q-value carries no evidence of the denominator behind it.

**Two further guards are now standing**, both because Phase 4 shipped an FDR
step-down that was correct-looking and wrong under ties: the correction is
computed by `statsmodels` *and* an independent hand-rolled step-up (agreeing to
1.1e-16), and monotonicity is asserted directly. **This family contains two tied
pairs**, so the guard was not hypothetical.

## What Phase 5 may and may not be used for

**May:** the phase is reportable in full as a pre-registered, corrected, null
result with an honest accounting of what was and was not testable — 18 cells
considered, 14 tested, 4 not assessable, 0 surviving.

**May not:** no Phase 5 result may be a headline claim. Only **q** is
reportable, never a raw p from P5-T3 or P5-T4. A cell below Phase 2's coverage
floor is **"not assessable"**, never a prognostic null (ADR 0008) — a score built
from genes at background measures background, whatever it is regressed against.
Every brain claim states **`TIME-B` n = 8** inline (hard constraint 8), and every
p-value carries n, effect size and a confidence interval (hard constraint 7).

## Consequences

- **Phase 5 is complete.** `config.yaml → phases.survival` stays `true`;
  `TARGETS_SURVIVAL` carries all of P5-T1…T5.
- **Phase 6 is next and is unblocked.** PROJECT_PLAN calls it the primary
  deliverable, and P6-T1's blocker was resolved at ADR 0023.
- `resources/` now has **six** sanctioned writers and P6-T2's RO-Crate must
  describe all six: `raw/`, `reference/`, `supplementary/`, `ligand_receptor/`,
  `clinical/`, `tcga/`.
- P6-T6's limitations document inherits three items from this phase: the n = 16
  cohort with 12/7 fitted, the bulk-versus-compartment mismatch, and the absence
  of any brain comparator.
- P6-T7's write-up must carry the asymmetry in §1 above. Stating "no signature
  was prognostic" without distinguishing the two cohorts' power would convert a
  measurement of this study's size into a claim about biology — the same error
  ADR 0022 guards against for Phase 4.
