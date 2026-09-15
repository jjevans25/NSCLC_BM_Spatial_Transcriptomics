# ADR 0026 — The multiplicity family is larger than ADR 0024 declared, and P5-T5 corrects it once

**Date:** 2026-09-15
**Status:** Accepted
**Task:** P5-T3 (before any Phase 5 p-value was computed)
**Corrects:** ADR 0024 §8, which is otherwise unchanged and still Accepted
**Constrained by:** ADR 0008 (a cell below the detection floor is "not
assessable", never a null), ADR 0012 (primary + sensitivity, both reported),
ADR 0016 §3 (a threshold chosen after seeing the result is not a threshold),
ADR 0018 (a restricted primary with a declared exploratory secondary),
ADR 0025 (the censored patient without follow-up), CLAUDE.md (changing a
multiplicity correction is a stop-and-ask)

## Context

ADR 0024 §8 pre-registered the Phase 5 multiplicity family as:

> **6 signatures × 2 cohorts = 12 tests**, corrected with Benjamini–Hochberg at
> α = 0.05

**That arithmetic is wrong, and it is wrong against ADR 0024's own §2.**

§8 counted *cohorts* — this study and TCGA LUAD — and implicitly assumed one
test per signature in this study. But §2 of the same ADR splits this study into
**two arms**, lung (`TIME-L`) and brain (`TIME-B`), and says they are estimated
within site and **never pooled**. Two arms × six signatures is **12 tests in this
study alone**, before TCGA contributes anything.

So the declared family undercounts, and an undercounted denominator makes every
q-value **too small**. That is the direction that matters: it manufactures
significance rather than hiding it.

The error was caught at P5-T3, before a single p-value was computed, by counting
the cells the task was about to fit. Nothing has been corrected against the wrong
number.

## Decision

### 1. The primary family is the ASSESSABLE cells, not every cell

ADR 0018 established a restricted primary with a declared exploratory secondary,
and P3-T3 implemented it: a gene enters a primary fit only if it clears the
detection floor, and the rest get a `not_assessable` row carrying no estimate.
Phase 5 inherits that shape unchanged.

From `results/tables/signature_coverage.tsv` (P2-T2, computed long before this
ADR):

| arm | assessable | not assessable |
|---|---|---|
| lung (`TIME-L`) | **5** — antigen_presentation, cytotoxicity, exhaustion, myeloid_m2, tls | myeloid_m1 |
| brain (`TIME-B`) | **3** — antigen_presentation, cytotoxicity, myeloid_m2 | exhaustion, myeloid_m1, tls |

**This study contributes 8 primary tests**, not 6 and not 12. TCGA LUAD (P5-T4)
contributes its own count, established there — TCGA is bulk RNA-seq with no
`NegProbe` background, so the detection floor is a property of *this assay* and
does not transfer; all six signatures are expected to be measurable there, but
that is counted at P5-T4 rather than assumed here.

**Not-assessable cells stay out of the denominator.** ADR 0008 forbids reporting
them as prognostic nulls at all, so including them would put un-interpretable
rows in the family and *deflate* every q — correcting for tests whose results may
not be stated. They are reported with `status = not_assessable`, an explanatory
note, and no estimate, exactly as `checkpoint_carrier_models.tsv` does.

**Sensitivity fits never enter the family.** `zscore` is a declared sensitivity
to `ssgsea`'s primary (ADR 0012), and P2-T3 and P3-T3 both keep sensitivities out
of their FDR families. Running a second scoring method does not double the search;
it measures whether the primary is method-dependent.

### 2. P5-T5 corrects ONCE, over the complete family

ADR 0024 §8's substantive decision — that this study and TCGA are **one family**,
because "does the same signature stratify an independent cohort" is a single
question asked twice — is correct and stands. Only the count was wrong.

The consequence is procedural: **the correction cannot happen until both cohorts
have been fitted.** So it happens at P5-T5, once, over the union.

### 3. P5-T3 therefore ships NO q-column, and the absence is deliberate

`survival_km_models.tsv` carries `p_raw` and no `q_bh`. This is unusual for this
project — every other model table has both — and it is the point: a q-value in
P5-T3 would be a corrected value that nothing corrected, and correcting there and
again at P5-T5 would be two corrections over overlapping families.

The rule **asserts that no `q_bh` column exists in any P5-T3 output**, so the
convention cannot be restored by reflex.

Every P5-T3 p-value is reported as **raw and uncorrected, in those words**,
beside its n, hazard ratio and confidence interval (hard constraint 7).

### 4. What this does not change

Nothing about ADR 0024 §7. A null remains **uninformative, not negative**, and
the family size does not bear on that: at 12 and 7 patients with splits of 6/6
and 3/4, the analysis detects only very large hazard ratios whatever the
denominator. A correctly-counted family makes the q-values honest; it does not
make the study powered.

**A separation at 3 versus 4 patients describes this cohort and is never an
inferential claim.** No Phase 5 result may be a headline claim.

## Consequences

- P5-T5's acceptance criterion is now specific: it corrects over the union of
  P5-T3's 8 assessable cells and P5-T4's TCGA cells, states the family size it
  used, and asserts it against both upstream tables rather than quoting a number.
- `docs/limitations.md` and any P6-T7 prose must carry the corrected family, not
  ADR 0024 §8's 12.
- **No `config/config.yaml` edit.** `survival.multiplicity.n_signatures: 6` and
  `n_cohorts: 2` remain true as written — they describe the signatures and the
  cohorts, and neither was the family size. The family is *derived* from those
  plus the per-arm assessability, and deriving it in the rule is what makes it
  checkable against the data rather than hand-maintained. A config edit would
  also cost a ~40 min Phases 1–4 rebuild for a number that must be measured.
- This ADR is a correction, not a reversal: ADR 0024 §8's one-family decision is
  reaffirmed, and only its arithmetic is superseded.
