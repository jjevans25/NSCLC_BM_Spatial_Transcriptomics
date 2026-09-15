# ADR 0025 — A censored patient with no follow-up time is excluded from fits, never imputed

**Date:** 2026-09-15
**Status:** Accepted
**Task:** P5-T2 (before any patient-level score was aggregated)
**Extends:** ADR 0024 (Phase 5 pre-registration), which did not anticipate this
**Constrained by:** CLAUDE.md ("flag, don't silently drop"; dropping a case is a
stop-and-ask), ADR 0008 (a near-background gene is "not assessable", never
"lower"), ADR 0012 (primary + sensitivity, both reported), ADR 0022 (an honest
null must be reachable rather than escapable)

## Context

ADR 0024 §3 fixed the censoring convention: the lung endpoint column spells
censoring as the literal `Alive`, the brain column as an **empty cell**, and the
two sets are identical — patients **6, 11 and 35**. P5-T1 asserted that identity
and it holds.

What ADR 0024 did **not** anticipate is the consequence. **Neither spelling
carries a duration.** `Alive` is a status, not a time, and the brain cell is
blank. Supplementary Data 1 has no last-contact date, no date of last follow-up
and no enrolment date — the full column list is in `docs/data-provenance.md` §Q5
and P5-T1 asserted it at 15 columns. So for a censored patient **this project
holds an event status and no time axis position at all.**

That is not a missing value to be handled. It is the absence of the second half
of a time-to-event observation.

**One of the three is in the Phase 5 cohort: patient 35, and it is in BOTH
arms** — `TIME-L` and `TIME-B`. So this is not a corner case. It is 1 of 13 lung
patients and 1 of 8 brain ones, and at this n one patient is between 8% and 13%
of an arm.

## Decision

### 1. Excluded from every fit, recorded with a reason, never imputed

A patient with an event status and no follow-up time does not enter a
Kaplan–Meier fit, a log-rank test, or any median split. The exclusion is
**recorded, not silent**: `results/tables/survival_excluded.tsv` names the
patient, the arm, the reason and the action, in the shape
`results/tables/qc_excluded.tsv` established at P0-T6.

It is a **separate file**, not an append to `qc_excluded.tsv`. That table is
P0-T6's output and one file must have one writer; a Phase 5 rule writing into a
Phase 0 table would make its provenance unreadable.

The patient is **kept in `survival_patient_scores.tsv`** with `has_followup =
False` and `in_primary_fit = False`. Flag-don't-drop applies here exactly as it
does to an AOI: the row exists, carries its scores, and says why it is not
fitted. A reader counting rows must be able to see the excluded case.

### 2. Administrative censoring at the cohort maximum is REJECTED

The standard alternative is to censor administratively at the last observed time
in the cohort — 143 months for lung, 100 for brain. Rejected, for two reasons,
and the second is the one that decides it:

- **It invents an observation.** This cohort has no last-contact date, so the
  imputed duration would not be a coarsened measurement of anything. It would be
  a number chosen by this project and then analysed as though the authors had
  recorded it.
- **It puts the invented value exactly where the estimator is most sensitive.**
  Censoring at the maximum makes patient 35 the longest-followed patient in both
  arms **by construction**. A median split at n = 13 and n = 8 turns on the
  ordering of the tail, so the fabricated point would be the single most
  influential observation in the analysis. That is the opposite of conservative.

Imputing a *shorter* time would be worse still — it invents an early loss to
follow-up that nothing supports.

### 3. The consequence, stated before any fit

| arm | cohort | enters a fit | median split |
|---|---|---|---|
| lung (`TIME-L`) | 13 | **12** | 6 / 6 |
| brain (`TIME-B`) | **8** | **7** | 3 / 4 |

**The brain arm lands exactly on `survival.model.min_arm_size: 3`.** It clears
the floor and does not have a patient to spare: one further exclusion, for any
reason, drops the smaller arm below the pre-registered minimum and the brain
analysis becomes not assessable rather than null. Recorded here so P5-T3 finds
this written down rather than discovering it.

**`TIME-B` n = 8 remains what every brain claim states inline** (hard constraint
8). Seven is the number entering a fit; eight is the cohort. Both belong in any
reported sentence, and reporting only the larger would overstate the analysis
while reporting only the smaller would understate the design.

### 4. This strengthens ADR 0024 §7 rather than qualifying it

ADR 0024 §7 fixed, before any curve, that a Phase 5 null is **uninformative, not
negative**. Twelve events over 12 lung patients and 7 over 7 brain ones makes
that more true, not less: every patient entering a fit has an event, so the
arms carry no censoring at all and the analysis is a comparison of ordered
survival times at n = 6/6 and n = 3/4.

**A separation at 3 versus 4 patients describes this cohort and is never an
inferential claim.** No Phase 5 result may be a headline claim.

## Consequences

- `survival_excluded.tsv` is a new Phase 5 output and P6-T2's RO-Crate must
  describe it alongside the other results tables.
- **No `config/config.yaml` edit.** This is a structural fact about the source
  data, not a tunable threshold — there is no value to set, because the choice is
  between recording an absence and fabricating a number. CLAUDE.md's "every
  threshold belongs in a sibling YAML" does not apply, and an edit would cost a
  ~40 min Phases 1–4 rebuild for a key nothing would ever vary. The arm-size
  floor that *is* tunable, `survival.model.min_arm_size`, already exists.
- P5-T3 fits complete cases and states both n's. P5-T5's multiplicity family is
  unchanged at 6 signatures × 2 cohorts.
- If a later cohort revision supplies follow-up times, this ADR is superseded
  rather than amended — the exclusion is a consequence of what the file contains,
  so a different file is a different decision.
