# ADR 0016 — Phase 3 pre-registration: the detection floor, the modelled compartments, and the background gradient

**Date:** 2026-09-02
**Status:** Accepted
**Task:** P3-T1 / P3-T2 (before any Phase 3 result was computed)
**Discharges:** ADR 0008 point 5
**Constrained by:** ADR 0009 §3 (batch is inseparable in three compartments),
ADR 0012 (primary + sensitivity, both reported)

## Context

Three thresholds decide most of what Phase 3 is allowed to say, and all three
have to be fixed before any result is looked at. PROJECT_PLAN §6 P3-T2 says
"set the floor in `config.yaml` before you look at the results"; §10's risk
register says "pre-register the floor"; `CLAUDE.md` says a threshold without an
ADR is a number someone made up. This is that ADR, and it was committed in the
same commit as the empty-panel fill, before `p3t2_checkpoint_detection` existed.

## Decision

### 1. The detection floor is 0.5 / 0.5 — Phase 2's values, not new ones

```yaml
checkpoints:
  detection:
    detected_in_aoi_fraction: 0.5    # gene present in a compartment
    min_assessable_fraction: 0.5     # panel assessable in a compartment
```

The per-value rule underneath is unchanged and is **not** re-derived: a gene is
detected in an AOI when `layers['q3']` exceeds
`qc.detection_background_multiple` (2.0) × that AOI's `NegProbe-WTX` level.
That is ADR 0007's definition, implemented in
`workflow/scripts/score_signatures.py`, and P3-T2 lifts it rather than writing a
second one. Two detection rules in one project is one too many.

`detected_in_aoi_fraction: 0.5` is chosen because it is **ADR 0008 point 4
verbatim**:

> A gene detected in fewer than half the `TIME-B` AOIs is reported as "not
> assessable in brain" — never as "lower in brain".

With `TIME-B` n = 8 that is 4 of 8. It is also `contexture.scoring`'s value, so
Phase 3 introduces no threshold Phase 2 did not already operate under, and the
two phases' coverage tables are directly comparable.

`min_assessable_fraction: 0.5` promotes the same rule to the panel: the panel is
assessable in a compartment when at least 5 of 9 genes are. It drives the
one-sentence summary in P3-T7 and nothing else; every model and every dot on the
P3-T5 figure is gated at gene level.

**Rejected: a stricter 0.6 gene floor** (5 of 8 in `TIME-B`). More conservative
against the background gradient, but it is a number with no precedent anywhere
in the project, and inventing a threshold at the moment it starts to matter is
what pre-registration is meant to prevent.

**Rejected: a gene-level floor only, with no panel-level fraction.** Closer to
ADR 0008, which only ever legislates per gene. Declined because P3-T5's caption
and P3-T7's prose then have no pre-registered basis for the sentence "the panel
is not assessable in brain", which is the phase's most likely headline.

### 2. Only four compartments are modelled; all seven are audited

```yaml
  audit_aoi_codes: ["L", "LB", "mLN", "TBME", "TIME-L", "TIME-B", "BC"]
  model_aoi_codes: ["L", "TIME-L", "LB", "TIME-B"]
```

From `results/tables/batch_crosstab.tsv` (P0-T3, Q3):

| | run A | run B | |
|---|---|---|---|
| `L` | 22 | 8 | spans both |
| `LB` | 21 | 6 | spans both |
| `TIME-L` | 9 | 6 | spans both |
| `TIME-B` | 6 | 2 | spans both, barely |
| **`TBME`** | **20** | **0** | **run A only** |
| `mLN` | 13 | 0 | run A only |
| `BC` | 0 | 7 | run B only |

The four modelled compartments all span both DSP runs, so **no Phase 3 contrast
confounds batch with biology**. Within lung, `L` vs `TIME-L` is clean; within
brain, `LB` vs `TIME-B` is clean.

`TBME` is excluded from the models. ADR 0009 §3 is explicit that batch is not
separable there by any covariate, and `TBME` is additionally the
lowest-detecting compartment in the project (21.4%), carries the second-highest
median background, and holds 4 of the project's 5 QC-flagged AOIs. A
glial-compartment checkpoint estimate would be a number whose confound cannot be
quantified, let alone removed. `mLN` is the only compartment at its site and so
yields no within-site contrast at all; `BC` is non-tumour.

**This narrows §6's P3-T3 from "tumour vs. immune vs. glial" to tumour vs.
immune, and that is a real loss.** It is recorded here rather than left implicit
so that a reader does not conclude the glial compartment was overlooked. P3-T7
must state it. The alternative — modelling `TBME` with an inline caveat — was
rejected on the grounds that ADR 0014 used to demote A5: a caveat that a number
cannot be interpreted is not a substitute for not producing the number, when the
number will be read off a figure regardless.

All seven compartments still appear in the P3-T2 audit and on the P3-T5 figure.
The audit is the phase's deliverable under Gate 3, so it must cover the
compartments no model can touch.

### 3. The background gradient is addressed by BOTH routes ADR 0008 offers

ADR 0008 point 5:

> Any brain-vs-lung checkpoint comparison must address the background gradient
> explicitly — by restricting to genes comfortably above background, or by
> modelling the negative-probe level — and say which it did.

Phase 3 does both, and says so here.

- **Restriction is the primary.** A gene enters a primary fit only if it clears
  the floor in *both* groups being compared. A gene that fails gets a row with
  `status = "not_assessable"` and no estimate — present in the table, never
  silently absent.
- **Modelling is a declared sensitivity.** `negprobe_log2` enters as a covariate
  in a refit whose estimate shift from the primary is reported, exactly the
  shape ADR 0012 established for `dsp_run`. This quantifies the gradient rather
  than asserting it away.

Both sensitivities (`dsp_run`, `negprobe_log2`) are run for both questions, and
all of it is reported. Nothing here is a tiebreaker for anything else.

**The background sensitivity is applied to the carrier model (P3-T3) as well as
the shift model (P3-T4), which goes beyond ADR 0008 point 5.** Point 5 mandates
it only for brain-vs-lung, but the measured gradient between `L` (median log2
4.96) and `TIME-L` (5.36) is *larger* than the one between `TIME-L` (5.36) and
`TIME-B` (5.30). The compartment contrast is therefore at least as exposed to
the artefact as the site contrast, and treating only the site contrast would
leave the phase's core claim unguarded.

**Rejected: `negprobe_log2` in the primary model.** It would use every gene
rather than only those above the floor, but it changes the primary
specification, and for a gene detected in 1 of 8 `TIME-B` AOIs no covariate can
rescue a quantity that was not measured — it would manufacture estimates for
exactly the genes ADR 0008 says are not assessable.

### 4. Multiplicity: two families, not one and not three

BH-FDR across genes, **within site, within model** for P3-T3; across genes,
**within compartment, within model** for P3-T4. T3 and T4 answer different
questions from overlapping AOIs; correcting across both at once would create a
third family nobody asked for, and correcting across the sensitivity refits
would treat the same hypothesis under a different adjustment as a new one
(ADR 0012).

### 5. The P3-T5 figure colours by the background-normalised ratio

*(Added at P3-T2, 2026-09-02, before the figure existed.)*

PROJECT_PLAN §6 P3-T5 says the dot plot's "colour = mean expression". It will
instead colour by **`median_negprobe_ratio`** — each gene's value divided by its
own AOI's NegProbe-WTX level, median across the compartment.

The reason is §3's artefact arriving at the deliverable. Q3 scaling does not
remove a background difference that tracks tissue type, so on a mean-log2 colour
scale an *undetected* gene still shows its background level — and because brain
background is higher, undetected genes would render **brighter in brain than in
lung**. A dot plot is exactly the artefact someone reads a compartment
difference off a colour from, which is the reasoning that made P2-T7 hatch its
not-assessable cells on the figure rather than in the caption.

The background-normalised ratio is immune to the gradient by construction, and
it is the statistic ADR 0008's own reconnaissance table reported, so the figure
speaks the same units as the ADR that governs the phase. Both `mean_log2_all`
and `mean_log2_detected` are still carried in
`results/tables/checkpoint_detection.tsv`; the deviation from §6's literal
wording is stated in the figure caption.

`mean_log2_detected` was rejected as the primary for a separate reason: its
denominator changes per cell, so two dots of the same colour could rest on 15
AOIs and on 4.

## Consequences

- `config/config.yaml` gains a `checkpoints:` block and `phases.checkpoints`
  goes true, in a single commit — that file is a declared input of
  `p0t7_assemble_h5ad`, so every edit rebuilds the `.h5ad` and re-runs Phases 1
  and 2. Phase 3 pays that cost exactly once, here.
- `workflow/schemas/config.schema.yaml` carries the `(1|patient_id)` regex guard
  on all six Phase 3 formula fields, so a config edit cannot silently drop the
  random intercept. Thirteen negative tests were run against the two schemas at
  P3-T1 and all thirteen were rejected as intended; the real files still
  validate.
- Changing any threshold in this ADR remains a stop-and-ask, and a result that
  makes one look wrong is a conversation, not an edit.
- **The likely outcome of Phase 3 under these choices is that most of the panel
  is not assessable in brain.** That is the anticipated result, not a failure of
  the design: the panel overlaps `exhaustion`, which cleared 1 of 6 genes in
  `TIME-B` (ADR 0014). Gate 3 licenses it explicitly. What these thresholds buy
  is the right to say it as a measured assay-sensitivity limit rather than as
  biology.
