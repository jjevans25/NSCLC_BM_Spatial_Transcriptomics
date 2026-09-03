# ADR 0019 — The detection-floor slider is a sensitivity display, not a threshold control

**Date:** 2026-09-03
**Status:** Accepted
**Task:** P3-T6 (written before `notebooks/apps/checkpoint_explorer.py` existed)
**Constrained by:** ADR 0016 §1 (the floor is pre-registered at 0.5),
ADR 0007 (the 2× background rule), ADR 0008 point 4 ("not assessable in
`<compartment>`", never "lower in"), ADR 0010 (app-tier notebooks ship their
data in `public/`), ADR 0018 (detection and expression move opposite ways)
**Applies to:** every interactive control this project ever puts over a
pre-registered number, not just this slider

## Context

P3-T6 is an app-tier marimo notebook whose central control is a slider over the
detection floor. `CLAUDE.md`'s stop-and-ask list names "changing a QC threshold"
as a scientific decision, and ADR 0016 §1 pre-registered the floor at 0.5 in the
same commit that created the panel, before `p3t2_checkpoint_detection` existed.
P3-T2's entire reason for existing is that setting a floor after seeing which
genes clear it is the failure mode the phase is shaped to prevent.

A slider that moves that number is therefore either a serious problem or the
most useful thing in the phase, depending on what it is *for*, and nothing in
the repository currently says which. The framing lived in a session note.

The tension is real, not bureaucratic. The Phase 3 deliverable is an audit whose
headline is that **40 of 63 gene × compartment cells are not assessable and the
panel clears the floor in `TIME-L` alone**. That claim rests on one number. A
reader who cannot see how the claim responds to that number has to take the
0.5 on trust — and "we pre-registered it" is an answer about *process*, not
about whether the conclusion is fragile. The honest way to present near-LOQ data
is to show the response surface, not to assert a point on it.

## Decision

**The app exposes the gene-level floor as a slider whose default is the
pre-registered value, and no position of that slider produces a result.**

Four parts:

1. **The floor of record is 0.5 and does not move.** ADR 0016 §1 stands
   unamended. `config/config.yaml`'s
   `checkpoints.detection.detected_in_aoi_fraction` is untouched by this task,
   and every table, figure, model and sentence in Phase 3 continues to be
   computed at 0.5. The slider changes what is *drawn*, never what is
   *recorded*.

2. **The default reproduces `checkpoint_dotplot.png` exactly.** The app opens at
   0.5, all seven compartments, both sites — and at those settings its hatch
   set, ring set, row order, column order and printed counts must equal P3-T5's.
   This is the acceptance criterion PROJECT_PLAN §6 states, and it is what makes
   the slider a *view onto* the deliverable rather than a second, competing one.
   If the two disagree, one of them is wrong.

3. **The app says what it is, on the page.** A caption naming 0.5 as
   pre-registered (ADR 0016) sits beside the slider permanently, and the instant
   the slider leaves 0.5 a banner appears saying the view is a sensitivity
   display and not the pre-registered floor. A reader who screenshots the app at
   0.75 gets an image that labels itself.

4. **The background multiple is deliberately NOT exposed.** `detection_rate` in
   `checkpoint_detection.tsv` is already computed at
   `qc.detection_background_multiple` = 2.0 (ADR 0007), so a control over it
   could only relabel an axis, not recompute anything — a slider that appears to
   move a threshold and does not is worse than no slider. Recomputing detection
   inside the notebook to make it real is forbidden twice: `CLAUDE.md` hard
   constraint 3 (no reported number originates in a notebook) and the fact that
   `p3t2_checkpoint_detection` owns that rule and there must be exactly one.

## Why this is not the failure P3-T2 exists to prevent

The failure mode is **choosing** a threshold after seeing the results — picking
0.4 because it admits a gene you wanted. Three properties separate this from
that:

- The choice was already made, in public, in a commit that predates the data
  (ADR 0016). The slider cannot unmake it.
- The slider has no write path. It touches no config, no table, no model, no
  `q`-value. Every number in `checkpoint_carrier_models.tsv` and
  `checkpoint_shift_models.tsv` was fitted at 0.5 and stays fitted at 0.5.
- The direction of use is the opposite one. The reader is not being invited to
  find a floor at which the panel looks better; they are being shown that
  `TIME-L` stays readable and `TBME` stays unreadable *across* the range, which
  is an argument that the audit's conclusion does not depend on the exact value.
  Where the conclusion *does* depend on it — and for the marginal cells it does
  — the slider makes that visible too, which is the point.

**Sensitivity analysis and threshold selection are the same arithmetic and
opposite epistemics**, and the difference lives entirely in whether the result
was fixed first. This project already relies on that distinction elsewhere:
ADR 0012 reports a batch term as a sensitivity rather than a primary, and P3-T3
reports the background adjustment the same way. The slider is that shape, made
continuous and handed to the reader.

## Consequences

- **`checkpoint_explorer.py` recomputes `assessable` from the slider** rather
  than reading the column — necessarily, since the column is fixed at 0.5. At
  the default it must equal the column, and the notebook asserts that. This is
  the one place the app derives anything, and it derives a *display* flag, not a
  result.
- **No model output reaches the dot plot canvas.** Detection and expression move
  opposite ways under the same background gradient (ADR 0018), so the carrier
  and shift estimates live in their own section below a rule, with the reason
  stated. PROJECT_PLAN §6's "app over the P3-T3/T4 model outputs" is satisfied
  by adjacency, not by superposition.
- **The exploratory tables carry a standing no-claim banner.** They have no
  `q_bh` column and ADR 0018 forbids resting anything on them; an interactive
  table is exactly the artefact from which someone copies a number without its
  provenance.
- **Gate 3's record moves to ADR 0020.** `docs/NEXT_STEPS.md` reserved 0019 for
  it; this ADR takes 0019 because it must exist before the notebook does.
- **The precedent generalises.** Any future control over a pre-registered number
  — a coverage floor, an FDR α, a power assumption — inherits these four parts:
  default at the registered value, no write path, self-labelling off-default,
  and nothing exposed that cannot actually be recomputed from declared inputs.
