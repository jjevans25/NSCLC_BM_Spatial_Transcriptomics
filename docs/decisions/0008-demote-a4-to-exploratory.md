# ADR 0008 — Aim A4 is demoted to exploratory at Gate 0

**Date:** 2026-08-24
**Status:** Accepted
**Task:** Gate 0
**Amends:** the aims table in `Markdowns/PROJECT_PLAN.md` §1.2

## Context

Gate 0 passed on all three conditions. This ADR records the one scope change
made at that gate.

**A4** is "map the checkpoint landscape by **compartment × site** — which cell
compartment carries each checkpoint gene, and whether that shifts brain vs
lung" (Phase 3). Two Phase 0 findings undercut it specifically.

**1. Most of the checkpoint panel sits at background, especially in the brain.**
Ratio of each gene to its own AOI's `NegProbe-WTX`, median across the immune
compartments, with the count of AOIs clearing 2× background:

| Gene | TIME-L ratio | TIME-B ratio | detected, TIME-L | detected, **TIME-B** |
|---|---|---|---|---|
| PDCD1 | 2.03 | 1.90 | 8/15 | **3/8** |
| CD274 | 2.89 | 2.32 | 13/15 | 5/8 |
| CTLA4 | 2.76 | 1.78 | 10/15 | **2/8** |
| LAG3 | 1.97 | 1.74 | 6/15 | 3/8 |
| HAVCR2 | 2.20 | 2.34 | 10/15 | 7/8 |
| TIGIT | 2.16 | 1.41 | 8/15 | **1/8** |
| IDO1 | 2.20 | 1.23 | 8/15 | **2/8** |
| PDCD1LG2 | 1.54 | 1.74 | 3/15 | 1/8 |
| BTLA | 0.76 | 0.88 | 0/15 | 0/8 |
| VSIR | 3.99 | 3.74 | 14/15 | 7/8 |

(Reconnaissance against `results/interim/aoi_normalised.h5ad`, using canonical
checkpoint genes. `config/checkpoints.yaml` is still empty — populating it is a
Phase 3 "stop and ask" and is not done here.)

For CTLA4, TIGIT and IDO1 there is essentially no brain-side measurement to
compare against: 1–2 of 8 `TIME-B` AOIs carry the gene above background. LAG3,
PDCD1LG2 and BTLA are at or below background at both sites.

**2. The background gradient biases the answer in A4's exact direction.**
P0-T6 found `NegProbe-WTX` rising monotonically from log2 4.96 in `L` to 5.67
in `BC` — brain background is higher. A gene sitting near background therefore
has an inflated denominator in brain and will read as **depleted in brain**
whether or not it is. "Reduced checkpoint expression in the brain TIME" is
precisely the headline A4 is shaped to produce, and it is the finding this
artefact manufactures.

Add the P0-T8 floor — 1.1–1.3 SD detectable at 80% power with `TIME-B` n = 8 —
and A4 is being asked to detect modest shifts, in genes that are mostly
undetected, against a bias pointing the same way as the expected result.

Nothing here is a defect in the data or the pipeline. It is the assay's
sensitivity meeting a low-abundance gene family.

## Decision

**A4 is exploratory, not a headline aim.** Concretely, and this is what the
label has to mean or it is only a word:

1. **Phase 3 still runs.** The rules, tables and figures are still built. The
   demotion is about what may be *claimed*, not about what is computed.
2. **No A4 result may be a headline finding**, in the report or in any summary.
3. **Every checkpoint gene reported states its detection**: how many AOIs
   exceeded background, per site, alongside any effect size. A gene's expression
   value is not reportable without it.
4. **A gene detected in fewer than half the `TIME-B` AOIs is reported as "not
   assessable in brain"** — never as "lower in brain". On the table above that
   is PDCD1, CTLA4, LAG3, TIGIT, IDO1 and PDCD1LG2, i.e. most of the panel.
5. **Any brain-vs-lung checkpoint comparison must address the background
   gradient explicitly** — by restricting to genes comfortably above background,
   or by modelling the negative-probe level — and say which it did.

**A5 is flagged, not demoted.** A5 (inferred crosstalk across adjacent
compartments, Phase 4) inherits `TBME`'s problems in full: `TBME` is wholly
inside one DSP run, so batch and biology are inseparable there, and it has the
lowest detection of any compartment at 21.4%. That is a real concern, but its
evidence arrives in Phase 2. **Re-decide A5 at Gate 2**, when there is something
to decide it on, rather than three phases early.

## Consequences

- **Phase 3's purpose changes from claim to reconnaissance.** Its honest output
  is "here is which checkpoints are measurable at all in this assay, and here is
  what that rules out", which is a genuinely useful negative result for a
  methods project — and a better one than a fragile positive.
- **The project's centre of gravity moves to A1–A3.** Those rest on
  compartments spanning both DSP runs, detecting normally, with validated
  labels. A6 (survival) is unaffected by this reasoning.
- **A "no checkpoint differences" result is now uninterpretable as biology**,
  and must be stated as an assay-sensitivity limit. Writing it up as a
  biological finding would be the specific error this ADR exists to prevent.
- **`config/checkpoints.yaml` is still a Phase 3 decision.** Whichever genes go
  in it, the detection reporting in point 3 applies. The table above is
  reconnaissance and confers no licence to populate the panel without asking.
- **Reversible.** If Phase 2 shows the background gradient is tractable and the
  detected subset (VSIR, CD274, HAVCR2) supports a real contrast, A4 can be
  restored by a superseding ADR. The demotion constrains claims, not effort.
