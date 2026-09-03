# ADR 0018 — The primary restriction stands; a declared exploratory secondary carries the rest

**Date:** 2026-09-02
**Status:** Accepted
**Task:** P3-T3 (decided before any Phase 3 model was fitted)
**Decided by:** the project owner, after reading the P3-T2 audit
**Does not supersede ADR 0016.** §3's primary restriction is unchanged, word for
word. This ADR adds a second table beside it.
**Shape borrowed from:** ADR 0012 (primary + declared sensitivity, both reported)

## Context

### The pre-registration bit harder than anyone expected

ADR 0016 §3 pre-registered that a gene enters a primary fit only if it clears
the detection floor in **both** groups being compared. That was written before
`p3t2_checkpoint_detection` existed. Applied to the measured audit it leaves:

| contrast | genes | which |
|---|---|---|
| P3-T3 lung, `L` vs `TIME-L` | **2** | CD274, CD276 |
| P3-T3 brain, `LB` vs `TIME-B` | **2** | CD274, CD276 |
| P3-T4 tumour, `L` vs `LB` | **2** | CD274, CD276 |
| P3-T4 immune, `TIME-L` vs `TIME-B` | **4** | CD274, CD276, HAVCR2, VSIR |

A nine-gene panel produces a two-gene answer in three of four contrasts.

### The restriction is symmetric; the artefact it defends against is directional

Measured background per compartment (median `negprobe_log2`, P0-T6):

| contrast | background shift | what the artefact does |
|---|---|---|
| `L` 4.96 → `TIME-L` 5.36 | **+0.40** | pushes checkpoints to look **lower** in immune |
| `LB` 4.99 → `TIME-B` 5.30 | **+0.31** | same |
| `L` 4.96 → `LB` 4.99 | +0.03 | nothing |
| `TIME-L` 5.36 → `TIME-B` 5.30 | −0.06 | slightly **higher** in brain |

Two facts ADR 0008 could not have known at Gate 0, because the measurement did
not exist:

1. **The gradient is a P3-T3 problem, not a P3-T4 problem.** The two
   within-compartment *site* contrasts — the ones ADR 0008 point 5 was written
   about — have essentially no gradient (+0.03, −0.06). The large gradients are
   in the *compartment* contrasts, which is why ADR 0016 §3 extended the
   background sensitivity to the carrier model.
2. **Where the gradient is large it runs WITH the observed direction — the
   opposite of what was argued when the options were drafted.** These genes are
   observed strongly enriched in the immune compartment:

   | gene | `L` → `TIME-L` | `LB` → `TIME-B` |
   |---|---|---|
   | VSIR | 8/30 → **14/15** | 11/27 → **7/8** |
   | HAVCR2 | 2/30 → **10/15** | 2/27 → **7/8** |
   | CTLA4 | 2/30 → **10/15** | 2/27 → 2/8 |
   | TIGIT | 0/30 → **8/15** | 2/27 → 1/8 |

   **The case for Option B rested on the claim that the artefact could only
   shrink these effects, never manufacture them. That claim was wrong, and
   P3-T3 measured it.** It conflated two quantities that move in opposite
   directions:

   - **Detection rate.** Detection is `q3 > 2 × negprobe`, so higher background
     raises the bar and a gene is detected *less* often. Higher background →
     looks *depleted*.
   - **Expression.** The models are on `log2(q3 + 1)`, where background
     contributes additively to the observed signal. Higher background →
     near-background values are *inflated* → looks *enriched*.

   Background is higher in the immune compartments (median `negprobe_log2`:
   `L` 4.96 → `TIME-L` 5.36). The effects are all positive (enriched in immune).
   So for the models the gradient is **permissive**: it could have contributed
   to the very effects Option B would have promoted on the grounds that it
   could not.

   `carrier_background_sensitivity` measures it directly. Adjusting for
   `negprobe_log2` **shrinks** the estimate for 8 of the 9 gene × stratum cells
   the primary excludes, and leaves the well-detected genes alone or slightly
   larger:

   | gene × stratum | primary | background-adjusted | shift |
   |---|---|---|---|
   | IDO1, lung | +0.524 | +0.277 | **−0.247** |
   | PDCD1, brain | +0.245 | +0.017 | **−0.227** |
   | TIGIT, lung | +1.036 | +0.840 | −0.196 |
   | HAVCR2, lung | +0.963 | +0.780 | −0.182 |
   | PDCD1, lung | +0.408 | +0.233 | −0.176 |
   | HAVCR2, brain | +1.217 | +1.051 | −0.166 |
   | VSIR, lung | +1.336 | +1.222 | −0.114 |
   | CTLA4, lung | +1.425 | +1.321 | −0.104 |
   | CD274, lung *(primary)* | +0.612 | +0.635 | +0.023 |
   | CD276, lung *(primary)* | +0.147 | +0.201 | +0.054 |
   | CD276, brain *(primary)* | +0.159 | +0.245 | +0.086 |

   IDO1 loses its nominal significance under the adjustment (p 0.059 → 0.329);
   PDCD1 in lung goes 0.0008 → 0.029.

   The restriction excludes these genes because they fail the floor in the
   *tumour* compartment — which is the finding, not a defect. **This ADR was
   written before the fits existed and originally repeated the inverted
   argument; the paragraph above replaces it.** The decision does not change —
   it is strengthened, because the rejected alternative turns out to have rested
   on a premise the data contradicts.

Also measured, and worth recording because it rules out the obvious alternative
explanation: `TIME-B`'s panel deficit is not a blanket sensitivity loss.
Genome-wide detection is `TIME-L` 35.0% vs `TIME-B` 35.2% — indistinguishable —
while the checkpoint panel goes 8/9 → 4/9.

### P3-T2b changed the terms of the argument

ADR 0017 recovered the count scale and validated it against the raw `.dcc`
counts for 120/120 AOIs. In `TIME-B`, background is ~30 counts and
CTLA4/TIGIT/IDO1 carry 37–44. An AOI holds a median of 645 distinct values
across 18,694 genes.

**That is the same fact as the detection floor, stated as resolution rather than
as background** — and it is the version that survives an argument about where
the threshold should sit. It also showed, in the published analysis, exactly what
a checkpoint claim looks like when detection is not audited first: nominally
significant effects (p = 0.0028–0.0111) at 0.39–0.89 log2 on genes carrying
31–74 counts, none of which reach that paper's own stated |log2FC| > 1.5 bar.

## Decision

**Option C. The primary restriction stands exactly as pre-registered, and a
declared exploratory secondary carries the genes it excludes.**

Three options were put; the two rejected are recorded below because the reason
for rejecting them is the substance of this ADR.

### The primary is untouched

`checkpoint_models.tsv` fits only genes clearing the floor in **both** groups —
2 / 2 / 2 / 4 genes. It is the **only** Phase 3 table carrying a q-value, and
ADR 0016 §3's specification is unchanged. Genes that fail get a row with
`status = "not_assessable"`, their per-group detection counts, and **no
estimate**.

### The secondary is declared, separate, and unadjusted

`checkpoint_exploratory_models.tsv` fits every gene assessable in **at least
one** group, using the identical formulas. It carries:

- **No q-value. No FDR of any kind.** Not a stricter correction — none. A
  q-value is a claim about a family of hypotheses, and this is not a family; it
  is a set of estimates chosen because the primary could not make them.
- **`in_primary`**, marking rows whose gene also passes the primary. Those
  estimates **must equal the primary's** and the script asserts it to 1e-10.
  This is what makes the two tables comparable rather than rival, and it is
  cheap: it re-derives the primary through the secondary's code path.
- **The background gradient, per row**: `background_delta` (median
  `negprobe_log2` of the second group minus the first, on the same orientation
  as the estimate) and `gradient_direction`, one of:
  - **`protective`** — the gradient runs *against* the observed effect, so the
    artefact can only shrink it. Whatever is left is real or larger.
  - **`permissive`** — the gradient runs *with* the observed effect, so the
    artefact could contribute to it. Read with suspicion.
  - **`negligible`** — |background_delta| below 0.1 log2.
- **`status = "exploratory"`** on every row, and the per-group detection counts
  that ADR 0008 requires of any reported checkpoint.

### The label is exploratory-within-exploratory, and it is not decoration

A4 is already exploratory (ADR 0008). This table is a second step down from
that. **No P3-T7 sentence may rest on it**, no figure may present it beside a
q-value, and the reading rule is unchanged: a gene below the floor in a
compartment is "not assessable in `<compartment>`", never "lower in
`<compartment>`". The secondary reports what a fit *would* say; it does not
promote a gene to assessable.

## What was rejected, and why it matters

**Option A — keep ADR 0016 §3 and report nothing else.** Defensible, and it is
the safest thing a reader could ask for. Rejected because it discards the
phase's most interesting measurement: the compartment signals where the bias is
demonstrably protective, on genes that go 2/30 → 10/15. Answering §6's "which
compartment carries each checkpoint" for two genes out of nine, while holding
evidence about the other seven that the artefact cannot explain, is a worse
report, not a more honest one.

**Option B — supersede ADR 0016 §3 with a directional rule**, admitting a gene
to the *primary* when it is assessable in one group and the gradient runs
against the observed direction. Rejected, and this is the load-bearing
rejection. **It is a post-hoc threshold change.** A reader could not distinguish
"we relaxed it for a good reason" from "we relaxed it until something appeared",
and this project's central claim is that it does not do the second. Option C
obtains the same information at the cost of one extra table and gives up nothing
but the q-value, which the excluded genes were never entitled to.

**And it would have admitted the wrong genes.** Option B's rule was
"assessable in at least one group *and* the gradient runs against the observed
direction". Per §Context above, the gradient runs *with* the direction for every
one of these genes on the expression scale, so the rule as written admits
nothing — or, if applied with the sign error it was drafted under, admits
exactly the genes whose effects the background is partly producing. Either way
it fails, and the failure was invisible until the models ran. **That is the
argument for pre-registration in one example: the reasoning offered for relaxing
a threshold was checkable only after the threshold had done its job.**

The asymmetry is the point: Option C is reversible by a reader who disagrees
(ignore the second table), Option B is not (the pre-registration is gone).

## Consequences

- **Two tables to explain, and a reader who quotes the wrong one.** This is
  Option C's real cost and it is not fully removable. Mitigations: the tables
  are separate files, only one has a `q_bh` column, every secondary row carries
  `status = "exploratory"`, and P3-T5's figure and P3-T7's prose take their
  numbers from the primary. The `in_primary` assertion means the overlap cannot
  silently disagree.
- **The P3-T5 figure is unaffected.** It is a detection figure (ADR 0016 §5,
  coloured by `median_negprobe_ratio`); neither table feeds it.
- **P3-T7 must state the split explicitly** — that the primary answers for two
  genes, that the secondary exists, why the excluded genes were excluded, and
  that the exclusion is itself the finding. Alongside ADR 0016 §2's requirement
  to state that the glial compartment was not modelled.
- **P3-T4 inherits this decision** unchanged; it is the same rule applied to the
  site contrast, where the gradient is near zero and therefore mostly
  `negligible` rather than `protective`.
- **Gate 3's record moves to ADR 0019.**
