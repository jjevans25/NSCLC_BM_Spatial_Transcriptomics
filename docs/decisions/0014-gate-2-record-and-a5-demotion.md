# ADR 0014 — Gate 2 record: direction reproduced, A5 demoted to exploratory

**Date:** 2026-08-28
**Status:** Accepted
**Task:** P2-T6 (records Phase 2)
**Gate:** this is the Gate 2 record
**Follows:** ADR 0008, which flagged A5 and deferred the decision to here

## 1. Gate 2's question, and the answer

> Do you recover the published direction (reduced antigen presentation and
> B/T function in brain)? **If yes:** pipeline validated. **If no:** stop and
> debug — an unexpected result at this n is far more likely to be a pipeline
> bug than new biology.

**Yes, with one arm unmeasurable rather than contradicted.** All estimates are
brain minus lung, within the `TIME` compartment, `score ~ site + (1|patient_id)`,
`TIME-B` n = 8.

| | coverage | estimate (z-score) | q | verdict |
|---|---|---|---|---|
| antigen presentation | **13/13 both sites** | −0.477 SD [−1.030, +0.075] | 0.098 | direction recovered, underpowered |
| cytotoxicity (T/NK) | 2/4 both sites | −0.946 SD [−1.590, −0.302] | **0.018** | recovered, significant |
| B-cell / TLS | **1/5 in brain** | — | — | **not assessable in brain** |

Antigen presentation is negative under both scoring methods and is the only set
fully detected at both sites — the set the gate turns on is also the cleanest
one. It does not clear FDR 0.05, and that is the predicted outcome rather than a
disappointment: −0.48 SD sits well below the **1.1–1.3 SD** floor P0-T8
established at this n, so the non-significance is uninformative, not negative.

Nothing surprising appeared, so the gate's debug branch was not triggered.
Supporting evidence: direction agrees with the 5-patient paired check in **12 of
12** signature × method cells; lme4 and `statsmodels.MixedLM` agree on the gate
signature to 5e-9; the two scoring methods correlate at ρ 0.79–0.97; and P2-T6's
one testable convergence — cytotoxicity — has the signature and the safeTME
deconvolution (CD8+NK 0.139 → 0.081) pointing the same way.

## 2. The finding the gate question did not ask about

**Three of six signatures failed their detection-coverage floor**, and two of
them show large, nominally significant "reductions in brain" that are exactly the
ADR 0008 artefact:

| set | `TIME-L` | `TIME-B` | |
|---|---|---|---|
| `exhaustion` | 4/6 | **1/6** | q = 0.005 — **not reportable** |
| `tls` | 3/5 | **1/5** | q = 0.018 — **not reportable** |
| `myeloid_m1` | **2/10** | **2/10** | not assessable at either site |

ADR 0008 predicted this for the checkpoint panel; Phase 2 shows it reaches
further. **Only antigen presentation, cytotoxicity and myeloid M2 are
interpretable.** For the other three the finding is "not assessable", never
"lower in brain".

The coverage gate demonstrably earned its keep: `myeloid_m1`'s signature falls in
brain while macrophage abundance rises, which reads as a methodological
contradiction between signatures and deconvolution. It is not one — 8 of its 10
genes are below background at both sites, so the score measures nothing. Without
the gate that would have been written up as a real conflict.

## 3. A5 is demoted to exploratory

ADR 0008 flagged A5 (inferred crosstalk between adjacent compartments, Phase 4)
rather than demoting it, on the grounds that the evidence would arrive in
Phase 2. It has, and it is worse than that ADR anticipated.

**A ligand–receptor analysis needs the ligand above background. Phase 2 measured
the ligand side directly, at both sites, and it is largely absent.** The genes
below background in `myeloid_m1` are CXCL10, CXCL11, IL1B, TNF, IL12B, CD80,
CD86 and NOS2 — 8 of 10, **at both sites, not just brain**. Add IL10 and CCL22
(`myeloid_m2`) and CCL19 and CCL21 (`tls`). Secreted ligands and chemokines are
systematically undetected in this assay, which is the class of gene A5 depends on
most.

The `TBME` compartment, one side of the brain adjacency A5 exists to ask about,
compounds it: **21.4% detection** (lowest of any compartment), median NegProbe
**48.2** (second highest), **20 of 20 AOIs in DSP run A** so batch is inseparable
from biology there (Q3), and 4 of the project's 5 QC-flagged AOIs.

**Decision: A5 is exploratory, on the same terms ADR 0008 set for A4.** Phase 4
runs — PROJECT_PLAN's Gate 4 already frames "no LR pair exceeded chance
expectation" as an honest and useful result, and the methods value of doing the
analysis properly is real. But **no A5 result may be a headline claim**, every
nomination states the detection status of both partners, and the anticipated
null is reported as an **assay-sensitivity limit**, never as evidence that the
crosstalk is absent.

Rejected: dropping Phase 4 entirely (loses the honest-null write-up Gate 4
envisages, and the methods demonstration) and restricting A5 to the lung
`L ↔ TIME-L` adjacency (better powered and batch-estimable, but it answers a
question this project is not about — A5 exists to ask about the brain metastasis).

## 4. Two defects in PROJECT_PLAN §6, recorded because the plan is gitignored

Found by the `/gate 2` audit. `Markdowns/PROJECT_PLAN.md` is not version
controlled, so this ADR is the durable record — the same reason ADR 0008 exists.

- **P2-T7's acceptance criterion is unfalsifiable.** "in Snakemake report,
  publication-grade, colourblind-safe" — the first and third are checkable and
  were checked; **"publication-grade" has no operational definition** and cannot
  be met or missed. §0 says a criterion that cannot be checked is not a
  criterion. It needs an operational form (renders in `snakemake --report`;
  caption states n and the power floor; no clipped text at 200 dpi) or removal.
- **The A5 re-decision was owned by nothing.** ADR 0008 mandates it "at Gate 2",
  but it appears in no Phase 2 task row and no Accept condition — it lived only
  in the aims table and `NEXT_STEPS.md`. A gate requirement that no task owns is
  one that gets skipped. Resolved here; the general fix is that any obligation
  an ADR defers to a later gate must land in that phase's task table.

## Consequences

- Phase 2 is complete; Phase 3 (A4, checkpoint landscape, already exploratory
  per ADR 0008) is next.
- **The detection ceiling is now the project's dominant limitation**, ahead of
  `TIME-B` n = 8. Phase 3 should expect it, and `docs/limitations.md` should say
  so — Phase 2 turned it from a checkpoint-panel worry into a measured,
  assay-wide property.
- Every Phase 2 claim carries its coverage status. The three not-assessable sets
  are reported as such wherever they appear, including on the P2-T7 heatmap,
  where they are hatched.
