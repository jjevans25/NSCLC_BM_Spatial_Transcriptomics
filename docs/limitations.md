# Limitations

What this reanalysis cannot support, and why. Written at the end of Phase 0;
every entry cites the rule output that establishes it. Provenance and the
answers to Q1–Q5 live in `docs/data-provenance.md`.

---

## 1. Power — the binding constraint (P0-T8)

**The study can only detect very large effects in its central comparison.**

`TIME-L` (n = 15) vs `TIME-B` (**n = 8**), simulated on the design exactly as
`config/samples.tsv` records it — 16 distinct patients, two contributing two
`TIME-L` AOIs, and five appearing in both groups — with a patient random
intercept and significance judged against t on 14 df:

| Between-patient ICC | Minimum detectable effect at 80% power |
|---|---|
| 0.00 | 1.30 SD |
| 0.25 | 1.21 SD |
| 0.50 | 1.11 SD |
| 0.75 | 0.85 SD |

**In one sentence, as §6 asks:** with `TIME-B` n = 8, the lung-vs-brain immune
contrast can detect a standardised difference of roughly **1.1–1.3 SD** at 80%
power for plausible ICCs, so anything short of a very large effect will be
missed rather than disproved.

Two consequences that are easy to get backwards:

- **A null result here is uninformative, not negative.** At this floor, failing
  to detect a 0.5 SD difference is the expected outcome whether or not one
  exists. Phase 2 must report confidence intervals and never phrase an absence
  of significance as evidence of similarity.
- **The MDE improves as ICC rises**, which looks wrong and is not. The five
  patients contributing to both groups turn part of the comparison into a
  within-subject contrast, and the more variance sits between patients, the more
  that pairing buys. It is also why the paired-immune subset is a genuine
  consistency check — but with n = 5 it stays a check, never a headline.

Estimated with `statsmodels` MixedLM; Phase 2 fits `lme4`. REML and the df
approximation differ slightly, so treat these as close estimates.
Evidence: `results/tables/power.tsv`, `results/figures/power_curves.png`.

## 2. Batch is partly inseparable from biology (P0-T3)

GEO exposes no batch field; the DSP run id survives only in the DCC filenames.
It splits 91 / 29, and **`mLN`, `TBME` and `BC` each sit entirely inside one
run**. For those compartments a batch term is not estimable, so any difference
involving them is also a run difference and no model can separate the two.

The core `TIME-L` / `TIME-B` contrast spans both runs, so batch is estimable
there — but on only **2 `TIME-B` AOIs** in the smaller run, which is thin enough
that a batch-adjusted brain-immune estimate should be treated as a sensitivity
check rather than the headline. Aim A4, resting on `TBME`, inherits the
confounding in full.
Evidence: `results/tables/batch_crosstab.tsv`.

## 3. Compartment labels rest on placement, not sorting (P0-T6)

Q2 resolved that the AOIs are genuinely antibody-segmented, but **PanCK was the
only collection mask**. CD45 and GFAP guided where a pathologist placed each
ROI; nothing was sorted on them. A `TIME` AOI is the PanCK-negative segment of
an ROI sited in a CD45-rich region — not a CD45-sorted population.

The marker check supports the labels (4/4 gating criteria, effects 2–6 log2), so
this is a constraint on *wording*, not a doubt about the data: never write
"CD45+ AOI", and never treat a compartment label as a cell-type label. Cell-type
fractions from `SpatialDecon` describe the mixture inside a PanCK-negative
segment; they do not validate the label.
Evidence: `results/tables/marker_sanity_verdict.tsv`.

## 4. No probe-level QC and no limit of quantitation (P0-T5)

GEO carries the DCCs but no PKC, and PKCs are distributed behind Bruker
registration rather than from a pinnable URL. Without one, RTS probe ids cannot
be mapped to genes, so:

- **Probe-level outlier removal never happened on our side.** We inherit the
  submitters' summarisation as given, with no visibility into it.
- **The standard GeoMx LOQ is not computable.** It needs
  `geomean(negative probes) × GeoSD(negative probes)²`; the deposited matrix has
  a single `NegProbe-WTX` row, so there is no SD. The 2× background multiple in
  ADR 0007 is a stand-in for a statistic we cannot calculate.
- **Raw gene-level counts do not exist for us at all**, so normalisation could
  only ever be verified, never redone (Q1).

Evidence: ADR 0006, ADR 0007.

## 5. Background is not constant across compartments (P0-T6)

`NegProbe-WTX` rises monotonically from log2 4.96 in `L` to 5.67 in `BC` — about
1.6×. Q3 scaling divides each AOI by its own third quartile and so does not
remove a background difference that tracks tissue type. **A gene expressed near
background will look differentially expressed between lung and brain
compartments on background alone.** Phase 2 should either restrict claims to
genes comfortably above the negative probe or model the background level
explicitly; the P0-T6 effects (2–6 log2) are far above this and unaffected.

## 6. `TBME` detects far less than the tumour cores (P0-T5)

Median detection at 2× background: `L` 42.5%, `LB` 42.8%, `mLN` 40.2%,
`TIME-L` 35.0%, `TIME-B` 35.2%, but **`TBME` 21.4% and `BC` 20.4%**. This is
brain parenchyma's lower transcriptional complexity plus its higher background,
not a defect — but it means a global QC cutoff reaches the glial compartment
first, and four `TBME` AOIs are flagged (retained) at the ADR 0007 thresholds.
`TBME` is therefore the compartment carrying both the batch confounding (§2) and
the detection deficit.

## 7. Per-AOI cellularity is unavailable (Q4)

Nuclei count and surface area do not survive into GEO — `cell type` is the only
per-sample characteristics field, and there is no DSP annotation worksheet. So
models **cannot** be weighted by AOI area or nuclei count, and there is no way
to ask whether library size tracks cellularity. QC rests on detection rate and
raw read depth alone.

## 8. Clinical metadata constraints (Q5)

Survival analysis is possible — Supplementary Data 1 carries two time-to-event
columns and a censoring value — but with three constraints. **Age is banded**
(`40s`, `60s`, …), so it enters a model as an ordered factor, never as a
continuous covariate. **Missingness is inconsistently coded** (`N/A`, `NA`,
`n/a`, `Unspecified`). And the table covers all 44 cohort patients while GEO
carries **35**, so nine rows have no expression data.

## 10. Detection is the dominant limitation, ahead of power (P2-T2)

Sections 4–6 anticipated this; Phase 2 measured it, and the result is larger in
scope than "`TBME` detects less than the tumour cores". **Three of six immune
signatures failed their per-site detection-coverage floor**, scored on the 23
`TIME` AOIs with a gene counting as detected only above
`qc.detection_background_multiple` × that AOI's NegProbe-WTX level:

| set | `TIME-L` | `TIME-B` | |
|---|---|---|---|
| `antigen_presentation` | 13/13 | 13/13 | the only fully-detected set |
| `myeloid_m2` | 6/9 | 6/9 | assessable |
| `cytotoxicity` | 2/4 | 2/4 | assessable, **at the floor** (rests on GZMB, NKG7) |
| `tls` | 3/5 | **1/5** | not assessable in brain |
| `exhaustion` | 4/6 | **1/6** | not assessable in brain |
| `myeloid_m1` | **2/10** | **2/10** | not assessable **at either site** |

Two consequences that outlive Phase 2:

- **`exhaustion` and `tls` returned large, nominally significant "reductions in
  brain"** (q = 0.005 and 0.018) that **are not reportable**. Brain background is
  higher (§5), so a near-background gene reads as depleted in brain
  artefactually — ADR 0008's argument, now observed. Anywhere those sets appear
  the finding is "not assessable in brain", never "lower in brain", and the
  P2-T7 heatmap hatches those cells on the figure itself.
- **Secreted ligands and chemokines are systematically undetected** — CXCL10,
  CXCL11, IL1B, TNF, IL12B, CD80, CD86, NOS2, IL10, CCL22, CCL19 and CCL21 all
  sit below background, several at *both* sites. A ligand–receptor analysis
  needs the ligand side above background, which is why A5 was demoted at Gate 2
  (ADR 0014).

The detection gate is also load-bearing rather than decorative: it is what
stopped `myeloid_m1`'s apparent signature-versus-deconvolution contradiction —
signature down in brain, macrophage abundance up — being written up as a real
methodological conflict. Eight of its ten genes measure nothing, so the score
cannot disagree with anything (P2-T6, `docs/analysis-notes.md`).

**Phase 3's checkpoint panel overlaps `exhaustion` almost entirely** (`PDCD1`,
`CTLA4`, `LAG3`, `HAVCR2`, `TIGIT`), so the same outcome should be expected
there and designed for, not discovered.

## 9. Scope

- **`mLN` is out of scope** for the core comparison by design (§A.2), and is
  carried only so the design table is complete.
- **One AOI more than the source paper.** GEO deposits 120; the paper analyses
  119, with `TBME` = 19 rather than 20. Which `TBME` AOI the authors excluded,
  and why, is not recoverable from GEO — so a figure-for-figure replication of
  their `TBME` results is not possible.
- **This is a learning and methods project, not a publication** (CLAUDE.md).
