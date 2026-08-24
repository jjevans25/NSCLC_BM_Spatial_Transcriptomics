# ADR 0007 — P0-T5 QC thresholds

**Date:** 2026-08-24
**Status:** Accepted
**Task:** P0-T5
**Chosen in:** `notebooks/review/qc_review.py` (review tier; run at the P0-T5 gate)

## Context

`config/config.yaml` carried four `null` QC thresholds, deliberately: the schema
permits null so an unmade decision cannot masquerade as a made one. This ADR
makes them.

Three facts shaped the choice, all of them from the data rather than from
convention.

**1. There are no zeros.** Every one of the 18,694 × 120 values is strictly
positive, minimum 2.12 (Q1, P0-T4). "Detected" therefore cannot mean
"non-zero"; it must mean "above this AOI's own background". The standard GeoMx
limit of quantitation is `geomean(negative probes) × GeoSD(negative probes)²`,
but the deposited matrix contains a **single** `NegProbe-WTX` row — no set, no
SD, no geomean. That statistic is not computable here (ADR 0006).

**2. Detection is strongly compartment-dependent.** Median detection at 2×
background: `L` 42.5%, `LB` 42.8%, `mLN` 40.2%, `TIME-L` 35.0%, `TIME-B` 35.2%,
but **`TBME` 21.4% and `BC` 20.4%**. Brain parenchyma is less transcriptionally
complex than tumour *and* carries higher background (P0-T6 found `NegProbe`
rising monotonically from 4.96 log2 in `L` to 5.67 in `BC`). Both push detection
down for reasons that are tissue, not defect. A global cutoff therefore reaches
the glial compartment first, and a cutoff set to "catch the low tail" would
mostly be deleting `TBME`.

**3. Library size from the matrix measures nothing.** The matrix is already
Q3-normalised, so its column sums vary at only CV 13.4% by construction. Real
depth is in the DCC headers, where raw reads span 0.66M–36M.

**Sequencing saturation turned out not to discriminate at all.** The minimum
across 120 AOIs is 0.614 — every AOI clears the 0.50 floor the source paper
claims. Independently computing this and reproducing the paper's claim is worth
more as a validation of the DCC parsing than as a filter.

What the candidate cutoffs would flag was evaluated in
`notebooks/review/qc_review.py`, whose sliders show the counts move live:

| detection @2× | raw reads | flagged | composition | `TIME-B` left |
|---|---|---|---|---|
| <2.5% | <0.75M | 2 | TBME 2 | 8 |
| **<5%** | **<1.0M** | **5** | **TBME 4, LB 1** | **8** |
| <10% | <1.5M | 9 | TBME 8, LB 1 | 8 |
| <10% | <2.0M | 12 | TBME 8, L 2, LB 1, **TIME-B 1** | **7** |

## Decision

```yaml
qc:
  detection_background_multiple: 2.0
  min_gene_detection_rate: 0.05
  min_library_size: 1000000        # raw reads, DCC header
  min_sequencing_saturation: 0.50
  flag_only: true
```

**2× background** as the detection definition: a stand-in for the uncomputable
LOQ, and close to where a real GeoMx LOQ typically lands.

**5% detection and 1.0M raw reads.** This marks five AOIs:

| AOI | Reason | Detection | Raw reads |
|---|---|---|---|
| TBME13 | detection | 2.1% | 2.79M |
| TBME20 | detection | 3.5% | 8.76M |
| TBME17 | detection | 4.0% | 13.84M |
| TBME09 | read depth | 34.3% | **0.66M** |
| LB35 | read depth | 56.5% | **0.90M** |

The two axes catch genuinely different failures, which is the argument for
keeping both: TBME17 has 13.8M reads and still detects 4% of genes (deep
sequencing of a low-complexity, high-background sample), while TBME09 detects a
normal 34.3% off only 0.66M reads (an under-sequenced but otherwise sound AOI).
Neither rule alone finds both.

**0.50 saturation** is documentation, not a filter — it flags nothing and
records the floor the source paper asserts, so that a future re-run on different
data would catch a violation.

### Why conservative rather than moderate

The <10% option flags 8 of 20 `TBME`. Since the `TBME` detection deficit is a
property of brain parenchyma, flagging 40% of the compartment states a verdict
on the tissue rather than on individual AOIs, and would invite a later reader to
exclude them. 5% marks only AOIs that are extreme *even against other `TBME`*
— the `TBME` median is 21.4%, so 2.1–4.0% is an order of magnitude below its own
compartment, not merely below the tumour cores.

### Why not the strict option

It costs `TIME-B12`, and `TIME-B` n = 8 is the binding constraint on the whole
project. `TIME-B12` detects at 59.8%, the *highest* of any `TIME-B` AOI; it
would be flagged purely on 1.90M reads. Trading a scarce brain-immune AOI for a
depth rule that its own detection rate contradicts is a bad trade.

## Consequences

- **Nothing is dropped.** `flag_only: true`, so all five appear in
  `results/tables/qc_excluded.tsv` with `action = flagged_retained` and stay in
  every downstream analysis. The filename is the plan's; the column is the truth.
- **`TIME-B` remains n = 8**, so P0-T8's power calculation stands on the
  design's stated numbers and hard constraint 8's inline `n = 8` is unqualified.
- **`TBME` carries 4 of its 20 AOIs flagged**, and `TBME` was already the
  compartment wholly confounded with DSP run (Q3). Aim A4 inherits both. This
  belongs in `docs/limitations.md` at P0-T8.
- **These thresholds are not portable.** They are calibrated against a single
  `NegProbe` row on an already-Q3-normalised matrix. Re-deriving them is
  mandatory if a PKC ever allows real LOQ computation (ADR 0006).
- **The compartment skew is documented, not corrected.** A defensible
  alternative — flagging within compartment — was rejected as unconventional,
  but it remains the right answer if a later phase needs `TBME` judged on its
  own terms.
