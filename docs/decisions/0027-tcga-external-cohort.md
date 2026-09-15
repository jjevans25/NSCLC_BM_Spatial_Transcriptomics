# ADR 0027 — The external cohort is UCSC Xena's TCGA LUAD, and the family is 14

**Date:** 2026-09-15
**Status:** Accepted
**Task:** P5-T4 (before any TCGA number was computed)
**Discharges:** ADR 0026 §1, which deferred TCGA's contribution to the
multiplicity family to be *measured here*
**Exercises:** ADR 0024 §1's pre-authorised fallback from cBioPortal to UCSC Xena
**Constrained by:** ADR 0005 (acquisition provenance, `resources/` read-only),
ADR 0008 (a cell below the detection floor is "not assessable", never a null),
ADR 0012 (primary + sensitivity, both reported), ADR 0024 §7 (a null is
uninformative, not negative), ADR 0025 (a case without follow-up is excluded,
never imputed), ADR 0026 (P5-T5 corrects once, over the whole family)

## Context

P5-T3 found **no immune signature stratified survival** in this cohort: across
eight assessable signature × arm cells, every hazard-ratio interval spans 1, at
lung n = 12 (6/6) and brain `TIME-B` n = 8 with 7 fitted (3/4).

PROJECT_PLAN §6 P5-T4 exists for exactly that outcome:

> Does the same signature stratify a properly powered independent cohort? **A
> negative here mostly tells you about your n, not about the biology.**

At ~500 patients TCGA LUAD can separate a power limit from an absent effect,
which is the one thing this project's own n cannot do. **A null in TCGA means
something a null here does not.**

## Decision

### 1. The source is UCSC Xena, because cBioPortal's bulk files are gone

ADR 0024 §1 named cBioPortal, because PMID 36216799's own TCGA analyses came from
there and a comparator should share provenance with the claim it contextualises.
It also pre-authorised the fallback in terms: *"UCSC Xena is the fallback if the
datahub tarball proves awkward; that substitution changes `config/tcga.yaml` only
and is free."*

It has proved awkward. Measured on 2026-09-15:

| URL | result |
|---|---|
| `cbioportal-datahub.s3.amazonaws.com/luad_tcga_pan_can_atlas_2018.tar.gz` | **HTTP 403** |
| `cbioportal-datahub.s3.amazonaws.com/luad_tcga.tar.gz` | **HTTP 403** |
| `media.githubusercontent.com/media/cBioPortal/datahub/.../luad_tcga_pan_can_atlas_2018.tar.gz` | **HTTP 404** |
| `www.cbioportal.org/api/studies/luad_tcga_pan_can_atlas_2018` | 200 |

Only the REST API answers, and it is **not usable under ADR 0005**: it serves
per-study/per-gene JSON from a live service, so there is no single artifact to
pin a SHA-256 to, and a comparator that can drift silently is not a comparator.
Paginating 20,000 genes × 500 patients out of a REST endpoint would also put
substantial new code between this project and its external check.

**Xena serves two static files, both pinned:**

| artifact | bytes | sha256 |
|---|---|---|
| `TCGA.LUAD.sampleMap/HiSeqV2.gz` | 30,884,695 | `5432e5a6…47aa00` |
| `survival/LUAD_survival.txt` | 34,045 | `224589c8…78e115` |

Both carry `last-modified: 2021-04-08` and stable ETags, so both are `fixed`
under ADR 0005 and a digest mismatch is a hard failure. `resources/tcga/` becomes
the project's **sixth** sanctioned writer.

The survival file is the **Liu et al. 2018 TCGA-CDR** standardised endpoint set,
which is the curated resource the field uses rather than raw TCGA follow-up.

### 2. The family is 14: this study's 8 plus TCGA's 6

ADR 0026 §1 fixed the primary family as the **assessable** cells and left TCGA's
count to be measured here. It is **6** — all six signatures.

**The detection floor does not transfer, and must not.** Phase 2's floor is
`q3 > 2 ×` that AOI's `NegProbe-WTX` (ADR 0007) — a property of *this* assay.
TCGA LUAD is bulk RNA-seq with no negative-probe channel, so nothing in it
corresponds to that threshold. Restricting TCGA to the five signatures assessable
in the GeoMx lung arm would import a limitation of this assay into a cohort that
does not have it.

`myeloid_m1` is therefore tested in TCGA and **its result is reportable**, even
though it is not assessable at either GeoMx site. It carries the label
**"TCGA-only, no GeoMx comparator"**, and it enters the family because a result
that may be reported belongs in the denominator that corrects it.

**Family = 8 + 6 = 14.** P5-T5 applies Benjamini–Hochberg once over all 14
(ADR 0026 §2). **P5-T4 computes no q-value**, and the rule asserts no `q_` column
reaches any output — the same guard P5-T3 carries, for the same reason.

A TCGA-side **gene-presence** table is still computed, because a signature gene
absent from the Xena matrix cannot be scored. It is labelled as *presence*, never
as detection: conflating "the gene is not in this matrix" with "the gene is at
background in this assay" would merge two unrelated limitations.

### 3. The endpoint is OS, and the day→month conversion is declared

`OS` / `OS.time`, the direct analogue of the GeoMx lung column *"Primary lung
cancer diagnosis to death (Months)"* — all-cause death from diagnosis — and what
the source paper itself used for its TCGA LUAD panel.

DSS was considered and rejected: it counts only cancer-attributed deaths, so it
is cleaner biologically but is **not what the GeoMx column measures**, and an
external check that swaps the endpoint compares two different questions.

**`OS.time` is in DAYS and the GeoMx endpoint is in months.** The conversion is
`365.25 / 12 = 30.4375` days per month, declared in `config/tcga.yaml` rather
than inlined, so it is reviewable beside the data it applies to.

### 4. Median split is primary; a continuous Cox is a declared sensitivity

**Comparability with P5-T3 is the entire point of this task.** A different
estimator would confound "different cohort" with "different method", which is
precisely the inference P5-T4 exists to license. So the primary is the same
median split, the same `score > median` tie side, and the same log-rank.

At n ≈ 500 a median split discards real information, and that is where TCGA's
power actually shows — so a **continuous Cox on the score** is run as a declared
sensitivity and reported beside the primary (ADR 0012's shape). Neither is the
tiebreak for the other. A primary null beside a significant continuous fit would
be a statement about the split, not about the biology, and it must be readable as
such.

### 5. Two limitations that are part of the deliverable, not footnotes

**This is not the same measurement.** A GeoMx `TIME-L` score comes from the
**PanCK-negative segment** of an ROI sited in a CD45-rich region — and per Q2's
live caveat that is not a CD45-sorted population either. A TCGA score comes from
**whole bulk tumour**: tumour cells, stroma and immune cells together. "The same
signature" is computed on two different things. **A disagreement between the
cohorts is therefore not necessarily a disagreement about biology**, and an
agreement is not necessarily a replication. This is the same class of caveat as
Phase 2's reference-matrix confound (`config.yaml → contexture.deconvolution`).

**There is no brain comparator at all.** TCGA LUAD is primary lung. The brain arm
— the one carrying `TIME-B` n = 8, the binding constraint on the entire project —
gets **no external check from this task**, and P6-T6's limitations document must
say so. P5-T4 contextualises the lung null only.

## Consequences

- P5-T5's family is **14**, asserted against `survival_km_models.tsv` and
  `tcga_km_models.tsv` rather than quoted from here.
- `resources/tcga/` is the sixth sanctioned writer and P6-T2's RO-Crate must
  describe it.
- **No `config/config.yaml` edit.** `survival.tcga_yaml` was named at P5-T1
  precisely so this task would cost no second Phases 1–4 rebuild, and the pin
  lives in `config/tcga.yaml`.
- **No conda environment edit.** `gseapy`, `scikit-survival` and `statsmodels`
  are all already present, so ADR 0023's pin-regeneration transition is avoided.
- ADR 0024 §1's cBioPortal preference is **superseded for this artifact only**.
  If the datahub returns, re-pinning is a `config/tcga.yaml` edit and costs
  nothing — which is why the pin was kept out of `config.yaml` in the first place.
