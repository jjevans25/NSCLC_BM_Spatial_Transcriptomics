---
name: geomx-dsp-analysis
description: >-
  NanoString GeoMx DSP assay vocabulary and analysis conventions — ROI vs.
  segment vs. AOI, Q3 normalisation, LOQ-based detection filtering, why AOIs
  within a patient are correlated, and the SpatialDecon reference-matrix
  caveat. Use whenever working with GeoMx DSP data, GeomxTools, standR, or
  SpatialDecon, or when writing about AOIs, compartments, or detection limits.
---

# GeoMx DSP analysis

Assay-specific knowledge that the general scientific skill libraries do not
cover. Written for `nsclc-brainmet-time` (GSE200563), but the assay facts are
general.

## Vocabulary: ROI, segment, AOI

These are routinely used interchangeably in papers and are **not** the same
thing. Getting them wrong changes what an n means.

| Term | What it is |
|---|---|
| **ROI** | Region of interest. The area an operator circles on the slide. |
| **Segment** | A sub-area *within* one ROI, isolated by masking on a morphology marker (PanCK⁺ vs. PanCK⁻, CD45⁺, GFAP⁺). Optional. |
| **AOI** | Area of illumination. The area actually exposed to UV so its oligos are cleaved and collected. **One AOI = one column of counts = the unit of analysis.** |

One ROI yields one AOI if unsegmented, or several AOIs if segmented — one per
segment. So **AOIs are not independent samples**: several can come from one
ROI, several ROIs from one slide, several slides from one patient.

**The distinction that matters most here:** antibody-based segmentation gives
you a genuine molecular *compartment*. Geometric ROIs with descriptive names
give you a *region* — a pathologist's judgement about where to draw a circle.
Both are legitimate; only the first supports a claim of the form "expression in
the immune compartment". Never assume which one a dataset used — check the
methods. In this project that is open question **Q2**, and it is unresolved.

## Q3 normalisation

Per-AOI counts scale with the area illuminated and the number of cells in it,
neither of which is biology. Normalisation must remove that.

**Q3 (upper-quartile):** divide each AOI's counts by the 75th percentile of its
own counts (conventionally over targets above LOQ), then rescale to the
geometric mean of those quartiles across AOIs.

Why Q3 rather than the alternatives, for WTA:

- **Not total counts.** In WTA most targets sit near background, so a library
  size is dominated by noise and tracks background more than signal.
- **Not housekeeping genes.** Assumes a stable HK set across compartments.
  Tumour, immune and glial AOIs are different cell types with different
  metabolic rates — that assumption is doing real work and is likely false.
- **Not background/negative-probe normalisation** *unless* signal-to-background
  is poor. It is the recommended fallback for low-signal AOIs, and it removes
  real signal when applied to healthy ones.

Q3 assumes the upper quartile is comparably composed across AOIs. Check that
assumption rather than asserting it: plot Q3 against library size and against
negative-probe geomean per AOI, coloured by compartment. A compartment whose Q3
sits systematically apart is a warning, not a detail.

**Verify before applying.** Public "processed" matrices are often already
normalised, sometimes already logged. Applying Q3 twice is silent and
destructive. Integer values ⇒ raw; a near-constant column-wise Q3 ⇒ already
normalised; a 0–20 range ⇒ already logged.

## LOQ and detection filtering

Each panel carries **negative probes** — sequences with no target — which
measure background *in that AOI*. Background varies by AOI, so the limit of
quantitation is per-AOI, not a global constant.

Conventionally:

```
LOQ_aoi = geomean(NegProbe_aoi) * geoSD(NegProbe_aoi) ^ n     # n = 2 typically
```

Two filters follow, and they are different:

- **AOI-level:** gene detection rate = fraction of targets above that AOI's LOQ.
  A low rate means the AOI under-performed — too few cells, poor tissue, failed
  collection.
- **Gene-level:** keep targets detected above LOQ in at least *X*% of AOIs.

**The censoring point, which matters more than it looks.** A value below LOQ is
**not zero and not absent** — it is censored: present but unquantifiable. Two
consequences:

1. Treating sub-LOQ as zero biases every downstream mean and fold change.
2. **Immune-checkpoint transcripts frequently sit near LOQ.** A "difference" in
   a checkpoint gene between compartments can be entirely a difference in how
   often it cleared background. Always report the detection rate alongside the
   effect, and test sensitivity to the LOQ multiplier — a result that survives
   only at n=1 is a result about the threshold.

**Flag, don't silently drop.** Every excluded AOI gets a row and a reason. In a
design where one compartment has n=8, losing two AOIs invisibly is the
difference between an analysis and a fiction.

## Why AOIs within a patient are correlated

Shared genetics, treatment history, fixation, staining batch, slide, scan. AOIs
from one patient are more alike than AOIs from different patients, for reasons
unrelated to the compartment being compared.

Treating them as independent inflates n and shrinks p-values. It is the most
common and most consequential error in GeoMx reanalyses.

**A random intercept for patient is mandatory, not stylistic** — `lme4`,
`lmerTest`, with `emmeans` for contrasts. This holds even where each patient
contributes one AOI per compartment: the paired structure is the point of the
design and discarding it wastes the very information that makes the comparison
worth doing. If a model will not converge, simplify the fixed effects before
you touch the random structure.

## SpatialDecon: the reference-matrix caveat

`SpatialDecon` estimates cell-type abundance per AOI against a **reference
profile matrix**. The estimates are only as good as that reference's match to
the tissue.

- **The reference must match the site.** `safeTME` is built for tumour
  microenvironment and has no CNS-resident populations. Applied to a brain AOI
  it cannot represent astrocytes, oligodendrocytes or microglia, so their signal
  is forced onto whatever cell types the matrix *does* contain — usually
  macrophages. Comparing a lung AOI to a brain AOI through a lung-appropriate
  matrix manufactures a difference.
- **Microglia and macrophages are barely separable** by transcriptome. In a
  brain-metastasis analysis this is precisely the distinction of interest, so
  state it as a limitation rather than a result.
- **Deconvolving an already-segmented AOI answers a narrower question.** In a
  CD45⁺ AOI, "immune fraction" is near 1 by construction. What deconvolution can
  offer is composition *within* the immune compartment — and the reference must
  be able to represent that composition.
- Abundance estimates are **relative and compositional**. Increases and
  decreases are coupled; treat them as such, not as independent measurements.

## Language discipline

- Adjacent compartments profiled separately do not show **colocalisation**.
  Write "inferred crosstalk between adjacent compartments."
- A ligand–receptor pair scored across two AOIs is a **nomination**, not an
  interaction.
- Sub-LOQ is "below the limit of quantitation", never "not expressed".
- Any statement about brain immune contexture in this project states **n = 8**
  inline.
