# Analysis notes

Prose that belongs with the results but not inside a script. Every number here
comes from a table a Snakemake rule produced; none of it originates in this
file or in a notebook (CLAUDE.md hard constraint 3).

---

## P2-T6 — do deconvolution and signature scoring agree? (2026-08-28)

**Source tables:** `results/tables/convergence_check.tsv`,
`convergence_summary.json`, `decon_composition.tsv`, `signature_models.tsv`,
`signature_coverage.tsv`.

**They agree, on the one signature where the question can be asked at all —
and that restriction is the finding, not a caveat attached to it.**

Of six signatures, **one** survives as a genuine convergence test. For
`cytotoxicity`, the mixed model puts the brain-vs-lung contrast at
**−0.946 SD [−1.590, −0.302]** (q = 0.018, `TIME-B` n = 8) and the safeTME
deconvolution puts the summed CD8 T-cell + NK proportion at **0.139 in `TIME-L`
against 0.081 in `TIME-B`** (Δ = −0.058). Two methods that share no machinery —
one a rank/z summary of four genes under a mixed model, the other a
constrained regression against a 906-gene reference matrix — point the same way.
That is the strongest single result in Phase 2.

The other five are disqualified, each for a reason worth stating separately:

| signature | why it is not a convergence test |
|---|---|
| `antigen_presentation` | **No counterpart.** It is a transcriptional program expressed by tumour cells and professional APCs alike, not a lineage `safeTME` resolves. Mapping it to "B + mDCs + macrophages" would silently swap a functional program for an abundance. |
| `exhaustion` | **Not assessable in brain** (1 of 6 genes above background in `TIME-B`). |
| `tls` | **Not assessable in brain** (1 of 5). |
| `myeloid_m1` | **Not assessable at either site** (2 of 10). |
| `myeloid_m2` | **Degenerate mapping.** It names the same cell type as `myeloid_m1`, because M1/M2 is a polarisation axis *within* a lineage and deconvolution resolves the lineage only. Macrophage abundance rising cannot confirm a polarisation score. |

**Which is trusted where they disagree: the signatures.** Deconvolution's answer
is conditional on the profile matrix chosen; a signature score is not. That is
not a preference, it is visible in this dataset — see the reference-confounding
note below.

**The one apparent disagreement is uninformative, and the coverage gate is what
tells us so.** `myeloid_m1` has the signature falling in brain (−0.292 SD) while
macrophage abundance *rises* (0.200 → 0.312). Read naively that is a
contradiction between methods. It is not: `myeloid_m1` clears its detection
floor at **neither** site — 8 of its 10 genes (CD80, CD86, CXCL10, CXCL11,
IL12B, IL1B, NOS2, TNF) sit below background — so its score is not a measurement
of M1 polarisation and cannot disagree with anything. Without the coverage gate
this would have been written up as a real methodological conflict.

`myeloid_m2` and macrophage abundance both rise in brain, which *looks* like
convergence and is not evidence of any: an M2 score and a macrophage proportion
are different quantities, and the same lineage underlies both arms of the
comparison.

### Reference confounding, measured rather than asserted

P2-T5 ran two arms. The primary applies `safeTME` to **both** sites, so one
reference underlies both halves and its brain-vs-lung comparison is a property
of the data. PROJECT_PLAN §6's two-matrix arm (`Lung_HCA` for `TIME-L`,
`Brain_Darmanis` for `TIME-B`) is reported as a sensitivity and **its cross-site
comparison must not be made** — for a sharper reason than "different
references":

> `Lung_HCA` resolves **9** lymphoid cell types. `Brain_Darmanis` — a
> glioblastoma infiltrating-front reference, not normal brain — resolves **0**,
> carrying only Neoplastic, Myeloid, Astrocyte, OPC, Oligodendrocyte, Neuron and
> Vascular types.

The two-matrix arm therefore **cannot structurally** show a lymphoid difference
between sites. A lymphoid gap between those panels would be the reference, not
the biology. Both counts are computed by `p2t5_deconvolution` and stored in
`decon_summary.json`, so the claim rests on a number rather than on prose. This
is why P2-T6 uses the safeTME arm alone.

### What the deconvolution was not given

`spatialdecon` was called **without** `raw`, so its Poisson error model and
per-observation weights are unused. Raw gene-level counts do not exist for this
dataset — Q1 established the DCCs are probe-level and GEO carries no PKC, so
they are not recoverable even in principle. The fit is unweighted, and every
composition number above should be read with that in mind. Deconvolution
proportions are **descriptive**: no p-value is attached to any of them, and the
contrast of record remains the mixed model.

A `TIME` AOI is the PanCK-negative segment of an ROI sited in a CD45-rich
region — **not** a CD45-sorted population (Q2's caveat). A deconvolved
proportion is an inference from bulk expression within that segment, never a
cell count.

---

## P3-T7 — which compartment carries PD-L1, and what would that mean? (2026-09-03)

**Source tables:** `results/tables/checkpoint_carrier_models.tsv`,
`checkpoint_shift_models.tsv`, `checkpoint_detection.tsv`,
`checkpoint_paired_concordance.tsv`, `decon_composition.tsv`.

**In the lung primary, PD-L1 transcript sits in the PanCK-negative compartment
rather than the tumour compartment — and that is the only checkpoint claim
Phase 3 supports.**

`CD274` is enriched in the lung immune compartment relative to the paired
tumour compartment by **+0.612 SD [+0.361, +0.863], q = 0.0001** (45 AOIs, 30
patients, `expression ~ compartment + (1|patient_id)`; detected 19/30 in `L`,
13/15 in `TIME-L`). It survives both pre-registered sensitivities — batch
+0.599, background +0.635 — and the background adjustment *raises* it, which
matters here: the immune compartments carry higher background, so that gradient
is permissive and could have manufactured an immune-side enrichment (ADR 0018).
It did not manufacture this one. `CD276` (B7-H3), the only other gene clearing
the primary restriction, shows no compartment preference at either site (+0.147
lung, +0.159 brain, both null), so this is not the panel-wide immune-side drift
an artefact would produce.

### What it would mean, and how far it goes

The routine NSCLC biomarker is PD-L1 protein on tumour cells, the basis of
KEYNOTE-024's TPS ≥ 50% selection (PMID:27718847). But tumour-cell and
immune-cell PD-L1 are separately scored and separately actionable: IMpower110
enrolled on PD-L1 expression in ≥ 1% of tumour cells **or** ≥ 1% of
tumour-infiltrating immune cells by the SP142 assay (Herbst et al. 2020,
*N Engl J Med* 383:1328–1339, PMID:32997907). A tumour whose PD-L1 is
predominantly stromal is one a tumour-cell-only assay can score as negative.

This dataset points at the second source in the lung primary. It does **not**
identify the cell type, and three limits say how far it goes, in order of what
they cost:

- **The compartment label is not a cell-type label.** A `TIME` AOI is the
  PanCK-negative segment of an ROI a pathologist sited in a CD45-rich region —
  **not** a CD45-sorted population (Q2). safeTME deconvolution puts `TIME-L` at
  a mean 0.200 macrophages, but also 0.198 fibroblasts and 0.060 endothelium:
  roughly a quarter of that compartment is not immune at all. **"Myeloid-derived"
  is a hypothesis this design cannot separate from stromal; what the contrast
  shows is "not tumour-cell-derived".**
- **Transcript is not protein, and neither is a TPS.** Nothing here is an IHC
  measurement, and the clinical biomarker is.
- **A4 is exploratory** (ADR 0008). This is reconnaissance that would justify a
  targeted study, not a claim that stands on its own.

### What Phase 3 cannot say, and why that silence is not evidence

**Brain is uninformative, not negative.** The same carrier contrast in brain is
+0.222 [−0.384, +0.829] with **`TIME-B` n = 8**, and the lung-vs-brain shift is
null in every assessable gene — largest `CD274` at −0.260 [−0.595, +0.075],
q = 0.462, again **`TIME-B` n = 8**. The design detects roughly **1.1–1.3 SD at
80% power** (P0-T8), so every one of these sits far below the floor. Writing them
up as "no difference between sites" would be a worse error than having no result.
The 5-patient paired check agrees with the unpaired direction for `CD274` (3 of 5
patients lung-higher) and deliberately carries no test statistic.

**Seven of nine genes never entered a primary fit.** A gene enters only if it
clears the detection floor in *both* groups compared (ADR 0016 §3, upheld as
Option C in ADR 0018); `CTLA4`, `PDCD1`, `LAG3`, `TIGIT`, `IDO1`, `HAVCR2` and
`VSIR` sit in a separate, unadjusted, exploratory table on which **no claim
rests**. Several show large apparent immune enrichments there, and every one is a
gene detected in 0–9 of 30 tumour AOIs. **The exclusion is itself the finding:**
40 of 63 gene × compartment cells are not assessable, and the panel clears the
floor in `TIME-L` alone.

**The glial compartment was not overlooked.** `TBME` — with `mLN` and `BC` —
sits wholly within one DSP run, so batch is inseparable from biology there
(ADR 0009 §3). Those three are audited and plotted but never modelled
(ADR 0016 §2), which is why §6's "tumour vs. immune vs. glial" is answered here
as tumour vs. immune only.
