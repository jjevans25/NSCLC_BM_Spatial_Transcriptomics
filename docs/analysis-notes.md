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
