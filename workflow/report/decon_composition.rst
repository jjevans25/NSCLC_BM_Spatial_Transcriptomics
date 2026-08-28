P2-T5 cell composition of the 23 ``TIME`` AOIs (Aim A3), estimated with
``SpatialDecon`` and renormalised to proportions **within each AOI**. Columns are
AOIs labelled by patient; colour is cell type.

**Panel A is the only panel whose brain-vs-lung comparison means anything.**
``safeTME`` (906 genes × 18 types, collapsed to 14) is applied to *both* sites,
so one reference underlies both halves and the difference between them is a
property of the data. **Panels B and C are PROJECT_PLAN §6's two-matrix
approach** — ``Lung_HCA`` for ``TIME-L``, ``Brain_Darmanis`` for ``TIME-B`` — and
they are drawn separately, with separate legends, on purpose.

**The reference-confounding caveat is specific and it is measured, not asserted.**
``Lung_HCA`` resolves **9** lymphoid cell types (CD8, CD4, Treg, B, plasma, NK,
pDC and others); ``Brain_Darmanis`` — a glioblastoma infiltrating-front
reference, not normal brain — resolves **0**, carrying only Neoplastic, Myeloid,
Astrocyte, OPC, Oligodendrocyte, Neuron and Vascular types. The two-matrix arm
therefore **cannot structurally** show a lymphoid difference between sites: a
lymphoid gap between panels B and C would be the reference, not the biology.
Both counts come from ``decon_summary.json``, computed by the rule, so the
caption cannot drift from what was fitted.

**What the deconvolution was not given.** ``spatialdecon`` was called **without**
``raw``, so its Poisson error model and per-observation weights are not in use.
That is not an oversight: raw gene-level counts do not exist for this dataset —
Q1 established the DCCs are probe-level and GEO carries no PKC, so gene-level
counts are not recoverable even in principle. The fit is unweighted.

Background is each AOI's NegProbe-WTX level broadcast across genes, the same
background quantity the project's detection rule uses
(``qc.detection_background_multiple`` × negprobe), so deconvolution and detection
are anchored to one definition rather than two.

**``TIME-B`` n = 8**, and the minimum detectable effect for the site contrast is
~1.1–1.3 SD at 80% power (P0-T8). Deconvolution proportions here are
**descriptive**: no p-value is attached to a composition difference, and the
statistical contrast of record remains the mixed model in
``signature_models.tsv``. P2-T6 compares the two.

A ``TIME`` AOI is the PanCK-negative segment of an ROI sited in a CD45-rich
region — **not** a CD45-sorted population (Q2's caveat). A deconvolved
proportion is an inference from bulk expression within that segment, not a cell
count.
