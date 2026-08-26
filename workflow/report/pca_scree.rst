P1-T1 feature selection and PCA. Highly variable genes are selected on
log2(Q3 + 1) values with scanpy's ``seurat`` flavour, which expects log input;
the count is set by ``landscape.n_top_genes`` in ``config/config.yaml``.

PCA is **zero-centred but not unit-scaled**. Scaling every gene to unit variance
is the single-cell convention, but it would give a gene detected in three AOIs
the same leverage as one detected in all 120 — and P0-T5 established that
detection varies systematically by compartment (``TBME`` 21.4% against ~42% in
the tumour cores). Unit scaling would therefore promote a detection artefact
into a principal component. The bulk-style zero-centred PCA is the more
conservative reading of the same data.

Left panel is the scree; right is the cumulative curve with the number of
components needed for 80% of variance. PC1–PC4 loadings are in
``results/tables/pca_loadings.tsv`` and the per-AOI coordinates in
``results/tables/pca_coords.tsv``.
