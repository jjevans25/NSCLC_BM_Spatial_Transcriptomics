P1-T4 hierarchical clustering of AOIs (Aim A2). AOI–AOI Spearman correlation
over the same 2000 highly variable genes the P1-T1/T2 ordination used, clustered
with average linkage on :math:`1-\rho`. Over all 18,694 genes every AOI pair
correlates at ~0.99 — the shared transcriptome dominates and no structure is
visible — so the HVG restriction is what makes this figure legible rather than a
uniform red square.

**Left — the heatmap**, ordered by the dendrogram above it, with three
annotation strips beneath: compartment, site, and patient. Patient is drawn
without a legend on purpose: with 42 levels the question is never *which*
patient a column is, only whether neighbouring columns share a colour.

**Right — the same question, answered numerically.** PROJECT_PLAN §6 expects
patients to cluster together across compartments ("*they usually do*"). Reading
a heatmap is not evidence for that either way, so the top panel sweeps the
dendrogram cut over k and reports the **adjusted Rand index** between the
resulting clusters and each factor, against a permutation null band (labels
shuffled, agreement recomputed). The bottom panel gives the simpler contrast:
mean correlation between AOIs **sharing** a factor level versus AOIs that do
not, with the difference annotated.

**The expectation is inverted in this dataset, and that is the result.**
Compartment structures the dendrogram; a patient's AOIs scatter across it. The
verdict sentence and every number behind it are in
``sample_clustering_summary.json``, with the full k sweep in
``cluster_agreement.tsv``. This is the clustering-side confirmation of the P1-T3
variance partition, which put compartment at 41% and patient at 22% by PVCA —
and it is why the ``(1|patient)`` random intercept in every downstream model has
to be argued from **non-independence** rather than from variance share.

Two things the figure cannot separate. ``site`` and ``compartment`` are
structurally confounded — ``lymph_node`` is tumour-only, ``glial_stroma`` and
``normal_control`` are brain-only — so their ARI curves are not independent of
each other. And the 7 ``BC`` controls are kept rather than removed, so that a
reader can see where they land instead of taking their exclusion on trust.

``TIME-B`` n = 8 throughout.
