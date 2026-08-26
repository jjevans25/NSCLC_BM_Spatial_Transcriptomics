P1-T3 variance partition — the number Gate 1 turns on (Aim A2). Each AOI factor
enters one mixed model per gene as a **crossed random intercept**
(``patient`` / ``site`` / ``compartment`` / ``dsp_run`` + residual, REML), so a
factor's share is shrunk toward zero in proportion to how little each of its
levels contributes.

**Left — per-gene partition.** The distribution of variance shares across all
converged genes, with **both null floors drawn on top**: a red cross for the
empirical permutation null (factor labels shuffled, model refit) and a black
tick for the analytic :math:`\eta^2` null, :math:`(k-1)/(n-1)`. The floors are
drawn rather than described because they are large: ``patient_id`` has 42 levels
across 120 AOIs, so its :math:`\eta^2` null alone is 0.345 and a naive read of
"patient explains half the variance" is close to an artefact of the level count.
A component is only interesting to the extent it clears its own floor.

**Centre — PVCA.** The same variance components fitted on the leading
eigenvectors of the sample–sample correlation matrix and averaged with
eigenvalue weights. This is the single dataset-level number the headline
sentence quotes; the left panel is its distributional backing.

**Right — what dropping** ``dsp_run`` **absorbs.** P1-T2 found PC1 shifting
consistently positive in DSP run B inside every compartment that spans both runs
(L +19.8, LB +9.7, TIME-L +4.7, TIME-B +7.6), so batch left out of the model does
not vanish — it lands somewhere. The paired bars show where, and the annotated
deltas are the amount. Two caveats travel with this term and cannot be fitted
away: ``dsp_run`` has only **two levels** and is **near-nested in patient** (just
1 of 42 patients spans both runs), so its own component is poorly estimated; and
it is estimable only in ``L``, ``LB``, ``TIME-L`` and ``TIME-B``, since ``mLN``,
``TBME`` and ``BC`` sit wholly within one run (Q3).

``site`` and ``compartment`` are structurally confounded — ``lymph_node`` is
tumour-only, ``glial_stroma`` and ``normal_control`` are brain-only — so a model
carrying both splits their shared variance by convention. Sensitivity fit ``S1``
in ``variance_partition_summary.json`` reports the alternative in which they are
a single 7-level ``aoi_code`` factor; ``S3`` reports the fit without the ``BC``
controls.

``TIME-B`` n = 8 throughout.
