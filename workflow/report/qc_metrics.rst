P0-T5 per-AOI quality control, by compartment. Metrics come from two places
because neither is sufficient alone.

The top row is **sequencing** quality, read from the DCC headers inside
``GSE200563_RAW.tar``: raw read depth, and saturation computed as
1 − deduplicated/aligned. These need no PKC. The dashed line on the saturation
panel is 0.50, the floor the source paper claims none of its ROIs fell below.

The bottom-left panel is **gene detection**, defined as the fraction of the
18,694 genes exceeding a multiple of that AOI's ``NegProbe-WTX`` level. It has
to be background-relative: the matrix contains no zeros at all (minimum value
2.12), so "non-zero" would mark every gene detected everywhere. The bottom-right
panel shows that background level itself, which is **not** constant across
compartments — it rises from lung tumour cores to brain parenchyma.

Library size is deliberately absent. The matrix is already Q3-normalised (Q1),
so its column sums are near-constant by construction and measure nothing; raw
reads in the top-left panel are the honest version of that metric.

Detection is markedly lower in ``TBME`` and ``BC`` than in tumour cores. That is
a real property of the tissue combined with its higher background, not a defect,
and it is why a single global cutoff would preferentially flag the glial
compartment. Policy is **flag, never drop** (``config qc.flag_only``): flagged
AOIs appear in ``results/tables/qc_excluded.tsv`` with a reason and remain in the
dataset. ``TIME-B`` n = 8 throughout, and losing 1–2 of those changes Phase 2.
