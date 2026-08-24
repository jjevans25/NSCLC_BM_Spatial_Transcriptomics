P0-T6 marker sanity check. Boxplots of the PROJECT_PLAN §6 marker panel across
the seven compartment labels, on log2(Q3-normalised + 1) values for all 120
AOIs. Boxes are coloured by the ``compartment`` field of
``config/compartment_map.yaml``; the dashed red line is the median
``NegProbe-WTX`` level, i.e. the background floor, so a group sitting on that
line has no signal rather than merely less of it.

This tests whether the compartment labels mean what they say. Q2 established
that the AOIs are genuinely antibody-segmented, but **PanCK was the only
collection mask** — CD45 and GFAP guided where a pathologist placed each ROI
rather than sorting cells — so the transcriptome is the only independent
evidence that placement worked. The gating criteria are §6's: tumour AOIs
epithelial-high, ``TIME`` AOIs ``PTPRC``-high, ``TBME`` AOIs ``GFAP``-high.

Effect sizes and confidence intervals for each criterion are in
``results/tables/marker_sanity_verdict.tsv``; per-group medians and quartiles
are in ``results/tables/marker_sanity.tsv``. Bootstrap CIs use the seed in
``config/config.yaml``.

Run **before** P0-T5, deliberately: this checks the labels, not the data
quality, and a wrong mapping should surface before effort goes into QC.
``TIME-B`` n = 8 throughout.
