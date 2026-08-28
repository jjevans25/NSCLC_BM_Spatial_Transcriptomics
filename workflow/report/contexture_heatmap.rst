P2-T7 immune contexture, brain metastasis vs. primary lung, resolved **within**
the ``TIME`` compartment (Aims A1/A3). Rows are the six versioned signatures of
``config/signatures.yaml`` (P2-T1, ADR 0011); columns are the 23 ``TIME`` AOIs,
blocked by site and ordered by patient within block. Colour is the **z-score
mean** — each gene standardised across these 23 AOIs, then averaged over the
set — on a symmetric, colourblind-safe diverging scale whose midpoint is a true
zero.

**What this figure is not.** It is a display of the scores, not the test. The
contrast of record is the mixed model in ``signature_models.tsv``:
``score ~ site + (1|patient_id)``, with the random intercept mandatory on
grounds of non-independence rather than variance share (ADR 0009), refit with
``dsp_run`` as a reported sensitivity (ADR 0012). Reading an effect off the
colours instead of off that table skips the patient structure entirely — 33 of
42 subjects contribute more than one AOI.

**Why ssGSEA is not the one plotted.** Both scoring methods were computed and
both are modelled; their agreement is the robustness evidence at this *n*
(``signature_scoring_summary.json`` carries the per-set Spearman ρ). Only the
z-score mean is comparable *between* signatures, because ssGSEA's NES is a rank
statistic whose magnitude depends on set size, and these sets run from 4 to 13
genes.

**Hatched cells mark a set that failed its detection-coverage floor at that
site, and the hatching is load-bearing.** Detection here is
background-relative — the matrix has no zeros at all (Q1), so a gene counts as
detected only above a multiple of that AOI's NegProbe-WTX level — and **brain
background is higher** (ADR 0008). A near-background gene therefore reads as
*depleted in brain* artefactually, which is precisely the misreading a heatmap
invites. For a hatched set the finding is **"not assessable in brain"**, never
"lower in brain". Per-set, per-site coverage is in ``signature_coverage.tsv``;
it is reported per site and never pooled, because P0-T5 established that
detection varies systematically by compartment.

**The binding constraint is stated on the figure and applies to every cell of
it: ``TIME-B`` n = 8.** The lung-vs-brain immune contrast detects roughly
1.1–1.3 SD at 80% power (P0-T8), so a signature showing no difference here is
**uninformative, not negative**. Per Gate 2, a *surprising* difference is more
likely a pipeline bug than biology.

A ``TIME`` AOI is the PanCK-negative segment of an ROI sited in a CD45-rich
region — **not** a CD45-sorted population (Q2's caveat). The compartment label
is not a cell-type label.
