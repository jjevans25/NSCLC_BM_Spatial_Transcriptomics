P3-T5 checkpoint panel **detection** by compartment (Aim A4, **exploratory** per
ADR 0008). Rows are the nine pre-registered genes of ``config/checkpoints.yaml``
(P3-T1, ADR 0015), ordered by how many compartments each clears — the sort key
is printed on the right-hand axis. Columns are all seven compartments, blocked
by site, covering all 120 AOIs. Dot area is the fraction of that compartment's
AOIs in which the gene exceeds **2.0 × that AOI's own** ``NegProbe-WTX`` level
(ADR 0007); colour is the median of that ratio across the compartment.

**What this figure is not.** It is a **detection** figure. It carries no effect
estimate, no confidence interval and no q-value; those are in
``checkpoint_carrier_models.tsv`` (P3-T3) and ``checkpoint_shift_models.tsv``
(P3-T4). The distinction is load-bearing, because **detection and expression
move in opposite directions under the same background gradient**: detection is
``q3 > 2 × NegProbe``, so higher background raises the bar and a gene is
detected *less* — it looks depleted; expression is ``log2(q3 + 1)``, where
background *adds* to the signal — it looks enriched. ADR 0018 records this
project making exactly that conflation once, and inverting an argument on it.
Reading a compartment difference off dot size is reading detection; reading it
off a model is reading expression. They are not the same claim.

**Why the colour is not "mean expression".** PROJECT_PLAN §6 specifies
"colour = mean expression". **ADR 0016 §5, written at P3-T2 before this figure
existed, changed it to** ``median_negprobe_ratio``. On a mean-log2 scale an
*undetected* gene still shows its own background level, and brain background is
higher — so undetected genes would render **brighter in brain than in lung**,
carrying the ADR 0008 artefact straight into the deliverable. The
background-normalised ratio is immune to the gradient by construction and is the
statistic ADR 0008's own reconnaissance table used. Both ``mean_log2_all`` and
``mean_log2_detected`` remain in ``checkpoint_detection.tsv``;
``mean_log2_detected`` was rejected separately because its denominator changes
per cell, so two dots of the same colour could rest on 15 AOIs and on 4.

**The colour midpoint is a real threshold, and the scale is piecewise-linear.**
The midpoint is not a display convenience: it is ``qc.detection_background_multiple``
read from the config, the same 2.0 that defines detection per AOI. Below it, the
median AOI in that compartment is indistinguishable from its own background.
Because the data are asymmetric about it (0.84 to 5.11), a ``TwoSlopeNorm`` maps
1.16 ratio units below the midpoint and 3.11 above onto the two halves of the
ramp independently — nothing is clipped, but **colour distance is not ratio
distance across the midpoint**. Read the colourbar's numeric ticks, not the
apparent gradient.

**In greyscale the colour is ambiguous by construction, and nothing depends on
it alone.** A diverging map sends both ends of the scale to dark grey, so a dot
at 0.84 and one at 5.11 look alike without colour. Every claim the figure makes
is carried redundantly by something achromatic — dot area, the hatch, the bold
edge, and the printed ``n_detected/n_aoi`` under every dot — so a greyscale
reader loses the ratio and nothing else.

**Colour side and hatching coincide on all 63 cells, and that is not a
coincidence.** ``ratio > 2.0`` and ``assessable`` are the same pre-registered
rule seen two ways: a gene is assessable when more than half its AOIs exceed 2×
background, and the median exceeds 2× under the same condition. They can
legitimately diverge only where ``n_aoi`` is even and the detection rate sits
*exactly* at the floor. No such cell exists here; the script logs the
concordance count every run as an observation, never as an assertion.

**Hatched cells mark a gene below its detection floor, and the hatch is
load-bearing.** The floor is ``detected_in_aoi_fraction: 0.5`` —
**pre-registered in ADR 0016 §1, committed before** ``p3t2_checkpoint_detection``
**existed**, and taken verbatim from ADR 0008 point 4 rather than invented for
Phase 3. For a hatched cell the finding is **"not assessable in that
compartment"**, never "lower in that compartment". **40 of 63 cells are
hatched; 23 clear the floor; the panel clears** ``min_assessable_fraction`` **in**
``TIME-L`` **alone.** Per Gate 3 that is a legitimate finding, not a failure —
and it is a finding only *because* the detection table exists to license it.
"Not assessable" is deliberately weaker than "absent": at these levels, with
hundreds of genes sharing an identical value in an AOI, absence and assay
insensitivity are not separable, and the mark is what stops a reader collapsing
them.

**Rings are a measured zero, not missing data.** Nine cells — ``TIGIT`` in
``L``, ``IDO1`` in ``TBME``, and seven of nine genes in ``BC`` — are detected in
zero AOIs. They are rendered distinctly rather than omitted, for the same reason
the hatch exists: the pre-registered panel reaches the figure whole, and which
genes turn out to be measurable is the Phase 3 *result*, not a Phase 3 input.

**The binding constraint is on the figure and applies to every dot:**
``TIME-B`` **n = 8.** The lung-vs-brain immune contrast detects roughly 1.1–1.3
SD at 80% power (P0-T8), so a small dot in ``TIME-B`` is an **assay-sensitivity
limit, not evidence of absence**. Per Gate 2, a surprising difference here is
more likely a pipeline bug than biology.

**mLN, TBME and BC (†) are shown for detection only.** Each sits **wholly within
one DSP run** — ``TBME`` 20/0, ``mLN`` 13/0, ``BC`` 0/7 (P0-T3, Q3) — so batch
is inseparable from biology there, and ADR 0009 §3 is explicit that no covariate
recovers it. They are audited and plotted but **never modelled** (ADR 0016 §2).
This narrows §6's "tumour vs. immune vs. glial" to tumour vs. immune, and it is
stated so a reader does not conclude the glial compartment was overlooked. **It
is not a statement about data quality.**

**QC, and what a compartment label is not.** Five QC-flagged AOIs are retained
under the flag-don't-drop policy (four in ``TBME``, one in ``LB``);
``detection_rate_qc_clean`` is carried in the table and **no cell changes
assessability under it**. And Q2's live caveat: a ``TIME`` AOI is the
PanCK-negative segment of an ROI sited in a CD45-rich region — **not** a
CD45-sorted population. The compartment label is not a cell-type label.
