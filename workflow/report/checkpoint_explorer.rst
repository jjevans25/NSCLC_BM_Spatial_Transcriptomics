P3-T6 checkpoint explorer — the Phase 3 detection audit, interactive. Move the
detection floor and watch which gene × compartment cells fall out of
"assessable"; narrow to a site or a subset of compartments; read the P3-T3 and
P3-T4 model estimates in their own section below the plot.

**Presentation only.** This is an app-tier notebook (``notebooks/apps/``), so it
carries no analysis logic (``CLAUDE.md`` hard constraint 4): every value it
shows was produced by ``p3t2_checkpoint_detection``, ``p3t3_carrier_models`` and
``p3t4_shift_models``, and all six files are declared inputs of the exporting
rule. An undeclared read would be a silent provenance hole. The one thing the
notebook derives is which cells to hatch at the floor you choose, and that is a
display flag, not a result.

**The slider is a sensitivity display, not a threshold control** (ADR 0019).
0.5 is pre-registered (ADR 0016 §1) and is the floor of record: every table,
model and q-value in Phase 3 was computed there, and at 0.5 this app reproduces
``checkpoint_dotplot.png`` exactly — same 40 hatched cells, same 9 rings, same
row order, same counts. The instant the slider leaves 0.5 the page says so in a
banner, so a screenshot taken at 0.75 carries its own caveat. Nothing the slider
reaches is a result, no model is refitted and no q-value changes. Setting a
floor *after* seeing which genes clear it is the specific failure P3-T2 exists
to prevent; the floor here was fixed in a commit that predates the data, and the
slider cannot unmake that choice.

The background multiple (2.0×, ADR 0007) is deliberately **not** exposed.
``detection_rate`` is already computed at it, so a control over it could only
relabel an axis — and recomputing detection inside a notebook is forbidden both
by hard constraint 3 and by ``p3t2_checkpoint_detection`` owning that rule.

**The thing to do with it** is the Phase 3 deliverable, which a static figure can
only assert: raise the floor and watch what survives. At the pre-registered 0.5,
**40 of 63 cells are not assessable and the panel clears the floor in**
``TIME-L`` **alone**. At 0.6 that is still true; by 0.75 no compartment carries
half the panel, and ``TIME-L`` is the last one standing throughout. ``TBME`` is
hatched at every value, and ``BC`` — the non-tumour brain control — is very
nearly empty, which is the negative control behaving. That stability is the
argument for the audit's conclusion, and this is where a reader can check it
rather than take it on trust.

**Dot area is detection, colour is background-normalised, and neither is an
effect size.** Colour is ``median q3 / NegProbe-WTX``, not mean expression: on a
mean-log2 scale an *undetected* gene still shows its background level, and brain
background is higher, so undetected genes would render brighter in brain than in
lung — the ADR 0008 artefact walking into the figure (ADR 0016 §5). The
colourbar's midpoint is read from ``qc.detection_background_multiple`` so it
cannot drift from the rule it depicts. ``TwoSlopeNorm`` makes the scale
piecewise-linear, so colour distance is not ratio distance across the midpoint;
every claim is carried redundantly by dot area, the hatch, the bold edge and the
printed counts, so a greyscale reader loses the ratio and nothing else.

**Model estimates sit in a separate section, and that separation is deliberate**
(ADR 0018). The plot is **detection** (``q3 > 2 × NegProbe``); the tables are
**expression** (``log2(q3 + 1)``). The two move *opposite* ways under the same
background gradient — higher background raises the detection bar so a gene looks
depleted, while it adds to signal so the same gene looks enriched — and this
project made that exact conflation once and reasoned from it. No effect size
appears on the dot plot and no dot appears in the tables. The four tabs are the
primary fits (the only Phase 3 output with a q-value), the two pre-registered
sensitivities, the genes the restriction *excluded* — the exclusion is itself
the finding — and the exploratory table, which carries no FDR and on which
**no claim may rest**.

Genes marked † — ``mLN``, ``TBME``, ``BC`` — are audited and plotted but never
modelled: each sits wholly within one DSP run, so batch is inseparable from
biology there (ADR 0009 §3, ADR 0016 §2). Not lesser data, a property of the
study design.

**Running it.** The export carries its own Python runtime (Pyodide) and its own
copy of the six pipeline outputs, so it needs no kernel, no server-side Python
and no network. It **must be served over HTTP** rather than opened as a
``file://`` URL — partly a browser restriction on WebAssembly, and partly
because the notebook fetches its data from ``public/`` alongside the page::

    python -m http.server --directory results/reports/checkpoint_explorer

That data path is ADR 0010, and it is the reason the tables are interactive
rather than static: bundling the frames through ``[tool.marimo.runtime]
cache_cells`` cannot work with ``mo.ui.table``, whose locally defined
``PandasTableManager`` cannot be pickled, so its cell output caches as an
``UnhashableStub`` and marimo re-runs the loader live — in the browser, against
no filesystem.

``TIME-B`` **n = 8**, and the lung-vs-brain immune contrast detects roughly
1.1–1.3 SD at 80% power (P0-T8). A small dot is an assay-sensitivity limit,
never evidence that the gene is absent. Aim A4 is exploratory (ADR 0008): no
result here may be a headline claim.
