P1-T5 landscape explorer — the Phase 1 ordination, interactive. Put any of the
30 principal components on either axis, switch to the UMAP, recolour by
compartment / site / patient / AOI code / DSP run, and search the AOI metadata
table.

**Presentation only.** This is an app-tier notebook (``notebooks/apps/``), so it
carries no analysis logic (``CLAUDE.md`` hard constraint 4): every value it
shows was produced by ``p1t1_pca_landscape`` and is read from
``pca_coords.tsv``, ``umap_coords.tsv`` and ``pca_summary.json``, all three
declared as inputs of the exporting rule. An undeclared read would be a silent
provenance hole.

**The thing to do with it** is the Phase 1 result in one gesture: colour by
**compartment**, then switch to **patient**. The first is structured, the second
is confetti. P1-T3 measured that split — compartment 41% of variance against
patient 22% by PVCA — and P1-T4 confirmed it from the clustering side
(compartment ARI 0.741 against patient 0.003, with same-patient AOIs landing in
the same cluster at *below* the chance rate). PROJECT_PLAN §6 expected patients
to cluster together across compartments; they do not, and this is the figure
where a reader can check that for themselves rather than take it on trust.

The **join each patient's AOIs** switch draws a line through every subject's
AOIs. It is the most direct way to see the same thing: the lines cross the whole
ordination instead of staying local.

**QC-flagged AOIs are ringed in red, not removed** (``qc.flag_only``). Whether
they sit apart from their compartment is exactly the question flagging exists to
raise, and it can only be answered by leaving them in.

**Running it.** The export carries its own Python runtime (Pyodide) and its own
copy of the three pipeline outputs, so it needs no kernel, no server-side Python
and no network. It **must be served over HTTP** rather than opened as a
``file://`` URL — partly a browser restriction on WebAssembly, and partly
because the notebook fetches its data from ``public/`` alongside the page::

    python -m http.server --directory results/reports/landscape_explorer

That data path is the reason the tables are interactive rather than static. The
exporting rule stages the notebook next to a ``public/`` folder, which
``marimo export html-wasm`` copies into the export, and ``mo.notebook_location()``
resolves to the served URL in the browser. The alternative — bundling the frames
through ``[tool.marimo.runtime] cache_cells`` — cannot work with ``mo.ui.table``:
that widget holds a locally defined ``PandasTableManager`` which cannot be
pickled, so its cell output caches as an ``UnhashableStub``, and an
unserializable output makes marimo re-run the cell's ancestors live — in the
browser, a loader with no filesystem to read. The tables then render as stub
text. See ADR 0010.

The axis dropdowns apply to the PCA; the UMAP has exactly two axes, and the
notebook says so rather than leaving the dropdowns looking live but ignored.

``TIME-B`` n = 8 throughout.
