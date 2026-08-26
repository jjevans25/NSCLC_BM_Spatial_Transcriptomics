# ADR 0010 — App-tier notebooks ship their data in `public/`, not in a cell cache

**Date:** 2026-08-26
**Status:** Accepted
**Task:** P1-T5
**Applies to:** every `notebooks/apps/` notebook exported by a rule, not just this one

## Context

P1-T5 requires an app-tier notebook exported to WASM HTML that "appears in the
Snakemake report under Interactive; opens offline with no kernel" (PROJECT_PLAN
§6), reading `pca_coords.tsv`, `umap_coords.tsv` and `pca_summary.json` as
**declared rule inputs** (§A.5).

A WASM export runs under Pyodide in the browser. There is no filesystem, so
`pd.read_csv("results/tables/pca_coords.tsv")` cannot work at view time. The
data has to reach the browser some other way, and the obvious mechanism is the
one marimo advertises for exactly this: `[tool.marimo.runtime] cache_cells =
true` with `marimo export html-wasm --execute`, which runs the notebook once at
export time and bundles cell outputs into the export.

**That mechanism was implemented, exported and tested in a browser. It does not
work for a notebook with tables.** Three findings, each of which cost a rebuild
and none of which is visible from `marimo check`:

1. **`mo.ui.table` cannot be cached, and failing to cache it is not benign.**
   The widget holds a `PandasTableManager` — a class defined *inside* a factory
   function, so unpicklable. Its cell output is therefore stored as an
   `UnhashableStub`. Per marimo's cached-lifecycle contract
   (`_runtime/executor/lifecycles/cached.py`), an unserializable output does not
   simply skip the cache: it forces a **live re-run of that cell's ancestors**.
   In the browser the ancestor is the loader, which has no files to read. Both
   tables rendered as `UnhashableStub` text.

   The plot and the controls survived, which is what made this confusing. A cell
   whose *defs* are live UI elements hits an explicit carve-out; a cell that
   *returns* UI it does not define does not.

2. **`mo.persistent_cache` is not a targeted substitute.** The natural fix —
   cache only the loader, leave the UI cells alone — fails silently. Only
   `cache_cells` writes the export manifest that the bundler reads
   (`dump_cache_manifests` is gated on `cache_cells_enabled` in
   `_runtime/callbacks/cache.py`). `persistent_cache` wrote to a local
   `notebooks/apps/__marimo__/` the exporter never looks at, and the export
   shipped **with no data at all** — no warning, no error.

3. **Pandas cannot read an `http://` URL under Pyodide.** There are no sockets,
   so urllib is non-functional and `read_csv(url)` raises. The failure is
   delayed and looks like success: the page renders the `--execute` preview,
   then throws an internal error seconds later when Pyodide boots and re-runs
   the loader.

## Decision

**App-tier notebooks ship their data as files in `public/`, and fetch it over
HTTP. No cell caching.**

- The exporting rule stages the notebook into a build directory under
  `results/interim/` beside a `public/` folder holding every input it reads.
  `marimo export html-wasm` copies `public/` into the export
  (`Exporter.export_public_folder`). Staging under `results/` rather than next
  to the notebook keeps generated files out of the source tree.
- The notebook resolves paths through `mo.notebook_location()`, which is the
  notebook's directory locally and the **served URL** in the browser, with a
  `mo.cli_args()` override and a `results/tables/` fallback so a bare
  `marimo edit` still opens it (§A.5).
- The loader branches on `sys.platform == "emscripten"` and reads through
  `pyodide.http.open_url()` there, plain paths everywhere else.
- `cache_cells` is **not** set, and `pyarrow` is not a dependency — it was only
  ever needed to deserialise the cached Arrow blobs.

## Consequences

- Tables stay fully interactive (search, sort, column summaries, sticky
  scrolling). This was the whole point: the cached alternative forced them to
  static HTML.
- The export **must be served over HTTP**. `file://` was already ruled out by
  the browser's WebAssembly restrictions; the same-origin data fetch is now a
  second, independent reason. `python -m http.server --directory
  results/reports/landscape_explorer`.
- Every file the app reads is duplicated into the export. At Phase 1 sizes
  (49 KB total) this is free. A future app over the full expression matrix must
  not copy the `.h5ad` in wholesale — it should read a rule-produced summary.
- The rule's `input:` remains the single declaration of what the app reads, so
  §A.5's provenance rule is preserved: the staged `public/` copy is derived from
  declared inputs, never a second source of truth.
- **`marimo check` does not protect against API skew.** The PEP 723 block
  installs *unpinned* `marimo` in both the export sandbox and Pyodide, which is
  not necessarily the version in `py-analysis`. `mo.ui.table(sortable=...,
  filterable=...)` — as documented in the `marimo-notebook` skill's `UI.md` —
  passed `marimo check --strict` and then failed at export runtime with
  `unexpected keyword argument`. Verify app-tier notebooks by running the export,
  not by linting them.
- The WASM lint rules (`marimo check --select MW`) passed throughout, including
  while the export was shipping no data at all. They check imports and packages,
  not whether the notebook's data actually arrives. A browser check is the only
  real acceptance test, and it is what caught all three findings above.
