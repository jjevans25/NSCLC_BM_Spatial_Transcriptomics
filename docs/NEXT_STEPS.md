# Next steps — Phase 1 (variance landscape)

**Last session:** 2026-08-26
**Branch:** `P1`, rebased onto `main` (`46466a0`, PR #2). Not pushed.
**Gate 0:** **PASSED.** Phase 0 merged. A4 demoted to exploratory (ADR 0008).
**Phase 1 status:** **COMPLETE.** P1-T1 … P1-T6 all done and committed.
Gate 1 met; what remains is the squash-merge that records it.

Phase 1's goal (§6): *the figure that justifies every mixed model that follows.*
Aim **A2** — what actually drives variance: patient, site, or compartment.

---

## Gate 1 — the answer

> **41% of variance is compartment-attributable.**

PVCA over the 14 leading PCs (60.4% of variance); per-gene median 21% over the
18,691 genes that converged of 18,694; permutation null 0.3%. `TIME-B` n = 8.
Table `results/tables/variance_partition.tsv`, headline and every sensitivity
fit in `results/tables/variance_partition_summary.json`.

**Compartment dominates. Patient is real but modest.** Full partition:

| Component | Per-gene median (all genes) | HVG-only median | PVCA weighted | Permutation null |
|---|---|---|---|---|
| compartment | 0.211 | 0.205 | **0.412** | 0.003 |
| patient | 0.124 | 0.141 | 0.218 | 0.026 |
| `dsp_run` | 0.078 | 0.034 | 0.096 | 0.002 |
| site | 0.006 | 0.005 | 0.023 | 0.003 |
| residual | 0.495 | 0.436 | 0.252 | 0.889 |

Four things in there matter more than the headline:

**a. The naive η² really was an artefact, and now there is a number for it.**
Median naive η² for patient is **0.447** against its own analytic floor of
**0.345** — barely clear of it. The variance-components fit puts patient at
0.124, against a permutation null of 0.026. So patient is genuinely present
(~5x its null) but nothing like "half the variance". Reporting the crude number
would have overstated it by roughly 3.5x.

**b. `dsp_run` absorption is confirmed but small.** Dropping the term (S2) moves
**+0.014** into compartment and **+0.023** into patient (HVG medians). So the
worry in the last session's note was directionally right — batch does inflate
compartment — but it is worth about one and a half percentage points, not the
distortion the PC1 shifts implied. Note also that most of what `dsp_run` holds
falls to *patient*, which is what you would expect from a term that is
near-nested in patient (1 of 42 patients spans both runs).

**c. Site does essentially no work.** 0.006 per-gene, 0.023 in PVCA — at or
barely above its null everywhere. Collapsing site and compartment into the
single 7-level `aoi_code` (S1) gives 0.190, against 0.205 for compartment alone.
Practically all the site signal is inseparable from compartment, which is what
the structural confounding (lymph_node tumour-only, glial_stroma and
normal_control brain-only) predicts.

**d. `BC` is not driving it.** Without the 7 controls (S3) compartment falls
0.205 → 0.184 and the ordering is unchanged.

**Consequence for P1-T6:** the ADR's argument **does** invert, exactly as
anticipated. `(1|patient)` cannot be justified on variance share — 12% is not a
dominant component. It is justified on **non-independence**: P12 and P24 each
have two `TIME-L`, P15 has two `TBME`, and five patients contribute paired
immune AOIs. Write that argument, not the one the plan expected.

---

## P1-T4 — hierarchical clustering: the expectation was inverted

Figure `results/figures/sample_correlation_heatmap.png`, k sweep in
`results/tables/cluster_agreement.tsv`, verdict in
`results/tables/sample_clustering_summary.json`. Spearman over the same 2000
HVGs, average linkage on 1 − ρ (median AOI–AOI ρ 0.400, range −0.045–0.872).

§6 said to "note whether patients cluster together across compartments (*they
usually do*)". **They do not**, and now there is a number rather than an
impression:

| Factor | Peak ARI (k) | Null p95 | Within vs between ρ |
|---|---|---|---|
| compartment | **0.741** (k=7) | 0.049 | 0.500 vs 0.320 (**+0.180**) |
| site | 0.034 (k=9) | 0.031 | 0.405 vs 0.384 (+0.021) |
| patient | **0.003** (k=3) | 0.004 | 0.517 vs 0.390 (+0.127) |

At k = 7 — where compartment agrees best — only **34%** of the 157 AOI pairs
sharing a patient land in the same cluster, against a **37%** baseline over all
pairs. Patient co-membership is at or slightly *below* chance.

**But do not stop at "patient does nothing", because that is the wrong read.**
Patient ARI is ~0 while within-patient correlation is +0.127, and those are not
in tension: compartment owns the top-level partition, and patient is a real
second-order effect that never becomes cluster structure. The controlled
contrast makes it explicit — restricted to the **102 same-patient AOI pairs
whose compartments differ** (so no part of it can be compartment agreement
wearing a patient label), mean ρ is **0.425 vs 0.318** for different-patient
pairs, **+0.107**.

That number is the one P1-T6 should quote. It is direct pairwise evidence of
non-independence that owes nothing to variance share, and it is measured across
compartments — exactly the structure `(1|patient)` exists to absorb.

---

## P1-T5 — landscape explorer: done

`notebooks/apps/landscape_explorer.py`, exported to self-contained WASM HTML by
`p1t5_landscape_explorer` into `results/reports/landscape_explorer/`. Appears in
the Snakemake report under **Interactive**. Dropdowns for embedding (PCA/UMAP),
x/y PC, colour-by, a switch that joins each patient's AOIs, and searchable,
sortable AOI metadata and design tables. **Verified working in a browser.**

**How the data reaches the browser, and why it is not cached.** The rule stages
the notebook into `results/interim/explorer_build/` beside a `public/` folder
holding the three pipeline outputs. `marimo export html-wasm` copies `public/`
into the export, and `mo.notebook_location()` resolves to the notebook's
directory locally and to the **served URL** in the browser, so the app fetches
its own data over HTTP.

The obvious-looking alternative — `[tool.marimo.runtime] cache_cells = true`
plus `--execute`, which bundles cell outputs into the export — **was tried and
does not work here.** Three findings worth keeping, because each cost a rebuild:

- `mo.ui.table` holds a locally defined `PandasTableManager` that **cannot be
  pickled**. Under `cache_cells` its cell output is stored as an
  `UnhashableStub`, and per marimo's cached-lifecycle contract an unserializable
  output forces a **live re-run of that cell's ancestors** — in the browser, a
  loader with no filesystem. The tables render as stub text. The plot and
  controls survive because a cell whose *defs* are live UI hits marimo's
  carve-out; a cell returning UI it does not define does not.
- `mo.persistent_cache` on just the loader is **not** a targeted substitute:
  only `cache_cells` writes the export manifest the bundler reads
  (`dump_cache_manifests` is gated on it in `_runtime/callbacks/cache.py`), so
  persistent_cache wrote to a local `__marimo__/` the exporter never saw and the
  export shipped with no data at all.
- **Pandas cannot read an `http://` URL under Pyodide** — no sockets, so urllib
  is non-functional. The loader branches on `sys.platform == "emscripten"` and
  uses `pyodide.http.open_url()`. Without this the page renders its embedded
  preview, then throws an internal error a few seconds later when Pyodide boots
  and re-runs the loader. That failure mode looks like a working app.

**The PEP 723 block installs *unpinned* `marimo`,** in the export sandbox and in
Pyodide alike — not necessarily the 0.24.0 in `py-analysis`. `mo.ui.table` was
called with `sortable=`/`filterable=` (as the `marimo-notebook` skill's `UI.md`
documents) and failed at **export runtime**, not at `marimo check`. The real
options are `show_search`, `show_column_summaries`, `max_height`,
`hover_template`. Assume any marimo API may skew between check time and run
time.

**It must be served over HTTP.** `file://` will not run WebAssembly, and the
data fetch is same-origin:
`python -m http.server --directory results/reports/landscape_explorer`.

**P1-T1's summary gained a field.** `pca_summary.json` now carries the full
30-element `variance_ratio`, not just `variance_ratio_pc1_4` — the explorer puts
any PC on either axis and needs the % for the axis label, and an app-tier
notebook may only read declared rule outputs.

---

## P1-T6 — the ADR: done

`docs/decisions/0009-mixed-models-and-the-phase-1-variance-landscape.md`. The
argument inverted exactly as anticipated: `(1|patient)` is mandatory but rests
on **non-independence**, not variance share — the design facts from P0-T3 plus
T4's cross-compartment contrast (+0.107 ρ over 102 same-patient pairs whose
compartments differ). It also records the naive-η² prohibition, `dsp_run`'s
place and its two permanent caveats, the site/compartment confounding, and the
method choices behind T3/T4.

`docs/decisions/0010-app-tier-notebooks-ship-data-in-public.md` records the
P1-T5 export mechanism.

---

## Phase 1 is complete — next is the gate, then Phase 2

**All of P1-T1 … P1-T6 are done and committed.** Nothing in Phase 1 remains.

To close the phase:

1. **Squash-merge `P1` into `main` with the Gate 1 result in the commit
   message** (ADR 0004). Gate 0's result went into `4e7a631` and PR #2 merged
   rather than squashed — do not repeat that. The message should carry:
   *Gate 1 PASSED — 41% of variance is compartment-attributable (PVCA, 14 PCs;
   per-gene median 21%; permutation null 0.3%). TIME-B n = 8.*
2. **Update `CLAUDE.md`** — it still says "current phase is Phase 1" and
   "Next: Phase 1 (A2 …)". Left untouched deliberately while Phase 1 was open.
3. **Re-decide A5 at Gate 2**, per ADR 0008. A4 stays exploratory.

Then Phase 2 (A1/A3 — the lung-vs-brain immune contrast, compartment-resolved).
Its models carry `(1|patient)` per ADR 0009, estimate the contrast **within**
compartment rather than pooled, and state `TIME-B` n = 8 inline. Remember the
power floor: ~1.1–1.3 SD at 80% (P0-T8), so a null there is uninformative, not
negative.

---

## State you'll have forgotten

- **`conda activate nsclc_bm_spatial` before every snakemake call.** The
  project standardises on that env's Snakemake **8.30**; base anaconda's 9.20
  writes incompatible provenance metadata and will provenance-trigger
  `p0t2_fetch_geo` into its protected outputs. Metadata was migrated to 8.30 on
  2026-08-24. `min_version("8.0")` will not warn you if you drift back.
- **`--use-conda` is mandatory, not optional.** Without it the software-env
  trigger fires, Snakemake tries to re-run P0-T2's downloads, and their
  `protected()` outputs abort the DAG build. A bare `snakemake -n` is not a
  clean dry run.
- **`py-analysis.conda-lock.yml` is stale.** P0-T5 added `xarray` to the env
  (see `docs/session-notes.md` for the sys.path fall-through that made it
  necessary) and `conda-lock` is not installed on this machine. A clean-room
  rebuild from the current lock fails at `p0t7_assemble_h5ad`. **P6-T1 blocker.**
- **`Markdowns/PROJECT_PLAN.md` is gitignored.** Its aims table was updated for
  A4/A5 but that edit is not version-controlled — **ADR 0008 is the record.**
- **`config/checkpoints.yaml` is still empty**, deliberately. Populating it is a
  Phase 3 stop-and-ask.
- **A4 is exploratory (ADR 0008); A5 is flagged, re-decide at Gate 2.**
- **`TIME-B` n = 8** and the power floor is **1.1–1.3 SD** at 80% (P0-T8). A null
  in Phase 2 is uninformative, not negative.
- **Gate 0's result lives in commit `4e7a631`**, not in the PR #2 merge commit —
  PR #2 was a merge rather than the squash ADR 0004 specifies. Harmless, but the
  gate result is one commit deeper than that ADR implies.
- **`.h5ad` contract:** `X` = log2(Q3+1), `layers['q3']` untransformed, 36 `obs`
  columns, git SHA + config hash in `uns`. 120 × 18,694.

## Useful commands

```bash
conda activate nsclc_bm_spatial           # first, always
snakemake -n --use-conda                  # must stay clean
snakemake --lint                          # must stay clean
snakemake --use-conda --cores 4           # build
snakemake --use-conda --report results/reports/workflow-report.html
uvx marimo check --strict notebooks/**/*.py
marimo edit notebooks/review/qc_review.py # sliders over the QC thresholds
/gate 1                                   # what Gate 1 still needs
```

If an env ever changes again:

```bash
snakemake --cleanup-metadata resources/raw/* results/interim/acquire/*.json
```

## Gate 1 — **MET** (2026-08-26), recorded in ADR 0009

> You can name the dominant variance component and point to the number.
> Downstream model specification is now defensible rather than conventional.

**Compartment, 41%** (`results/tables/variance_partition_summary.json`).
Satisfied by P1-T3, confirmed independently by P1-T4, and the model
specification it licenses is written down in ADR 0009 — so "defensible rather
than conventional" holds in the literal sense the criterion asks for.

Per ADR 0004 the gate result belongs in the **squash-merge message** for the
`P1` branch, not in a commit of its own — Gate 0's went into `4e7a631` and PR #2
merged rather than squashed, so do not repeat that here.
