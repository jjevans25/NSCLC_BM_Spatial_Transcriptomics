# Next steps — Phase 2 (immune contexture, brain vs lung)

**Last session:** 2026-08-26
**Branch:** `main` at `0d3fa9c`, **pushed**. Phase 2 needs a new `P2` branch.
**Gate 0:** PASSED. **Gate 1:** PASSED — squash-merged with the result in the
commit message, as ADR 0004 asks.
**Phase 1:** **COMPLETE** (P1-T1 … P1-T6). Nothing outstanding.

Phase 2's goal (§6): Aims **A1/A3** — is the brain-metastasis immune
microenvironment different from the primary lung tumour, and in which
direction, resolved by compartment.

---

## Gate 2 — read this before writing any code

> **Do you recover the published direction (reduced antigen presentation and
> B/T function in brain)?** If yes: pipeline validated, proceed with confidence.
> If no: **stop and debug** — an unexpected result at this n is far more likely
> to be a pipeline bug than new biology. **Resist the opposite instinct.**

This is the only gate in the project whose criterion is "reproduce someone
else's answer". A1–A3 are the **positive control on the pipeline**, not a
novelty claim (§1.3). Treat a surprising Phase 2 result as a bug report against
yourself.

**A5 is re-decided at this gate** (ADR 0008). It inherits `TBME`'s batch
confounding and 21.4% detection. Decide it explicitly and write it down; do not
let it drift into Phase 4 undecided.

---

## What Phase 1 decided that binds Phase 2

**ADR 0009 is the model contract.** Three things follow, and none is optional:

1. **Every model carries `(1|patient)`** — justified on non-independence, never
   on variance share. Patient is only 12% of variance and its clustering ARI is
   0.003; the load-bearing evidence is the +0.107 ρ contrast across the 102
   same-patient AOI pairs whose compartments differ. Changing the
   random-effects structure is a **stop and ask**.
2. **Estimate within compartment, never pooled across compartments.**
   Compartment is 41% of variance — pooling `TIME` with `LB`/`L` would swamp the
   site effect you are trying to measure.
3. **No variance share or p-value without its null / n / CI / effect size.**
   Hard constraint 7, and §2 of ADR 0009 for why a naive η² is unusable here.

---

## Three things measured this session that change Phase 2's shape

**a. The `TIME` AOIs are QC-clean.** Zero of the 23 `TIME-L` + `TIME-B` AOIs are
QC-flagged. The flagging in P0-T5 hit 3 `TBME` AOIs, none of them here — so
`TIME-B` n = 8 is intact and Phase 2 starts at full strength. Design every
analysis so it still degrades gracefully, but the feared degradation did not
happen.

**b. Batch is estimable in the `TIME` contrast, but barely.**

| | run A | run B |
|---|---|---|
| `TIME-L` | 9 | 6 |
| `TIME-B` | **6** | **2** |

Unlike `mLN`/`TBME`/`BC` — each wholly inside one run — the `TIME` comparison
does span both. But `TIME-B` has **2 AOIs in run B**, so `score ~ site +
dsp_run + (1|patient)` is fragile and will likely fit singularly.

**P2-T3 must decide this explicitly, and it is a stop-and-ask.** The plan
specifies `score ~ site + (1|patient)` with no batch term. P1-T3 found `dsp_run`
carries ~10% of variance, which argues for adjusting; the 6/2 split argues
against. Suggested resolution: fit the plan's model as primary, refit with
`dsp_run` as a **sensitivity**, and report both — the same shape P1-T3's S2 used,
which is what let it quantify batch absorption at +0.014 rather than assert it.
Record the choice in an ADR (next free number: **0011**).

**c. Five patients have both `TIME-L` and `TIME-B`** — confirmed against the
data, matching §2.3. That is P2-T4's paired check, and it is a
direction-of-effect consistency check only. **No p-values from n = 5.**

---

## Start here

### 1. P2-T1 — version the signatures (~2 h)

`config/signatures.yaml` is currently `signatures: {}`. Populate: cytotoxicity
(`GZMB, PRF1, GNLY, NKG7`), exhaustion (`PDCD1, LAG3, HAVCR2, TIGIT, CTLA4,
TOX`), antigen presentation (HLA-I/II, `B2M, TAP1, TAP2, NLRC5`), M1/M2 myeloid,
TLS (`CXCL13, CCL19, CCL21, CR2, MS4A1`).

- **Every set carries a `source:`** (MSigDB ID or PMID). Un-sourced gene sets are
  how reanalyses become unreproducible. **Accept:** YAML validates; every set
  sourced.
- **There is no signatures schema yet** — `workflow/schemas/` holds only
  `config.schema.yaml` and `samples.schema.yaml`. Write
  `signatures.schema.yaml` and validate at Snakefile load, the way `samples.tsv`
  is validated in `common.smk`.
- **Adding or removing a gene from a signature is a stop-and-ask**
  (`CLAUDE.md`), including while first populating the file.
- Antigen presentation is the set Gate 2 turns on — the published direction is
  *reduced* antigen presentation in brain. Get its membership and source right
  before anything else.

### 2. P2-T2 — score signatures (~3 h)

ssGSEA (`gseapy`, already in `py-analysis`) **and** a z-score mean, on the
`TIME-L` + `TIME-B` AOIs (23 AOIs). Two methods because **agreement between them
is the robustness evidence at this n**.

**Record per-set detection coverage.** A set with 2 of 6 genes detected is not a
signature, and P0-T5 established that detection varies systematically by
compartment — so coverage must be reported per site, not pooled. Reuse the
`detection_at_2x` columns already in `obs` and `config.qc.detection_background_multiple`
rather than inventing a second detection rule.

### 3. P2-T3 — mixed models (~3 h)

`score ~ site + (1|patient)` per signature in `lme4`/`lmerTest`, cross-checked on
one signature in `statsmodels.MixedLM`. Report estimate, 95% CI, df method, raw
p, BH-FDR across signatures. **Singular fits surfaced, not hidden.**

See §b above — the `dsp_run` decision belongs here, with an ADR.

### 4. P2-T4 — paired sensitivity (~1.5 h)

The 5 patients with both `TIME-L` and `TIME-B`. **Direction of effect only.**
Per-signature concordance table against the unpaired result.

### 5. P2-T5 — deconvolution (~3 h)

`SpatialDecon`, two-matrix: lung reference for `TIME-L`, brain reference for
`TIME-B`. Renormalise composition **within** compartment. **Accept:**
stacked-bar figure **plus a written caveat that cross-site comparison is
reference-confounded** — that caveat is part of the deliverable, not a footnote.

### 6. P2-T6 — convergence check (~1 h)

Do deconvolution and signatures agree on direction? Where they disagree, write
which you trust and why (usually signatures — they do not depend on reference
choice). **Out:** paragraph in `docs/analysis-notes.md`.

### 7. P2-T7 — deliverable figure (~2 h)

Heatmap: signatures × AOIs, columns grouped by site, annotated by patient.
Publication-grade, colourblind-safe, in the Snakemake report.

---

## State you'll have forgotten

- **Neither R environment has ever been built.** Only two `py-analysis` envs
  exist under `.snakemake/conda/`. The first Phase 2 rule with
  `conda: "../envs/r-stats.yaml"` triggers a fresh solve, and `r-geomx.yaml`
  pulls **five Bioconductor packages** (`GeomxTools`, `standR`, `SpatialDecon`,
  `limma`, `edgeR`). Budget real time for that solve, and do it *early* rather
  than discovering it at P2-T5. It is the largest un-derisked step in Phase 2.
- **`conda activate nsclc_bm_spatial` before every snakemake call.** That env's
  Snakemake is **8.30**; base anaconda's 9.20 writes provenance metadata 8.30
  reads as stale and will re-trigger `p0t2_fetch_geo` into its protected
  outputs. `min_version("8.0")` will not warn you.
- **`--use-conda` is mandatory, not optional.** Without it the software-env
  trigger fires and P0-T2's `protected()` outputs abort the DAG build. A bare
  `snakemake -n` is not a clean dry run.
- **Editing a `script:` rule can leave `snakemake -n` dirty** in a way only a
  re-run clears. Snakemake 8.30 stores `code: None` for every `script:` rule
  here, yet after this session's edits three rules reported "code has changed".
  A full re-run cleared it and reproduced every output byte-identically. Not
  understood; do not reach for `--cleanup-metadata` (it makes snakemake *forget*
  provenance, which is wrong for a FAIR target) unless an env genuinely changed.
- **`py-analysis.conda-lock.yml` is stale** and `conda-lock` is not installed on
  this machine. P0-T5 added `xarray`; a clean-room rebuild fails at
  `p0t7_assemble_h5ad`. **Still the P6-T1 blocker.** Nothing this session made it
  worse — P1-T4 implements ARI in numpy rather than adding scikit-learn, and
  P1-T5 dropped its pyarrow dependency.
- **App-tier notebooks: read ADR 0010 before writing P3-T6's
  `checkpoint_explorer.py`.** It ships data in `public/` and fetches over HTTP;
  `cache_cells` breaks `mo.ui.table`, `mo.persistent_cache` silently ships an
  export with no data, and pandas cannot read an `http://` URL under Pyodide.
  Also: `marimo check` passes while the export is broken, and the PEP 723 block
  installs *unpinned* marimo, so APIs skew between check time and run time.
  **Verify app-tier notebooks by exporting and opening them, never by linting.**
- **`Markdowns/PROJECT_PLAN.md` is gitignored.** Its aims table was updated for
  A4/A5 but that edit is not version-controlled — **ADR 0008 is the record.**
- **`config/checkpoints.yaml` is still empty**, deliberately. Populating it is a
  Phase 3 stop-and-ask.
- **A4 is exploratory (ADR 0008); A5 is re-decided at Gate 2** — i.e. now.
- **`TIME-B` n = 8** and the power floor is **1.1–1.3 SD** at 80% (P0-T8). **A
  null in Phase 2 is uninformative, not negative** — and per Gate 2, a
  *surprising* result is more likely a bug than biology.
- **`.h5ad` contract:** `X` = log2(Q3+1), `layers['q3']` untransformed, 36 `obs`
  columns, git SHA + config hash in `uns`. 120 × 18,694.
- **`config/config.yaml` is a declared input of `p0t7_assemble_h5ad`**, so any
  config edit rebuilds the `.h5ad` and re-runs everything downstream — including
  P1-T3's ~25-minute fit. Batch config changes rather than making them one at a
  time.

## Conventions that bite if forgotten

- **Turn `phases.contexture: true` in the same commit that adds the rules**, and
  populate `TARGETS_CONTEXTURE` in `04_contexture.smk` (currently `[]`) — a
  phase contributes to `rule all` only when both are true.
- One `.smk` per phase; every rule declares `log:`, `benchmark:`, `conda:` and
  `threads:`. No exceptions.
- **R owns** `lme4`/`lmerTest`/`emmeans` and `SpatialDecon`; **Python owns**
  AnnData, plotting, enrichment. The handoff is plain TSV with `digits = 17`
  and an asserted round-trip — **not** `zellkonverter` (ADR 0001).
- Commits: `P2-T<n>: imperative summary`. One phase per branch, squash-merged
  with the gate result in the message (ADR 0004).
- Never write "colocalisation" — **"inferred crosstalk between adjacent
  compartments"**, no exceptions.
- Never write "CD45+ AOI". A `TIME` AOI is the PanCK-negative segment of an ROI
  sited in a CD45-rich region — **not** a CD45-sorted population (Q2's caveat).

## Useful commands

```bash
conda activate nsclc_bm_spatial           # first, always
git checkout -b P2                        # one branch per phase (ADR 0004)
snakemake -n --use-conda                  # must stay clean
snakemake --lint                          # must stay clean
snakemake --use-conda --cores 4           # build
snakemake --use-conda --conda-create-envs-only   # build the R envs early
snakemake --use-conda --report results/reports/workflow-report.html
uvx marimo check --strict notebooks/**/*.py
/gate 2                                   # what Gate 2 still needs
```

Serving the Phase 1 explorer (it needs HTTP, not `file://`):

```bash
python -m http.server --directory results/reports/landscape_explorer
```

---

## Phase 1, for reference

Gate 1: **41% of variance is compartment-attributable** (PVCA over 14 PCs
spanning 60.4%; per-gene median 21% over 18,691 converged genes; permutation
null 0.3%). `TIME-B` n = 8.

| Component | Per-gene median | PVCA | Permutation null |
|---|---|---|---|
| compartment | 0.211 | **0.412** | 0.003 |
| patient | 0.124 | 0.218 | 0.026 |
| `dsp_run` | 0.078 | 0.096 | 0.002 |
| site | 0.006 | 0.023 | 0.003 |
| residual | 0.495 | 0.252 | 0.889 |

P1-T4 confirmed it from the clustering side: compartment ARI **0.741**, patient
**0.003**, same-patient AOIs sharing a cluster at *below* the chance rate —
inverting PROJECT_PLAN §6's stated expectation, which is recorded as measured
false in ADR 0009 so it is not reinstated from the plan.

Artefacts: `variance_partition.tsv`, `variance_partition_summary.json`,
`cluster_agreement.tsv`, `sample_clustering_summary.json`,
`results/figures/variance_partition.png`,
`results/figures/sample_correlation_heatmap.png`,
`results/reports/landscape_explorer/`. ADRs **0009** (models) and **0010**
(app-tier notebook exports).
