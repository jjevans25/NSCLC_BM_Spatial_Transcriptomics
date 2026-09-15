# Compartment-Resolved Spatial Transcriptomics of the NSCLC Tumour Immune Microenvironment

**Primary lung tumour vs. brain metastasis.** A reproducible reanalysis of
NanoString GeoMx DSP whole-transcriptome profiles.

| | |
|---|---|
| **Project code** | `nsclc-brainmet-time` |
| **Data** | GEO [GSE200563](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200563) — 120 AOIs, 35 NSCLC patients + 7 non-tumour brain controls |
| **Source publication** | PMID [36216799](https://pubmed.ncbi.nlm.nih.gov/36216799/) |
| **Orchestration** | Snakemake 8.x |
| **Notebooks** | [marimo](https://marimo.io) only — no `.ipynb` anywhere |
| **Compliance target** | FAIR (F1–F4, A1–A2, I1–I3, R1.1–R1.3) |
| **Status** | **Phases 0–5 complete, Gates 0–5 all passed.** Phase 6 (FAIR packaging, reproduction, write-up) is next |
| **Licence** | MIT (code) · CC-BY-4.0 (derived data and figures) |

> **This is a learning and methods project, not a publication.** Every analytical
> aim below is reported with the constraint that bounds it, and several are
> explicitly exploratory. Read [How to read these results](#how-to-read-these-results)
> before quoting any number.

## What

Brain metastases are a common and poorly-served outcome in NSCLC, and the immune
microenvironment they sit in is not the one the primary tumour sits in.
GSE200563 profiled tumour, immune and glial compartments separately across
matched lung and brain sites, which makes it possible to ask the compartment
question directly rather than inferring it from bulk tissue.

This repository asks how the tumour **immune** microenvironment differs between
primary lung tumour and brain metastasis, compartment by compartment — and does
so as a fully reproducible, FAIR-packaged workflow.

## What it found

Every number below is produced by a Snakemake rule and carries its constraint.
Full detail in [`docs/analysis-notes.md`](docs/analysis-notes.md) and the
[ADRs](docs/decisions/).

**The inputs are demonstrably the published inputs.** Against the source paper's
own deposited Source Data, **all 2,243,280 values agree exactly** (max absolute
difference 0.0) and the published `NegProbe-WTX` row equals this project's
background column. The paper's 119-vs-120 AOI gap is closed by measurement: it
dropped `TBME15b`. *(P3-T2b, [ADR 0017](docs/decisions/0017-external-validation-against-the-source-publication.md))*

**Compartment, not patient or site, is what structures this data — 41% of
variance.** PVCA over 14 PCs spanning 60% of variance; per-gene median 21% across
18,691 genes; permutation null 0.3%. Independently confirmed by clustering:
compartment ARI **0.741**, patient **0.003**. This is why every downstream
contrast is estimated *within* compartment and never pooled across them.
*(Gate 1, [ADR 0009](docs/decisions/0009-mixed-models-and-the-phase-1-variance-landscape.md))*

**The published lung-vs-brain immune direction reproduces.** Brain minus lung
within the immune compartment, **`TIME-B` n = 8**: cytotoxicity **−0.95 SD
[−1.59, −0.30]** (q = 0.018); antigen presentation **−0.48 SD [−1.03, +0.08]**
(q = 0.098, and the only set detected 13/13 at *both* sites). *(Gate 2,
[ADR 0014](docs/decisions/0014-gate-2-record-and-a5-demotion.md))*

**PD-L1 is carried by the lung immune compartment**, +0.61 SD [+0.36, +0.86],
q = 0.0001, across 45 AOIs and 30 patients — the compartment-resolved analysis
bulk RNA-seq structurally cannot do. *(P3-T3,
[ADR 0018](docs/decisions/0018-primary-restriction-stands-with-a-declared-exploratory-secondary.md))*

**The dominant limitation is detection, and measuring it is the larger result.**
40 of 63 checkpoint gene × compartment cells sit below the pre-registered
detection floor; the panel clears it in `TIME-L` alone. Across 2,239
ligand–receptor interactions the assay resolves the matrix axis roughly
**eightfold** more often than the soluble one — ECM-Receptor 45.3% brain / 39.6%
lung admitted, Cell-Cell Contact 13.7% / 17.3%, Secreted Signaling **6.4% /
5.3%**. *(Gates 3–4)*

**Two honest nulls, both pre-registered and both reportable as nulls.** No
inferred ligand–receptor pair exceeded chance expectation in either adjacency
(0 of 268 lung, 0 of 260 brain). No immune signature stratified survival:
**0 of 14** tests survive FDR 0.05, smallest q = 0.26. *(Gates 4–5,
[ADR 0022](docs/decisions/0022-gate-4-record.md),
[ADR 0028](docs/decisions/0028-gate-5-record.md))*

## How to read these results

Four rules bind every number this project reports. They are enforced in code and
checked on every build by a repository-wide language audit.

- **A null here is usually uninformative, not negative.** With `TIME-B` **n = 8**
  the lung-vs-brain immune contrast detects roughly **1.1–1.3 SD** at 80% power.
  The ligand–receptor and survival nulls are **assay-sensitivity and cohort-size
  limits**, never evidence that the biology is absent. Phase 5's two nulls are
  not even the same kind: this cohort fitted 12 and 7 patients, TCGA fitted 502.
- **"Not assessable" is not "lower".** Background is higher in brain, so a
  near-background gene reads as depleted there *artefactually*. A gene or
  signature below the detection floor is reported as **not assessable in that
  compartment** and never as a difference.
  *([ADR 0008](docs/decisions/0008-demote-a4-to-exploratory.md))*
- **A compartment label is not a cell-type label.** PanCK was the *only*
  collection mask. A `TIME` AOI is the PanCK-negative segment of an ROI sited in
  a CD45-rich region — **not** a CD45-sorted population. Never "CD45+ AOI".
- **This project does not claim colocalisation.** The assay has no coordinates,
  so the claim it makes is **inferred crosstalk between adjacent compartments**,
  and nothing stronger.

Two aims are **exploratory by decision, not by outcome**: the checkpoint
landscape ([ADR 0008](docs/decisions/0008-demote-a4-to-exploratory.md)) and the
crosstalk analysis ([ADR 0014](docs/decisions/0014-gate-2-record-and-a5-demotion.md)).
No result from either may be a headline claim.

## Why it is built this way

Three architectural commitments, each a decision rather than a default:

- **Snakemake is the only source of truth.** Every reported number, table and
  figure is produced by a rule. Deleting `results/` and re-running must
  reproduce everything; if that is not true, the workflow is lying about
  provenance.
- **Notebooks are marimo, and only marimo.** Pure `.py`, diffable, no hidden
  state, deterministic execution order, and runnable as a script by a rule — so
  the exploratory artifact and the pipeline step never silently diverge. See
  [`notebooks/README.md`](notebooks/README.md).
- **The R/Python boundary is a file boundary.** R owns the assay-specific work
  (`GeomxTools`, `standR`, `SpatialDecon`, mixed models); Python owns AnnData,
  plotting, enrichment, survival and packaging. They exchange plain TSV + JSON,
  not a bridge library
  ([ADR 0001](docs/decisions/0001-anndata-handoff-without-zellkonverter.md)).

A fourth commitment is procedural: **every threshold is pre-registered in an ADR
committed before the result it governs.** A threshold chosen after seeing which
genes clear it is not a threshold.

## How to run

```bash
# 1. Bootstrap the driver environment (Snakemake + Python)
conda env create -f environment.yml
conda activate nsclc_bm_spatial

# 2. Point the per-rule conda envs at a fixed prefix (once per machine)
ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"

# 3. Check the workflow parses and the DAG is sane
snakemake --lint
snakemake -n --use-conda --conda-prefix "$HOME/nsclc-envs"

# 4. Run — per-rule conda envs are created on first use
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" --cores 4

# 5. Build the provenance report
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" \
  --report results/reports/workflow-report.html
```

> **`--use-conda` and `--conda-prefix` are not optional.** Without them the
> software-environment rerun trigger fires, Snakemake tries to re-run the
> acquisition rules, and their `protected()` outputs abort the DAG. **A bare
> `snakemake -n` is not a clean dry run.** Run from the `nsclc_bm_spatial`
> environment specifically — a newer Snakemake writes `.snakemake/metadata` in a
> format 8.30 reads as stale, with the same result.
> See [ADR 0013](docs/decisions/0013-r-on-a-path-with-spaces-and-the-float-boundary.md).

Target rules, so you can run a slice without editing the Snakefile:

| Target | Runs |
|---|---|
| `all` | everything enabled in `config/config.yaml → phases` |
| `qc_only` | acquisition → annotation → QC (Phase 0) |
| `phase2` | through the immune-contexture analysis |
| `phase3` | through the checkpoint landscape |
| `phase4` | through the inferred-crosstalk analysis |
| `fair` | everything, including RO-Crate packaging |

A phase only contributes targets when **both** its config switch is on and its
`TARGETS_*` list in `workflow/rules/*.smk` is populated. That is deliberate: a
partially-built workflow reports an empty DAG rather than a missing-input error.

## Layout

```
config/          paths, thresholds, switches; samples.tsv; per-source manifests
workflow/        Snakefile, rules/, scripts/, envs/, schemas/, report/
resources/       READ-ONLY. Acquired inputs + checksums. Never written by a rule
results/         everything derived: interim/, tables/, figures/, logs/, reports/
notebooks/       marimo only, three tiers — explore/, review/, apps/
metadata/        RO-Crate, ontology terms, dataset description
docs/            data provenance, ADRs, analysis notes, limitations
```

`resources/` is read-only and `results/` is fully derived. Those two facts
together are what make the reproduction claim checkable rather than aspirational.

It has **six sanctioned writers**, each with a pinned digest and a
`*_checksums.sha256` / `*_provenance.tsv` pair, and nothing else may write there:

| Directory | Source | Task |
|---|---|---|
| `raw/` | GEO GSE200563 | P0-T2 |
| `reference/` | SpatialDecon cell profiles | P2-T5 |
| `supplementary/` | the source publication's supplementary files | P3-T2b |
| `ligand_receptor/` | CellChatDB v2, commit-pinned | P4-T2 |
| `clinical/` | Supplementary Data 1 (clinical + survival) | P5-T1 |
| `tcga/` | TCGA LUAD via UCSC Xena | P5-T4 |

## Status

**Phases 0–5 are complete and Gates 0–5 have all passed.** Phase 6 — FAIR
packaging, clean-room reproduction and write-up — is next; PROJECT_PLAN calls it
the primary deliverable.

All five open questions from the plan are **resolved**, with citations, in
[`docs/data-provenance.md`](docs/data-provenance.md): the GEO matrix is
Q3-normalised (Q1), the AOIs are antibody-segmented (Q2), slide/batch is
recoverable from the DCC filenames (Q3), per-AOI nuclei counts did not survive
into GEO (Q4), and patient-level survival metadata is joinable (Q5).

What the design cannot support is recorded in
[`docs/limitations.md`](docs/limitations.md), and it is worth reading before the
results.

## Citing

Cite this workflow via [`CITATION.cff`](CITATION.cff). Cite the **data**
separately: GEO GSE200563, and PMID 36216799 for the original study. This
repository makes no claim over the upstream data — see
[`LICENSE-DATA`](LICENSE-DATA).

## Plan and decision record

The full six-phase plan with per-task acceptance criteria and phase gates lives
in `Markdowns/PROJECT_PLAN.md`, which is **gitignored and therefore not in this
repository**. The durable, public record is
[`docs/decisions/`](docs/decisions/) — every threshold, scope change and gate
verdict, each committed before the result it governs — together with
[`CHANGELOG.md`](CHANGELOG.md) and
[`docs/NEXT_STEPS.md`](docs/NEXT_STEPS.md).
