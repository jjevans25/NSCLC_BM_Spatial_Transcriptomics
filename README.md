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
| **Status** | Phase 0. Learning / methods project, not intended for peer-reviewed publication |
| **Licence** | MIT (code) · CC-BY-4.0 (derived data and figures) |

## What

Brain metastases are a common and poorly-served outcome in NSCLC, and the
immune microenvironment they sit in is not the one the primary tumour sits in.
GSE200563 profiled tumour, immune and glial compartments separately across
matched lung and brain sites, which makes it possible to ask the compartment
question directly rather than inferring it from bulk tissue.

This repository asks how the tumour **immune** microenvironment differs between
primary lung tumour and brain metastasis, compartment by compartment — and does
so as a fully reproducible, FAIR-packaged workflow.

## Why it is built this way

Three architectural commitments, each of which is a decision rather than a
default:

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
  plotting, enrichment and packaging. They exchange plain TSV + JSON, not a
  bridge library ([ADR 0001](docs/decisions/0001-anndata-handoff-without-zellkonverter.md)).

## How to run

```bash
# 1. Bootstrap environment (Snakemake + Python)
conda env create -f environment.yml
conda activate nsclc_bm_spatial

# 2. Check the workflow parses and the DAG is sane
snakemake --lint
snakemake -n

# 3. Run — per-rule conda envs are created on first use
snakemake --software-deployment-method conda --cores 4

# 4. Build the provenance report
snakemake --report results/reports/workflow-report.html
```

Target rules, so you can run a slice without editing the Snakefile:

| Target | Runs |
|---|---|
| `all` | everything enabled in `config/config.yaml → phases` |
| `qc_only` | acquisition → annotation → QC (Phase 0) |
| `phase2` | through the immune-contexture analysis |
| `phase3` | through the checkpoint landscape |
| `fair` | everything, including RO-Crate packaging |

A phase only contributes targets when **both** its config switch is on and its
`TARGETS_*` list in `workflow/rules/*.smk` is populated. That is deliberate: a
partially-built workflow reports an empty DAG rather than a missing-input error.

## Layout

```
config/          paths, thresholds, switches; samples.tsv; compartment map
workflow/        Snakefile, rules/, scripts/, envs/, schemas/, report/
resources/       READ-ONLY. Raw GEO downloads + checksums. Never written by a rule
results/         everything derived: interim/, tables/, figures/, logs/, reports/
notebooks/       marimo only, three tiers — explore/, review/, apps/
metadata/        RO-Crate, ontology terms, dataset description
docs/            data provenance, ADRs, analysis notes, limitations
.claude/         agent hooks, commands, project-local skills
```

`resources/` is read-only and `results/` is fully derived. Those two facts
together are what make the reproduction claim checkable rather than aspirational.

## Status and open questions

Phase 0 is in progress. Five questions from the plan (§3) are unresolved and two
of them can still change Phase 0 — most importantly **Q2**, whether the AOIs
were antibody-segmented or are geometric ROIs with descriptive names. Answers
land in `docs/data-provenance.md`.

## Citing

Cite this workflow via [`CITATION.cff`](CITATION.cff). Cite the **data**
separately: GEO GSE200563, and PMID 36216799 for the original study. This
repository makes no claim over the upstream data — see
[`LICENSE-DATA`](LICENSE-DATA).

## Plan

The full six-phase plan, with per-task acceptance criteria and phase gates,
is in `Markdowns/PROJECT_PLAN.md`.
