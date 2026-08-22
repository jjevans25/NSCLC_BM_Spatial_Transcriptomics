# 0001 — The R→Python handoff writes TSV, not `.h5ad` via zellkonverter

- **Status:** accepted
- **Date:** 2026-08-22
- **Affects:** PROJECT_PLAN.md §4.1 (boundary object), §4.4 (`r-geomx.yaml`),
  §4.2 (`workflow/scripts/to_anndata.R`), P0-T7
- **Supersedes:** the §4.1 sentence "write a single `.h5ad` (AnnData) via `zellkonverter`"

## Context

§4.1 makes a single `.h5ad` the boundary object between the R half (Phase 0 QC and
normalisation) and the Python half (Phases 1, 3–5), and nominates `zellkonverter`
as the writer. §4.4 therefore lists `bioconductor-zellkonverter` in `r-geomx.yaml`.

A dry-run solve of `r-geomx.yaml` on `spike/env-feasibility` showed what that
actually pulls in, and what it does not:

```
bioconductor-basilisk        1.18.0     present
bioconductor-basilisk.utils  1.18.0     present
r-reticulate                 1.46.0     present

anndata                                 ABSENT
h5py                                    ABSENT
scipy                                   ABSENT
pandas                                  ABSENT
```

zellkonverter does not talk to a Python interpreter that conda installed. It goes
through basilisk, which provisions its **own private Python environment at run
time**, on first use.

That breaks the workflow's central provenance claim three ways:

1. **conda-lock cannot pin it.** §4.4's whole argument is that unpinned
   environments are the most common reason a reanalysis stops reproducing. The
   basilisk env is invisible to the lock file, so the locks would be lying.
2. **It needs the network mid-run.** P6-T1 is a clean-room reproduction treated
   as a gate, not a formality. A rule that downloads a Python environment while
   executing can pass on this laptop in August and fail on a fresh checkout in
   October, for reasons no committed file records.
3. **Snakemake cannot see it.** §4.3 requires every rule to declare its
   environment via `conda:`. A runtime-provisioned interpreter is an undeclared
   dependency — the same class of provenance hole that §4.5 rejects `.ipynb` for.

## Decision

**Drop `zellkonverter`. The boundary object stays `.h5ad`; only the writer changes.**

The handoff becomes two rules across the language boundary:

| Step | Env | Writes |
|---|---|---|
| R export | `r-geomx.yaml` | `results/interim/expr_normalised.tsv`, `obs.tsv`, `var.tsv`, `uns.json` |
| Python assembly | `py-analysis.yaml` | `results/interim/geomx.h5ad` |

`workflow/scripts/to_anndata.R` is replaced by an R export step (plain
`write.table`, no new dependency) plus `workflow/scripts/to_anndata.py`, which
builds the AnnData with the `anndata` package already listed in §4.4's
`py-analysis.yaml`.

P0-T7's acceptance criteria are unchanged and still met: the `uns.json` carries
the git SHA, config hash, and run date; `obs.tsv` carries the full annotation and
the §7.2 ontology IDs, and must still match `config/samples.tsv`.

## Consequences

**Good**

- `r-geomx.yaml` loses basilisk, basilisk.utils, and reticulate. Every remaining
  dependency is declared and lockable. Verified: re-solved clean at 387 packages
  (was 395), with all six remaining Bioconductor packages at identical versions.
- The boundary is a plain text table. It can be inspected with `head`, diffed,
  checksummed into `resources/checksums.sha256`, and read by any tool — which is
  a better FAIR-I story than a binary written by an R shim around a hidden Python.
- The failure mode moves to the left. A malformed TSV fails immediately and
  legibly; a basilisk provisioning failure surfaces as an opaque reticulate error
  deep in a Snakemake job.

**Bad**

- One extra rule and one intermediate artifact in the DAG.
- TSV round-trips floats as text. Write with `digits = 17` so the values are
  exactly recoverable, and assert equality in the Python step rather than
  assuming it.
- Sparse matrices densify in TSV. Irrelevant here — 120 AOIs × ~18k genes of
  Q3-normalised WTA data is dense and small (~50 MB uncompressed, gzip it).

**Neutral**

- `bioconductor-spatialdecon` still pulls Python 3.14, numpy, and yq into
  `r-geomx.yaml`. That is a bioconda recipe dependency, fully declared in the
  solve and pinned by conda-lock — it is not the basilisk problem and needs no
  action.

## Alternatives rejected

- **Keep zellkonverter, add the Python stack to `r-geomx.yaml` and force
  reticulate onto the ambient interpreter.** Preserves §4.1 as written, but the
  basilisk override mechanism would have to be verified and then re-verified on
  every zellkonverter bump. It buys a binary boundary object at the cost of a
  standing reproducibility risk.
- **Hand off as `.rds` and read with `rpy2`.** Trades an R-side Python dependency
  for a Python-side R dependency. Strictly worse.
- **Stay in R for everything.** Contradicts §4.1, and discards the K-Dense Python
  skills that motivate the hybrid split in the first place.
