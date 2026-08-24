# ADR 0006 — P0-T5 runs in Python, because GeomxTools cannot run at all

**Date:** 2026-08-24
**Status:** Accepted
**Task:** P0-T5
**Supersedes:** the tool assignment for P0-T5 in `Markdowns/PROJECT_PLAN.md` §6

## Context

`CLAUDE.md` draws a language boundary: "R owns the assay-specific work:
`GeomxTools`, `standR`, `SpatialDecon`, `limma`/`edgeR`, and the mixed models
… **Do not reimplement `standR` or `GeomxTools` in Python.**" §6 accordingly
lists P0-T5's tools as `standR`, `GeomxTools`, `marimo-pair`.

Both packages are built around a `NanoStringGeoMxSet`, and
`readNanoStringGeoMxSet()` requires three things: the DCC files, a **PKC**
probe-kit configuration, and a phenotype worksheet.

We have the DCCs. We do not have a PKC, and cannot get one on terms this
project can accept:

- **GEO does not carry it.** GSE200563's supplementary files are the processed
  matrix, `GSE200563_RAW.tar` (120 `.dcc.gz`) and `filelist.txt`. There is no
  `.pkc` and no annotation worksheet (`docs/data-provenance.md`, "Known gaps").
- **It is not publicly fetchable.** `Hs_R00000001_WTA.pkc` is distributed
  through NanoString/Bruker's resource portal behind registration, not from a
  stable URL that P0-T2 could pin with a SHA-256. Adding an unpinnable,
  login-gated file to a read-only `resources/` would put a hole straight
  through the provenance chain that ADR 0005 exists to close.

Without a PKC the DCC `<Code_Summary>` block — raw counts keyed by RTS probe
ids — cannot be mapped to genes. So there is no GeoMxSet, no probe-level
outlier removal, no per-target LOQ, and no raw gene-level count matrix.

The alternatives considered:

1. **Plain R in the `r-geomx` env.** Honours the boundary literally. But it
   builds a large Bioconductor environment for packages that are then never
   called, and the metrics computed would be identical arithmetic in a
   different language.
2. **Pause P0-T5** until a PKC is obtained by hand. Blocks Gate 0 and Phase 1
   on an external dependency with no timeline, to recover analyses (probe
   outlier removal) that the already-Q3-normalised deposited matrix cannot
   benefit from anyway.
3. **Python.** Chosen.

## Decision

**P0-T5 is a Python rule** (`workflow/scripts/qc_metrics.py`, `py-analysis`).

The boundary in `CLAUDE.md` is re-read as what it was always for: *do not
rewrite the assay packages' logic in Python*. That prohibition is not engaged
here, because none of that logic is available to rewrite. What P0-T5 computes
is arithmetic on values already in hand:

- **Sequencing metrics from the DCC headers** — `Raw`, `Trimmed`, `Stitched`,
  `Aligned`, `umiQ30`, `rtsQ30`, and the deduplicated total implied by
  `<Code_Summary>`. These need no PKC. Saturation is `1 − dedup/aligned`.
- **Detection from the processed matrix**, relative to `NegProbe-WTX`.
- **Normalisation verification, not application** (Q1; ADR-free because the
  data forces it — raw gene-level counts do not exist for us).

**R is not abandoned.** It re-enters at Phase 2 for `lme4`/`lmerTest`/`emmeans`,
where the patient random intercept is mandatory (§2.3) and R genuinely has no
Python equivalent worth using. `workflow/envs/r-stats.yaml` stands unchanged.
`workflow/envs/r-geomx.yaml` is kept but currently unused; it becomes live only
if a PKC is ever obtained.

## Consequences

What this costs:

- **No probe-level QC.** `GeomxTools` would flag outlier probes within a target
  and drop them before summarisation. We inherit the submitters' summarisation
  as given, with no visibility into it. State this in `docs/limitations.md`
  (P0-T8).
- **No per-target LOQ.** The standard GeoMx limit of quantitation is
  `geomean(negative probes) × GeoSD(negative probes)²`. The deposited matrix has
  **one** `NegProbe-WTX` row, so there is no SD and no geomean across probes.
  The 2× background multiple in ADR 0007 is a stand-in for a statistic we
  cannot compute.
- **A stated project constraint is now weaker.** "R owns the assay-specific
  work" survives only in Phase 2. A future reader will find `r-geomx.yaml` in
  the tree with nothing using it; this ADR is why.
- **The boundary is now a judgement call, not a rule.** That is a real loss of
  crispness, and the mitigation is that this ADR names the exact condition that
  released it — no PKC, therefore no GeoMxSet — rather than a general licence.

What it buys:

- Gate 0 and Phase 1 are not blocked on a login-gated file.
- `resources/` keeps its property that every byte has a pinned checksum and a
  public URL (ADR 0005).
- One conda environment instead of two for Phase 0, so the clean-room
  reproduction at P6-T1 has less to rebuild.

**Reversible.** If a PKC is obtained, add it to the P0-T2 manifest as a `fixed`
artifact, write the R rule, and supersede this ADR. Nothing here forecloses it.
