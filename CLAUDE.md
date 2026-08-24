# CLAUDE.md — operating instructions for this repository

## The project, in a paragraph

`nsclc-brainmet-time` reanalyses NanoString GeoMx DSP whole-transcriptome
profiles from GEO **GSE200563** (120 AOIs; 35 NSCLC patients + 7 non-tumour
brain controls; source publication **PMID 36216799**) to ask how the tumour
**immune** microenvironment differs between a primary lung tumour and a brain
metastasis, resolved by tissue compartment. Orchestration is **Snakemake 8.x**;
notebooks are **marimo only**; the compliance target is **FAIR**. This is a
learning and methods project, not a publication. The full plan lives in
`Markdowns/PROJECT_PLAN.md` — when this file and the plan disagree, the plan
wins, and the disagreement is a bug in this file.

**GATE 0 PASSED (2026-08-24). Phase 0 complete (P0-T1 … P0-T8); current phase
is Phase 1 (variance landscape).** All five open questions are answered
(`docs/data-provenance.md`), and `docs/limitations.md` records what the design
cannot support.

Gate 0 evidence: marker sanity check 4/4
(`results/tables/marker_sanity_verdict.tsv`), design table matches §2.1 exactly
(120 AOIs), Q2 resolved to antibody segmentation.

**Scope change at the gate: Aim A4 is demoted to exploratory (ADR 0008).** Most
of the checkpoint panel sits at background — CTLA4 is above background in 2 of
8 `TIME-B` AOIs, TIGIT in 1, IDO1 in 2 — and brain background is *higher*, so a
near-background gene reads as depleted in brain artefactually. Phase 3 still
runs; **no A4 result may be a headline claim**, every checkpoint reported states
its per-site detection count, and a gene detected in fewer than half the
`TIME-B` AOIs is "not assessable in brain", never "lower in brain".
**A5 is flagged, not demoted** — re-decide at Gate 2.

Phase 0 output: `results/interim/aoi_normalised.h5ad` — 120 AOIs × 18,694
genes, `X` = log2(Q3 + 1), `layers['q3']` the untransformed values, 36 `obs`
columns, git SHA and config hash in `uns`.

**The number to keep in view:** with `TIME-B` n = 8, the lung-vs-brain immune
contrast detects roughly **1.1–1.3 SD** at 80% power (P0-T8). A null result in
Phase 2 is uninformative, not negative.

Next: Phase 1 (A2 — what drives variance: patient vs site vs compartment).

## The design table — use these numbers, never estimates

| Code | Site | Compartment (as labelled) | AOIs | Patients |
|---|---|---|---|---|
| `L` | Lung | Primary tumour | 30 | 30 |
| `LB` | Brain | Brain-metastasis tumour | 27 | 27 |
| `mLN` | Lymph node | Nodal metastasis tumour | 13 | 13 |
| `TBME` | Brain | Tumour–brain microenvironment (glial) | 20 | 19 |
| `TIME-L` | Lung | Tumour immune microenvironment | 15 | 13 |
| `TIME-B` | Brain | Tumour immune microenvironment | 8 | 8 |
| `BC` | Brain | Non-tumour brain control | 7 | — |
| | | **Total** | **120** | **35 + controls** |

Structural facts with statistical consequences:

- **AOIs within a patient are not independent.** P12 and P24 each have two
  `TIME-L`; P15 has two `TBME`. A random intercept for patient is **mandatory**,
  not stylistic.
- **`TIME-B` n = 8** is the binding constraint on the entire project. Design
  every analysis so it degrades gracefully if 1–2 of those AOIs fail QC.
- Only **5 patients** (5, 12, 15, 19, 35) have paired immune AOIs at both sites.
  That is a consistency check, never a headline result.

## Hard constraints

1. **`resources/` is read-only.** No rule, script or notebook writes into it.
2. **Notebooks are marimo (`.py`) only. Never create a `.ipynb`.** If you
   encounter one, convert it first: `marimo convert old.ipynb -o notebooks/explore/<name>.py`.
   A `PreToolUse` hook blocks `.ipynb` writes (ADR 0002) — do not try to route
   around it.
3. **No reported number, table, or figure may originate in a notebook.**
   Notebooks read pipeline outputs; `workflow/scripts/` produces them. If a
   notebook exploration matters, promote it to a script + rule *first*, then
   have the notebook read that rule's output.
4. **No analysis code outside `workflow/scripts/`.** Notebooks under
   `notebooks/apps/` may contain presentation logic only.
5. **Every new output is produced by a Snakemake rule**, never by an ad-hoc
   script run. If you ran something by hand to get a file, it is not a result
   yet.
6. **Never write "colocalisation."** Write "inferred crosstalk between adjacent
   compartments." The assay does not support the stronger claim.
7. **Never report a p-value without n, effect size, and a confidence interval.**
8. **Any claim about brain immune contexture states `TIME-B` n = 8 inline.**

## The R/Python boundary

R owns the assay-specific work: `GeomxTools`, `standR`, `SpatialDecon`,
`limma`/`edgeR`, and the mixed models (`lme4`, `lmerTest`, `emmeans`). Python
owns AnnData, `scanpy`, plotting, enrichment, survival and packaging. **Do not
reimplement `standR` or `GeomxTools` in Python.** The handoff crosses the
boundary as plain files, not `zellkonverter` (ADR 0001): R writes
`expr_normalised.tsv` + `obs.tsv` + `var.tsv` + `uns.json`, and a Python rule
assembles the `.h5ad`. TSVs are written with `digits = 17` and the round-trip is
asserted on read.

## Conventions

- **Rules:** one `.smk` per phase in `workflow/rules/`, named `NN_phase.smk`.
  Rule names are prefixed with their task ID where one exists. Every rule
  declares `log:`, `benchmark:`, `conda:` and `threads:` — no exceptions.
- **Targets:** a phase contributes to `rule all` only when its switch in
  `config/config.yaml → phases` is on **and** its `TARGETS_*` list in the `.smk`
  is populated. Turn the switch on in the same commit that adds the rules.
- **Logs:** `results/logs/<rule>.log`. **Benchmarks:** `results/benchmarks/<rule>.tsv`.
- **Seeds:** every stochastic step (UMAP, permutation tests) takes
  `config["seed"]` and passes it explicitly. No implicit RNG anywhere.
- **Figures:** wrapped in `report(..., category=...)` with a caption, so
  `snakemake --report` is populated for free.
- **Config:** validated against `workflow/schemas/config.schema.yaml` at load;
  `samples.tsv` against `workflow/schemas/samples.schema.yaml` when it exists.
- **Commits:** `P<phase>-T<task>: imperative summary`. One phase per branch,
  squash-merged with the gate result in the message — because gates are
  phase-level, so a per-task branch would have no gate result to record
  (ADR 0004). Infrastructure tasks may branch and merge early; nothing that
  produces a number, table, figure or threshold may.
- **ADRs:** real decisions go in `docs/decisions/NNNN-slug.md`. A threshold
  without an ADR is a number someone made up.

## Installed skills

Pinned in `skills-lock.json`; restore with `npx skills experimental_install`.
The marimo set is deliberately a subset — the seven skills in
`marimo-team/skills` with no task behind them in this plan (paper-demo,
streamlit conversion, anywidget, molab badge, batch scheduling) are not
installed, on the same context-cost and supply-chain reasoning §5.1.2 applies
to the K-Dense library.

| Skill | Source | Use |
|---|---|---|
| `marimo-notebook` | `marimo-team/skills` | Correct notebook structure, `@app.cell` patterns, reactivity, `mo.ui` |
| `jupyter-to-marimo` | `marimo-team/skills` | The §4.5 conversion route for incoming Jupyter material |
| `wasm-compatibility` | `marimo-team/skills` | Checks an app-tier notebook before `marimo export html-wasm` |
| `marimo-pair` | `marimo-team/marimo-pair` | Drives a live kernel: run cells, inspect real state, iterate on plots |
| `retro-marimo-pair` | `marimo-team/marimo-pair` | Session retrospective on pairing friction |

**`marimo-pair` executes code in the user's live kernel.** That does not create
an exemption from hard constraint 5: anything discovered in a pairing session is
still not a result until a Snakemake rule produces it. Use the scratchpad to
*find* the answer, then write the script and the rule. Requires `bash`, `curl`
and `jq` on PATH (all present), and a notebook running under
`marimo edit --watch`.

## Stop and ask — these are scientific decisions, not implementation details

- Changing a QC threshold.
- Adding or removing a gene from a signature or the checkpoint panel.
- Dropping an AOI (and note: the policy is **flag, don't silently drop** —
  exclusions are listed in `results/tables/qc_excluded.tsv` with a reason).
- Changing a statistical model, its random-effects structure, or the
  multiplicity correction.
- Anything that would change the answer to an open question in §3 of the plan.

## Open questions (PROJECT_PLAN §3) — status as of 2026-08-24

Full answers, with citations, in `docs/data-provenance.md`. **Provisional** means
the evidence is metadata, not the data — P0-T4 confirms it against the file
before anything is hardened against it.

| # | Question | Status |
|---|---|---|
| Q1 | Is the GEO matrix raw counts, Q3-normalised, or already log-transformed? | **Resolved — Q3-normalised.** Column Q3 CV 0.059% vs library-size CV 13.4%. No zeros (min 2.12), so detection is background-relative |
| Q2 | Antibody-segmented compartments or geometric ROIs? | **Resolved — segmented.** PanCK+/PanCK− UV-cleavage within marker-guided ROIs |
| Q3 | Is slide / TMA / batch identifiable per AOI? | **Resolved — yes, from DCC filenames.** Two DSP runs (91/29), in `samples.tsv`. **`mLN`, `TBME`, `BC` are each wholly within one run** — batch inseparable from biology there |
| Q4 | Does per-AOI nuclei count / surface area survive into GEO metadata? | **Resolved — no.** `cell type` is GEO's only characteristics field |
| Q5 | Is patient-level clinical/survival metadata extractable and joinable? | **Resolved — yes.** Two time-to-event columns in Supplementary Data 1; Phase 5 is viable |

**Q2's caveat is the live constraint, and it outlives the question.** PanCK was
the *only* collection mask. CD45 and GFAP guided where a pathologist placed the
ROI; nothing was collected on a CD45 or GFAP mask. A `TIME` AOI is the
PanCK-negative segment of an ROI sited in a CD45-rich region — **not** a
CD45-sorted population; same for `TBME` and GFAP. Never write "CD45+ AOI", and
never treat a compartment label as a cell-type label. This is the same overclaim
hard constraint 6 forbids for "colocalisation", and it is why P0-T6's marker
sanity check is evidence rather than ceremony.

## Useful commands

```bash
snakemake -n --use-conda          # dry run (must stay clean)
snakemake --lint                  # must stay clean
# ALWAYS pass --use-conda. Without it the software-env rerun trigger fires,
# Snakemake tries to re-run P0-T2's downloads, and their protected() outputs
# raise ProtectedOutputException. A bare `snakemake -n` is not a clean dry run.
snakemake --report results/reports/workflow-report.html
uvx marimo check notebooks/**/*.py
marimo edit --watch notebooks/review/qc_review.py   # live pairing
```
