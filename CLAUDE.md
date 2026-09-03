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

**GATE 2 PASSED (2026-08-28). Phases 0, 1 and 2 are complete (P0-T1 … P0-T8,
P1-T1 … P1-T6, P2-T1 … P2-T7); current phase is Phase 3 (A4 — the checkpoint
landscape by compartment, exploratory per ADR 0008).** All five open questions
are answered (`docs/data-provenance.md`), and `docs/limitations.md` records what
the design cannot support.

Gate 2 evidence (ADR 0014): the published direction is reproduced. Brain minus
lung, within `TIME`, `TIME-B` n = 8 — antigen presentation **−0.477 SD
[−1.030, +0.075]** (q = 0.098, and the only set detected 13/13 at *both* sites),
cytotoxicity **−0.946 SD [−1.590, −0.302]** (q = 0.018). Antigen presentation
does not clear FDR 0.05 and is not expected to: −0.48 SD sits below the
1.1–1.3 SD power floor, so that is uninformative, not negative. Direction agrees
with the paired check in 12 of 12 cells; lme4 and statsmodels agree to 5e-9.

**The dominant limitation is now detection, ahead of `TIME-B` n = 8.** Three of
six Phase 2 signatures failed their coverage floor — `exhaustion` 1/6 and `tls`
1/5 in brain, `myeloid_m1` 2/10 at *both* sites — and two of those show large,
nominally significant "reductions in brain" that are the ADR 0008 artefact and
**are not reportable**. Only antigen presentation, cytotoxicity and myeloid M2
are interpretable. Secreted ligands and chemokines are systematically
undetected.

Gate 1 evidence: **41% of variance is compartment-attributable** — PVCA over 14
PCs spanning 60.4% of variance, per-gene median 21% across the 18,691 of 18,694
genes that converged, permutation null 0.3%
(`results/tables/variance_partition_summary.json`, P1-T3). Confirmed
independently by P1-T4's clustering: compartment ARI **0.741**, patient
**0.003**. The model specification this licenses is **ADR 0009**.

Gate 0 evidence: marker sanity check 4/4
(`results/tables/marker_sanity_verdict.tsv`), design table matches §2.1 exactly
(120 AOIs), Q2 resolved to antibody segmentation.

**Scope change at Gate 2: Aim A5 is demoted to exploratory (ADR 0014)** — a
ligand–receptor analysis needs the ligand above background, and Phase 2 measured
the ligand side directly at both sites and found it largely absent. Phase 4 runs;
no A5 result may be a headline claim, and the anticipated null is an
assay-sensitivity limit, never evidence that the crosstalk is absent.

**Scope change at Gate 0: Aim A4 is demoted to exploratory (ADR 0008).** Most
of the checkpoint panel sits at background — CTLA4 is above background in 2 of
8 `TIME-B` AOIs, TIGIT in 1, IDO1 in 2 — and brain background is *higher*, so a
near-background gene reads as depleted in brain artefactually. Phase 3 still
runs; **no A4 result may be a headline claim**, every checkpoint reported states
its per-site detection count, and a gene detected in fewer than half the
`TIME-B` AOIs is "not assessable in brain", never "lower in brain".
A5 was flagged here and **re-decided at Gate 2 — demoted (ADR 0014)**.

Phase 0 output: `results/interim/aoi_normalised.h5ad` — 120 AOIs × 18,694
genes, `X` = log2(Q3 + 1), `layers['q3']` the untransformed values, 36 `obs`
columns, git SHA and config hash in `uns`.

**The number to keep in view:** with `TIME-B` n = 8, the lung-vs-brain immune
contrast detects roughly **1.1–1.3 SD** at 80% power (P0-T8). A null result in
Phase 2 is uninformative, not negative.

**A2 is answered: compartment drives the variance, not patient and not site.**
Phase 2 therefore estimates the lung-vs-brain contrast **within** compartment,
never pooled across compartments, and every model carries `(1|patient)` for the
reason ADR 0009 gives — non-independence, not variance share.

Phase 2 output: `signature_scores.tsv`, `signature_coverage.tsv`,
`signature_models.tsv`, `paired_concordance.tsv`, `decon_composition.tsv`,
`convergence_check.tsv`; figures `contexture_heatmap.png`,
`decon_composition.png`; prose in `docs/analysis-notes.md`.

**P3-T2b added an external check (ADR 0017), and it is the only one the project
has.** Everything else is internal: P3-T2's three "independent implementations"
all read the same `obs['negprobe']` and test the same hypothesis. Against the
source paper's own deposited Source Data, **all 2,243,280 values of
`layers['q3']` are identical** and the published `NegProbe-WTX` row equals
`obs['negprobe']` exactly — so this project's inputs are demonstrably the
published inputs. Two consequences worth carrying:

- **The 119-vs-120 gap is closed: the paper dropped `TBME15b`** (agreement
  1.0000 vs `TBME15a`, 0.0003 vs `TBME15b`). No AOI failed QC; P15 contributes
  two `TBME` AOIs and the paper used one. All 20 are analysed here.
- **The count scale is measured**, validated against the `.dcc` raw counts for
  120/120 AOIs. In `TIME-B`, background is ~30 counts and CTLA4/TIGIT/IDO1 carry
  37–44 — the same fact as the detection floor, stated as resolution. An AOI
  holds a median of 645 distinct values across 18,694 genes.

**P3-T3 (ADR 0018): the primary restriction stands, a declared exploratory
secondary carries the rest.** A gene enters a primary fit only if it clears the
detection floor in *both* groups — 2 of 9 genes per site. `CD274` (PD-L1) is
enriched in the lung immune compartment, **+0.612 SD [+0.361, +0.863], q =
0.0001**, surviving both sensitivities. Everything else is a separate,
unadjusted, exploratory-within-exploratory table on which no claim may rest.

**Detection and expression move OPPOSITE ways under the same background
gradient, and conflating them inverts the argument.** Detection is
`q3 > 2 × negprobe`, so higher background → detected less → looks *depleted*.
Expression is `log2(q3 + 1)`, where background adds to signal, so higher
background → looks *enriched*. The models are on expression, background is
higher in the immune compartments, so the gradient is **permissive** there —
adjusting for `negprobe_log2` shifts the excluded genes by a mean of −0.155 and
the admitted ones by +0.010. The genes the restriction excludes are the genes
the adjustment moves.

**P3-T4: the lung-vs-brain checkpoint shift is null in every assessable gene,
and uninformative rather than negative.** Primary is 4 genes in the immune
compartment and 2 in tumour; largest is CD274 at −0.260 SD [−0.595, +0.075],
q = 0.462, brain minus lung, `TIME-B` n = 8. All far below the 1.1–1.3 SD power
floor. Every P3-T4 gradient is `negligible` (site background deltas +0.03 and
−0.06) — **the background artefact is a P3-T3 problem, not a P3-T4 problem**,
now measured rather than assumed. The paired check is direction only in BOTH
compartments and agrees 9 of 11; the 23-patient tumour set deliberately gets no
p-value either.

**P3-T5 is the deliverable figure** (`checkpoint_dotplot.png`): genes ×
compartment × site, all 7 compartments, all 120 AOIs. Dot area = detection rate,
colour = `median_negprobe_ratio` with its midpoint read from
`qc.detection_background_multiple` so it cannot drift from ADR 0007's rule,
hatch = below the pre-registered floor. **40 of 63 cells are hatched and that is
the finding** — the panel clears the floor in `TIME-L` alone.

**P3-T6 is the same audit made interactive** (`checkpoint_explorer.py`, exported
to WASM at `results/reports/checkpoint_explorer`). **Its detection-floor slider
is a sensitivity display, not a threshold control** (ADR 0019): the floor of
record stays 0.5, the app defaults there and reproduces P3-T5's figure there,
it labels itself the moment it leaves, and nothing it reaches is a result. The
2× background multiple is deliberately not exposed — `detection_rate` is already
computed at it. Model estimates sit in their own section **below** the plot and
never on it: the plot is detection, they are expression, and the two move
opposite ways under the same gradient (ADR 0018).

Next: Phase 3 P3-T7. A4's honest output is reconnaissance, not a claim
(ADR 0008).

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
  `TIME-L`; P15 has two `TBME`; 33 of 42 subjects contribute more than one AOI.
  A random intercept for patient is **mandatory**, not stylistic — but justify
  it on **non-independence, never on variance share** (ADR 0009). Patient is
  only 12% of variance per-gene and its clustering ARI is 0.003; the evidence
  that matters is the +0.107 ρ contrast across the 102 same-patient AOI pairs
  whose *compartments differ*. A reader who sees only the ARI will conclude the
  random intercept can be dropped. It cannot.
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
conda activate nsclc_bm_spatial   # ALWAYS first — see below
# ALWAYS pass --conda-prefix "$HOME/nsclc-envs" too. That symlink points AT
# .snakemake/conda, so the envs do not move — but conda embeds the path string
# it is given, and this working directory contains a space ("Biomedical Data
# Science"). Three separate conda/bioconda scripts interpolate that prefix
# unquoted: conda-forge's R wrapper, Snakemake's post-deploy hook, and
# bioconda's installBiocDataPackage.sh. The last one makes r-geomx impossible
# to CREATE. Recreate the symlink with:
#     ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"
# Do NOT relocate the envs to a genuinely different directory: that fires the
# software-env rerun trigger on p0t2_fetch_geo and its protected() outputs
# abort the DAG. See ADR 0013.
snakemake -n --use-conda --conda-prefix "$HOME/nsclc-envs"   # must stay clean
snakemake --lint                  # must stay clean
# ALWAYS pass --use-conda. Without it the software-env rerun trigger fires,
# Snakemake tries to re-run P0-T2's downloads, and their protected() outputs
# raise ProtectedOutputException. A bare `snakemake -n` is not a clean dry run.
#
# ALWAYS run from the nsclc_bm_spatial env (Snakemake 8.30). Base anaconda has
# 9.20, which writes .snakemake/metadata in a format 8.30 reads as stale — it
# then provenance-triggers P0-T2 into the same protected-output abort.
snakemake --report results/reports/workflow-report.html
uvx marimo check notebooks/**/*.py
marimo edit --watch notebooks/review/qc_review.py   # live pairing
```
