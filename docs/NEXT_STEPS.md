# Next steps — Phase 3 (checkpoint landscape by compartment)

**Last session:** 2026-08-28
**Branch:** `main` at `af12d48` (Phase 2 squash-merge), **pushed**. Phase 3 needs
a new `P3` branch.
**Gate 0:** PASSED. **Gate 1:** PASSED. **Gate 2:** PASSED — ADR 0014 is the record.
**Phase 2:** **COMPLETE** (P2-T1 … P2-T7). Nothing outstanding.

Phase 3's goal (§6): Aim **A4** — which compartment carries each checkpoint gene,
and whether that shifts brain vs. lung. §6 calls it "the most differentiable
analysis in the project"; ADR 0008 demoted it to **exploratory** at Gate 0 and
nothing since has changed that.

---

## Read this before writing any code

**Phase 2 turned the detection ceiling from a checkpoint-panel worry into a
measured, assay-wide property, and it is now the project's dominant limitation —
ahead of `TIME-B` n = 8.**

Three of six Phase 2 signatures failed their coverage floor:

| set | `TIME-L` | `TIME-B` | |
|---|---|---|---|
| `antigen_presentation` | 13/13 | 13/13 | the only fully-detected set |
| `myeloid_m2` | 6/9 | 6/9 | assessable |
| `cytotoxicity` | 2/4 | 2/4 | assessable, **at the floor** |
| `tls` | 3/5 | **1/5** | not assessable in brain |
| `exhaustion` | 4/6 | **1/6** | not assessable in brain |
| `myeloid_m1` | **2/10** | **2/10** | not assessable **at either site** |

`exhaustion` and `tls` returned large, nominally significant "reductions in
brain" (q = 0.005 and 0.018) that are **not reportable** — they are the ADR 0008
artefact, exactly as predicted. **P3-T1's panel overlaps `exhaustion` almost
entirely** (`PDCD1`, `CTLA4`, `LAG3`, `HAVCR2`, `TIGIT`), so expect the same
outcome and design P3-T2 to produce it cleanly rather than to be surprised by it.

**Secreted ligands and chemokines are systematically undetected** — CXCL10,
CXCL11, IL1B, TNF, IL12B, IL10, CCL22, CCL19, CCL21 all sit below background.
That is what demoted A5 (ADR 0014).

**Gate 3 already anticipates the null:** "If everything is not-assessable, that's
still a legitimate finding." Phase 2 makes that the *likely* outcome. Write
P3-T2 as the deliverable it is, not as a gate to get past.

---

## What Phase 2 decided that binds Phase 3

1. **`(1|patient_id)` in every model**, on non-independence grounds (ADR 0009).
   The config schema now enforces it with a regex on the formula fields — a
   config edit that drops it fails at Snakefile load. Changing the
   random-effects structure remains a stop-and-ask.
2. **Primary + batch sensitivity, both reported** (ADR 0012). `dsp_run` moved
   estimates by a median of 0.0076 (max 0.0465) and changed no direction or
   ordering. Half the fits (12/24) were singular and are surfaced, not hidden.
3. **The coverage gate is not optional.** It is what stopped `myeloid_m1`'s
   apparent signature/deconvolution contradiction being written up as a real
   methodological conflict. P3-T2 is the same mechanism; give it the same weight.
4. **Never "lower in brain" for a set below its floor** — "not assessable in
   brain" (ADR 0008). The P2-T7 heatmap hatches such cells *on the figure*,
   because a heatmap is exactly what someone reads "lower" off a colour from.
   P3-T5 must do the equivalent: §6 says not-assessable genes are "rendered
   distinctly, not omitted".

---

## Before you write anything — the one `config.yaml` edit

**Phase 2 made this more expensive than it was.** `config/config.yaml` is a
declared input of `p0t7_assemble_h5ad` — its SHA-256 goes into `uns`, which is
correct provenance — so *any* edit, a comment included, rebuilds the `.h5ad` and
now re-runs **Phase 1 and Phase 2**: P1-T3's ~16-minute variance fit plus the
whole contexture chain. **Batch every config change into a single commit.**

Everything Phase 3 needs from `config.yaml`, in one edit:

```yaml
phases:
  checkpoints: true          # same commit that populates TARGETS_CHECKPOINTS

checkpoints:
  panel_yaml: config/checkpoints.yaml
  # PRE-REGISTERED. Chosen and committed BEFORE looking at any result
  # (§6 P3-T2, and ADR 0008 for why). Choosing it after is the failure this
  # task exists to prevent.
  detected_in_aoi_fraction: <choose>   # gene detected in >= this share of a group's AOIs
  min_assessable_fraction: <choose>    # below this, the gene is "not assessable" there
  compartments: [...]                  # which aoi_codes enter the models
  model:
    carrier: "expression ~ compartment + (1|patient_id)"  # P3-T3, within site
    shift:   "expression ~ site + (1|patient_id)"         # P3-T4, within compartment
    sensitivity_covariate: dsp_run                        # ADR 0012's shape
    df_method: satterthwaite
    fdr_method: fdr_bh
    fdr_alpha: 0.05
```

Mirror the `contexture` block's schema style in
`workflow/schemas/config.schema.yaml` — `additionalProperties: false`, and
**reuse the `(1|patient_id)` regex guard on both formula fields**. That guard is
what makes a config edit unable to silently drop the random intercept, and it
was negative-tested when it went in.

**The panel itself goes in `config/checkpoints.yaml`, which is deliberately NOT
an input of `p0t7_assemble_h5ad`** — same reasoning as `signatures.yaml`.
Editing the panel costs seconds; editing `config.yaml` costs a Phase 1 + 2
rerun. Keep gene membership out of `config.yaml`.

`workflow/rules/05_checkpoints.smk` is still a stub with
`TARGETS_CHECKPOINTS = []`. Populate the list and flip `phases.checkpoints` in
the same commit — a phase contributes to `rule all` only when both are true.

**Next free ADR number: 0015.**

---

## What to reuse — do not write any of this twice

Phase 3 is structurally Phase 2 again with a different panel and one extra
question. Almost everything has a working template:

| need | reuse | note |
|---|---|---|
| detection per gene × group | `workflow/scripts/score_signatures.py` | the `q3 > multiple × negprobe` block. **Do not invent a second detection rule** |
| load + validate a panel YAML | the `SIGNATURES` block in `workflow/rules/common.smk` | validates at Snakefile load, so a malformed panel fails at parse |
| panel schema | `workflow/schemas/signatures.schema.yaml` | including the `PMID:` / `MSigDB:` source pattern that rejects free text |
| resolve genes against the matrix | `workflow/scripts/resolve_signatures.py` | fails loudly on a symbol that stops resolving |
| mixed models | `workflow/scripts/fit_signature_models.R` | primary + sensitivity, Satterthwaite df, BH within family, `isSingular()`, and **both** the warning and message handlers |
| implementation cross-check | `workflow/scripts/crosscheck_mixedlm.py` | compares estimate and SE only; df/p are not comparable |
| paired direction check | `workflow/scripts/paired_check.py` | asserts no test statistic reaches the table |
| rendering not-assessable | `workflow/scripts/contexture_heatmap.py` | hatches those cells *on the figure*; §6 asks P3-T5 for the same |
| fetching a pinned artifact | `p2t5_fetch_reference` + `p2t5_record_reference_provenance` | reuses `acquire_geo.py` / `record_provenance.py` unchanged — but note the wildcard **must** be named `artifact`, since the script reads it by name |
| app-tier WASM export | `p1t5_landscape_explorer` in `03_landscape.smk` | stages the notebook beside a `public/` dir; **read ADR 0010 first** |

---

## Start here

### 1. P3-T1 — define the panel (~1.5 h)

`config/checkpoints.yaml` is still empty, deliberately. Populate it with
`CD274`, `PDCD1`, `CTLA4`, `LAG3`, `HAVCR2`, `TIGIT`, `IDO1`, `VSIR`, `CD276`,
one line of sourced clinical rationale each. **Adding or removing a gene is a
stop-and-ask.**

**Verify every citation against the record.** Two of the four PMIDs proposed for
`signatures.yaml` in Phase 2 pointed at unrelated papers — right journal, right
year, adjacent subject — and were caught only by querying PubMed. A schema can
enforce that a source is *present*, never that it is *correct* (ADR 0011).

Write `workflow/schemas/checkpoints.schema.yaml` and validate it at Snakefile
load, the way `signatures.yaml` is validated in `common.smk`.

### 2. P3-T2 — detection audit FIRST (~2 h)

**Set the detection floor in `config.yaml` before looking at any result.** §6
says pre-registered; ADR 0008 says the same thing for a different reason.

Reuse the existing rule rather than inventing a second one: a gene is detected in
an AOI when `layers['q3'] > qc.detection_background_multiple × obs['negprobe']`.
`score_signatures.py` already implements exactly this — lift it rather than
rewrite it.

Report per gene × compartment, **never pooled**. Output the panel split into
assessable / not assessable.

### 3. P3-T3 — which compartment carries each checkpoint? (~3 h)

`expression ~ compartment + (1|patient_id)`, **within site**. This is the
analysis bulk RNA-seq structurally cannot do, and §6 calls it the project's core
claim. `fit_signature_models.R` is the template — primary/sensitivity pair,
Satterthwaite df, BH-FDR within family, singular-fit reporting, and the
`site`-relevel assertion (here it becomes a *compartment* reference level, so
assert whichever one you pick).

**Phase 2 never had to face this and Phase 3 does: one of the compartment
contrasts is batch-confounded and cannot be fixed by modelling.** Run split by
compartment, from `results/tables/batch_crosstab.tsv` (P0-T3, Q3):

| | run A | run B | |
|---|---|---|---|
| `L` | 22 | 8 | spans both |
| `LB` | 21 | 6 | spans both |
| `TIME-L` | 9 | 6 | spans both |
| `TIME-B` | 6 | 2 | spans both, barely |
| **`TBME`** | **20** | **0** | **run A only** |
| `mLN` | 13 | 0 | run A only |
| `BC` | 0 | 7 | run B only |

So **within lung**, `L` vs `TIME-L` is clean — both span runs. **Within brain**,
`LB` vs `TIME-B` is clean, but **any contrast involving `TBME` confounds batch
with biology and no covariate recovers it** (ADR 0009 §3). `TBME` is also the
lowest-detecting compartment at 21.4% and carries 4 of the project's 5
QC-flagged AOIs. Every glial-compartment result must say this inline; it is the
same reasoning that demoted A5.

### 4. P3-T4 — does the profile shift lung → brain? (~3 h)

Within compartment, never pooled — §6 is explicit that a tumour-compartment
`CD274` shift and an immune-compartment `CD274` shift mean different things
therapeutically. Unpaired primary + the 5-patient paired sensitivity (direction
only, no p-values — `paired_check.py` asserts none reach the table). **Every
brain claim states `TIME-B` n = 8 inline** (hard constraint 8).

Note the two tasks answer different questions with the same fits, so decide the
FDR family deliberately: T3 corrects across genes within a site, T4 across genes
within a compartment. Correcting across both at once would be a third family
nobody asked for.

### 5. P3-T5 — deliverable dot plot (~2.5 h)

Genes × (compartment × site); dot size = detection rate, colour = mean
expression. **Not-assessable genes are rendered distinctly, not omitted** — §6
says so, and absence of evidence has to be visible. `contexture_heatmap.py`
hatches exactly this case on the figure and its caption explains why; copy the
approach.

The caption carries `TIME-B` n = 8 and the power floor, states the detection
floor and where it was pre-registered, and — if any contrast involves `TBME` —
says that batch is inseparable from biology there.

### 6. P3-T6 — `checkpoint_explorer.py`, the flagship artifact (~4 h)

A marimo app over the T3/T4 outputs: a slider for the detection floor, a
compartment multiselect, a site toggle, driving the dot plot live. §6's point is
that a reader can *watch* genes fall out of "assessable" as the floor rises —
the most honest presentation available of near-LOQ WTA data, and the one place
where this project's central limitation becomes something you can feel rather
than read.

**Read ADR 0010 first.** It ships data in `public/` and fetches over HTTP;
`cache_cells` breaks `mo.ui.table`, `mo.persistent_cache` silently ships an
export with no data, and pandas cannot read an `http://` URL under Pyodide. The
PEP 723 block installs *unpinned* marimo, so APIs skew between check time and
run time. **`marimo check` passes while the export is broken — verify by
exporting and opening it, never by linting.** `p1t5_landscape_explorer` is the
working rule to copy.

Accept criterion worth honouring literally: **the app's defaults must exactly
reproduce the static P3-T5 figure.** If they don't, one of the two is wrong, and
finding out which is the point of the criterion.

### 7. P3-T7 — clinical interpretation (~2 h)

~600 words appended to `docs/analysis-notes.md` (the file now exists; P2-T6's
section is the format). Hedge proportionally to n.

---

## Stop and ask — Phase 3's are not implementation details

`CLAUDE.md` lists these generally; these are the ones Phase 3 will actually hit.

1. **The panel membership** (P3-T1). Adding or removing a checkpoint gene is a
   scientific decision, including while first populating the empty file.
2. **The detection floor** (P3-T2). It is a QC threshold, it must be
   pre-registered, and it decides which genes are reportable at all — which at
   this detection level is most of the result. Needs an ADR citing where the
   value came from; a threshold without an ADR is a number someone made up.
3. **Any change to the random-effects structure.** `(1|patient_id)` is
   mandatory on non-independence grounds (ADR 0009) and the config schema
   enforces it. Singular fits are *not* a reason to drop it — half of Phase 2's
   fits were singular and it stayed.
4. **Dropping an AOI.** The policy is flag, don't drop; exclusions go in
   `results/tables/qc_excluded.tsv` with a reason.

## Gate 3, and why the likely outcome is fine

> **GATE 3** — At least one checkpoint gene shows an interpretable
> compartment-resolved pattern with a CI you're willing to show a stranger. If
> everything is not-assessable, that's still a legitimate finding.

Phase 2 makes "everything not-assessable" a genuinely plausible result, not a
worst case: the panel overlaps the `exhaustion` signature almost entirely, and
`exhaustion` cleared 1 of 6 genes in brain. **Write P3-T2 as the deliverable it
is.** A clean, well-quantified statement of what a WTA assay cannot see in
8 brain immune AOIs is a better output than a fragile positive — and per ADR
0008 a fragile positive here could not have been reported as biology anyway.

What would make it a weak phase is not a null; it is a null reported without
the detection table that licenses it.

---

## State you'll have forgotten

- **`--conda-prefix "$HOME/nsclc-envs"` is MANDATORY on this machine**, alongside
  `--use-conda`. That symlink points *at* `.snakemake/conda`, so nothing moves —
  but this working directory contains a space, and **three** separate scripts
  interpolate the conda prefix unquoted: conda-forge's R wrapper, Snakemake's
  post-deploy hook, and bioconda's `installBiocDataPackage.sh`. The last makes
  `r-geomx` impossible to *create*. Recreate the symlink with
  `ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"`. **Do not** relocate the
  envs to a genuinely different path — that fires the software-env trigger on
  `p0t2_fetch_geo` and its `protected()` outputs abort the DAG. ADR 0013.
- **`rule p2t0_repair_r_env` must run before any R rule.** It quotes conda-R's
  relocated path assignments; declare `results/interim/env_repair/<env>.ok` as an
  input of every R rule. The wildcard already accepts `r-stats|r-geomx`.
- **Both R envs are now built and verified.** `r-stats` (lme4 2.0.6 / Matrix
  1.7.5) and `r-geomx` (SpatialDecon 1.16.0, GeomxTools, standR). `r-geomx.yaml`
  carries `r-lme4>=2.0` + `r-reformulas` pins without which SpatialDecon solves
  green and then fails to load.
- **A green `--conda-create-envs-only` is not evidence an env works.** It
  happened twice in one session. Load the libraries.
- **ADR 0001's exact float round-trip does not hold Python → R.** R's parser is
  not correctly rounded (`R_strtod` reads Python's `0.44525532065683371` one ULP
  low, on 34 of 276 values). Exactness is kept R → Python, where it does hold;
  Python → R asserts structure, and numeric agreement is checked end-to-end by
  `p2t3_model_crosscheck` instead. ADR 0013.
- **`config/config.yaml` is a declared input of `p0t7_assemble_h5ad`** — its
  SHA-256 goes into `uns`, so *any* edit, comment included, correctly rebuilds
  the `.h5ad` and re-runs Phase 1 (~16 min). Batch config changes. Signature-like
  content belongs in its own file for this reason; `config/signatures.yaml` and
  `config/checkpoints.yaml` are not inputs of that rule.
- **`py-analysis.conda-lock.yml` is stale** and `conda-lock` is not installed.
  Still the P6-T1 blocker. Phase 2 added no new Python dependency.
- **P6-T1 must be run on a path containing a space**, or it will not exercise
  `p2t0_repair_r_env` or the symlink and will report a false pass.
- **`Markdowns/PROJECT_PLAN.md` is gitignored** — ADRs are the durable record.
  ADR 0014 §4 lists two defects found in §6 at the Gate 2 audit: P2-T7's
  "publication-grade" criterion is unfalsifiable, and the A5 re-decision was
  owned by no task row.

## Useful commands

```bash
conda activate nsclc_bm_spatial
ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"   # once per machine
git checkout -b P3
snakemake -n --use-conda --conda-prefix "$HOME/nsclc-envs"   # must stay clean
snakemake --lint                                             # must stay clean
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" --cores 4
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" \
  --report results/reports/workflow-report.html
/gate 3
```

---

## Phase 2, for reference

Gate 2 (ADR 0014): the published direction is reproduced. Brain minus lung,
within `TIME`, `TIME-B` n = 8:

| signature | coverage | estimate (z-score) | q | |
|---|---|---|---|---|
| antigen presentation | 13/13 both | −0.477 SD [−1.030, +0.075] | 0.098 | direction recovered, underpowered |
| cytotoxicity | 2/4 both | −0.946 SD [−1.590, −0.302] | **0.018** | recovered, significant |
| myeloid M2 | 6/9 both | +0.091 SD [−0.333, +0.514] | 0.660 | null |

−0.48 SD sits below the 1.1–1.3 SD power floor, so antigen presentation's
non-significance is uninformative, not negative. Paired direction agreed 12/12;
lme4 vs statsmodels 5e-9; the two scoring methods correlate at ρ 0.79–0.97.

P2-T5 deconvolution ran two arms — `safeTME` on both sites (comparable) and
§6's two-matrix approach (**not** comparable: `Lung_HCA` resolves 9 lymphoid
types, `Brain_Darmanis` resolves **0**). P2-T6 found exactly **one** signature
genuinely testable for convergence — cytotoxicity — and it agrees: −0.946 SD
against CD8+NK proportion 0.139 → 0.081.

Artefacts: `signature_membership.tsv`, `signature_scores.tsv`,
`signature_coverage.tsv`, `signature_models.tsv`, `model_crosscheck.json`,
`paired_concordance.tsv`, `decon_composition.tsv`, `convergence_check.tsv`,
`results/figures/contexture_heatmap.png`, `results/figures/decon_composition.png`,
`docs/analysis-notes.md`. ADRs **0011** (signature membership and sourcing),
**0012** (batch term), **0013** (conda R on a path with spaces; the float
boundary), **0014** (Gate 2 record; A5 demotion).
