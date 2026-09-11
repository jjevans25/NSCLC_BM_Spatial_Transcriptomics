# Next steps — Phase 4 (inferred cross-compartment crosstalk)

**Last session:** 2026-09-11
**Branch:** `main` at `bce2bde` — **Phase 3 is MERGED and PUSHED.**
`origin/main` is the same commit, working tree clean. **Phase 4 starts by
branching `P4` from `bce2bde`** (ADR 0004, one branch per phase).
**Gate 0:** PASSED. **Gate 1:** PASSED. **Gate 2:** PASSED (ADR 0014).
**Gate 3:** PASSED (ADR 0020). **Gate 4:** open.
**Phase 4:** nothing started. `workflow/rules/06_crosstalk.smk` is a stub with
`TARGETS_CROSSTALK = []`, and `config/config.yaml → phases.crosstalk` is
`false`. Both flip in the same commit as the first rule.
**Next free ADR number: 0021.** P4-T1's pairing rule is ADR 0021 and it must be
committed **before anything runs** — the plan says so explicitly.

**Phase 3 is fully closed.** The last open acceptance criterion — ADR 0010's
browser check on the P3-T6 WASM export — was **verified on 2026-09-11**: served
over HTTP from `results/reports/checkpoint_explorer`, the dot plot renders and
the four model tabs render as *tables* rather than `UnhashableStub` text. That
stub-text failure is the one ADR 0010 says only a browser can catch — `marimo
check --select MW` passed throughout while P1-T5's export was shipping no data
at all — so this was the real test, not a formality. **Nothing from Phase 3 is
outstanding.**

---

## START HERE — read this before writing a line of Phase 4

Phase 4's goal (§6) is Aim **A5**: nominate inferred ligand–receptor crosstalk
between adjacent compartments. **ADR 0014 demoted A5 to exploratory at Gate 2,
and it did so on evidence, not caution.** The terms are binding:

- **No A5 result may be a headline claim.**
- **Every nomination states the detection status of both partners**, per site.
- **The anticipated null is an assay-sensitivity limit, NEVER evidence that the
  crosstalk is absent.** This is the same distinction Phases 2 and 3 both turned
  on, and it is the one a reader gets wrong from silence.
- **Never write "colocalisation" or "spatially adjacent"** (hard constraint 6).
  The phrase is "inferred crosstalk between adjacent compartments." You have no
  coordinates. P4-T7 greps the whole repo for violations.

**The ground has shifted once since ADR 0014, in Phase 4's favour, and the
handoff would be misleading without it.** ADR 0014's case for demotion was that
"secreted ligands and chemokines are systematically undetected" — measured on
`myeloid_m1`, `myeloid_m2` and `tls`, which are *cytokine* signatures. That
finding holds and is if anything stronger than it looked. But it does not
generalise to the whole ligand side: membrane and matrix ligands are detected in
nearly every AOI. See **Reconnaissance** below. Phase 4 is therefore not vacuous
— it is *restricted*, and naming which class survives is itself part of the
deliverable.

---

## The tasks (PROJECT_PLAN §6, Phase 4 — est. 10–14 h)

### P4-T1 — Build paired sets, and write ADR 0021 FIRST (~1.5 h)

**Verified against `results/tables/design_matrix.tsv`, so use these, not the
plan's estimates — they agree, and now they are checked:**

| adjacency | patients | who |
|---|---|---|
| lung, `L` ↔ `TIME-L` | **13** | P5, P12, P15, P18, P19, P24, P26, P29, P30, P32, P35, P40, P43 |
| brain, `LB` ↔ `TIME-B` | **8** | P5, P12, P14, P15, P19, P20, P31, P35 |

**Only two patients have a duplicate AOI on either side, and both duplicates are
on the same side of the same adjacency:** P12 and P24 each contribute **two
`TIME-L`** AOIs and one `L`. There is no duplicate anywhere in the brain set.
`P15`'s duplicate is `TBME`, which Phase 4 does not pair.

So ADR 0021 decides exactly one thing: **average the two `TIME-L` AOIs, or take
the higher-QC one.** Recommendation — **average on the log2 scale, pre-registered
before looking at any rho** — because "higher-QC" needs a tiebreak rule that
does not exist yet, and inventing one after seeing the correlations is the
failure mode ADR 0016 §3 was written to prevent. Record the rejected option and
state that the choice moves 2 of 13 lung patients and 0 of 8 brain ones, so it
cannot be load-bearing. Then run the alternative as a **sensitivity** and report
both, the way ADR 0012 handles batch.

**Accept:** two pairing tables; ADR 0021 committed before the first rule runs.

### P4-T2 — Ingest the LR database, pinned (~2 h)

CellChatDB or CellPhoneDB. **Pin the version** and put it under
`resources/` — which is **read-only to rules** (hard constraint 1), so it is
acquired the way P0-T2 and P3-T2b acquired theirs: a sanctioned writer with a
provenance row and a checksum (`resources/supplementary_provenance.tsv` and
`resources/supplementary_checksums.sha256` are the pattern; ADR 0005 has the
rules and P3-T2b is the most recent worked example — it was the third sanctioned
writer, this is the fourth).

**The detection filter is a pre-registered threshold and therefore needs an ADR
line, not a number chosen at the keyboard.** Pre-set it *before* seeing the
ranked table: a pair enters only if **both** members clear detection in ≥ a
stated fraction of AOIs **in the compartment that member is measured in** —
ligand in the tumour AOIs, receptor in the immune AOIs — at the project's one
detection rule, `q3 > 2.0 × negprobe` (ADR 0007). Reuse the 0.5 floor Phase 3
pre-registered unless there is a stated reason not to; a *different* floor here
than in Phase 3 is a decision, not a detail.

**Do not add a dependency to `py-analysis.yaml`** — that is `p0t2_fetch_geo`'s
env and its outputs are `protected()`; the software-env trigger aborts the DAG.
P3-T2b parsed `.xlsx` with `zipfile` + `ElementTree` rather than add `openpyxl`.
A CellPhoneDB/CellChatDB release is CSV or TSV, so this should not arise — but
check before reaching for a package.

**Accept:** filtered LR table with counts before *and* after filtering. The
before/after gap is a reportable number, not bookkeeping — it is the quantitative
form of ADR 0014's claim.

### P4-T3 — Correlate across patients (~3 h)

Ligand in the PanCK⁺ tumour AOI vs. receptor in the paired immune AOI, across
patients, **within site**. **Spearman, not Pearson** — at n = 13 and n = 8
Pearson is a coin flip on one outlier. Rank by |rho|.

**Both directions of each pair must be tested and labelled**, since "ligand on
tumour, receptor on immune" and the reverse are different biological claims and
CellChatDB is directed. Report n per pair explicitly: a pair can be filtered in
overall and still lose patients to per-AOI detection.

**Accept:** ranked table with rho, raw p, and n per pair.

### P4-T4 — Null calibration — THE task the gate turns on (~2.5 h)

With ~1000 pairs at n = 13, |rho| > 0.7 arises by chance. **Permute patient
labels** — the pairing, not the expression — to get an empirical null, and report
an **empirical FDR** column. Without this the ranking is uninterpretable, and
Gate 4 is written against exactly this number.

Seed from `config["seed"]`, passed explicitly (no implicit RNG anywhere). Permute
*within* site so the null preserves the two n's. `p1t3_variance_partition`'s
permutation null is the in-repo precedent for how this project structures one.

**Accept:** null distribution figure + empirical FDR column. **Skill:**
`statistical-analysis`.

### P4-T5 — Nominate 5–10 pairs (~1.5 h)

Top-ranked, detection-filtered, biologically coherent, **framed explicitly as
hypothesis generation**. Each row carries rho, empirical FDR, n, the per-site
detection of both partners (ADR 0014), and one line of biological rationale.

**If nothing survives the empirical FDR, that is the deliverable** — write it as
such and do not go looking for a softer threshold. See Gate 4.

### P4-T6 — External cross-reference (~2 h)

Against the Nat Commun 2024 NSCLC Visium result: `NRP1–VEGFA`,
`NECTIN2–TIGIT`, `LGALS9–HAVCR2` reported adjacent; `PD1–PDL1` not. Informal
concordance, **not validation**. **Skills:** `paper-lookup`, `literature-review`.

**Two of the four comparators are already compromised, and you should know which
before you write the task rather than after:**

- **`LGALS9` is not on the panel at all.** The 18,694 measured genes include
  `LGALS1`, `LGALS3`, `LGALS8`, `LGALS9C` and others — **`LGALS9` itself is
  absent.** So `LGALS9–HAVCR2` cannot be checked, and **substituting `LGALS9C`
  is a stop-and-ask** (adding a gene to an analysis panel — see below). Report it
  as not measurable; do not quietly swap the paralog.
- **`TIGIT` is detected in 0 of 30 `L` AOIs**, so `NECTIN2–TIGIT` cannot be
  evaluated in the direction that matters in lung. `NRP1–VEGFA` is the one
  comparator where both partners are detected nearly everywhere and the
  concordance check is genuinely informative.

The paragraph on what agreement and disagreement each mean is the actual
deliverable here. Disagreement with a Visium study is **not** evidence against
either result — different assay, different resolution, no coordinates on this
side. Say so.

### P4-T7 — Network figure + language audit (~2 h)

LR network, nodes coloured by compartment, edges by rho. **Skill:** `networkx`.
Then:

```bash
grep -rniE "colocali[sz]|spatial(ly)? (proximit|adjacen|close)" \
  --exclude-dir=.git --exclude-dir=.snakemake .
```

Fix every hit. **Accept:** figure + a clean grep.

Figure conventions that are already settled and must not be re-litigated: wrap
in `report(..., category=...)` with a caption; **`RdBu_r` means a signed effect**
(P2-T7) and **`PRGn` centred on the detection multiple means a background ratio**
(P3-T5, ADR 0016 §5) — a third quantity needs a third colour language, not a
borrowed one. Open the PNG. P3-T5's title clipped at both canvas edges on first
render *despite* the plan warning about it.

---

## Reconnaissance — NOT A RESULT (hard constraint 5)

Computed in a scratchpad against `results/interim/expr_q3.tsv` and
`obs['negprobe']` at the ADR 0007 rule, to size Phase 4 before committing to it.
**None of it is a result, none of it may be cited, and P4-T2's rule must
reproduce any of it that ends up mattering.** ADR 0008's Gate 0 reconnaissance
is the precedent for this being legitimate *and* for it being non-citable.

**Genes clearing 50% detection in both members of an adjacency:**
~5,560 (lung `L` + `TIME-L`) and ~5,507 (brain `LB` + `TIME-B`), out of 18,694.
**Phase 4 is not starved of genes in general.** It is starved of a specific
class:

| gene | `L` 30 | `TIME-L` 15 | `LB` 27 | `TIME-B` 8 | |
|---|---|---|---|---|---|
| VEGFA | 30 | 15 | 26 | 7 | membrane/matrix — fine |
| CD47 | 29 | 14 | 24 | 7 | fine |
| CXCL12 | 29 | 15 | 22 | 7 | fine |
| NECTIN2 | 28 | 13 | 25 | 8 | fine |
| TGFB1 | 27 | 15 | 22 | 8 | fine |
| NRP1 | 24 | 15 | 20 | 8 | fine |
| CXCL9 | 17 | 13 | 5 | 5 | brain-limited |
| SIRPA | 6 | 12 | 5 | 6 | immune-side only |
| TNF | 7 | 6 | 4 | 3 | **below floor everywhere** |
| CCL21 | 5 | 10 | 1 | 0 | **below floor** |
| CCL19 | 3 | 9 | 2 | 1 | **below floor** |
| CD80 | 0 | 6 | 3 | 3 | **below floor** |
| CXCL10 | 1 | 4 | 1 | 1 | **below floor** |
| CXCL11 | 3 | 3 | 2 | 1 | **below floor** |
| IL1B | 1 | 3 | 2 | 0 | **below floor** |
| IFNG | 0 | 3 | 2 | 0 | **below floor** |
| CCR7 | 0 | 5 | 0 | 2 | **below floor** |
| CXCR3 | 0 | 2 | 1 | 1 | **below floor** |
| IL10 | 2 | 0 | 3 | 0 | **below floor** |
| IL12B | 0 | 0 | 1 | 0 | **below floor** |
| LGALS9 | — | — | — | — | **not on the panel** |

**Read this the way ADR 0014 would.** The soluble-immune axis — interferons,
interleukins, the inflammatory chemokines and their receptors — is gone at both
sites, exactly as ADR 0014 said. What survives is the **membrane-bound and
matrix-associated** axis: VEGFA–NRP1, CD47–SIRPA, TGFB1, NECTIN2, CXCL12. Those
are real ligand–receptor biology and they are measurable here.

**Two consequences for how Phase 4 gets written:**

1. **The expected Gate 4 outcome is no longer "probably nothing survives."** It
   is "a structurally biased subset survives, and the bias is the finding" —
   which is the same shape as Gate 3's answer, where the audit was the larger
   result. Plan the write-up for that, not for a bare null.
2. **A nomination list dominated by VEGFA/TGFB1/CD47 is a detection artefact
   waiting to be over-read.** These are high-expression, broadly-expressed genes;
   their correlation across 13 patients may be driven by shared AOI quality, not
   by crosstalk. **Consider pre-registering a negative control in P4-T4** — e.g.
   permuted pairs restricted to matched-expression deciles — so "high-expression
   genes correlate with high-expression genes" is separable from a nomination.
   That is a model decision, so it belongs in ADR 0021 or its own ADR, decided
   before the ranking is seen.

---

## What to reuse

| need | reuse |
|---|---|
| detection per gene × group | the ADR 0007 rule as lifted in `workflow/scripts/checkpoint_detection.py` and `score_signatures.py` — **do not write a third copy** |
| keying AOIs | **`aoi_code`, never `compartment`.** `compartment` is degenerate: `L`, `LB` and `mLN` are all `tumour`, so grouping by it silently pools lung and brain |
| a permutation null | `workflow/scripts/variance_partition.py` (P1-T3) — structure and the `config["seed"]` convention |
| paired patient sets | `results/tables/design_matrix.tsv` (counts) and `obs.tsv` (`patient_id` × `aoi_code`); `checkpoint_paired_check.py` is the *cross-site* pairing and is the wrong axis for Phase 4, but its no-test-statistic assertions are the pattern worth copying |
| acquiring an external file under `resources/` | `workflow/scripts/external_validation.py` + `resources/supplementary_provenance.tsv` (ADR 0005, P3-T2b — the third sanctioned writer) |
| exact TSV reads in Python | `pd.read_csv(..., float_precision="round_trip")` — **mandatory**, ADR 0013 postscript |
| rendering not-assessable | `checkpoint_dotplot.py` (dots) or `contexture_heatmap.py` (hatched cells) — absence of evidence must be visible, never omitted |
| a network figure | nothing in-repo yet; `networkx` skill, and P3-T5's layout lessons (title on short lines; `main \| key \| colourbar`; open the PNG) |
| an external comparator at all | `results/tables/external_validation_summary.json` (ADR 0017) |
| mixed models, if Phase 4 needs one | `fit_checkpoint_models.R` is parameterised and drives two contrasts already — extend, do not copy |

---

## Stop and ask — Phase 4's live ones

1. **The P4-T1 pairing rule** — average vs. higher-QC for P12 and P24's duplicate
   `TIME-L`. **ADR 0021, before anything runs.** This is the plan's own explicit
   instruction, not an inference.
2. **The P4-T2 detection filter** — a pre-set threshold. Choosing it after seeing
   the ranked table is the failure ADR 0016 §3 exists to prevent.
3. **Substituting `LGALS9C` for the absent `LGALS9`** (P4-T6), or adding any gene
   to any analysis panel. Gene membership is a scientific decision — ADR 0011 and
   ADR 0015, and ADR 0015 exists *because* citation checking caught four errors
   the schema could not.
4. **Any change to the random-effects structure** if a model appears.
   `(1|patient_id)` is mandatory on non-independence grounds (ADR 0009); the
   config schema enforces it with a regex. Singular fits are not a reason to drop
   it.
5. **Dropping an AOI.** Flag, don't drop — `results/tables/qc_excluded.tsv` with
   a reason.
6. **Re-scoping A5.** ADR 0014 §3 explicitly rejected restricting Phase 4 to the
   lung `L ↔ TIME-L` adjacency, on the grounds that it answers a question this
   project is not about. If the brain arm proves unworkable at n = 8, that is a
   new decision with a new ADR, not a quiet narrowing.

---

## Gate 4

> **GATE 4** — Empirical FDR computed and nominations survive it. If nothing
> survives, report that: "no LR pair exceeded chance expectation at n=13" is an
> honest, useful, publishable-to-blog result, and it's a better outcome than a
> ranked list you can't defend.

**The gate turns on P4-T4, not on P4-T5.** The pass condition is *the empirical
FDR being computed and honoured* — a surviving list and an empty list both pass;
a ranked list without a null does not. Do not let P4-T5 drift into the gate's
subject.

Write the Gate 4 record as an ADR (0022 or later, after ADR 0021), the way
ADR 0014 and ADR 0020 record Gates 2 and 3, and squash-merge `P4` into `main`
with the verdict in the message (ADR 0004). `Markdowns/PROJECT_PLAN.md` is
gitignored, so the ADR is the durable record of both the gate and of any defect
found in the plan's own §6 text — ADR 0014 §4 is the precedent for recording
those.

---

## State you'll have forgotten

- **`--conda-prefix "$HOME/nsclc-envs"` is MANDATORY on this machine**, alongside
  `--use-conda`. The symlink points *at* `.snakemake/conda`, so nothing moves —
  but this working directory contains a space, and **three** separate
  conda/bioconda scripts interpolate the prefix unquoted. Recreate with
  `ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"`. **Do not** relocate the
  envs: that fires the software-env trigger on `p0t2_fetch_geo` and its
  `protected()` outputs abort the DAG. ADR 0013.
- **A bare `snakemake -n` is not a clean dry run** — without `--use-conda` the
  same trigger fires. And always run from the `nsclc_bm_spatial` env
  (Snakemake 8.30); base anaconda's 9.20 writes `.snakemake/metadata` in a format
  8.30 reads as stale, into the same protected-output abort.
- **`config/config.yaml` may now be edited — once, deliberately, and early.**
  Phase 4 needs `phases.crosstalk: true` plus a `crosstalk:` block, and its
  SHA-256 goes into `uns`, so the edit rebuilds the `.h5ad` and re-runs Phases
  1–3. **Do it in the first Phase 4 commit, with all of Phase 4's config keys
  present**, the way P3-T1 paid it once in `f1f55b8` — and check the rebuild the
  same way, by confirming the existing tables come back byte-identical except
  `h5ad_summary.json`'s `git_sha` and `config_sha256`. P3-T1's rebuild returned
  36 of 37 identical; a second diff is a bug, not noise.
- **A clean Phase 4 dry run is `p4tN… + all` only.** If `snakemake -n` proposes
  Phase 0/1/2/3 work *after* the config edit has been absorbed, something else
  was touched — stop rather than spend 20 minutes.
- **`rule p2t0_repair_r_env` must run before any R rule.** Declare
  `results/interim/env_repair/<env>.ok` as an input of every R rule. Both R envs
  are built: `r-stats` (lme4 2.0.6 / Matrix 1.7.5) and `r-geomx`. `r-stats.yaml`'s
  header says "Phases 2-3" — **if Phase 4 uses R, that header is now wrong** and
  the env may need a Phase 4 line. Phase 4 is mostly Python (correlation,
  permutation, `networkx`), so it may need no R at all.
- **A green `--conda-create-envs-only` is not evidence an env works.** Load the
  libraries.
- **Neither R's nor pandas' default float parser is correctly rounded**
  (ADR 0013 §3 and its P3-T2 postscript). Exactness is a property of the
  *reader*. Assert structure across the boundary; check numbers end-to-end.
- **`py-analysis.conda-lock.yml` is stale** and `conda-lock` is not installed.
  Still the P6-T1 blocker. **If P4-T2 adds a Python dependency, it must not go in
  `py-analysis.yaml`** (see P4-T2) — a new env is cheaper than a rebuilt Phase 0.
- **P6-T1 must be run on a path containing a space**, or it will not exercise
  `p2t0_repair_r_env` and will report a false pass.
- **`results/reports/landscape_explorer` renders on a black background** — P1-T5's
  notebook, deliberately left alone because editing it fires that rule's rerun
  trigger. The fix is the `plt.rcParams` reset in
  `notebooks/apps/checkpoint_explorer.py`. Worth doing when something else
  touches Phase 1.
- **`Markdowns/PROJECT_PLAN.md` is gitignored** — ADRs are the durable record.

## Useful commands

```bash
conda activate nsclc_bm_spatial
ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"   # once per machine
git switch -c P4                                     # from bce2bde, ADR 0004
snakemake -n --use-conda --conda-prefix "$HOME/nsclc-envs"   # must stay clean
snakemake --lint                                             # must stay clean
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" --cores 4
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" \
  --report results/reports/workflow-report.html
grep -rniE "colocali[sz]|spatial(ly)? (proximit|adjacen|close)" \
  --exclude-dir=.git --exclude-dir=.snakemake .          # P4-T7, must come back clean
/gate 4
```

---

## Phase 3, for reference

Gate 3 (ADR 0020) passed on two clauses satisfied by different tasks.

**An interpretable compartment-resolved pattern:** `CD274` (PD-L1) is enriched in
the **lung immune** compartment, **+0.612 SD [+0.361, +0.863], q = 0.0001**,
`expression ~ compartment + (1|patient_id)`, 45 AOIs across 30 patients, detected
19/30 in `L` and 13/15 in `TIME-L`. Survives both pre-registered sensitivities
(batch +0.599, background +0.635); lme4 vs. statsmodels agree to 1.1e-6.

**The audit, and it is the larger result:** **40 of 63 gene × compartment cells
are below the pre-registered floor** and the panel clears it in **`TIME-L`
alone** (`TIME-L` 8/9 genes, `TIME-B` 4/9, `LB` 3/9, `mLN` 3/9, `L` 2/9, `TBME`
2/9, `BC` 1/9).

**What Phase 3 deliberately did not say:** the lung-vs-brain shift is null in
every assessable gene — largest `CD274` at −0.260 [−0.595, +0.075], q = 0.462,
**`TIME-B` n = 8** — and against the 1.1–1.3 SD power floor those nulls are
**uninformative, not negative**.

**Three facts Phase 3 established that Phase 4 inherits:**

1. **This project's inputs are demonstrably the published inputs** (ADR 0017) —
   all 2,243,280 values of `layers['q3']` identical to the paper's deposited
   Source Data, and the 119-vs-120 gap closed (the paper dropped `TBME15b`).
2. **Detection and expression move OPPOSITE ways under the same background
   gradient** (ADR 0018): detection is `q3 > 2 × negprobe`, so higher background
   → looks depleted; expression is `log2(q3 + 1)`, so higher background → looks
   enriched. **Phase 4 correlates expression across patients, so the relevant
   question is whether AOI-level background varies with patient** — if it does,
   two background-tracking genes correlate for no biological reason. Measure it;
   do not assume it either way. P3-T4 measured its own version and found the site
   gradient negligible, which is not the same gradient.
3. **`TBME`, `mLN` and `BC` each sit wholly within one DSP run**, so batch is
   inseparable from biology there (Q3, ADR 0009 §3). Phase 4's adjacencies are
   `L ↔ TIME-L` and `LB ↔ TIME-B`, which are batch-estimable — **`TBME` is not
   one of Phase 4's compartments**, despite being the glial side of a brain
   adjacency, and ADR 0014 §3 records why.

## Phase 2, for reference

Brain minus lung, within `TIME`, `TIME-B` n = 8: antigen presentation
**−0.477 SD [−1.030, +0.075]** (q = 0.098, 13/13 detected at both sites),
cytotoxicity **−0.946 SD [−1.590, −0.302]** (q = 0.018), myeloid M2 null.
Three of six signatures failed their coverage floor — `exhaustion` 1/6 and `tls`
1/5 in brain, `myeloid_m1` 2/10 at *both* sites — and the large nominally
significant "reductions in brain" among them are **not reportable**. The genes
that failed are the cytokine genes, which is the ADR 0014 finding Phase 4
inherits and the Reconnaissance section above refines.
