# Next steps — Phase 5 (prognostic association), setup complete

**Last session:** 2026-09-15
**Branch:** `P5`, branched from `main` at ADR 0004's one-branch-per-phase rule.
*(No SHA is quoted here on purpose: a commit hash in this file is stale the
moment the next commit lands, and a stale one reads as authoritative. Run
`git log --oneline -1`.)*
**Gate 0:** PASSED. **Gate 1:** PASSED. **Gate 2:** PASSED (ADR 0014).
**Gate 3:** PASSED (ADR 0020). **Gate 4:** PASSED (ADR 0022).
**Gate 5:** not yet reached — it is *"timebox respected"*, not "a result".
**Next free ADR number: 0025.**

**The Phase 5 setup block is DONE.** All three items that gated it are closed,
and the three constraints on the join are now enforced by the manifest rather
than remembered:

| was | now |
|---|---|
| MOESM4 absent from the repo | fetched, pinned, **fifth sanctioned writer** at `resources/clinical/` |
| `phases.survival: false` | `true`, with a complete `survival:` block, rebuild paid once |
| no schema for a `survival:` block | `config.schema.yaml` has one; `clinical.schema.yaml` is new |

**The join gate PASSES, and a rule says so.** `p5t1_clinical_join` re-derives
both `docs/data-provenance.md` §Q5 assertions on every run: all 35 GEO patients
appear in Supplementary Data 1, and all 113 non-control AOI-code suffixes equal
their patient number. It carries **no fuzzy matching and no fallback path** —
the shape `lr_nominations.py` has, so an honest abandonment stayed reachable
rather than escapable. 9 of 9 assertions pass. Outputs:
`results/tables/clinical_cohort.tsv` (44 rows, 35 `in_geo`, 16 `in_phase5`) and
`clinical_join_summary.json`.

---

## START HERE — Phase 6, the primary deliverable

**Phase 5 is COMPLETE and Gate 5 PASSED (ADR 0028).** Nothing is outstanding.

**P6-T1 is the clean-room reproduction** and its blocker was resolved at
ADR 0023: fresh clone, fresh envs, `snakemake --cores 4 all` from nothing,
timed. Three things carry into it, all already recorded further down this file:
the pins are `osx-arm64` only, regenerating a pin is mandatory after any env
YAML edit, and **the run must happen on a path containing a space** or it will
not exercise `p2t0_repair_r_env` and will report a false pass.

`resources/` now has **six** sanctioned writers — `raw/`, `reference/`,
`supplementary/`, `ligand_receptor/`, `clinical/`, `tcga/` — and **P6-T2's
RO-Crate must describe all six**, each with its `*_checksums.sha256` and
`*_provenance.tsv` pair.

---

## What Phase 5 found, stated correctly

**0 of 14 pre-registered tests survive BH at FDR 0.05.** Smallest q = **0.2563**
(TCGA `antigen_presentation`, HR 0.70 [0.52, 0.94], raw p 0.018, n = 502 with
182 events).

| | |
|---|---|
| cells considered | **18** (12 GeoMx signature × arm, 6 TCGA) |
| tested | **14** (8 GeoMx assessable + 6 TCGA) |
| not assessable | **4** — `exhaustion` brain, `tls` brain, `myeloid_m1` both sites |
| survive FDR 0.05 | **0** |

**Only q is reportable — never a raw p from P5-T3 or P5-T4.** The TCGA
`antigen_presentation` signal is the strongest in the phase and it lands at
q = 0.26: a hypothesis, not a finding. P5-T5 exists so that raw p = 0.018 cannot
be quoted alone.

### The one distinction P6-T7 must not collapse

**The two nulls are not the same kind of null.**

- **This study**: 12 lung and 7 brain patients fitted, splits 6/6 and 3/4,
  intervals spanning eight- to thirty-fold. A null is **uninformative, not
  negative** (ADR 0024 §7).
- **TCGA**: 502 patients, 182 events, intervals ~0.5–1.2 wide. A null is a
  **much stronger** statement — though still not proof of absence.

**That asymmetry is the deliverable, not any q-value.** Writing "no signature was
prognostic" without distinguishing the two cohorts' power converts a measurement
of this study's *size* into a claim about *biology* — the same error ADR 0022
guards against for Phase 4.

And the comparison is **bounded structurally** (ADR 0027 §5): GeoMx measures the
**PanCK-negative segment**, TCGA **whole bulk tumour**, so a disagreement is not
necessarily about biology — and **there is no brain comparator at all**, so
nothing in Phase 5 externally checks `TIME-B` n = 8.

One live example: the cohorts' `antigen_presentation` point estimates point
**opposite ways** (GeoMx HR 2.19 [0.65, 7.35] vs TCGA 0.70 [0.52, 0.94]). **Not a
contradiction** — the GeoMx interval spans 1 and is consistent with both
directions, on a different measurement. "The cohorts disagree" would be wrong
twice over.

### Three things Phase 5 learned that Phase 6 should not relearn

1. **A declared family size is a number someone wrote down.** ADR 0024 §8's "12"
   forgot that §2 of the same ADR splits this study into two arms. The true
   family is 14 and P5-T5 **counts it from the upstream tables**. The direction
   matters: an undercounted denominator makes every q **too small**, and a
   q-value carries no evidence of the denominator behind it.
2. **Ties are where an FDR implementation goes wrong quietly** — Phase 4's
   lesson, and this family has **two tied pairs**. The correction is computed by
   `statsmodels` *and* an independent hand-rolled step-up (agreeing to 1.1e-16)
   with monotonicity asserted directly.
3. **Open the PNG — four occurrences now.** P5-T3's annotations collided with
   the curves (an interval read as a different number); P5-T4's suptitle ran off
   **both** edges, because matplotlib neither wraps nor clips one. Both were
   caught only by looking.

---

## The number that governs Phase 5

**n = 16, not ~35 — and 19 events across both arms.**

`signature_scores.tsv` covers 13 lung + 8 brain patients, 5 shared. PROJECT_PLAN
§6 P5-T3's "at ~35 patients this is descriptive" counts patients with a *tumour*
AOI, **none of which carries a signature score**. ADR 0024 §2 records why Phase 5
does not chase the larger number: reaching it means scoring `L`/`LB`, a
stop-and-ask scope change that re-opens the detection floor where the panel
clears far less (P3-T2: `L` 2/9 genes, `LB` 3/9, against `TIME-L` 8/9). **A
larger n bought with a less measurable compartment is not a larger n.**

Measured by `p5t1_clinical_join`, not quoted: **lung 13 patients / 12 events / 1
censored; brain `TIME-B` n = 8 / 7 events / 1 censored.** After ADR 0025's
exclusion of the one patient with no follow-up time, **12 and 7 enter a fit** and
the median split is **6/6 and 3/4** — the brain arm exactly on
`min_arm_size = 3`, with no patient to spare.

**So a Phase 5 null is uninformative, not negative** — as the 1.1–1.3 SD power
floor makes a Phase 2 null uninformative and ADR 0022 makes Phase 4's empty table
an assay limit. Fixed in ADR 0024 §7 **before the first curve**, so it cannot be
softened later if something happens to separate. A separation at 6 versus 7
patients describes this cohort; it is never an inferential claim.

**Three of the six signatures failed their Phase 2 detection floor** —
`exhaustion` 1/6 and `tls` 1/5 in brain, `myeloid_m1` 2/10 at *both* sites. They
are scored and corrected within the family and reported as **"not assessable"**
wherever they failed, never as a prognostic null. ADR 0008's rule is not
suspended by a change of outcome variable: a score built from genes at background
measures background, whatever it is regressed against.

---

## The one sentence Phase 4 earned

> No inferred ligand–receptor pair exceeded chance expectation in either
> adjacency — lung n = 13, brain **`TIME-B` n = 8** — and the assay resolves
> matrix-associated interactions roughly **eight times** as often as secreted
> ones.

Both halves matter and the second is the larger one. **The null is an
assay-sensitivity limit, never evidence that the crosstalk is absent.** Anyone
writing this up (P6-T7) who states the first half without the second has
converted a measurement of the assay into a claim about biology.

---

## Which phase is next — decide this before writing code

PROJECT_PLAN's slip rule: *"if you're behind at the end of week 4, cut Phase 5
entirely and protect Phase 6. A finished, reproducible, well-documented
four-phase project beats a sprawling six-phase one that nobody can run."*

**Phase 5 is a stretch goal with a one-week timebox and Gate 5 is "timebox
respected".** Q5 resolved **yes**, and its setup and join are now done (ADR 0024,
P5-T1) — what remains is P5-T2..T5, at the top of this file. The binding
constraint is no longer availability; it is **n = 16 with 19 events**.

**Phase 6 is arguably the primary deliverable** (PROJECT_PLAN §6's own words).
**P6-T1's blocker is resolved** (ADR 0023) — the environments are pinned and
verified, so the task is no longer gated on anything. What remains of it is the
clean-room run itself: a fresh clone, fresh envs, `snakemake --cores 4 all` from
nothing, timed.

**One caveat carries forward into that run:** the pins are `osx-arm64` only, so
a clean-room reproduction on Linux solves the YAMLs unpinned. That is no worse
than before, but it is not a guarantee — see below.

---

## The P6-T1 blocker is RESOLVED (ADR 0023), and the note about it was wrong

The long-standing note read "`py-analysis.conda-lock.yml` is stale and
`conda-lock` is not installed". **Both halves were false and the truth was
worse:** `conda-lock` was installed all along, and Snakemake never read
`*.conda-lock.yml` at all — it reads `<env>.<platform>.pin.txt`
(`deployment/conda.py:128`). **Nothing was pinned.** Two of the three locks were
also stale in content, missing `xarray` and `r-reformulas`, both of which were
added to fix real failures — so wiring them up would have reproduced two known
bugs in the very run meant to prove the pipeline works.

**Now fixed.** `workflow/envs/<env>.osx-arm64.pin.txt` exists for all three
envs, generated with `conda list --explicit --md5` from the environments that
actually produced every committed result. Verified by forcing one rule per env
to re-execute: all three outputs byte-identical, including the R mixed-model
fit. The three dead `*.conda-lock.yml` files are removed.

**What remains for P6-T1:**

- **The pins are `osx-arm64` only.** On `linux-64` Snakemake finds no pin and
  solves the YAML, exactly as it did before — nothing is lost, but nothing is
  guaranteed either. **P6-T6's limitations document must say so.** If
  cross-platform reproduction is wanted, render a `linux-64.pin.txt` from a
  `conda-lock` solve; `conda-lock` is still in `environment.yml` for that.
- **Regenerating a pin is now mandatory after any env YAML edit**, because a
  stale pin installs the wrong environment silently instead of failing.
- **P6-T1 must be run on a path containing a space**, or it will not exercise
  `p2t0_repair_r_env` and will report a false pass. This working directory has
  one ("Biomedical Data Science"); a clean clone elsewhere may not.
- **`r-geomx.yaml`'s header is stale**: it says "Phase 0" but
  `p4t2c_export_lr_database` now runs in it too. Editing it fires the
  software-env trigger on every rule in that env AND invalidates the pin, so it
  needs the ADR 0023 transition (rebuild, re-repair, `--touch`) — do it in a
  commit where that is the point, not casually.

---

## What Phase 4 left for Phase 6 to use

| need | where |
|---|---|
| the Gate 4 record and every Phase 4 decision | `docs/decisions/0022-gate-4-record.md`, `0021-phase-4-pre-registration.md` |
| the detection-by-class table (the headline) | `results/tables/crosstalk_lr_summary.json` → `by_site` |
| the honest null | `results/tables/crosstalk_fdr.tsv`, `crosstalk_null_summary.json` |
| the empty nomination table, which is a result | `results/tables/crosstalk_nominations.tsv` |
| the concordance prose for P6-T7 | `docs/analysis-notes.md`, P4-T6 section |
| figures, both `report()`-wrapped with captions | `crosstalk_null.png`, `crosstalk_network.png` |
| a build-time language check | `p4t7_language_audit` — **runs on every build**, so P6-T7's prose is checked automatically |
| a permutation-null implementation | `workflow/scripts/lr_null_calibration.py` — identity permutation reproduces the observed rho to 1e-16 |

---

## Things Phase 4 learned that the next phase should not relearn

1. **Open the PNG.** Both Phase 4 figures clipped on first render — the network
   showed "CEACAM" for both CEACAM1 and CEACAM5, and the colourbar label went
   off-canvas. P3-T5 had the same failure and the plan warned about it. Use the
   `main | colourbar` gridspec from `checkpoint_dotplot.py`; **do not use
   `tight_layout` with an attached colourbar** — it warns and then clips.
2. **Assertions catch what review does not.** The FDR step-down was wrong under
   ties (Spearman on 8 patients takes 42 distinct values across 260 rows) and
   would have shipped a plausible ranked table. The monotonicity assertion
   caught it. Write the assertion that would fail if you were wrong.
3. **A hand-transcribed orientation has nothing checking it.** `CD96–NECTIN1`
   was entered backwards into the comparator and reported as absent from a
   database that contains it. Where the underlying measure is undirected, match
   both orientations.
4. **A grep that cannot tell a prohibition from a claim is unrunnable.**
   PROJECT_PLAN §6 P4-T7's literal instruction fails on `CLAUDE.md`'s own hard
   constraint 6, because the files that state a rule must name what they forbid.
5. **An explicit-looking seed can still be an implicit one.** The permutation
   null was seeded `default_rng([seed, hash(site) % 2**31])`, which reads as
   correct and is not: **Python salts string hashing per process**, so two runs
   of the same commit gave different nulls (min FDR 0.228 vs 0.225). Derive
   stream ids from a stable digest, never from `hash()` of a string. It survived
   review and was caught only because a branch operation forced a re-run.
6. **A check with no rerun trigger reports a stale pass.** `p4t7_language_audit`
   declared no inputs, so it would have run once and never again — and its green
   output would have kept asserting something nobody re-tested. It now keys on a
   digest of every file it scans.
7. **Every THRESHOLD belongs in a sibling YAML, not in `config/config.yaml`.**
   The Phase 4 edit cost a full Phases 1–3 rebuild, verified at 58 of 59 tables
   byte-identical. Note this is *not* the same as "config.yaml will never be
   edited again": a phase switch lives there, so **Phase 5 must edit it and will
   pay the rebuild once**, exactly as Phases 3 and 4 each did. Pay it once, with
   every key present, and verify it the same way.

---

## State you'll have forgotten

- **`--conda-prefix "$HOME/nsclc-envs"` is MANDATORY on this machine**, alongside
  `--use-conda`. Recreate with `ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"`.
  **Do not** relocate the envs: that fires the software-env trigger on
  `p0t2_fetch_geo` and its `protected()` outputs abort the DAG. ADR 0013.
- **A bare `snakemake -n` is not a clean dry run.** Always run from the
  `nsclc_bm_spatial` env (Snakemake 8.30); base anaconda's 9.20 writes
  `.snakemake/metadata` in a format 8.30 reads as stale.
- **`resources/` now has FIVE sanctioned writers** — `raw/` (P0-T2),
  `reference/` (P2-T5), `supplementary/` (P3-T2b), `ligand_receptor/` (P4-T2)
  and **`clinical/` (P5-T1)**. Each has a `*_checksums.sha256` and
  `*_provenance.tsv` pair, and P6-T2's RO-Crate must describe all five. A
  **sixth**, `tcga/` (P5-T4), is named in `config.yaml` but not yet created.
- **`results/reports/landscape_explorer` still renders on a black background** —
  P1-T5's notebook, left alone because editing it fires that rule's rerun
  trigger. The fix is the `plt.rcParams` reset in
  `notebooks/apps/checkpoint_explorer.py`. Worth doing when something else
  touches Phase 1, and P6-T4b is that occasion.
- **Regenerate the pin after ANY edit to an env YAML** (ADR 0023). A stale pin
  installs the wrong environment *silently* rather than failing, which is the
  same failure shape as the blocker just cleaned up. The edit also invalidates
  the env hash, so it needs the full transition: `--conda-create-envs-only`, a
  **scoped** `-f` on the two `env_repair` sentinels (a bare `--forcerun` aborts
  on protected outputs), then `--touch`.
- **`Markdowns/PROJECT_PLAN.md` is gitignored** — ADRs are the durable record.
  ADR 0022 §5 records two further defects in §6's Phase 4 text.

## Useful commands

```bash
conda activate nsclc_bm_spatial
ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"   # once per machine
snakemake -n --use-conda --conda-prefix "$HOME/nsclc-envs"   # must stay clean
snakemake --lint                                             # must stay clean
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" --cores 4
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" \
  --report results/reports/workflow-report.html
/gate 5
```

---

## Phase 4, for reference

**Gate 4 (ADR 0022) passed on its second clause.** The empirical FDR was
computed and honoured; **0 of 268 lung and 0 of 260 brain** primary
direction-rows clear FDR 0.05, smallest 0.499 and 0.224. Observed |rho| sits
**on** the permuted-pairing null. The gate turns on P4-T4, not P4-T5:
`lr_nominations.py` was written with no fallback and no second threshold, so an
empty table was reachable rather than escapable.

**The larger result is the detection gap by class**, ADR 0014's claim made
numeric on 2,239 interactions: ECM-Receptor **45.3%** brain / **39.6%** lung
admitted, Cell-Cell Contact 13.7% / 17.3%, Secreted Signaling **6.4%** / **5.3%**.

**Background is a patient-level property** — across a patient's two paired AOIs,
`negprobe_log2` ρ +0.548 brain / +0.352 lung, `gene_detection_rate` ρ +0.667 /
+0.231. Three of four above threshold, which is why ADR 0021 §6's
abundance-matched Null B was pre-registered. P3-T4's *site* gradient was
negligible; that is a different gradient.

**`VEGFA–NRP1`, which this file previously called the one informative
comparator, is absent from CellChatDB v2 entirely** — it carries VEGFA to FLT1
and KDR, with NRP1 only as a SEMA3x co-receptor. A database difference, not a
detection one. Four of the five comparators are unevaluable for four different
reasons.

## Phase 3, for reference

`CD274` (PD-L1) is enriched in the **lung immune** compartment, **+0.612 SD
[+0.361, +0.863], q = 0.0001**, 45 AOIs across 30 patients. **40 of 63 gene ×
compartment cells are below the pre-registered floor** and the panel clears it in
`TIME-L` alone. The lung-vs-brain shift is null in every assessable gene —
largest `CD274` at −0.260 [−0.595, +0.075], q = 0.462, **`TIME-B` n = 8** — and
against the 1.1–1.3 SD power floor those nulls are **uninformative, not
negative**.

## Phase 2, for reference

Brain minus lung, within `TIME`, `TIME-B` n = 8: antigen presentation
**−0.477 SD [−1.030, +0.075]** (q = 0.098, 13/13 detected at both sites),
cytotoxicity **−0.946 SD [−1.590, −0.302]** (q = 0.018), myeloid M2 null. Three
of six signatures failed their coverage floor; the genes that failed are the
cytokine genes, which Phase 4 has now measured across a whole LR database.
