# Changelog

All notable changes to this project are recorded here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added — P5-T5: multiplicity honesty, and GATE 5 PASSED (ADR 0028) (2026-09-15)
- **GATE 5 PASSED.** Gate 5 is *"timebox respected"*, not "a result was found".
  Phase 5 ran P5-T1 through P5-T5, produced every output through a rule, and
  stopped where it was designed to stop. **Phase 5 is complete.**
- **0 of 14 pre-registered tests survive Benjamini–Hochberg at FDR 0.05.**
  Smallest q = **0.2563** (TCGA `antigen_presentation`, HR 0.70 [0.52, 0.94],
  raw p = 0.018, n = 502 with 182 events). **Only q is reportable** — the raw
  p cannot be quoted alone, which is the entire reason P5-T5 exists.
- **The honest accounting: 18 cells considered, 14 tested, 4 not assessable,
  0 surviving.** The four are `exhaustion` brain, `tls` brain and `myeloid_m1`
  at both sites — kept **out of the denominator**, because ADR 0008 forbids
  reporting them as prognostic nulls at all and correcting over them would
  deflate every q.
- **The family is COUNTED from the upstream tables, never quoted.** ADR 0026
  exists because ADR 0024 §8's "12" forgot that §2 of the same ADR splits this
  study into two arms that are never pooled. The direction of that error is what
  matters: an undercounted denominator makes every q **too small**, and a
  q-value carries no evidence of the denominator behind it.
- **The correction is computed twice** — `statsmodels.multipletests` and an
  independent hand-rolled BH step-up, agreeing to **1.1e-16** — with
  monotonicity asserted directly. Phase 4 shipped a step-down that was
  correct-looking and wrong under ties, and **this family contains two tied
  pairs**, so the guard was not hypothetical. 6 of 6 assertions pass.
- Sensitivities (`zscore` in both cohorts, TCGA's continuous Cox) stay **out of
  the family** (ADR 0012): a second estimator does not double the search, it
  measures whether the primary is estimator-dependent.

#### The distinction P6-T7 must not collapse
- **The two nulls are not the same kind of null.** This study: 12 and 7 patients
  fitted, intervals spanning eight- to thirty-fold — **uninformative, not
  negative** (ADR 0024 §7). TCGA: 502 patients, 182 events, intervals ~0.5–1.2
  wide — a **much stronger** statement. **That asymmetry is the deliverable**,
  not any q-value. Stating "no signature was prognostic" without distinguishing
  the cohorts' power converts a measurement of this study's size into a claim
  about biology.
- The comparison is bounded structurally (ADR 0027 §5) — PanCK-negative segment
  vs whole bulk tumour, and **no brain comparator exists at all**. The cohorts'
  `antigen_presentation` estimates point opposite ways and that is **not** a
  contradiction: the GeoMx interval spans 1, on a different measurement.

### Added — P5-T4: the TCGA LUAD external check (ADR 0027) (2026-09-15)
- **`resources/tcga/` is the project's SIXTH sanctioned writer**, on ADR 0005's
  terms: `protected()` outputs, pinned digests, `pin_status: verified`.
- **`p5t4c_tcga_survival` scores and fits TCGA LUAD at n = 502 with 182 events**,
  median split 251/251, same split rule and tie side as P5-T3 because
  comparability is the point. 7 of 7 assertions pass.
- **The result that makes Phase 5 interpretable:** one of six signatures clears
  its own interval — `antigen_presentation`, **HR 0.70 [0.52, 0.94], raw
  p = 0.018, uncorrected** (higher score, better survival). The other five span
  1. So **this study's null is a power limit, not evidence of absence** —
  exactly what PROJECT_PLAN §6 said P5-T4 was for.
- **The family is now fully determined at 14**: 8 (P5-T3) + 6 (TCGA), which
  discharges ADR 0026 §1's deferral. P5-T5 corrects once over the union; no
  `q_` column exists in either upstream table and both assert it.
- **Source switched to UCSC Xena, and it was forced.** cBioPortal's bulk files
  return **403** (datahub) and **404** (LFS mirror); only its live REST API
  answers, which serves per-study JSON and cannot carry a pinned SHA-256 under
  ADR 0005. ADR 0024 §1 pre-authorised this fallback in terms; ADR 0027 records
  the substitution with the HEAD evidence.
- **All six signatures are tested, and that is a decision.** Phase 2's detection
  floor is `q3 > 2 × NegProbe-WTX` — a property of the *GeoMx assay* — and TCGA
  has no negative-probe channel, so restricting would import a limitation the
  external cohort does not have. `myeloid_m1` is labelled **TCGA-only, no GeoMx
  comparator**. A **gene-presence** table is computed instead and labelled as
  presence, never detection: all 47 genes are present.
- The cohort is **constructed, not taken as given**: only sample type `01`
  survives (the matrix carries **59 `-11` normals and 2 `-02` recurrences**),
  9 cases with missing follow-up and 4 with non-positive time are excluded and
  counted — ADR 0025's shape applied to an external cohort — and `OS.time` is
  converted from days at `365.25 / 12`, declared in the manifest.

#### Two limitations that bound every number P5-T4 produces (ADR 0027 §5)
- **It is not the same measurement.** GeoMx `TIME-L` is the PanCK-negative
  *segment* of an ROI; TCGA is *whole bulk tumour*. A disagreement between the
  cohorts is not necessarily a disagreement about biology, and an agreement is
  not necessarily a replication.
- **There is no brain comparator at all.** TCGA LUAD is primary lung, so this
  bounds the **lung** null only. Nothing here speaks to `TIME-B` n = 8, the
  binding constraint on the project. P6-T6 must say so.

#### Caught by looking at the output
- **The median split and the continuous Cox disagree in both directions** — 2 of
  6 cells. `antigen_presentation` clears on the split and not continuously;
  `tls` the reverse (HR/unit 0.38 [0.16, 0.89], p 0.026). ADR 0027 §4 required
  this be readable, so it is now a flagged column and a summary field rather
  than only visible in the log.
- **The figure's suptitle ran off BOTH edges of the canvas** on first render —
  a ~150-character line, clipped to "n split…" and "…No brain comparator e".
  matplotlib neither wraps nor clips a suptitle; it simply draws past the
  figure. Split into three short lines. **Phase 4 lesson 1, fourth occurrence.**

### Added — P5-T3: Kaplan-Meier and log-rank, and a corrected family (ADR 0026) (2026-09-15)
- **`p5t3_survival_km` fits the 8 assessable signature × arm cells** — lung 5,
  brain 3 — by median split within arm, and emits models, an exploratory table,
  the KM curve points and a 12-panel figure. 6 of 6 assertions pass.
- **The result: no immune signature stratified survival.** Every hazard-ratio
  interval spans 1. Lung n = 12 (6/6), brain **`TIME-B` n = 8, 7 fitted** (3/4).
  The one cell with raw p < 0.05 — `tls` in lung, **HR 4.51 [0.88, 23.16],
  raw p = 0.042** — has an interval spanning 1 and is **not a detected effect**.
- **ADR 0026 corrects ADR 0024 §8's multiplicity arithmetic.** §8 declared
  "6 signatures × 2 cohorts = 12", which counted cohorts and forgot that §2 of
  the same ADR splits this cohort into two arms that are never pooled. The
  primary family is the **assessable** cells — **8 here** plus TCGA's, counted at
  P5-T4. An undercounted denominator makes every q too small, which is the
  direction that manufactures significance.
- **P5-T3 therefore emits NO q-values, and the rule asserts none exist.** The
  family spans this study and TCGA, so it cannot be corrected until both are
  fitted; P5-T5 corrects once over the union. Every p here is raw and
  uncorrected, in those words.
- Four of twelve cells are **not assessable** (`exhaustion` brain, `tls` brain,
  `myeloid_m1` at *both* sites) and get a row with no estimate — the
  `checkpoint_carrier_models.tsv` shape. They are **hatched on the figure rather
  than omitted**, so a reader can tell "not tested" from "tested, null".
- `sksurv.compare_survival` for the log-rank, `statsmodels` PHReg for the hazard
  ratio and CI. `lifelines` deliberately unused — adding it would fire the
  software-env trigger on `p0t2_fetch_geo`'s protected outputs (ADR 0023).

#### Three things caught by looking at the output, not by the assertions
- **The log-rank and the hazard-ratio interval disagree on `tls` lung** —
  score-test p = 0.042 against a Wald interval of [0.88, 23.16]. **Where they
  disagree the interval governs.** Now flagged in the table, log, summary and
  **on the figure panel**, so the p cannot be read off the plot alone. Hard
  constraint 7 is what made it visible.
- **The p-agreement tolerance was the wrong check.** It asserted the two
  p-values agree within 0.05 and passed at **0.046** — nearly failing on real
  score-vs-Wald divergence rather than on a bug. Replaced with an **exact** one:
  inverting the group indicator must negate the Cox coefficient (residual
  3.3e-16). Assert what must hold exactly, not what usually holds.
- **The figure's annotations collided with the curves on first render** — a step
  crossed the decimal point of "[0.88, 23.16]" so it read "23 16", with "(n=12)"
  struck through. Fixed with an opaque bbox and headroom. **Phase 4 lesson 1,
  third occurrence: open the PNG.**
- Identical results across different signatures are **not** a bug: brain
  `cytotoxicity` and `myeloid_m2` share a partition, and lung
  `antigen_presentation` and `exhaustion` differ only by swapping two patients
  with **identical survival times** (115 months, both events).

### Added — P5-T2: patient-level signature scores (ADR 0025) (2026-09-15)
- **`p5t2_patient_scores` collapses Phase 2's 276 per-AOI scores to 252
  patient-level rows** — 21 patient-sites × 6 signatures × 2 methods — and emits
  the analysis-ready table P5-T3 fits directly: score, coverage verdict and
  outcome in one row, so no joining happens inside a modelling script. 8 of 8
  assertions pass.
- Aggregation is `mean`, pre-registered in ADR 0024 and matching ADR 0021 §1's
  duplicate rule so Phase 5 invents no second convention. It touches **2 of 21
  patient-sites — P12 and P24, lung only** — asserted against ADR 0021 §1 rather
  than assumed; the other 19 pass through **bitwise unchanged**, also asserted.
- `ssgsea` is primary and `zscore` a declared sensitivity (ADR 0012's shape),
  both carried so P5-T3 can report whether a split is method-dependent.
- **Coverage is carried, never applied.** Four of twelve signature × site cells
  are below Phase 2's floor (`exhaustion` brain, `tls` brain, `myeloid_m1` at
  *both* sites). Filtering them here would hide a Phase 2 result; they are
  "not assessable", never a prognostic null (ADR 0008).

#### ADR 0025 — a censored patient with no follow-up time
- **ADR 0024 fixed the censoring convention but not its consequence.** Neither
  spelling carries a duration: `Alive` is a status and the brain cell is blank,
  and Supplementary Data 1 has no last-contact date. So this project holds an
  event status and **no time-axis position** for the three censored patients.
- **Patient 35 is the one in the Phase 5 cohort, and it is in BOTH arms.**
  Excluded from fits, recorded in a new `results/tables/survival_excluded.tsv`
  with an arm-specific reason, retained in the scores table with
  `has_followup = False` — flag-don't-drop, **never imputed**.
- Administrative censoring at the cohort maximum was **rejected**: it invents an
  observation *and* places it exactly where a median split is most sensitive,
  making the fabricated point the most influential one in the analysis.
- **Consequence, fixed before any fit: lung 12 and brain 7 enter a fit, splits
  6/6 and 3/4.** The brain arm lands **exactly on `min_arm_size = 3`** and has no
  patient to spare — one further exclusion makes the brain analysis not
  assessable rather than null. The rule fails rather than proceeds if it drops
  below.
- **No `config.yaml` edit** — this is a structural fact about the source data,
  not a tunable threshold, so no second Phases 1–4 rebuild was paid.

### Added — P5-T1: the clinical join, and its hard gate (2026-09-15)
- **`p5t1_clinical_join` re-derives both `docs/data-provenance.md` §Q5
  assertions on every run** and passes: all 35 GEO patients appear in
  Supplementary Data 1, all 113 non-control AOI-code suffixes equal their
  patient number. **No fuzzy matching and no fallback path** — the shape
  `lr_nominations.py` has, so an honest abandonment stayed reachable rather than
  escapable. 9 of 9 assertions pass.
- **Added** `results/tables/clinical_cohort.tsv` (44 rows, 35 `in_geo`, 16
  `in_phase5`) and `clinical_join_summary.json`. The nine supplementary rows
  without expression data are recorded as `in_geo = False`, not dropped
  (flag-don't-drop). Raw columns are preserved beside the derived
  `{arm}_months` / `{arm}_event` pair.
- The `.xlsx` reader is **lifted** from `external_validation.py` rather than
  imported, the way `checkpoint_detection.py` lifts its detection rule.
  `openpyxl` stays out of `py-analysis.yaml`.
- **Caught by assertion, not by review:** `pd.DataFrame` turns the brain
  endpoint's empty cells into `NaN`, and `NaN is None` is `False`, so a
  censoring test that works on the lung column (where censoring is the string
  `Alive`) silently counted the censored brain patient as a **death with a
  missing time**. The first run reported brain as 8 events / 0 censored against
  ADR 0024's pre-registered 7 / 1. The fix is `pd.isna`; the guard is a new
  per-arm assertion that event flags agree with the censoring set in **both**
  arms. The two columns spelling censoring differently is what made the bug
  asymmetric and survivable.
- Confirmed counts: **lung 13 patients / 12 events / 1 censored; brain
  `TIME-B` n = 8 / 7 events / 1 censored**, patient 35 censored in both. A
  censored patient has **no follow-up time in either column**, so the time is
  left missing rather than imputed — P5-T3 decides how to handle three patients
  with an event flag and no time.

### Added — Phase 5 setup: the clinical table, pre-registered (ADR 0024) (2026-09-15)
- **ADR 0024 pre-registers Phase 5 in eight sections**, committed before any
  Phase 5 number existed: the manifest placement, the cohort, the endpoints and
  censoring convention, the normalisation rule, the join gate, the split, what a
  null may say, and the multiplicity family. A median-split Kaplan–Meier at
  n = 13 and n = 8 is the easiest quantity in this project to manufacture after
  the fact, and Phase 5 has strictly more exposure than Phase 4 did.
- **Added `resources/clinical/` as the project's FIFTH sanctioned writer** —
  Supplementary Data 1 (`41467_2022_33365_MOESM4_ESM.xlsx`, 13,210 bytes),
  pinned in a new `config/clinical.yaml` with a new
  `workflow/schemas/clinical.schema.yaml`. Filed there rather than under
  `supplementary/`: those files are an external *comparator*, this is analysis
  *input*, and reusing `config/external_validation.yaml` would also have changed
  `p3t2b_external_validation`'s params hash and re-parsed a 45 MB workbook for
  no reason connected to Phase 5.
- **Opened `config/config.yaml → phases.survival` and the `survival:` block in a
  single edit**, every P5-T1..T5 key present including P5-T4's `tcga_yaml`, so
  the Phases 1–4 rebuild that a `config.yaml` edit forces is paid exactly once —
  as Phases 3 and 4 each did. `survival` added to `config.schema.yaml`, whose
  root is `additionalProperties: false`.
- **No conda environment was edited.** `scikit-survival` was already in
  `py-analysis.yaml`, so Phase 5 avoids ADR 0023's pin-regeneration transition
  entirely — and must keep avoiding it.

#### Three things the deposited file does that the docs did not record
- **The header cells carry embedded newlines** (`Age at \nNSCLC diagnosis`,
  `Primary lung cancer \ndiagnosis to death (Months)`) and ` Location of BrM `
  has a leading *and* trailing space. `docs/data-provenance.md` §Q5 transcribed
  the column list already-normalised, so it must not be matched verbatim — a
  caveat that was nowhere until now.
- **The two endpoint columns spell censoring differently.** The lung column
  carries the literal `Alive`; the brain column carries an **empty cell** for the
  same three patients (6, 11, 35). The sets are identical, and that identity is
  what licenses reading a blank as *censored* rather than *missing* — asserted,
  never assumed. Reading them as missing would drop the only censored
  observations in the study.
- **The sheet ends in a legend, not in data**: three blank rows then three
  treatment-abbreviation lines, which a read-to-end-of-sheet ingests as patients
  with a null ID. The manifest bounds the data rows explicitly.
- `Gender` also carries a lowercase `m` and `Grade` an `n/a`, neither of them in
  §Q5's list of inconsistencies.

#### Corrected
- **`docs/data-provenance.md` §Q5 assigned the join assertions to P0-T3, which
  never ran them** — Supplementary Data 1 was not in the repository until P5-T1
  fetched it, so there was nothing to assert against. They are P5-T1's gate.
- **PROJECT_PLAN §6 P5-T3's "~35 patients" is wrong**, and `docs/limitations.md`
  §8 inherits it. Phase 2 scored `TIME-L`/`TIME-B` only, so the patients carrying
  a signature score number **16** — 13 lung + 8 brain, 5 shared. The 35 are the
  patients with a *tumour* AOI and none of them has a score. ADR 0024 §2 records
  why Phase 5 does not chase the larger number.
- PROJECT_PLAN §6 P5-T5 says "5 signatures"; `config/signatures.yaml` resolves
  **six**. The multiplicity family is 6 × 2 = 12.

### Fixed — P6-T1: environments are actually pinned now (ADR 0023) (2026-09-14)
- **The recorded blocker was wrong on both counts.** `conda-lock` was installed
  all along, and Snakemake 8.30 never read `*.conda-lock.yml` — it reads
  `<env>.<platform>.pin.txt` (`deployment/conda.py:128`). **No environment was
  pinned at all**; the three committed lock files were read by nothing.
- Two of the three locks were also stale in content, missing `xarray`
  (py-analysis) and `r-reformulas` (r-geomx). Both had been added to fix real
  failures — a numpy ABI fall-through and the lme4 2.0 / SpatialDecon breakage —
  so installing from those locks would have reproduced two known bugs.
- **Added** `workflow/envs/{py-analysis,r-stats,r-geomx}.osx-arm64.pin.txt`,
  generated with `conda list --explicit --md5` from the environments that
  produced every committed result, so the pin describes what happened rather
  than predicting what would happen.
- **Removed** the three dead `*.conda-lock.yml` files.
- **Corrected** `environment.yml`'s header, which documented the conda-lock
  workflow and is what let the misreading persist.
- Transition needed three steps, recorded in ADR 0023 §3: suppressing the
  software-env trigger alone leaves the repo dirty, because no job runs and
  Snakemake never updates its provenance record. `--touch` is what absorbs it,
  and is truthful only because the pins came from the envs that ran.
- Verified: all three envs load their libraries (including
  `SpatialDecon` + `lme4 2.0.6` + `reformulas 0.4.4`, the combination ADR 0013
  documents as broken), one rule per env forced to re-execute reproduced
  byte-identical output, and **97 of 97 tables and figures are unchanged**.


### Added — Phase 4: inferred crosstalk between adjacent compartments (2026-09-14)
- **GATE 4 PASSED (ADR 0022).** The empirical FDR was computed and honoured and
  nothing survived it: 0 of 268 lung and 0 of 260 brain primary direction-rows
  clear FDR 0.05. Gate 4 licenses this outcome in terms. It is an
  **assay-sensitivity limit, never evidence that the crosstalk is absent**.
- **ADR 0021** pre-registers Phase 4 in nine sections, committed before the
  first rule ran: the duplicate-AOI collapse rule, the CellChatDB v2 pin, the
  0.5 detection floor, the 1:1-primary/complex-exploratory split, keeping every
  paired patient in the primary, two nulls with one multiplicity family per
  site, the nomination cap, the `LGALS9` non-substitution, and testing both
  directions.
- `config/config.yaml` gains a `crosstalk:` block and `phases.crosstalk: true`
  in one edit; 12 negative tests against the new schema block were rejected as
  intended. The rebuild returned **58 of 59 tables byte-identical**, differing
  only in `h5ad_summary.json`'s `git_sha` and `config_sha256`.
- `workflow/rules/06_crosstalk.smk`: 12 rules (`p4t1_build_pairs` …
  `p4t7_language_audit`) and 10 new scripts. **Pure Python plus one 
  narrow R rule** that reads the pinned `.rda` with base `load()` in the
  existing `r-geomx` env — no conda env YAML edited, no dependency added.
- `config/ligand_receptor.yaml` + schema: CellChatDB v2 pinned at commit
  `75253cd`, the project's **fourth sanctioned writer** under ADR 0005, with an
  `expect` block asserted against the parsed database and the P4-T6 comparator
  declared with the source paper's own sentence quoted.
- **The detection gap by interaction class is the phase's larger result** and is
  ADR 0014's claim made numeric on 2,239 interactions: ECM-Receptor 45.3% brain
  / 39.6% lung admitted, Cell-Cell Contact 13.7% / 17.3%, Secreted Signaling
  6.4% / 5.3%. The matrix axis is measurable roughly eightfold more often.
- **Background is a patient-level property** (P4-T3a): across a patient's two
  paired AOIs, `negprobe_log2` ρ +0.548 brain / +0.352 lung and
  `gene_detection_rate` ρ +0.667 / +0.231. Measured rather than inherited from
  P3-T4's negligible *site* gradient, which is a different gradient.
- `p4t7_language_audit` makes PROJECT_PLAN §6's "grep the repo" a **build-time
  assertion**, and classifies a match inside a prohibition separately from a
  claim — the literal instruction fails on CLAUDE.md's own hard constraint 6.
  0 violations, 14 prohibitions across 8 files.
- Figures `crosstalk_null.png` and `crosstalk_network.png`, both wrapped in
  `report()` with captions.

### Fixed — Phase 4 defects caught by checks rather than by eye (2026-09-14)
- **The empirical-FDR step-down was wrong under ties.** At n = 8 Spearman takes
  42 distinct values across 260 rows; per-row ranks gave tied rows different
  denominators. Now computed once per distinct |rho| and mapped back by value.
  Caught by the monotonicity assertion.
- **`CD96–NECTIN1` was entered into the P4-T6 comparator backwards** and
  reported as absent from a database that contains it. The comparator's measure
  is undirected, so the lookup now matches either orientation.
- **The permutation null was not reproducible.** The per-site RNG stream was
  seeded from `hash(site)`, and Python salts string hashing per process — so it
  looked like an explicit seed and was not one. Two runs of the same commit gave
  min empirical FDR 0.228 and 0.225. Now derived from a sha256 of the site name;
  two independent runs verified byte-identical. The conclusion is unchanged
  (0 rows clear FDR 0.05 in either arm); the quoted minima are now 0.499 lung
  and 0.224 brain.
- **Both figures clipped on first render** — the network's right-hand labels
  ("CEACAM" for CEACAM1 and CEACAM5) and the colourbar label. Fixed with the
  `main | colourbar` gridspec. Caught only by opening the PNGs.


### Added — P0-T1: repo scaffold and workflow skeleton (2026-08-23)
- Full repository tree per PROJECT_PLAN §4.2.
- `workflow/Snakefile` with target rules `all`, `qc_only`, `phase2`, `phase3`,
  `fair`; config validated at load against `workflow/schemas/config.schema.yaml`.
- Nine phase rule stubs plus `rules/common.smk`; each phase contributes to a
  target only when both its config switch is on and its `TARGETS_*` list is
  populated.
- `config/config.yaml`, `compartment_map.yaml` (§A.2), and empty
  `signatures.yaml` / `checkpoints.yaml`.
- JSON schemas for config and `samples.tsv`.
- `workflow/envs/py-analysis.yaml` (joins the two validated R env specs).
- `.claude/settings.json` registering two hooks: `no-ipynb` (PreToolUse) and
  `marimo-check` (PostToolUse). Hook commands are **quoted** —
  `"$CLAUDE_PROJECT_DIR/..."` — because this project's absolute path contains
  spaces, and `/bin/sh` word-splits an unquoted expansion.
- `notebooks/README.md` stating the three-tier rules, and a stub explore-tier
  notebook `notebooks/explore/00_first_look.py` for Q1.
- MIT (code) and CC-BY-4.0 (derived data) licences, `CITATION.cff`,
  `codemeta.json`, `CLAUDE.md`.
- ADR 0002 recording the deviation from §5.6's hook registration.

### Added — P0-T1: marimo skills (2026-08-23)
- Installed and pinned five marimo skills (`skills-lock.json` 163 -> 168, all
  K-Dense pins unchanged): `marimo-notebook`, `jupyter-to-marimo`,
  `wasm-compatibility` from `marimo-team/skills`; `marimo-pair`,
  `retro-marimo-pair` from `marimo-team/marimo-pair`.
- Subset rather than full install, per §5.1.2's context-cost and supply-chain
  reasoning: the seven remaining skills in `marimo-team/skills` map to no task
  in this plan.
- Reviewed both bundled `marimo-pair` scripts before use; verified `bash`,
  `curl`, `jq` prerequisites and smoke-tested `discover-servers.sh`.
- Recorded the skill inventory and the pair-mode/constraint-5 interaction in
  `CLAUDE.md`.

### Added — P0-T1: commands, project skill, env locks, security scan (2026-08-23)
- Seven slash commands in `.claude/commands/` per §5.3: `/newrule`, `/dryrun`,
  `/gate`, `/adr`, `/fairscan`, `/marimo-check`, `/newnb`. marimo publishes no
  ready-made `/marimo-check` command in the installed skill trees, so it is
  written here rather than copied.
- Project-local skill `.claude/skills/geomx-dsp-analysis/` (149 lines, §5.4):
  ROI/segment/AOI vocabulary, Q3 rationale, LOQ censoring, within-patient
  correlation, SpatialDecon reference-matrix caveat.
- `.gitignore`: `.claude/skills/` narrowed to `.claude/skills/*` so the
  project-local skill can be re-included. It is source, not an installed
  dependency — nothing in `skills-lock.json` can restore it.
- conda-lock lockfiles for all four env specs, both `osx-arm64` and `linux-64`
  (§4.4): `environment.conda-lock.yml` and `workflow/envs/*.conda-lock.yml`.
  Verified no `zellkonverter` / `basilisk` / `reticulate` entered via a
  transitive dependency (ADR 0001).
- Security-scanned all 169 installed skills (§5.1.2): 148 clean, 21 flagged.
  Findings triaged in ADR 0003 — four of five Phase-0-relevant HIGH/MEDIUM
  findings were false positives on source inspection.

### Earlier
- ADR 0001: AnnData handoff without zellkonverter; validated `r-geomx` and
  `r-stats` env specs on osx-arm64 and linux-64.
