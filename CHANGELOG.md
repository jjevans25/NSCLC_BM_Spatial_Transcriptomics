# Changelog

All notable changes to this project are recorded here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

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
