# Changelog

All notable changes to this project are recorded here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

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
