# Session notes — environment and tooling gotchas

Things that cost time once and would cost it again. Not provenance (that is
`docs/data-provenance.md`), not decisions (those are `docs/decisions/`) — just
state about this working environment that is not derivable from the code.

Carried forward from `docs/NEXT_STEPS.md`, which was retired on 2026-08-24 once
its three action items were done or superseded.

## Deliberate-looking-wrong things in the repo

- **`config.yaml` QC thresholds are `null` on purpose.** The schema permits
  `null` so that an unmade decision cannot masquerade as a made one. They get
  real values at the P0-T5 gate, with an ADR citing
  `notebooks/review/qc_review.py`. Do not "fix" them by inventing a number.
- **`phases.*` switches start `false`.** Each is flipped on in the *same commit*
  that adds its rules and populates the matching `TARGETS_*` list. Either alone
  contributes nothing — `targets()` in `workflow/rules/common.smk` requires both.
- **Env locks are provisional.** `environment.yml`'s header intends locking at
  P6-T1, but `workflow/envs/*.conda-lock.yml` were generated early. Regenerate
  them if any `workflow/envs/*.yaml` changes; treat them as unpinned until P6-T1.

## Tooling

- **`gh` is not installed** (confirmed 2026-08-24). PRs must be opened in the
  browser, or `brew install gh`.
- **`marimo-pair` needs a browser-attached kernel.** `marimo edit --headless`
  starts a server but creates no session, so `execute-code.sh` fails with
  "No active sessions". Open the notebook in a browser first.
- **The `no-ipynb` hook matches literal command text.** A `Bash` command that
  merely *quotes* an `.ipynb` creation string gets blocked too. Harmless, but
  surprising the first time. See ADR 0002 for why it is `PreToolUse` on `Bash`.
- **System `python3` has no `pyyaml`.** For quick YAML checks outside a rule,
  `uvx --with pyyaml python -c ...` works without touching the project envs.

## Snakemake

- **Always pass `--use-conda`.** Every rule declares `conda:`, so without the
  flag the recorded software stack does not match and the `software-env` rerun
  trigger fires. For P0-T2 that means Snakemake tries to re-download artifacts
  whose outputs are `protected()`, and the DAG build dies with
  `ProtectedOutputException`. `snakemake -n` alone is not a clean dry run.
- **To re-fetch a protected artifact** you must defeat the protection on
  purpose: `chmod u+w resources/raw/<file> && rm resources/raw/<file>`. Editing
  `workflow/scripts/acquire_geo.py` also marks those jobs out of date (the
  `code` trigger) and needs the same recovery.
- **`--list-params-changes` over-reports.** It will name a file immediately
  after a clean run with no edits, while `snakemake -n` correctly reports
  nothing to do. Do not use it to diagnose a rerun; read the `reason:` line in
  the dry-run output instead.

## Skills

**15 flagged skills are still uncleared** (ADR 0003). `literature-review` (P4),
`citation-management` and `scientific-schematics` (P6) need the same
read-the-source triage at the phase that first uses them. The other ~12 map to
no task in this plan and are removal candidates.

`paper-lookup` was cleared for use in ADR 0003 and is the skill behind the Q2/Q5
answers in `docs/data-provenance.md`.
