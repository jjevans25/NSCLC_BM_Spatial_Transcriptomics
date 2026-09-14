# ADR 0023 — Environments are pinned with `.pin.txt` from the envs that ran, not with conda-lock

**Date:** 2026-09-14
**Status:** Accepted
**Task:** P6-T1 (unblocking it)
**Supersedes:** PROJECT_PLAN §4.4's "commit the locks", in mechanism though not
in intent
**Related:** ADR 0005 (pinned digests for data), ADR 0013 (conda on a path with
spaces)

## Context

`docs/NEXT_STEPS.md` had carried this blocker since Phase 3:

> **`py-analysis.conda-lock.yml` is stale** and `conda-lock` is not installed.
> Still the P6-T1 blocker.

**Both halves are wrong, and the truth is worse than either.**

**`conda-lock` is installed** — 4.0.2, in the `nsclc_bm_spatial` env, exactly
where `environment.yml` puts it. That half was never true, or stopped being
true and nobody re-checked.

**Snakemake never read the lock files.** Snakemake 8.30 looks for an auxiliary
file named `<env>.<platform>.pin.txt`
(`snakemake/deployment/conda.py:128`); it has no support for a
`<env>.conda-lock.yml` sitting beside the environment file. No `.pin.txt`
existed anywhere in `workflow/envs/`. **So nothing was pinned at all**, and the
three committed `*.conda-lock.yml` files — 775 KB of them — were read by
nothing. A clean-room P6-T1 would have solved the three YAMLs fresh against
whatever conda-forge served that day.

**And two of the three locks were stale in a way that would have bitten if they
ever had been wired up.** `py-analysis.conda-lock.yml` has no `xarray` and
`r-geomx.conda-lock.yml` has no `r-reformulas`. Both packages were added to
their YAMLs *after* the locks were generated, and both were added to fix a
failure that had already happened:

- `xarray`, because `anndata >= 0.11` does a conditional `import xarray` and
  Snakemake appends its own interpreter's `site-packages` *after* the job env —
  so a missing optional import silently falls through to the user's base
  install and fails as a numpy binary-incompatibility error that names nothing
  relevant.
- `r-reformulas`, because `lme4 >= 1.1_36` split its formula handling into that
  package, and the unconstrained solve produced an `lme4` that neither provides
  nor requires it. `library(SpatialDecon)` died on an ABI mismatch plus
  "there is no package called 'reformulas'".

**Installing from those locks would have reproduced both bugs exactly**, in a
clean-room run whose whole purpose is to prove the pipeline works. A lock file
that is both unread and wrong is worse than no lock file, because it looks like
the problem is solved.

## Decision

### 1. Pin with `<env>.<platform>.pin.txt`, which is what Snakemake actually reads

```
workflow/envs/py-analysis.osx-arm64.pin.txt   245 packages
workflow/envs/r-stats.osx-arm64.pin.txt       170 packages
workflow/envs/r-geomx.osx-arm64.pin.txt       389 packages
```

Generated with `conda list --explicit --md5`, so every line carries a URL and a
checksum. Snakemake logs `Using pinnings from …` and installs from the pin;
if the pin fails it falls back to solving the YAML and says so, which is the
right degradation.

### 2. The pins come from the environments that actually ran, not from a fresh solve

This is the load-bearing half of the decision.

`conda-lock` would have re-solved the YAMLs and produced a lock for *today's*
conda-forge. That is a valid environment; it is **not the environment that
produced the committed results.** Every number in `results/tables/`, every ADR
that quotes one, and the Gate 0–4 record were produced by the three env
directories under `.snakemake/conda`. Pinning from those, with
`conda list --explicit`, makes the pin a **description of what happened** rather
than a prediction about what would happen.

The pins carry, demonstrably, the two packages the stale locks lacked —
`xarray 2026.7.0`, `r-reformulas 0.4.4` — and `r-lme4 2.0_6`, satisfying the
`r-lme4>=2.0` constraint ADR 0013's investigation established.

**Rejected: re-solving with `conda-lock` for `osx-arm64` and `linux-64`.** It
would give Linux reproducibility, which this does not (§4 below). It was
declined because a pin that does not match the run it claims to pin is a
provenance claim nobody checked, and this project has spent four phases
declining exactly that trade — ADR 0017 checked the inputs against the
publication's own file rather than against ourselves, and this is the same
instinct applied to the software.

### 3. The one-time env-hash change is absorbed with `--touch`, not by suppressing the trigger

`content_pin` feeds the environment hash (`conda.py:257`), so adding pin files
changes the hash of all three envs, which fires the software-env rerun trigger
on **every** rule — including `p0t2_fetch_geo`, whose outputs are `protected()`
and which would abort the DAG with `ProtectedOutputException`. Dry run before
the transition: **57 jobs, every one of them "Software environment definition
has changed".**

**Suppressing the trigger is not sufficient on its own, and finding that out
is part of this record.** The obvious procedure —

```
snakemake ... --rerun-triggers mtime params input code
```

— reports "Nothing to be done" and leaves the repository **dirty**: because no
job runs, Snakemake never updates its recorded provenance, so a plain dry run
still asks for all 55 jobs. Skipping a trigger is not the same as absorbing it.

The transition that works is three steps:

```
# 1. build the envs from the pins
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" --conda-create-envs-only

# 2. re-repair the NEW R envs for the path with spaces (ADR 0013). Scoped with
#    -f to the two sentinels: a bare --forcerun pulls in the whole DAG and
#    aborts on p3t2b_fetch_supplementary's protected outputs.
snakemake ... --rerun-triggers mtime params input code \
    -f results/interim/env_repair/r-stats.ok results/interim/env_repair/r-geomx.ok

# 3. refresh the provenance record without recomputing
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" --cores 4 --touch
```

`--touch` is the honest instrument here **precisely because the pins were taken
from the envs that produced these outputs**: the software did not change, only
Snakemake's record of which directory it lived in. Touching asserts "these
outputs are current for this environment", and that assertion is true by
construction. It would be a lie under any other pinning source, which is a second
reason §2 chose this one.

### 4. The pinned envs are verified by running the pipeline in them, not by building them

ADR 0013's standing lesson is that a green `--conda-create-envs-only` is not
evidence an env works. So all three were loaded —

- `r-geomx`: `SpatialDecon`, `GeomxTools` and `lme4 2.0.6` with
  `reformulas 0.4.4`, which is the exact combination ADR 0013 documents as
  failing when solved unconstrained;
- `r-stats`: `lmerTest` + `lme4 2.0.6`;
- `py-analysis`: `anndata 0.13.2`, `xarray 2026.7.0`, `numpy 2.5.2`,
  `networkx 3.6.1`

— and then one rule per environment was **forced to re-execute**:
`p4t7_language_audit` (py-analysis), `p4t2c_export_lr_database` (r-geomx) and
`p3t3_carrier_models` (r-stats). **All three produced byte-identical output**,
including the mixed-model fit, which is the most numerically sensitive artifact
in the project. Across the whole transition, **97 of 97 tables and figures are
byte-identical** to the pre-pin snapshot.

**Rejected: unprotecting `resources/` and re-running everything including the
GEO downloads.** A cleaner provenance chain in principle, but it re-downloads
inputs whose digests are already pinned and verified under ADR 0005, to
re-derive outputs from software that is by construction identical. It buys a
stronger-sounding claim and no more actual evidence.

**Rejected: deferring the pins to P6-T1 itself.** P6-T1 is a fresh clone with
fresh envs, so it would be the first real use — but this repository's dry run
would stay permanently dirty until then, and "the dry run must stay clean" is
the standard every gate in this project has been checked against.

## Consequences

- **P6-T1 is unblocked in the sense that matters**: a clean-room run on
  `osx-arm64` now installs the exact packages that produced the results, rather
  than solving fresh.
- **This does not make the project reproducible on Linux.** The pins are
  `osx-arm64` only, because `conda list --explicit` can only describe an
  environment that exists and these three only exist for this platform. On
  `linux-64` Snakemake finds no pin file, says so, and solves the YAML — which
  is what it did on every platform until now, so nothing is lost. **P6-T6's
  limitations document must state this**, and a multi-platform lock rendered to
  `linux-64.pin.txt` is the obvious Phase 6 follow-up if cross-platform
  reproduction is wanted.
- **The three `*.conda-lock.yml` files are removed.** They were read by nothing,
  two of three described environments this project would actively reject, and
  their presence is what produced a wrong blocker note that survived two phases
  and two gates. A file that looks authoritative, is never read, and is wrong is
  worse than an absent one — the next reader sees "conda-lock.yml" and concludes
  pinning is handled. They are recoverable from git history if a multi-platform
  lock is ever wanted; `conda-lock` stays in `environment.yml` for that.
- **`environment.yml`'s header comment still documents the `conda-lock`
  workflow** and should be corrected to describe the `.pin.txt` mechanism, or
  the same misreading will recur.
- Regenerating a pin after any env YAML edit is now a required step, not an
  optional one, because a stale pin installs the wrong environment silently
  rather than failing.
