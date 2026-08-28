# ADR 0013 — conda R on a path with spaces, and the limit of the float round-trip

**Date:** 2026-08-28
**Status:** Accepted
**Task:** P2-T3, extended at P2-T5 (infrastructure)
**Amends:** ADR 0001 (the R/Python file boundary)

## Context

Phase 2 contains the first R rule in the project. Two problems surfaced on the
first attempt to run it, both invisible until then, and neither specific to the
script that hit them.

## 1. conda-forge R does not start on a path containing a space

The project lives under `.../Biomedical Data Science/...`. conda relocates R at
install time by rewriting the build-time prefix into `lib/R/bin/R`, a `/bin/sh`
wrapper, and it writes the assignments **unquoted**:

```sh
R_HOME_DIR=/Users/.../Biomedical Data Science/.../lib/R
```

The shell word-splits that, and R dies before it starts with errors that name
neither R nor the path:

```
.../lib/R/bin/R: line 4: Data: command not found
.../lib/R/bin/R: line 250: /etc/ldpaths: No such file or directory
```

`Rscript` does not avoid it — the `Rscript` binary execs this same wrapper. This
affects **every** R rule, present and future, in both `r-stats` and `r-geomx`.

### What was rejected, and why

- **Snakemake's post-deploy hook** (`workflow/envs/<env>.post-deploy.sh`) is the
  designed place for exactly this repair. It cannot be used: Snakemake invokes
  the hook itself unquoted (`{interpreter} {deploy_file}`), so on this path
  env creation fails with
  `bash: /Users/jarrettevans/Documents/Biomedical: No such file or directory`
  and takes the whole environment down with it. **Same bug class, one level up.**
- **Relocating the envs** with `--conda-prefix` to a space-free directory does
  fix R properly, at the root. It was rejected because it changes the recorded
  software environment and fires the provenance rerun trigger on
  `p0t2_fetch_geo`, whose outputs are `protected()` — a dry run confirms the DAG
  then wants to re-download the GEO artifacts and would abort with
  `ProtectedOutputException`. Buying an R fix with a broken acquisition
  provenance chain is a bad trade.

### Decision

`rule p2t0_repair_r_env` runs `workflow/scripts/repair_r_env.sh` inside the
target env, quoting the four relocated assignments, and touches a sentinel that
every R rule declares as an input. It is idempotent, it **verifies that R
actually starts** before touching the sentinel, and it is version-controlled —
so a clean-room rebuild repairs itself rather than requiring a manual edit
inside `.snakemake/conda/`.

The rule is parameterised by env (`r-stats|r-geomx`), so P2-T5's deconvolution
inherits the fix without a second mechanism.

## 1b. The same bug blocks env CREATION, and the fix is a symlink

Repairing R's wrapper after the fact is not sufficient. Building `r-geomx` fails
outright:

```
installBiocDataPackage.sh: line 26: $TARBALL: ambiguous redirect
```

bioconda's post-link script for `bioconductor-genomeinfodbdata` interpolates the
env prefix into a shell redirect unquoted. It runs **during** env creation, long
before any rule of ours could intervene, so `p2t0_repair_r_env` cannot help and
the environment simply does not exist. This is the **third** instance of one root
cause — conda-forge's R wrapper, Snakemake's post-deploy hook, and now bioconda's
post-link script all interpolate the prefix unquoted.

**Decision: build the envs through a space-free symlink to the same directory.**

```bash
ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" ...
```

`~/nsclc-envs` *is* `.snakemake/conda` — same inode, so the envs stay where they
already are, inside the project and gitignored. conda embeds the path string it
was given, which now contains no space, and every post-link script works.

The critical property, and the reason this beats a genuine relocation: because
the symlink resolves to the directory Snakemake already recorded, **it does not
fire the software-environment rerun trigger**. A real `--conda-prefix` move to
`~/.snakemake-conda/...` does — a dry run marks `p0t2_fetch_geo` out of date,
and its `protected()` outputs then abort the DAG with `ProtectedOutputException`
before anything runs. Verified both ways: the symlinked prefix leaves the DAG
untouched; the relocated one does not.

**This flag is now required for every invocation on this machine**, and it is
machine-local rather than in-repo because `$HOME` is not knowable from the
repository. `docs/NEXT_STEPS.md` and `CLAUDE.md` carry it. A clean-room
reproduction (P6-T1) on a path *without* spaces needs neither the symlink nor
`p2t0_repair_r_env` — which is exactly why P6-T1 must be run on a path *with*
one, or it will report a false pass.

## 2. R's parser is not correctly rounded, so the round-trip cannot be exact

ADR 0001 says the R/Python boundary is plain files written with `digits = 17`
and an **asserted exact round-trip**. Asserting that on the Python → R direction
fails, and the reason is a property of R rather than of this pipeline:

| | value |
|---|---|
| Python wrote | `0.44525532065683371` = `0x1.c7f102c25d47cp-2` |
| R parsed | `0x1.c7f102c25d47bp-2` |

**One ULP low.** `as.numeric`, `scan` and `read.table` all go through
`R_strtod` and all three agree on the wrong value; Python's own parser recovers
the correct one. It affected **34 of 276** values on the first run.

Nor is it rescuable by comparing the two 17-digit renderings as strings: `%g`
strips trailing zeros, so a 1-ULP difference can change the string *length*
(`0.09206579933276339` vs `0.092065799332763404`) and a prefix rule then fails
on formatting rather than on precision.

### Decision

**The exactness assertion holds in the direction where it is achievable and is
replaced by a stronger check in the direction where it is not.**

- **R → Python stays exact.** Python's parser is correctly rounded, so
  `assemble_h5ad.py`'s existing assertion is sound and is unchanged. ADR 0001's
  contract is intact where ADR 0001 actually applies.
- **Python → R asserts structure, not bytes.** The R script requires every value
  to parse to a finite double and the frame to carry its expected columns. That
  catches the failures the check exists for — a truncated file, a shifted
  column, a locale decimal comma, text where a double belongs.
- **Numeric agreement across the boundary is verified end-to-end instead**, by
  `rule p2t3_model_crosscheck`: statsmodels refits the Gate 2 signature from the
  same TSV and must agree with lme4 to 1e-4. Observed agreement is ~5e-9 and
  ~2e-8. This is not circular, unlike any comparison R can make against its own
  parse, and it tests what actually matters — that both paths reach the same
  answer.

A 1-ULP difference is ~1e-16 relative, against an `lmer` convergence tolerance
of ~1e-8.

## Consequences

- Every R rule declares `results/interim/env_repair/<env>.ok` as an input.
- **Every invocation on this machine needs
  `--conda-prefix "$HOME/nsclc-envs"`**, with that symlink in place (§1b).
- **A green `--conda-create-envs-only` is not evidence that an env works.** Both
  R envs solved cleanly and neither could run R. Load the libraries. That was
  the second instance in one session — the `P2-T0` commit found `r-geomx`
  solving green with an lme4/Matrix ABI mismatch that appeared only on
  `library(SpatialDecon)` — and §1b was the third.
- **The general lesson: a path containing a space is not a supported
  configuration for the conda/bioconda toolchain.** Three independent scripts in
  three different projects got the same quoting wrong. Expect a fourth; the
  symptom is always a shell error naming a fragment of the path
  ("Biomedical", "Data") or an "ambiguous redirect", never the real cause.
- P6-T1's clean-room reproduction must run on a path with a space, or it will
  not exercise the repair rule and will report a false pass.
