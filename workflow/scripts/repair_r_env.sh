#!/usr/bin/env bash
# Repair conda-forge R inside an activated conda env whose prefix contains a space.
#
# conda relocates R at install time by rewriting the build-time prefix into
# lib/R/bin/R, a /bin/sh wrapper, and it writes the assignments UNQUOTED:
#
#     R_HOME_DIR=/Users/.../Biomedical Data Science/.../lib/R
#
# On a prefix containing a space the shell word-splits that and R dies before it
# starts, with errors naming neither R nor the path:
#
#     .../lib/R/bin/R: line 4: Data: command not found
#     .../lib/R/bin/R: line 250: /etc/ldpaths: No such file or directory
#
# This project lives under "Biomedical Data Science", so EVERY R rule hits it.
# Rscript does not sidestep it — the Rscript binary execs this same wrapper.
#
# Why this is a rule and not a post-deploy script or a --conda-prefix change:
#   * Snakemake's post-deploy hook invokes the script UNQUOTED itself
#     (`{interpreter} {deploy_file}`), so on this path the hook fails with
#     "bash: /Users/jarrettevans/Documents/Biomedical: No such file or
#     directory" and takes env creation down with it. Same bug, one level up.
#   * Relocating the envs with --conda-prefix to a space-free path does work,
#     but it changes the recorded software environment and fires the provenance
#     rerun trigger on p0t2_fetch_geo, whose outputs are protected() — the DAG
#     then aborts with ProtectedOutputException before doing anything.
# See ADR 0013.
#
# Idempotent: quoting an already-quoted assignment is a no-op.
set -euo pipefail

wrapper="${CONDA_PREFIX:?repair_r_env: CONDA_PREFIX is not set — run me from a conda-enabled rule}/lib/R/bin/R"

if [ ! -f "$wrapper" ]; then
    echo "repair_r_env: no R wrapper at $wrapper — nothing to repair" >&2
    exit 1
fi

perl -pi -e 's{^(R_HOME_DIR|R_SHARE_DIR|R_INCLUDE_DIR|R_DOC_DIR)=(?!")(.*)$}{$1="$2"}' "$wrapper"

# Prove it starts, so a still-broken wrapper fails here rather than inside the
# first analysis rule, where the error names neither R nor the path.
"${CONDA_PREFIX}/bin/Rscript" --vanilla -e 'cat("repair_r_env: R starts —", R.version.string, "\n")' >&2

echo "repair_r_env: quoted the relocated path assignments in $wrapper" >&2
