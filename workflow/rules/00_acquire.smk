# ---------------------------------------------------------------------------
# Phase 0 — acquisition: GEO download, checksums, idempotent re-run.
#
# Owner task: P0-T2 (see Markdowns/PROJECT_PLAN.md §6).
#
# resources/ is read-only (CLAUDE.md hard constraint 1). This is the one rule
# that writes into it, and the constraint is then enforced by the filesystem
# rather than by comment: artifact outputs are protected(), so Snakemake chmods
# them 0444 on success and every later rule is mechanically unable to modify
# them. ADR 0005.
#
# TO RE-FETCH an artifact you must therefore defeat that on purpose:
#     chmod u+w resources/raw/<file> && rm resources/raw/<file>
# That friction is the point. It also means editing acquire_geo.py marks these
# jobs out of date (--rerun-triggers code) and Snakemake will refuse to
# overwrite; same recovery.
#
# AND: always invoke with --use-conda. Without it the software-env rerun trigger
# fires against the recorded conda stack, these jobs are marked out of date, and
# the DAG build dies on the protected outputs before doing anything. A bare
# `snakemake -n` is not a clean dry run once this phase is on.
#
# Note what is NOT a param of p0t2_fetch_geo: the pinned SHA-256. Verification
# lives in p0t2_record_provenance so that adding a pin re-runs only that cheap
# rule, instead of invalidating a protected download. See acquire_geo.py.
# ---------------------------------------------------------------------------

import re

ACQUIRE = config["acquire"]
ARTIFACTS = ACQUIRE["artifacts"]

# dest filename -> manifest key. The wildcard is the filename, because that is
# what the output path actually is; the key is carried through for the record.
KEY_BY_DEST = {spec["dest"]: key for key, spec in ARTIFACTS.items()}

if len(KEY_BY_DEST) != len(ARTIFACTS):
    raise WorkflowError(
        "config acquire.artifacts: two artifacts share a `dest` filename. "
        "Each must land at its own path under paths.raw."
    )

RAW_DIR = PATHS["raw"]
ACQUIRE_META = f"{PATHS['interim']}/acquire"


rule p0t2_fetch_geo:
    """Download one GEO artifact and record what arrived."""
    wildcard_constraints:
        # Only the four manifested filenames. Without this the {artifact}
        # wildcard would offer to produce any path under resources/raw/, and a
        # later rule wanting a file there would collide with this one.
        artifact="|".join(re.escape(dest) for dest in sorted(KEY_BY_DEST)),
    output:
        artifact=protected(f"{RAW_DIR}/{{artifact}}"),
        meta=f"{ACQUIRE_META}/{{artifact}}.json",
    params:
        url=lambda w: f"{ACQUIRE['base_url']}/{ARTIFACTS[KEY_BY_DEST[w.artifact]]['remote']}",
        key=lambda w: KEY_BY_DEST[w.artifact],
        stability=lambda w: ARTIFACTS[KEY_BY_DEST[w.artifact]]["stability"],
    log:
        f"{PATHS['logs']}/p0t2_fetch_geo_{{artifact}}.log",
    benchmark:
        f"{PATHS['benchmarks']}/p0t2_fetch_geo_{{artifact}}.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/acquire_geo.py"


rule p0t2_record_provenance:
    """Verify the acquired artifacts and write the provenance records."""
    input:
        artifacts=[f"{RAW_DIR}/{dest}" for dest in sorted(KEY_BY_DEST)],
        meta=[f"{ACQUIRE_META}/{dest}.json" for dest in sorted(KEY_BY_DEST)],
    output:
        checksums=f"{PATHS['resources']}/checksums.sha256",
        provenance=f"{PATHS['resources']}/provenance.tsv",
    params:
        artifacts=ARTIFACTS,
    log:
        f"{PATHS['logs']}/p0t2_record_provenance.log",
    benchmark:
        f"{PATHS['benchmarks']}/p0t2_record_provenance.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/record_provenance.py"


# The phase contributes to `all` only because config phases.acquire is true AND
# this list is populated — both, by design (targets() in common.smk).
TARGETS_ACQUIRE = [
    f"{PATHS['resources']}/checksums.sha256",
    f"{PATHS['resources']}/provenance.tsv",
]
