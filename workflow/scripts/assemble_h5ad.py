"""P0-T7b — assemble the .h5ad from the TSV handoff, asserting round-trip.

Owner task: P0-T7. Driven by rule p0t7_assemble_h5ad. See ADR 0001.

The acceptance criterion from PROJECT_PLAN §6 is that the float round-trip is
*exact*: the TSVs are written with 17 significant digits, and this side asserts
bit-for-bit equality after re-reading rather than an approximate closeness. A
`np.allclose` here would pass on a lossy writer and hide the bug the digits
setting exists to prevent.
"""

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)


def git_sha():
    """Best-effort commit id. Computed here rather than passed as a param so a
    new commit does not mark every downstream job out of date."""
    try:
        sha = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10
        )
        if sha.returncode != 0:
            return "unknown"
        dirty = subprocess.run(
            ["git", "status", "--porcelain"], capture_output=True, text=True, timeout=10
        )
        return sha.stdout.strip() + ("-dirty" if dirty.stdout.strip() else "")
    except (OSError, subprocess.SubprocessError):
        return "unknown"


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    x_log = np.loadtxt(snakemake.input.expr, delimiter="\t", dtype=np.float64)
    x_q3 = np.loadtxt(snakemake.input.expr_q3, delimiter="\t", dtype=np.float64)
    obs = pd.read_csv(snakemake.input.obs, sep="\t", index_col=0)
    var = pd.read_csv(snakemake.input.var, sep="\t", index_col=0)
    uns = json.loads(Path(snakemake.input.uns).read_text(encoding="utf-8"))

    if x_log.shape != (obs.shape[0], var.shape[0]):
        raise RuntimeError(
            f"shape mismatch: X {x_log.shape} vs obs {obs.shape[0]} x var {var.shape[0]}"
        )

    # --- the §6 acceptance criterion: exact, not approximate -----------------
    pseudocount = uns["normalisation"]["pseudocount"]
    recomputed = np.log2(x_q3 + pseudocount)
    exact = bool(np.array_equal(x_log, recomputed))
    max_ulp_gap = float(np.max(np.abs(x_log - recomputed))) if not exact else 0.0
    emit(f"round-trip log2(q3+{pseudocount:g}) reproduces X exactly: {exact}")
    if not exact:
        raise RuntimeError(
            "float round-trip is not exact — X does not reproduce from the q3 layer "
            f"(max absolute gap {max_ulp_gap:.3g}). The TSV writer is lossy; check "
            "the %.17g format in export_tsv.py before trusting anything downstream."
        )

    # obs must still be the samples.tsv contract (§6).
    samples = pd.read_csv(snakemake.input.samples, sep="\t", dtype=str)
    missing = sorted(set(samples.columns) - set(obs.columns) - {"aoi_label"})
    if missing:
        raise RuntimeError(f"obs is missing samples.tsv column(s): {missing}")
    if sorted(obs.index) != sorted(samples["aoi_label"]):
        raise RuntimeError("obs index does not match samples.tsv aoi_label")
    emit(f"obs carries all {len(samples.columns)} samples.tsv columns, {obs.shape[1]} total")

    config_hash = hashlib.sha256(
        Path(snakemake.input.config).read_bytes()
    ).hexdigest()

    adata = ad.AnnData(
        X=x_log,
        obs=obs,
        var=var,
        layers={"q3": x_q3},
        uns={
            **uns,
            "project": snakemake.params.project_code,
            "geo_accession": snakemake.params.geo_accession,
            "source_pmid": snakemake.params.source_pmid,
            "git_sha": git_sha(),
            "config_sha256": config_hash,
            "run_date_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            "round_trip_exact": exact,
        },
    )
    adata.write_h5ad(snakemake.output.h5ad, compression="gzip")
    emit(f"\nwrote {snakemake.output.h5ad}")
    emit(f"  X       {adata.X.shape}  log2(q3+{pseudocount:g})")
    # anndata 0.13 lists a None key in layers as an alias for X; drop it.
    emit(f"  layers  {[k for k in adata.layers if k is not None]}")
    emit(f"  obs     {adata.obs.shape[1]} columns")
    emit(f"  var     {adata.var.shape[1]} columns")
    emit(f"  git     {adata.uns['git_sha']}")
    emit(f"  config  sha256 {config_hash[:16]}...")

    # Prove it loads, from disk, as a separate object (§6: "loads in Python").
    check = ad.read_h5ad(snakemake.output.h5ad)
    if not np.array_equal(check.X, x_log):
        raise RuntimeError("h5ad X does not match what was written")
    if not np.array_equal(check.layers["q3"], x_q3):
        raise RuntimeError("h5ad q3 layer does not match what was written")
    emit("\nre-read from disk: X and layers['q3'] both bit-identical")

    Path(snakemake.output.summary).write_text(
        json.dumps(
            {
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "layers": [k for k in adata.layers if k is not None],
                "obs_columns": list(adata.obs.columns),
                "round_trip_exact": exact,
                "git_sha": adata.uns["git_sha"],
                "config_sha256": config_hash,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
