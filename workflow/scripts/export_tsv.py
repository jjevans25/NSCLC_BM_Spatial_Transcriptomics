"""P0-T7a — write the AnnData handoff as plain TSVs.

Owner task: P0-T7. Driven by rule p0t7_export_tsv.

ADR 0001 puts the R/Python boundary at plain files rather than zellkonverter:
one side writes `expr_normalised.tsv` + `obs.tsv` + `var.tsv` + `uns.json`, the
other assembles the `.h5ad`. ADR 0006 then moved P0-T5 to Python, so today both
sides are Python and the file boundary buys nothing at runtime.

It is kept anyway, for two reasons. The interface is the part ADR 0001 actually
decided, and it is what a future R rule would write into unchanged if a PKC
ever makes `GeomxTools` usable. And the round-trip assertion in the assembler
is a real check on `%.17g` formatting that would otherwise go untested until
the day it matters.

Floats are written with 17 significant digits, which is what IEEE-754 double
round-tripping requires. The assembler asserts exact equality on read.

X is log2(Q3 + 1) — Q1 established the matrix arrives Q3-normalised, so this
step is the log, not the normalisation. The untransformed Q3 values are carried
alongside so nothing is lost.
"""

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

NEG_PROBE = "NegProbe-WTX"
FMT = "%.17g"

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    log_base = snakemake.params.log_base
    pseudocount = snakemake.params.pseudocount
    method = snakemake.params.norm_method
    if method != "verify":
        raise RuntimeError(
            f"normalisation.method is {method!r}; P0-T7 expects 'verify' because "
            "the GEO matrix arrives already Q3-normalised (Q1). Re-check P0-T4."
        )
    if log_base != 2:
        raise RuntimeError(f"only log2 is implemented; config asks for base {log_base}")

    # ------------------------------------------------------------------ load
    genes, values = [], []
    with gzip.open(snakemake.input.matrix, "rt", encoding="utf-8") as handle:
        aois = handle.readline().rstrip("\n").split("\t")[1:]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[0])
            values.append(fields[1:])
    full = pd.DataFrame(np.asarray(values, dtype=np.float64), index=genes, columns=aois)

    # The negative probe is a control, not a gene. It leaves the gene matrix and
    # survives as a per-AOI obs column, where it belongs.
    negprobe = full.loc[NEG_PROBE]
    q3 = full.drop(index=NEG_PROBE)
    emit(f"genes: {q3.shape[0]}  AOIs: {q3.shape[1]}  ({NEG_PROBE} moved to obs)")

    qc = pd.read_csv(snakemake.input.qc_metrics, sep="\t").set_index("aoi_label")
    if set(qc.index) != set(q3.columns):
        raise RuntimeError("qc_metrics.tsv and the matrix disagree on AOI labels")

    # obs must carry the whole samples.tsv contract (§6: "obs columns match
    # samples.tsv"), not just the subset qc_metrics.tsv happens to echo. Start
    # from samples.tsv and add the QC columns it does not already have.
    samples = pd.read_csv(snakemake.input.samples, sep="\t", dtype=str).set_index("aoi_label")
    if set(samples.index) != set(q3.columns):
        raise RuntimeError("samples.tsv and the matrix disagree on AOI labels")
    extra = [c for c in qc.columns if c not in samples.columns]

    # AOIs (obs) as rows, genes (var) as columns — AnnData's orientation.
    obs = samples.loc[q3.columns].join(qc[extra])
    obs["negprobe"] = negprobe.loc[q3.columns]
    obs.index.name = "aoi_label"

    x_q3 = q3.T.to_numpy(dtype=np.float64)
    x_log = np.log2(x_q3 + pseudocount)

    var = pd.DataFrame(index=pd.Index(q3.index, name="gene"))
    var["mean_log2"] = x_log.mean(axis=0)
    var["detected_in_n_aoi"] = (
        x_q3 > np.asarray(obs["negprobe"])[:, None] * snakemake.params.background_multiple
    ).sum(axis=0)

    # ----------------------------------------------------------------- write
    np.savetxt(snakemake.output.expr, x_log, delimiter="\t", fmt=FMT)
    emit(f"wrote {snakemake.output.expr}  shape {x_log.shape}  fmt {FMT}")

    np.savetxt(snakemake.output.expr_q3, x_q3, delimiter="\t", fmt=FMT)
    emit(f"wrote {snakemake.output.expr_q3}")

    obs.to_csv(snakemake.output.obs, sep="\t", float_format=FMT)
    var.to_csv(snakemake.output.var, sep="\t", float_format=FMT)
    emit(f"wrote {snakemake.output.obs}  ({obs.shape[1]} columns)")
    emit(f"wrote {snakemake.output.var}  ({var.shape[1]} columns)")

    Path(snakemake.output.uns).write_text(
        json.dumps(
            {
                "layer_order": ["X = log2(q3 + %g)" % pseudocount, "layers['q3']"],
                "normalisation": {
                    "method": method,
                    "applied_by": "GEO submitters (Q3, verified at P0-T4)",
                    "log_base": log_base,
                    "pseudocount": pseudocount,
                },
                "negprobe_in_obs": True,
                "detection_background_multiple": snakemake.params.background_multiple,
                "float_format": FMT,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    emit(f"wrote {snakemake.output.uns}")
