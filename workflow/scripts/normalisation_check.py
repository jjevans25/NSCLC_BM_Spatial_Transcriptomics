"""P0-T4 — settle Q1 empirically: what state is the GEO matrix in?

Owner task: P0-T4. Driven by rule p0t4_normalisation_check.

`docs/data-provenance.md` already carries a *provisional* answer from the GEO
metadata (every GSM says counts "were scaled to the 75th percentile"). This
rule is the check against the file itself, which is what §3 actually asks for.

The four discriminating tests, from PROJECT_PLAN §3:

    | Check                       | Raw    | Q3-normalised   | Log |
    |-----------------------------|--------|-----------------|-----|
    | all values integer          | yes    | no              | no  |
    | value range                 | 0-10^4 | 0-10^4 rescaled | ~20 |
    | column-wise Q3 across AOIs  | wide   | ~constant       | wide|
    | minimum non-zero value      | 1      | fractional      | frac|

The third is the decisive one: Q3 normalisation is precisely the operation that
makes each AOI's third quartile equal, so a near-zero spread in column-wise Q3
is the signature, and a wide spread rules it out.

`NegProbe-WTX` is excluded from the gene statistics — it is a control, not a
gene, and leaving it in would contaminate the detection and library-size
figures that P0-T5 builds on.
"""

import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NEG_PROBE = "NegProbe-WTX"
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    # ------------------------------------------------------------------ load
    genes, values = [], []
    with gzip.open(snakemake.input.matrix, "rt", encoding="utf-8") as handle:
        aois = handle.readline().rstrip("\n").split("\t")[1:]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[0])
            values.append(fields[1:])

    full = pd.DataFrame(np.asarray(values, dtype=float), index=genes, columns=aois)
    emit(f"matrix: {full.shape[0]} rows x {full.shape[1]} AOIs")

    neg = full.loc[NEG_PROBE] if NEG_PROBE in full.index else None
    expr = full.drop(index=NEG_PROBE) if neg is not None else full
    emit(f"genes (excluding {NEG_PROBE}): {expr.shape[0]}")

    flat = expr.to_numpy().ravel()

    # ----------------------------------------------------------- the 4 tests
    n_values = flat.size
    n_integer = int(np.sum(flat == np.floor(flat)))
    frac_integer = n_integer / n_values
    n_negative = int(np.sum(flat < 0))
    n_zero = int(np.sum(flat == 0))
    nonzero = flat[flat > 0]
    min_nonzero = float(nonzero.min())
    vmax = float(flat.max())

    col_q3 = expr.quantile(0.75)
    q3_cv = float(col_q3.std(ddof=1) / col_q3.mean())
    col_sum = expr.sum()
    sum_cv = float(col_sum.std(ddof=1) / col_sum.mean())

    emit("\n--- the four §3 tests ---")
    emit(f"  all values integer      : {frac_integer:.4%} integer  -> {'yes' if frac_integer > 0.999 else 'NO'}")
    emit(f"  value range             : {flat.min():.4g} .. {vmax:.4g}")
    emit(f"  column-wise Q3 spread   : mean {col_q3.mean():.4f}, sd {col_q3.std(ddof=1):.4g}, CV {q3_cv:.3%}")
    emit(f"  minimum non-zero value  : {min_nonzero:.6g}")
    emit(f"  (library-size CV, for contrast: {sum_cv:.2%})")
    emit(f"  negative values         : {n_negative}")
    emit(f"  exact zeros             : {n_zero}")

    # ------------------------------------------------------------- the verdict
    is_integer = frac_integer > 0.999
    q3_constant = q3_cv < 0.01          # Q3 scaling makes this ~0 by construction
    log_like = vmax < 100               # log2 of a 10^4-scale count tops out ~14

    if is_integer:
        verdict, reason = "raw_counts", "values are integer"
    elif log_like:
        verdict, reason = "log_transformed", f"maximum value {vmax:.3g} is on a log scale"
    elif q3_constant:
        verdict, reason = (
            "q3_normalised",
            f"column-wise Q3 is constant to within CV {q3_cv:.3%} while library "
            f"size varies by CV {sum_cv:.1%} — the signature of per-AOI Q3 scaling",
        )
    else:
        verdict, reason = "unknown", "no test was decisive — inspect by hand before proceeding"

    expected = snakemake.params.expected_state
    agrees = verdict == expected

    emit(f"\nVERDICT: {verdict}\n  because {reason}")
    emit(f"  GEO metadata predicted: {expected} -> {'AGREES' if agrees else 'DISAGREES'}")

    # ------------------------------------------------------------ per-AOI table
    per_aoi = pd.DataFrame(
        {
            "aoi_label": expr.columns,
            "library_size": col_sum.to_numpy().round(4),
            "q3": col_q3.to_numpy().round(6),
            "median": expr.median().to_numpy().round(4),
            "max": expr.max().to_numpy().round(4),
            "n_zero_genes": (expr == 0).sum().to_numpy(),
        }
    )
    if neg is not None:
        per_aoi["negprobe"] = neg.to_numpy().round(4)
    per_aoi.to_csv(snakemake.output.per_aoi, sep="\t", index=False)
    emit(f"\nwrote {snakemake.output.per_aoi}")

    Path(snakemake.output.summary).write_text(
        json.dumps(
            {
                "n_genes": int(expr.shape[0]),
                "n_aoi": int(expr.shape[1]),
                "negprobe_present": neg is not None,
                "fraction_integer": round(frac_integer, 6),
                "n_negative_values": n_negative,
                "n_exact_zeros": n_zero,
                "min_nonzero": min_nonzero,
                "max": vmax,
                "col_q3_mean": round(float(col_q3.mean()), 6),
                "col_q3_cv": round(q3_cv, 8),
                "library_size_cv": round(sum_cv, 6),
                "verdict": verdict,
                "reason": reason,
                "geo_metadata_expected": expected,
                "agrees_with_metadata": agrees,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    emit(f"wrote {snakemake.output.summary}")

    # ----------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

    ax = axes[0]
    ax.plot(np.sort(col_q3.to_numpy()), marker="o", ms=3, lw=1.0, color="#3E7CB1")
    ax.set_title(f"Column-wise Q3 across AOIs\nCV = {q3_cv:.3%} — the decisive test", fontsize=10)
    ax.set_xlabel("AOI (sorted)", fontsize=9)
    ax.set_ylabel("third quartile", fontsize=9)
    span = max(col_q3.max() - col_q3.min(), 1e-9)
    ax.set_ylim(col_q3.min() - 4 * span, col_q3.max() + 4 * span)
    ax.grid(alpha=0.25, lw=0.6)

    ax = axes[1]
    ax.plot(np.sort(col_sum.to_numpy()), marker="o", ms=3, lw=1.0, color="#B0763F")
    ax.set_title(f"Library size across AOIs\nCV = {sum_cv:.1%} — varies, as it should", fontsize=10)
    ax.set_xlabel("AOI (sorted)", fontsize=9)
    ax.set_ylabel("sum of gene values", fontsize=9)
    ax.grid(alpha=0.25, lw=0.6)

    ax = axes[2]
    sub = flat[flat > 0]
    ax.hist(np.log10(sub), bins=120, color="#5B8C5A", alpha=0.85)
    ax.set_title(f"Value distribution (non-zero)\nmax {vmax:.4g}, min non-zero {min_nonzero:.3g}", fontsize=10)
    ax.set_xlabel("log10(value)", fontsize=9)
    ax.set_ylabel("count", fontsize=9)
    ax.grid(alpha=0.25, lw=0.6)

    fig.suptitle(
        f"P0-T4 / Q1 — GEO matrix normalisation state: {verdict.replace('_', ' ').upper()}"
        f"   ({expr.shape[0]} genes x {expr.shape[1]} AOIs)",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(snakemake.output.figure, dpi=200, bbox_inches="tight")
    plt.close(fig)
    emit(f"wrote {snakemake.output.figure}")

    if not agrees:
        emit(
            "\nWARNING: the file disagrees with what the GEO metadata claimed. "
            "The file wins. Update docs/data-provenance.md Q1 and re-check "
            "config normalisation.method before P0-T5."
        )
