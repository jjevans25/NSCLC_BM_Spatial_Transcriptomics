"""P0-T5 — per-AOI QC metrics, flags, and normalisation verification.

Owner task: P0-T5. Driven by rule p0t5_qc_metrics.

Python rather than R, deliberately: PROJECT_PLAN §6 names GeomxTools/standR,
but both need a PKC to map the DCC's RTS probe ids to genes and GEO carries
none. No GeoMxSet can be built, so there is nothing here to reimplement.
See ADR 0006. R re-enters at Phase 2 for the lme4 mixed models.

Two sources, because neither alone is enough:

  DCC headers (inside GSE200563_RAW.tar) give the *sequencing* metrics — raw,
  trimmed, stitched and aligned reads, Q30 rates, and the deduplicated count
  implied by the Code_Summary block. These need no PKC. This is where real
  library size lives.

  The processed matrix gives *detection*, relative to NegProbe-WTX. It cannot
  give library size: it is already Q3-normalised (Q1), so its column sums are
  near-constant by construction and measure nothing.

Thresholds are read from config and may be null. Null means the decision has
not been made yet: metrics are still computed and written, nothing is flagged,
and the log says so. The value arrives from notebooks/review/qc_review.py with
an ADR, per §6.

Policy is flag, never drop (config qc.flag_only). Flagged AOIs are listed in
qc_excluded.tsv with a reason and stay in the dataset.
"""

import gzip
import io
import json
import re
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NEG_PROBE = "NegProbe-WTX"
CODE_ORDER = ["L", "LB", "mLN", "TIME-L", "TIME-B", "TBME", "BC"]
COMPARTMENT_COLOUR = {
    "tumour": "#B0763F", "immune": "#3E7CB1",
    "glial_stroma": "#5B8C5A", "normal_control": "#8C8C8C",
}

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)


def parse_dcc(text):
    """Pull the header metrics and the deduplicated count out of one DCC."""
    out = {}
    for key in ("Raw", "Trimmed", "Stitched", "Aligned"):
        m = re.search(rf"^{key},(\d+)\s*$", text, re.M)
        out[key.lower()] = int(m.group(1)) if m else np.nan
    for key in ("umiQ30", "rtsQ30"):
        m = re.search(rf"^{key},([\d.]+)\s*$", text, re.M)
        out[key.lower()] = float(m.group(1)) if m else np.nan
    m = re.search(r"^ID,(\S+)\s*$", text, re.M)
    out["scan_id"] = m.group(1) if m else None

    block = re.search(r"<Code_Summary>(.*?)</Code_Summary>", text, re.S)
    counts = [
        int(line.split(",")[1])
        for line in block.group(1).strip().splitlines()
        if "," in line
    ] if block else []
    out["deduplicated"] = int(sum(counts))
    out["n_probes"] = len(counts)
    return out


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    qc_cfg = snakemake.params.qc
    thresholds = {
        "min_gene_detection_rate": qc_cfg["min_gene_detection_rate"],
        "min_library_size": qc_cfg["min_library_size"],
        "min_sequencing_saturation": qc_cfg["min_sequencing_saturation"],
    }
    background_multiple = qc_cfg["detection_background_multiple"]
    flag_only = qc_cfg["flag_only"]

    samples = pd.read_csv(snakemake.input.samples, sep="\t", dtype=str)

    # ------------------------------------------------------- DCC header metrics
    # Read straight out of the tar: no intermediate extraction, and resources/
    # stays untouched (hard constraint 1).
    dcc_rows = {}
    with tarfile.open(snakemake.input.raw_tar) as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".dcc.gz"):
                continue
            gsm = member.name.split("_", 1)[0].replace("./", "")
            with gzip.open(io.BytesIO(tar.extractfile(member).read()), "rt") as handle:
                dcc_rows[gsm] = parse_dcc(handle.read())
    emit(f"parsed {len(dcc_rows)} DCC headers from {snakemake.input.raw_tar}")

    dcc = pd.DataFrame.from_dict(dcc_rows, orient="index")
    missing = sorted(set(samples["gsm_id"]) - set(dcc.index))
    if missing:
        raise RuntimeError(f"no DCC for GSM(s): {missing}")

    qc = samples.set_index("gsm_id").join(dcc, how="left").reset_index()

    # GeoMx sequencing saturation: the fraction of aligned reads that collapse
    # onto already-seen molecules. 1 - dedup/aligned.
    qc["pct_aligned"] = qc["aligned"] / qc["raw"]
    qc["pct_trimmed_retained"] = qc["trimmed"] / qc["raw"]
    qc["sequencing_saturation"] = 1.0 - (qc["deduplicated"] / qc["aligned"])

    # ------------------------------------------------------ detection from matrix
    genes, values = [], []
    with gzip.open(snakemake.input.matrix, "rt", encoding="utf-8") as handle:
        aois = handle.readline().rstrip("\n").split("\t")[1:]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[0])
            values.append(fields[1:])
    full = pd.DataFrame(np.asarray(values, dtype=float), index=genes, columns=aois)
    neg = full.loc[NEG_PROBE]
    expr = full.drop(index=NEG_PROBE)

    qc = qc.set_index("aoi_label")
    qc["negprobe"] = neg
    qc["q3_norm_sum"] = expr.sum()  # post-normalisation; reported, never a threshold

    # The whole detection curve, not one number. The background multiple is
    # itself a threshold decision, so the notebook needs to slide over it.
    curve_multiples = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0]
    for mult in curve_multiples:
        qc[f"detection_at_{mult:g}x"] = (expr > neg * mult).mean()

    if background_multiple is None:
        emit("\ndetection_background_multiple is null — no single detection rate "
             "is chosen yet; the curve is written for the review notebook.")
        qc["gene_detection_rate"] = np.nan
    else:
        qc["gene_detection_rate"] = (expr > neg * background_multiple).mean()
        emit(f"\ndetection defined at >{background_multiple:g}x NegProbe-WTX")

    qc = qc.reset_index()

    # --------------------------------------------------------------- flagging
    reasons = {label: [] for label in qc["aoi_label"]}
    unset = [name for name, value in thresholds.items() if value is None]
    if background_multiple is None:
        unset.append("detection_background_multiple")

    if unset:
        emit(
            "\nTHRESHOLDS NOT SET: " + ", ".join(sorted(unset)) + "\n"
            "  Nothing is flagged. Metrics are written so the P0-T5 gate can "
            "choose values in notebooks/review/qc_review.py; they then go into\n"
            "  config/config.yaml with an ADR citing that notebook."
        )
    else:
        checks = [
            ("gene_detection_rate", thresholds["min_gene_detection_rate"],
             "gene detection rate below threshold"),
            ("raw", thresholds["min_library_size"],
             "raw sequencing reads below threshold"),
            ("sequencing_saturation", thresholds["min_sequencing_saturation"],
             "sequencing saturation below threshold"),
        ]
        for column, limit, why in checks:
            for label, value in zip(qc["aoi_label"], qc[column]):
                if pd.notna(value) and value < limit:
                    reasons[label].append(f"{why} ({value:.4g} < {limit:g})")

    qc["qc_flag"] = [len(reasons[label]) > 0 for label in qc["aoi_label"]]
    qc["qc_reason"] = ["; ".join(reasons[label]) for label in qc["aoi_label"]]

    columns = [
        "aoi_label", "gsm_id", "aoi_code", "patient_id", "compartment", "dsp_run",
        "raw", "trimmed", "stitched", "aligned", "deduplicated", "n_probes",
        "pct_aligned", "pct_trimmed_retained", "sequencing_saturation",
        "umiq30", "rtsq30", "negprobe", "q3_norm_sum",
        *[f"detection_at_{m:g}x" for m in curve_multiples],
        "gene_detection_rate", "qc_flag", "qc_reason",
    ]
    qc[columns].round(6).to_csv(snakemake.output.metrics, sep="\t", index=False)
    emit(f"\nwrote {snakemake.output.metrics}")

    # Flag, don't drop (§6). The file is named qc_excluded.tsv by the plan, but
    # while flag_only is true nothing in it is actually removed — the column
    # says so explicitly rather than leaving the filename to imply otherwise.
    flagged = qc[qc["qc_flag"]].copy()
    flagged["action"] = "flagged_retained" if flag_only else "excluded"
    flagged[["aoi_label", "gsm_id", "aoi_code", "patient_id", "qc_reason", "action"]].to_csv(
        snakemake.output.excluded, sep="\t", index=False
    )
    emit(f"wrote {snakemake.output.excluded} ({len(flagged)} flagged, "
         f"action={'flagged_retained' if flag_only else 'excluded'})")

    # ------------------------------------------- normalisation verification (Q1)
    col_q3 = expr.quantile(0.75)
    q3_cv = float(col_q3.std(ddof=1) / col_q3.mean())
    verified = bool(snakemake.params.norm_method == "verify" and q3_cv < 0.01)
    emit(f"\nnormalisation: method={snakemake.params.norm_method}, "
         f"column-Q3 CV={q3_cv:.4%} -> {'VERIFIED' if verified else 'NOT VERIFIED'}")

    # ----------------------------------------------------------- the 4-panel figure
    colours = [
        COMPARTMENT_COLOUR[
            qc.loc[qc["aoi_code"] == code, "compartment"].iloc[0]
        ] for code in CODE_ORDER
    ]

    def by_code(column):
        return [qc.loc[qc["aoi_code"] == code, column].to_numpy() for code in CODE_ORDER]

    def style(ax, bp):
        for patch, colour in zip(bp["boxes"], colours):
            patch.set_facecolor(colour)
            patch.set_alpha(0.75)
            patch.set_edgecolor("black")
            patch.set_linewidth(0.8)
        ax.set_xticks(range(1, len(CODE_ORDER) + 1))
        ax.set_xticklabels(CODE_ORDER, rotation=45, ha="right", fontsize=9)
        ax.grid(axis="y", alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4))
    box = dict(patch_artist=True, widths=0.62,
               medianprops=dict(color="black", lw=1.5),
               flierprops=dict(marker="o", ms=3, mfc="black", mec="none", alpha=0.6))

    ax = axes[0, 0]
    style(ax, ax.boxplot([v / 1e6 for v in by_code("raw")], **box))
    ax.set_title("Raw sequencing reads (DCC header)", fontsize=11, fontweight="bold")
    ax.set_ylabel("million reads", fontsize=9)

    ax = axes[0, 1]
    style(ax, ax.boxplot(by_code("sequencing_saturation"), **box))
    ax.axhline(0.5, ls="--", lw=1.0, color="crimson")
    ax.set_title("Sequencing saturation  (1 − dedup/aligned)", fontsize=11, fontweight="bold")
    ax.set_ylabel("saturation", fontsize=9)

    ax = axes[1, 0]
    detection_column = (
        "gene_detection_rate" if background_multiple is not None else "detection_at_2x"
    )
    shown_multiple = background_multiple if background_multiple is not None else 2.0
    style(ax, ax.boxplot([v * 100 for v in by_code(detection_column)], **box))
    ax.set_title(
        f"Gene detection rate  (> {shown_multiple:g}× NegProbe)"
        + ("" if background_multiple is not None else "  [provisional]"),
        fontsize=11, fontweight="bold",
    )
    ax.set_ylabel("% of 18,694 genes detected", fontsize=9)

    ax = axes[1, 1]
    style(ax, ax.boxplot(by_code("negprobe"), **box))
    ax.set_title("NegProbe-WTX level (background)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Q3-normalised value", fontsize=9)

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.75, ec="black", lw=0.8)
        for c in COMPARTMENT_COLOUR.values()
    ]
    fig.legend(handles, list(COMPARTMENT_COLOUR), loc="lower center", ncol=4,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.005))
    state = "thresholds unset" if unset else f"{int(qc['qc_flag'].sum())} AOI(s) flagged, none dropped"
    fig.suptitle(
        f"P0-T5 QC — 120 AOIs, TIME-B n=8 — {state}",
        fontsize=12.5, fontweight="bold", y=0.985,
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.96))
    fig.savefig(snakemake.output.figure, dpi=200, bbox_inches="tight")
    plt.close(fig)
    emit(f"wrote {snakemake.output.figure}")

    Path(snakemake.output.summary).write_text(
        json.dumps(
            {
                "n_aoi": int(len(qc)),
                "thresholds": thresholds,
                "detection_background_multiple": background_multiple,
                "thresholds_unset": sorted(unset),
                "flag_only": bool(flag_only),
                "n_flagged": int(qc["qc_flag"].sum()),
                "normalisation_method": snakemake.params.norm_method,
                "column_q3_cv": round(q3_cv, 8),
                "normalisation_verified": verified,
                "raw_reads_median": float(qc["raw"].median()),
                "raw_reads_min": float(qc["raw"].min()),
                "saturation_median": round(float(qc["sequencing_saturation"].median()), 6),
                "saturation_min": round(float(qc["sequencing_saturation"].min()), 6),
                "detection_curve_multiples": curve_multiples,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    emit(f"wrote {snakemake.output.summary}")
