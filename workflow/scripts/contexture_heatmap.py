"""P2-T7 — the Phase 2 deliverable figure: signatures x AOIs, grouped by site.

Owner task: P2-T7. Driven by rule p2t7_contexture_heatmap.

Choices worth stating, because each could reasonably go the other way:

  * **The z-score mean is plotted, not ssGSEA.** Both are computed and both are
    modelled; only one can be a heatmap. The z-score mean is on a common,
    centred, interpretable scale across signatures (SD units about the
    23-AOI mean), whereas ssGSEA's NES is a rank-enrichment statistic whose
    magnitude is not comparable between sets of different size. Plotting the
    latter would put a 4-gene and a 13-gene signature on one colour scale and
    invite a comparison the statistic does not support. The ssGSEA result is in
    `signature_scores.tsv` and in every model row.

  * **Columns are grouped by site, then ordered by patient**, rather than
    clustered. A dendrogram here would let the display choose the grouping,
    which is the question the figure is supposed to answer.

  * **Colour is symmetric about zero** on a ColorBrewer RdBu diverging map,
    which is colourblind-safe. The limit is the largest absolute score, so no
    value is clipped and the midpoint is a real zero rather than the data mean.

  * **Sets that failed their coverage floor are marked on the figure itself**,
    not just in the caption. ADR 0008's rule is that such a set is "not
    assessable" at that site, never "lower" there — and a heatmap is exactly
    the artefact from which someone reads "lower in brain" off a colour. The
    mark is the defence against its own most likely misreading.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

SITE_COLOUR = {"TIME-L": "#4F6D7A", "TIME-B": "#C1666B"}

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    aoi_codes = list(snakemake.params.aoi_codes)

    scores = pd.read_csv(snakemake.input.scores, sep="\t")
    coverage = pd.read_csv(snakemake.input.coverage, sep="\t")

    d = scores.loc[scores["method"] == "zscore"].copy()
    emit(f"plotting method=zscore, {len(d)} rows")

    # Columns: site block, then patient within block.
    order = (
        d[["aoi_label", "aoi_code", "patient_id"]]
        .drop_duplicates()
        .assign(_k=lambda x: x["aoi_code"].map({c: i for i, c in enumerate(aoi_codes)}))
        .sort_values(["_k", "patient_id", "aoi_label"])
    )
    cols = order["aoi_label"].tolist()

    mat = (
        d.pivot_table(index="signature", columns="aoi_label", values="score")
        .loc[sorted(d["signature"].unique()), cols]
    )
    emit(f"matrix: {mat.shape[0]} signatures x {mat.shape[1]} AOIs")

    counts = order["aoi_code"].value_counts().to_dict()
    emit(f"AOIs per site: {counts}")

    vmax = float(np.nanmax(np.abs(mat.to_numpy())))
    emit(f"symmetric colour limit: +/-{vmax:.3f}")

    not_assessable = {
        (r["signature"], r["aoi_code"])
        for _, r in coverage.iterrows()
        if not r["assessable"]
    }
    emit(f"not-assessable cells to mark: {sorted(not_assessable)}")

    # ------------------------------------------------------------------ figure
    fig = plt.figure(figsize=(13.0, 6.0))
    gs = fig.add_gridspec(
        2, 2, height_ratios=[0.055, 1], width_ratios=[1, 0.018],
        hspace=0.06, wspace=0.02, left=0.17, right=0.93, top=0.82, bottom=0.19,
    )
    ax_site = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0])
    ax_cbar = fig.add_subplot(gs[1, 1])

    # site strip
    for i, label in enumerate(cols):
        code = order.loc[order["aoi_label"] == label, "aoi_code"].iloc[0]
        ax_site.add_patch(Rectangle((i, 0), 1, 1, facecolor=SITE_COLOUR[code], lw=0))
    ax_site.set_xlim(0, len(cols))
    ax_site.set_ylim(0, 1)
    ax_site.set_xticks([])
    ax_site.set_yticks([])
    for spine in ax_site.spines.values():
        spine.set_visible(False)

    boundaries = np.cumsum([counts[c] for c in aoi_codes])[:-1]
    start = 0
    for code in aoi_codes:
        n = counts[code]
        ax_site.text(
            start + n / 2, 1.55, f"{code}  (n = {n})",
            ha="center", va="bottom", fontsize=10, color=SITE_COLOUR[code],
            fontweight="bold",
        )
        start += n

    im = ax.imshow(
        mat.to_numpy(), aspect="auto", cmap="RdBu_r",
        vmin=-vmax, vmax=vmax, interpolation="nearest",
    )
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(order["patient_id"].tolist(), rotation=90, fontsize=7)
    ax.set_xlabel("patient", fontsize=10)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels([s.replace("_", " ") for s in mat.index], fontsize=10)

    for b in boundaries:
        ax.axvline(b - 0.5, color="black", lw=2.0)
        ax_site.axvline(b, color="white", lw=2.0)

    # Cross-hatch the cells whose set is not assessable at that site.
    for r, sig in enumerate(mat.index):
        for c, label in enumerate(cols):
            code = order.loc[order["aoi_label"] == label, "aoi_code"].iloc[0]
            if (sig, code) in not_assessable:
                ax.add_patch(
                    Rectangle(
                        (c - 0.5, r - 0.5), 1, 1, fill=False, hatch="xxx",
                        edgecolor="black", lw=0.0, alpha=0.55,
                    )
                )

    cb = fig.colorbar(im, cax=ax_cbar)
    cb.set_label("signature score\n(mean gene z, SD units)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    lines = [
        "Immune contexture, brain metastasis vs. primary lung — within the TIME compartment",
        f"TIME-L n = {counts['TIME-L']}   ·   TIME-B n = {counts['TIME-B']}"
        "   ·   minimum detectable effect ~1.1–1.3 SD at 80% power (P0-T8)",
    ]
    if not_assessable:
        # Wrapped onto its own line rather than appended: with six signatures
        # the single-line form ran past both edges of the canvas and clipped
        # the n = 8 the figure is required to state.
        by_set = {}
        for sig, code in sorted(not_assessable):
            by_set.setdefault(sig, []).append(code)
        marked = "; ".join(
            f"{s.replace('_', ' ')} ({', '.join(c)})" for s, c in sorted(by_set.items())
        )
        lines.append(
            f"hatched = below detection-coverage floor, not assessable there: {marked}"
        )
    fig.suptitle("\n".join(lines), fontsize=10, y=0.985)

    fig.savefig(snakemake.output.figure, dpi=200)
    plt.close(fig)

    emit(f"wrote {snakemake.output.figure}")
