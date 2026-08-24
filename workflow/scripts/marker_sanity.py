"""P0-T6 — marker sanity check: do the compartment labels mean what they say?

Owner task: P0-T6. Driven by rule p0t6_marker_sanity.

PROJECT_PLAN §6 calls this "the 20-minute check that occasionally saves a
month", and Q2's resolution is what makes it load-bearing rather than
ceremonial. The compartments are real antibody segmentation, but **PanCK was
the only collection mask** (docs/data-provenance.md): CD45 and GFAP guided
where a pathologist placed the ROI, they did not sort cells. So `TIME` and
`TBME` labels rest on human placement plus a PanCK-negative mask, and the
transcriptome is the only independent evidence that the placement worked.

Acceptance criteria, verbatim from §6:

    tumour AOIs epithelial-high, `TIME` AOIs `PTPRC`-high, `TBME` AOIs
    `GFAP`-high. If it fails, stop the project and re-derive the compartment
    mapping.

Markers are exactly the five §6 names (EPCAM/KRT19, PTPRC, GFAP/AQP4). Nothing
is added to that panel here — that would be a "stop and ask" decision.
`NegProbe-WTX` is carried alongside as the **background floor**, not as a
marker: without it "lower" and "absent" are indistinguishable, and the
interesting claim is that epithelial signal in non-tumour compartments sits at
noise rather than merely below tumour.

**This runs before P0-T5.** That is deliberate: it tests the *labels*, not the
data quality, and a wrong mapping should be found before effort goes into QC.
Re-run it after P0-T5 to confirm QC did not change the picture.

Statistics: the criteria are ordinal (group > rest), so the log2 transform —
being monotone and global — cannot change any verdict. The per-AOI Q3 scaling
that the comparison *does* depend on was applied by the submitters (Q1). No
p-values: effect sizes are reported as a difference in median log2 with a
seeded bootstrap CI and an n, per hard constraint 7.
"""

import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MARKERS = ["EPCAM", "KRT19", "PTPRC", "GFAP", "AQP4"]
NEG_PROBE = "NegProbe-WTX"
PANEL = MARKERS + [NEG_PROBE]

# Display order: tumour cores, then immune, then glial, then control.
CODE_ORDER = ["L", "LB", "mLN", "TIME-L", "TIME-B", "TBME", "BC"]

# The three §6 criteria, as (marker, codes that should be high, label).
CRITERIA = [
    ("EPCAM", ["L", "LB", "mLN"], "tumour AOIs epithelial-high (EPCAM)"),
    ("KRT19", ["L", "LB", "mLN"], "tumour AOIs epithelial-high (KRT19)"),
    ("PTPRC", ["TIME-L", "TIME-B"], "TIME AOIs PTPRC-high"),
    ("GFAP", ["TBME"], "TBME AOIs GFAP-high"),
]
# Supporting, not gating: §6 names AQP4 alongside GFAP but states the criterion
# in terms of GFAP only.
SUPPORTING = [("AQP4", ["TBME", "BC"], "brain parenchyma AQP4-high (supporting)")]

COMPARTMENT_COLOUR = {
    "tumour": "#B0763F",
    "immune": "#3E7CB1",
    "glial_stroma": "#5B8C5A",
    "normal_control": "#8C8C8C",
}

N_BOOT = 2000
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)


def median_diff_ci(rng, high, rest, n_boot=N_BOOT):
    """Difference in median log2 (high - rest), with a bootstrap 95% CI."""
    point = float(np.median(high) - np.median(rest))
    draws = np.empty(n_boot)
    for i in range(n_boot):
        draws[i] = np.median(
            rng.choice(high, size=high.size, replace=True)
        ) - np.median(rng.choice(rest, size=rest.size, replace=True))
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return point, float(lo), float(hi)


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    seed = snakemake.params.seed
    rng = np.random.default_rng(seed)
    emit(f"seed: {seed} (config['seed'], passed explicitly)")

    # ------------------------------------------------------------ load
    samples = pd.read_csv(snakemake.input.samples, sep="\t", dtype=str)
    samples = samples.set_index("aoi_label")

    wanted = set(PANEL)
    found = {}
    with gzip.open(snakemake.input.matrix, "rt", encoding="utf-8") as handle:
        columns = handle.readline().rstrip("\n").split("\t")[1:]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in wanted:
                found[fields[0]] = np.asarray(fields[1:], dtype=float)

    absent = sorted(wanted - set(found))
    if absent:
        raise RuntimeError(
            f"marker(s) absent from the expression matrix: {absent}. "
            "The panel is fixed by PROJECT_PLAN §6 — do not substitute."
        )

    expr = pd.DataFrame(found, index=columns)
    if set(expr.index) != set(samples.index):
        raise RuntimeError(
            "matrix columns and samples.tsv aoi_labels disagree — P0-T3 asserts "
            "this join, so something has changed underneath it."
        )
    expr = expr.loc[samples.index]

    # log2(x+1). Monotone and global, so it cannot flip any verdict; it only
    # makes the boxplots legible.
    lg = np.log2(expr + 1.0)
    lg["aoi_code"] = samples["aoi_code"]
    lg["compartment"] = samples["compartment"]

    floor = float(np.median(lg[NEG_PROBE]))
    emit(f"background floor: median log2 {NEG_PROBE} = {floor:.2f}\n")

    # -------------------------------------------------- per-group statistics
    stats = []
    for marker in PANEL:
        for code in CODE_ORDER:
            values = lg.loc[lg["aoi_code"] == code, marker].to_numpy()
            q1, med, q3 = np.percentile(values, [25, 50, 75])
            stats.append(
                {
                    "marker": marker,
                    "aoi_code": code,
                    "n": int(values.size),
                    "median_log2": round(float(med), 4),
                    "q1_log2": round(float(q1), 4),
                    "q3_log2": round(float(q3), 4),
                    "median_over_background": round(float(med - floor), 4),
                }
            )
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(snakemake.output.stats, sep="\t", index=False)
    emit(f"wrote {snakemake.output.stats}")

    # ------------------------------------------------------------- verdicts
    verdicts = []
    for marker, high_codes, label in CRITERIA + SUPPORTING:
        gating = (marker, high_codes, label) in CRITERIA
        mask = lg["aoi_code"].isin(high_codes)
        high = lg.loc[mask, marker].to_numpy()
        rest = lg.loc[~mask, marker].to_numpy()
        delta, lo, hi = median_diff_ci(rng, high, rest)
        # "High" means both directions agree: higher than the other AOIs, and
        # above the negative-probe floor. A marker can be relatively higher and
        # still be background, which would not support the label.
        passed = bool(lo > 0 and np.median(high) > floor)
        verdicts.append(
            {
                "criterion": label,
                "gating": gating,
                "marker": marker,
                "high_in": "+".join(high_codes),
                "n_high": int(high.size),
                "n_rest": int(rest.size),
                "median_high_log2": round(float(np.median(high)), 4),
                "median_rest_log2": round(float(np.median(rest)), 4),
                "delta_median_log2": round(delta, 4),
                "ci95_low": round(lo, 4),
                "ci95_high": round(hi, 4),
                "above_background": bool(np.median(high) > floor),
                "verdict": "PASS" if passed else "FAIL",
            }
        )
    verdict_df = pd.DataFrame(verdicts)
    verdict_df.to_csv(snakemake.output.verdict, sep="\t", index=False)
    emit(f"wrote {snakemake.output.verdict}\n")

    for row in verdicts:
        tag = "GATING " if row["gating"] else "support"
        emit(
            f"  [{row['verdict']}] {tag} {row['criterion']}\n"
            f"           median {row['median_high_log2']:.2f} vs {row['median_rest_log2']:.2f}"
            f"  Δ {row['delta_median_log2']:+.2f} log2"
            f"  95% CI [{row['ci95_low']:+.2f}, {row['ci95_high']:+.2f}]"
            f"  n {row['n_high']}/{row['n_rest']}"
        )

    gating_failed = [r for r in verdicts if r["gating"] and r["verdict"] == "FAIL"]

    # ---------------------------------------------------------------- figure
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.6), sharex=True)
    for ax, marker in zip(axes.ravel(), PANEL):
        groups = [lg.loc[lg["aoi_code"] == c, marker].to_numpy() for c in CODE_ORDER]
        bp = ax.boxplot(
            groups, patch_artist=True, widths=0.62, medianprops=dict(color="black", lw=1.6),
            flierprops=dict(marker="o", ms=2.5, mfc="black", mec="none", alpha=0.55),
        )
        for patch, code in zip(bp["boxes"], CODE_ORDER):
            compartment = samples.loc[samples["aoi_code"] == code, "compartment"].iloc[0]
            patch.set_facecolor(COMPARTMENT_COLOUR[compartment])
            patch.set_alpha(0.75)
            patch.set_edgecolor("black")
            patch.set_linewidth(0.8)
        ax.axhline(floor, ls="--", lw=1.0, color="crimson", zorder=0)
        is_neg = marker == NEG_PROBE
        ax.set_title(
            f"{marker}{'  (background control)' if is_neg else ''}",
            fontsize=11, fontweight="bold" if not is_neg else "normal",
        )
        ax.set_ylabel("log2(Q3-normalised + 1)", fontsize=9)
        ax.grid(axis="y", alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_xticks(range(1, len(CODE_ORDER) + 1))
        ax.set_xticklabels(CODE_ORDER, rotation=45, ha="right", fontsize=9)

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=colour, alpha=0.75, ec="black", lw=0.8)
        for colour in COMPARTMENT_COLOUR.values()
    ]
    handles.append(plt.Line2D([0], [0], ls="--", color="crimson", lw=1.0))
    fig.legend(
        handles,
        list(COMPARTMENT_COLOUR) + [f"median {NEG_PROBE} (background)"],
        loc="lower center", ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.005),
    )
    banner = "PASS" if not gating_failed else "FAIL"
    fig.suptitle(
        f"P0-T6 marker sanity check — {banner}    "
        f"(n=120 AOIs; TIME-B n=8; pre-QC, seed {seed})",
        fontsize=12.5, fontweight="bold", y=0.985,
    )
    fig.tight_layout(rect=(0, 0.055, 1, 0.96))
    fig.savefig(snakemake.output.figure, dpi=200, bbox_inches="tight")
    plt.close(fig)
    emit(f"\nwrote {snakemake.output.figure}")

    Path(snakemake.output.summary).write_text(
        json.dumps(
            {
                "seed": seed,
                "n_boot": N_BOOT,
                "background_floor_log2": round(floor, 4),
                "gating_criteria": len(CRITERIA),
                "gating_passed": len(CRITERIA) - len(gating_failed),
                "verdict": banner,
                "failed": [r["criterion"] for r in gating_failed],
                "pre_qc": True,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    # This rule does not raise on failure. Snakemake deletes the outputs of a
    # failed job, which would destroy the very evidence needed to decide what to
    # do — and §6's instruction ("stop the project and re-derive the compartment
    # mapping") is a human call at the gate, not something a rule should take.
    # The verdict is in the outputs, in the figure title, and in this banner.
    rule_line = "=" * 74
    if gating_failed:
        message = (
            f"\n{rule_line}\nSTOP — MARKER SANITY CHECK FAILED\n{rule_line}\n"
            + "\n".join(f"  FAILED: {r['criterion']}" for r in gating_failed)
            + "\n\nPROJECT_PLAN §6: 'If it fails, stop the project and re-derive the\n"
            "compartment mapping.' Do not proceed to Phase 1 on these labels.\n"
            f"{rule_line}"
        )
    else:
        message = (
            f"\n{rule_line}\n"
            f"MARKER SANITY CHECK PASSED — {len(CRITERIA)}/{len(CRITERIA)} gating criteria.\n"
            "The compartment labels are supported by the transcriptome.\n"
            f"{rule_line}"
        )
    emit(message)
    print(message)
