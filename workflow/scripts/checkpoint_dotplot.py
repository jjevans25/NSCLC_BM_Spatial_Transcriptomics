"""P3-T5 — the Phase 3 deliverable figure: genes x (compartment x site).

Owner task: P3-T5. Driven by rule p3t5_checkpoint_dotplot.

Choices worth stating, because each could reasonably go the other way:

  * **This is a DETECTION figure and no model output reaches it.** P3-T3 and
    P3-T4 estimate effects on EXPRESSION; this plots DETECTION. The two move
    OPPOSITE ways under the same background gradient -- detection is
    `q3 > 2 x NegProbe-WTX`, so higher background raises the bar and a gene
    looks depleted; expression is `log2(q3+1)`, where background ADDS to the
    signal and a gene looks enriched. ADR 0018 records this project making that
    exact conflation once already, and inverting an argument on it. Putting an
    effect size on the same canvas as a detection rate is how it happens again.

  * **Colour is `median_negprobe_ratio`, not mean expression.** PROJECT_PLAN
    section 6 says "colour = mean expression"; ADR 0016 section 5, written
    before this figure existed, changed it. On a mean-log2 scale an UNDETECTED
    gene still shows its background level, and brain background is higher, so
    undetected genes would render BRIGHTER IN BRAIN than in lung -- the ADR 0008
    artefact walking into the deliverable. Both mean-log2 columns stay in
    checkpoint_detection.tsv. `mean_log2_detected` was rejected separately: its
    denominator changes per cell, so two dots of one colour could rest on 15
    AOIs and on 4.

  * **The colour midpoint IS the detection rule, read from config.** vcenter is
    `qc.detection_background_multiple` (ADR 0007), not a hardcoded 2.0, so the
    midpoint cannot drift from the rule it depicts. A sequential map's
    perceptual middle would land near 2.97 -- a number with no referent. The one
    meaningful number on this scale is the threshold.

  * **TwoSlopeNorm makes the colour scale PIECEWISE-LINEAR, and that is a real
    cost.** Roughly 1.16 ratio units below the centre and 3.11 above are each
    mapped to half the ramp, so colour distance is not ratio distance across the
    midpoint. Declared here, in the caption, and mitigated with explicit
    colourbar ticks plus a drawn rule at the threshold. The alternative --
    a sequential map with one tick at 2.0 -- puts a categorical boundary
    somewhere the eye does not read it.

  * **PRGn, deliberately not RdBu_r.** P2-T7 uses RdBu_r for a SIGNED EFFECT on
    a symmetric scale. Two phase deliverables must not share a colour language
    for two different quantities, which is the whole point of ADR 0016 section 5.
    PuOr_r was the considered alternative; its orange end sits close to the
    project's SITE_COLOUR palette.

  * **In greyscale the colour is ambiguous, and that is why nothing depends on
    it alone.** A diverging map sends BOTH ends to dark grey, so a dot at 0.84
    and one at 5.1 are indistinguishable without colour. Verified by
    desaturating the PNG. Every claim the figure makes is carried redundantly by
    something achromatic -- dot area, the hatch, the bold edge, and the printed
    counts -- so the greyscale reader loses the ratio and nothing else. That
    redundancy is the reason a diverging map is affordable here at all.

  * **Zero detection is an open RING, not a shrunken dot.** Area is proportional
    to detection rate with no minimum floor: a floor would make 0/30 and 1/20
    indistinguishable and turn the size key into a lie, in a figure whose entire
    subject is the counts. A ring reads as "an empty dot" -- a MEASURED ZERO,
    which is the strongest statement on the figure. An `x` would read as
    excluded, and would collide with the hatch's meaning.

  * **Not-assessable cells are hatched, and the hatch is load-bearing.** Same
    mark and same reasoning as contexture_heatmap.py: ADR 0008 point 4 says the
    finding is "not assessable in <compartment>", never "lower in
    <compartment>", and a dot plot is exactly the artefact someone reads a
    compartment difference off. 40 of 63 cells are hatched and that DOMINATES
    the figure, which is correct -- the panel clears the floor in TIME-L alone,
    and that is the phase's deliverable under Gate 3. Marking only the
    assessable cells was rejected: an unmarked dot reads as a confident zero in
    exactly the cells where absence and invisibility are not separable.

  * **Every count is printed under its dot.** Makes the figure self-verifying
    against checkpoint_detection.tsv, greyscale-safe, and unambiguous at the
    zeros -- the defence against the worst misreading, "the gene is not on the
    plot".

  * **Columns are grouped by site, tumour before immune within site.** This
    DEVIATES from the project's compartment-major CODE_ORDER (pca_landscape.py,
    qc_metrics.py); those are per-AOI QC figures with no site contrast to align.
    Here the axis IS compartment x site, and this order puts P3-T3's carrier
    contrast adjacent within each site and P3-T4's shift contrast positionally
    parallel between them.

  * **Single-DSP-run compartments carry a dagger, not a fade.** Every graphical
    encoding of "modelled vs not" -- hollow, faded, desaturated -- implies rank.
    A dagger reads as "see note", and the note carries the reason: batch is
    inseparable from biology there (ADR 0009 section 3, ADR 0016 section 2).
    That is a property of the study design, not a judgement about the data.

  * **Rows are ordered by how many compartments clear the pre-registered floor,
    never by a pooled detection rate.** Pooling a rate across compartments is
    the one operation ADR 0016 and checkpoint_detection.py both forbid, and it
    would be forbidden as a sort key too. The key used averages nothing, is
    built from `assessable`, and is PRINTED on the right-hand axis so the order
    explains itself and claims nothing.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

# Site-major, tumour before immune within site. See the docstring.
CODE_ORDER = ["L", "TIME-L", "mLN", "LB", "TIME-B", "TBME", "BC"]

# Reused verbatim from pca_landscape.py so one study has one site palette.
SITE_COLOUR = {"lung": "#C1666B", "brain": "#4F6D7A", "lymph_node": "#D4B483"}
SITE_LABEL = {"lung": "lung", "brain": "brain", "lymph_node": "lymph node"}

CMAP = "PRGn"
S_MAX = 260.0        # pt^2 at detection_rate = 1.0
S_ZERO = 26.0        # the ring; smaller than every drawn dot
DOT_EDGE = "#2B2B2B"
HATCH_EDGE = "#8A8A8A"

plt.rcParams["hatch.linewidth"] = 0.6

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    audit_codes = list(snakemake.params.audit_aoi_codes)
    model_codes = set(snakemake.params.model_aoi_codes)
    floor = float(snakemake.params.detected_in_aoi_fraction)
    vcenter = float(snakemake.params.background_multiple)
    expected_n_aoi = int(snakemake.params.expected_n_aoi)

    emit("P3-T5 — checkpoint detection dot plot")
    emit(f"  detection rule   q3 > {vcenter} x that AOI's NegProbe-WTX (ADR 0007)")
    emit(f"  gene floor       detected in >= {floor} of a compartment's AOIs")
    emit("                   (PRE-REGISTERED, ADR 0016 — from ADR 0008 point 4)")
    emit(f"  colour basis     median_negprobe_ratio (ADR 0016 §5), midpoint {vcenter}")
    emit()

    # ------------------------------------------------------------------ read
    # round_trip: ADR 0013's postscript — pandas' default C parser is not
    # correctly rounded, and assertion 3 compares detection_rate to the floor.
    det = pd.read_csv(
        snakemake.input.detection, sep="\t", float_precision="round_trip"
    )
    with open(snakemake.input.summary, encoding="utf-8") as handle:
        summary = json.load(handle)

    # -------------------------------------------------------------- assertions
    if len(det) != 63:
        raise RuntimeError(
            f"expected 63 rows (9 genes x 7 compartments), got {len(det)}. "
            "The audit changed shape; the figure's layout is no longer valid."
        )
    if set(det["aoi_code"]) != set(audit_codes):
        raise RuntimeError(
            f"compartments in the table {sorted(set(det['aoi_code']))} do not "
            f"match config audit_aoi_codes {sorted(audit_codes)}."
        )
    if set(CODE_ORDER) != set(audit_codes):
        raise RuntimeError(
            f"this script's CODE_ORDER {CODE_ORDER} does not cover "
            f"audit_aoi_codes {sorted(audit_codes)} — a compartment would be "
            "silently dropped from the deliverable figure."
        )

    per_gene_total = det.groupby("gene")["n_aoi"].sum()
    if not (per_gene_total == expected_n_aoi).all():
        raise RuntimeError(
            f"per-gene AOI totals {per_gene_total.to_dict()} do not all equal "
            f"{expected_n_aoi}. The figure must account for every AOI in the "
            "study; PROJECT_PLAN §2.1 is authoritative."
        )

    # Re-derive the pre-registered flag at plot time rather than trusting the
    # column. Catches a stale table against a changed floor.
    derived_assessable = det["detection_rate"] >= floor
    if not (derived_assessable == det["assessable"]).all():
        bad = det.loc[derived_assessable != det["assessable"], ["gene", "aoi_code"]]
        raise RuntimeError(
            f"`assessable` disagrees with detection_rate >= {floor} on "
            f"{bad.to_dict('records')}. The table and the pre-registered floor "
            "have diverged."
        )
    derived_modelled = det["aoi_code"].isin(model_codes)
    if not (derived_modelled == det["modelled"]).all():
        raise RuntimeError(
            "`modelled` disagrees with config model_aoi_codes."
        )

    if summary["detected_in_aoi_fraction"] != floor:
        raise RuntimeError(
            f"summary floor {summary['detected_in_aoi_fraction']} != config "
            f"{floor}: figure and table are under different numbers."
        )
    if summary["detection_background_multiple"] != vcenter:
        raise RuntimeError(
            f"summary background multiple {summary['detection_background_multiple']}"
            f" != config {vcenter}."
        )
    if not str(summary["figure_colour_basis"]).startswith("median_negprobe_ratio"):
        raise RuntimeError(
            "checkpoint_detection_summary.json does not declare "
            "median_negprobe_ratio as the figure colour basis. P3-T2 wrote that "
            "contract forward to this rule (ADR 0016 §5); do not plot something "
            "else because the string moved."
        )
    emit("assertions: 6 of 6 passed (shape, AOI totals, assessable, modelled,")
    emit("            floor agreement, colour-basis contract)")
    emit()

    # ------------------------------------------------------------- row order
    # Compartments clearing the pre-registered floor. Built from `assessable`,
    # averages nothing, and is printed on the figure.
    clears = det.groupby("gene")["assessable"].sum().astype(int)
    genes = sorted(clears.index, key=lambda g: (-clears[g], g))
    alias = det.drop_duplicates("gene").set_index("gene")["alias"]
    n_comp = len(CODE_ORDER)

    emit("row order — compartments clearing the floor (descending, ties alphabetical):")
    for g in genes:
        emit(f"  {g:8s} ({alias[g]:7s})  {clears[g]}/{n_comp}")
    emit()

    emit("column order — site-major, tumour before immune within site:")
    code_site = det.drop_duplicates("aoi_code").set_index("aoi_code")
    for c in CODE_ORDER:
        emit(
            f"  {c:8s} {code_site.at[c, 'site']:11s} "
            f"{code_site.at[c, 'compartment']:14s} n = {int(code_site.at[c, 'n_aoi']):3d}"
            f"   {'modelled' if c in model_codes else 'audit only (single DSP run) †'}"
        )
    emit()

    # ------------------------------------------------------------- the grid
    idx = det.set_index(["gene", "aoi_code"])
    emit("detection per gene x compartment (n_detected/n_aoi), as plotted:")
    emit()
    header = "  gene     alias   " + "".join(f"{c:>9s}" for c in CODE_ORDER)
    emit(header)
    emit("  " + "-" * (len(header) - 2))
    for g in genes:
        cells = "".join(
            f"{str(int(idx.at[(g, c), 'n_detected'])) + '/' + str(int(idx.at[(g, c), 'n_aoi'])):>9s}"
            for c in CODE_ORDER
        )
        emit(f"  {g:8s} {alias[g]:7s} {cells}")
    emit()

    ratio = det["median_negprobe_ratio"]
    # Guard the norm: TwoSlopeNorm raises if vmin >= vcenter, and a rerun on a
    # different subset could put every cell on one side.
    vmin = min(float(ratio.min()), 0.95 * vcenter)
    vmax = max(float(ratio.max()), 1.05 * vcenter)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    emit(f"colour: vmin {vmin:.4f}  vcenter {vcenter}  vmax {vmax:.4f}")
    near = det.assign(_d=(ratio - vcenter).abs()).nsmallest(3, "_d")
    emit("  nearest the midpoint (these need the dot edge to stay visible):")
    for r in near.itertuples(index=False):
        emit(f"    {r.gene:8s} {r.aoi_code:8s} ratio {r.median_negprobe_ratio:.4f}")
    emit("  NOTE TwoSlopeNorm is piecewise-linear: "
         f"{vcenter - vmin:.2f} units below the midpoint and {vmax - vcenter:.2f} "
         "above each map to half the ramp, so colour distance is NOT ratio "
         "distance across it.")
    emit()

    hatched = det.loc[~det["assessable"], ["gene", "aoi_code"]]
    rings = det.loc[det["n_detected"] == 0, ["gene", "aoi_code"]]
    emit(f"not assessable (hatched): {len(hatched)} of {len(det)} cells")
    emit("  " + ", ".join(f"{r.gene}/{r.aoi_code}" for r in hatched.itertuples(index=False)))
    emit(f"detected in zero AOIs (ring): {len(rings)} cells")
    emit("  " + ", ".join(f"{r.gene}/{r.aoi_code}" for r in rings.itertuples(index=False)))
    emit()

    # Observation, never an assertion: the colour side and the hatch are the
    # same pre-registered rule seen twice, and they agree on every cell today.
    # They can legitimately diverge when n_aoi is even and the rate sits exactly
    # at the floor -- the two middle values then straddle the threshold.
    discordant = det.loc[(det["median_negprobe_ratio"] > vcenter) != det["assessable"]]
    emit(f"colour/hatch concordance: {len(discordant)} discordant cells of {len(det)}")
    emit("  (ratio > midpoint) and (assessable) are the same rule seen twice.")
    emit("  They CAN diverge where n_aoi is even and the rate is exactly at the")
    emit("  floor; this is logged as an observation, never asserted.")
    if len(discordant):
        for r in discordant.itertuples(index=False):
            emit(f"    {r.gene}/{r.aoi_code} ratio {r.median_negprobe_ratio:.4f} "
                 f"assessable={r.assessable}")

    flips = det.loc[
        det["detection_rate_qc_clean"].notna()
        & ((det["detection_rate_qc_clean"] >= floor) != det["assessable"])
    ]
    emit(f"QC-clean stability: {len(flips)} cells change assessability when the")
    emit("  5 QC-flagged AOIs are excluded (flag-don't-drop; the figure uses all).")
    emit()

    # ------------------------------------------------------------------ figure
    x_of = {c: i for i, c in enumerate(CODE_ORDER)}
    y_of = {g: i for i, g in enumerate(genes)}

    # Column order is main | key | colourbar. The twinx that carries the sort
    # key attaches to the main axes' right spine, so anything immediately beside
    # it collides with those labels; wspace buys them room, and the colourbar
    # goes last where its own label has a clear margin.
    fig = plt.figure(figsize=(13.8, 7.8))
    gs = fig.add_gridspec(
        2, 3, height_ratios=[0.05, 1], width_ratios=[1, 0.30, 0.026],
        hspace=0.06, wspace=0.34,
        left=0.115, right=0.878, top=0.735, bottom=0.125,
    )
    ax_site = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0])
    ax_key = fig.add_subplot(gs[1, 1])
    ax_cbar = fig.add_subplot(gs[1, 2])
    ax_key.set_axis_off()

    # site strip
    site_of = {c: code_site.at[c, "site"] for c in CODE_ORDER}
    for i, c in enumerate(CODE_ORDER):
        ax_site.add_patch(
            Rectangle((i, 0), 1, 1, facecolor=SITE_COLOUR[site_of[c]], lw=0)
        )
    ax_site.set_xlim(0, n_comp)
    ax_site.set_ylim(0, 1)
    ax_site.set_xticks([])
    ax_site.set_yticks([])
    for spine in ax_site.spines.values():
        spine.set_visible(False)

    # site labels + boundaries
    blocks = []
    start = 0
    for i, c in enumerate(CODE_ORDER):
        if i + 1 == n_comp or site_of[CODE_ORDER[i + 1]] != site_of[c]:
            blocks.append((site_of[c], start, i + 1))
            start = i + 1
    site_n = det.drop_duplicates("aoi_code").groupby("site")["n_aoi"].sum()
    for site, a, b in blocks:
        ax_site.text(
            (a + b) / 2, 1.5, f"{SITE_LABEL[site]}  (n = {site_n[site]})",
            ha="center", va="bottom", fontsize=9.5,
            color=SITE_COLOUR[site], fontweight="bold",
        )
        if b < n_comp:
            ax_site.axvline(b, color="white", lw=2.0)
            ax.axvline(b - 0.5, color="black", lw=1.6, zorder=2)

    # row bands, so the eye tracks a gene across seven columns
    for g in genes:
        if y_of[g] % 2 == 0:
            ax.add_patch(
                Rectangle((-0.5, y_of[g] - 0.5), n_comp, 1,
                          facecolor="#F7F7F7", lw=0, zorder=0)
            )

    # hatched cells: not assessable
    for r in hatched.itertuples(index=False):
        ax.add_patch(
            Rectangle(
                (x_of[r.aoi_code] - 0.5, y_of[r.gene] - 0.5), 1, 1,
                fill=False, hatch="///", edgecolor=HATCH_EDGE, lw=0.0,
                alpha=0.45, zorder=1,
            )
        )

    halo = [pe.withStroke(linewidth=2.2, foreground="white")]

    # zero-detection rings
    zx = [x_of[r.aoi_code] for r in rings.itertuples(index=False)]
    zy = [y_of[r.gene] for r in rings.itertuples(index=False)]
    if zx:
        ring = ax.scatter(
            zx, zy, s=S_ZERO, marker="o", facecolors="none",
            edgecolors="#3A3A3A", linewidths=1.0, zorder=3,
        )
        ring.set_path_effects(halo)

    # the dots
    drawn = det.loc[det["n_detected"] > 0]
    dx = [x_of[c] for c in drawn["aoi_code"]]
    dy = [y_of[g] for g in drawn["gene"]]
    sizes = S_MAX * drawn["detection_rate"].to_numpy()
    edges = np.where(drawn["assessable"].to_numpy(), 1.30, 0.55)
    sc = ax.scatter(
        dx, dy, s=sizes, c=drawn["median_negprobe_ratio"].to_numpy(),
        cmap=CMAP, norm=norm, edgecolors=DOT_EDGE, linewidths=edges, zorder=3,
    )
    sc.set_path_effects(halo)

    # printed counts
    for r in det.itertuples(index=False):
        t = ax.text(
            x_of[r.aoi_code], y_of[r.gene] + 0.33,
            f"{int(r.n_detected)}/{int(r.n_aoi)}",
            ha="center", va="center", fontsize=5.4, color="#4A4A4A", zorder=5,
        )
        t.set_path_effects(halo)

    ax.set_xlim(-0.5, n_comp - 0.5)
    ax.set_ylim(len(genes) - 0.5, -0.5)
    ax.set_xticks(range(n_comp))
    ax.set_xticklabels(
        [
            f"{c}{'' if c in model_codes else ' †'}\nn = {int(code_site.at[c, 'n_aoi'])}"
            for c in CODE_ORDER
        ],
        fontsize=9,
    )
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([f"{g}  ({alias[g]})" for g in genes], fontsize=9.5)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_color("#BBBBBB")

    # right-hand axis: the sort key, printed
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(range(len(genes)))
    ax2.set_yticklabels([f"{clears[g]}/{n_comp}" for g in genes], fontsize=8)
    # A short header above the column rather than a rotated ylabel: the ylabel
    # ran straight into the key panel.
    ax2.text(
        1.004, 1.035, f"clears the\n{floor:g} floor",
        transform=ax.transAxes, fontsize=7.4, ha="left", va="bottom",
        color="#444444", linespacing=1.2,
    )
    ax2.tick_params(length=0, pad=2)
    for spine in ax2.spines.values():
        spine.set_visible(False)

    # colourbar
    cb = fig.colorbar(sc, cax=ax_cbar)
    ticks = sorted({
        round(vmin, 2), 1.0, vcenter, 3.0, 4.0, 5.0, round(vmax, 2)
    })
    cb.ax.tick_params(labelsize=7.4)
    cb.ax.plot(
        [0, 1], [vcenter, vcenter], color="black", lw=1.4,
        transform=cb.ax.get_yaxis_transform(), clip_on=False,
    )
    cb.ax.plot(
        [0, 1], [1.0, 1.0], color="0.35", lw=0.9, ls=":",
        transform=cb.ax.get_yaxis_transform(),
    )
    # The threshold and the background level are named IN THE TICK LABELS.
    # Side annotations were tried and collided with the key panel; a tick label
    # cannot drift away from the value it names.
    shown = [t for t in ticks if vmin <= t <= vmax]
    cb.set_ticks(shown)
    cb.set_ticklabels([
        f"{t:g} ← threshold" if t == vcenter
        else (f"{t:g} ← background" if t == 1.0 else f"{t:g}")
        for t in shown
    ])
    cb.set_label("median  q3 / NegProbe-WTX  (ADR 0016 §5)", fontsize=8.0,
                 labelpad=4)

    # ------------------------------------------------------------------- key
    # Drawn with the SAME sizes and the SAME scatter call as the field, so the
    # key cannot drift from what it describes.
    ax_key.set_xlim(0, 1)
    ax_key.set_ylim(0, 1)
    key_rates = [0.25, 0.50, 0.75, 1.00]
    ys = np.linspace(0.80, 0.56, len(key_rates))
    ax_key.scatter(
        [0.18] * len(key_rates), ys, s=[S_MAX * r for r in key_rates],
        facecolor="#DDDDDD", edgecolors=DOT_EDGE, linewidths=0.55,
    )
    for r, y in zip(key_rates, ys):
        ax_key.text(0.36, y, f"{r:.0%} of AOIs", fontsize=7.5, va="center")
    ax_key.text(0.02, 0.89, "dot area = detection rate", fontsize=8,
                fontweight="bold")
    ax_key.scatter([0.18], [0.44], s=S_ZERO, marker="o", facecolors="none",
                   edgecolors="#3A3A3A", linewidths=1.0)
    ax_key.text(0.36, 0.44, "detected in 0 AOIs", fontsize=7.5, va="center")
    ax_key.add_patch(
        Rectangle((0.10, 0.26), 0.16, 0.10, fill=False, hatch="///",
                  edgecolor=HATCH_EDGE, lw=0.0, alpha=0.45)
    )
    ax_key.text(0.36, 0.31, f"below the {floor:g} floor —\n\"not assessable\" there",
                fontsize=7.5, va="center")
    ax_key.scatter([0.18], [0.14], s=S_MAX * 0.6, facecolor="#DDDDDD",
                   edgecolors=DOT_EDGE, linewidths=1.30)
    ax_key.text(0.36, 0.14, "bold edge = clears\nthe floor", fontsize=7.5,
                va="center")

    # ------------------------------------------------------------------ title
    # KEEP EVERY LINE SHORT. P2-T7's docstring records a caption line that ran
    # past both canvas edges and clipped the n = 8 the figure is required to
    # state; the first draft of this figure reproduced it exactly. Verify by
    # opening the PNG, never by assuming.
    lines = [
        "Checkpoint panel detection by compartment — 9 pre-registered genes × "
        f"{n_comp} compartments, all {expected_n_aoi} AOIs",
        "(P3-T5, Aim A4 — exploratory per ADR 0008)",
        "dot area = fraction of that compartment's AOIs above background   ·   "
        f"colour = median q3 / NegProbe-WTX, midpoint at the {vcenter:g}× rule",
        "ring = detected in 0 AOIs   ·   hatched = below the pre-registered "
        f"{floor:g} floor, \"not assessable\" there (ADR 0016)",
        "TIME-B n = 8; the lung-vs-brain immune contrast detects ~1.1–1.3 SD at "
        "80% power (P0-T8)",
        "— a small dot is an assay-sensitivity limit, NOT evidence the gene is "
        "absent",
        "† detection only: each sits wholly within one DSP run, so batch is "
        "inseparable from biology there (ADR 0009 §3)",
    ]
    fig.suptitle("\n".join(lines), fontsize=8.4, y=0.995,
                 va="top", linespacing=1.45)

    fig.savefig(snakemake.output.figure, dpi=200)
    plt.close(fig)

    emit(f"wrote {snakemake.output.figure}")
