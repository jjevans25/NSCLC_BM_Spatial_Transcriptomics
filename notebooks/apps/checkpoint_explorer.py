# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo",
#     "pandas",
#     "matplotlib",
#     "numpy",
# ]
# ///

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.md(r"""
    # P3-T6 — checkpoint explorer

    **Tier:** app (`notebooks/apps/`) — **presentation logic only** (`CLAUDE.md`
    hard constraint 4). Every number here was produced by a Snakemake rule.
    This notebook reads `p3t2_checkpoint_detection`'s and
    `p3t3`/`p3t4`'s outputs and computes nothing of its own except which cells
    to hatch at the floor you choose.

    **What the plot is.** The nine pre-registered checkpoint genes (ADR 0015)
    × all seven compartments, over all 120 AOIs. Dot area is the fraction of
    that compartment's AOIs in which the gene sits above background — the one
    detection rule this project has, `q3 > 2.0 ×` that AOI's `NegProbe-WTX`
    (ADR 0007). Colour is the median of that ratio, **not** mean expression: on
    a mean-log2 scale an *undetected* gene still shows its background level, and
    brain background is higher, so undetected genes would render brighter in
    brain than in lung — the ADR 0008 artefact walking into the figure
    (ADR 0016 §5).

    **What to do with it**, and it is the Phase 3 deliverable: raise the
    detection floor and watch what survives. At the pre-registered 0.5,
    **40 of 63 cells are not assessable and the panel clears the floor in
    `TIME-L` alone**. Raise the floor and `TIME-L` is the last compartment
    standing; `TBME` is hatched at every value. That stability is the argument
    the static figure can only assert.

    **The slider is a sensitivity display, not a threshold control**
    (ADR 0019). 0.5 is pre-registered (ADR 0016) and is the floor of record for
    every table, model and sentence in this phase. Nothing you can reach with
    the slider is a result.

    `TIME-B` **n = 8**. The lung-vs-brain immune contrast detects roughly
    **1.1–1.3 SD at 80% power** (P0-T8), so a small dot is an assay-sensitivity
    limit and **never** evidence that the gene is absent.
    """)
    return


@app.cell
def _(mo):
    # Parameters (PROJECT_PLAN §A.5), resolved for three different homes.
    #
    # Identical mechanism to notebooks/apps/landscape_explorer.py, and for the
    # reason ADR 0010 gives: the exporting rule stages this notebook next to a
    # `public/` directory holding every pipeline output it reads, because
    # `marimo export html-wasm` copies `public/` into the export.
    # `mo.notebook_location()` resolves to that directory locally and to the
    # served URL in the browser, so the WASM build fetches its data over HTTP
    # instead of needing a filesystem.
    #
    # Bundling the frames via `[tool.marimo.runtime] cache_cells` cannot work
    # here: an `mo.ui.table` holds a locally defined `PandasTableManager` that
    # cannot be pickled, so its cell output caches as an `UnhashableStub`, and
    # an unserializable output makes marimo re-run the cell's ancestors live.
    # In the browser that ancestor is the loader, which has no files to read.
    args = mo.cli_args()
    location = mo.notebook_location()

    def resolve(key, filename):
        override = args.get(key)
        if override:
            return str(override)
        if location is not None:
            staged = location / "public" / filename
            # In WASM this is a URL and always used; locally it is a real path,
            # so fall back when the notebook is opened outside the build dir.
            if not hasattr(staged, "exists") or staged.exists():
                return str(staged)
        return f"results/tables/{filename}"

    detection_path = resolve("detection", "checkpoint_detection.tsv")
    summary_path = resolve("summary", "checkpoint_detection_summary.json")
    carrier_path = resolve("carrier", "checkpoint_carrier_models.tsv")
    carrier_expl_path = resolve(
        "carrier_exploratory", "checkpoint_carrier_exploratory.tsv"
    )
    shift_path = resolve("shift", "checkpoint_shift_models.tsv")
    shift_expl_path = resolve(
        "shift_exploratory", "checkpoint_shift_exploratory.tsv"
    )
    return (
        carrier_expl_path,
        carrier_path,
        detection_path,
        shift_expl_path,
        shift_path,
        summary_path,
    )


@app.cell
def _(
    carrier_expl_path,
    carrier_path,
    detection_path,
    shift_expl_path,
    shift_path,
    summary_path,
):
    import json
    import sys
    from pathlib import Path

    import pandas as pd

    # Under Pyodide the six paths above are http(s) URLs, and pandas cannot open
    # those: the browser has no sockets, so urllib is not functional and
    # `read_csv` on a URL raises. `pyodide.http.open_url` performs a synchronous
    # same-origin fetch and hands back a file-like object, which is the
    # supported way to read a URL in WASM. Outside Pyodide these are ordinary
    # paths and open normally. The failure this avoids is DELAYED — the
    # `--execute` preview renders first and Pyodide throws seconds later, which
    # is why it survived a lint (ADR 0010).
    if sys.platform == "emscripten":
        from pyodide.http import open_url

        def read_source(path):
            return open_url(str(path))
    else:

        def read_source(path):
            return path

    def read_table(path):
        # float_precision="round_trip": pandas' default C parser is not
        # correctly rounded (ADR 0013 postscript), and this notebook compares
        # `detection_rate` against a floor.
        return pd.read_csv(
            read_source(path), sep="\t", float_precision="round_trip"
        )

    det = read_table(detection_path)
    carrier = read_table(carrier_path)
    carrier_expl = read_table(carrier_expl_path)
    shift = read_table(shift_path)
    shift_expl = read_table(shift_expl_path)

    if sys.platform == "emscripten":
        summary = json.load(read_source(summary_path))
    else:
        summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    return carrier, carrier_expl, det, pd, shift, shift_expl, summary


@app.cell
def _(det, summary):
    # ------------------------------------------------------------------------
    # Encodings. workflow/scripts/checkpoint_dotplot.py is THE REFERENCE for
    # every one of these, and at the default controls this notebook must
    # reproduce results/figures/checkpoint_dotplot.png exactly. Two copies of
    # the encoding is a real cost, and it is forced: `marimo export html-wasm`
    # bundles the notebook plus `public/` and nothing else, so a shared module
    # sitting beside the notebook would not be importable in the browser.
    # ------------------------------------------------------------------------

    # Site-major, tumour before immune within site. This DEVIATES from the
    # project's compartment-major CODE_ORDER because here the axis IS
    # compartment x site: it puts P3-T3's carrier contrast adjacent within each
    # site and P3-T4's shift contrast positionally parallel between them.
    CODE_ORDER = ["L", "TIME-L", "mLN", "LB", "TIME-B", "TBME", "BC"]

    # Reused verbatim from pca_landscape.py so one study has one site palette.
    SITE_COLOUR = {"lung": "#C1666B", "brain": "#4F6D7A", "lymph_node": "#D4B483"}
    SITE_LABEL = {"lung": "lung", "brain": "brain", "lymph_node": "lymph node"}

    CMAP = "PRGn"          # deliberately NOT P2-T7's RdBu_r, which means a
    S_MAX = 260.0          # SIGNED EFFECT. Two deliverables must not share a
    S_ZERO = 26.0          # colour language for two different quantities.
    DOT_EDGE = "#2B2B2B"
    HATCH_EDGE = "#8A8A8A"

    # Both read from the pipeline, never hardcoded, so neither can drift from
    # the rule it depicts.
    PREREG_FLOOR = float(summary["detected_in_aoi_fraction"])   # ADR 0016 §1
    MIN_ASSESSABLE = float(summary["min_assessable_fraction"])  # ADR 0016 §1
    VCENTER = float(summary["detection_background_multiple"])   # ADR 0007

    ALIAS = det.drop_duplicates("gene").set_index("gene")["alias"]
    CODE_INFO = det.drop_duplicates("aoi_code").set_index("aoi_code")

    # The colour-basis contract P3-T2 wrote forward to every figure over this
    # table (ADR 0016 §5). Checked here for the same reason the static script
    # checks it: so a changed string cannot silently license plotting something
    # else.
    COLOUR_BASIS_OK = str(summary["figure_colour_basis"]).startswith(
        "median_negprobe_ratio"
    )

    # The pre-registered flag, re-derived. At PREREG_FLOOR this must equal the
    # `assessable` column the pipeline wrote; the static figure asserts the same
    # thing from the other side.
    DEFAULT_MATCHES_TABLE = bool(
        ((det["detection_rate"] >= PREREG_FLOOR) == det["assessable"]).all()
    )
    return (
        ALIAS,
        CMAP,
        CODE_INFO,
        CODE_ORDER,
        COLOUR_BASIS_OK,
        DEFAULT_MATCHES_TABLE,
        DOT_EDGE,
        HATCH_EDGE,
        MIN_ASSESSABLE,
        PREREG_FLOOR,
        SITE_COLOUR,
        SITE_LABEL,
        S_MAX,
        S_ZERO,
        VCENTER,
    )


@app.cell
def _(CODE_ORDER, PREREG_FLOOR, mo):
    floor = mo.ui.slider(
        start=0.0,
        stop=1.0,
        step=0.05,
        value=PREREG_FLOOR,
        show_value=True,
        label="detection floor",
    )
    compartments = mo.ui.multiselect(
        options=CODE_ORDER, value=CODE_ORDER, label="compartments"
    )
    site_filter = mo.ui.radio(
        options=["all sites", "lung", "brain", "lymph node"],
        value="all sites",
        inline=True,
        label="site",
    )
    mo.vstack([
        mo.hstack([floor, site_filter], justify="start", gap=2, wrap=True),
        compartments,
    ])
    return compartments, floor, site_filter


@app.cell
def _(floor):
    # Snapped to the slider's own step. `start + i * step` accumulates binary
    # error -- 0.05 x 10 is 0.5000000000000001, not 0.5 -- and an unsnapped
    # value would both miss the equality test against the pre-registered floor
    # and flip any cell whose detection rate sits EXACTLY on it. No such cell
    # exists in the current table (0.5 needs an even n_aoi and none of 30, 8 or
    # 20 lands there), and this must not depend on that staying true.
    floor_value = round(float(floor.value), 2)
    return (floor_value,)


@app.cell
def _(PREREG_FLOOR, floor_value, mo):
    # The app labels itself the moment it leaves the pre-registered value, so a
    # screenshot taken at 0.75 carries its own caveat (ADR 0019, decision 3).
    floor_banner = (
        mo.callout(
            mo.md(
                f"**Pre-registered floor — {PREREG_FLOOR:g}** (ADR 0016 §1). "
                "This is the floor of record: every table, model and q-value in "
                "Phase 3 was computed here, and this view reproduces "
                "`checkpoint_dotplot.png`."
            ),
            kind="neutral",
        )
        if floor_value == PREREG_FLOOR
        else mo.callout(
            mo.md(
                f"**Sensitivity view — floor {floor_value:g}, not the "
                f"pre-registered {PREREG_FLOOR:g}.** Nothing shown at this "
                "setting is a result (ADR 0019). The floor of record is "
                f"{PREREG_FLOOR:g} (ADR 0016 §1); no model was refitted, no "
                "q-value changed, and the tables below are unaffected."
            ),
            kind="warn",
        )
    )
    floor_banner
    return


@app.cell
def _(CODE_INFO, CODE_ORDER, compartments, site_filter):
    # Column order is CODE_ORDER applied to the survivors, never rebuilt from
    # the filtered frame — filtering must hide columns, never reorder them.
    SITE_OF_CHOICE = {"lung": "lung", "brain": "brain", "lymph node": "lymph_node"}
    chosen = set(compartments.value)
    if site_filter.value != "all sites":
        wanted = SITE_OF_CHOICE[site_filter.value]
        chosen = {c for c in chosen if CODE_INFO.at[c, "site"] == wanted}
    visible_codes = [c for c in CODE_ORDER if c in chosen]
    return (visible_codes,)


@app.cell
def _():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Rectangle

    plt.rcParams["hatch.linewidth"] = 0.6
    return Rectangle, TwoSlopeNorm, np, pe, plt


@app.cell
def _(CODE_ORDER, det, floor_value):
    # `assessable` at the CURRENT floor. This is the one thing the notebook
    # derives, and it derives a display flag, not a result (ADR 0019).
    view = det.assign(assessable_at_floor=det["detection_rate"] >= floor_value)

    # The sort key counts ALL SEVEN compartments at the current floor, never the
    # visible ones: the slider is a sensitivity control and should re-rank the
    # rows, while the compartment filter is a pure view and should not. The
    # denominator is printed on the right-hand axis so the key explains itself.
    #
    # It averages nothing. Pooling a detection rate across compartments is the
    # one operation ADR 0016 and checkpoint_detection.py both forbid, and it
    # would be forbidden as a sort key too.
    clears = view.groupby("gene")["assessable_at_floor"].sum().astype(int)
    gene_order = sorted(clears.index, key=lambda g: (-clears[g], g))
    n_all_codes = len(CODE_ORDER)
    return clears, gene_order, n_all_codes, view


@app.cell
def _(
    ALIAS,
    CMAP,
    CODE_INFO,
    DOT_EDGE,
    HATCH_EDGE,
    Rectangle,
    SITE_COLOUR,
    SITE_LABEL,
    S_MAX,
    S_ZERO,
    TwoSlopeNorm,
    VCENTER,
    clears,
    floor_value,
    gene_order,
    mo,
    n_all_codes,
    np,
    pe,
    plt,
    view,
    visible_codes,
):
    # Every encoding below is checkpoint_dotplot.py's, and the reasoning for each
    # is in that script's docstring. In one line apiece:
    #   * detection only, no effect size on this canvas -- ADR 0018
    #   * colour = median_negprobe_ratio, midpoint AT the detection rule
    #   * TwoSlopeNorm is piecewise-linear, so colour distance != ratio distance
    #   * zero detection is an open RING, never a shrunken dot
    #   * below the floor is HATCHED, and the hatch is load-bearing
    #   * every count is printed, so the figure is greyscale-safe and verifiable
    shown = view[view["aoi_code"].isin(visible_codes)]
    n_comp = len(visible_codes)

    x_of = {c: i for i, c in enumerate(visible_codes)}
    y_of = {g: i for i, g in enumerate(gene_order)}

    # PIN THE STYLE. The export sandbox renders under a DARK matplotlib theme --
    # verified by decoding the figure out of a built index.html, where the
    # unbanded rows came back black and the default-coloured tick labels white.
    # checkpoint_dotplot.py runs under matplotlib's defaults and this app has to
    # reproduce it exactly, so ambient rcParams cannot be allowed to decide.
    # Resetting here rather than at import time because the theme is applied to
    # the kernel, not to this module.
    #
    # results/reports/landscape_explorer has the same black background and is
    # NOT fixed here: editing that notebook fires p1t5_landscape_explorer's
    # rerun trigger, and it is P1-T5's file.
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams["hatch.linewidth"] = 0.6

    fig = plt.figure(figsize=(13.8, 7.8), facecolor="white")
    gs = fig.add_gridspec(
        2, 3, height_ratios=[0.05, 1], width_ratios=[1, 0.30, 0.026],
        hspace=0.06, wspace=0.34,
        left=0.115, right=0.878, top=0.93, bottom=0.125,
    )
    ax_site = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0])
    ax_key = fig.add_subplot(gs[1, 1])
    ax_cbar = fig.add_subplot(gs[1, 2])
    ax_key.set_axis_off()
    for panel in (ax_site, ax, ax_key):
        panel.set_facecolor("white")

    # Guard the norm the same way the static script does: TwoSlopeNorm raises if
    # vmin >= vcenter, and a filtered subset could put every cell on one side.
    ratio = shown["median_negprobe_ratio"] if n_comp else view["median_negprobe_ratio"]
    vmin = min(float(ratio.min()), 0.95 * VCENTER)
    vmax = max(float(ratio.max()), 1.05 * VCENTER)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=VCENTER, vmax=vmax)

    site_of = {c: CODE_INFO.at[c, "site"] for c in visible_codes}
    for i, code in enumerate(visible_codes):
        ax_site.add_patch(
            Rectangle((i, 0), 1, 1, facecolor=SITE_COLOUR[site_of[code]], lw=0)
        )
    ax_site.set_xlim(0, max(n_comp, 1))
    ax_site.set_ylim(0, 1)
    ax_site.set_xticks([])
    ax_site.set_yticks([])
    for spine in ax_site.spines.values():
        spine.set_visible(False)

    # site labels + the black rule between sites
    blocks = []
    block_start = 0
    for i, code in enumerate(visible_codes):
        if i + 1 == n_comp or site_of[visible_codes[i + 1]] != site_of[code]:
            blocks.append((site_of[code], block_start, i + 1))
            block_start = i + 1
    site_n = (
        shown.drop_duplicates("aoi_code").groupby("site")["n_aoi"].sum()
        if n_comp
        else {}
    )
    for site, a, b in blocks:
        ax_site.text(
            (a + b) / 2, 1.5, f"{SITE_LABEL[site]}  (n = {site_n[site]})",
            ha="center", va="bottom", fontsize=9.5,
            color=SITE_COLOUR[site], fontweight="bold",
        )
        if b < n_comp:
            ax_site.axvline(b, color="white", lw=2.0)
            ax.axvline(b - 0.5, color="black", lw=1.6, zorder=2)

    # row bands, so the eye tracks a gene across the columns
    for gene in gene_order:
        if y_of[gene] % 2 == 0:
            ax.add_patch(
                Rectangle((-0.5, y_of[gene] - 0.5), max(n_comp, 1), 1,
                          facecolor="#F7F7F7", lw=0, zorder=0)
            )

    hatched = shown.loc[~shown["assessable_at_floor"], ["gene", "aoi_code"]]
    rings = shown.loc[shown["n_detected"] == 0, ["gene", "aoi_code"]]

    for row in hatched.itertuples(index=False):
        ax.add_patch(
            Rectangle(
                (x_of[row.aoi_code] - 0.5, y_of[row.gene] - 0.5), 1, 1,
                fill=False, hatch="///", edgecolor=HATCH_EDGE, lw=0.0,
                alpha=0.45, zorder=1,
            )
        )

    halo = [pe.withStroke(linewidth=2.2, foreground="white")]

    zx = [x_of[r.aoi_code] for r in rings.itertuples(index=False)]
    zy = [y_of[r.gene] for r in rings.itertuples(index=False)]
    if zx:
        ring = ax.scatter(
            zx, zy, s=S_ZERO, marker="o", facecolors="none",
            edgecolors="#3A3A3A", linewidths=1.0, zorder=3,
        )
        ring.set_path_effects(halo)

    drawn = shown.loc[shown["n_detected"] > 0]
    sc = ax.scatter(
        [x_of[c] for c in drawn["aoi_code"]],
        [y_of[g] for g in drawn["gene"]],
        s=S_MAX * drawn["detection_rate"].to_numpy(),
        c=drawn["median_negprobe_ratio"].to_numpy(),
        cmap=CMAP, norm=norm, edgecolors=DOT_EDGE,
        linewidths=np.where(drawn["assessable_at_floor"].to_numpy(), 1.30, 0.55),
        zorder=3,
    )
    sc.set_path_effects(halo)

    for row in shown.itertuples(index=False):
        label = ax.text(
            x_of[row.aoi_code], y_of[row.gene] + 0.33,
            f"{int(row.n_detected)}/{int(row.n_aoi)}",
            ha="center", va="center", fontsize=5.4, color="#4A4A4A", zorder=5,
        )
        label.set_path_effects(halo)

    ax.set_xlim(-0.5, max(n_comp, 1) - 0.5)
    ax.set_ylim(len(gene_order) - 0.5, -0.5)
    ax.set_xticks(range(n_comp))
    ax.set_xticklabels(
        [
            f"{c}{'' if CODE_INFO.at[c, 'modelled'] else ' †'}"
            f"\nn = {int(CODE_INFO.at[c, 'n_aoi'])}"
            for c in visible_codes
        ],
        fontsize=9,
    )
    ax.set_yticks(range(len(gene_order)))
    ax.set_yticklabels(
        [f"{g}  ({ALIAS[g]})" for g in gene_order], fontsize=9.5
    )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_color("#BBBBBB")

    # right-hand axis: the sort key, printed. Counts all seven compartments at
    # the current floor, which is why the denominator is spelled out.
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(range(len(gene_order)))
    ax2.set_yticklabels(
        [f"{clears[g]}/{n_all_codes}" for g in gene_order], fontsize=8
    )
    ax2.text(
        1.004, 1.035, f"clears the\n{floor_value:g} floor\n(of {n_all_codes})",
        transform=ax.transAxes, fontsize=7.4, ha="left", va="bottom",
        color="#444444", linespacing=1.2,
    )
    ax2.tick_params(length=0, pad=2)
    for spine in ax2.spines.values():
        spine.set_visible(False)

    cb = fig.colorbar(sc, cax=ax_cbar)
    cb.ax.tick_params(labelsize=7.4)
    cb.ax.plot(
        [0, 1], [VCENTER, VCENTER], color="black", lw=1.4,
        transform=cb.ax.get_yaxis_transform(), clip_on=False,
    )
    cb.ax.plot(
        [0, 1], [1.0, 1.0], color="0.35", lw=0.9, ls=":",
        transform=cb.ax.get_yaxis_transform(),
    )
    # The threshold and the background level are named IN THE TICK LABELS. Side
    # annotations were tried in P3-T5 and collided with the key panel; a tick
    # label cannot drift away from the value it names.
    cb_ticks = [
        t for t in sorted({round(vmin, 2), 1.0, VCENTER, 3.0, 4.0, 5.0,
                           round(vmax, 2)})
        if vmin <= t <= vmax
    ]
    cb.set_ticks(cb_ticks)
    cb.set_ticklabels([
        f"{t:g} ← threshold" if t == VCENTER
        else (f"{t:g} ← background" if t == 1.0 else f"{t:g}")
        for t in cb_ticks
    ])
    cb.set_label("median  q3 / NegProbe-WTX  (ADR 0016 §5)", fontsize=8.0,
                 labelpad=4)

    # ------------------------------------------------------------------- key
    # Drawn with the SAME sizes as the field, so the key cannot drift from what
    # it describes.
    ax_key.set_xlim(0, 1)
    ax_key.set_ylim(0, 1)
    key_rates = [0.25, 0.50, 0.75, 1.00]
    key_ys = np.linspace(0.80, 0.56, len(key_rates))
    ax_key.scatter(
        [0.18] * len(key_rates), key_ys, s=[S_MAX * r for r in key_rates],
        facecolor="#DDDDDD", edgecolors=DOT_EDGE, linewidths=0.55,
    )
    for rate, y in zip(key_rates, key_ys):
        ax_key.text(0.36, y, f"{rate:.0%} of AOIs", fontsize=7.5, va="center")
    ax_key.text(0.02, 0.89, "dot area = detection rate", fontsize=8,
                fontweight="bold")
    ax_key.scatter([0.18], [0.44], s=S_ZERO, marker="o", facecolors="none",
                   edgecolors="#3A3A3A", linewidths=1.0)
    ax_key.text(0.36, 0.44, "detected in 0 AOIs", fontsize=7.5, va="center")
    ax_key.add_patch(
        Rectangle((0.10, 0.26), 0.16, 0.10, fill=False, hatch="///",
                  edgecolor=HATCH_EDGE, lw=0.0, alpha=0.45)
    )
    ax_key.text(
        0.36, 0.31,
        f"below the {floor_value:g} floor —\n\"not assessable\" there",
        fontsize=7.5, va="center",
    )
    ax_key.scatter([0.18], [0.14], s=S_MAX * 0.6, facecolor="#DDDDDD",
                   edgecolors=DOT_EDGE, linewidths=1.30)
    ax_key.text(0.36, 0.14, "bold edge = clears\nthe floor", fontsize=7.5,
                va="center")

    # No suptitle: in a web page the caption belongs in markdown, where it
    # cannot clip. P2-T7 and P3-T5 both lost title text off both canvas edges.
    plot_output = (
        fig
        if n_comp
        else mo.callout(
            mo.md(
                "**No compartments selected.** The compartment multiselect and "
                "the site filter intersect; widen one of them."
            ),
            kind="warn",
        )
    )
    plot_output
    return hatched, rings, shown


@app.cell
def _(
    COLOUR_BASIS_OK,
    DEFAULT_MATCHES_TABLE,
    MIN_ASSESSABLE,
    PREREG_FLOOR,
    floor_value,
    hatched,
    mo,
    rings,
    shown,
    visible_codes,
):
    # A running count, so the headline number is checkable against P3-T5's log
    # rather than squinted at. At the default this reads 40 / 23 / 9.
    n_cells = len(shown)
    n_hatched = len(hatched)
    n_rings = len(rings)
    # The panel-level floor, also pre-registered and also read from the summary
    # rather than written down here (ADR 0016 §1).
    assessable_codes = [
        code
        for code in visible_codes
        if shown.loc[shown["aoi_code"] == code, "assessable_at_floor"].mean()
        >= MIN_ASSESSABLE
    ]
    assessable_text = (
        ", ".join("`" + code + "`" for code in assessable_codes)
        or "no compartment"
    )
    mo.md(
        f"""
        **At floor {floor_value:g}, over {len(visible_codes)} compartment(s):
        {n_hatched} of {n_cells} cells are not assessable**, {n_cells - n_hatched}
        clear the floor, and {n_rings} are detected in zero AOIs.
        The panel is assessable (≥ {MIN_ASSESSABLE:g} of its nine genes) in:
        **{assessable_text}**.

        *Self-checks:* the `assessable` column reproduces
        `detection_rate ≥ {PREREG_FLOOR:g}` — **{DEFAULT_MATCHES_TABLE}**;
        the colour-basis contract in `checkpoint_detection_summary.json` still
        declares `median_negprobe_ratio` — **{COLOUR_BASIS_OK}**. Both must be
        `True`; they are the same assertions `p3t5_checkpoint_dotplot` makes.

        † detection only: `mLN`, `TBME` and `BC` each sit wholly within one DSP
        run, so batch is inseparable from biology there (ADR 0009 §3) and they
        are audited and plotted but never modelled (ADR 0016 §2). Not lesser
        data — a property of the study design.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### The detection table

    Every value on the plot above, searchable. `detection_rate` is the dot area,
    `median_negprobe_ratio` the colour, `assessable` the pre-registered
    verdict at 0.5 — **the column is fixed at the floor of record and does not
    move with the slider**; `assessable_at_floor` is what the plot draws.

    `detection_rate_qc_clean` excludes the five QC-flagged AOIs, which are
    retained and scored (flag-don't-drop). Excluding them flips no
    assessability verdict, so `TBME`'s deficit is the compartment and not those
    AOIs.
    """)
    return


@app.cell
def _(mo, shown):
    detection_columns = [
        "gene", "alias", "aoi_code", "site", "compartment", "n_aoi",
        "n_detected", "detection_rate", "detection_rate_qc_clean",
        "median_negprobe_ratio", "assessable", "assessable_at_floor",
        "modelled",
    ]
    detection_view = shown[detection_columns].copy()
    for column in [
        "detection_rate", "detection_rate_qc_clean", "median_negprobe_ratio"
    ]:
        detection_view[column] = detection_view[column].round(4)
    mo.ui.table(
        detection_view,
        page_size=12,
        selection=None,
        show_search=True,
        show_column_summaries="stats",
        max_height=420,
        label=f"{len(detection_view)} gene × compartment cells",
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---

    ## Model estimates — a different quantity, and deliberately a different canvas

    Everything above is **detection**: is the gene above background at all.
    Everything below is **expression**, `log2(q3 + 1)`, estimated by
    `expression ~ … + (1|patient_id)` in `lme4`/`lmerTest`.

    **They move OPPOSITE ways under the same background gradient, and conflating
    them inverts the argument.** Detection is `q3 > 2 × negprobe`, so higher
    background raises the bar and a gene looks *depleted*; expression is
    `log2(q3 + 1)`, where background adds to signal, so higher background makes
    a gene look *enriched*. This project made that exact conflation once and
    reasoned from it (ADR 0018), which is why no effect size appears on the dot
    plot and no dot appears in these tables.

    The random intercept is mandatory on **non-independence** grounds, not
    variance share (ADR 0009): 33 of 42 subjects contribute more than one AOI.
    Singular fits are surfaced in `is_singular` and are not a reason to drop it.

    **`TIME-B` n = 8**, and the lung-vs-brain immune contrast detects roughly
    1.1–1.3 SD at 80% power (P0-T8). Every P3-T4 shift estimate is far below
    that floor, so those nulls are **uninformative, not negative**.

    A `TIME` AOI is the **PanCK-negative segment** of an ROI sited in a
    CD45-rich region — not a CD45-sorted population (Q2). The compartment label
    is not a cell-type label.
    """)
    return


@app.cell
def _(pd):
    PRIMARY_COLUMNS = [
        "stratum", "gene", "estimate", "ci_low", "ci_high", "p_raw", "q_bh",
        "n_obs", "n_patient", "detection_rate_reference",
        "detection_rate_test", "is_singular",
    ]
    SENSITIVITY_COLUMNS = [
        "stratum", "gene", "model", "estimate", "ci_low", "ci_high", "p_raw",
        "delta_vs_primary", "n_obs", "n_patient", "is_singular",
    ]
    EXCLUDED_COLUMNS = [
        "stratum", "gene", "status", "detection_rate_reference",
        "detection_rate_test", "assessable_reference", "assessable_test",
        "note",
    ]
    # No q_bh, by construction: the exploratory table carries none (ADR 0018).
    EXPLORATORY_COLUMNS = [
        "stratum", "gene", "estimate", "ci_low", "ci_high", "p_raw",
        "delta_vs_primary", "gradient_direction", "n_obs", "n_patient",
        "detection_rate_reference", "detection_rate_test", "is_singular",
    ]

    def model_view(frame, columns):
        """Curate and round. Selects only columns that exist, so a table
        without `q_bh` simply does not show one rather than raising."""
        present = [c for c in columns if c in frame.columns]
        out = frame[present].copy()
        for column in present:
            if pd.api.types.is_float_dtype(out[column]):
                out[column] = out[column].round(4)
        return out.reset_index(drop=True)

    # checkpoint_carrier_models.tsv and checkpoint_shift_models.tsv each hold
    # three kinds of row in one file: the primary fits, the two pre-registered
    # sensitivity refits, and a placeholder per gene the restriction excluded.
    # `model` is empty on that last kind, so `status` is what separates them.
    def primary_rows(frame):
        return frame[frame["model"] == "primary"]

    def sensitivity_rows(frame):
        return frame[
            frame["model"].isin(["batch_sensitivity", "background_sensitivity"])
        ]

    def excluded_rows(frame):
        return frame[frame["status"] == "not_assessable"]

    return (
        EXCLUDED_COLUMNS,
        EXPLORATORY_COLUMNS,
        PRIMARY_COLUMNS,
        SENSITIVITY_COLUMNS,
        excluded_rows,
        model_view,
        primary_rows,
        sensitivity_rows,
    )


@app.cell
def _(
    EXCLUDED_COLUMNS,
    EXPLORATORY_COLUMNS,
    PRIMARY_COLUMNS,
    SENSITIVITY_COLUMNS,
    carrier,
    carrier_expl,
    excluded_rows,
    mo,
    model_view,
    primary_rows,
    sensitivity_rows,
    shift,
    shift_expl,
):
    def model_table(frame, columns, label):
        return mo.ui.table(
            model_view(frame, columns),
            page_size=10,
            selection=None,
            show_search=True,
            max_height=340,
            label=label,
        )

    def pair(frame_carrier, frame_shift, columns, note):
        return mo.vstack([
            note,
            mo.md(
                "**P3-T3 — the carrier contrast.** `immune − tumour`, within "
                "site. `L` vs `TIME-L` in lung, `LB` vs `TIME-B` in brain; "
                "never pooled, because `compartment` is degenerate across "
                "sites."
            ),
            model_table(frame_carrier, columns, "carrier"),
            mo.md(
                "**P3-T4 — the shift contrast.** `brain − lung`, within "
                "compartment. A tumour-compartment shift and an "
                "immune-compartment shift mean different things "
                "therapeutically, so they are never pooled either."
            ),
            model_table(frame_shift, columns, "shift"),
        ])

    primary_note = mo.callout(
        mo.md(
            "**The primary tables — the only Phase 3 output with a q-value.** A "
            "gene enters a primary fit only if it clears the detection floor in "
            "**both** groups compared (ADR 0016 §3, upheld as Option C in "
            "ADR 0018), which admits 2 of 9 genes per site here and 4 of 9 / 2 "
            "of 9 for the shift. `CD274` (PD-L1) is enriched in the lung immune "
            "compartment, **+0.612 SD [+0.361, +0.863], q = 0.0001**, and "
            "survives both pre-registered sensitivities."
        ),
        kind="success",
    )
    sensitivity_note = mo.callout(
        mo.md(
            "**Sensitivities, reported not adjudicated** (ADR 0012's shape). "
            "`batch_sensitivity` adds `dsp_run`; `background_sensitivity` adds "
            "`negprobe_log2`, which ADR 0008 point 5 requires. "
            "`delta_vs_primary` is the shift each induces. The background "
            "adjustment moves the genes the primary *excludes* by a mean of "
            "−0.155 and the ones it *admits* by +0.010 — which is the "
            "pre-registration working."
        ),
        kind="info",
    )
    excluded_note = mo.callout(
        mo.md(
            "**The exclusions are themselves the finding.** These rows were not "
            "fitted because the gene is below the floor in at least one group — "
            "\"not assessable in that compartment\", **never** \"lower in that "
            "compartment\" (ADR 0008 point 4). Brain background is *higher*, so "
            "a near-background gene reads as depleted in brain artefactually, "
            "which is precisely the headline this phase is shaped to produce and "
            "must not."
        ),
        kind="neutral",
    )
    exploratory_note = mo.callout(
        mo.md(
            "**No claim may rest on this table** (ADR 0018). It carries the "
            "genes the primary restriction excludes: **no FDR, no q-value, "
            "unadjusted**, declared exploratory-within-exploratory. The large "
            "apparent immune enrichments here — `CTLA4` +1.425, `VSIR` +1.336, "
            "`TIGIT` +1.036, `HAVCR2` +0.963 in lung — are every one of them a "
            "gene detected in 0–9 of 30 tumour AOIs. The detection counts are "
            "reportable; the estimates are not."
        ),
        kind="danger",
    )

    mo.ui.tabs({
        "primary": pair(
            primary_rows(carrier), primary_rows(shift),
            PRIMARY_COLUMNS, primary_note,
        ),
        "sensitivities": pair(
            sensitivity_rows(carrier), sensitivity_rows(shift),
            SENSITIVITY_COLUMNS, sensitivity_note,
        ),
        "not assessable": pair(
            excluded_rows(carrier), excluded_rows(shift),
            EXCLUDED_COLUMNS, excluded_note,
        ),
        "exploratory": pair(
            carrier_expl, shift_expl, EXPLORATORY_COLUMNS, exploratory_note
        ),
    })
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---

    **Provenance.** Every file this app reads is a declared input of
    `p3t6_checkpoint_explorer` and is staged into `public/` beside the export
    (ADR 0010); an undeclared read would be a silent provenance hole
    (PROJECT_PLAN §A.5). Aim A4 is **exploratory** (ADR 0008): no result here
    may be a headline claim.
    """)
    return


if __name__ == "__main__":
    app.run()
