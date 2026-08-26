# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo",
#     "pandas",
#     "matplotlib",
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
    # P1-T5 — landscape explorer

    **Tier:** app (`notebooks/apps/`) — **presentation logic only** (`CLAUDE.md`
    hard constraint 4). Every number here was produced by a Snakemake rule; this
    notebook reads `p1t1_pca_landscape`'s outputs and computes nothing of its
    own. Put any PC on either axis, recolour, and read the metadata table.

    **What the ordination is.** PCA on the top 2000 highly variable genes of the
    120 AOIs, zero-centred but **not** unit-scaled — scaling would give a gene
    detected in three AOIs the same leverage as one detected in all 120, which
    P0-T5 showed would convert a detection artefact into a principal component.

    **What to look for**, and it is the Phase 1 result: colour by
    **compartment** and the structure is obvious; colour by **patient** and it
    dissolves. P1-T3 put the number on that — compartment 41% of variance
    against patient 22% (PVCA) — and P1-T4 confirmed it from the other side
    (compartment ARI 0.741, patient 0.003). PROJECT_PLAN §6 expected patients to
    cluster together. They do not.

    **QC-flagged AOIs are ringed, never dropped** (`qc.flag_only`). Whether they
    sit apart from their compartment is the question flagging exists to raise.

    `TIME-B` **n = 8** — the binding constraint on the whole project.
    """)
    return


@app.cell
def _(mo):
    # Parameters (PROJECT_PLAN §A.5), resolved for three different homes.
    #
    # The exporting rule stages this notebook next to a `public/` directory
    # holding the three pipeline outputs, because `marimo export html-wasm`
    # copies `public/` into the export. `mo.notebook_location()` then resolves
    # to that directory locally and to the served URL in the browser, so the
    # WASM build fetches its data over HTTP instead of needing a filesystem.
    #
    # That is what keeps the tables interactive. The alternative — bundling the
    # frames via `[tool.marimo.runtime] cache_cells` — cannot work here: a
    # `mo.ui.table` holds a locally defined `PandasTableManager` that cannot be
    # pickled, so its cell output is cached as an `UnhashableStub`, and an
    # unserializable output makes marimo re-run the cell's ancestors live. In
    # the browser that ancestor is the loader, which has no files to read.
    args = mo.cli_args()
    location = mo.notebook_location()

    def resolve(key, filename, fallback):
        override = args.get(key)
        if override:
            return str(override)
        if location is not None:
            staged = location / "public" / filename
            # In WASM this is a URL and always used; locally it is a real path,
            # so fall back when the notebook is opened outside the build dir.
            if not hasattr(staged, "exists") or staged.exists():
                return str(staged)
        return fallback

    coords_path = resolve(
        "coords", "pca_coords.tsv", "results/tables/pca_coords.tsv"
    )
    umap_path = resolve("umap", "umap_coords.tsv", "results/tables/umap_coords.tsv")
    summary_path = resolve(
        "summary", "pca_summary.json", "results/tables/pca_summary.json"
    )
    return coords_path, summary_path, umap_path


@app.cell
def _(coords_path, summary_path, umap_path):
    import json
    import sys
    from pathlib import Path

    import pandas as pd

    # Under Pyodide the three paths above are http(s) URLs, and pandas cannot
    # open those: the browser has no sockets, so urllib is not functional and
    # `read_csv` on a URL raises. `pyodide.http.open_url` performs a synchronous
    # same-origin fetch and hands back a file-like object, which is the
    # supported way to read a URL in WASM. Outside Pyodide these are ordinary
    # paths and open normally.
    if sys.platform == "emscripten":
        from pyodide.http import open_url

        def read_source(path):
            return open_url(str(path))
    else:

        def read_source(path):
            return path

    coords = pd.read_csv(read_source(coords_path), sep="\t", index_col="aoi_label")
    umap = pd.read_csv(read_source(umap_path), sep="\t", index_col="aoi_label")
    if sys.platform == "emscripten":
        summary = json.load(read_source(summary_path))
    else:
        summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    landscape = coords.join(umap, how="left")
    return landscape, summary


@app.cell
def _(landscape, summary):
    COMPARTMENT_COLOUR = {
        "tumour": "#B0763F",
        "immune": "#3E7CB1",
        "glial_stroma": "#5B8C5A",
        "normal_control": "#8C8C8C",
    }
    SITE_COLOUR = {"lung": "#C1666B", "brain": "#4F6D7A", "lymph_node": "#D4B483"}
    CODE_COLOUR = {
        "L": "#B0763F",
        "LB": "#8C5A2B",
        "mLN": "#D4B483",
        "TIME-L": "#3E7CB1",
        "TIME-B": "#2A5A87",
        "TBME": "#5B8C5A",
        "BC": "#8C8C8C",
    }
    RUN_COLOUR = {r: c for r, c in zip(
        sorted(landscape["dsp_run"].astype(str).unique()), ["#7B4B94", "#C1666B"]
    )}
    PALETTE = {
        "compartment": COMPARTMENT_COLOUR,
        "site": SITE_COLOUR,
        "aoi_code": CODE_COLOUR,
        "dsp_run": RUN_COLOUR,
    }

    pc_columns = [c for c in landscape.columns if c.startswith("PC")]
    variance_ratio = summary["variance_ratio"]

    def hue_colours(values, cmap):
        """One colour per value for a high-cardinality factor.

        Built here rather than looped over in the plot cell: a loop variable
        defined inside only one branch of that cell becomes a conditionally
        defined cell output, which marimo cannot cache and which would be a
        genuine bug the moment anything downstream referenced it.
        """
        levels = sorted(set(values))
        lookup = {
            level: cmap(i / max(len(levels) - 1, 1))
            for i, level in enumerate(levels)
        }
        return [lookup[v] for v in values], len(levels)

    def axis_label(column):
        """PC labels carry their % variance; UMAP axes have no such quantity."""
        if column.startswith("PC"):
            index = int(column[2:]) - 1
            if index < len(variance_ratio):
                return f"{column} ({100 * variance_ratio[index]:.1f}%)"
        return column

    return PALETTE, axis_label, hue_colours, pc_columns


@app.cell
def _(mo, pc_columns):
    embedding = mo.ui.dropdown(
        options=["PCA", "UMAP"], value="PCA", label="embedding"
    )
    x_axis = mo.ui.dropdown(options=pc_columns, value="PC1", label="x axis")
    y_axis = mo.ui.dropdown(options=pc_columns, value="PC2", label="y axis")
    colour_by = mo.ui.dropdown(
        options=["compartment", "site", "patient_id", "aoi_code", "dsp_run"],
        value="compartment",
        label="colour by",
    )
    join_patients = mo.ui.switch(value=False, label="join each patient's AOIs")
    mo.hstack(
        [embedding, x_axis, y_axis, colour_by, join_patients],
        justify="start",
        gap=1,
        wrap=True,
    )
    return colour_by, embedding, join_patients, x_axis, y_axis


@app.cell
def _(embedding, mo, x_axis, y_axis):
    # UMAP has exactly two axes, so the PC dropdowns do not apply to it. Say so
    # rather than leaving them looking live but ignored.
    axis_note = (
        mo.md(
            f"Showing **UMAP1 vs UMAP2**. The axis dropdowns "
            f"(`{x_axis.value}` / `{y_axis.value}`) apply to the PCA only."
        )
        if embedding.value == "UMAP"
        else mo.md(
            f"Showing **{x_axis.value} vs {y_axis.value}** of the PCA. "
            "Switch the embedding dropdown for the UMAP."
        )
    )
    axis_note
    return


@app.cell
def _():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    return np, plt


@app.cell
def _(
    PALETTE,
    axis_label,
    colour_by,
    embedding,
    hue_colours,
    join_patients,
    landscape,
    np,
    plt,
    x_axis,
    y_axis,
):
    if embedding.value == "UMAP":
        xcol, ycol = "UMAP1", "UMAP2"
    else:
        xcol, ycol = x_axis.value, y_axis.value

    xs = landscape[xcol].to_numpy()
    ys = landscape[ycol].to_numpy()
    key = colour_by.value
    levels = landscape[key].astype(str)
    flagged = landscape["qc_flag"].astype(str).isin(["True", "true"]).to_numpy()

    fig, ax = plt.subplots(figsize=(8.2, 6.2))

    # Computed unconditionally: a name defined in only one branch becomes a
    # conditionally defined cell output, which marimo cannot cache.
    point_colours, n_levels = hue_colours(levels, plt.get_cmap("hsv"))

    if key in PALETTE:
        for level, colour in PALETTE[key].items():
            mask = (levels == level).to_numpy()
            if not mask.any():
                continue
            ax.scatter(xs[mask], ys[mask], s=52, c=colour, edgecolor="black",
                       linewidth=0.5, alpha=0.88, label=level)
        ax.legend(frameon=False, fontsize=8, loc="best")
    else:
        # patient_id: 42 levels, so a categorical legend carries no information.
        # Colour by hue and let the question be whether neighbours share one.
        ax.scatter(xs, ys, s=52, color=point_colours, edgecolor="black",
                   linewidth=0.5, alpha=0.9)
        ax.set_title(f"{n_levels} subjects, no legend", fontsize=9,
                     loc="right", color="#666666")

    patients = landscape["patient_id"].astype(str)
    if join_patients.value:
        for _subject in patients.unique():
            _mask = (patients == _subject).to_numpy()
            if _mask.sum() < 2:
                continue
            _order = np.argsort(xs[_mask])
            ax.plot(xs[_mask][_order], ys[_mask][_order], color="#555555",
                    lw=0.7, alpha=0.45, zorder=1)

    ax.scatter(xs[flagged], ys[flagged], s=190, facecolors="none",
               edgecolors="crimson", linewidth=1.4, zorder=5)

    ax.set_xlabel(axis_label(xcol), fontsize=9)
    ax.set_ylabel(axis_label(ycol), fontsize=9)
    ax.grid(alpha=0.22, lw=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### AOI metadata

    Sortable and filterable — narrow to a compartment, a patient or a DSP run to
    see which points above they are. The selected axes are included as columns.
    `qc_flag` marks the AOIs ringed in red; they are retained, never dropped.
    """)
    return


@app.cell
def _(landscape, mo, x_axis, y_axis):
    metadata_columns = [
        "aoi_code", "patient_id", "site", "compartment", "dsp_run", "qc_flag",
    ]
    metadata = landscape[metadata_columns].copy()
    for axis_column in dict.fromkeys([x_axis.value, y_axis.value, "UMAP1", "UMAP2"]):
        metadata[axis_column] = landscape[axis_column].round(3)
    metadata = metadata.reset_index()
    aoi_table = mo.ui.table(
        metadata,
        page_size=12,
        selection=None,
        show_search=True,
        show_column_summaries="stats",
        max_height=420,
        label=f"{len(metadata)} AOIs",
    )
    aoi_table
    return


@app.cell
def _(landscape, mo):
    counts = (
        landscape.groupby(["aoi_code", "site", "compartment"], observed=True)
        .size()
        .rename("n_aoi")
        .reset_index()
        .sort_values("n_aoi", ascending=False)
    )
    design_table = mo.ui.table(
        counts, page_size=8, selection=None, show_search=False
    )
    mo.vstack([
        mo.md("### Design table, as loaded — compare against PROJECT_PLAN §2.1"),
        design_table,
        mo.md(
            f"**{len(landscape)} AOIs**, "
            f"{landscape['patient_id'].nunique()} subjects, "
            f"`TIME-B` n = {int((landscape['aoi_code'] == 'TIME-B').sum())}."
        ),
    ])
    return


if __name__ == "__main__":
    app.run()
