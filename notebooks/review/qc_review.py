# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo",
#     "pandas",
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
    # P0-T5 — QC threshold review

    **Tier:** review (`notebooks/review/`) — not in the Snakemake DAG. A human
    runs this at the P0-T5 gate. It reads `results/tables/qc_metrics.tsv` and
    writes nothing.

    **Choosing a QC threshold is a "stop and ask" decision** (`CLAUDE.md`). This
    notebook is where the choice is *made*; `config/config.yaml` is where it
    *lives*, and the R— sorry, the pipeline rule reads config, never this file.
    The chosen values need an ADR citing this notebook by path.

    **Policy is flag, never drop** (`qc.flag_only: true`). Nothing below is
    removed from the dataset — flagged AOIs go to `qc_excluded.tsv` with a
    reason and stay in. So these sliders decide what gets *marked*, not what
    gets deleted.

    ### Watch `TIME-B`

    It has **n = 8**, and it is the binding constraint on the whole project.
    Losing 1–2 changes what Phase 2 can support.
    """)
    return


@app.cell
def _(mo):
    # Parameters (PROJECT_PLAN §A.5).
    args = mo.cli_args()
    metrics_path = args.get("metrics", "results/tables/qc_metrics.tsv")
    return (metrics_path,)


@app.cell
def _(metrics_path):
    import pandas as pd

    qc = pd.read_csv(metrics_path, sep="\t")
    CODE_ORDER = ["L", "LB", "mLN", "TIME-L", "TIME-B", "TBME", "BC"]
    return CODE_ORDER, pd, qc


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. Detection is background-relative, and the background is not flat

    The matrix has **no zeros** (minimum 2.12), so "detected" cannot mean
    "non-zero" — it means "above this AOI's own `NegProbe-WTX` level, by some
    multiple". That multiple is itself a threshold decision, so it gets a
    slider too.
    """)
    return


@app.cell
def _(mo):
    background_multiple = mo.ui.dropdown(
        options={f"{m:g}x": m for m in [1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0]},
        value="2x",
        label="Detection threshold: gene value must exceed this multiple of NegProbe",
    )
    background_multiple
    return (background_multiple,)


@app.cell
def _(CODE_ORDER, background_multiple, mo, qc):
    detection_column = f"detection_at_{background_multiple.value:g}x"
    by_code = (
        qc.groupby("aoi_code")[detection_column]
        .agg(["size", "median", "min", "max"])
        .reindex(CODE_ORDER)
    )
    detection_rows = "\n".join(
        f"| `{code}` | {int(r['size'])} | {r['median'] * 100:.1f}% | "
        f"{r['min'] * 100:.1f}% | {r['max'] * 100:.1f}% |"
        for code, r in by_code.iterrows()
    )
    mo.md(f"""
    **Detection rate at >{background_multiple.value:g}× background, by compartment**

    | Compartment | n | median | min | max |
    |---|---|---|---|---|
    {detection_rows}

    `TBME` and `BC` detect far less than the tumour cores at any multiple. That
    is the tissue plus its higher background, not bad data — and it is exactly
    why a single global cutoff flags the glial compartment first.
    """)
    return (detection_column,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. The thresholds

    Move these and watch the table at the bottom. `min_library_size` is on
    **raw sequencing reads from the DCC header**, not the matrix column sum —
    the matrix is already Q3-normalised, so its sums are near-constant by
    construction and measure nothing.
    """)
    return


@app.cell
def _(mo):
    min_detection = mo.ui.slider(
        start=0.0, stop=50.0, step=0.5, value=0.0,
        label="min gene detection rate (%)", show_value=True, full_width=True,
    )
    min_reads = mo.ui.slider(
        start=0.0, stop=5.0, step=0.05, value=0.0,
        label="min raw sequencing reads (millions)", show_value=True, full_width=True,
    )
    min_saturation = mo.ui.slider(
        start=0.0, stop=0.95, step=0.01, value=0.0,
        label="min sequencing saturation", show_value=True, full_width=True,
    )
    mo.vstack([min_detection, min_reads, min_saturation])
    return min_detection, min_reads, min_saturation


@app.cell
def _(detection_column, min_detection, min_reads, min_saturation, qc):
    flagged_mask = (
        (qc[detection_column] * 100 < min_detection.value)
        | (qc["raw"] / 1e6 < min_reads.value)
        | (qc["sequencing_saturation"] < min_saturation.value)
    )
    flagged = qc[flagged_mask]
    surviving = qc[~flagged_mask]
    return flagged, surviving


@app.cell
def _(CODE_ORDER, mo, qc, surviving):
    full_counts = qc.groupby("aoi_code").size().reindex(CODE_ORDER, fill_value=0)
    kept_counts = surviving.groupby("aoi_code").size().reindex(CODE_ORDER, fill_value=0)
    count_rows = "\n".join(
        f"| `{code}` | {full_counts[code]} | **{kept_counts[code]}** | "
        f"{kept_counts[code] - full_counts[code]:+d} |"
        for code in CODE_ORDER
    )
    time_b_kept = int(kept_counts["TIME-B"])
    if time_b_kept == 8:
        time_b_note = "`TIME-B` still n = 8. The binding constraint holds at its stated value."
    elif time_b_kept >= 6:
        time_b_note = (
            f"**`TIME-B` down to n = {time_b_kept}.** Every claim about brain immune "
            "contexture must state this inline (hard constraint 8), and Phase 2's "
            "power drops with it — check P0-T8 before accepting."
        )
    else:
        time_b_note = (
            f"**`TIME-B` down to n = {time_b_kept}. This is too few.** At this "
            "threshold the brain immune comparison is not worth running; either "
            "relax the cutoff or accept that Aim A2 is out of reach."
        )
    mo.md(f"""
    ### What the thresholds would flag

    **{len(qc) - len(surviving)} of {len(qc)} AOIs flagged** (still retained —
    `flag_only` is true).

    | Compartment | Full | Unflagged | Δ |
    |---|---|---|---|
    {count_rows}

    {time_b_note}
    """)
    return


@app.cell
def _(detection_column, flagged, mo, pd):
    mo.md("### The flagged AOIs")
    columns = ["aoi_label", "aoi_code", "patient_id", "raw", "sequencing_saturation", detection_column]
    table = (
        flagged[columns].sort_values(detection_column)
        if len(flagged)
        else pd.DataFrame({"note": ["nothing flagged at these thresholds"]})
    )
    table
    return


@app.cell
def _(mo, qc):
    mo.md("### Full metrics, sortable")
    qc[
        ["aoi_label", "aoi_code", "patient_id", "dsp_run", "raw", "aligned",
         "sequencing_saturation", "negprobe", "detection_at_2x"]
    ].sort_values("detection_at_2x")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---

    ## 3. Writing the decision down

    Once the values are chosen, they go into `config/config.yaml`:

    ```yaml
    qc:
      detection_background_multiple: <multiple>
      min_gene_detection_rate: <fraction, not percent>
      min_library_size: <raw reads>
      min_sequencing_saturation: <fraction>
    ```

    then `snakemake --use-conda` re-runs `p0t5_qc_metrics`, and the choice gets
    an ADR in `docs/decisions/` **citing this notebook by path** — that citation
    is P0-T5's acceptance criterion.

    A threshold without an ADR is a number someone made up.
    """)
    return


if __name__ == "__main__":
    app.run()
