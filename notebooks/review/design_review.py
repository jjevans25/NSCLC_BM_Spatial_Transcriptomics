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
    # P0-T3 — Design review

    **Tier:** review (`notebooks/review/`) — not in the Snakemake DAG. A human
    runs this at a gate. It **reads** `results/` and `config/samples.tsv` and
    writes nothing.

    **Every number below was produced by `p0t3_parse_samples`.** This notebook
    re-renders and cross-filters them; it does not compute them. If something
    here needs to become a reported result, it gets promoted to the script
    first (`notebooks/README.md`, hard constraint 3).

    **What to decide here:** whether the parsed design matches PROJECT_PLAN
    §2.1, and **which contrasts actually survive** once you filter to the
    subset an analysis would use. The binding constraint is `TIME-B` n = 8.
    """)
    return


@app.cell
def _(mo):
    # Parameters (PROJECT_PLAN §A.5). Defaults make the notebook openable with a
    # bare `marimo edit`; anything automated passes explicit paths.
    args = mo.cli_args()
    samples_path = args.get("samples", "config/samples.tsv")
    design_path = args.get("design", "results/tables/design_matrix.tsv")
    batch_path = args.get("batch", "results/tables/batch_crosstab.tsv")
    summary_path = args.get("summary", "results/tables/design_summary.json")
    return batch_path, design_path, samples_path, summary_path


@app.cell
def _(batch_path, design_path, samples_path, summary_path):
    import json
    from pathlib import Path

    import pandas as pd

    samples = pd.read_csv(samples_path, sep="\t", dtype=str)
    samples["replicate_flag"] = samples["replicate_flag"] == "True"
    samples["is_control"] = samples["is_control"] == "True"

    design = pd.read_csv(design_path, sep="\t", index_col=0)
    batch = pd.read_csv(batch_path, sep="\t", index_col=0)
    summary = json.loads(Path(summary_path).read_text())
    return batch, design, pd, samples, summary


@app.cell
def _(mo, summary):
    mo.md(f"""
    ## 1. Does the parse match §2.1?

    | | |
    |---|---|
    | AOIs | **{summary["n_aoi"]}** (expected 120) |
    | Subjects | {summary["n_subjects"]} — {summary["n_patients"]} patients + {summary["n_controls"]} controls |
    | Replicate AOIs | {summary["replicate_aois"]} |

    The per-compartment counts are asserted against `compartment_map.yaml` in
    the rule, so if this notebook opens at all, they matched — the pipeline
    fails loudly rather than reporting a mismatch here.

    ### Design matrix — subjects x compartments
    """)
    return


@app.cell
def _(design):
    design
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Cross-filter: which contrasts survive?

    Pick a subset the way an analysis would, and watch the per-compartment
    counts move. **Watch `TIME-B`.** It starts at 8, and any filter that costs
    it 1-2 AOIs changes what Phase 2 can support.
    """)
    return


@app.cell
def _(mo, samples):
    site_filter = mo.ui.multiselect(
        options=sorted(samples["site"].unique()),
        value=sorted(samples["site"].unique()),
        label="Site",
    )
    compartment_filter = mo.ui.multiselect(
        options=sorted(samples["compartment"].unique()),
        value=sorted(samples["compartment"].unique()),
        label="Compartment",
    )
    run_filter = mo.ui.multiselect(
        options=sorted(samples["dsp_run"].unique()),
        value=sorted(samples["dsp_run"].unique()),
        label="DSP run (batch)",
    )
    drop_replicates = mo.ui.checkbox(
        value=False, label="Keep only the first AOI per (subject, compartment)"
    )
    exclude_controls = mo.ui.checkbox(value=False, label="Exclude BC controls")

    mo.vstack([
        mo.hstack([site_filter, compartment_filter], justify="start"),
        run_filter,
        mo.hstack([drop_replicates, exclude_controls], justify="start"),
    ])
    return (
        compartment_filter,
        drop_replicates,
        exclude_controls,
        run_filter,
        site_filter,
    )


@app.cell
def _(
    compartment_filter,
    drop_replicates,
    exclude_controls,
    run_filter,
    samples,
    site_filter,
):
    subset = samples[
        samples["site"].isin(site_filter.value)
        & samples["compartment"].isin(compartment_filter.value)
        & samples["dsp_run"].isin(run_filter.value)
    ]
    if exclude_controls.value:
        subset = subset[~subset["is_control"]]
    if drop_replicates.value:
        subset = subset.drop_duplicates(subset=["patient_id", "aoi_code"], keep="first")
    return (subset,)


@app.cell
def _(mo, samples, subset):
    surviving = (
        subset.groupby("aoi_code")
        .size()
        .reindex(["L", "LB", "mLN", "TBME", "TIME-L", "TIME-B", "BC"], fill_value=0)
    )
    full = (
        samples.groupby("aoi_code")
        .size()
        .reindex(["L", "LB", "mLN", "TBME", "TIME-L", "TIME-B", "BC"], fill_value=0)
    )
    lines = "\n".join(
        f"| `{code}` | {full[code]} | **{surviving[code]}** | {surviving[code] - full[code]:+d} |"
        for code in full.index
    )
    time_b = int(surviving["TIME-B"])
    verdict = (
        f"**`TIME-B` n = {time_b}.**"
        + (
            " Unchanged — the binding constraint still holds at its stated value."
            if time_b == 8
            else f" Down {8 - time_b} from 8. Anything resting on the brain immune"
            " compartment must be re-checked against this number, and every claim"
            " about brain immune contexture states it inline (hard constraint 8)."
        )
        if time_b
        else "**`TIME-B` is empty in this subset — no brain immune contrast exists here.**"
    )
    mo.md(f"""
    **{len(subset)} of {len(samples)} AOIs, {subset["patient_id"].nunique()} subjects.**

    | Compartment | Full | Subset | Δ |
    |---|---|---|---|
    {lines}

    {verdict}
    """)
    return


@app.cell
def _(subset):
    subset[
        ["aoi_label", "patient_id", "aoi_code", "site", "compartment",
         "replicate_flag", "dsp_run"]
    ]
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 3. Paired immune AOIs

    `CLAUDE.md` states that only five patients carry immune AOIs at **both**
    sites, and that this is a consistency check, never a headline result. This
    is where that claim is visible.
    """)
    return


@app.cell
def _(mo, samples):
    immune = samples[samples["aoi_code"].isin(["TIME-L", "TIME-B"])]
    paired = sorted(
        {
            patient
            for patient, group in immune.groupby("patient_id")
            if set(group["aoi_code"]) == {"TIME-L", "TIME-B"}
        },
        key=lambda p: int(p.lstrip("P")),
    )
    mo.md(f"""
    Patients with both `TIME-L` and `TIME-B`: **{", ".join(paired)}** (n = {len(paired)}).

    Expected per `CLAUDE.md`: P5, P12, P15, P19, P35 (n = 5).
    **{"Match." if paired == ["P5", "P12", "P15", "P19", "P35"] else "DOES NOT MATCH — investigate before Phase 2."}**
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 4. Batch — and where it is inseparable from biology

    GEO exposes no batch field. The DSP run id survives only in the DCC
    filenames (`docs/data-provenance.md`, Q3), which is why P0-T2 acquires
    `filelist.txt`.
    """)
    return


@app.cell
def _(batch):
    batch
    return


@app.cell
def _(mo, summary):
    confounded = summary["compartments_confounded_with_dsp_run"]
    mo.md(f"""
    **Compartments sitting entirely inside one DSP run: {", ".join(f"`{c}`" for c in confounded)}.**

    For these, batch and biology cannot be separated — a difference between
    `BC` and anything else is a difference between run B and run A at the same
    time. That is a limitation to state (P0-T8), not a defect to correct.

    The core comparison is safe: `TIME-L` and `TIME-B` both span the two runs,
    so a batch term is estimable there. Note it rests on only 2 `TIME-B` AOIs
    in the smaller run, so the batch estimate for the brain immune compartment
    is thin.
    """)
    return


@app.cell
def _(mo, pd, samples):
    mo.md("### Subjects contributing more than one AOI of the same compartment")
    reps = samples[samples["replicate_flag"]][
        ["patient_id", "aoi_code", "aoi_label", "dsp_run"]
    ].sort_values(["patient_id", "aoi_code"])
    display = (
        reps
        if len(reps)
        else pd.DataFrame({"note": ["none — every (subject, compartment) is unique"]})
    )
    display
    return


@app.cell
def _(mo):
    mo.md(r"""
    These are why the random intercept for patient is **mandatory, not
    stylistic** (`CLAUDE.md`). Treating these AOIs as independent inflates the
    effective n and narrows every interval that depends on it.

    ---

    ### Gate check

    - [ ] 120 AOIs, per-compartment counts match §2.1
    - [ ] Paired-immune patients are the expected five
    - [ ] `TIME-B` n = 8 before any QC filtering
    - [ ] Batch confounding understood and carried to P0-T5 / P0-T8

    Remaining for Gate 0: **P0-T6** (marker sanity check). Q2 is answered
    (`docs/data-provenance.md`).
    """)
    return


if __name__ == "__main__":
    app.run()
