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
    # 00 — First look at GSE200563

    **Tier:** explore (`notebooks/explore/`) — not in the Snakemake DAG, may write
    to `results/interim/` only, and **never** the source of a reported number.

    **Purpose:** answer **Q1** from PROJECT_PLAN §3 — is
    `GSE200563_processed_data.txt.gz` raw counts, Q3-normalised, or already
    log-transformed? The checks that settle it:

    | Check | Raw counts | Q3-normalised | Log-transformed |
    |---|---|---|---|
    | All values integer | yes | no | no |
    | Value range | 0 – 10⁴⁺ | 0 – 10⁴⁺, rescaled | roughly 0 – 20 |
    | Column-wise Q3 spread across AOIs | wide | ~constant | wide |
    | Minimum non-zero value | 1 | fractional | fractional |

    Whatever this notebook concludes gets written up in `docs/data-provenance.md`
    by **P0-T4** and drives `normalisation.method` in `config/config.yaml`. The
    conclusion lives there; this notebook is only how it was reached.

    **Status: stub.** The matrix arrives with P0-T2 — nothing to inspect yet.
    """)
    return


@app.cell
def _(mo):
    # Parameters (PROJECT_PLAN §A.5). Defaults make the notebook openable with a
    # bare `marimo edit`; anything automated passes explicit paths.
    args = mo.cli_args()
    matrix_path = args.get(
        "matrix", "resources/raw/GSE200563_processed_data.txt.gz"
    )
    config_path = args.get("config", "config/config.yaml")
    return config_path, matrix_path


@app.cell
def _(config_path, matrix_path, mo):
    from pathlib import Path

    have_matrix = Path(matrix_path).exists()

    mo.md(
        f"""
    - matrix: `{matrix_path}` — **{"present" if have_matrix else "not downloaded yet (run P0-T2)"}**
    - config: `{config_path}`
    """
    )
    return (have_matrix,)


@app.cell
def _(have_matrix, mo):
    mo.stop(
        not have_matrix,
        mo.md(
            "**Waiting on P0-T2.** Download the GEO matrix with "
            "`snakemake acquire`, then re-run this notebook to answer Q1."
        ),
    )

    mo.md("TODO (P0-T4): integer check, value range, column-wise Q3 spread.")
    return


if __name__ == "__main__":
    app.run()
