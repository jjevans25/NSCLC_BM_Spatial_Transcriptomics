# ---------------------------------------------------------------------------
# Phase 0 — annotation: sample-title parsing, design matrix, .h5ad export.
#
# Owner task: P0-T3, P0-T7 (see Markdowns/PROJECT_PLAN.md §6).
#
# P0-T3 writes config/samples.tsv, which common.smk then reads and validates on
# every subsequent Snakefile load. That ordering is deliberate, not circular:
# before the rule runs the file is absent and SAMPLES is None (empty DAG, by
# design); after it runs the file is committed and becomes the contract every
# later phase is checked against.
# ---------------------------------------------------------------------------


rule p0t3_parse_samples:
    """Parse GEO sample titles into samples.tsv + the design and batch tables."""
    input:
        soft=f"{PATHS['raw']}/GSE200563_family.soft.gz",
        filelist=f"{PATHS['raw']}/GSE200563_filelist.txt",
        matrix=f"{PATHS['raw']}/GSE200563_processed_data.txt.gz",
        compartment_map=config["samples"]["compartment_map"],
    output:
        samples=config["samples"]["tsv"],
        design=f"{PATHS['tables']}/design_matrix.tsv",
        batch=f"{PATHS['tables']}/batch_crosstab.tsv",
        summary=f"{PATHS['tables']}/design_summary.json",
    params:
        expected_n_aoi=config["samples"]["expected_n_aoi"],
    log:
        f"{PATHS['logs']}/p0t3_parse_samples.log",
    benchmark:
        f"{PATHS['benchmarks']}/p0t3_parse_samples.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/parse_samples.py"


TARGETS_ANNOTATE = [
    config["samples"]["tsv"],
    f"{PATHS['tables']}/design_matrix.tsv",
    f"{PATHS['tables']}/batch_crosstab.tsv",
    f"{PATHS['tables']}/design_summary.json",
]
