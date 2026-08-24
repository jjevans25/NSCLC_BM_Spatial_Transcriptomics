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


rule p0t7_export_tsv:
    """R/Python boundary as plain files (ADR 0001): write the handoff TSVs."""
    input:
        matrix=f"{PATHS['raw']}/GSE200563_processed_data.txt.gz",
        qc_metrics=f"{PATHS['tables']}/qc_metrics.tsv",
        samples=config["samples"]["tsv"],
    output:
        expr=f"{PATHS['interim']}/expr_normalised.tsv",
        expr_q3=f"{PATHS['interim']}/expr_q3.tsv",
        obs=f"{PATHS['interim']}/obs.tsv",
        var=f"{PATHS['interim']}/var.tsv",
        uns=f"{PATHS['interim']}/uns.json",
    params:
        norm_method=config["normalisation"]["method"],
        log_base=config["normalisation"]["log_base"],
        pseudocount=config["normalisation"]["pseudocount"],
        background_multiple=config["qc"]["detection_background_multiple"],
    log:
        f"{PATHS['logs']}/p0t7_export_tsv.log",
    benchmark:
        f"{PATHS['benchmarks']}/p0t7_export_tsv.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/export_tsv.py"


rule p0t7_assemble_h5ad:
    """Assemble the .h5ad and assert the float round-trip is exact (§6)."""
    input:
        expr=f"{PATHS['interim']}/expr_normalised.tsv",
        expr_q3=f"{PATHS['interim']}/expr_q3.tsv",
        obs=f"{PATHS['interim']}/obs.tsv",
        var=f"{PATHS['interim']}/var.tsv",
        uns=f"{PATHS['interim']}/uns.json",
        samples=config["samples"]["tsv"],
        config="config/config.yaml",
    output:
        h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
        summary=f"{PATHS['tables']}/h5ad_summary.json",
    params:
        project_code=config["project"]["code"],
        geo_accession=config["project"]["geo_accession"],
        source_pmid=config["project"]["source_pmid"],
    log:
        f"{PATHS['logs']}/p0t7_assemble_h5ad.log",
    benchmark:
        f"{PATHS['benchmarks']}/p0t7_assemble_h5ad.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/assemble_h5ad.py"


TARGETS_ANNOTATE = [
    f"{PATHS['interim']}/aoi_normalised.h5ad",
    f"{PATHS['tables']}/h5ad_summary.json",
    config["samples"]["tsv"],
    f"{PATHS['tables']}/design_matrix.tsv",
    f"{PATHS['tables']}/batch_crosstab.tsv",
    f"{PATHS['tables']}/design_summary.json",
]
