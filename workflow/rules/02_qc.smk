# ---------------------------------------------------------------------------
# Phase 0 — QC and normalisation (R), marker sanity check, power check.
#
# Owner task: P0-T4, P0-T5, P0-T6, P0-T8 (see Markdowns/PROJECT_PLAN.md §6).
#
# P0-T4 lives here rather than with acquisition because Q1's answer is what
# decides `normalisation.method` — verify vs apply — which is this file's
# business.
#
# P0-T6 lands before P0-T5 on purpose. It tests whether the compartment
# *labels* are supported by the transcriptome, which is a question about the
# annotation rather than about data quality — and §6 is explicit that a failure
# stops the project. Finding that out before investing in QC is the cheaper
# order. It reads the acquired matrix directly and should be re-run once P0-T5
# produces the QC'd, normalised object.
# ---------------------------------------------------------------------------


rule p0t4_normalisation_check:
    """Q1 — is the GEO matrix raw counts, Q3-normalised, or log-transformed?"""
    input:
        matrix=f"{PATHS['raw']}/GSE200563_processed_data.txt.gz",
    output:
        figure=report(
            f"{PATHS['figures']}/normalisation_check.png",
            caption="../report/normalisation_check.rst",
            category="Phase 0 — acquisition, annotation, QC",
            labels={"task": "P0-T4", "check": "normalisation state (Q1)"},
        ),
        per_aoi=f"{PATHS['tables']}/normalisation_check.tsv",
        summary=f"{PATHS['tables']}/normalisation_check_summary.json",
    params:
        # What the GEO metadata claims. The rule reports whether the file agrees;
        # if it does not, the file wins.
        expected_state="q3_normalised",
    log:
        f"{PATHS['logs']}/p0t4_normalisation_check.log",
    benchmark:
        f"{PATHS['benchmarks']}/p0t4_normalisation_check.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/normalisation_check.py"


rule p0t6_marker_sanity:
    """Do the compartment labels mean what they say? (PROJECT_PLAN §6, Gate 0)"""
    input:
        matrix=f"{PATHS['raw']}/GSE200563_processed_data.txt.gz",
        samples=config["samples"]["tsv"],
    output:
        figure=report(
            f"{PATHS['figures']}/marker_sanity.png",
            caption="../report/marker_sanity.rst",
            category="Phase 0 — acquisition, annotation, QC",
            labels={"task": "P0-T6", "check": "marker sanity"},
        ),
        stats=f"{PATHS['tables']}/marker_sanity.tsv",
        verdict=f"{PATHS['tables']}/marker_sanity_verdict.tsv",
        summary=f"{PATHS['tables']}/marker_sanity_summary.json",
    params:
        seed=SEED,
    log:
        f"{PATHS['logs']}/p0t6_marker_sanity.log",
    benchmark:
        f"{PATHS['benchmarks']}/p0t6_marker_sanity.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/marker_sanity.py"


rule p0t5_qc_metrics:
    """Per-AOI QC metrics and flags; verify (not apply) normalisation."""
    input:
        raw_tar=f"{PATHS['raw']}/GSE200563_RAW.tar",
        matrix=f"{PATHS['raw']}/GSE200563_processed_data.txt.gz",
        samples=config["samples"]["tsv"],
    output:
        figure=report(
            f"{PATHS['figures']}/qc_metrics.png",
            caption="../report/qc_metrics.rst",
            category="Phase 0 — acquisition, annotation, QC",
            labels={"task": "P0-T5", "check": "per-AOI QC"},
        ),
        metrics=f"{PATHS['tables']}/qc_metrics.tsv",
        excluded=f"{PATHS['tables']}/qc_excluded.tsv",
        summary=f"{PATHS['tables']}/qc_summary.json",
    params:
        # Thresholds come from config and may be null — see the script. The rule
        # reads config, never the notebook (§6).
        qc=config["qc"],
        norm_method=config["normalisation"]["method"],
    log:
        f"{PATHS['logs']}/p0t5_qc_metrics.log",
    benchmark:
        f"{PATHS['benchmarks']}/p0t5_qc_metrics.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/qc_metrics.py"


rule p0t8_power_check:
    """What effect is detectable at 80% power with TIME-B n=8? (§6)"""
    input:
        samples=config["samples"]["tsv"],
    output:
        figure=report(
            f"{PATHS['figures']}/power_curves.png",
            caption="../report/power_check.rst",
            category="Phase 0 — acquisition, annotation, QC",
            labels={"task": "P0-T8", "check": "power for TIME-L vs TIME-B"},
        ),
        power=f"{PATHS['tables']}/power.tsv",
        summary=f"{PATHS['tables']}/power_summary.json",
    params:
        seed=SEED,
        n_sim=1000,
    log:
        f"{PATHS['logs']}/p0t8_power_check.log",
    benchmark:
        f"{PATHS['benchmarks']}/p0t8_power_check.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 4
    script:
        "../scripts/power_check.py"


TARGETS_QC = [
    f"{PATHS['figures']}/power_curves.png",
    f"{PATHS['tables']}/power.tsv",
    f"{PATHS['tables']}/power_summary.json",
    f"{PATHS['figures']}/qc_metrics.png",
    f"{PATHS['tables']}/qc_metrics.tsv",
    f"{PATHS['tables']}/qc_excluded.tsv",
    f"{PATHS['tables']}/qc_summary.json",
    f"{PATHS['figures']}/normalisation_check.png",
    f"{PATHS['tables']}/normalisation_check.tsv",
    f"{PATHS['tables']}/normalisation_check_summary.json",
    f"{PATHS['figures']}/marker_sanity.png",
    f"{PATHS['tables']}/marker_sanity.tsv",
    f"{PATHS['tables']}/marker_sanity_verdict.tsv",
    f"{PATHS['tables']}/marker_sanity_summary.json",
]
