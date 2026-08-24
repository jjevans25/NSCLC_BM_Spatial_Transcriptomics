# ---------------------------------------------------------------------------
# Phase 0 — QC and normalisation (R), marker sanity check, power check.
#
# Owner task: P0-T5, P0-T6, P0-T8 (see Markdowns/PROJECT_PLAN.md §6).
#
# P0-T6 lands before P0-T5 on purpose. It tests whether the compartment
# *labels* are supported by the transcriptome, which is a question about the
# annotation rather than about data quality — and §6 is explicit that a failure
# stops the project. Finding that out before investing in QC is the cheaper
# order. It reads the acquired matrix directly and should be re-run once P0-T5
# produces the QC'd, normalised object.
# ---------------------------------------------------------------------------


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


TARGETS_QC = [
    f"{PATHS['figures']}/marker_sanity.png",
    f"{PATHS['tables']}/marker_sanity.tsv",
    f"{PATHS['tables']}/marker_sanity_verdict.tsv",
    f"{PATHS['tables']}/marker_sanity_summary.json",
]
