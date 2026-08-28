# ---------------------------------------------------------------------------
# Phase 2 — immune contexture, brain vs. lung (A1/A3).
#
# Owner task: Phase 2 (see Markdowns/PROJECT_PLAN.md §6).
#
# Everything in this file is estimated WITHIN the TIME compartment. Gate 1
# measured compartment at 41% of variance, so pooling TIME with L/LB would
# swamp the site effect this phase exists to measure (ADR 0009).
#
# The binding constraint is TIME-B n = 8, with a power floor of ~1.1-1.3 SD at
# 80% (P0-T8). A null here is uninformative, not negative — and per Gate 2, a
# SURPRISING result is more likely a pipeline bug than biology.
# ---------------------------------------------------------------------------

# Explicit rather than relying on 00_acquire.smk's import leaking through the
# shared include namespace: this file's wildcard constraint needs it, and an
# include-order change should not be able to break that silently.
import re

CONTEXTURE = config["contexture"]

# --- P2-T5 reference matrices ----------------------------------------------
# Same acquisition contract as P0-T2 (ADR 0005): pinned digests, protected()
# outputs, verification in a SEPARATE rule so that adding a pin re-runs only the
# cheap check instead of invalidating a write-protected download.
#
# base_url pins a commit, not a branch -- the config schema enforces the
# trailing 40-hex SHA, because CellProfileLibrary is a live GitHub repo and a
# branch ref would let the reference matrix change under a rerun.
REFERENCE = config["reference"]
REF_ARTIFACTS = REFERENCE["artifacts"]
REF_KEY_BY_DEST = {spec["dest"]: key for key, spec in REF_ARTIFACTS.items()}

if len(REF_KEY_BY_DEST) != len(REF_ARTIFACTS):
    raise WorkflowError(
        "config reference.artifacts: two artifacts share a `dest` filename."
    )

REF_DIR = PATHS["reference"]
REF_META = f"{PATHS['interim']}/reference"

# Sentinel marking that an R conda env has had its relocated-path assignments
# quoted. This project's working directory contains a space ("Biomedical Data
# Science") and conda-forge's R writes those assignments unquoted, so R does not
# start at all in a freshly created env. workflow/scripts/repair_r_env.sh has
# the full diagnosis and says why this is a rule rather than a post-deploy hook
# or a --conda-prefix change; ADR 0013 records the decision.
R_ENV_OK = f"{PATHS['interim']}/env_repair/{{r_env}}.ok"


rule p2t0_repair_r_env:
    """Quote conda-forge R's relocated paths so R starts on a path with spaces.

    Runs inside the target env, so $CONDA_PREFIX resolves to the env being
    repaired. Idempotent, and it verifies R actually starts before touching the
    sentinel — a still-broken wrapper fails here rather than inside an analysis
    rule, where the error message names neither R nor the path.
    """
    output:
        ok=touch(R_ENV_OK),
    wildcard_constraints:
        r_env="r-stats|r-geomx",
    log:
        f"{PATHS['logs']}/p2t0_repair_{{r_env}}.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t0_repair_{{r_env}}.tsv"
    conda:
        "../envs/{r_env}.yaml"
    threads: 1
    shell:
        'bash workflow/scripts/repair_r_env.sh > "{log}" 2>&1'


rule p2t1_resolve_signatures:
    """Resolve config/signatures.yaml against the measured genes (P2-T1).

    config/signatures.yaml is hand-written; hard constraint 5 says nothing is a
    result until a rule produces it. This is that rule. A symbol that does not
    resolve stops the workflow rather than shrinking its set silently.
    """
    input:
        h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
        signatures=CONTEXTURE["signatures_yaml"],
    output:
        membership=f"{PATHS['tables']}/signature_membership.tsv",
        summary=f"{PATHS['tables']}/signature_membership_summary.json",
    params:
        signatures=SIGNATURES,
    log:
        f"{PATHS['logs']}/p2t1_resolve_signatures.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t1_resolve_signatures.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/resolve_signatures.py"


rule p2t2_score_signatures:
    """Score the six signatures on the 23 TIME AOIs, two ways, with coverage (P2-T2).

    Two methods because agreement between them is the robustness evidence at
    TIME-B n = 8, not because one checks the other.

    Coverage is reported per set PER SITE and never pooled: detection varies
    systematically by compartment (P0-T5) and brain background is higher
    (ADR 0008), so a near-background gene reads as depleted in brain
    artefactually. A set below the coverage floor at a site is "not assessable"
    there.
    """
    input:
        h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
        membership=f"{PATHS['tables']}/signature_membership.tsv",
    output:
        scores=f"{PATHS['tables']}/signature_scores.tsv",
        coverage=f"{PATHS['tables']}/signature_coverage.tsv",
        summary=f"{PATHS['tables']}/signature_scoring_summary.json",
    params:
        signatures=SIGNATURES,
        aoi_codes=CONTEXTURE["aoi_codes"],
        scoring=CONTEXTURE["scoring"],
        background_multiple=config["qc"]["detection_background_multiple"],
        seed=SEED,
    log:
        f"{PATHS['logs']}/p2t2_score_signatures.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t2_score_signatures.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 4
    script:
        "../scripts/score_signatures.py"


rule p2t3_signature_models:
    """The lung-vs-brain contrast per signature, lme4/lmerTest (P2-T3).

    Primary `score ~ site + (1|patient_id)`; sensitivity adds `dsp_run`.
    Both are reported, with singular fits surfaced and the estimate shift
    between them quantified, because TIME-B splits 6/2 across the two DSP runs
    and the adjusted fit is barely identified. ADR 0012.

    site is releveled to lung so the coefficient is BRAIN MINUS LUNG. Gate 2
    expects a negative antigen-presentation estimate; getting the reference
    level backwards would invert the gate's verdict, so the script asserts it.
    """
    input:
        scores=f"{PATHS['tables']}/signature_scores.tsv",
        env_ok=f"{PATHS['interim']}/env_repair/r-stats.ok",
    output:
        table=f"{PATHS['tables']}/signature_models.tsv",
        summary=f"{PATHS['tables']}/signature_models_summary.json",
    params:
        primary=CONTEXTURE["model"]["primary"],
        sensitivity=CONTEXTURE["model"]["sensitivity"],
        fdr_alpha=CONTEXTURE["model"]["fdr_alpha"],
        df_method=CONTEXTURE["model"]["df_method"],
        crosscheck_signature=CONTEXTURE["model"]["crosscheck_signature"],
    log:
        f"{PATHS['logs']}/p2t3_signature_models.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t3_signature_models.tsv"
    conda:
        "../envs/r-stats.yaml"
    threads: 2
    script:
        "../scripts/fit_signature_models.R"


rule p2t3_model_crosscheck:
    """Refit the Gate 2 signature in statsmodels and compare to lme4 (P2-T3).

    An implementation check, not a statistical one: estimates and SEs are
    compared, df and p-values are not, because lmerTest uses Satterthwaite and
    statsmodels a Wald z with no finite-sample correction.
    """
    input:
        scores=f"{PATHS['tables']}/signature_scores.tsv",
        models=f"{PATHS['tables']}/signature_models.tsv",
    output:
        crosscheck=f"{PATHS['tables']}/model_crosscheck.json",
    params:
        crosscheck_signature=CONTEXTURE["model"]["crosscheck_signature"],
        # Loose enough to absorb optimiser differences between two independent
        # REML implementations, tight enough that a real disagreement in the
        # third decimal of a signature score fails.
        tolerance=1e-4,
    log:
        f"{PATHS['logs']}/p2t3_model_crosscheck.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t3_model_crosscheck.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/crosscheck_mixedlm.py"


rule p2t4_paired_check:
    """Within-patient direction check on the 5 doubly-sampled patients (P2-T4).

    DIRECTION OF EFFECT ONLY. No p-values from n = 5 (PROJECT_PLAN §2.3): the
    script emits no test statistic at all, and asserts that none reached the
    table. The paired and unpaired results share AOIs and are not independent
    evidence.
    """
    input:
        scores=f"{PATHS['tables']}/signature_scores.tsv",
        models=f"{PATHS['tables']}/signature_models.tsv",
    output:
        deltas=f"{PATHS['tables']}/paired_deltas.tsv",
        concordance=f"{PATHS['tables']}/paired_concordance.tsv",
        summary=f"{PATHS['tables']}/paired_summary.json",
    log:
        f"{PATHS['logs']}/p2t4_paired_check.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t4_paired_check.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/paired_check.py"


rule p2t7_contexture_heatmap:
    """The Phase 2 deliverable figure: signatures x AOIs, grouped by site (P2-T7).

    Plots the z-score mean, because it is the method whose magnitude is
    comparable between signatures of different size; the ssGSEA result is in
    the score table and in every model row.

    Sets that failed their coverage floor are hatched ON THE FIGURE. ADR 0008's
    rule is "not assessable at that site", never "lower there", and a heatmap
    is exactly the artefact someone reads "lower in brain" off a colour from.
    """
    input:
        scores=f"{PATHS['tables']}/signature_scores.tsv",
        coverage=f"{PATHS['tables']}/signature_coverage.tsv",
    output:
        figure=report(
            f"{PATHS['figures']}/contexture_heatmap.png",
            caption="../report/contexture_heatmap.rst",
            category="Phase 2 — immune contexture",
            labels={"task": "P2-T7", "figure": "contexture heatmap"},
        ),
    params:
        aoi_codes=CONTEXTURE["aoi_codes"],
    log:
        f"{PATHS['logs']}/p2t7_contexture_heatmap.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t7_contexture_heatmap.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/contexture_heatmap.py"


rule p2t5_fetch_reference:
    """Download one SpatialDecon profile matrix into resources/reference (P2-T5).

    Reuses workflow/scripts/acquire_geo.py unchanged -- nothing in it is
    GEO-specific; it takes a url, a key and a stability tier and writes the
    artifact plus a provenance sidecar.
    """
    wildcard_constraints:
        # Named `artifact`, not `reference`, because acquire_geo.py reads
        # snakemake.wildcards.artifact by name. Renaming the wildcard here is
        # the cheap side of that coupling: editing the Phase 0 script instead
        # would fire the code rerun-trigger on p0t2_fetch_geo and abort the DAG
        # on its protected() outputs.
        artifact="|".join(re.escape(dest) for dest in sorted(REF_KEY_BY_DEST)),
    output:
        artifact=protected(f"{REF_DIR}/{{artifact}}"),
        meta=f"{REF_META}/{{artifact}}.json",
    params:
        url=lambda w: (
            f"{REFERENCE['base_url']}/"
            f"{REF_ARTIFACTS[REF_KEY_BY_DEST[w.artifact]]['remote']}"
        ),
        key=lambda w: REF_KEY_BY_DEST[w.artifact],
        stability=lambda w: REF_ARTIFACTS[REF_KEY_BY_DEST[w.artifact]]["stability"],
    log:
        f"{PATHS['logs']}/p2t5_fetch_reference_{{artifact}}.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t5_fetch_reference_{{artifact}}.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/acquire_geo.py"


rule p2t5_record_reference_provenance:
    """Verify the reference matrices against their pinned digests (P2-T5).

    Reuses workflow/scripts/record_provenance.py unchanged -- it resolves paths
    from its inputs rather than a hardcoded root, so it writes any pair of
    checksum/provenance outputs.
    """
    input:
        artifacts=[f"{REF_DIR}/{dest}" for dest in sorted(REF_KEY_BY_DEST)],
        meta=[f"{REF_META}/{dest}.json" for dest in sorted(REF_KEY_BY_DEST)],
    output:
        checksums=f"{PATHS['resources']}/reference_checksums.sha256",
        provenance=f"{PATHS['resources']}/reference_provenance.tsv",
    params:
        artifacts=REF_ARTIFACTS,
    log:
        f"{PATHS['logs']}/p2t5_record_reference_provenance.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t5_record_reference_provenance.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/record_provenance.py"


rule p2t5_deconvolution:
    """Cell composition of the TIME AOIs, two reference arms (P2-T5).

    safeTME on BOTH sites -- one reference, so that arm's composition is
    comparable brain vs lung. Plus PROJECT_PLAN §6's two-matrix approach as the
    sensitivity, which is reference-confounded across sites and says so.

    Reads the plain TSVs p0t7_export_tsv already writes, which is the ADR 0001
    boundary used as designed -- no zellkonverter, no second export.
    """
    input:
        expr_q3=f"{PATHS['interim']}/expr_q3.tsv",
        obs=f"{PATHS['interim']}/obs.tsv",
        var=f"{PATHS['interim']}/var.tsv",
        references=[f"{REF_DIR}/{dest}" for dest in sorted(REF_KEY_BY_DEST)],
        provenance=f"{PATHS['resources']}/reference_checksums.sha256",
        env_ok=f"{PATHS['interim']}/env_repair/r-geomx.ok",
    output:
        table=f"{PATHS['tables']}/decon_composition.tsv",
        summary=f"{PATHS['tables']}/decon_summary.json",
    params:
        aoi_codes=CONTEXTURE["aoi_codes"],
        primary_reference=CONTEXTURE["deconvolution"]["primary_reference"],
        ref_by_code={
            code: f"{REF_DIR}/{ref}.RData"
            for code, ref in CONTEXTURE["deconvolution"][
                "sensitivity_references"
            ].items()
        },
    log:
        f"{PATHS['logs']}/p2t5_deconvolution.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t5_deconvolution.tsv"
    conda:
        "../envs/r-geomx.yaml"
    threads: 2
    script:
        "../scripts/deconvolve.R"


rule p2t5_decon_figure:
    """Stacked composition, one panel per reference arm (P2-T5).

    Three panels rather than two: the two tissue references resolve different
    cell types (9 lymphoid vs 0), so drawing them together with a shared legend
    would invite exactly the cross-site comparison that arm cannot support.
    """
    input:
        table=f"{PATHS['tables']}/decon_composition.tsv",
        summary=f"{PATHS['tables']}/decon_summary.json",
    output:
        figure=report(
            f"{PATHS['figures']}/decon_composition.png",
            caption="../report/decon_composition.rst",
            category="Phase 2 — immune contexture",
            labels={"task": "P2-T5", "figure": "cell composition"},
        ),
    params:
        aoi_codes=CONTEXTURE["aoi_codes"],
    log:
        f"{PATHS['logs']}/p2t5_decon_figure.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t5_decon_figure.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/decon_figure.py"


rule p2t6_convergence_check:
    """Do deconvolution and signature scoring agree on direction? (P2-T6)

    Answering that honestly first requires saying for which signatures the
    question is askable at all. Three things disqualify one -- no lineage
    counterpart, below its detection-coverage floor, or a mapping shared with
    another signature -- and the script tests all three rather than reporting a
    disagreement as though it were about biology.

    safeTME arm only: the two-matrix arm uses a different reference per site and
    Brain_Darmanis resolves no lymphoid type, so a cross-site direction from it
    would be the reference rather than the biology.
    """
    input:
        models=f"{PATHS['tables']}/signature_models.tsv",
        coverage=f"{PATHS['tables']}/signature_coverage.tsv",
        composition=f"{PATHS['tables']}/decon_composition.tsv",
    output:
        table=f"{PATHS['tables']}/convergence_check.tsv",
        summary=f"{PATHS['tables']}/convergence_summary.json",
    params:
        signatures=SIGNATURES,
        aoi_codes=CONTEXTURE["aoi_codes"],
    log:
        f"{PATHS['logs']}/p2t6_convergence_check.log",
    benchmark:
        f"{PATHS['benchmarks']}/p2t6_convergence_check.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/convergence_check.py"


TARGETS_CONTEXTURE = [
    f"{PATHS['tables']}/signature_membership.tsv",
    f"{PATHS['tables']}/signature_membership_summary.json",
    f"{PATHS['tables']}/signature_scores.tsv",
    f"{PATHS['tables']}/signature_coverage.tsv",
    f"{PATHS['tables']}/signature_scoring_summary.json",
    f"{PATHS['tables']}/signature_models.tsv",
    f"{PATHS['tables']}/signature_models_summary.json",
    f"{PATHS['tables']}/model_crosscheck.json",
    f"{PATHS['tables']}/paired_deltas.tsv",
    f"{PATHS['tables']}/paired_concordance.tsv",
    f"{PATHS['tables']}/paired_summary.json",
    f"{PATHS['figures']}/contexture_heatmap.png",
    f"{PATHS['resources']}/reference_checksums.sha256",
    f"{PATHS['resources']}/reference_provenance.tsv",
    f"{PATHS['tables']}/decon_composition.tsv",
    f"{PATHS['tables']}/decon_summary.json",
    f"{PATHS['figures']}/decon_composition.png",
    f"{PATHS['tables']}/convergence_check.tsv",
    f"{PATHS['tables']}/convergence_summary.json",
]
