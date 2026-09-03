# ---------------------------------------------------------------------------
# Phase 3 — checkpoint landscape by compartment (A4).
#
# Owner task: Phase 3 (see Markdowns/PROJECT_PLAN.md §6).
#
# A4 IS EXPLORATORY, NOT A HEADLINE AIM (ADR 0008). Everything here is still
# computed — the demotion is about what may be CLAIMED, not about what is built.
# Three rules bind every output in this file:
#
#   1. Every checkpoint gene reported states its per-compartment detection.
#      An expression value without it is not reportable.
#   2. A gene below the pre-registered floor is "not assessable in
#      <compartment>", NEVER "lower in <compartment>". Brain background is
#      HIGHER (log2 4.96 in L to 5.67 in BC), so a near-background gene reads as
#      depleted in brain artefactually — which is precisely the headline A4 is
#      shaped to produce.
#   3. Every brain claim states TIME-B n = 8 inline (hard constraint 8).
#
# Phase 2 made outcome 2 the LIKELY one, not the worst case: the panel overlaps
# the `exhaustion` signature almost entirely and `exhaustion` cleared 1 of 6
# genes in brain. Gate 3 already licenses that — "if everything is
# not-assessable, that's still a legitimate finding". What would make this phase
# weak is not a null; it is a null reported without the detection table that
# licenses it.
#
# Models run ONLY on the four compartments spanning both DSP runs. TBME, mLN and
# BC sit wholly inside one run (Q3) and ADR 0009 §3 is explicit that no
# covariate recovers that; they are audited and plotted, never modelled.
# ---------------------------------------------------------------------------

# Imported explicitly rather than inherited from 04_contexture.smk's namespace.
# Includes share one global namespace, so `re` is in scope either way — but a
# rule whose wildcard constraint depends on which .smk was included first is a
# reordering away from a confusing NameError.
import re

CHECKPOINT = config["checkpoints"]


rule p3t1_resolve_checkpoints:
    """Resolve config/checkpoints.yaml against the measured genes (P3-T1).

    config/checkpoints.yaml is hand-written; hard constraint 5 says nothing is a
    result until a rule produces it. This is that rule. A symbol that does not
    resolve stops the workflow rather than shrinking the panel silently.

    It does NOT filter by detection. The panel is pre-registered (ADR 0015) and
    every gene reaches every downstream table, including the ones that will read
    "not assessable" everywhere — which genes are measurable is the Phase 3
    result, not a Phase 3 input.
    """
    input:
        h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
        panel=CHECKPOINT["panel_yaml"],
    output:
        membership=f"{PATHS['tables']}/checkpoint_membership.tsv",
        summary=f"{PATHS['tables']}/checkpoint_membership_summary.json",
    params:
        checkpoints=CHECKPOINTS,
        signature_overlap=CHECKPOINT_SIGNATURE_OVERLAP,
    log:
        f"{PATHS['logs']}/p3t1_resolve_checkpoints.log",
    benchmark:
        f"{PATHS['benchmarks']}/p3t1_resolve_checkpoints.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/resolve_checkpoints.py"


rule p3t2_checkpoint_detection:
    """Detection per checkpoint gene per compartment — the Phase 3 deliverable (P3-T2).

    Not a gate to get past. ADR 0008 demoted A4 because most of the panel sits
    at background in brain, and Gate 3 says "if everything is not-assessable,
    that's still a legitimate finding". A null here is only reportable BECAUSE
    of this table, so the table is the result.

    Reuses the project's one detection rule — q3 > qc.detection_background_multiple
    x that AOI's NegProbe-WTX (ADR 0007) — lifted from score_signatures.py rather
    than rewritten.

    Keyed on aoi_code, never on compartment: `compartment` is degenerate across
    sites (L, LB and mLN are all `tumour`), so grouping by it would pool lung and
    brain tumour AOIs and destroy the contrast this phase measures.

    Audits all seven compartments; the long expression table it emits carries
    only the four that span both DSP runs, which makes accidental pooling across
    the single-run compartments impossible rather than merely discouraged
    (ADR 0016).

    Three assertions carry this rule, two of them against INDEPENDENT
    implementations: per-gene totals must equal var['detected_in_n_aoi'] from
    P0-T7's export_tsv.py, and the TIME-L/TIME-B counts must equal ADR 0008's
    reconnaissance table exactly. A disagreement stops the phase.
    """
    input:
        h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
        membership=f"{PATHS['tables']}/checkpoint_membership.tsv",
    output:
        detection=f"{PATHS['tables']}/checkpoint_detection.tsv",
        expression=f"{PATHS['tables']}/checkpoint_expression.tsv",
        summary=f"{PATHS['tables']}/checkpoint_detection_summary.json",
    params:
        checkpoints=CHECKPOINTS,
        audit_aoi_codes=CHECKPOINT["audit_aoi_codes"],
        model_aoi_codes=CHECKPOINT["model_aoi_codes"],
        detection=CHECKPOINT["detection"],
        background_multiple=config["qc"]["detection_background_multiple"],
    log:
        f"{PATHS['logs']}/p3t2_checkpoint_detection.log",
    benchmark:
        f"{PATHS['benchmarks']}/p3t2_checkpoint_detection.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 2
    script:
        "../scripts/checkpoint_detection.py"


# --- P3-T2b external validation --------------------------------------------
# The source publication deposits its Source Data. That makes this project's
# matrix checkable against the authors' own file rather than against itself,
# which is the one thing P3-T2's 25 cross-checks cannot do: they run three
# implementations of ADR 0007's rule over one input column and one hypothesis.
#
# Fetched under ADR 0005's rules into a THIRD resources subdirectory, alongside
# raw/ (P0-T2, GEO) and reference/ (P2-T5, SpatialDecon). Same two-rule split
# for the same reason — the digests cannot be params of the fetch rule, whose
# outputs are protected().
SUPP_DIR = f"{PATHS['resources']}/supplementary"
SUPP_META = f"{PATHS['interim']}/supplementary"

if EXTERNAL_VALIDATION is not None:
    EXTERNAL_ARTIFACTS = EXTERNAL_VALIDATION["artifacts"]

    rule p3t2b_fetch_supplementary:
        """Download one supplementary file from the source publication (P3-T2b).

        Reuses workflow/scripts/acquire_geo.py unchanged, exactly as
        p2t5_fetch_reference does — nothing in it is GEO-specific. The wildcard
        must be named `artifact` because that script reads
        snakemake.wildcards.artifact by name; editing it instead would fire the
        code rerun-trigger on p0t2_fetch_geo and abort the DAG on its
        protected() outputs (ADR 0005).

        Springer static-content, not nature.com: the article URL 303s to an IdP
        authorize endpoint and returns no file.
        """
        wildcard_constraints:
            artifact="|".join(
                re.escape(dest) for dest in sorted(EXTERNAL_KEY_BY_DEST)
            ),
        output:
            artifact=protected(f"{SUPP_DIR}/{{artifact}}"),
            meta=f"{SUPP_META}/{{artifact}}.json",
        params:
            url=lambda w: (
                f"{EXTERNAL_VALIDATION['base_url']}/"
                f"{EXTERNAL_ARTIFACTS[EXTERNAL_KEY_BY_DEST[w.artifact]]['remote']}"
            ),
            key=lambda w: EXTERNAL_KEY_BY_DEST[w.artifact],
            stability=lambda w: EXTERNAL_ARTIFACTS[
                EXTERNAL_KEY_BY_DEST[w.artifact]
            ]["stability"],
        log:
            f"{PATHS['logs']}/p3t2b_fetch_supplementary_{{artifact}}.log",
        benchmark:
            f"{PATHS['benchmarks']}/p3t2b_fetch_supplementary_{{artifact}}.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/acquire_geo.py"

    rule p3t2b_record_supplementary_provenance:
        """Verify the supplementary files against their pinned digests (P3-T2b).

        Reuses workflow/scripts/record_provenance.py unchanged. Every artifact
        here is `fixed`: these are author-uploaded static-content files, and
        unlike GEO's SOFT family file nothing regenerates them server-side, so
        drift is a failure rather than an observation.
        """
        input:
            artifacts=[
                f"{SUPP_DIR}/{dest}" for dest in sorted(EXTERNAL_KEY_BY_DEST)
            ],
            meta=[
                f"{SUPP_META}/{dest}.json" for dest in sorted(EXTERNAL_KEY_BY_DEST)
            ],
        output:
            checksums=f"{PATHS['resources']}/supplementary_checksums.sha256",
            provenance=f"{PATHS['resources']}/supplementary_provenance.tsv",
        params:
            artifacts=EXTERNAL_ARTIFACTS,
        log:
            f"{PATHS['logs']}/p3t2b_record_supplementary_provenance.log",
        benchmark:
            f"{PATHS['benchmarks']}/p3t2b_record_supplementary_provenance.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/record_provenance.py"

    rule p3t2b_external_validation:
        """layers['q3'] against the publication's own Source Data (P3-T2b).

        The first check in this project that can fail for a reason other than
        its own arithmetic. Four assertions, each licensing the next:

          1. Every one of the 120 x 18,694 values equals the deposited value
             exactly, and the published NegProbe-WTX row equals obs['negprobe'].
             Anything else stops the phase.
          2. The published "TBME15" column equals TBME15a and differs from
             TBME15b — which is the whole of the paper's 119-vs-120 gap, and
             answers the open thread in docs/limitations.md by measurement.
          3. Per-AOI count quantisation, recovered from the value lattice. A
             detection floor is a claim about background; implied counts are a
             claim about resolution.
          4. All 99 published log2FC values and p-values reproduce, which
             establishes the estimator (mean-of-log2, Student's t) and is what
             licenses any statement about which genes the published list omits.

        It changes NO threshold and NO panel. The genes the paper names are
        reported next to ours and are not adopted: config/checkpoints.yaml is
        pre-registered (ADR 0015) and adding to it after seeing an audit is the
        failure P3-T2 exists to prevent.

        openpyxl is deliberately absent — see the script docstring. Adding it to
        py-analysis.yaml would fire the software-env trigger on p0t2_fetch_geo.
        """
        input:
            h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
            membership=f"{PATHS['tables']}/checkpoint_membership.tsv",
            detection=f"{PATHS['tables']}/checkpoint_detection.tsv",
            supplementary=[
                f"{SUPP_DIR}/{dest}" for dest in sorted(EXTERNAL_KEY_BY_DEST)
            ],
            checksums=f"{PATHS['resources']}/supplementary_checksums.sha256",
            # P0-T2's archive, read here for the ONLY thing the processed matrix
            # cannot supply: integer raw counts. Read-only and already
            # protected(); this rule never untars it to disk.
            raw_archive=f"{PATHS['raw']}/GSE200563_RAW.tar",
            samples=config["samples"]["tsv"],
        output:
            matrix=f"{PATHS['tables']}/external_validation_matrix.tsv",
            quantisation=f"{PATHS['tables']}/aoi_quantisation.tsv",
            counts=f"{PATHS['tables']}/checkpoint_count_quantisation.tsv",
            reproduction=f"{PATHS['tables']}/published_deg_reproduction.tsv",
            contrast=f"{PATHS['tables']}/published_fibrosis_contrast.tsv",
            summary=f"{PATHS['tables']}/external_validation_summary.json",
        params:
            external=EXTERNAL_VALIDATION,
            checkpoints=CHECKPOINTS,
        log:
            f"{PATHS['logs']}/p3t2b_external_validation.log",
        benchmark:
            f"{PATHS['benchmarks']}/p3t2b_external_validation.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 2
        script:
            "../scripts/external_validation.py"


rule p3t3_carrier_models:
    """Which compartment carries each checkpoint, within site (P3-T3).

    `expression ~ compartment + (1|patient_id)`, fitted separately within lung
    (L vs TIME-L) and within brain (LB vs TIME-B). Never pooled across sites:
    `compartment` is degenerate there (L, LB and mLN are all `tumour`), so a
    pooled fit would silently mix the two contrasts this phase exists to
    separate.

    TWO TABLES (ADR 0018, decided after the P3-T2 audit and before any fit):

      primary      genes clearing the floor in BOTH compartments — ADR 0016 §3
                   exactly as pre-registered. The ONLY Phase 3 table with a
                   q-value. Genes that fail get a not_assessable row carrying
                   their detection counts and no estimate.
      exploratory  genes clearing the floor in AT LEAST ONE compartment, same
                   formulas, NO FDR, every row labelled and carrying the
                   background-gradient direction. Exploratory-within-exploratory:
                   no P3-T7 sentence may rest on it.

    The exploratory fit re-derives the primary for genes in both, and the script
    asserts they agree to 1e-10 — the two tables cannot silently diverge.

    Reference level is the TUMOUR compartment and the script asserts it, so the
    coefficient is immune minus tumour. Getting that backwards would invert the
    sign of the phase's central number.
    """
    input:
        expression=f"{PATHS['tables']}/checkpoint_expression.tsv",
        env_ok=f"{PATHS['interim']}/env_repair/r-stats.ok",
    output:
        primary=f"{PATHS['tables']}/checkpoint_carrier_models.tsv",
        exploratory=f"{PATHS['tables']}/checkpoint_carrier_exploratory.tsv",
        summary=f"{PATHS['tables']}/checkpoint_carrier_summary.json",
    params:
        task="P3-T3",
        contrast_var="compartment",
        reference_level="tumour",
        test_level="immune",
        strata_var="site",
        primary=CHECKPOINT["model"]["carrier"],
        batch_sensitivity=CHECKPOINT["model"]["carrier_batch_sensitivity"],
        background_sensitivity=CHECKPOINT["model"][
            "carrier_background_sensitivity"
        ],
        detected_in_aoi_fraction=CHECKPOINT["detection"][
            "detected_in_aoi_fraction"
        ],
        df_method=CHECKPOINT["model"]["df_method"],
        fdr_alpha=CHECKPOINT["model"]["fdr_alpha"],
    log:
        f"{PATHS['logs']}/p3t3_carrier_models.log",
    benchmark:
        f"{PATHS['benchmarks']}/p3t3_carrier_models.tsv"
    conda:
        "../envs/r-stats.yaml"
    threads: 2
    script:
        "../scripts/fit_checkpoint_models.R"


rule p3t3_checkpoint_crosscheck:
    """Refit the pre-registered gene in statsmodels and compare to lme4 (P3-T3).

    An implementation check, and the end-to-end numeric guarantee across the
    R/Python boundary: the R script asserts only structure on read, because
    R_strtod is not correctly rounded (ADR 0013 §3). Estimates and SEs are
    compared; df and p-values are not, because lmerTest uses Satterthwaite and
    statsmodels a Wald z with no finite-sample correction.
    """
    input:
        expression=f"{PATHS['tables']}/checkpoint_expression.tsv",
        models=f"{PATHS['tables']}/checkpoint_carrier_models.tsv",
    output:
        crosscheck=f"{PATHS['tables']}/checkpoint_carrier_crosscheck.json",
    params:
        task="P3-T3",
        crosscheck_gene=CHECKPOINT["model"]["crosscheck_gene"],
        contrast_var="compartment",
        reference_level="tumour",
        test_level="immune",
        strata_var="site",
        # Same tolerance as p2t3_model_crosscheck: loose enough to absorb
        # optimiser differences between two independent REML implementations,
        # tight enough that a real disagreement fails.
        tolerance=1e-4,
    log:
        f"{PATHS['logs']}/p3t3_checkpoint_crosscheck.log",
    benchmark:
        f"{PATHS['benchmarks']}/p3t3_checkpoint_crosscheck.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/crosscheck_checkpoint_mixedlm.py"


rule p3t4_shift_models:
    """Does the checkpoint profile shift lung -> brain, within compartment (P3-T4).

    `expression ~ site + (1|patient_id)`, fitted separately within the tumour
    compartment (L vs LB) and within the immune compartment (TIME-L vs TIME-B).
    NEVER pooled across compartments: a tumour-compartment CD274 shift and an
    immune-compartment CD274 shift mean different things therapeutically, and
    ADR 0009 licenses the within-compartment contrast specifically.

    Drives fit_checkpoint_models.R UNCHANGED — the same script P3-T3 uses, with
    a different contrast variable. The arithmetic is identical and a second copy
    would be a second place for the reference level to drift.

    `site` is releveled to lung and asserted, so estimates are BRAIN MINUS LUNG.

    FDR is a DIFFERENT FAMILY from P3-T3's: BH across genes within compartment
    within model, where T3 corrects within site. T3 and T4 answer different
    questions from the same AOIs; correcting across both at once would be a
    third family nobody asked for (config.yaml, checkpoints.model).

    Every brain claim states TIME-B n = 8 inline (hard constraint 8), and the
    1.1-1.3 SD power floor means a null here is uninformative, not negative.
    """
    input:
        expression=f"{PATHS['tables']}/checkpoint_expression.tsv",
        env_ok=f"{PATHS['interim']}/env_repair/r-stats.ok",
    output:
        primary=f"{PATHS['tables']}/checkpoint_shift_models.tsv",
        exploratory=f"{PATHS['tables']}/checkpoint_shift_exploratory.tsv",
        summary=f"{PATHS['tables']}/checkpoint_shift_summary.json",
    params:
        task="P3-T4",
        contrast_var="site",
        reference_level="lung",
        test_level="brain",
        strata_var="compartment",
        primary=CHECKPOINT["model"]["shift"],
        batch_sensitivity=CHECKPOINT["model"]["shift_batch_sensitivity"],
        background_sensitivity=CHECKPOINT["model"]["shift_background_sensitivity"],
        detected_in_aoi_fraction=CHECKPOINT["detection"][
            "detected_in_aoi_fraction"
        ],
        df_method=CHECKPOINT["model"]["df_method"],
        fdr_alpha=CHECKPOINT["model"]["fdr_alpha"],
    log:
        f"{PATHS['logs']}/p3t4_shift_models.log",
    benchmark:
        f"{PATHS['benchmarks']}/p3t4_shift_models.tsv"
    conda:
        "../envs/r-stats.yaml"
    threads: 2
    script:
        "../scripts/fit_checkpoint_models.R"


rule p3t4_checkpoint_crosscheck:
    """Refit the pre-registered gene in statsmodels for the shift models (P3-T4)."""
    input:
        expression=f"{PATHS['tables']}/checkpoint_expression.tsv",
        models=f"{PATHS['tables']}/checkpoint_shift_models.tsv",
    output:
        crosscheck=f"{PATHS['tables']}/checkpoint_shift_crosscheck.json",
    params:
        task="P3-T4",
        crosscheck_gene=CHECKPOINT["model"]["crosscheck_gene"],
        contrast_var="site",
        reference_level="lung",
        test_level="brain",
        strata_var="compartment",
        tolerance=1e-4,
    log:
        f"{PATHS['logs']}/p3t4_checkpoint_crosscheck.log",
    benchmark:
        f"{PATHS['benchmarks']}/p3t4_checkpoint_crosscheck.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/crosscheck_checkpoint_mixedlm.py"


rule p3t4_checkpoint_paired_check:
    """Within-patient direction check on the paired patients (P3-T4).

    A NEW FILE rather than an edit to paired_check.py: editing that would fire
    the code rerun-trigger on p2t4_paired_check and re-run Phase 2 for no
    numerical gain.

    DIRECTION ONLY. NO TEST STATISTIC, in either compartment. The immune
    compartment has 5 paired patients and PROJECT_PLAN §2.3 forbids a p-value
    there in terms. The tumour compartment has 23, which would support one --
    and it still does not get one, because that would be a test this phase never
    pre-registered, introduced at the point it became available. The script
    asserts no p-value column can reach either table.
    """
    input:
        expression=f"{PATHS['tables']}/checkpoint_expression.tsv",
        primary=f"{PATHS['tables']}/checkpoint_shift_models.tsv",
        exploratory=f"{PATHS['tables']}/checkpoint_shift_exploratory.tsv",
    output:
        deltas=f"{PATHS['tables']}/checkpoint_paired_deltas.tsv",
        concordance=f"{PATHS['tables']}/checkpoint_paired_concordance.tsv",
        summary=f"{PATHS['tables']}/checkpoint_paired_summary.json",
    params:
        # PROJECT_PLAN §2.3. Asserted, not discovered: if this moves, the subset
        # or the design changed and the phase should stop.
        expected_immune_paired=5,
    log:
        f"{PATHS['logs']}/p3t4_checkpoint_paired_check.log",
    benchmark:
        f"{PATHS['benchmarks']}/p3t4_checkpoint_paired_check.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/checkpoint_paired_check.py"


rule p3t5_checkpoint_dotplot:
    """The Phase 3 deliverable figure: genes x compartment x site (P3-T5).

    A DETECTION figure. No model output reaches it, deliberately: P3-T3 and
    P3-T4 estimate effects on expression, and detection and expression move
    OPPOSITE ways under the same background gradient (ADR 0018). Detection is
    q3 above a multiple of background, so higher background raises the bar and a
    gene looks depleted; expression is log2 of that value, where background adds
    to the signal and it looks enriched. Putting an effect size on the same
    canvas as a detection rate is how those two get conflated, and that
    conflation already cost this phase one inverted argument.

    Colour is median_negprobe_ratio, NOT PROJECT_PLAN section 6's literal "mean
    expression" -- ADR 0016 section 5, decided before this figure existed. On a
    mean-log2 scale an undetected gene still shows its background level, and
    brain background is higher, so undetected genes would render BRIGHTER IN
    BRAIN than in lung: the ADR 0008 artefact arriving in the deliverable. Both
    mean-log2 columns stay in the table. The colour midpoint is read from
    qc.detection_background_multiple so it cannot drift from the rule it depicts.

    Genes below the pre-registered floor are marked ON THE FIGURE, the same mark
    and the same reasoning as p2t7_contexture_heatmap: the finding is "not
    assessable in that compartment", never "lower in that compartment". 40 of 63
    cells are marked, which dominates the figure and is correct -- the panel
    clears the floor in TIME-L alone, and Gate 3 licenses that as a finding.

    All seven compartments, because the audit is the phase's deliverable under
    Gate 3. The three that sit wholly within one DSP run carry a dagger and a
    footnote rather than a fade -- they are not lesser data (ADR 0016 section 2).
    """
    input:
        detection=f"{PATHS['tables']}/checkpoint_detection.tsv",
        summary=f"{PATHS['tables']}/checkpoint_detection_summary.json",
    output:
        figure=report(
            f"{PATHS['figures']}/checkpoint_dotplot.png",
            caption="../report/checkpoint_dotplot.rst",
            category="Phase 3 — checkpoint landscape",
            labels={"task": "P3-T5", "figure": "checkpoint detection dot plot"},
        ),
    params:
        audit_aoi_codes=CHECKPOINT["audit_aoi_codes"],
        model_aoi_codes=CHECKPOINT["model_aoi_codes"],
        detected_in_aoi_fraction=CHECKPOINT["detection"][
            "detected_in_aoi_fraction"
        ],
        background_multiple=config["qc"]["detection_background_multiple"],
        expected_n_aoi=config["samples"]["expected_n_aoi"],
    log:
        f"{PATHS['logs']}/p3t5_checkpoint_dotplot.log",
    benchmark:
        f"{PATHS['benchmarks']}/p3t5_checkpoint_dotplot.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/checkpoint_dotplot.py"


rule p3t6_checkpoint_explorer:
    """Export the app-tier checkpoint explorer to self-contained WASM HTML (P3-T6).

    Every path the notebook reads is declared here. An undeclared read is a
    silent provenance hole (PROJECT_PLAN section A.5).

    The notebook is staged into a build directory beside a `public/` folder
    holding the six pipeline outputs, because `marimo export html-wasm` copies
    `public/` into the export and `mo.notebook_location()` resolves to the
    served URL in the browser. That is what lets the exported app fetch its data
    over HTTP, which is what keeps `mo.ui.table` interactive: bundling the
    frames through `cache_cells` instead would cache each table cell's output as
    an `UnhashableStub` and render stub text in place of the tables (ADR 0010).

    THE SLIDER IS A SENSITIVITY DISPLAY, NOT A THRESHOLD CONTROL (ADR 0019).
    The floor of record stays 0.5 (ADR 0016 section 1) and this rule changes no
    config: the app defaults to it, reproduces p3t5_checkpoint_dotplot's figure
    there, and labels itself the moment it leaves it. Nothing the slider reaches
    is a result. qc.detection_background_multiple is deliberately NOT exposed --
    detection_rate is already computed at it, and recomputing detection in a
    notebook is forbidden by hard constraint 3 and by p3t2_checkpoint_detection
    owning that rule.

    NO MODEL OUTPUT REACHES THE DOT PLOT. The carrier and shift estimates are
    inputs, but they render in their own section below a rule, because the plot
    is DETECTION and they are EXPRESSION, and the two move opposite ways under
    the same background gradient (ADR 0018). Adjacency satisfies PROJECT_PLAN
    section 6's "app over the P3-T3/T4 model outputs"; superposition is the
    conflation that already cost this phase one inverted argument.

    Staged under results/ rather than next to the notebook so the source tree
    stays free of generated files.
    """
    input:
        notebook="notebooks/apps/checkpoint_explorer.py",
        detection=f"{PATHS['tables']}/checkpoint_detection.tsv",
        summary=f"{PATHS['tables']}/checkpoint_detection_summary.json",
        carrier=f"{PATHS['tables']}/checkpoint_carrier_models.tsv",
        carrier_exploratory=f"{PATHS['tables']}/checkpoint_carrier_exploratory.tsv",
        shift=f"{PATHS['tables']}/checkpoint_shift_models.tsv",
        shift_exploratory=f"{PATHS['tables']}/checkpoint_shift_exploratory.tsv",
    output:
        app=report(
            directory(f"{PATHS['reports']}/checkpoint_explorer"),
            htmlindex="index.html",
            caption="../report/checkpoint_explorer.rst",
            category="Interactive",
            labels={"task": "P3-T6", "app": "checkpoint explorer"},
        ),
    params:
        build=f"{PATHS['interim']}/checkpoint_explorer_build",
    log:
        f"{PATHS['logs']}/p3t6_checkpoint_explorer.log",
    benchmark:
        f"{PATHS['benchmarks']}/p3t6_checkpoint_explorer.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    shell:
        "rm -rf {params.build} && mkdir -p {params.build}/public && "
        "cp {input.notebook} {params.build}/checkpoint_explorer.py && "
        "cp {input.detection} {input.summary} {input.carrier} "
        "{input.carrier_exploratory} {input.shift} {input.shift_exploratory} "
        "{params.build}/public/ && "
        "marimo export html-wasm --execute --mode run -f "
        "{params.build}/checkpoint_explorer.py -o {output.app} > {log} 2>&1"


TARGETS_CHECKPOINTS = [
    f"{PATHS['tables']}/checkpoint_membership.tsv",
    f"{PATHS['tables']}/checkpoint_membership_summary.json",
    f"{PATHS['tables']}/checkpoint_detection.tsv",
    f"{PATHS['tables']}/checkpoint_expression.tsv",
    f"{PATHS['tables']}/checkpoint_detection_summary.json",
    f"{PATHS['tables']}/checkpoint_carrier_models.tsv",
    f"{PATHS['tables']}/checkpoint_carrier_exploratory.tsv",
    f"{PATHS['tables']}/checkpoint_carrier_summary.json",
    f"{PATHS['tables']}/checkpoint_carrier_crosscheck.json",
    f"{PATHS['tables']}/checkpoint_shift_models.tsv",
    f"{PATHS['tables']}/checkpoint_shift_exploratory.tsv",
    f"{PATHS['tables']}/checkpoint_shift_summary.json",
    f"{PATHS['tables']}/checkpoint_shift_crosscheck.json",
    f"{PATHS['tables']}/checkpoint_paired_deltas.tsv",
    f"{PATHS['tables']}/checkpoint_paired_concordance.tsv",
    f"{PATHS['tables']}/checkpoint_paired_summary.json",
    f"{PATHS['figures']}/checkpoint_dotplot.png",
    f"{PATHS['reports']}/checkpoint_explorer",
]

if EXTERNAL_VALIDATION is not None:
    TARGETS_CHECKPOINTS += [
        f"{PATHS['resources']}/supplementary_checksums.sha256",
        f"{PATHS['resources']}/supplementary_provenance.tsv",
        f"{PATHS['tables']}/external_validation_matrix.tsv",
        f"{PATHS['tables']}/aoi_quantisation.tsv",
        f"{PATHS['tables']}/checkpoint_count_quantisation.tsv",
        f"{PATHS['tables']}/published_deg_reproduction.tsv",
        f"{PATHS['tables']}/published_fibrosis_contrast.tsv",
        f"{PATHS['tables']}/external_validation_summary.json",
    ]
