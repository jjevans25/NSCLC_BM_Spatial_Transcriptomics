# ---------------------------------------------------------------------------
# Phase 4 — inferred crosstalk between adjacent compartments (A5).
#
# Owner task: Phase 4 (see Markdowns/PROJECT_PLAN.md §6).
#
# NEVER "colocalisation" and never "spatially adjacent" (hard constraint 6).
# There are no coordinates. The phrase is "inferred crosstalk between adjacent
# compartments", everywhere, no exceptions — and p4t7_language_audit enforces it
# over the repository rather than trusting anyone to remember.
#
# A5 IS EXPLORATORY, NOT A HEADLINE AIM (ADR 0014). Everything here is still
# computed — the demotion is about what may be CLAIMED, not about what is built.
# Three rules bind every output in this file:
#
#   1. Every nomination states the per-site detection status of BOTH partners.
#      A rho without it is not reportable.
#   2. The anticipated null is an ASSAY-SENSITIVITY LIMIT, never evidence that
#      the crosstalk is absent. ADR 0014 demoted A5 because the ligand side is
#      largely below background — CXCL10, CXCL11, IL1B, TNF, IL12B, CD80, CD86,
#      NOS2, IL10, CCL22, CCL19, CCL21 — and a correlation that cannot be
#      computed is not a correlation of zero.
#   3. Every brain claim states TIME-B n = 8 inline (hard constraint 8).
#
# GATE 4 TURNS ON p4t4_null_calibration, NOT ON p4t5_nominations. The pass
# condition is the empirical FDR being computed and honoured: a surviving list
# and an empty list both pass, a ranked list without a null does not. With ~1000
# pairs at n = 13, |rho| > 0.7 arises by chance.
#
# Every threshold this file reads was pre-registered in ADR 0021 and committed
# before the first rule ran. A ranked correlation is the easiest quantity in
# this project to manufacture by choosing a threshold late, which is why
# PROJECT_PLAN §6's own P4-T1 acceptance criterion is "ADR committed".
#
# Adjacencies are L <-> TIME-L and LB <-> TIME-B, keyed on aoi_code and never on
# compartment (`compartment` is degenerate: L, LB and mLN are all `tumour`).
# TBME is absent despite being the glial side of a brain adjacency: 20 of 20 of
# its AOIs sit in one DSP run, so batch is inseparable from biology there
# (ADR 0009 §3). ADR 0014 §3 records it; re-scoping is a stop-and-ask.
# ---------------------------------------------------------------------------

# Imported explicitly rather than inherited from another .smk's namespace.
# Includes share one global namespace, so `re` is in scope either way — but
# a rule whose wildcard constraint depends on include order is a reordering
# away from a confusing NameError (05_checkpoints.smk makes the same note).
import re

CROSSTALK = config["crosstalk"]


rule p4t1_build_pairs:
    """The two paired adjacency sets, and how duplicate AOIs collapse (P4-T1).

    Lung `L` <-> `TIME-L` across 13 patients; brain `LB` <-> `TIME-B` across
    **8**. Both counts are ASSERTED against crosstalk.pairing.expected_n_patients
    and independently against design_matrix.tsv — if either moves, the design
    changed and the phase stops rather than adapting to it.

    P12 and P24 each contribute two `TIME-L` AOIs and one `L`; nothing in the
    brain set is duplicated, and P15's duplicate is `TBME`, which Phase 4 does
    not pair. How those two collapse is the one thing ADR 0021 §1 decides:

      primary      mean of log2(q3 + 1) across the duplicates
      sensitivity  the AOI with the highest gene_detection_rate

    Both tables are written and neither is a tiebreaker for the other
    (ADR 0012's shape). The choice moves 2 of 13 lung patients and 0 of 8 brain
    ones, so it cannot be load-bearing — and the sensitivity exists so that is
    checkable rather than asserted.

    No AOI is dropped. This is flag-don't-drop applied to a collapse rather than
    an exclusion, and qc_excluded.tsv gains nothing from Phase 4.
    """
    input:
        h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
        qc=f"{PATHS['tables']}/qc_metrics.tsv",
        design=f"{PATHS['tables']}/design_matrix.tsv",
    output:
        pairs=f"{PATHS['tables']}/crosstalk_pairs.tsv",
        sensitivity=f"{PATHS['tables']}/crosstalk_pairs_sensitivity.tsv",
        summary=f"{PATHS['tables']}/crosstalk_pairs_summary.json",
    params:
        adjacencies=CROSSTALK["adjacencies"],
        pairing=CROSSTALK["pairing"],
    log:
        f"{PATHS['logs']}/p4t1_build_pairs.log",
    benchmark:
        f"{PATHS['benchmarks']}/p4t1_build_pairs.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/crosstalk_pairs.py"


# --- P4-T2 the ligand–receptor database -------------------------------------
# CellChatDB v2 (jinworks), pinned at a commit, PRE-REGISTERED in ADR 0021 §2.
# Which interactions exist is Phase 4's panel, and a panel chosen after seeing
# which pairs correlate is not a panel.
#
# Fetched under ADR 0005's rules into a FOURTH resources subdirectory, alongside
# raw/ (P0-T2, GEO), reference/ (P2-T5, SpatialDecon) and supplementary/
# (P3-T2b, the source publication). Same two-rule split for the same reason —
# the digests cannot be params of the fetch rule, whose outputs are protected().
LR_DIR = f"{PATHS['resources']}/ligand_receptor"
LR_META = f"{PATHS['interim']}/ligand_receptor"

if LIGAND_RECEPTOR is not None:
    LR_ARTIFACTS = LIGAND_RECEPTOR["artifacts"]

    rule p4t2a_fetch_lr_db:
        """Download the pinned CellChatDB release (P4-T2).

        Reuses workflow/scripts/acquire_geo.py UNCHANGED, exactly as
        p2t5_fetch_reference and p3t2b_fetch_supplementary do — nothing in it is
        GEO-specific. The wildcard must be named `artifact` because that script
        reads snakemake.wildcards.artifact by name; editing it instead would
        fire the code rerun-trigger on p0t2_fetch_geo and abort the DAG on its
        protected() outputs (ADR 0005).

        base_url pins a COMMIT, enforced by pattern in
        ligand_receptor.schema.yaml. A branch ref would silently change the
        database under a rerun and the pinned digest would become a failure
        rather than a control.
        """
        wildcard_constraints:
            artifact="|".join(re.escape(dest) for dest in sorted(LR_KEY_BY_DEST)),
        output:
            artifact=protected(f"{LR_DIR}/{{artifact}}"),
            meta=f"{LR_META}/{{artifact}}.json",
        params:
            url=lambda w: (
                f"{LIGAND_RECEPTOR['base_url']}/"
                f"{LR_ARTIFACTS[LR_KEY_BY_DEST[w.artifact]]['remote']}"
            ),
            key=lambda w: LR_KEY_BY_DEST[w.artifact],
            stability=lambda w: LR_ARTIFACTS[LR_KEY_BY_DEST[w.artifact]]["stability"],
        log:
            f"{PATHS['logs']}/p4t2a_fetch_lr_db_{{artifact}}.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t2a_fetch_lr_db_{{artifact}}.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/acquire_geo.py"

    rule p4t2b_record_lr_provenance:
        """Verify the database against its pinned digest (P4-T2).

        Reuses workflow/scripts/record_provenance.py UNCHANGED. Every artifact
        here is `fixed`, enforced by enum in the schema: a commit-pinned raw
        file cannot move, so drift is a failure rather than an observation —
        unlike GEO's SOFT family file, which GEO regenerates server-side.
        """
        input:
            artifacts=[f"{LR_DIR}/{dest}" for dest in sorted(LR_KEY_BY_DEST)],
            meta=[f"{LR_META}/{dest}.json" for dest in sorted(LR_KEY_BY_DEST)],
        output:
            checksums=f"{PATHS['resources']}/ligand_receptor_checksums.sha256",
            provenance=f"{PATHS['resources']}/ligand_receptor_provenance.tsv",
        params:
            artifacts=LR_ARTIFACTS,
        log:
            f"{PATHS['logs']}/p4t2b_record_lr_provenance.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t2b_record_lr_provenance.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/record_provenance.py"

    rule p4t2c_export_lr_database:
        """Read the .rda and cross the R/Python boundary as TSV (P4-T2).

        The narrowest R rule in the project, and it exists for one reason:
        CellChatDB ships as a base-R .rda and nothing in the Python env reads
        one. `load()` is a BASE function, so this needs no CellChat package and
        no new dependency — it runs in the EXISTING, UNMODIFIED r-geomx env,
        which already reads .RData deconvolution references in deconvolve.R.
        Editing any env YAML would needlessly re-run Phase 2's and Phase 3's
        fits; editing py-analysis.yaml would abort the DAG outright.

        It does NO science: no symbol is resolved, no detection floor applied,
        nothing ranked. The only thing it drops is the pre-registered
        `exclude_annotations` classes — v2's 994 Non-protein Signaling
        interactions, whose ligands are metabolites with no transcript to
        measure. Counting those as "filtered out for low detection" would
        inflate the before/after gap P4-T2 exists to report honestly.

        Every value in config/ligand_receptor.yaml's `expect` block is asserted
        against the parsed file. A database that parses to different numbers
        means the pin moved or the parse is wrong, and both stop the phase.

        The handoff is plain TSV, never zellkonverter (ADR 0001). Nothing
        exported is a float — the database is symbols and labels — which is the
        one case where the round-trip is exact by construction rather than by
        the %.17g convention ADR 0013 §3 requires elsewhere.
        """
        input:
            rda=f"{LR_DIR}/{LR_ARTIFACTS['cellchatdb_human']['dest']}",
            checksums=f"{PATHS['resources']}/ligand_receptor_checksums.sha256",
            env_ok=f"{PATHS['interim']}/env_repair/r-geomx.ok",
        output:
            interactions=f"{PATHS['tables']}/crosstalk_lr_interactions.tsv",
            complexes=f"{PATHS['tables']}/crosstalk_lr_complexes.tsv",
        params:
            database=LIGAND_RECEPTOR["database"],
            version=LIGAND_RECEPTOR["version"],
            exclude_annotations=LIGAND_RECEPTOR["exclude_annotations"],
            expect=LIGAND_RECEPTOR["expect"],
        log:
            f"{PATHS['logs']}/p4t2c_export_lr_database.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t2c_export_lr_database.tsv"
        conda:
            "../envs/r-geomx.yaml"
        threads: 1
        script:
            "../scripts/export_lr_database.R"

    rule p4t2d_resolve_lr_pairs:
        """Resolve the database against the panel and apply the detection filter (P4-T2).

        THE BEFORE/AFTER GAP IS THE DELIVERABLE, not bookkeeping. It is the
        quantitative form of ADR 0014's claim that "secreted ligands and
        chemokines are systematically undetected", which is why it is reported
        BY INTERACTION CLASS — and why CellChatDB was chosen over CellPhoneDB at
        all (ADR 0021 §2). A count that lumps Secreted Signaling in with
        Cell-Cell Contact hides the finding inside the summary of it.

        Reuses the project's one detection rule — q3 > 2 x that AOI's
        NegProbe-WTX (ADR 0007) — lifted from checkpoint_detection.py rather
        than rewritten. Two detection rules in one project is one too many.

        Detection is computed over ALL AOIs of an aoi_code (L 30, TIME-L 15,
        LB 27, TIME-B 8), not over the paired subset, so Phase 4's coverage
        numbers compare directly with Phase 2's and Phase 3's — the only
        coverage numbers a reader has to calibrate against.

        Admission is evaluated PER DIRECTION (ADR 0021 §9). CellChatDB is
        directed, and a pair can be admissible with the ligand measured in the
        tumour compartment and not with it measured in the immune one, because
        the two partners are being asked about in different places.

        THREE TABLES, pre-registered in ADR 0021 §4:
          membership   every usable interaction, admitted or not, with
                       per-compartment detection for both partners. An
                       interaction that fails is a ROW SAYING SO, never an
                       absence.
          filtered     one-to-one, admitted. The only Phase 4 table that will
                       carry an empirical FDR.
          exploratory  complex, every subunit clearing the floor. NO FDR of any
                       kind, ever — the script asserts no test statistic can
                       reach it, the way checkpoint_paired_check.py does.
        """
        input:
            h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
            interactions=f"{PATHS['tables']}/crosstalk_lr_interactions.tsv",
            complexes=f"{PATHS['tables']}/crosstalk_lr_complexes.tsv",
            pairs=f"{PATHS['tables']}/crosstalk_pairs.tsv",
        output:
            membership=f"{PATHS['tables']}/crosstalk_lr_membership.tsv",
            filtered=f"{PATHS['tables']}/crosstalk_lr_filtered.tsv",
            exploratory=f"{PATHS['tables']}/crosstalk_lr_exploratory.tsv",
            summary=f"{PATHS['tables']}/crosstalk_lr_summary.json",
        params:
            adjacencies=CROSSTALK["adjacencies"],
            detection=CROSSTALK["detection"],
            background_multiple=config["qc"]["detection_background_multiple"],
            database=LIGAND_RECEPTOR["database"],
            version=LIGAND_RECEPTOR["version"],
            expect=LIGAND_RECEPTOR["expect"],
        log:
            f"{PATHS['logs']}/p4t2d_resolve_lr_pairs.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t2d_resolve_lr_pairs.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 2
        script:
            "../scripts/resolve_lr_pairs.py"


    rule p4t3a_crosstalk_expression:
        """The paired expression matrix, and the background it rests on (P4-T3).

        Collapses each patient's AOIs per ADR 0021 §1 — mean of log2(q3 + 1) for
        the primary, highest gene_detection_rate for the sensitivity — and
        writes both, so the claim that the choice is not load-bearing stays
        checkable rather than asserted.

        A gene counts as detected for a patient only if it clears the floor in
        EVERY constituent AOI. The conservative reading, chosen because
        n_detected_both is what tells a reader how much of a rho rests on values
        at background.

        IT ALSO MEASURES THE GRADIENT THAT COULD INVALIDATE THE PHASE. ADR 0018
        established that detection and expression move opposite ways under the
        same background gradient; Phases 2 and 3 met that as a GROUP-CONTRAST
        problem. Phase 4 correlates ACROSS PATIENTS, so the question here is a
        different one — does a patient's two paired AOIs SHARE their background?
        If they do, two genes that merely track background correlate across
        patients for no biological reason at all. P3-T4 measured the SITE
        gradient and found it negligible; that is NOT this gradient, and
        inheriting the verdict would assume exactly what needs checking.

        The rule REPORTS that number. It does not adjust for it and does not
        gate on it: adjusting would change the pre-registered specification, and
        gating would make a threshold out of a measurement taken to inform one.
        P4-T4's abundance-matched Null B is the pre-registered response and runs
        whatever the number turns out to be.
        """
        input:
            h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
            pairs=f"{PATHS['tables']}/crosstalk_pairs.tsv",
            pairs_sensitivity=f"{PATHS['tables']}/crosstalk_pairs_sensitivity.tsv",
            filtered=f"{PATHS['tables']}/crosstalk_lr_filtered.tsv",
            exploratory=f"{PATHS['tables']}/crosstalk_lr_exploratory.tsv",
            qc=f"{PATHS['tables']}/qc_metrics.tsv",
        output:
            expression=f"{PATHS['tables']}/crosstalk_expression.tsv",
            sensitivity=f"{PATHS['tables']}/crosstalk_expression_sensitivity.tsv",
            background=f"{PATHS['tables']}/crosstalk_background_check.tsv",
            summary=f"{PATHS['tables']}/crosstalk_expression_summary.json",
        params:
            adjacencies=CROSSTALK["adjacencies"],
            pairing=CROSSTALK["pairing"],
            detection=CROSSTALK["detection"],
            background_multiple=config["qc"]["detection_background_multiple"],
        log:
            f"{PATHS['logs']}/p4t3a_crosstalk_expression.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t3a_crosstalk_expression.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 2
        script:
            "../scripts/crosstalk_expression.py"

    rule p4t3b_lr_correlation:
        """Spearman, ligand against paired receptor, across patients (P4-T3).

        SPEARMAN, NEVER PEARSON — PROJECT_PLAN §6 P4-T3 is explicit, and the
        reason is the n: at 13 lung patients and 8 brain ones a Pearson
        coefficient is a coin flip on one outlier. The schema makes
        crosstalk.correlation.method an enum of one and the script re-asserts
        it, so a config edit cannot quietly change the estimator.

        BOTH DIRECTIONS, LABELLED (ADR 0021 §9). CellChatDB is directed, and
        "ligand on the tumour AOI, receptor on the paired immune AOI" is a
        different biological claim from its reverse. Neither is a control for
        the other.

        The primary keeps EVERY paired patient with n_detected_both as a column
        (ADR 0021 §5). Dropping patients whose values sit at background is
        expression-dependent selection — it removes the low values of BOTH
        partners preferentially, which induces a correlation rather than
        removing one — and it would make n vary per pair, so P4-T4's null would
        need recalibrating at every distinct n. The detected-both refit is a
        declared sensitivity with a delta column (ADR 0012's shape).

        NOTHING THIS RULE PRODUCES IS INTERPRETABLE ON ITS OWN. With ~1000
        interactions at n = 13, |rho| > 0.7 arises by chance, and PROJECT_PLAN
        §6 says the ranking is uninterpretable without the empirical null. This
        table is an INPUT to p4t4_null_calibration, which is what Gate 4 turns
        on.
        """
        input:
            expression=f"{PATHS['tables']}/crosstalk_expression.tsv",
            filtered=f"{PATHS['tables']}/crosstalk_lr_filtered.tsv",
            exploratory=f"{PATHS['tables']}/crosstalk_lr_exploratory.tsv",
        output:
            primary=f"{PATHS['tables']}/crosstalk_correlation.tsv",
            exploratory=f"{PATHS['tables']}/crosstalk_correlation_exploratory.tsv",
            summary=f"{PATHS['tables']}/crosstalk_correlation_summary.json",
        params:
            adjacencies=CROSSTALK["adjacencies"],
            correlation=CROSSTALK["correlation"],
        log:
            f"{PATHS['logs']}/p4t3b_lr_correlation.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t3b_lr_correlation.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 2
        script:
            "../scripts/lr_correlation.py"



    rule p4t4_null_calibration:
        """The empirical null and the FDR. THE TASK GATE 4 TURNS ON (P4-T4).

        PROJECT_PLAN §6 P4-T4: "With ~1000 LR pairs at n=13, some |rho|>0.7
        arises by chance. Permute patient labels to get an empirical null and
        report an empirical FDR. WITHOUT THIS, THE RANKING IS UNINTERPRETABLE."

        Gate 4's pass condition is this number being computed and HONOURED. A
        surviving list and an empty list both pass; a ranked list without a null
        does not. Nothing downstream may reach past it for a softer threshold.

        TWO NULLS (ADR 0021 §6), and only the first is the gate.

        NULL A — the gate. Permutes THE PAIRING, not the expression: which
        patient's immune AOI is matched to which patient's tumour AOI, WITHIN
        site, so both n's are preserved exactly (13 and 8). Each gene's values,
        distribution, abundance and detection are untouched; what is destroyed
        is precisely the claim the analysis makes, that the ligand and the
        receptor were measured in the same person.

        NULL B — the abundance-matched control. A nomination list dominated by
        VEGFA, TGFB1 and CD47 is a detection artefact waiting to be over-read:
        those are broadly expressed, and across 13 patients their correlation
        may track shared AOI quality rather than crosstalk. p4t3a MEASURES
        whether that quality is shared across a patient's paired AOIs; this
        control answers the question whatever that came out at, because a
        pre-registered control is not contingent on the confound turning out to
        be large. It QUALIFIES a nomination — not a second gate, no second FDR
        family, and a row clearing A but not B is reported as such, never
        removed.

        FAMILY: all primary direction-rows WITHIN SITE (ADR 0021 §6). A null
        built at n = 13 cannot be applied to the brain's n = 8, so the two sites
        are corrected separately. Both directions of an interaction sit in the
        SAME family — they are distinct hypotheses asked in one sweep of one
        database, and splitting them would make two families where one was
        asked.

        Seeds come from config["seed"], passed explicitly. No implicit RNG
        anywhere; p1t3_variance_partition is the in-repo precedent.

        Spearman is Pearson on ranks and a permutation only reorders values, so
        ranks are computed ONCE and a whole permutation is one elementwise
        multiply. The identity permutation is asserted to reproduce p4t3b's rho
        before any permuted value is trusted.
        """
        input:
            expression=f"{PATHS['tables']}/crosstalk_expression.tsv",
            correlation=f"{PATHS['tables']}/crosstalk_correlation.tsv",
        output:
            fdr=f"{PATHS['tables']}/crosstalk_fdr.tsv",
            null=f"{PATHS['tables']}/crosstalk_null.tsv",
            summary=f"{PATHS['tables']}/crosstalk_null_summary.json",
            figure=report(
                f"{PATHS['figures']}/crosstalk_null.png",
                caption="../report/crosstalk_null.rst",
                category="Phase 4 — inferred crosstalk",
                labels={"task": "P4-T4", "figure": "empirical null distribution"},
            ),
        params:
            seed=SEED,
            adjacencies=CROSSTALK["adjacencies"],
            null_calibration=CROSSTALK["null_calibration"],
        log:
            f"{PATHS['logs']}/p4t4_null_calibration.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t4_null_calibration.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 2
        script:
            "../scripts/lr_null_calibration.py"



    rule p4t5_nominations:
        """Nominate 5-10 pairs, framed as hypothesis generation (P4-T5).

        Top-ranked, detection-filtered, surviving the empirical FDR. Each row
        carries rho, the empirical FDR, the abundance-matched p, n,
        n_detected_both and THE PER-SITE DETECTION COUNTS OF BOTH PARTNERS —
        ADR 0014's binding term, not a nicety. A nomination without them is not
        reportable, and the script joins them from
        crosstalk_lr_membership.tsv and asserts the two sources agree.

        IF NOTHING SURVIVES, THAT IS THE DELIVERABLE. Gate 4 licenses it in
        terms: "no LR pair exceeded chance expectation at n=13" is an honest,
        useful result and a better outcome than a ranked list you can't defend.
        The script has NO fallback path and no second threshold to reach for.

        GATE 4 DOES NOT TURN ON THIS TABLE. It turns on p4t4_null_calibration —
        on the empirical FDR being computed and honoured. A surviving list and
        an empty list both pass. `max_nominations` truncates a long list; it
        never licenses a longer one.

        The rationale is LOOKED UP from the database's own pathway_name,
        annotation and evidence string, never generated. Prose about why a
        correlation makes biological sense, written after seeing that it
        correlated, is how a ranked list becomes a story.
        """
        input:
            fdr=f"{PATHS['tables']}/crosstalk_fdr.tsv",
            membership=f"{PATHS['tables']}/crosstalk_lr_membership.tsv",
        output:
            nominations=f"{PATHS['tables']}/crosstalk_nominations.tsv",
            summary=f"{PATHS['tables']}/crosstalk_nominations_summary.json",
        params:
            adjacencies=CROSSTALK["adjacencies"],
            null_calibration=CROSSTALK["null_calibration"],
            # PROJECT_PLAN §6 P4-T5, "Nominate 5-10 pairs". A rule param rather
            # than a config key, exactly as p3t4_checkpoint_paired_check carries
            # expected_immune_paired=5 from PROJECT_PLAN §2.3: this is the
            # PLAN'S number, not a threshold this phase invented, so it does not
            # belong in the pre-registered config block — and config.yaml is a
            # declared input of p0t7_assemble_h5ad, so a key here would rebuild
            # the .h5ad and re-run Phases 1-3 for a display cap.
            max_nominations=10,
        log:
            f"{PATHS['logs']}/p4t5_nominations.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t5_nominations.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/lr_nominations.py"



    rule p4t6_external_crossref:
        """Informal concordance against an external NSCLC Visium study (P4-T6).

        De Zuani et al. 2024, Nat Commun 15:4388 (PMID 38782901), pinned in
        config/ligand_receptor.yaml with the paper's OWN SENTENCE quoted,
        pre-registered before any concordance was computed.

        NOT VALIDATION, and the three reasons are structural rather than
        cautionary — they are carried on every output row, because a concordance
        table read without them is a replication claim:

          1. They used CellPhoneDB; this project uses CellChatDB. A pair missing
             from our panel is a DATABASE difference, never a biological one.
          2. They measured co-expression WITHIN A VISIUM SPOT, tumour section vs
             background section. This project measures a CROSS-PATIENT
             CORRELATION between two compartments and has NO COORDINATES AT ALL.
             Different quantities: agreement is encouraging, disagreement is
             evidence against neither.
          3. They profiled primary lung (LUAD/LUSC). THE BRAIN ARM HAS NO
             COMPARATOR, and TIME-B n = 8.

        Two comparators were known compromised before the script existed:
        LGALS9 is absent from the 18,694 measured genes (reported not_on_panel;
        LGALS9C is NOT substituted, and the script asserts it never appears —
        that would be a membership change made to keep a comparator alive), and
        TIGIT sits below the floor in the lung tumour compartment. VEGFA-NRP1 is
        the one genuinely informative check.

        `not_on_panel`, `not_in_database`, `below_detection_floor` and
        `not_recovered` are kept as four separate statuses so an assay limit
        cannot collapse into "did not replicate".
        """
        input:
            h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
            membership=f"{PATHS['tables']}/crosstalk_lr_membership.tsv",
            fdr=f"{PATHS['tables']}/crosstalk_fdr.tsv",
        output:
            concordance=f"{PATHS['tables']}/crosstalk_external_concordance.tsv",
            summary=f"{PATHS['tables']}/crosstalk_external_summary.json",
        params:
            external_comparator=LIGAND_RECEPTOR["external_comparator"],
            null_calibration=CROSSTALK["null_calibration"],
            adjacencies=CROSSTALK["adjacencies"],
        log:
            f"{PATHS['logs']}/p4t6_external_crossref.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t6_external_crossref.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/lr_external_crossref.py"

    rule p4t7_network_figure:
        """The Phase 4 deliverable figure: the inferred-crosstalk network (P4-T7).

        Nodes are genes coloured by the compartment they are measured in; edges
        are inferred ligand-receptor relationships coloured by rho. Laid out
        BIPARTITE — tumour-side left, immune-side right — because that is what
        the design is; a spring layout would invent a topology the data does not
        have.

        NEVER colocalisation and never spatially adjacent. There are no
        coordinates. An edge means the two compartments' expression correlated
        ACROSS PATIENTS — inference from a database plus a correlation, not a
        measured proximity (hard constraint 6).

        EDGE COLOUR IS RdBu_r, AND THAT IS A DELIBERATE REUSE RATHER THAN A
        BORROW. The project's colour language is fixed: RdBu_r means a signed
        effect (P2-T7), PRGn centred on the detection multiple means a
        background ratio (P3-T5, ADR 0016 §5). A Spearman rho IS a signed effect
        on [-1, 1] — the same quantity class — so speaking that language is
        consistent, and inventing a third would imply rho is a third kind of
        quantity. The caption states it.

        EDGE STYLE, NOT EDGE PRESENCE, CARRIES THE FDR. Solid clears it, dashed
        does not. Drawing only survivors would produce a blank panel for an
        outcome Gate 4 explicitly licenses — "no LR pair exceeded chance
        expectation at n=13" — and a blank panel reads as a broken rule rather
        than a finding.

        Complex interactions are absent: they carry no FDR of any kind
        (ADR 0021 §4), so they have no style to draw.
        """
        input:
            fdr=f"{PATHS['tables']}/crosstalk_fdr.tsv",
            nominations=f"{PATHS['tables']}/crosstalk_nominations.tsv",
        output:
            figure=report(
                f"{PATHS['figures']}/crosstalk_network.png",
                caption="../report/crosstalk_network.rst",
                category="Phase 4 — inferred crosstalk",
                labels={"task": "P4-T7", "figure": "inferred crosstalk network"},
            ),
            summary=f"{PATHS['tables']}/crosstalk_network_summary.json",
        params:
            adjacencies=CROSSTALK["adjacencies"],
            null_calibration=CROSSTALK["null_calibration"],
            # A display cap, not a threshold: it decides how many of the
            # strongest edges are legible on one canvas, and nothing about
            # which are real — that is the edge STYLE, from the pre-registered
            # FDR. Same standing as p4t5's max_nominations.
            top_n_edges=25,
        log:
            f"{PATHS['logs']}/p4t7_network_figure.log",
        benchmark:
            f"{PATHS['benchmarks']}/p4t7_network_figure.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/lr_network.py"


rule p4t7_language_audit:
    """Audit the repository for claims the assay cannot support (P4-T7).

    PROJECT_PLAN §6 P4-T7: "grep the entire repo for 'colocali' and 'spatial
    proximity' and fix every hit. You have no coordinates."

    A RULE rather than a command someone remembers to run. Hard constraint 5
    says every output is produced by a rule, and a check that lives in a handoff
    note survives exactly as long as the person who wrote the note. As a rule, a
    violation FAILS THE BUILD.

    Two families of forbidden claim:

      * Never a claim of colocalisation, spatial adjacency or spatial
        proximity — this is GeoMx DSP. An AOI is a region a pathologist placed;
        there is no spot grid and no coordinate. The supported phrase is
        "inferred crosstalk between adjacent compartments" (hard constraint 6).
      * Never a CD45+ or GFAP+ AOI — PanCK was the ONLY collection mask (Q2).
        CD45 and GFAP guided ROI PLACEMENT; nothing was collected on either, so
        a compartment label is never a cell-type label.

    The script assembles its patterns from fragments so it does not match its
    own source. Excluding it by name would be the obvious alternative and would
    create the one place in the repository where the forbidden phrase could be
    written without consequence.

    Deliberately NOT inside the `if LIGAND_RECEPTOR is not None` block: the
    language rule is the project's, not the database's, and it must still run on
    a checkout that has not fetched CellChatDB.
    """
    output:
        report=f"{PATHS['tables']}/crosstalk_language_audit.json",
    params:
        # Not read by the script — it re-derives the file list itself. This is
        # the RERUN TRIGGER: it changes whenever any scanned file changes, so
        # the audit cannot report a stale pass.
        tracked_digest=_tracked_text_digest(),
    log:
        f"{PATHS['logs']}/p4t7_language_audit.log",
    benchmark:
        f"{PATHS['benchmarks']}/p4t7_language_audit.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/language_audit.py"



TARGETS_CROSSTALK = [
    f"{PATHS['tables']}/crosstalk_language_audit.json",
    f"{PATHS['tables']}/crosstalk_pairs.tsv",
    f"{PATHS['tables']}/crosstalk_pairs_sensitivity.tsv",
    f"{PATHS['tables']}/crosstalk_pairs_summary.json",
]

if LIGAND_RECEPTOR is not None:
    TARGETS_CROSSTALK += [
        f"{PATHS['resources']}/ligand_receptor_checksums.sha256",
        f"{PATHS['resources']}/ligand_receptor_provenance.tsv",
        f"{PATHS['tables']}/crosstalk_lr_interactions.tsv",
        f"{PATHS['tables']}/crosstalk_lr_complexes.tsv",
        f"{PATHS['tables']}/crosstalk_lr_membership.tsv",
        f"{PATHS['tables']}/crosstalk_lr_filtered.tsv",
        f"{PATHS['tables']}/crosstalk_lr_exploratory.tsv",
        f"{PATHS['tables']}/crosstalk_lr_summary.json",
        f"{PATHS['tables']}/crosstalk_expression.tsv",
        f"{PATHS['tables']}/crosstalk_expression_sensitivity.tsv",
        f"{PATHS['tables']}/crosstalk_background_check.tsv",
        f"{PATHS['tables']}/crosstalk_expression_summary.json",
        f"{PATHS['tables']}/crosstalk_correlation.tsv",
        f"{PATHS['tables']}/crosstalk_correlation_exploratory.tsv",
        f"{PATHS['tables']}/crosstalk_correlation_summary.json",
        f"{PATHS['tables']}/crosstalk_fdr.tsv",
        f"{PATHS['tables']}/crosstalk_null.tsv",
        f"{PATHS['tables']}/crosstalk_null_summary.json",
        f"{PATHS['figures']}/crosstalk_null.png",
        f"{PATHS['tables']}/crosstalk_nominations.tsv",
        f"{PATHS['tables']}/crosstalk_nominations_summary.json",
        f"{PATHS['tables']}/crosstalk_external_concordance.tsv",
        f"{PATHS['tables']}/crosstalk_external_summary.json",
        f"{PATHS['figures']}/crosstalk_network.png",
        f"{PATHS['tables']}/crosstalk_network_summary.json",
    ]
