# ---------------------------------------------------------------------------
# Phase 5 — prognostic association (stretch). Gated on Q5, which resolved YES.
#
# Owner task: Phase 5 (see Markdowns/PROJECT_PLAN.md §6).
#
# GATE 5 IS "TIMEBOX RESPECTED", NOT "A RESULT WAS FOUND". PROJECT_PLAN calls
# this a stretch goal on a one-week timebox, the slip rule says cut it to protect
# Phase 6, and P5-T1's own acceptance criterion permits "a documented
# abandonment". The project is complete and coherent without this phase.
#
# Every threshold this file reads was pre-registered in ADR 0024 and committed
# before the first rule ran. Phase 5's central quantity is a median-split
# Kaplan-Meier at n = 13 (lung) and n = 8 (brain), and that is the easiest
# quantity in this project to manufacture after the fact: the split point, the
# endpoint, the censoring convention and the covariate set are each a degree of
# freedom, and each can be chosen to make a curve separate.
#
# Four rules bind every output in this file:
#
#   1. NO PHASE 5 RESULT MAY BE A HEADLINE CLAIM. With 12 events in lung and 7
#      in brain, this detects only very large hazard ratios. A null is
#      UNINFORMATIVE, NOT NEGATIVE — the same status the 1.1–1.3 SD power floor
#      gives a Phase 2 null (P0-T8) and ADR 0022 gives Phase 4's empty
#      nomination table.
#   2. A separation at 6 versus 7 patients DESCRIBES THIS COHORT. It is never an
#      inferential claim, however small its p-value.
#   3. Three of the six signatures failed their Phase 2 detection floor —
#      `exhaustion` 1/6 and `tls` 1/5 in brain, `myeloid_m1` 2/10 at BOTH sites.
#      Wherever they failed they are "not assessable", NEVER a prognostic null.
#      ADR 0008's rule is not suspended by a change of outcome variable: a score
#      built from genes at background measures background, whatever it is
#      regressed against.
#   4. Every brain claim states `TIME-B` n = 8 inline (hard constraint 8), and
#      every p-value carries n, effect size and a confidence interval (hard
#      constraint 7) — for a log-rank that means the hazard ratio and its CI
#      beside the test, never the p-value alone.
#
# THE JOIN IS A HARD GATE AND ITS FAILURE ACTION IS TO STOP THE PHASE. Two
# assertions from docs/data-provenance.md §Q5: every GEO patient number must
# appear in Supplementary Data 1, and the AOI-code numeric suffix must equal the
# patient number. NO FUZZY MATCHING, AT ALL — PROJECT_PLAN §6 is explicit ("do
# not spend three days on fuzzy ID matching for a stretch goal"). They are
# re-derived on every run; having been checked once by hand satisfies nothing.
#
# The cohort is the 16 distinct `TIME` patients (13 lung + 8 brain, 5 shared).
# NOT extended to L/LB to reach PROJECT_PLAN's "~35 patients": that 35 counts
# patients with a TUMOUR AOI, none of which carries a signature score, and
# reaching it would re-open the detection floor where the panel clears far less
# (P3-T2: L 2/9 genes, LB 3/9, against TIME-L 8/9). ADR 0024 §2; re-scoping is a
# stop-and-ask.
# ---------------------------------------------------------------------------

# Imported explicitly rather than inherited from another .smk's namespace.
# Includes share one global namespace, so `re` is in scope either way — but a
# rule whose wildcard constraint depends on include order is a reordering away
# from a confusing NameError (05_checkpoints.smk and 06_crosstalk.smk both make
# this note).
import re

SURVIVAL = config["survival"]


# --- P5-T1 the clinical table ----------------------------------------------
# Supplementary Data 1, pinned in config/clinical.yaml and PRE-REGISTERED in
# ADR 0024 §1.
#
# Fetched under ADR 0005's rules into a FIFTH resources subdirectory, alongside
# raw/ (P0-T2, GEO), reference/ (P2-T5, SpatialDecon), supplementary/ (P3-T2b,
# the source publication) and ligand_receptor/ (P4-T2, CellChatDB). Same
# two-rule split for the same reason — the digests cannot be params of the fetch
# rule, whose outputs are protected().
#
# NOT filed under supplementary/ with the other Springer files, and that is a
# decision rather than an oversight (ADR 0024 §1): those are an external
# COMPARATOR this project is checked against, this is analysis INPUT the project
# models. Reusing config/external_validation.yaml would also have changed
# p3t2b_external_validation's params hash and re-parsed a 45 MB workbook for no
# reason connected to Phase 5.
CLINICAL_DIR = f"{PATHS['resources']}/clinical"
CLINICAL_META = f"{PATHS['interim']}/clinical"

if CLINICAL is not None:
    CLINICAL_ARTIFACTS = CLINICAL["artifacts"]

    rule p5t1_fetch_clinical:
        """Download Supplementary Data 1 from the source publication (P5-T1).

        Reuses workflow/scripts/acquire_geo.py UNCHANGED, exactly as
        p2t5_fetch_reference, p3t2b_fetch_supplementary and p4t2a_fetch_lr_db
        do — nothing in it is GEO-specific. The wildcard must be named
        `artifact` because that script reads snakemake.wildcards.artifact by
        name; editing it instead would fire the code rerun-trigger on
        p0t2_fetch_geo and abort the DAG on its protected() outputs (ADR 0005).

        Springer static-content, not nature.com and not the PMC instance path:
        the article URL 303s to an IdP authorize endpoint and returns no file,
        and the PMC path returns an HTML interstitial rather than the xlsx
        (docs/data-provenance.md §Q5).
        """
        wildcard_constraints:
            artifact="|".join(
                re.escape(dest) for dest in sorted(CLINICAL_KEY_BY_DEST)
            ),
        output:
            artifact=protected(f"{CLINICAL_DIR}/{{artifact}}"),
            meta=f"{CLINICAL_META}/{{artifact}}.json",
        params:
            url=lambda w: (
                f"{CLINICAL['base_url']}/"
                f"{CLINICAL_ARTIFACTS[CLINICAL_KEY_BY_DEST[w.artifact]]['remote']}"
            ),
            key=lambda w: CLINICAL_KEY_BY_DEST[w.artifact],
            stability=lambda w: CLINICAL_ARTIFACTS[
                CLINICAL_KEY_BY_DEST[w.artifact]
            ]["stability"],
        log:
            f"{PATHS['logs']}/p5t1_fetch_clinical_{{artifact}}.log",
        benchmark:
            f"{PATHS['benchmarks']}/p5t1_fetch_clinical_{{artifact}}.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/acquire_geo.py"

    rule p5t1_record_clinical_provenance:
        """Verify the clinical table against its pinned digest (P5-T1).

        Reuses workflow/scripts/record_provenance.py UNCHANGED — it resolves
        paths from its inputs rather than a hardcoded root, so it writes any
        pair of checksum/provenance outputs.

        `fixed`, so a mismatch is a hard failure rather than an observation.
        This one file is the sole source of every Phase 5 outcome; letting it
        drift silently would invalidate the phase without any rule noticing.
        """
        input:
            artifacts=[
                f"{CLINICAL_DIR}/{dest}" for dest in sorted(CLINICAL_KEY_BY_DEST)
            ],
            meta=[
                f"{CLINICAL_META}/{dest}.json"
                for dest in sorted(CLINICAL_KEY_BY_DEST)
            ],
        output:
            checksums=f"{PATHS['resources']}/clinical_checksums.sha256",
            provenance=f"{PATHS['resources']}/clinical_provenance.tsv",
        params:
            artifacts=CLINICAL_ARTIFACTS,
        log:
            f"{PATHS['logs']}/p5t1_record_clinical_provenance.log",
        benchmark:
            f"{PATHS['benchmarks']}/p5t1_record_clinical_provenance.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/record_provenance.py"


    rule p5t1_clinical_join:
        """Parse Supplementary Data 1 and join it to the GEO design (P5-T1).

        THE HARD GATE. PROJECT_PLAN §6: "if it doesn't join cleanly on patient
        ID, stop the phase. Do not spend three days on fuzzy ID matching for a
        stretch goal." Its acceptance criterion is "a joined table **or** a
        documented abandonment", and both are legitimate outcomes.

        So the script has NO FUZZY MATCHING, NO FALLBACK PATH and no second
        threshold — the shape `lr_nominations.py` has, for the reason ADR 0022
        gives: an honest failure has to be REACHABLE rather than escapable, or
        the gate is decorative. Both §Q5 assertions raise rather than warn.

        Four quirks of the deposited file are handled explicitly, each measured
        rather than assumed (ADR 0024 §3-§4):

          1. Header cells carry EMBEDDED NEWLINES, so docs/data-provenance.md
             §Q5's column list is normalised and must not be matched verbatim.
          2. The sheet ends in a three-line legend an unbounded read ingests as
             patients with a null ID; the manifest bounds the data rows.
          3. The two endpoint columns spell censoring DIFFERENTLY — `Alive` in
             lung, an EMPTY CELL in brain. The two sets are asserted identical,
             which is what licenses reading a blank as censored rather than
             missing.
          4. Missingness has five spellings and `Gender` carries a lowercase
             `m`; a token outside the declared set is a hard failure.

        It joins and asserts. It fits NOTHING — P5-T2 aggregates, P5-T3 models —
        and it drops no patient: the nine supplementary rows without expression
        data are recorded as `in_geo = False` (flag-don't-drop).

        The .xlsx reader is LIFTED from external_validation.py rather than
        imported, the way checkpoint_detection.py lifts its detection rule from
        score_signatures.py. This project has no cross-script imports, Snakemake
        script-mode sibling imports are not reliable across versions, and
        editing external_validation.py would fire its code rerun-trigger.
        `openpyxl` stays out of py-analysis.yaml (ADR 0005, ADR 0013, ADR 0023).
        """
        input:
            clinical=[
                f"{CLINICAL_DIR}/{dest}" for dest in sorted(CLINICAL_KEY_BY_DEST)
            ],
            checksums=f"{PATHS['resources']}/clinical_checksums.sha256",
            samples=config["samples"]["tsv"],
        output:
            clinical=f"{PATHS['tables']}/clinical_cohort.tsv",
            summary=f"{PATHS['tables']}/clinical_join_summary.json",
        params:
            clinical=CLINICAL,
            survival=SURVIVAL,
        log:
            f"{PATHS['logs']}/p5t1_clinical_join.log",
        benchmark:
            f"{PATHS['benchmarks']}/p5t1_clinical_join.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/clinical_join.py"


    rule p5t2_patient_scores:
        """Collapse Phase 2's per-AOI signature scores to one per patient (P5-T2).

        PROJECT_PLAN §6 P5-T2 says "specify the aggregation rule". It is
        specified in `config.yaml → survival.aggregation` and pre-registered in
        ADR 0024, so this rule IMPLEMENTS it rather than choosing it:

          mean of the per-AOI scores  — matching ADR 0021 §1's duplicate_rule,
                                        so Phase 5 invents no second convention
          ssgsea PRIMARY              — zscore is a declared sensitivity
                                        (ADR 0012's shape); Phase 2 treated
                                        their AGREEMENT as the robustness
                                        evidence and neither is the other's
                                        tiebreaker

        The collapse touches 2 of 21 patient-sites — P12 and P24, both `TIME-L`,
        nothing in the brain set — so it cannot be load-bearing, and that is
        asserted rather than assumed.

        Output is the ANALYSIS-READY table P5-T3 fits directly: score, coverage
        verdict and outcome in one row, so no joining happens inside a modelling
        script.

        IT DOES NOT FILTER ON `assessable`. Three of the six signatures failed
        their Phase 2 coverage floor (`exhaustion` 1/6 and `tls` 1/5 in brain,
        `myeloid_m1` 2/10 at BOTH sites). They are carried with the flag
        attached, because which signatures are measurable is a Phase 2 RESULT
        and dropping them here would hide it. Wherever the flag is False the
        finding is "not assessable", NEVER a prognostic null (ADR 0008) — a
        score built from genes at background measures background, whatever it is
        regressed against.

        THE CENSORED PATIENT (ADR 0025). A censored patient carries an event
        status and NO TIME: `Alive` is a status, not a duration, the brain cell
        is blank, and Supplementary Data 1 has no last-contact date. **Patient
        35 is the one in the Phase 5 cohort and it is in BOTH arms.** It is kept
        in the table with `has_followup = False`, recorded in
        `survival_excluded.tsv` with a reason, and NEVER imputed — administrative
        censoring at the cohort maximum would invent an observation and place it
        exactly where a median split is most sensitive.

        So **12 lung and 7 brain patients enter a fit**, splits 6/6 and 3/4, and
        the brain arm lands EXACTLY on `survival.model.min_arm_size`. The rule
        fails rather than proceeds if it ever drops below.

        It fits, splits and tests NOTHING — P5-T3 owns all three.
        """
        input:
            scores=f"{PATHS['tables']}/signature_scores.tsv",
            coverage=f"{PATHS['tables']}/signature_coverage.tsv",
            clinical=f"{PATHS['tables']}/clinical_cohort.tsv",
        output:
            scores=f"{PATHS['tables']}/survival_patient_scores.tsv",
            excluded=f"{PATHS['tables']}/survival_excluded.tsv",
            summary=f"{PATHS['tables']}/survival_patient_scores_summary.json",
        params:
            survival=SURVIVAL,
        log:
            f"{PATHS['logs']}/p5t2_patient_scores.log",
        benchmark:
            f"{PATHS['benchmarks']}/p5t2_patient_scores.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/patient_scores.py"


    rule p5t3_survival_km:
        """Median split, Kaplan-Meier and log-rank, within arm (P5-T3).

        PROJECT_PLAN §6 P5-T3: "Median split per signature. Report n per arm; at
        ~35 patients this is descriptive." The real n is 12 lung and 7 brain
        entering a fit (ADR 0025), so it is MORE descriptive than the plan
        assumed, and this rule's job is to produce that description without
        letting it read as an inference.

        THE SPLIT IS `score > median`, COMPUTED WITHIN THE ARM, AND THE TIE SIDE
        IS DECLARED — at n = 7 the median IS an observation and which side it
        falls on changes the arm sizes. That gives 6/6 and 3/4, which ADR 0025 §3
        recorded before any curve existed, and which this rule ASSERTS rather
        than discovers. Not a tertile, not an optimal cutpoint, not a
        maximally-selected rank statistic (ADR 0024 §6).

        TWO ESTIMATORS, AND THEIR AGREEMENT IS THE CHECK. The log-rank is
        `sksurv.compare.compare_survival` — the `scikit-survival` skill
        PROJECT_PLAN §6 names — and the hazard ratio, its CI and a second
        p-value come from statsmodels PHReg, because hard constraint 7 forbids a
        p-value without n, effect size and a confidence interval, and a bare
        log-rank p is exactly that. `lifelines` is deliberately unused: it is
        absent from py-analysis.yaml and adding it would fire the software-env
        trigger on p0t2_fetch_geo's protected() outputs (ADR 0023).

        NO q-VALUES, BY DECISION (ADR 0026 §3). The family spans this study AND
        TCGA (P5-T4), so it cannot be corrected until both are fitted; P5-T5
        corrects ONCE over the union. ADR 0024 §8's declared family of 12 was
        arithmetically wrong — it counted cohorts and forgot this cohort has two
        arms — and ADR 0026 corrects it to the ASSESSABLE cells, 8 here plus
        TCGA's. The rule asserts no `q_` column reaches any output, so the
        convention cannot be restored by reflex.

        FOUR OF TWELVE CELLS ARE NOT FITTED, AND THAT IS PART OF THE RESULT.
        `exhaustion` and `tls` in brain and `myeloid_m1` at BOTH sites sit below
        Phase 2's coverage floor. They get a `not_assessable` row with a note and
        no estimate — the checkpoint_carrier_models.tsv shape — and are HATCHED
        on the figure rather than omitted, so a reader can tell "not tested" from
        "tested, null". ADR 0008: a score built from genes at background measures
        background, whatever it is regressed against.

        A declared EXPLORATORY table fits all twelve with no FDR, every row
        labelled (ADR 0018's shape). No claim may rest on it; it exists so the
        restriction is visible as a choice rather than as an absence.
        """
        input:
            scores=f"{PATHS['tables']}/survival_patient_scores.tsv",
        output:
            models=f"{PATHS['tables']}/survival_km_models.tsv",
            exploratory=f"{PATHS['tables']}/survival_km_exploratory.tsv",
            curves=f"{PATHS['tables']}/survival_km_curves.tsv",
            summary=f"{PATHS['tables']}/survival_km_summary.json",
            figure=report(
                f"{PATHS['figures']}/survival_km.png",
                caption="../report/survival_km.rst",
                category="Phase 5 — prognostic association",
                labels={"task": "P5-T3", "figure": "Kaplan-Meier by median signature split"},
            ),
        params:
            survival=SURVIVAL,
        log:
            f"{PATHS['logs']}/p5t3_survival_km.log",
        benchmark:
            f"{PATHS['benchmarks']}/p5t3_survival_km.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/survival_km.py"


# --- P5-T4 the TCGA LUAD external cohort ------------------------------------
# PRE-REGISTERED in ADR 0027. Fetched under ADR 0005's rules into a SIXTH
# resources subdirectory, after raw/ (P0-T2), reference/ (P2-T5), supplementary/
# (P3-T2b), ligand_receptor/ (P4-T2) and clinical/ (P5-T1). Same two-rule split
# for the same reason — the digests cannot be params of the fetch rule, whose
# outputs are protected().
#
# UCSC Xena, NOT cBioPortal. ADR 0024 §1 preferred cBioPortal because the source
# paper's own TCGA analyses came from there, and pre-authorised this fallback.
# Measured 2026-09-15 the bulk files are gone — 403 on the datahub, 404 on the
# LFS mirror — and only the live REST API answers, which serves per-study JSON
# and so cannot be digest-pinned under ADR 0005.
TCGA_DIR = f"{PATHS['resources']}/tcga"
TCGA_META = f"{PATHS['interim']}/tcga"

if TCGA is not None:
    TCGA_ARTIFACTS = TCGA["artifacts"]

    rule p5t4a_fetch_tcga:
        """Download one TCGA LUAD file from UCSC Xena (P5-T4).

        Reuses workflow/scripts/acquire_geo.py UNCHANGED, as every other
        acquisition rule in this project does. The wildcard must be named
        `artifact` because that script reads snakemake.wildcards.artifact by
        name; editing it instead would fire the code rerun-trigger on
        p0t2_fetch_geo and abort the DAG on its protected() outputs (ADR 0005).

        The remote key contains a PERCENT-ENCODED SLASH — Xena names the object
        "TCGA.LUAD.sampleMap/HiSeqV2.gz" and the %2F is part of the path. It must
        not be normalised away or the request 404s. `dest` is a bare filename by
        schema, so the local name deliberately differs from the remote key.
        """
        wildcard_constraints:
            artifact="|".join(
                re.escape(dest) for dest in sorted(TCGA_KEY_BY_DEST)
            ),
        output:
            artifact=protected(f"{TCGA_DIR}/{{artifact}}"),
            meta=f"{TCGA_META}/{{artifact}}.json",
        params:
            url=lambda w: (
                f"{TCGA['base_url']}/"
                f"{TCGA_ARTIFACTS[TCGA_KEY_BY_DEST[w.artifact]]['remote']}"
            ),
            key=lambda w: TCGA_KEY_BY_DEST[w.artifact],
            stability=lambda w: TCGA_ARTIFACTS[
                TCGA_KEY_BY_DEST[w.artifact]
            ]["stability"],
        log:
            f"{PATHS['logs']}/p5t4a_fetch_tcga_{{artifact}}.log",
        benchmark:
            f"{PATHS['benchmarks']}/p5t4a_fetch_tcga_{{artifact}}.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/acquire_geo.py"

    rule p5t4b_record_tcga_provenance:
        """Verify the TCGA files against their pinned digests (P5-T4).

        Reuses workflow/scripts/record_provenance.py UNCHANGED. Both artifacts
        are `fixed` — static since 2021 with stable ETags — so a mismatch is a
        hard failure rather than an observation. An external comparator that can
        drift silently is not a comparator.
        """
        input:
            artifacts=[
                f"{TCGA_DIR}/{dest}" for dest in sorted(TCGA_KEY_BY_DEST)
            ],
            meta=[
                f"{TCGA_META}/{dest}.json" for dest in sorted(TCGA_KEY_BY_DEST)
            ],
        output:
            checksums=f"{PATHS['resources']}/tcga_checksums.sha256",
            provenance=f"{PATHS['resources']}/tcga_provenance.tsv",
        params:
            artifacts=TCGA_ARTIFACTS,
        log:
            f"{PATHS['logs']}/p5t4b_record_tcga_provenance.log",
        benchmark:
            f"{PATHS['benchmarks']}/p5t4b_record_tcga_provenance.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/record_provenance.py"


    rule p5t4c_tcga_survival:
        """Score and fit the TCGA LUAD external cohort (P5-T4).

        P5-T3 found no immune signature stratified survival here. PROJECT_PLAN
        §6 P5-T4 exists for that outcome: "a negative here mostly tells you about
        your n, not about the biology." At ~500 patients this cohort separates a
        power limit from an absent effect — the one thing this project's own n
        cannot do. A NULL IN TCGA MEANS SOMETHING A NULL HERE DOES NOT.

        Reuses gseapy.ssgsea exactly as score_signatures.py calls it (min_size
        explicit — gseapy defaults it to 15 and every set here is 4-13 genes) and
        the estimator pair from survival_km.py, including its EXACT wiring check:
        inverting the group indicator must negate the Cox coefficient.

        ALL SIX SIGNATURES ARE TESTED, so TCGA contributes 6 and the Phase 5
        family is 8 + 6 = 14 — which discharges ADR 0026 §1's deferral. Phase 2's
        detection floor is q3 above a multiple of NegProbe-WTX (ADR 0007), a
        property of the GEOMX ASSAY; TCGA has no negative-probe channel, so
        restricting to the five assessable in the GeoMx lung arm would import a
        limitation the external cohort does not have. `myeloid_m1` is labelled
        TCGA-only, no GeoMx comparator.

        What it computes instead is GENE PRESENCE — whether a signature's genes
        are in the matrix at all. PRESENCE IS NOT DETECTION, and the two are
        labelled separately so they are never conflated.

        NO q-VALUES (ADR 0026 §3): P5-T5 corrects once over all 14. Asserted.

        TWO LIMITATIONS THAT BOUND EVERY NUMBER IT PRODUCES (ADR 0027 §5):
          - THIS IS NOT THE SAME MEASUREMENT. GeoMx TIME-L is the PanCK-negative
            SEGMENT of an ROI; TCGA is WHOLE BULK TUMOUR. A disagreement between
            cohorts is not necessarily a disagreement about biology.
          - THERE IS NO BRAIN COMPARATOR. TCGA LUAD is primary lung, so this
            bounds the LUNG null only — nothing here speaks to TIME-B n = 8.

        The cohort is constructed, not taken as given: only sample type `01` is
        kept (the matrix carries `11` normals and `02` recurrences), and cases
        with missing or non-positive follow-up are excluded and counted, never
        imputed — ADR 0025's shape applied to an external cohort.
        """
        input:
            expression=f"{TCGA_DIR}/{TCGA_ARTIFACTS['expression']['dest']}",
            survival=f"{TCGA_DIR}/{TCGA_ARTIFACTS['survival']['dest']}",
            checksums=f"{PATHS['resources']}/tcga_checksums.sha256",
        output:
            scores=f"{PATHS['tables']}/tcga_patient_scores.tsv",
            presence=f"{PATHS['tables']}/tcga_gene_presence.tsv",
            models=f"{PATHS['tables']}/tcga_km_models.tsv",
            curves=f"{PATHS['tables']}/tcga_km_curves.tsv",
            summary=f"{PATHS['tables']}/tcga_survival_summary.json",
            figure=report(
                f"{PATHS['figures']}/tcga_km.png",
                caption="../report/tcga_km.rst",
                category="Phase 5 — prognostic association",
                labels={"task": "P5-T4", "figure": "TCGA LUAD external check"},
            ),
        params:
            tcga=TCGA,
            survival=SURVIVAL,
            signatures=SIGNATURES,
            scoring=config["contexture"]["scoring"],
            seed=SEED,
        log:
            f"{PATHS['logs']}/p5t4c_tcga_survival.log",
        benchmark:
            f"{PATHS['benchmarks']}/p5t4c_tcga_survival.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 4
        script:
            "../scripts/tcga_survival.py"


    rule p5t5_multiplicity:
        """Correct the Phase 5 family once, and account for everything (P5-T5).

        PROJECT_PLAN §6 P5-T5 in full: "Multiplicity honesty. You tested 5
        signatures x 2 cohorts. Say so, with correction." This is the ONLY place
        in Phase 5 where a q-value is computed, and only q is reportable — a raw
        p from P5-T3 or P5-T4 is not.

        THE FAMILY IS COUNTED, NEVER QUOTED. ADR 0026 exists because ADR 0024
        §8's family size was written down rather than counted and was wrong: it
        forgot that ADR 0024 §2 splits this study into two arms that are never
        pooled. An undercounted denominator makes every q TOO SMALL — the
        direction that manufactures significance. So this rule derives the
        family by counting rows in the two upstream tables and asserts the
        reconciliation.

        IN THE FAMILY: the 8 assessable GeoMx signature x arm cells (P5-T3) and
        the 6 TCGA signatures (P5-T4), primary scoring method. ONE family, not
        two — "does the same signature stratify an independent cohort" is a
        single question asked twice (ADR 0024 §8, ADR 0026 §2).

        OUT: `zscore` and TCGA's continuous Cox, which are declared sensitivities
        (ADR 0012) — a second estimator does not double the search. And the 4
        GeoMx cells below Phase 2's coverage floor, because ADR 0008 forbids
        reporting them as prognostic nulls at all, so correcting over them would
        put un-interpretable rows in the denominator and DEFLATE every q.

        TIES ARE THE FAILURE MODE AND THIS FAMILY HAS THEM. Phase 4 shipped a
        step-down that was correct-looking and wrong under ties, caught only by a
        monotonicity assertion (NEXT_STEPS lesson 2). Here the correction is
        computed TWICE — statsmodels and an independent hand-rolled step-up —
        and they must agree exactly, with monotonicity asserted directly.

        It re-fits NOTHING: every estimate, interval and raw p is lifted from the
        upstream tables unchanged.
        """
        input:
            geomx=f"{PATHS['tables']}/survival_km_models.tsv",
            tcga=f"{PATHS['tables']}/tcga_km_models.tsv",
        output:
            family=f"{PATHS['tables']}/survival_multiplicity.tsv",
            accounting=f"{PATHS['tables']}/survival_multiplicity_accounting.tsv",
            not_assessable=f"{PATHS['tables']}/survival_not_assessable.tsv",
            summary=f"{PATHS['tables']}/survival_multiplicity_summary.json",
        params:
            survival=SURVIVAL,
        log:
            f"{PATHS['logs']}/p5t5_multiplicity.log",
        benchmark:
            f"{PATHS['benchmarks']}/p5t5_multiplicity.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/survival_multiplicity.py"


# Phase 5 is complete: P5-T1..T5 all have rules. TARGETS_SURVIVAL stays limited to what exists: the
# switch in config.yaml is the user's INTENT, this list is what actually exists,
# and a phase contributes to `rule all` only when both agree.
#
# Populate the rules and this list in the same commit, as every other phase did.
TARGETS_SURVIVAL = []

if CLINICAL is not None:
    TARGETS_SURVIVAL += [
        f"{PATHS['resources']}/clinical_checksums.sha256",
        f"{PATHS['resources']}/clinical_provenance.tsv",
        f"{PATHS['tables']}/clinical_cohort.tsv",
        f"{PATHS['tables']}/clinical_join_summary.json",
        f"{PATHS['tables']}/survival_patient_scores.tsv",
        f"{PATHS['tables']}/survival_excluded.tsv",
        f"{PATHS['tables']}/survival_patient_scores_summary.json",
        f"{PATHS['tables']}/survival_km_models.tsv",
        f"{PATHS['tables']}/survival_km_exploratory.tsv",
        f"{PATHS['tables']}/survival_km_curves.tsv",
        f"{PATHS['tables']}/survival_km_summary.json",
        f"{PATHS['figures']}/survival_km.png",
    ]

if TCGA is not None:
    TARGETS_SURVIVAL += [
        f"{PATHS['resources']}/tcga_checksums.sha256",
        f"{PATHS['resources']}/tcga_provenance.tsv",
        f"{PATHS['tables']}/tcga_patient_scores.tsv",
        f"{PATHS['tables']}/tcga_gene_presence.tsv",
        f"{PATHS['tables']}/tcga_km_models.tsv",
        f"{PATHS['tables']}/tcga_km_curves.tsv",
        f"{PATHS['tables']}/tcga_survival_summary.json",
        f"{PATHS['figures']}/tcga_km.png",
        f"{PATHS['tables']}/survival_multiplicity.tsv",
        f"{PATHS['tables']}/survival_multiplicity_accounting.tsv",
        f"{PATHS['tables']}/survival_not_assessable.tsv",
        f"{PATHS['tables']}/survival_multiplicity_summary.json",
    ]
