# ---------------------------------------------------------------------------
# Phase 6 — FAIR packaging: ontology terms, RO-Crate, citation metadata,
# report audit. The phase PROJECT_PLAN §6 calls "arguably the primary
# deliverable": the artifact someone else can actually run.
#
# Owner task: Phase 6 (see Markdowns/PROJECT_PLAN.md §6).
#
# FOUR THINGS BIND EVERY RULE IN THIS FILE.
#
#   1. NOTHING HERE PRODUCES A RESULT. Phase 6 packages what Phases 0-5 found;
#      it does not compute, revise or reinterpret any of it. A rule in this file
#      that changed a number would be a bug, and the clean-room comparison
#      (P6-T1) is what would catch it.
#
#   2. THE CRATE IS MANIFEST-DRIVEN, NOT A FILESYSTEM WALK (ADR 0030). Every
#      described entity is named in config/ro_crate.yaml and declared as a rule
#      input, so the DAG knows what the crate depends on and a described file
#      that no rule produces fails the build. A crate built by walking results/
#      would describe whatever happened to be on disk that afternoon.
#
#   3. A COMPARTMENT LABEL IS NEVER A CELL-TYPE LABEL, AND THIS IS THE PHASE
#      WHERE THAT STOPS BEING PROSE. PanCK was the only collection mask; CD45
#      and GFAP guided ROI placement only (Q2). p4t7_language_audit can catch
#      the phrase in a document; it cannot stop a JSON-LD graph from asserting
#      the claim in RDF to a machine that will not read it sceptically. So the
#      ontology mapping carries an explicit `relation`, compartments resolve as
#      `enriched_for`, and common.smk refuses the manifest at parse time if one
#      of them says `is_a`.
#
#   4. `snakemake --report` IS A COMMAND, NOT A RULE, and that is the one
#      sanctioned exception to hard constraint 5 (ADR 0030 §3). It consumes the
#      completed DAG, so it cannot be a node inside it. What CAN be a rule is
#      the check that the report will be worth opening — p6t4_report_audit —
#      and that is where the acceptance criterion lives.
# ---------------------------------------------------------------------------

METADATA_DIR = "metadata"


# --- P6-T2a the ontology mapping -------------------------------------------
# PROJECT_PLAN §7.2, five phases late. See the script docstring for why that is
# cheaper here than it sounds, and expensive in a way that is worth recording.
if ONTOLOGY_TERMS is not None:

    rule p6t2a_resolve_ontology_terms:
        """Resolve every metadata field against EBI OLS4 (P6-T2a, FAIR I2).

        §7.2: "Resolve each against EBI OLS4 — don't hand-type IDs from memory,
        that's exactly the error class the skill exists to prevent."

        THE MANIFEST IS A PIN, NOT A LOOKUP TABLE. Every CURIE in
        config/ontology_terms.yaml came from OLS4; this rule asks OLS4 again on
        every build and asserts the answer is unchanged. That is the ADR 0005
        contract applied to a vocabulary instead of a download, and it buys what
        a one-off curation cannot: ONTOLOGY DRIFT FAILS THE BUILD. A term that
        is obsoleted, merged, or relabelled comes back as a different id under
        the same label, and a curated TSV would absorb that silently forever.

        NO FUZZY MATCHING AND NO FALLBACK PATH. Exact primary label, restricted
        to the declared ontology, restricted to hits where it is the DEFINING
        ontology, and the resolved CURIE must equal the pinned one. Retries
        exist for a flaky network and never widen a match. An unresolved term
        stops the phase and is a stop-and-ask — what a field means is a
        scientific question, not an implementation detail.

        The network is a real dependency here, and that is not new: every fetch
        rule in Phases 0, 2, 3, 4 and 5 has one. What is new is that this one
        asks a question rather than downloading an answer.
        """
        input:
            manifest=config["fair"]["ontology_yaml"],
        output:
            tsv=f"{METADATA_DIR}/ontology-terms.tsv",
            summary=f"{METADATA_DIR}/ontology_resolution_summary.json",
        params:
            terms=ONTOLOGY_TERMS,
            ols=config["fair"]["ols"],
        log:
            f"{PATHS['logs']}/p6t2a_resolve_ontology_terms.log",
        benchmark:
            f"{PATHS['benchmarks']}/p6t2a_resolve_ontology_terms.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/resolve_ontology_terms.py"



# --- P6-T2 the RO-Crate -----------------------------------------------------
# The crate describes THIS repository: only ro-crate-metadata.json is written,
# every @id is a relative path into the working tree, and nothing is copied.
if RO_CRATE is not None:

    # The fetched inputs, resolved from the SIX manifests rather than globbed.
    # Every one of these is also a rule output, so naming them here makes the
    # crate wait for the fetch rather than describe an empty directory.
    RO_CRATE_MANIFEST_PARTS = []
    if RO_CRATE_MANIFEST_GROUPS:
        _fetched = []
        _fetched += [
            f"{PATHS['raw']}/{_spec['dest']}"
            for _spec in config["acquire"]["artifacts"].values()
        ]
        _fetched += [
            f"{PATHS['reference']}/{_spec['dest']}"
            for _spec in config["reference"]["artifacts"].values()
        ]
        if EXTERNAL_VALIDATION is not None:
            _fetched += [
                f"{PATHS['resources']}/supplementary/{_spec['dest']}"
                for _spec in EXTERNAL_VALIDATION["artifacts"].values()
            ]
        if LIGAND_RECEPTOR is not None:
            _fetched += [
                f"{PATHS['resources']}/ligand_receptor/{_spec['dest']}"
                for _spec in LIGAND_RECEPTOR["artifacts"].values()
            ]
        if CLINICAL is not None:
            _fetched += [
                f"{PATHS['resources']}/clinical/{_spec['dest']}"
                for _spec in CLINICAL["artifacts"].values()
            ]
        if TCGA is not None:
            _fetched += [
                f"{PATHS['resources']}/tcga/{_spec['dest']}"
                for _spec in TCGA["artifacts"].values()
            ]
        _fetched = sorted(set(_fetched))

        for _group in RO_CRATE_MANIFEST_GROUPS:
            if len(_fetched) < _group["expect_min"]:
                raise WorkflowError(
                    f"config ro_crate.parts[{_group['group']}]: the six fetch "
                    f"manifests declare {len(_fetched)} artifacts, fewer than "
                    f"the {_group['expect_min']} expected. A manifest has "
                    "emptied, or one is absent from this checkout."
                )
            RO_CRATE_MANIFEST_PARTS += [
                {"path": _path, "group": _group["group"]} for _path in _fetched
            ]


    # The derived half of the crate's parts, composed HERE rather than in
    # common.smk because this file is included last and is the only place where
    # every phase's TARGETS list exists.
    #
    # Taken from the targets rather than from a glob over results/, for the
    # reason the schema gives at length: a glob returns NOTHING on a clean
    # checkout — exactly what P6-T1 creates — so the crate would describe zero
    # figures and still validate. Sourcing from the targets makes the described
    # outputs identical to the promised outputs by construction.
    #
    # Phase 6's own targets are excluded: the crate cannot be a part of itself,
    # and RO-Crate 1.1 keeps its metadata descriptor outside hasPart for the
    # same reason.
    RO_CRATE_DERIVED_PARTS = []
    if RO_CRATE_DERIVED_GROUPS:
        _described = {part["path"] for part in RO_CRATE_SOURCE_PARTS} | {
            part["path"] for part in RO_CRATE_MANIFEST_PARTS
        }
        _derived_targets = sorted(
            {
                target
                for target in targets(*[p for p in ALL_PHASES if p != "fair"])
                if target not in _described
            }
        )
        for _group in RO_CRATE_DERIVED_GROUPS:
            if len(_derived_targets) < _group["expect_min"]:
                raise WorkflowError(
                    f"config ro_crate.parts[{_group['group']}]: the workflow "
                    f"declares {len(_derived_targets)} targets, fewer than the "
                    f"{_group['expect_min']} expected. Either a phase switch is "
                    "off or a TARGETS list has emptied — a crate describing no "
                    "outputs would otherwise validate."
                )
            RO_CRATE_DERIVED_PARTS += [
                {"path": _path, "group": _group["group"]}
                for _path in _derived_targets
            ]

    RO_CRATE_ALL_PARTS = (
        RO_CRATE_SOURCE_PARTS + RO_CRATE_MANIFEST_PARTS + RO_CRATE_DERIVED_PARTS
    )


    rule p6t2_build_ro_crate:
        """RO-Crate 1.1 metadata for the whole workflow (P6-T2, FAIR F2/F3/A2).

        Inputs carry the GEO accession as a PID, outputs carry checksums, the
        licences are distinct for code and derived data, and the ontology terms
        come from p6t2a with their `relation` intact.

        EVERY DESCRIBED FILE IS A DECLARED INPUT. That is what makes the crate a
        node in the DAG rather than a snapshot of whatever was on disk: a
        described file that no rule produces fails the build, and a rule that
        reruns rebuilds the crate.

        RO-Crate 1.1 EXPLICITLY — ro-crate-py 0.15 emits 1.2 by default and the
        acceptance criterion names 1.1. A crate that claims a version nobody
        checked it against is an unfalsifiable claim.

        `datePublished` is the HEAD COMMIT DATE, not the run date: the crate
        describes a commit, and a wall-clock stamp would make this the one
        output in the project that never reproduces byte-for-byte.
        """
        input:
            manifest=config["fair"]["ro_crate_yaml"],
            ontology=f"{METADATA_DIR}/ontology-terms.tsv",
            # CITATION.cff is read for the crate's `version` and is already in
            # the `documentation` group, but naming it here makes the
            # dependency visible rather than incidental.
            citation="CITATION.cff",
            parts=[part["path"] for part in RO_CRATE_ALL_PARTS],
        output:
            # AT THE REPOSITORY ROOT, not under metadata/. RO-Crate 1.1 requires
            # the metadata file to sit at the root of the crate it describes,
            # because every data entity @id is a path relative to it. Under
            # metadata/ the @ids would resolve one directory too deep and the
            # crate would not validate — and "validates against the RO-Crate 1.1
            # spec" is P6-T2's acceptance criterion. PROJECT_PLAN §4.2's tree
            # puts it at metadata/ro-crate-metadata.json; that is a defect in
            # the plan, recorded in ADR 0030 §7. The crate's companions — which
            # are this project's own reports about it, not part of the spec —
            # stay in metadata/.
            crate="ro-crate-metadata.json",
            summary=f"{METADATA_DIR}/ro_crate_summary.json",
        params:
            fair=config["fair"],
            manifest=RO_CRATE,
            parts=RO_CRATE_ALL_PARTS,
        log:
            f"{PATHS['logs']}/p6t2_build_ro_crate.log",
        benchmark:
            f"{PATHS['benchmarks']}/p6t2_build_ro_crate.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/build_ro_crate.py"

    rule p6t2b_validate_ro_crate:
        """Check the crate against RO-Crate 1.1 and against the disk (P6-T2).

        A SEPARATE RULE because a builder that grades its own output grades it
        kindly, and because P6-T2's acceptance criterion is validation, not
        generation.

        Four families of check, each written so a wrong crate FAILS it:
        1.1 structure (descriptor, conformsTo, root data entity, required
        properties); referential integrity (no dangling @id, every hasPart both
        described and present); agreement with the bytes on disk (every File
        entity re-hashed INDEPENDENTLY, not compared against what the builder
        recorded); and the Q2 rule — every compartment annotation pointing at a
        cell-type term must say `enriched_for`, because the crate is a second
        file from the one common.smk guards and only this one is read by
        machines downstream.

        Internal validation only, and ADR 0017's lesson applies: checking
        yourself against yourself measures consistency, not correctness. The
        external check is out of band (ADR 0030) because no validator is in any
        pinned environment, and adding one would re-open ADR 0029's hole.
        """
        input:
            crate="ro-crate-metadata.json",
        output:
            summary=f"{METADATA_DIR}/ro_crate_validation.json",
        log:
            f"{PATHS['logs']}/p6t2b_validate_ro_crate.log",
        benchmark:
            f"{PATHS['benchmarks']}/p6t2b_validate_ro_crate.tsv"
        conda:
            "../envs/py-analysis.yaml"
        threads: 1
        script:
            "../scripts/validate_ro_crate.py"



# --- P6-T3 citation and software metadata ----------------------------------
# CITATION.cff and codemeta.json are hand-maintained project metadata, not
# results, so this rule CHECKS them rather than generating them. Both carried a
# repository URL that 404s for five phases and five gates, because nothing ever
# compared them to `git remote`.
rule p6t3_citation_audit:
    """Audit CITATION.cff and codemeta.json against reality (P6-T3).

    Five families of check, each written so a wrong file fails it:

      1. Both files carry every required field, at the right schema version.
      2. They agree with each other AND with `git remote get-url origin` — the
         check that was missing while both pointed at a URL that does not exist.
      3. Every cited package version matches the version in the pin that
         installs it. ADR 0026's lesson in a different costume: a number quoted
         rather than counted from its source goes stale silently, and a reader
         cannot tell a current version string from a stale one. Snakemake is the
         one exception, checked against environment.yml, because it drives the
         workflow from the bootstrap env and is deliberately in no pin.
      4. The DOI policy. Until P6-T5 mints one, neither file may claim a DOI —
         a persistent identifier that does not resolve is worse than none.
      5. ADR 0029, made durable: no env YAML may carry a `pip:` section and
         every env must have a pin. A decision that lives only in an ADR is a
         decision that gets re-made.
    """
    input:
        cff="CITATION.cff",
        codemeta="codemeta.json",
        environment="environment.yml",
        envs=sorted(str(p) for p in Path("workflow/envs").glob("*.yaml")),
        pins=sorted(str(p) for p in Path("workflow/envs").glob("*.pin.txt")),
    output:
        summary=f"{METADATA_DIR}/citation_audit.json",
    params:
        fair=config["fair"],
    log:
        f"{PATHS['logs']}/p6t3_citation_audit.log",
    benchmark:
        f"{PATHS['benchmarks']}/p6t3_citation_audit.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/citation_audit.py"


# --- P6-T4 the report ------------------------------------------------------
# `snakemake --report` is a COMMAND, not a rule, and that is the one sanctioned
# exception to hard constraint 5 (ADR 0030 §3): it consumes the completed DAG,
# so it cannot be a node inside it. What can be a rule is the condition that
# makes the report worth generating, and this is it.
rule p6t4_report_audit:
    """Check that `snakemake --report` will be worth opening (P6-T4).

    Every figure wrapped in report(); every caption file present and non-empty;
    every category phase-shaped, with the app notebooks under "Interactive";
    no orphaned caption file; and every reported figure asked for by a target,
    because a figure no target requests is a hole in the report.

    An uncaptioned figure still renders. It just renders uselessly, and nobody
    notices until a reader asks what they are looking at.

    EVERY SCANNED FILE IS A DECLARED INPUT. p4t7_language_audit shipped with no
    rerun trigger and would have run once and then asserted a green result
    forever (NEXT_STEPS lesson 6); a check that can go stale is worse than no
    check, because its output keeps making a claim nobody re-tested.
    """
    input:
        rules=sorted(str(p) for p in Path("workflow/rules").glob("*.smk")),
        captions=sorted(str(p) for p in Path("workflow/report").glob("*.rst")),
    output:
        summary=f"{METADATA_DIR}/report_audit.json",
    params:
        paths=PATHS,
        # The DAG's own answer to "what does this workflow promise to produce",
        # so a reported figure that nothing asks for is caught here rather than
        # discovered as a gap in the rendered report.
        #
        # Phase 6 is excluded explicitly rather than by accident of ordering:
        # TARGETS_FAIR is populated at the bottom of this file, so a bare
        # targets(*ALL_PHASES) here would evaluate to the same thing while
        # LOOKING like it covered Phase 6. It produces no figures, but a reader
        # should not have to know the include order to know that.
        targets=targets(*[_phase for _phase in ALL_PHASES if _phase != "fair"]),
    log:
        f"{PATHS['logs']}/p6t4_report_audit.log",
    benchmark:
        f"{PATHS['benchmarks']}/p6t4_report_audit.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    script:
        "../scripts/report_audit.py"


# Phase 6 is under construction: the switch in config.yaml is the user's INTENT,
# this list is what actually EXISTS, and a phase contributes to `rule all` only
# when both agree. Populate the rules and this list in the same commit, as every
# other phase did.
TARGETS_FAIR = []

if ONTOLOGY_TERMS is not None:
    TARGETS_FAIR += [
        f"{METADATA_DIR}/ontology-terms.tsv",
        f"{METADATA_DIR}/ontology_resolution_summary.json",
    ]

if RO_CRATE is not None:
    TARGETS_FAIR += [
        "ro-crate-metadata.json",
        f"{METADATA_DIR}/ro_crate_summary.json",
        f"{METADATA_DIR}/ro_crate_validation.json",
    ]

TARGETS_FAIR += [
    f"{METADATA_DIR}/citation_audit.json",
    f"{METADATA_DIR}/report_audit.json",
]


# A LITERAL path in the ro_crate manifest must either exist now or be something
# the workflow promises to produce. This sits at the BOTTOM of the file rather
# than beside the rule because TARGETS_FAIR is populated above and Phase 6
# produces two of the files the crate describes — checked any earlier, the
# ontology mapping would look like a typo.
#
# common.smk cannot make this check at all: it does not know the targets.
# Without it a typo becomes a missing-input error several minutes into a run.
# `Snakefile` was exactly that typo — the file is `workflow/Snakefile`, and
# nothing caught it until the builder tried to hash it.
if RO_CRATE is not None:
    _promised = set(targets(*ALL_PHASES))
    _unaccounted = sorted(
        _part["path"]
        for _part in RO_CRATE_SOURCE_PARTS
        if not Path(_part["path"]).exists() and _part["path"] not in _promised
    )
    if _unaccounted:
        raise WorkflowError(
            f"config ro_crate.parts names {_unaccounted}, which do not exist "
            "and which no rule produces. A literal path in the manifest is a "
            "promise that the crate will describe that file."
        )
