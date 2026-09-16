"""P6-T2 — build the RO-Crate 1.1 metadata describing this workflow.

Owner task: P6-T2. Driven by rule p6t2_build_ro_crate.

PROJECT_PLAN §6: "`metadata/ro-crate-metadata.json` describing inputs (GEO
accession as PID), workflow, outputs, licences, authorship, ontology terms.
**Accept:** validates against the RO-Crate 1.1 spec."

**The crate is written to the repository ROOT, not to `metadata/`, and the plan
is wrong about that** (ADR 0030 §7). RO-Crate 1.1 requires the metadata file at
the root of the crate it describes, because every data entity `@id` is resolved
relative to it: under `metadata/` every path would resolve one directory too
deep and the root data entity `./` would denote `metadata/` rather than the
repository. The plan's own acceptance criterion is the spec, and the spec and
the plan's directory tree disagree.

**RO-Crate 1.1, explicitly, not the library default.** `ro-crate-py` 0.15.1
emits 1.2 unless told otherwise; the acceptance criterion names 1.1, so the
crate is built with `ROCrate(version="1.1")` and its `conformsTo` says so. A
crate that claims a different version from the one it was checked against is
the same defect as an uncounted multiplicity denominator — the claim is
unfalsifiable by the person reading it.

**Only the metadata file is written.** `crate.write()` would COPY every
described file into a crate directory; `crate.metadata.write()` writes
`ro-crate-metadata.json` alone, leaving the repository as the crate payload and
every `@id` a relative path into it. That is what makes the crate describe *this
repository* rather than a duplicate of it.

**`datePublished` is the HEAD commit date, not today.** RO-Crate 1.1 requires
the property, and `datetime.now()` would make the output differ on every run —
which would destroy the byte-comparison P6-T1's clean room rests on, and would
be the only non-reproducible file the workflow produces. The commit date is the
honest answer to "when was this published": the crate describes a commit.

**The ontology mapping keeps its `relation`, and that is the point of carrying
it into RDF at all.** A compartment label is never a cell-type label (Q2): PanCK
was the only collection mask, and CD45 and GFAP guided where a pathologist
placed an ROI. So `immune` is emitted as a `PropertyValue` whose
`additionalType` is `enriched_for` and whose `valueReference` points at
CL:0000738 — never as a bare `is_a` triple. Prose can be read sceptically; a
JSON-LD graph is consumed by machines that cannot be, and the machine-readable
overclaim would outlive every document that qualifies it.

**What the crate says about access is the FAIR A1.2 statement**, in words,
because "N/A" is not an answer a reader can check: no controlled-access data,
no authentication anywhere in the workflow, every input retrievable over HTTPS.
"""

import hashlib
import json
import mimetypes
import subprocess
from pathlib import Path

import yaml
from rocrate.model.contextentity import ContextEntity
from rocrate.model.person import Person
from rocrate.rocrate import ROCrate

# Extensions the stdlib does not know, or knows differently from how this
# project uses them. Everything else falls through to mimetypes.
ENCODING_FORMAT = {
    ".smk": "text/x-snakemake",
    ".h5ad": "application/x-hdf5",
    ".rda": "application/x-r-data",
    ".RData": "application/x-r-data",
    ".R": "text/x-r-source",
    ".tsv": "text/tab-separated-values",
    ".pin.txt": "text/plain",
    ".sha256": "text/plain",
    ".cff": "application/x-yaml",
    ".yaml": "application/x-yaml",
    ".yml": "application/x-yaml",
    ".md": "text/markdown",
}

# Which licence covers what (R1.1). Code is MIT, derived data and figures are
# CC-BY-4.0, and both files sit at the repository root. The distinction is not
# decorative: the inputs are someone else's data and the outputs are derived
# from them, so a single repository-wide licence would be wrong in one
# direction or the other.
DATA_GROUPS = {"results", "inputs"}

crate_config = snakemake.params.fair
manifest = snakemake.params.manifest
parts = snakemake.params.parts
identifiers = crate_config["identifiers"]
licences = crate_config["licences"]

# Read from p6t2a's OUTPUT, not from the manifest. The manifest carries the pin;
# the TSV carries what OLS4 actually resolved, including the IRI — and taking it
# from the file makes the crate depend on the resolution having succeeded rather
# than on the pin having been written.
with open(snakemake.input.ontology, encoding="utf-8") as _handle:
    _header = _handle.readline().rstrip("\n").split("\t")
    ONTOLOGY_ROWS = [
        dict(zip(_header, _line.rstrip("\n").split("\t")))
        for _line in _handle
        if _line.strip()
    ]

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)


def digest(path):
    """sha256 + size, streamed — the source data workbook alone is 45 MB."""
    hasher = hashlib.sha256()
    size = 0
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(chunk)
            size += len(chunk)
    return hasher.hexdigest(), size


def encoding_format(path):
    # The Snakefile has no extension, and `application/octet-stream` on the
    # workflow's own entry point is the least useful answer the crate could
    # give about it.
    if Path(path).name in ("Snakefile", "LICENSE", "LICENSE-DATA"):
        return "text/x-snakemake" if Path(path).name == "Snakefile" else "text/plain"
    for suffix, fmt in ENCODING_FORMAT.items():
        if path.endswith(suffix):
            return fmt
    guessed, _ = mimetypes.guess_type(path)
    return guessed or "application/octet-stream"


def entry_description(group):
    """The manifest's own words for a group, reused on its directory entities."""
    for candidate in manifest["parts"]:
        if candidate["group"] == group:
            return candidate["description"]
    return ""


def git(*args):
    return subprocess.run(
        ["git", *args], capture_output=True, text=True, check=True
    ).stdout.strip()


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    emit("P6-T2 — RO-Crate 1.1 metadata")
    emit("=" * 70)
    emit()

    commit = git("rev-parse", "HEAD")
    commit_date = git("show", "-s", "--format=%cs", "HEAD")
    emit(f"  commit         {commit}  (the commit the crate was BUILT FROM;\n                 the crate carries CITATION.cff's version, not a SHA —\n                 a file cannot contain the hash of the commit that\n                 contains it)")
    emit(f"  datePublished  {commit_date}  (the COMMIT date, not today — the")
    emit("                 crate describes a commit, and a wall-clock stamp")
    emit("                 would make this the one output that never")
    emit("                 reproduces byte-for-byte)")
    emit(f"  described      {len(parts)} files in "
         f"{len({part['group'] for part in parts})} groups")
    emit()

    crate = ROCrate(version="1.1")

    # `sha256` IS NOT IN THE RO-CRATE 1.1 CONTEXT, and every File entity here
    # carries one. A compacted JSON-LD document may only use keys its @context
    # defines, so 284 `sha256` keys made the crate fail RO-Crate 1.1 REQUIRED
    # 2.1 — while ro-crate-py wrote them happily and this project's own
    # validator passed it 23/23. The external validator found it (ADR 0030 §4,
    # and ADR 0017's lesson exactly: checking yourself against yourself measures
    # consistency, not correctness).
    #
    # The term is declared against the URI RO-Crate 1.2 later standardised, so
    # this says the same thing 1.2 says rather than inventing a private term.
    crate.metadata.extra_terms = {"sha256": "http://schema.org/sha256"}

    dataset = manifest["dataset"]

    crate.name = " ".join(dataset["name"].split())
    crate.description = dataset["description"].strip()
    crate.datePublished = commit_date
    crate.keywords = list(dataset["keywords"])

    root = crate.root_dataset
    root["identifier"] = identifiers["repository"]
    root["url"] = identifiers["repository"]

    # THE VERSION IS THE RELEASE VERSION, NOT THE COMMIT SHA, and the difference
    # is not cosmetic. A file cannot contain the hash of the commit that
    # contains it: writing the SHA here would mean the crate always named its
    # own PARENT commit, and adding a rerun trigger to fix that would put the
    # repository on a treadmill — commit, crate goes stale, rebuild, commit.
    #
    # So the crate carries the version CITATION.cff carries, READ FROM IT rather
    # than repeated, which also makes the two agree by construction.
    # p6t3_citation_audit already holds CITATION.cff and codemeta.json to the
    # same number. The commit the crate was BUILT FROM is recorded in the log
    # and in ro_crate_summary.json, where it is a build record rather than a
    # claim about identity.
    citation = yaml.safe_load(Path("CITATION.cff").read_text(encoding="utf-8"))
    root["version"] = str(citation["version"])

    # FAIR A1.2, in words rather than "N/A". A reader can check this sentence
    # against the workflow; they cannot check an abbreviation.
    root["conditionsOfAccess"] = (
        "Open. No controlled-access data and no authentication anywhere in the "
        "workflow: every input is retrieved over HTTPS from a public archive "
        "(NCBI GEO, the source publication's supplementary files, the "
        "NanoString CellProfileLibrary, CellChatDB and UCSC Xena), each pinned "
        "by SHA-256 and re-verified on every run."
    )

    # --- licences (R1.1) ---------------------------------------------------
    code_licence = crate.add(
        ContextEntity(
            crate,
            licences["code"],
            {
                "@type": "CreativeWork",
                "name": "MIT License",
                "description": "Covers the workflow code (LICENSE).",
            },
        )
    )
    data_licence = crate.add(
        ContextEntity(
            crate,
            licences["data"],
            {
                "@type": "CreativeWork",
                "name": "Creative Commons Attribution 4.0 International",
                "description": (
                    "Covers derived data and figures (LICENSE-DATA). The inputs "
                    "are third-party data under their own terms; only what this "
                    "workflow derives is licensed here."
                ),
            },
        )
    )
    root["license"] = code_licence

    # --- authorship --------------------------------------------------------
    authors = []
    for index, author in enumerate(manifest["authors"]):
        properties = {"name": author["name"], "email": author["email"]}
        # An ORCID is to an author what a DOI is to the work. Absent rather
        # than invented; the crate simply does not assert one.
        author_id = author.get("orcid") or f"#author-{index + 1}"
        authors.append(crate.add(Person(crate, author_id, properties)))
    root["author"] = authors

    # --- the input dataset, by PID (F1) ------------------------------------
    geo = crate.add(
        ContextEntity(
            crate,
            identifiers["geo_url"],
            {
                "@type": "Dataset",
                "name": (
                    f"GEO {identifiers['geo_accession']} — NanoString GeoMx DSP "
                    "whole transcriptome profiles, 120 AOIs"
                ),
                "identifier": identifiers["geo_accession"],
                "description": (
                    "The primary data this workflow reanalyses. Deposited by "
                    "the authors of the source publication; retrieved and "
                    "verified against pinned SHA-256 digests by P0-T2. The "
                    "matrix arrives Q3-normalised, not raw (Q1), so detection "
                    "is necessarily background-relative."
                ),
            },
        )
    )
    root["isBasedOn"] = geo

    # AND listed in hasPart, which is not optional. RO-Crate 1.1 REQUIRED 14.1:
    # a data entity MUST be linked to the root data entity through hasPart, and
    # a `Dataset` at an absolute URI is a WEB-BASED DATA ENTITY, not a
    # contextual one. This project's validator originally asserted the
    # opposite — that an external dataset is not a part — on the reasoning that
    # a crate cannot contain GEO. That reasoning is wrong: hasPart here means
    # "this crate is about this", not "this crate ships this", and the external
    # validator rejected the crate over it.
    root.append_to("hasPart", geo)

    article = crate.add(
        ContextEntity(
            crate,
            identifiers["pubmed_url"],
            {
                "@type": "ScholarlyArticle",
                "identifier": f"PMID:{identifiers['source_pmid']}",
                "name": (
                    "Source publication for GEO "
                    f"{identifiers['geo_accession']} (PMID "
                    f"{identifiers['source_pmid']})"
                ),
                "description": (
                    "The study that generated and published these data. This "
                    "workflow is an independent reanalysis, not a reproduction "
                    "of its pipeline; P3-T2b checks this project's inputs "
                    "against the article's own deposited Source Data and finds "
                    "all 2,243,280 normalised values identical."
                ),
            },
        )
    )
    root["citation"] = article

    if identifiers.get("doi"):
        root["sameAs"] = crate.add(
            ContextEntity(
                crate,
                identifiers["doi"],
                {"@type": "CreativeWork", "name": "Archived release"},
            )
        )
        emit(f"  doi            {identifiers['doi']}")
    else:
        # Written only when populated. There is no DOI and that is a decision,
        # not an omission (ADR 0032): P6-T5 was abandoned rather than deferred.
        # The branch stays because the mechanism stays — setting the config
        # field is all it takes — but nothing here guesses a value, because a
        # DOI that does not resolve is worse than an absent one.
        emit("  doi            none — P6-T5 abandoned by decision (ADR 0032)")

    # --- the ontology mapping (I2), with its relations intact --------------
    annotations = []
    for row in ONTOLOGY_ROWS:
        term = crate.add(
            ContextEntity(
                crate,
                row["iri"],
                {
                    "@type": "DefinedTerm",
                    "name": row["label"],
                    "termCode": row["curie"],
                    "inDefinedTermSet": row["ontology"].upper(),
                },
            )
        )
        annotation = crate.add(
            ContextEntity(
                crate,
                f"#annotation-{row['field']}-{row['value']}",
                {
                    "@type": "PropertyValue",
                    "propertyID": row["field"],
                    "name": row["value"],
                    "value": row["label"],
                    "valueReference": {"@id": term.id},
                    # THE Q2 CAVEAT, MADE MACHINE-READABLE. `enriched_for` says
                    # the compartment is enriched for this cell type by ROI
                    # placement — not sorted, not pure, not a composition
                    # claim. Emitting `is_a` here would put the overclaim into
                    # RDF, where it outlives every document that qualifies it.
                    "additionalType": row["relation"],
                    "description": row["notes"],
                },
            )
        )
        annotations.append(annotation)
    root["variableMeasured"] = annotations
    emit(f"  ontology       {len(annotations)} terms, "
         f"{sum(1 for r in ONTOLOGY_ROWS if r['relation'] == 'enriched_for')}"
         " of them enriched_for (never is_a — Q2)")
    emit()

    # --- the workflow itself -----------------------------------------------
    snakemake_language = crate.add(
        ContextEntity(
            crate,
            "#snakemake",
            {
                "@type": "ComputerLanguage",
                "name": "Snakemake",
                "url": "https://snakemake.readthedocs.io/",
                "version": "8.30.0",
            },
        )
    )

    # --- the described files -----------------------------------------------
    group_entities = {}
    for group in manifest["parts"]:
        group_entities[group["group"]] = crate.add(
            ContextEntity(
                crate,
                f"#part-{group['group']}",
                {
                    "@type": "Collection",
                    "name": group["group"],
                    "description": " ".join(group["description"].split()),
                },
            )
        )

    n_bytes = 0
    n_directories = 0
    for part in parts:
        path = part["path"]

        # TWO TARGETS ARE DIRECTORIES — the WASM exports of the app notebooks,
        # declared with Snakemake's directory(). A directory has no digest, so
        # it becomes a Dataset entity and its ENTRY POINT is described as a File
        # with one. That keeps the integrity claim to something this crate's
        # validator can actually re-check: an aggregate digest over a tree would
        # be a convention nobody verifies, which is the kind of claim this
        # project spends its ADRs refusing.
        if Path(path).is_dir():
            n_directories += 1
            contents = sorted(f for f in Path(path).rglob("*") if f.is_file())
            index = Path(path) / "index.html"
            dataset_properties = {
                "name": Path(path).name,
                "description": (
                    f"{' '.join(entry_description(part['group']).split())} "
                    f"{len(contents)} files; the entry point is index.html and "
                    "is described separately with its digest."
                ),
                "isPartOf": {"@id": f"#part-{part['group']}"},
                "license": {"@id": data_licence.id},
            }
            if index.exists():
                index_sha, index_size = digest(str(index))
                n_bytes += index_size
                crate.add_file(
                    str(index),
                    str(index),
                    properties={
                        "name": "index.html",
                        "sha256": index_sha,
                        "contentSize": str(index_size),
                        "encodingFormat": "text/html",
                        "isPartOf": {"@id": f"#part-{part['group']}"},
                        "license": {"@id": data_licence.id},
                    },
                )
                dataset_properties["hasPart"] = {"@id": str(index)}
            crate.add_dataset(path, path, properties=dataset_properties)
            continue

        sha256, size = digest(path)
        n_bytes += size

        properties = {
            "name": Path(path).name,
            "sha256": sha256,
            "contentSize": str(size),
            "encodingFormat": encoding_format(path),
            "isPartOf": {"@id": f"#part-{part['group']}"},
        }
        if part["group"] in DATA_GROUPS:
            properties["license"] = {"@id": data_licence.id}
        else:
            properties["license"] = {"@id": code_licence.id}

        if path == "workflow/Snakefile":
            # The Workflow RO-Crate profile's central entity. The profile is
            # NOT claimed in conformsTo — this crate is checked against
            # RO-Crate 1.1 and nothing else, and claiming a profile nobody
            # validated would be the same unfalsifiable assertion the
            # version-pinning above avoids.
            properties["@type"] = ["File", "SoftwareSourceCode", "ComputationalWorkflow"]
            properties["programmingLanguage"] = {"@id": snakemake_language.id}
            properties["description"] = (
                "Entry point. Declares the target rules and includes one .smk "
                "per phase; each phase contributes to `all` only when its "
                "config switch is on and its TARGETS list is populated."
            )

        entity = crate.add_file(path, path, properties=properties)
        if path == "workflow/Snakefile":
            root["mainEntity"] = entity

    emit(f"  payload        {n_bytes / 1e6:.1f} MB across {len(parts)} entries "
         f"({n_directories} of them directories, described as Datasets)")
    emit()

    output = Path(snakemake.output.crate)
    output.parent.mkdir(parents=True, exist_ok=True)
    crate.metadata.write(str(output.parent))
    emit(f"wrote {output}")

    # The summary is what /fairscan and a reader check without parsing JSON-LD.
    graph = json.loads(output.read_text(encoding="utf-8"))["@graph"]
    by_group = {}
    for part in parts:
        by_group[part["group"]] = by_group.get(part["group"], 0) + 1

    summary = {
        "task": "P6-T2",
        "conforms_to": "https://w3id.org/ro/crate/1.1",
        "built_from_commit": commit,
        "date_published": commit_date,
        "n_graph_entities": len(graph),
        "n_parts": len(parts),
        "n_directories": n_directories,
        "n_bytes": n_bytes,
        "n_files_by_group": dict(sorted(by_group.items())),
        "n_ontology_terms": len(annotations),
        "doi": identifiers.get("doi"),
        "licence_code": licences["code"],
        "licence_data": licences["data"],
        "note": (
            "datePublished is the HEAD commit date, not the run date: the crate "
            "describes a commit, and a wall-clock stamp would make this the one "
            "workflow output that never reproduces byte-for-byte."
        ),
    }
    Path(snakemake.output.summary).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    emit(f"wrote {snakemake.output.summary}")
    emit()
    emit(f"{len(graph)} entities in the graph. Validation is p6t2b's job, not")
    emit("this rule's — a builder that grades its own output grades it kindly.")
