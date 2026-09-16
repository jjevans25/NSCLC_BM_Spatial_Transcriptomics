"""P6-T2b — check the RO-Crate against RO-Crate 1.1, and against the disk.

Owner task: P6-T2. Driven by rule p6t2b_validate_ro_crate.

PROJECT_PLAN §6 P6-T2's acceptance criterion is "validates against the RO-Crate
1.1 spec", and a builder that grades its own output grades it kindly. So the
check is a separate rule, reading only the written file, and it fails the build.

**Every assertion here is one a wrong crate would fail.** That is the standard
Phase 4 set when its FDR step-down shipped correct-looking and wrong under ties
(NEXT_STEPS lesson 2): write the assertion that would fail if you were wrong,
not the one that confirms you were right.

The checks are in four families:

  **1. RO-Crate 1.1 structure.** The metadata descriptor exists, is typed
  `CreativeWork`, conforms to `https://w3id.org/ro/crate/1.1` and is `about` the
  root; the root data entity is `./`, is a `Dataset`, and carries `name`,
  `description`, `datePublished` and `license`. Every `@id` is unique.

  **2. Referential integrity.** Every `{"@id": ...}` reference in the graph
  resolves — to another entity, to an absolute URI, or to a file that exists.
  A dangling `@id` is the characteristic failure of a hand-assembled JSON-LD
  graph and it is invisible to anything that only parses the file.

  **3. The crate describes the bytes actually on disk.** Every `File` entity's
  `sha256` and `contentSize` are recomputed from the file and compared. This is
  the check that makes the crate evidence rather than an assertion — and it is
  deliberately an INDEPENDENT recomputation rather than a comparison against
  what the builder recorded, because those would agree even if both were wrong.

  **4. The Q2 caveat survived into RDF.** Every compartment annotation pointing
  at a cell-type term must carry `additionalType: enriched_for`. PanCK was the
  only collection mask; CD45 and GFAP guided ROI placement (Q2). `common.smk`
  refuses a manifest that says otherwise and this refuses a crate that says
  otherwise, because the manifest and the crate are two different files and only
  one of them is read by machines downstream.

Plus a round-trip: `ro-crate-py` re-reads the written file. Parsing is a weak
check on its own — it is here because a file that this project's own library
cannot read back is not a crate whatever else is true of it.
"""

import hashlib
import json
import re
from pathlib import Path

from rocrate.rocrate import ROCrate

CRATE_PROFILE = "https://w3id.org/ro/crate/1.1"
CRATE_CONTEXT = "https://w3id.org/ro/crate/1.1/context"
ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")

crate_path = Path(snakemake.input.crate)
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

failures = []
checks = []


def check(name, condition, detail=""):
    checks.append({"check": name, "passed": bool(condition), "detail": detail})
    if not condition:
        failures.append(f"{name}: {detail}")
    return bool(condition)


def references(node):
    """Every {'@id': ...} reachable from a node, at any depth."""
    if isinstance(node, dict):
        if set(node) == {"@id"}:
            yield node["@id"]
            return
        for key, value in node.items():
            if key == "@id":
                continue
            yield from references(value)
    elif isinstance(node, list):
        for item in node:
            yield from references(item)


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    emit("P6-T2b — RO-Crate 1.1 conformance and on-disk agreement")
    emit("=" * 70)
    emit()
    emit(f"  crate    {crate_path}")
    emit(f"  profile  {CRATE_PROFILE}")
    emit()

    document = json.loads(crate_path.read_text(encoding="utf-8"))
    graph = document.get("@graph", [])
    by_id = {entity["@id"]: entity for entity in graph if "@id" in entity}

    # --- 1. structure ------------------------------------------------------
    emit("1. RO-Crate 1.1 structure")

    check(
        "filename_is_ro_crate_metadata_json",
        crate_path.name == "ro-crate-metadata.json",
        f"named {crate_path.name}",
    )
    # The context may be the bare 1.1 URL or a LIST whose first element is it,
    # followed by term declarations. Both are legal; the list form is what a
    # crate needs the moment it uses a key 1.1 does not define.
    context = document.get("@context")
    context_list = context if isinstance(context, list) else [context]
    check(
        "context_is_1_1",
        context_list and context_list[0] == CRATE_CONTEXT,
        f"@context is {context!r}, expected {CRATE_CONTEXT!r} first",
    )

    # EVERY KEY MUST BE DEFINED, and `sha256` is the one this crate uses that
    # RO-Crate 1.1 does not define — 1.2 later added it as schema.org/sha256.
    # A compacted JSON-LD document using an undeclared key is not valid
    # JSON-LD, which is RO-Crate 1.1 REQUIRED 2.1, and neither ro-crate-py nor
    # this validator's first version noticed: an external validator did.
    declared = {}
    for item in context_list:
        if isinstance(item, dict):
            declared.update(item)
    uses_sha256 = any("sha256" in entity for entity in graph)
    check(
        "sha256_is_declared_in_the_context",
        not uses_sha256 or "sha256" in declared,
        f"{sum('sha256' in e for e in graph)} entities carry `sha256`, which "
        "RO-Crate 1.1's context does not define. Declare it in an extra term "
        "map rather than dropping the checksums.",
    )
    check(
        "ids_are_unique",
        len(by_id) == len([e for e in graph if "@id" in e]),
        "two entities share an @id",
    )
    check(
        "every_entity_has_id_and_type",
        all("@id" in e and "@type" in e for e in graph),
        f"{sum(1 for e in graph if '@id' not in e or '@type' not in e)} entities lack one",
    )

    descriptor = by_id.get("ro-crate-metadata.json", {})
    check(
        "descriptor_present",
        bool(descriptor),
        "no entity with @id 'ro-crate-metadata.json'",
    )
    check(
        "descriptor_is_creative_work",
        descriptor.get("@type") == "CreativeWork",
        f"@type is {descriptor.get('@type')!r}",
    )
    check(
        "descriptor_conforms_to_1_1",
        descriptor.get("conformsTo", {}).get("@id") == CRATE_PROFILE,
        f"conformsTo is {descriptor.get('conformsTo')!r}",
    )
    check(
        "descriptor_is_about_root",
        descriptor.get("about", {}).get("@id") == "./",
        f"about is {descriptor.get('about')!r}",
    )

    root = by_id.get("./", {})
    check("root_present", bool(root), "no entity with @id './'")
    check(
        "root_is_dataset",
        root.get("@type") == "Dataset",
        f"@type is {root.get('@type')!r}",
    )
    for required in ("name", "description", "datePublished", "license"):
        check(
            f"root_has_{required}",
            bool(root.get(required)),
            f"root data entity has no {required}",
        )
    check(
        "date_published_is_iso_8601",
        bool(ISO_DATE.match(str(root.get("datePublished", "")))),
        f"datePublished is {root.get('datePublished')!r}",
    )
    emit(f"   {sum(1 for c in checks if c['passed'])}/{len(checks)} structural checks passed")
    emit()

    # --- 2. referential integrity -----------------------------------------
    emit("2. Referential integrity")

    dangling = []
    for entity in graph:
        for reference in references(entity):
            if reference in by_id:
                continue
            if reference.startswith(("http://", "https://", "#")):
                # An absolute URI is a legitimate external reference; a '#'
                # id that is not in the graph is not, and is caught below.
                if reference.startswith("#"):
                    dangling.append((entity["@id"], reference))
                continue
            if not Path(reference).exists():
                dangling.append((entity["@id"], reference))

    check(
        "no_dangling_references",
        not dangling,
        f"{len(dangling)} unresolved: {dangling[:5]}",
    )

    # An EXTERNAL dataset — GEO, referenced by absolute URI through isBasedOn —
    # is a contextual entity, not a part of this crate. ro-crate-py warns about
    # it on read ("looks like a data entity but it's not listed in the root
    # dataset's hasPart"), and the warning is a heuristic rather than a defect:
    # a Dataset at an absolute URI is not something this crate contains.
    # Asserted rather than ignored, so the intent is checked instead of assumed.
    parts = [ref["@id"] for ref in root.get("hasPart", [])]
    external = [
        entity["@id"]
        for entity in graph
        if entity.get("@type") == "Dataset"
        and entity["@id"].startswith(("http://", "https://"))
    ]
    check(
        "web_based_data_entities_are_parts",
        all(entity in parts for entity in external),
        f"{[e for e in external if e not in parts]} is a Dataset at an absolute "
        "URI — a WEB-BASED DATA ENTITY — and RO-Crate 1.1 REQUIRED 14.1 says "
        "every data entity MUST be linked from the root through hasPart.",
    )

    missing_parts = [p for p in parts if p not in by_id]
    check(
        "every_haspart_is_described",
        not missing_parts,
        f"{len(missing_parts)} hasPart ids have no entity: {missing_parts[:5]}",
    )
    # Web-based data entities are excluded: GEO is a part of this crate in the
    # RO-Crate sense (the crate is ABOUT it) and has no local file, which is
    # what makes it a different kind of part rather than a missing one.
    absent_parts = [
        p
        for p in parts
        if not p.startswith(("http://", "https://")) and not Path(p).exists()
    ]
    check(
        "every_local_haspart_exists_on_disk",
        not absent_parts,
        f"{len(absent_parts)} described files are absent: {absent_parts[:5]}",
    )
    emit(f"   {len(parts)} parts, {len(by_id)} entities, {len(dangling)} dangling")
    emit()

    # --- 3. the crate describes the bytes on disk -------------------------
    emit("3. Agreement with the files on disk")

    files = [
        entity
        for entity in graph
        if "File" in (
            entity["@type"] if isinstance(entity["@type"], list) else [entity["@type"]]
        )
    ]
    mismatches = []
    for entity in files:
        path = Path(entity["@id"])
        if not path.exists():
            mismatches.append(f"{entity['@id']}: absent")
            continue
        hasher = hashlib.sha256()
        size = 0
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                hasher.update(chunk)
                size += len(chunk)
        if entity.get("sha256") != hasher.hexdigest():
            mismatches.append(f"{entity['@id']}: sha256 differs")
        if entity.get("contentSize") != str(size):
            mismatches.append(
                f"{entity['@id']}: contentSize {entity.get('contentSize')} != {size}"
            )

    check(
        "every_file_entity_matches_disk",
        not mismatches,
        f"{len(mismatches)} mismatched: {mismatches[:5]}",
    )
    check(
        "every_file_entity_is_complete",
        all(
            e.get("sha256") and e.get("contentSize") and e.get("encodingFormat")
            for e in files
        ),
        "a File entity lacks sha256, contentSize or encodingFormat",
    )
    emit(f"   {len(files)} file entities re-hashed, {len(mismatches)} mismatches")
    emit()

    # --- 4. the Q2 caveat survived into RDF -------------------------------
    emit("4. A compartment label is never a cell-type label (Q2)")

    overclaims = []
    compartment_terms = 0
    for entity in graph:
        if entity.get("@type") != "PropertyValue":
            continue
        if entity.get("propertyID") != "compartment":
            continue
        term_id = entity.get("valueReference", {}).get("@id", "")
        term = by_id.get(term_id, {})
        if term.get("inDefinedTermSet") != "CL":
            continue
        compartment_terms += 1
        if entity.get("additionalType") != "enriched_for":
            overclaims.append(
                f"{entity['@id']} -> {term.get('termCode')} as "
                f"'{entity.get('additionalType')}'"
            )

    check(
        "compartment_cell_type_terms_are_enriched_for",
        not overclaims,
        (
            f"{len(overclaims)} compartment(s) assert a cell type as fact: "
            f"{overclaims}. PanCK was the only collection mask; CD45 and GFAP "
            "guided ROI placement only."
        ),
    )
    emit(f"   {compartment_terms} compartment→cell-type annotations, all enriched_for"
         if not overclaims
         else f"   {len(overclaims)} OVERCLAIM(S)")
    emit()

    # --- round-trip --------------------------------------------------------
    emit("5. Round-trip through ro-crate-py")
    try:
        reread = ROCrate(str(crate_path.parent))
        n_data_entities = len(list(reread.data_entities))
        check(
            "round_trip_reads_back",
            n_data_entities == len(parts),
            f"read back {n_data_entities} data entities, crate lists {len(parts)}",
        )
        emit(f"   read back {n_data_entities} data entities")
    except Exception as exc:  # noqa: BLE001 — any failure here is a failure
        check("round_trip_reads_back", False, f"{type(exc).__name__}: {exc}")
        emit(f"   FAILED: {exc}")
    emit()

    passed = sum(1 for c in checks if c["passed"])
    summary = {
        "task": "P6-T2b",
        "crate": str(crate_path),
        "profile": CRATE_PROFILE,
        "n_checks": len(checks),
        "n_passed": passed,
        "n_failed": len(checks) - passed,
        "n_graph_entities": len(by_id),
        "n_parts": len(parts),
        "n_file_entities": len(files),
        "n_compartment_cell_type_annotations": compartment_terms,
        "checks": checks,
        "valid": not failures,
        "note": (
            "Internal validation. The crate is also checked out of band against "
            "an external validator (P6-T2, ADR 0030) — this project's standing "
            "lesson from ADR 0017 is that checking yourself against yourself "
            "measures consistency, not correctness."
        ),
    }
    Path(snakemake.output.summary).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    emit(f"wrote {snakemake.output.summary}")
    emit()

    if failures:
        for failure in failures:
            emit(f"  FAIL  {failure}")
        raise RuntimeError(
            f"{len(failures)} of {len(checks)} RO-Crate checks failed — see "
            f"{snakemake.log[0]}. Fix the crate; do not relax the check."
        )

    emit(f"VALID. {passed}/{len(checks)} checks passed: the crate conforms to")
    emit("RO-Crate 1.1, every reference resolves, every described file matches")
    emit("its bytes on disk, and no compartment is asserted as a cell type.")
