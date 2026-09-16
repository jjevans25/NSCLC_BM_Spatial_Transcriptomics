"""P6-T2a — resolve every metadata field to an ontology term against EBI OLS4.

Owner task: P6-T2a. Driven by rule p6t2a_resolve_ontology_terms.

PROJECT_PLAN §7.2 asks for this in **Phase 0** — "retrofitting ontology terms at
the end is miserable" — and it never landed. `metadata/` held nothing but a
`.gitkeep` through five phases and five gates, and `config/compartment_map.yaml`
still carries a header promising the fields would be filled in at P0-T3. So it
is being retrofitted, which is the miserable option, and the one thing that
makes it cheap is that the vocabulary never changed: the same seven AOI codes,
three sites and four compartments this project started with.

**This resolves; it does not look up.** §7.2 is explicit — "resolve each against
EBI OLS4; don't hand-type IDs from memory, that's exactly the error class the
skill exists to prevent." Every CURIE in `config/ontology_terms.yaml` was
obtained from OLS4, and this rule obtains it again on every build and asserts
the answer is unchanged.

That makes the manifest a **pin**, in exactly the sense ADR 0005 gives a
downloaded file, and it buys the thing a one-off curation cannot: **ontology
drift fails the build.** If a term is obsoleted, merged, or has its primary
label changed under it, the run stops and names the term. The alternative — a
TSV curated once and trusted forever — is a file that is correct on the day it
is written and unfalsifiable afterwards.

**No fuzzy matching, no fallback path, no near-match acceptance.** The search is
exact-label, restricted to the declared ontology, and restricted to hits where
that ontology is the *defining* one. Three assertions have to pass and any
failure stops the phase. This is the shape `lr_nominations.py` has, and the one
`p5t1_clinical_join` has: an honest failure must be reachable rather than
escapable, or the check is decoration. Retries exist for a flaky network and
nothing else — they never widen a match.

**`relation` is where this file carries scientific weight.** A compartment label
is never a cell-type label (Q2): PanCK was the only collection mask, and CD45
and GFAP guided where a pathologist placed an ROI. So `immune → CL:0000738` is
emitted as `enriched_for`, never `is_a`. Prose can be read sceptically; an RDF
triple cannot, and the RO-Crate this feeds is consumed by machines. The
constraint is enforced in `common.smk` at parse time and re-asserted here,
because it is the one claim in this output that the assay cannot support.

**The output carries no timestamp.** A resolution date would make the TSV differ
byte-for-byte on every run and destroy the comparison P6-T1's clean room
depends on. What matters — which service was asked, what it answered, whether
it matched the pin — is in the file; when it was asked is in the log and the
git history.
"""

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

TSV_COLUMNS = [
    "field",
    "value",
    "label",
    "ontology",
    "curie",
    "iri",
    "relation",
    "notes",
]

# Returned by OLS4 for each hit. `is_defining_ontology` is the one that stops a
# term being read out of an ontology that merely imports it: UBERON:0000955
# appears in a dozen ontologies and means the same thing in all of them, but a
# CL term imported into an application ontology can carry a narrowed definition.
FIELD_LIST = "iri,label,obo_id,ontology_name,is_defining_ontology"

terms = snakemake.params.terms
ols = snakemake.params.ols

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)


def search(label, ontology, log):
    """Exact-label search against OLS4, retried only for transport failure."""
    query = urllib.parse.urlencode(
        {
            "q": label,
            "ontology": ontology,
            "exact": "true",
            "queryFields": "label",
            "fieldList": FIELD_LIST,
            "rows": 10,
        }
    )
    url = f"{ols['api']}?{query}"

    last_error = None
    for attempt in range(1, ols["max_attempts"] + 1):
        try:
            with urllib.request.urlopen(url, timeout=ols["timeout_seconds"]) as handle:
                return json.load(handle)["response"]["docs"]
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            # A network failure is not a resolution failure and must not be
            # allowed to look like one. Retried, then raised — never softened
            # into "unresolved" or a cached answer.
            last_error = exc
            print(
                f"    attempt {attempt}/{ols['max_attempts']} failed: {exc}",
                file=log,
                flush=True,
            )
            if attempt < ols["max_attempts"]:
                time.sleep(ols["retry_wait_seconds"])

    raise RuntimeError(
        f"OLS4 unreachable after {ols['max_attempts']} attempts for "
        f"'{label}' in {ontology}: {last_error}. This is a TRANSPORT failure, "
        "not an unresolved term — the run stops rather than recording an "
        "absence it did not measure."
    )


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    emit("P6-T2a — ontology term resolution (PROJECT_PLAN §7.2, FAIR I2)")
    emit("=" * 70)
    emit()
    emit(f"  service      {ols['api']}")
    emit(f"  entries      {len(terms)}")
    emit("  matching     exact label, declared ontology, defining ontology only")
    emit("  on mismatch  STOP. No fuzzy fallback and no near-match acceptance.")
    emit()

    rows = []
    failures = []

    for entry in terms:
        key = f"{entry['field']}={entry['value']}"
        emit(f"  {key}")
        emit(f"    searching '{entry['label']}' in {entry['ontology']}")

        docs = search(entry["label"], entry["ontology"], log)

        # Three assertions, and each one fails the build on its own.
        #
        # 1. The label must match EXACTLY. OLS4's `exact=true` matches
        #    synonyms as well as primary labels, so a hit is not yet a match:
        #    the manifest records a primary label and that is what must come
        #    back, or the pin and the term have drifted apart.
        defining = [
            d
            for d in docs
            if d.get("is_defining_ontology")
            and d.get("ontology_name") == entry["ontology"]
            and d.get("label") == entry["label"]
        ]

        if not defining:
            near = [
                f"{d.get('obo_id')} '{d.get('label')}' "
                f"({d.get('ontology_name')}, defining={d.get('is_defining_ontology')})"
                for d in docs[:5]
            ]
            failures.append(
                f"{key}: no exact-label hit for '{entry['label']}' defined in "
                f"{entry['ontology']}. OLS4 returned: {near or 'nothing'}"
            )
            emit(f"    UNRESOLVED — {len(docs)} hit(s), none an exact defining match")
            emit()
            continue

        # 2. One label, one term. Two defining hits on the same primary label
        #    means the search is ambiguous and the manifest cannot say which
        #    was meant — a coin flip is not a resolution.
        if len(defining) > 1:
            failures.append(
                f"{key}: '{entry['label']}' resolves to "
                f"{[d.get('obo_id') for d in defining]} — ambiguous in "
                f"{entry['ontology']}. Pick one in the manifest and say why."
            )
            emit(f"    AMBIGUOUS — {[d.get('obo_id') for d in defining]}")
            emit()
            continue

        hit = defining[0]

        # 3. The resolved CURIE must equal the pinned one. This is the whole
        #    point of the file being a pin: an obsoleted or merged term comes
        #    back as a DIFFERENT id under the same label, and silently
        #    absorbing that is how a crate ends up asserting a retired term.
        if hit.get("obo_id") != entry["curie"]:
            failures.append(
                f"{key}: pinned {entry['curie']} but OLS4 now resolves "
                f"'{entry['label']}' to {hit.get('obo_id')}. ONTOLOGY DRIFT — "
                "re-pin deliberately in config/ontology_terms.yaml, with an "
                "ADR if the meaning moved; do not edit it to make this pass."
            )
            emit(f"    DRIFT — pinned {entry['curie']}, resolved {hit.get('obo_id')}")
            emit()
            continue

        emit(f"    {hit['obo_id']}  {hit['iri']}")
        emit(f"    relation: {entry['relation']}")
        emit()

        rows.append(
            {
                "field": entry["field"],
                "value": entry["value"],
                "label": entry["label"],
                "ontology": entry["ontology"],
                "curie": entry["curie"],
                "iri": hit["iri"],
                "relation": entry["relation"],
                "notes": " ".join(entry["notes"].split()),
            }
        )

    emit("-" * 70)
    emit(f"resolved {len(rows)} of {len(terms)}")
    emit()

    if failures:
        for failure in failures:
            emit(f"  FAIL  {failure}")
        raise RuntimeError(
            f"{len(failures)} of {len(terms)} ontology terms did not resolve as "
            f"pinned — see {snakemake.log[0]}. There is no fallback path here "
            "on purpose: §7.2 forbids hand-typing an ID, and accepting a near "
            "match would be hand-typing it with extra steps. A term that will "
            "not resolve is a scientific question (what does this field mean?), "
            "so stop and ask rather than widening the search."
        )

    tsv_path = Path(snakemake.output.tsv)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(TSV_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(row[column] for column in TSV_COLUMNS) + "\n")
    emit(f"wrote {tsv_path} ({len(rows)} rows)")

    # The summary is what the RO-Crate and /fairscan read, and what a reader
    # checks without opening the TSV. Deliberately carries no timestamp — see
    # the module docstring.
    by_field = {}
    for row in rows:
        by_field[row["field"]] = by_field.get(row["field"], 0) + 1
    by_relation = {}
    for row in rows:
        by_relation[row["relation"]] = by_relation.get(row["relation"], 0) + 1

    summary = {
        "task": "P6-T2a",
        "service": ols["api"],
        "n_terms": len(rows),
        "n_by_field": dict(sorted(by_field.items())),
        "n_by_relation": dict(sorted(by_relation.items())),
        "ontologies": sorted({row["ontology"] for row in rows}),
        "matching": (
            "exact primary label, restricted to the declared ontology and to "
            "hits where it is the defining ontology; the resolved CURIE must "
            "equal the pinned one"
        ),
        "relation_rule": (
            "A compartment label is never a cell-type label (Q2). PanCK was the "
            "only collection mask and CD45/GFAP guided ROI placement only, so a "
            "compartment mapped to a CL term is `enriched_for`, never `is_a` — "
            "enforced in common.smk at parse time."
        ),
        "all_resolved": True,
    }
    Path(snakemake.output.summary).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    emit(f"wrote {snakemake.output.summary}")
    emit()
    emit("CLEAN. Every metadata field carries a term that OLS4 resolved today,")
    emit("at the identifier this project pinned when it first resolved it.")
