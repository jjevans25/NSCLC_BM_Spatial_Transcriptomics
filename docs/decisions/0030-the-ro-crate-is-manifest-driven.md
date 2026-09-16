# ADR 0030 — The RO-Crate is manifest-driven and rule-generated; the report is a command

**Date:** 2026-09-16
**Status:** Accepted
**Task:** P6-T2, P6-T4
**Related:** ADR 0005 (pinned digests and read-only resources), ADR 0017
(external validation), ADR 0029 (a pin must carry every dependency)

## Context

PROJECT_PLAN §6 asks Phase 6 for `metadata/ro-crate-metadata.json` describing
"inputs (GEO accession as PID), workflow, outputs, licences, authorship,
ontology terms", accepted when it "validates against the RO-Crate 1.1 spec". It
also asks for a `snakemake --report`, accepted when the report "opens standalone
with no broken links".

Both raise the same question and answer it differently, and both answers are
worth recording because both are easy to get wrong in a way that still looks
finished.

## Decision

### 1. The crate is built from a manifest, not from a filesystem walk

`config/ro_crate.yaml` names what the crate describes. Every entity resolves to
a concrete path at parse time and those paths are **declared inputs** of
`p6t2_build_ro_crate`.

The obvious alternative — walk `results/` and `workflow/` and describe what is
there — fails three ways at once. It cannot declare its inputs, so the DAG does
not know the crate depends on them and Snakemake cannot order it last. It
describes whatever happened to be on disk, so a stale figure from a deleted rule
enters the crate as though it were current. And it is **silently empty on a
clean checkout**: a glob over `results/figures/` matches nothing before anything
has run, and the crate would describe zero figures and still validate — in P6-T1,
the one run whose entire purpose is to be a clean checkout.

Two kinds of group exist for that last reason:

- **`include` (source).** Committed files, as globs expanded against the working
  tree. Deterministic because git decides what is there, with an `expect_min`
  per group because **a glob that stops matching returns an empty list rather
  than an error.**
- **`derived: true`.** Everything the workflow declares it produces, taken from
  the phase `TARGETS_*` lists rather than from a glob. This makes the set of
  described outputs **identical to the set of promised outputs by
  construction** — the crate cannot describe an output the workflow does not
  produce, and cannot omit one it does.

### 2. RO-Crate 1.1 explicitly, because that is what the acceptance criterion names

`ro-crate-py` 0.15.1 emits **1.2** by default. The crate is built with
`ROCrate(version="1.1")`, and `conformsTo` says `https://w3id.org/ro/crate/1.1`.

A crate that claims one version while having been checked against another is the
same defect as an uncounted multiplicity denominator (ADR 0026): the claim is
unfalsifiable by the reader, and the file carries no evidence of which standard
anyone actually applied.

### 3. `datePublished` is the HEAD commit date, not the run date

RO-Crate 1.1 requires the property. `datetime.now()` would make the crate the
**one output in this project that never reproduces byte-for-byte**, which would
break P6-T1's comparison of a clean-room run against this working tree — and
that comparison is the acceptance criterion for the whole phase.

The commit date is also the more honest answer to "when was this published":
**the crate describes a commit.** Rebuilding the same commit tomorrow describes
the same thing and should say so.

### 4. Validation is a separate rule, and the external check is out of band

`p6t2b_validate_ro_crate` reads only the written file. A builder that grades its
own output grades it kindly, and P6-T2's acceptance criterion is validation, not
generation.

It asserts four families, each written so a wrong crate fails it: 1.1 structure;
referential integrity (no dangling `@id`, every `hasPart` both described and
present on disk); **agreement with the bytes on disk**, by re-hashing every File
entity *independently* rather than comparing against what the builder recorded;
and the Q2 rule below.

**The external check stays out of band**, run by hand with
`uvx --from roc-validator rocrate-validator -y validate --profile-identifier
ro-crate-1.1 .` and recorded here. It is not a rule because no validator is in
any pinned environment, and adding one through a `pip:` section is exactly the
hole ADR 0029 just closed.

**It earned its place on the first run.** The internal validator passed the
crate **23 of 23**; the external one failed it on **two REQUIRED checks**, and
both were real:

1. **`sha256` is not in the RO-Crate 1.1 context.** A compacted JSON-LD document
   may only use keys its `@context` defines, and every one of the 284 File
   entities carried a checksum — so the crate violated REQUIRED 2.1 while
   ro-crate-py wrote it happily. Fixed by declaring the term against
   `http://schema.org/sha256`, the URI **RO-Crate 1.2 later standardised**, so
   the crate says what 1.2 says rather than inventing a private term.
2. **The GEO dataset was not reachable from the root through `hasPart`.** A
   `Dataset` at an absolute URI is a **web-based data entity**, and REQUIRED
   14.1 says every data entity must be linked from the root. This project's
   validator had asserted the *opposite* — an explicit check named
   `external_datasets_are_not_parts`, on the reasoning that a crate cannot
   contain GEO. The reasoning was wrong: `hasPart` here means "this crate is
   about this", not "this crate ships this". ro-crate-py had warned about it on
   read and the warning was dismissed as a heuristic.

That is ADR 0017's lesson, demonstrated rather than asserted: **checking
yourself against yourself measures consistency, not correctness** — and the
second defect is the sharper case, because the internal check was not merely
silent, it confidently encoded the wrong rule. Both checks are now inverted or
added in `validate_ro_crate.py`, so the internal validator would catch either
recurrence. After the fixes: **24/24 internal, 38/38 external REQUIRED.**

### 5. The Q2 caveat is carried into RDF as `enriched_for`, never `is_a`

A compartment label is never a cell-type label. PanCK was the only collection
mask; CD45 and GFAP guided where a pathologist placed an ROI, and nothing was
collected on either.

`p4t7_language_audit` can catch the forbidden phrase in a document. **It cannot
stop a JSON-LD graph from asserting the claim in RDF**, and a graph is consumed
by machines that cannot read it sceptically — the machine-readable overclaim
would outlive every document that qualifies it.

So `config/ontology_terms.yaml` carries an explicit `relation`, a compartment
mapped to a CL term may only be `enriched_for`, **`common.smk` refuses the
manifest at parse time if one says `is_a`**, and `p6t2b_validate_ro_crate`
refuses the crate if one does. Two files, two checks, deliberately: the manifest
and the crate are different artifacts and only one of them is read downstream.

### 6. `snakemake --report` is a command, not a rule — the one sanctioned exception to hard constraint 5

Hard constraint 5 says every output is produced by a Snakemake rule. The report
**cannot** be, and this is a structural fact rather than a convenience: it
consumes the completed DAG. A rule producing it would be a rule whose input is
the result of running every rule, including itself.

Stating the exception is the point. An unstated exception is indistinguishable
from a violation, and the next reader would either "fix" it into a cycle or
conclude the constraint is negotiable.

What **can** be a rule is the condition that makes the report worth generating,
and `p6t4_report_audit` is it: every figure wrapped in `report()`, every caption
file present and non-empty, every category phase-shaped with the app notebooks
under "Interactive", no orphaned caption, and every reported figure asked for by
a target. An uncaptioned figure still renders — it just renders uselessly, and
nobody notices until a reader asks what they are looking at.

Every scanned file is a declared input, so the audit cannot report a stale pass.
That is `p4t7_language_audit`'s lesson (NEXT_STEPS lesson 6), applied on the
first attempt this time rather than after review.

### 7. The crate sits at the repository root, and PROJECT_PLAN's tree is wrong about that

PROJECT_PLAN §4.2 shows `metadata/ro-crate-metadata.json`. **That location cannot
validate.** RO-Crate 1.1 requires the metadata file at the root of the crate it
describes, because every data entity `@id` is a path resolved relative to it —
under `metadata/` every path would resolve one directory too deep, and the root
data entity `./` would denote `metadata/` rather than the repository.

So the crate is written to `ro-crate-metadata.json` at the repository root. Its
companions — `metadata/ro_crate_summary.json` and
`metadata/ro_crate_validation.json` — stay under `metadata/`, because those are
this project's reports *about* the crate and the spec says nothing about them.

This is a defect in the plan rather than a decision against it: §6 P6-T2's own
acceptance criterion is "validates against the RO-Crate 1.1 spec", and §4.2's
tree contradicts it. Recorded here the way ADR 0022 §5 recorded two defects in
§6's Phase 4 text, because `Markdowns/PROJECT_PLAN.md` is gitignored and the
ADRs are the durable record.

## Consequences

- **The crate is reproducible.** Same commit, same inputs, byte-identical
  `ro-crate-metadata.json` — so it participates in P6-T1's comparison like every
  other output instead of being excluded from it.
- **A described file that no rule produces fails the build**, and a rule that
  reruns rebuilds the crate. The crate cannot drift from the repository.
- **`expect_min` is now the thing standing between a renamed directory and a
  crate that describes nothing.** It is a floor, not a count: it catches
  collapse, not omission of a single file.
- **The report's acceptance criterion is split.** The automatable half is a rule
  and fails the build; "opens standalone with no broken links" is checked by
  opening it, and P6-T1's record says who opened it and when.
- **Phase 6's own outputs are not described by the crate.** It cannot be a part
  of itself, and RO-Crate 1.1 keeps its metadata descriptor outside the root
  dataset's `hasPart` for the same reason.
