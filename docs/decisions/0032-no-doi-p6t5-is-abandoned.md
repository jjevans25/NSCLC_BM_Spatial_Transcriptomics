# ADR 0032 — P6-T5 is abandoned: no DOI is minted, and the mechanism stays

**Date:** 2026-09-16
**Status:** Accepted
**Task:** P6-T5 (abandoning it)
**Related:** ADR 0031 (Gate 6 record), ADR 0030 (the crate), ADR 0005 (pinned
identifiers)

## Context

P6-T5 asks for a Zenodo release: tag `v1.0.0`, archive through the
GitHub–Zenodo integration, mint a DOI, backfill it into `README.md` and
`CITATION.cff`. Its acceptance criterion is "a resolvable DOI".

Everything up to the mint is done. Gate 6 passed (ADR 0031), the clean room is
green, and `fair.identifiers.doi` sits at `null` waiting for a value.

**The remaining step is not one this project can take.** It needs a Zenodo
account, the integration switched on for this repository in a browser, and a
published GitHub release — and it has an ordering constraint that makes it worse
than it looks: **the integration must be enabled before the release is
published**, or Zenodo never sees it and a second release has to be cut.

**The decision, taken by the repository owner: do not mint one.**

## Decision

### 1. No DOI. P6-T5 is abandoned, not deferred

"Deferred" is what `docs/NEXT_STEPS.md` said an hour ago, and it is the wrong
word for a task nobody intends to do. A deferred task is an open item that
accrues; an abandoned one is closed, and closing it is what lets the next reader
stop wondering.

This is the same move Phase 5's own pre-registration permitted for itself
("a documented abandonment") and the same move ADR 0014 made for Aim A5. **The
project's position is that stopping deliberately, in writing, beats an open item
that never closes.**

### 2. The mechanism stays, fully wired and fully enforced

**Nothing is removed.** `fair.identifiers.doi` stays in `config/config.yaml` at
`null`; `build_ro_crate.py` still writes the crate's `sameAs` the moment the
field is populated; and `p6t3_citation_audit` still enforces the policy **in
both directions** —

- with no DOI in config, it **fails** if `CITATION.cff` or `codemeta.json`
  claims one;
- with a DOI in config, it **fails** if either file is missing it.

So minting one later is a one-line edit plus a re-run, and a half-finished
backfill still cannot ship. **The capability is not the same thing as the
artifact**, and only the artifact is being declined.

### 3. No placeholder, ever

`doi` stays `null` rather than becoming `10.5281/zenodo.XXXXXXX` or a
"pending" string. **A persistent identifier that does not resolve is worse than
an absent one**: it looks like provenance, survives copy-paste into someone
else's bibliography, and fails silently years later. The audit refuses one by
construction, and this ADR is the reason not to relax it.

### 4. What this costs, stated rather than glossed

Two FAIR indicators are now **permanently partial rather than pending**:

| | Was | Now |
|---|---|---|
| **F1** globally unique PIDs | partial — "no DOI until P6-T5" | **partial, closed.** The inputs have PIDs (GEO accession, PMID); the *workflow* has a repository URL and a git SHA, which are globally unique but not persistent in the archival sense |
| **F4** indexed and searchable | partial — "Zenodo record pending" | **partial, closed.** Public GitHub repository, no archival index |

The honest reading: **the data this project analyses is FAIR at the source, and
the reanalysis is reproducible, documented and openly licensed but not
archivally citable.** That is a real gap and it is smaller than it sounds for a
learning project whose own README says it is not a publication.

What is *not* affected: the RO-Crate, the checksums, the pinned environments, the
provenance record, the clean-room reproduction, and every other FAIR indicator —
**12 pass, 3 partial, 0 fail** is unchanged by this decision, because both
affected indicators were already partial.

## Consequences

- **Phase 6 is complete**, and with it the project. `docs/NEXT_STEPS.md` no
  longer carries an open item.
- **`CITATION.cff` cites the workflow by repository URL and version**, which is
  what a citation of it can honestly be.
- **If a DOI is ever wanted**, the route is unchanged and takes about an hour:
  enable the Zenodo integration **first**, tag `v1.0.0`, publish a GitHub
  release, then set `fair.identifiers.doi`, add `doi:` to `CITATION.cff` and
  `identifier` to `codemeta.json`, and re-run. The audit will refuse anything
  incomplete, and the crate will pick the value up on its own.
- **The version number stays `1.0.0`** in both citation files. It marks a
  finished piece of work, which this is; it was never a claim that an archive
  exists.
