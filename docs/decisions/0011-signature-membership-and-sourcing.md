# ADR 0011 — signature membership is the plan's; the source is cited for the concept

**Date:** 2026-08-28
**Status:** Accepted
**Task:** P2-T1
**Applies to:** `config/signatures.yaml`, `workflow/schemas/signatures.schema.yaml`

## Context

PROJECT_PLAN §6's P2-T1 gives both a gene list and an accept criterion, and the
two pull in different directions:

> cytotoxicity (`GZMB, PRF1, GNLY, NKG7`), exhaustion (…), antigen presentation
> (HLA-I/II, `B2M, TAP1, TAP2, NLRC5`), M1/M2 myeloid modules, TLS (…).
> **Each gene set carries a `source:` field** (MSigDB ID or PMID). Un-sourced
> gene sets are how reanalyses become unreproducible.

The plan's memberships are **not** any published set. Rooney's cytolytic
activity metric is two genes (GZMA, PRF1), not four. Cabrita's TLS signature is
eight genes and contains neither CR2 nor MS4A1. No single publication defines
the six-gene exhaustion list, and TOX post-dates the review that covers the
other five. So "carry a source" cannot mean "the source defines this list"
without changing the lists.

Three options were considered and the choice was put to the project owner:

1. **Plan-literal, enumerated** — keep the plan's membership, enumerate what it
   abbreviates, cite the source for the *concept*.
2. **Published sets verbatim** — replace membership with canonical MSigDB /
   published sets, so membership and source are the same object.
3. **Hybrid** — extend each set to a canonical superset where one exists.

## Decision

**Option 1.** Membership is PROJECT_PLAN §6's, enumerated where the plan
abbreviates. `source:` cites the publication or MSigDB set for the construct,
and a `source_note:` states, per set, how the membership differs from the cited
source's own gene list.

Two reasons. Sets here are scored on 23 AOIs against a background-relative
detection rule at `TIME-B` n = 8, and the larger published sets (REACTOME MHC-I
+ MHC-II is ~40 genes against the plan's 13) dilute a small signal with genes
that are mostly below background in brain — where the background is *higher*
(ADR 0008). And option 3 would have produced two different provenance stories
inside one file, since only some sets have a canonical superset containing the
plan's genes; `tls` does not.

`HLA-I/II` is enumerated as HLA-A/B/C + B2M + TAP1/TAP2 (class I) and
HLA-DRA/DRB1/DPA1/DPB1/DMA + CIITA (class II). NLRC5 and CIITA are
transactivators — regulators of the pathway rather than members of it, so they
appear in neither Reactome set; they are in the plan's list and are kept.

### The citations were wrong and this is why the task exists

Four PMIDs were proposed for this file at the start of Phase 2. Queried against
PubMed, **two pointed at unrelated papers**:

| proposed | actually is | corrected to |
|---|---|---|
| PMID 30635236 | Kurtulus et al. 2019, *Immunity* — PD-1⁻CD8⁺ TILs | **29634943** (Thommen & Schumacher 2018, *Cancer Cell*) |
| PMID 17082599 | He et al. 2006, *J Immunol* — CD8⁺ IL-17 T cells | **17082649** (Martinez et al. 2006, *J Immunol*) |

Both were plausible — right journal, right year, adjacent subject — and neither
would have been caught by reading the file. A schema can enforce that a source
is *present* and well-formed; it cannot enforce that it is *correct*. Citations
are verified against the record, and the schema's `PMID:<digits>|MSigDB:M<digits>`
pattern exists only to stop free text creeping in.

Verified and kept: PMID 25594174 (Rooney 2015, *Cell*), PMID 31942071 (Cabrita
2020, *Nature*), PMID 31207603 (Khan 2019, *Nature*, cited for TOX),
MSigDB M1062 (Reactome R-HSA-983170) and M705 (Reactome R-HSA-2132295).

### Membership is resolved by a rule, not trusted

`p2t1_resolve_signatures` re-resolves every symbol against the 18,694 genes in
the `.h5ad` on each run and **fails** on one that does not appear. All 47 gene
slots (44 unique) resolve today. A signature that silently shrank from four
genes to three is not a smaller signature, it is a different one, and the
scoring step must never be the place that discovers it.

The rule deliberately does **not** filter by detection. Whether a gene clears
background is per-AOI and per-site, and it is P2-T2's job; conflating "absent
from the panel" with "below background in brain" is the exact error ADR 0008
warns about.

## Consequences

- **Adding or removing a gene from any set here remains a stop-and-ask**
  (`CLAUDE.md`), including a rename resolved against HGNC.
- **`antigen_presentation` is the set Gate 2 turns on** — the published
  direction is *reduced* antigen presentation in brain. Its membership was
  settled before any score was computed.
- **`exhaustion` is expected to fail its coverage floor in brain.** ADR 0008
  measured CTLA4 above background in 2 of 8 `TIME-B` AOIs and TIGIT in 1. When
  it does, it is reported as **"not assessable in brain"**, never as "lower in
  brain".
- **No gene is shared between two sets, and that is now measured rather than
  assumed.** `p2t1_resolve_signatures` computes pairwise overlap and records it
  in `signature_membership_summary.json`; it reports overlap rather than
  rejecting it, since overlap would be legitimate. Today `overlaps` is empty —
  47 gene slots, 47 unique genes — so the six scores are not coupled through
  shared membership. This will NOT hold across files: `exhaustion` carries
  CTLA4 and TIGIT, which the Phase 3 checkpoint panel will also carry, and the
  overlap check here does not look outside `config/signatures.yaml`. Any Phase 3
  statement relating a checkpoint to the exhaustion score must account for that.
- `config/signatures.yaml` is a separate file from `config/config.yaml` on
  purpose: the latter is a declared input of `p0t7_assemble_h5ad`, so editing it
  rebuilds the `.h5ad` and re-runs Phase 1. A gene-set edit must not cost that.
