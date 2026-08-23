---
description: Run the FAIR checklist and report which items fail
allowed-tools: Bash(ls:*), Bash(test:*), Bash(git log:*), Read, Grep
---

Audit this repository against the FAIR self-assessment in
`Markdowns/PROJECT_PLAN.md` §7, principle by principle (F1–F4, A1–A2, I1–I3,
R1.1–R1.3).

For each principle report: **pass / partial / fail**, the evidence you checked
(an actual path or command output), and — for anything not passing — the
smallest concrete change that would fix it.

Things that are cheap to check and easy to forget:

- `CITATION.cff` and `codemeta.json` present, parseable, and *current*.
- Licences present and distinct for code (MIT) vs. derived data (CC-BY-4.0).
- `metadata/ontology-terms.tsv` — every metadata field mapped to
  UBERON / CL / NCIT / EFO / OBI (§7.2). Unmapped fields are an **I** failure.
- `resources/checksums.sha256` covering every raw input.
- Env specs pinned, with locks committed.
- `uns` block in the `.h5ad` carrying git SHA, config hash, run date.
- Every figure carrying a `report(..., category=...)` caption.
- RO-Crate present and valid (Phase 6).

End with a count: N pass, N partial, N fail. Do not grade generously — the
point of this command is to find the gaps before a reviewer does.
