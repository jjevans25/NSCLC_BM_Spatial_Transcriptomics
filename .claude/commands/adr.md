---
description: Write an architecture decision record into docs/decisions/
argument-hint: <short title of the decision>
allowed-tools: Bash(ls:*), Read, Write
---

Write an ADR for: **$ARGUMENTS**

1. `ls docs/decisions/` to find the next number (zero-padded, four digits).
2. Create `docs/decisions/NNNN-kebab-case-slug.md` with exactly these sections:

```markdown
# ADR NNNN — <Title>

**Date:** <today, absolute>
**Status:** Accepted
**Task:** <P0-T5 etc., if there is one>

## Context

What situation forced a decision. Include the numbers that mattered — an ADR
that says "we chose a threshold" without saying which alternatives were on the
table and what they cost is not a record, it is a note.

## Decision

What was chosen, stated so someone can act on it without reading the context.

## Consequences

What this makes easy, what it makes hard, and what it rules out. Include the
consequences you dislike — those are the ones future-you needs.
```

Rules:

- **An ADR records a decision that was actually made**, with its real reasoning.
  Do not invent a rationale to make a choice look more considered than it was.
- If the decision came out of a review notebook (a QC threshold, say), **cite
  the notebook by path** — that link is the acceptance criterion in P0-T5.
- If it deviates from `Markdowns/PROJECT_PLAN.md`, say so explicitly and say
  which section, so the deviation is visible rather than silent.
- Do not edit an existing ADR to change a decision. Supersede it: new ADR, and
  mark the old one `**Status:** Superseded by ADR NNNN`.
