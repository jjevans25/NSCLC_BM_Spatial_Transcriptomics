# ADR 0003 — Skill security scan: results and accepted risk

**Date:** 2026-08-23
**Status:** Accepted
**Task:** P0-T1
**Implements:** PROJECT_PLAN §5.1.2 ("Security step — do it, don't skip it")

## Context

169 agent skills are installed: 163 from `K-Dense-AI/scientific-agent-skills`,
5 marimo skills, and 1 project-local (`geomx-dsp-analysis`). Skills run with
full agent permissions, so §5.1.2 requires a scan. It had not been run against
the K-Dense set at any point before now.

Command (reproducible):

```bash
uvx --from cisco-ai-skill-scanner skill-scanner scan-all .claude/skills \
    --recursive --use-behavioral
```

Scanner: `skill-scanner 2.0.13`. `--use-behavioral` is static AST + taint dataflow; it needs
no API key. The optional LLM, VirusTotal and Cisco AI Defense analyzers were
**not** run — they require external API keys and would send skill contents to a
third party.

## Result

**148 clean, 21 flagged.** The flagged set:

| Severity | Skills |
|---|---|
| CRITICAL | `autoskill`, `citation-management`, `infographics`, `latex-posters`, `literature-review`, `pacsomatic`, `pi-agent`, `research-lookup`, `scientific-schematics`, `scientific-slides`, `xlsx` |
| HIGH | `docx`, `generate-image`, `geomaster`, `histolab`, `marimo-pair`, `modal`, `nextflow`, `paper-lookup`, `statistical-power`, `waypoint-bio` |

Six skills are load-bearing for this project right now. Five are clean (INFO
only): `marimo-notebook`, `jupyter-to-marimo`, `wasm-compatibility`,
`retro-marimo-pair`, `geomx-dsp-analysis`. The sixth, `marimo-pair`, is
flagged HIGH.

## Triage of the findings that touch Phase 0

Each was checked by reading the flagged source, not by trusting the label.

| Skill | Finding | Verdict |
|---|---|---|
| `marimo-pair` | HIGH — "executes bash but Bash tool not in allowed-tools" | **False positive.** Its frontmatter declares `allowed-tools: Bash(bash **/scripts/discover-servers.sh *), Bash(bash **/scripts/execute-code.sh *), Read`. The rule looks for a bare `Bash` token and cannot parse the scoped `Bash(...)` form, so it misreads a *tighter* declaration as a missing one. |
| `marimo-pair` | MEDIUM — "hardcoded password or secret in variable" | **False positive.** The line is `token="\${MARIMO_TOKEN:-}"` — reading from the environment, the opposite of hardcoding. |
| `statistical-power` | HIGH — "infinite loop without clear exit condition" | **False positive.** `while True:` at `simulate_power.py:81` breaks two lines later on `est_hi.power >= target_power or hi >= 1_000_000`. A bounded doubling search. |
| `paper-lookup` | MEDIUM — "environment variable harvesting" | **False positive.** Reads three named config vars (`OPENALEX_EMAIL`, `OPENALEX_API_KEY`, `CROSSREF_MAILTO`), not the environment at large. |
| `paper-lookup` | HIGH — "dangerous data flow in command pipeline" | **Real pattern, accepted.** `curl … | python3 scripts/arxiv_atom.py -` pipes untrusted remote XML into a local parser. The parser is bundled and inspectable; the risk is parser robustness, not code execution. Accepted for read-only literature queries. |

Both `marimo-pair` scripts were additionally read in full (567 lines) before
first use. `discover-servers.sh` touches only the local registry and loopback;
`execute-code.sh` POSTs to a marimo kernel, warns on non-local hosts, and
prefers `MARIMO_TOKEN` over `--token` (which would be visible in `ps`).

## Decision

- The six skills Phase 0 depends on are cleared for use.
- The remaining 15 flagged skills are **not cleared**. Several appear in the
  §5.1.2 install list for later phases — `literature-review` (P4),
  `citation-management` and `scientific-schematics` (P6). Each gets the same
  read-the-source triage **at the phase that first uses it**, not before.
- Skills flagged and never used by this plan (`autoskill`, `infographics`,
  `latex-posters`, `pacsomatic`, `pi-agent`, `modal`, `waypoint-bio`,
  `geomaster`, `generate-image`, `histolab`, `nextflow`, `scientific-slides`,
  `xlsx`, `docx`, `research-lookup`) are candidates for removal. §5.1.2 argues
  for a topical subset precisely to shrink this surface; installing all 163
  went the other way. Pruning is deferred, not resolved.

## Consequences

- A HIGH or CRITICAL label from this tool is a **prompt to read the source**,
  not a verdict. Four of the five Phase-0-relevant findings were false
  positives on inspection; treating the labels as authoritative would have
  removed working, correctly-scoped skills.
- The scan is a point-in-time result. `skills-lock.json` pins content hashes,
  so an unreviewed change to an installed skill is detectable — re-run the scan
  whenever those pins move.
- No LLM-based analyzer was run, so semantic prompt-injection in skill *prose*
  is not covered by this scan.
