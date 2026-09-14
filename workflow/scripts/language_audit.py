"""P4-T7 — audit the repository for claims the assay cannot support.

Owner task: P4-T7. Driven by rule p4t7_language_audit.

PROJECT_PLAN §6 P4-T7: "grep the entire repo for 'colocali' and 'spatial
proximity' and fix every hit. You have no coordinates."

**This is a rule rather than a shell command someone remembers to run**, because
hard constraint 5 says every output is produced by a Snakemake rule, and because
a check that lives in a handoff note is a check that survives exactly as long as
the person who wrote the note. As a rule, a violation **fails the build**.

The forbidden claims are hard constraint 6 and Q2's caveat:

  * **Never a claim of colocalisation, spatial adjacency or spatial
    proximity.** This is GeoMx DSP. An AOI is a region a pathologist placed;
    there are no coordinates and no spot grid. The supported phrase is
    **"inferred crosstalk between adjacent compartments"**, and the difference
    is not stylistic — one claims a measured spatial relationship and the other
    claims an inference from two compartments' expression.

  * **Never a CD45+ or GFAP+ AOI.** PanCK was the ONLY collection mask (Q2).
    CD45 and GFAP guided where a pathologist placed an ROI; nothing was
    collected on either. A `TIME` AOI is the PanCK-negative segment of an ROI
    sited in a CD45-rich region — **not** a CD45-sorted population — so a
    compartment label is never a cell-type label.

**A match is not automatically a violation, and this is the part PROJECT_PLAN
§6's literal instruction gets wrong.** The plan says to grep the entire repo for the
two forbidden phrases and fix every hit. Run literally, that fails on
`CLAUDE.md`'s own hard constraint 6, on `config/compartment_map.yaml`'s
header, and on `docs/limitations.md` §3 — because **the files that state the
rule have to name the thing they forbid.** A check that cannot tell a
prohibition from a claim would force the rule's own statement to be deleted to
make the check pass, which is the exact inversion of what it is for.

So every match is classified. A match inside a prohibition — a line, or the line
before it, carrying `never`, `not`, `do not`, `forbid`, `avoid`, `rather than`,
`instead of`, `cannot` — is counted as a **prohibition** and reported
separately. Everything else is a **violation** and fails the build. Both counts
are in the report, because a repository whose prohibition count drops to zero
has probably lost its rules rather than fixed its language.

**A quotation is a third category, and it is not a loophole.**
`config/ligand_receptor.yaml` quotes the P4-T6 comparator's own sentence, which
says "colocalized" — because that is what the authors wrote, and ADR 0021 §2
quotes rather than paraphrases precisely so the verdict is not written by us.
**Editing a quotation to satisfy our style rule would falsify the source.** So a
match inside a declared quotation field is counted separately and does not fail
the build. The cue is narrow on purpose — the field must be named for quoting
(`source_claim`), or the surrounding lines must say the text is quoted — so it
cannot be reached by writing a claim and calling it a quote.

The classifier is deliberately generous to prohibitions and quotations and
strict about nothing else: a false prohibition costs one unflagged sentence,
while a false violation makes the rule unrunnable and the rule then gets
deleted.

**The patterns are assembled from fragments so this file does not match itself.**
Writing the literal here would make the audit report its own source on every
run, and the obvious fix — excluding this file by name — would create a place
where the forbidden phrase could be written without consequence. The
scratchpad-quoted terms in the report are split the same way.

Scope is `git ls-files`: what is committed is what makes a claim. Generated
output under `results/` is excluded because it is derived from code this audit
already covers, and a phrase there is a symptom whose cause is in a script.
"""

import json
import re
import subprocess
from pathlib import Path

# Assembled from fragments — see the module docstring. This file must not match
# its own patterns, and excluding it by name would create a blind spot.
PATTERNS = [
    (
        "spatial_claim",
        re.compile(
            r"coloc" + r"ali[sz]" + r"|"
            r"spatial(?:ly)?[\s_-]+(?:proximit|adjacen|close|neighbour|neighbor)",
            re.IGNORECASE,
        ),
        "no coordinates — say 'inferred cross" + "talk between adjacent compartments'",
    ),
    (
        "cell_type_label",
        re.compile(r"\b(?:CD45|GFAP)\s*\+\s*AOI", re.IGNORECASE),
        "PanCK was the only collection mask (Q2) — a compartment label is not a "
        "cell-type label",
    ),
]

# A match inside one of these is the rule being STATED, not broken. The files
# that forbid a phrase have to name it. Checked against the matched line and the
# one before it, because this project's prose wraps — compartment_map.yaml puts
# "Do not write" on one line and the forbidden term on the next.
PROHIBITION_CUES = re.compile(
    r"\b(?:never|not|forbid(?:s|den)?|avoid|rather than|instead of|cannot|"
    r"don'?t|do not|no longer|must not|may not|prohibit(?:s|ed)?)\b",
    re.IGNORECASE,
)

# A match inside a QUOTATION is someone else's claim, reported verbatim. The
# P4-T6 comparator's own sentence says "colocalized"; editing it to satisfy this
# project's style rule would falsify the source. Narrow by design — the field
# must be named for quoting, or the surrounding lines must say so — so it cannot
# be reached by writing a claim and calling it a quote. Six lines of lookback,
# because a YAML folded block puts its key well above the text.
QUOTATION_CUES = re.compile(
    r"source_claim|\bquoted\b|\bverbatim\b|their own words|own sentence",
    re.IGNORECASE,
)
QUOTATION_LOOKBACK = 6

# Binary and generated paths carry no prose claims.
SKIP_SUFFIX = {
    ".png", ".pdf", ".h5ad", ".rda", ".RData", ".gz", ".tar", ".xlsx",
    ".sha256", ".lock", ".ipynb",
}

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    emit("P4-T7 — repository language audit")
    emit("=" * 70)
    emit()
    emit("  PROJECT_PLAN §6 P4-T7: grep the entire repo and fix every hit.")
    emit("  A RULE rather than a remembered command (hard constraint 5), so a")
    emit("  violation fails the build instead of surviving until someone looks.")
    emit()

    tracked = subprocess.run(
        ["git", "ls-files"], capture_output=True, text=True, check=True
    ).stdout.split("\n")
    files = [
        f
        for f in tracked
        if f
        and Path(f).suffix not in SKIP_SUFFIX
        and not f.startswith("results/")
    ]
    emit(f"tracked files scanned: {len(files)} "
         f"(of {len([f for f in tracked if f])} tracked)")
    emit("  excluded: binary/generated suffixes, and results/ — a phrase there")
    emit("  is a symptom whose cause is in a script this audit already covers.")
    emit()

    hits = []
    prohibitions = []
    quotations = []
    for path in files:
        try:
            text = Path(path).read_text(encoding="utf-8")
        except (UnicodeDecodeError, FileNotFoundError):
            continue
        lines = text.split("\n")
        for lineno, line in enumerate(lines, start=1):
            for name, pattern, why in PATTERNS:
                m = pattern.search(line)
                if not m:
                    continue
                # Window of the matched line plus the one before it: this
                # project's prose wraps, and the negation often sits on the
                # previous line.
                window = (lines[lineno - 2] + " " if lineno > 1 else "") + line
                quote_window = " ".join(
                    lines[max(0, lineno - 1 - QUOTATION_LOOKBACK) : lineno]
                )
                record = {
                    "file": path,
                    "line": lineno,
                    "pattern": name,
                    "match": m.group(0),
                    "why": why,
                    "text": line.strip()[:200],
                }
                if QUOTATION_CUES.search(quote_window):
                    record["category"] = "quotation"
                    quotations.append(record)
                elif PROHIBITION_CUES.search(window):
                    record["category"] = "prohibition"
                    prohibitions.append(record)
                else:
                    record["category"] = "violation"
                    hits.append(record)

    emit(f"patterns checked: {[n for n, _, _ in PATTERNS]}")
    emit(f"matches: {len(hits) + len(prohibitions) + len(quotations)}")
    emit(f"  prohibitions (the rule being STATED):  {len(prohibitions)}")
    emit(f"  quotations   (someone ELSE'S claim):    {len(quotations)}")
    emit(f"  violations   (the rule being BROKEN):   {len(hits)}")
    emit()
    emit("  PROJECT_PLAN §6 P4-T7 says 'fix every hit'. Run literally that")
    emit("  fails on CLAUDE.md's own hard constraint 6 and on")
    emit("  config/compartment_map.yaml's header — the files that state a rule")
    emit("  have to name what they forbid. A check that cannot tell the two")
    emit("  apart would force the rule's statement to be deleted to make the")
    emit("  check pass, which inverts what it is for.")
    emit()
    if prohibitions:
        emit("prohibitions, by file (expected, and not an error):")
        by_file = {}
        for r in prohibitions:
            by_file.setdefault(r["file"], 0)
            by_file[r["file"]] += 1
        for f, c in sorted(by_file.items()):
            emit(f"  {c:3d}  {f}")
        emit()
        emit("  A repository whose prohibition count falls to zero has probably")
        emit("  lost its rules rather than fixed its language.")
        emit()

    if hits:
        emit("VIOLATIONS — every one must be fixed, not suppressed")
        emit("-" * 70)
        for h in hits:
            emit(f"  {h['file']}:{h['line']}  [{h['pattern']}]  "
                 f"matched {h['match']!r}")
            emit(f"    {h['text']}")
            emit(f"    -> {h['why']}")
            emit()

    report = {
        "task": "P4-T7",
        "n_files_scanned": len(files),
        "n_violations": len(hits),
        "n_prohibitions": len(prohibitions),
        "n_quotations": len(quotations),
        "patterns": [n for n, _, _ in PATTERNS],
        "violations": hits,
        "prohibitions": prohibitions,
        "quotations": quotations,
        "classification_rule": (
            "A match inside a prohibition — the line, or the one before it, "
            "carrying never / not / forbid / avoid / rather than / instead of / "
            "cannot — is the rule being STATED and is counted separately. "
            "Everything else fails the build. PROJECT_PLAN §6's literal 'fix "
            "every hit' would fail on CLAUDE.md's own hard constraint 6: the "
            "files that state a rule have to name what they forbid. A match "
            "inside a declared quotation field is a third category: the P4-T6 "
            "comparator's own sentence says 'colocalized', and editing a "
            "quotation to satisfy this project's style rule would falsify the "
            "source (ADR 0021 §2)."
        ),
        "scope": (
            "git ls-files, excluding binary/generated suffixes and results/. "
            "What is committed is what makes a claim."
        ),
        "rule": (
            "This assay carries no coordinates. The supported phrase is "
            "'inferred crosstalk between adjacent compartments' (hard "
            "constraint 6). PanCK was the only collection mask, so a "
            "compartment label is never a cell-type label (Q2)."
        ),
        "clean": not hits,
    }
    Path(snakemake.output.report).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    emit(f"wrote {snakemake.output.report}")

    if hits:
        raise RuntimeError(
            f"{len(hits)} language violation(s) — see {snakemake.log[0]}. "
            "This assay has no coordinates: the phrase is 'inferred crosstalk "
            "between adjacent compartments' (hard constraint 6), and a "
            "compartment label is never a cell-type label (Q2). Fix each hit; "
            "do not widen the exclusions."
        )

    emit()
    emit("CLEAN. No committed file claims a spatial relationship this assay")
    emit("cannot support, and no compartment label is used as a cell-type label.")
    emit(f"({len(prohibitions)} matches are the rule being stated and "
         f"{len(quotations)} are quoted from an external source.)")
