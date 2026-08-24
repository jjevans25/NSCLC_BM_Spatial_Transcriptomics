"""Parse GEO sample titles into samples.tsv, plus the design and batch tables.

Owner task: P0-T3. Driven by rule p0t3_parse_samples.

The title grammar is PROJECT_PLAN §A.1 and it has two forms:

    Patient 5 [TIME-B05]          tumour-bearing cases, numbered 1..44
    Non-tumor brain #1 [BC03]     the seven non-tumour brain controls

Replicate AOIs carry a trailing letter (`TIME-L12a`, `TBME15b`).

`L` MUST come last in the code alternation. Put it first and it matches the
leading `L` of `LB`, silently relabelling 27 brain-metastasis AOIs as primary
lung tumour — §A.1 calls this the single most likely bug in Phase 0, and it
would be invisible downstream because the counts would still sum to 120.

Everything this script asserts is a fail-loud check, not a warning. A design
table that silently disagrees with PROJECT_PLAN §2.1 is worse than a crash.
"""

import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import yaml

# --- title grammar (PROJECT_PLAN §A.1, verbatim) ---------------------------
PATTERN = (
    r"^(?:Patient\s+(?P<patient>\d+)"
    r"|Non-tumor\s+brain\s+#(?P<control>\d+))"
    r"\s+\[(?P<code>TIME-L|TIME-B|TBME|mLN|LB|BC|L)"
    r"(?P<num>\d+)(?P<rep>[a-z])?\]$"
)
TITLE_RX = re.compile(PATTERN)

# GSM6573697_DSP-1012300141221-A-A02.dcc.gz
DCC_RX = re.compile(r"^(GSM\d+)_(DSP-\d+)-([A-Z])-([A-Z]\d+)\.dcc\.gz$")

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)


class DesignError(RuntimeError):
    """A structural disagreement with PROJECT_PLAN §2.1."""


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    # ---------------------------------------------------------------- inputs
    compartment_map = yaml.safe_load(
        Path(snakemake.input.compartment_map).read_text(encoding="utf-8")
    )

    with gzip.open(snakemake.input.soft, "rt", encoding="utf-8") as handle:
        soft_text = handle.read()

    # ------------------------------------------------- parse the SOFT records
    records = []
    for block in soft_text.split("^SAMPLE = ")[1:]:
        gsm = block.split("\n", 1)[0].strip()
        title = re.search(r"!Sample_title = (.*)", block).group(1).strip()
        cell_type = re.search(r"!Sample_characteristics_ch1 = (.*)", block).group(1).strip()
        records.append((gsm, title, cell_type))

    emit(f"SOFT records: {len(records)}")

    unmatched = [(g, t) for g, t, _ in records if not TITLE_RX.match(t)]
    if unmatched:
        detail = "\n".join(f"  {g}  {t!r}" for g, t in unmatched)
        raise DesignError(
            f"{len(unmatched)} sample title(s) do not match the §A.1 grammar:\n{detail}\n"
            "Fix the pattern, do not special-case the row."
        )

    rows = []
    for gsm, title, cell_type in records:
        f = TITLE_RX.match(title).groupdict()
        is_control = f["control"] is not None

        # Namespace the grouping key. `Non-tumor brain #1 [BC03]` would otherwise
        # yield patient_id "3" from the code suffix and collide with NSCLC
        # Patient 3 — the random intercept would then pool an unrelated control
        # and a tumour case as one subject. The BC suffix is a TMA/core number,
        # not a subject id; the `#n` is the subject.
        patient_id = f"BC{int(f['control'])}" if is_control else f"P{int(f['patient'])}"

        rows.append(
            {
                "gsm_id": gsm,
                "sample_title": title,
                # The bracketed label is the expression-matrix column header
                # (`L01`, `TIME-L12a`), so it is the join key to the data.
                "aoi_label": f"{f['code']}{f['num']}{f['rep'] or ''}",
                "aoi_code": f["code"],
                "patient_id": patient_id,
                "is_control": is_control,
                "cell_type": cell_type,
                "_num": int(f["num"]),
                "_rep": f["rep"],
            }
        )

    samples = pd.DataFrame(rows)

    # ------------------------------------------- assert: suffix == patient id
    tumour = samples[~samples["is_control"]]
    mism = tumour[tumour["_num"] != tumour["patient_id"].str.lstrip("P").astype(int)]
    if len(mism):
        detail = "\n".join(f"  {r.gsm_id}  {r.sample_title}" for r in mism.itertuples())
        raise DesignError(
            f"{len(mism)} AOI code suffix(es) disagree with the patient number in "
            f"the same title — a GEO transcription error:\n{detail}"
        )
    emit(f"suffix == patient number: OK ({len(tumour)} non-control rows)")

    # -------------------------------------------- annotate from the compartment map
    unknown = sorted(set(samples["aoi_code"]) - set(compartment_map))
    if unknown:
        raise DesignError(f"AOI code(s) absent from compartment_map.yaml: {unknown}")

    samples["site"] = samples["aoi_code"].map(lambda c: compartment_map[c]["site"])
    samples["compartment"] = samples["aoi_code"].map(
        lambda c: compartment_map[c]["compartment"]
    )

    # --------------------------------------------------------- replicate flag
    # Schema definition: true where a patient contributes more than one AOI of
    # this code. The trailing letter in the title is the authors' own marking of
    # the same thing, so the two must agree — if they ever disagree, one of them
    # is wrong and silently picking either would hide it.
    pair_counts = Counter(zip(samples["patient_id"], samples["aoi_code"]))
    samples["replicate_flag"] = [
        pair_counts[(p, c)] > 1
        for p, c in zip(samples["patient_id"], samples["aoi_code"])
    ]
    by_letter = samples["_rep"].notna()
    if samples["replicate_flag"].tolist() != by_letter.tolist():
        disagree = samples[samples["replicate_flag"] != by_letter]
        detail = "\n".join(f"  {r.gsm_id}  {r.sample_title}" for r in disagree.itertuples())
        raise DesignError(
            "replicate_flag disagrees with the trailing-letter marking:\n"
            f"{detail}\n"
            "Either a replicate is unlabelled or a label is spurious."
        )
    n_rep = int(samples["replicate_flag"].sum())
    emit(f"replicate rows: {n_rep} (letter marking and duplicate-pair rule agree)")

    # ---------------------------------------------------- DSP run from filelist
    # GEO exposes no batch field (Q3); the DSP run id survives only in the DCC
    # filenames. See docs/data-provenance.md.
    dsp = {}
    for line in Path(snakemake.input.filelist).read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) < 5 or fields[-1] != "DCC":
            continue
        m = DCC_RX.match(fields[1])
        if not m:
            raise DesignError(f"unparsable DCC filename in filelist: {fields[1]!r}")
        dsp[m.group(1)] = (m.group(2), m.group(4))

    missing = sorted(set(samples["gsm_id"]) - set(dsp))
    extra = sorted(set(dsp) - set(samples["gsm_id"]))
    if missing or extra:
        raise DesignError(
            f"GSM sets differ between SOFT and filelist — missing {missing}, extra {extra}"
        )
    samples["dsp_run"] = samples["gsm_id"].map(lambda g: dsp[g][0])
    samples["dsp_well"] = samples["gsm_id"].map(lambda g: dsp[g][1])
    emit(f"DSP runs: {dict(Counter(samples['dsp_run']))}")

    # ------------------------------------------------------ assert: the counts
    observed = Counter(samples["aoi_code"])
    expected = {code: spec["expected_aoi"] for code, spec in compartment_map.items()}
    wrong = {c: (observed[c], expected[c]) for c in expected if observed[c] != expected[c]}
    if wrong:
        detail = "\n".join(f"  {c}: got {g}, expected {e}" for c, (g, e) in sorted(wrong.items()))
        raise DesignError(
            "per-compartment counts disagree with PROJECT_PLAN §2.1:\n"
            f"{detail}\n"
            "§2.1 is authoritative — fix the parser, do not adjust the expectation."
        )

    n_expected = snakemake.params.expected_n_aoi
    if len(samples) != n_expected:
        raise DesignError(f"expected {n_expected} AOIs, parsed {len(samples)}")

    if samples["aoi_label"].duplicated().any():
        dups = sorted(samples.loc[samples["aoi_label"].duplicated(), "aoi_label"])
        raise DesignError(f"duplicate AOI labels (not a usable join key): {dups}")

    # --------------------------------------- assert: the join key actually joins
    # The expression matrix is keyed by AOI label (`L01`, `TIME-L12a`), not by
    # GSM id, so aoi_label is the only bridge between annotation and data. Read
    # just the header and prove the correspondence is exact and one-to-one — a
    # partial match would surface much later as unexplained NaN columns.
    with gzip.open(snakemake.input.matrix, "rt", encoding="utf-8") as handle:
        matrix_columns = handle.readline().rstrip("\n").split("\t")[1:]

    only_matrix = sorted(set(matrix_columns) - set(samples["aoi_label"]))
    only_samples = sorted(set(samples["aoi_label"]) - set(matrix_columns))
    repeated = sorted({c for c in matrix_columns if matrix_columns.count(c) > 1})
    if only_matrix or only_samples or repeated:
        raise DesignError(
            "aoi_label does not join one-to-one onto the expression matrix "
            "columns:\n"
            f"  only in matrix : {only_matrix}\n"
            f"  only in samples: {only_samples}\n"
            f"  duplicated     : {repeated}"
        )
    emit(f"aoi_label <-> matrix columns: 1:1 OK ({len(matrix_columns)} columns)")

    emit("\nper-compartment counts (all match §2.1):")
    for code in sorted(expected):
        emit(f"  {code:<7} {observed[code]:>3}")
    emit(f"  {'TOTAL':<7} {len(samples):>3}")

    # ------------------------------------------------------------ samples.tsv
    columns = [
        "gsm_id", "sample_title", "aoi_label", "aoi_code", "patient_id",
        "site", "compartment", "replicate_flag", "is_control",
        "dsp_run", "dsp_well", "cell_type",
    ]
    out = samples[columns].sort_values(["aoi_code", "patient_id", "aoi_label"])
    out.to_csv(snakemake.output.samples, sep="\t", index=False)
    emit(f"\nwrote {snakemake.output.samples} ({len(out)} rows)")

    # ----------------------------------------------- design matrix: patients x codes
    code_order = [c for c in ["L", "LB", "mLN", "TBME", "TIME-L", "TIME-B", "BC"]
                  if c in expected]
    design = (
        pd.crosstab(samples["patient_id"], samples["aoi_code"])
        .reindex(columns=code_order, fill_value=0)
    )
    # Patients first (numeric order), then controls — not lexicographic, which
    # would read P1, P10, P11, P12, P2.
    design = design.reindex(
        sorted(design.index, key=lambda p: (p.startswith("BC"), int(re.sub(r"\D", "", p))))
    )
    design["total"] = design.sum(axis=1)
    design.loc["TOTAL"] = design.sum(axis=0)
    design.to_csv(snakemake.output.design, sep="\t")
    emit(f"wrote {snakemake.output.design} ({len(design) - 1} subjects)")

    # ------------------------------------------- batch crosstab: codes x DSP run
    batch = pd.crosstab(samples["aoi_code"], samples["dsp_run"]).reindex(
        index=code_order, fill_value=0
    )
    batch["total"] = batch.sum(axis=1)
    batch.loc["TOTAL"] = batch.sum(axis=0)
    batch.to_csv(snakemake.output.batch, sep="\t")
    emit(f"wrote {snakemake.output.batch}")

    # Which compartments sit entirely in one run cannot have batch separated
    # from biology. Reported, not enforced — it is a limitation, not an error.
    runs = [c for c in batch.columns if c != "total"]
    confounded = [
        code for code in code_order
        if sum(batch.loc[code, r] > 0 for r in runs) == 1
    ]
    emit("\nDSP run x compartment:")
    emit(batch.to_string())
    if confounded:
        emit(
            "\nNOTE: these compartments sit entirely within one DSP run, so batch "
            f"and biology cannot be separated for them: {', '.join(confounded)}. "
            "Carry this into the P0-T5 QC model and docs/limitations.md (P0-T8)."
        )

    # ------------------------------------------------------------- provenance
    Path(snakemake.output.summary).write_text(
        json.dumps(
            {
                "n_aoi": int(len(samples)),
                "n_subjects": int(samples["patient_id"].nunique()),
                "n_patients": int(samples.loc[~samples["is_control"], "patient_id"].nunique()),
                "n_controls": int(samples.loc[samples["is_control"], "patient_id"].nunique()),
                "counts_by_code": {c: int(observed[c]) for c in code_order},
                "replicate_aois": int(n_rep),
                "dsp_runs": {k: int(v) for k, v in Counter(samples["dsp_run"]).items()},
                "compartments_confounded_with_dsp_run": confounded,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    emit(f"wrote {snakemake.output.summary}")
