"""P4-T1 — build the two paired adjacency sets, and pre-register how duplicates collapse.

Owner task: P4-T1. Driven by rule p4t1_build_pairs.

Phase 4 correlates a ligand measured in a tumour AOI against a receptor measured
in the SAME PATIENT'S immune AOI, across patients, within site. That requires a
patient-level pairing, and two of the thirteen lung patients contribute two
`TIME-L` AOIs rather than one. How those collapse is a scientific decision, not
an implementation detail — PROJECT_PLAN §6 P4-T1 says "write it in an ADR before
running anything" — and this script implements ADR 0021 §1 rather than choosing.

  PRIMARY      mean of log2(q3 + 1) across a patient's duplicate AOIs
  SENSITIVITY  the AOI with the highest gene_detection_rate (qc_metrics.tsv)

Both are built here, both are written, and neither is a tiebreaker for the
other. The choice moves 2 of 13 lung patients and 0 of 8 brain ones, so it
cannot be load-bearing — and the sensitivity table exists so that claim is
checkable rather than asserted.

**This is flag-don't-drop applied to a collapse rather than an exclusion.** No
AOI is discarded: both of P12's and both of P24's `TIME-L` AOIs contribute to
the primary, and `results/tables/qc_excluded.tsv` gains nothing from Phase 4.

Three things it deliberately does NOT do:

  * **Key on `compartment`.** `compartment` is degenerate across sites — `L`,
    `LB` and `mLN` are all `tumour` — so grouping by it would pool lung and
    brain tumour AOIs and destroy the contrast Phase 4 measures. Everything
    here keys on `aoi_code`.

  * **Pair `TBME`.** It is the glial side of a brain adjacency and it is not one
    of Phase 4's compartments: 20 of 20 of its AOIs sit in one DSP run, so batch
    is inseparable from biology there (Q3, ADR 0009 §3). ADR 0014 §3 records the
    decision and re-scoping it is a stop-and-ask. P15's duplicate AOI is `TBME`,
    which is why it does not appear below as a duplicate.

  * **Discover the cohort size.** `crosstalk.pairing.expected_n_patients` is
    ASSERTED against the data. Lung 13 and brain 8 come from
    `results/tables/design_matrix.tsv` and PROJECT_PLAN §2.1; if either count
    moves, the design changed and the phase stops rather than adapting to it.
    `TIME-B` n = 8 is the binding constraint on the entire project and every
    brain claim states it inline (hard constraint 8).
"""

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    adjacencies = snakemake.params.adjacencies
    pairing = snakemake.params.pairing
    primary_rule = pairing["duplicate_rule"]
    sensitivity_rule = pairing["duplicate_sensitivity"]
    expected = pairing["expected_n_patients"]

    emit("P4-T1 — paired adjacency sets (ADR 0021 §1, pre-registered)")
    emit("=" * 70)
    emit()
    emit(f"  primary collapse rule      {primary_rule}")
    emit(f"  sensitivity collapse rule  {sensitivity_rule}")
    emit(f"  expected patients          {expected}")
    emit()
    emit("  Aim A5 is EXPLORATORY (ADR 0014). No Phase 4 result may be a")
    emit("  headline claim; the anticipated null is an assay-sensitivity limit,")
    emit("  never evidence that the crosstalk is absent.")
    emit()

    if primary_rule == sensitivity_rule:
        raise RuntimeError(
            "crosstalk.pairing.duplicate_rule and duplicate_sensitivity are "
            f"both '{primary_rule}'. The sensitivity exists to run the REJECTED "
            "option (ADR 0021 §1); making it the same as the primary would "
            "report one choice twice and quietly retire the check."
        )

    adata = ad.read_h5ad(snakemake.input.h5ad)
    obs = adata.obs.copy()
    obs["aoi_label"] = obs.index.astype(str)
    emit(f"AOIs: {len(obs)} (from {snakemake.input.h5ad})")

    qc = pd.read_csv(snakemake.input.qc, sep="\t", float_precision="round_trip")
    qc = qc.set_index("aoi_label")
    emit(f"QC metrics: {len(qc)} AOIs (from {snakemake.input.qc})")
    emit()

    # --- build the pairing, one adjacency at a time -------------------------
    rows = []
    for site in sorted(adjacencies):
        codes = adjacencies[site]
        tumour_code, immune_code = codes["tumour"], codes["immune"]
        emit(f"{site}: {tumour_code} <-> {immune_code}")

        side_aois = {}
        for side, code in (("tumour", tumour_code), ("immune", immune_code)):
            sub = obs.loc[obs["aoi_code"] == code]
            if sub.empty:
                raise RuntimeError(
                    f"aoi_code '{code}' matched no AOI. "
                    "crosstalk.adjacencies keys on aoi_code, never on "
                    "compartment — check config/config.yaml against the design "
                    "table in CLAUDE.md."
                )
            side_aois[side] = sub
            emit(f"  {side:7s} {code:7s} {len(sub):3d} AOIs, "
                 f"{sub['patient_id'].nunique():2d} patients")

        shared = sorted(
            set(side_aois["tumour"]["patient_id"])
            & set(side_aois["immune"]["patient_id"]),
            key=lambda p: int(str(p).lstrip("P")) if str(p)[1:].isdigit() else 0,
        )
        emit(f"  paired patients: {len(shared)} — {shared}")

        want = expected[site]
        if len(shared) != want:
            raise RuntimeError(
                f"{site}: found {len(shared)} paired patients, expected {want} "
                f"(crosstalk.pairing.expected_n_patients). The count is "
                "ASSERTED, not discovered: if the design changed, that is a "
                "stop-and-ask, not a number to edit down. Paired set was "
                f"{shared}."
            )

        for patient in shared:
            row = {"site": site, "patient_id": patient}
            for side, code in (("tumour", tumour_code), ("immune", immune_code)):
                sub = side_aois[side]
                mine = sub.loc[sub["patient_id"] == patient]
                labels = sorted(mine["aoi_label"])
                rates = qc.loc[labels, "gene_detection_rate"]
                best = str(rates.idxmax())
                row[f"{side}_aoi_code"] = code
                row[f"{side}_aoi_labels"] = ";".join(labels)
                row[f"{side}_n_aoi"] = len(labels)
                row[f"{side}_is_duplicate"] = len(labels) > 1
                row[f"{side}_highest_detection_aoi"] = best
                row[f"{side}_gene_detection_rate_min"] = float(rates.min())
                row[f"{side}_gene_detection_rate_max"] = float(rates.max())
                row[f"{side}_dsp_run"] = ";".join(sorted(set(mine["dsp_run"])))
                row[f"{side}_qc_flag_any"] = bool(mine["qc_flag"].any())
            rows.append(row)
        emit()

    pairs = pd.DataFrame(rows)

    # --- duplicates, named rather than counted ------------------------------
    dup_mask = pairs["tumour_is_duplicate"] | pairs["immune_is_duplicate"]
    dups = pairs.loc[dup_mask]
    emit("duplicate AOIs — the only thing ADR 0021 §1 decides:")
    if dups.empty:
        emit("  none")
    for r in dups.itertuples(index=False):
        for side in ("tumour", "immune"):
            labels = getattr(r, f"{side}_aoi_labels").split(";")
            if len(labels) > 1:
                rates = qc.loc[labels, "gene_detection_rate"]
                emit(
                    f"  {r.site:5s} {r.patient_id:4s} {side:7s} "
                    f"{getattr(r, f'{side}_aoi_code'):7s} {labels}  "
                    f"detection {[round(float(x), 4) for x in rates]}  "
                    f"-> sensitivity picks {rates.idxmax()}"
                )
    n_dup_lung = int(dup_mask[pairs["site"] == "lung"].sum())
    n_dup_brain = int(dup_mask[pairs["site"] == "brain"].sum())
    emit()
    emit(
        f"  the collapse rule moves {n_dup_lung} of "
        f"{int((pairs['site'] == 'lung').sum())} lung patients and "
        f"{n_dup_brain} of {int((pairs['site'] == 'brain').sum())} brain ones, "
        "so it cannot be load-bearing (ADR 0021 §1)."
    )
    emit()

    # --- the sensitivity table: same patients, one AOI per side -------------
    sensitivity = pairs.copy()
    for side in ("tumour", "immune"):
        sensitivity[f"{side}_aoi_labels"] = sensitivity[
            f"{side}_highest_detection_aoi"
        ]
        sensitivity[f"{side}_n_aoi"] = 1
    sensitivity["collapse_rule"] = sensitivity_rule
    pairs["collapse_rule"] = primary_rule

    # An AOI reachable from the primary but not from the sensitivity is the
    # whole of what the two tables differ by; state it as a number rather than
    # leaving a reader to diff two files.
    def _aoi_set(frame):
        out = set()
        for side in ("tumour", "immune"):
            for cell in frame[f"{side}_aoi_labels"]:
                out.update(cell.split(";"))
        return out

    primary_aois = _aoi_set(pairs)
    sensitivity_aois = _aoi_set(sensitivity)
    dropped = sorted(primary_aois - sensitivity_aois)
    emit(f"AOIs used by the primary:     {len(primary_aois)}")
    emit(f"AOIs used by the sensitivity: {len(sensitivity_aois)}")
    emit(f"  reachable from the primary only: {dropped}")
    emit(
        "  Those AOIs are NOT excluded from the project — the primary averages "
        "them in. Flag-don't-drop; qc_excluded.tsv is untouched."
    )
    emit()

    # --- assertions ---------------------------------------------------------
    n_checks = 0

    for frame, label in ((pairs, "primary"), (sensitivity, "sensitivity")):
        if frame.duplicated(["site", "patient_id"]).any():
            raise RuntimeError(
                f"{label}: a patient appears twice within one adjacency. The "
                "pairing is one row per patient per site by construction."
            )
        n_checks += 1

    for side in ("tumour", "immune"):
        if (sensitivity[f"{side}_n_aoi"] != 1).any():
            raise RuntimeError(
                f"sensitivity: {side} side has a row with more than one AOI. "
                "highest_detection_rate must resolve to exactly one AOI."
            )
        n_checks += 1

    # Every AOI named must exist, and must carry the aoi_code claimed for it.
    known = set(obs["aoi_label"])
    for frame, label in ((pairs, "primary"), (sensitivity, "sensitivity")):
        for side in ("tumour", "immune"):
            for _, r in frame.iterrows():
                for aoi in r[f"{side}_aoi_labels"].split(";"):
                    if aoi not in known:
                        raise RuntimeError(f"{label}: unknown AOI '{aoi}'.")
                    got = str(obs.loc[obs["aoi_label"] == aoi, "aoi_code"].iloc[0])
                    if got != r[f"{side}_aoi_code"]:
                        raise RuntimeError(
                            f"{label}: AOI {aoi} is aoi_code '{got}', not "
                            f"'{r[f'{side}_aoi_code']}'."
                        )
        n_checks += 1

    # Cross-check against design_matrix.tsv — an INDEPENDENT count, written by
    # p0t3_parse_samples, not by anything in this file.
    design = pd.read_csv(snakemake.input.design, sep="\t", float_precision="round_trip")
    design = design.loc[design["patient_id"] != "TOTAL"].set_index("patient_id")
    for site in sorted(adjacencies):
        codes = adjacencies[site]
        want = set(
            design.index[
                (design[codes["tumour"]] > 0) & (design[codes["immune"]] > 0)
            ]
        )
        got = set(pairs.loc[pairs["site"] == site, "patient_id"])
        if want != got:
            raise RuntimeError(
                f"{site}: pairing disagrees with design_matrix.tsv. "
                f"only in design: {sorted(want - got)}; "
                f"only here: {sorted(got - want)}. Two independent counts of "
                "the same design must agree."
            )
    n_checks += 1

    emit(f"assertions: {n_checks} of {n_checks} passed "
         "(one row per patient per site; sensitivity resolves to one AOI; "
         "AOI labels and codes exist; design_matrix.tsv agreement)")
    emit()

    pairs.to_csv(snakemake.output.pairs, sep="\t", index=False)
    sensitivity.to_csv(snakemake.output.sensitivity, sep="\t", index=False)

    summary = {
        "task": "P4-T1",
        "aim": "A5",
        "exploratory": True,
        "exploratory_adr": "ADR 0014",
        "pre_registered": "ADR 0021 §1",
        "duplicate_rule": primary_rule,
        "duplicate_sensitivity": sensitivity_rule,
        "adjacencies": {
            site: dict(adjacencies[site]) for site in sorted(adjacencies)
        },
        "n_patients": {
            site: int((pairs["site"] == site).sum()) for site in sorted(adjacencies)
        },
        "expected_n_patients": dict(expected),
        "patients": {
            site: list(pairs.loc[pairs["site"] == site, "patient_id"])
            for site in sorted(adjacencies)
        },
        "n_patients_with_duplicate": {"lung": n_dup_lung, "brain": n_dup_brain},
        "duplicate_patients": sorted(set(dups["patient_id"])) if not dups.empty else [],
        "n_aoi_primary": len(primary_aois),
        "n_aoi_sensitivity": len(sensitivity_aois),
        "aoi_primary_only": dropped,
        "n_time_b": int((pairs["site"] == "brain").sum()),
        "power_floor_sd": "1.1-1.3",
        "reporting_rule": (
            "TIME-B n = 8. No A5 result may be a headline claim (ADR 0014); "
            "every nomination states the per-site detection of BOTH partners; "
            "a null is an assay-sensitivity limit, never evidence that the "
            "crosstalk is absent. Never write 'colocalisation' or 'spatially "
            "adjacent' — the phrase is 'inferred crosstalk between adjacent "
            "compartments' (hard constraint 6)."
        ),
        "n_assertions": n_checks,
        "h5ad": str(snakemake.input.h5ad),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit(
        f"wrote {snakemake.output.pairs}, {snakemake.output.sensitivity}, "
        f"{snakemake.output.summary}"
    )
    emit()
    emit(
        "REPORTING RULE — brain is TIME-B n = 8, the binding constraint on the "
        "project (hard constraint 8). A5 is exploratory (ADR 0014)."
    )
