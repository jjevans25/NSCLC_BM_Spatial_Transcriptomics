"""P3-T4 — the within-patient direction check for the checkpoint panel.

Owner task: P3-T4. Driven by rule p3t4_checkpoint_paired_check.

A NEW FILE rather than an edit to paired_check.py, deliberately: editing that
would fire the code rerun-trigger on p2t4_paired_check and re-run Phase 2 for no
numerical gain. Same reason crosscheck_checkpoint_mixedlm.py is its own file.

Restricted to patients contributing an AOI at BOTH sites within a compartment,
the lung-vs-brain difference is measured *within* a patient, so it is free of
between-patient confounding by construction. That is the whole of what this task
is for.

**DIRECTION OF EFFECT ONLY. NO TEST STATISTIC, IN EITHER COMPARTMENT.**

The immune compartment has 5 paired patients and PROJECT_PLAN §2.3 forbids a
p-value there in terms — "a consistency check, never a headline result".

The tumour compartment has 23, which would support one. **It still does not get
one.** A test introduced at the moment it becomes available, on a contrast the
phase never pre-registered a test for, is the same move ADR 0018 rejected in
Option B — and it would arrive with no multiplicity family, no power statement,
and no ADR. If a paired tumour test is wanted it is a stop-and-ask, decided
before it is computed, not a column that appears because the n happened to be
large enough.

So this script computes no test statistic at all: there is no p-value column to
be tempted by, and no significance to read into a sign count. The assertion at
the end enforces that mechanically rather than by intention.

A patient contributing more than one AOI at a site has those AOIs averaged
before differencing, so each patient contributes exactly one delta and nobody is
weighted by how many AOIs they happen to have. P12 contributes two `TIME-L`.

Concordance is measured against the unpaired fit for the same gene x
compartment: the PRIMARY estimate where one exists, otherwise the EXPLORATORY
one, with `compared_against` saying which. A gene the primary declined to model
is not silently compared against a q-valued number it never had. Paired and
unpaired results share AOIs and are not independent evidence — agreement is
reassurance about direction, not confirmation.

Every row carries its per-site detection counts (ADR 0008): a delta computed
from values that are mostly background is a number, not a measurement, and the
counts are what let a reader tell the difference.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

FMT = "%.17g"

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    expected_immune = snakemake.params.expected_immune_paired

    expression = pd.read_csv(
        snakemake.input.expression, sep="\t", float_precision="round_trip"
    )
    primary = pd.read_csv(
        snakemake.input.primary, sep="\t", float_precision="round_trip"
    )
    exploratory = pd.read_csv(
        snakemake.input.exploratory, sep="\t", float_precision="round_trip"
    )

    emit("P3-T4 — within-patient direction check, checkpoint panel")
    emit("DIRECTION ONLY. No test statistic in either compartment — see the")
    emit("module docstring for why the 23-patient tumour set does not get one.")
    emit()

    # --------------------------------------------------------- paired patients
    paired_by_compartment = {}
    for compartment in sorted(expression["compartment"].unique()):
        block = expression.loc[expression["compartment"] == compartment]
        per_patient = block.groupby("patient_id")["site"].nunique()
        paired = sorted(per_patient[per_patient == 2].index)
        paired_by_compartment[compartment] = paired
        emit(f"{compartment:7s}: {len(paired)} patients with AOIs at both sites")
        emit(f"         {paired}")

    n_immune = len(paired_by_compartment.get("immune", []))
    if n_immune != expected_immune:
        raise RuntimeError(
            f"expected {expected_immune} paired immune patients (PROJECT_PLAN "
            f"§2.3), found {n_immune}: {paired_by_compartment.get('immune')}. "
            "Do not adjust the expectation to match the data."
        )
    emit()
    emit(f"immune paired n = {n_immune}, as PROJECT_PLAN §2.3 states.")
    emit("TIME-B n = 8 overall; the paired subset is smaller still.")
    emit()

    # ------------------------------------------------------------------ deltas
    delta_frames = []
    for compartment, paired in paired_by_compartment.items():
        block = expression.loc[
            (expression["compartment"] == compartment)
            & (expression["patient_id"].isin(paired))
        ]
        # One value per patient x site x gene: a patient with two AOIs at a site
        # must not count twice.
        per_site = (
            block.groupby(["gene", "patient_id", "site"], observed=True)["expression"]
            .mean()
            .unstack("site")
        )
        n_collapsed = len(block) - len(per_site) * 2
        emit(f"{compartment}: {len(block)} AOI rows -> {len(per_site)} patient x "
             f"gene cells ({n_collapsed} collapsed by averaging)")

        # Detection counts for the same cells, so a delta is never read without
        # knowing how much of it is above background (ADR 0008).
        detected = (
            block.groupby(["gene", "patient_id", "site"], observed=True)["detected"]
            .max()
            .unstack("site")
        )

        frame = per_site.reset_index()
        frame["compartment"] = compartment
        frame["delta_brain_minus_lung"] = frame["brain"] - frame["lung"]
        det = detected.reset_index()
        frame["detected_lung"] = det["lung"].to_numpy()
        frame["detected_brain"] = det["brain"].to_numpy()
        delta_frames.append(frame)

    deltas = pd.concat(delta_frames, ignore_index=True)
    deltas = deltas[
        [
            "compartment", "gene", "patient_id", "lung", "brain",
            "delta_brain_minus_lung", "detected_lung", "detected_brain",
        ]
    ].sort_values(["compartment", "gene", "patient_id"])
    deltas.to_csv(snakemake.output.deltas, sep="\t", index=False, float_format=FMT)
    emit()

    # ------------------------------------------------------------- concordance
    prim = primary.loc[primary["model"] == "primary"].set_index(
        ["stratum", "gene"]
    )
    expl = exploratory.loc[exploratory["model"] == "primary"].set_index(
        ["stratum", "gene"]
    )

    rows = []
    for (compartment, gene), g in deltas.groupby(
        ["compartment", "gene"], observed=True
    ):
        delta = g["delta_brain_minus_lung"].to_numpy()
        mean_delta = float(np.mean(delta))

        key = (compartment, gene)
        estimate, source = np.nan, "none"
        if key in prim.index and pd.notna(prim.loc[key, "estimate"]):
            estimate = float(prim.loc[key, "estimate"])
            source = "primary"
        elif key in expl.index and pd.notna(expl.loc[key, "estimate"]):
            estimate = float(expl.loc[key, "estimate"])
            source = "exploratory"

        # np.sign(0.0) is 0.0, which would silently read as "agrees" against a
        # zero estimate; require both non-zero and matching.
        agrees = (
            bool(
                np.sign(mean_delta) == np.sign(estimate)
                and mean_delta != 0
                and estimate != 0
            )
            if source != "none"
            else None
        )

        rows.append(
            {
                "compartment": compartment,
                "gene": gene,
                "n_patients": int(len(delta)),
                "n_brain_higher": int((delta > 0).sum()),
                "n_lung_higher": int((delta < 0).sum()),
                "mean_delta_brain_minus_lung": mean_delta,
                "median_delta_brain_minus_lung": float(np.median(delta)),
                "paired_direction": "brain_higher" if mean_delta > 0 else "lung_higher",
                "n_patients_detected_lung": int(g["detected_lung"].sum()),
                "n_patients_detected_brain": int(g["detected_brain"].sum()),
                "unpaired_estimate": estimate,
                "unpaired_direction": (
                    "brain_higher" if estimate > 0 else "lung_higher"
                    if source != "none" else ""
                ),
                "compared_against": source,
                "direction_agrees": agrees,
            }
        )

    concordance = pd.DataFrame(rows).sort_values(["compartment", "gene"])

    forbidden = {"p", "p_raw", "p_value", "pvalue", "q", "q_bh", "statistic",
                 "t_value", "ci_low", "ci_high"}
    assert not (forbidden & set(concordance.columns)), (
        "a test statistic or interval reached the paired table; this check "
        "reports direction only (PROJECT_PLAN §2.3), in both compartments"
    )

    concordance.to_csv(
        snakemake.output.concordance, sep="\t", index=False, float_format=FMT
    )

    # ---------------------------------------------------------------- report
    emit("within-patient direction vs. the unpaired fit:")
    for compartment in sorted(concordance["compartment"].unique()):
        sub = concordance.loc[concordance["compartment"] == compartment]
        n_pat = int(sub["n_patients"].iloc[0])
        emit()
        emit(f"  compartment = {compartment}  (n = {n_pat} paired patients)")
        emit(f"    {'gene':8s} {'paired':>8s} {'brain>lung':>11s} "
             f"{'unpaired':>9s} {'source':>12s}  detected L/B  direction")
        for _, r in sub.iterrows():
            verdict = (
                "agrees" if r["direction_agrees"] is True
                else "DISAGREES" if r["direction_agrees"] is False
                else "no unpaired fit"
            )
            est = (
                f"{r['unpaired_estimate']:+9.3f}"
                if r["compared_against"] != "none" else f"{'—':>9s}"
            )
            emit(
                f"    {r['gene']:8s} {r['mean_delta_brain_minus_lung']:+8.3f} "
                f"{r['n_brain_higher']:>6d}/{r['n_patients']:<4d} {est} "
                f"{r['compared_against']:>12s}  "
                f"{r['n_patients_detected_lung']:>2d}/{r['n_patients_detected_brain']:<2d}"
                f"          {verdict}"
            )

    comparable = concordance.loc[concordance["direction_agrees"].notna()]
    n_agree = int(comparable["direction_agrees"].sum())
    emit()
    emit(f"direction agrees in {n_agree} of {len(comparable)} comparable "
         f"gene x compartment cells "
         f"({len(concordance) - len(comparable)} had no unpaired fit)")
    for compartment in sorted(comparable["compartment"].unique()):
        sub = comparable.loc[comparable["compartment"] == compartment]
        emit(f"  {compartment}: {int(sub['direction_agrees'].sum())} of {len(sub)}")

    emit()
    emit("NO P-VALUE IS REPORTED HERE AND NONE SHOULD BE ADDED — not for the")
    emit("5-patient immune set (PROJECT_PLAN §2.3 forbids it) and not for the")
    emit("23-patient tumour set either, where the n would support one. A test")
    emit("introduced at the moment it becomes available, on a contrast nobody")
    emit("pre-registered a test for, is the move ADR 0018 rejected. It would be")
    emit("a stop-and-ask, decided before computing it.")
    emit()
    emit("The paired and unpaired results SHARE AOIs and are not independent")
    emit("evidence. Agreement is reassurance about direction, not confirmation.")
    emit("Every brain claim states TIME-B n = 8 inline (hard constraint 8).")

    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "paired_patients": {
                    k: v for k, v in sorted(paired_by_compartment.items())
                },
                "n_paired_patients": {
                    k: len(v) for k, v in sorted(paired_by_compartment.items())
                },
                "n_cells": int(len(concordance)),
                "n_comparable": int(len(comparable)),
                "n_direction_agrees": n_agree,
                "reports_p_values": False,
                "reports_intervals": False,
                "note": (
                    "Direction of effect only, in BOTH compartments. The immune "
                    "set is n = 5 and PROJECT_PLAN §2.3 forbids a p-value; the "
                    "tumour set is n = 23 and does not get one either, because "
                    "that test was never pre-registered and adding it on "
                    "availability is the move ADR 0018 rejected. Paired and "
                    "unpaired results share AOIs and are not independent."
                ),
                "concordance": concordance.to_dict(orient="records"),
            },
            handle,
            indent=2,
            sort_keys=True,
            default=str,
        )
        handle.write("\n")

    emit()
    emit(f"wrote {len(deltas)} delta rows and {len(concordance)} concordance rows")
