"""P2-T4 — the within-patient direction check, on the 5 paired patients.

Owner task: P2-T4. Driven by rule p2t4_paired_check.

Five patients contribute both a TIME-L and a TIME-B AOI. Restricted to them,
the lung-vs-brain difference is measured *within* a patient, so it is free of
between-patient confounding by construction — which is worth having, and is the
only thing this task is for.

**DIRECTION OF EFFECT ONLY. NO P-VALUES FROM n = 5.** PROJECT_PLAN §6 and
CLAUDE.md both say so, and §2.3 calls this "a consistency check, never a
headline result". This script therefore does not compute a test statistic at
all: there is no p-value column to be tempted by, and no significance to read
into a sign count. Five paired observations agreeing in direction is weak
evidence pointing the same way as the unpaired fit; it is not confirmation of
it, and the two are not independent — the same AOIs are in both.

A patient contributing more than one AOI at a site has those AOIs averaged
before differencing, so each patient contributes exactly one delta and no
patient is weighted by how many AOIs they happen to have.
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

    scores = pd.read_csv(snakemake.input.scores, sep="\t")
    models = pd.read_csv(snakemake.input.models, sep="\t")

    # --------------------------------------------------------- paired patients
    per_patient_sites = scores.groupby("patient_id")["site"].nunique()
    paired = sorted(per_patient_sites[per_patient_sites == 2].index)
    emit(f"patients with TIME AOIs at both sites: {len(paired)} -> {paired}")

    # §2.3 says five. If this moves, the subset or the design changed.
    if len(paired) != 5:
        raise RuntimeError(
            f"expected 5 paired patients (PROJECT_PLAN §2.3), found {len(paired)}: "
            f"{paired}. Do not adjust the expectation to match the data."
        )

    d = scores.loc[scores["patient_id"].isin(paired)].copy()

    # One value per patient x site x signature x method: a patient with two
    # AOIs at a site must not count twice.
    per_site = (
        d.groupby(["signature", "method", "patient_id", "site"], observed=True)["score"]
        .mean()
        .unstack("site")
    )
    n_collapsed = len(d) - len(per_site) * 2
    emit(f"AOI rows collapsed to patient means: {n_collapsed}")

    per_site["delta_brain_minus_lung"] = per_site["brain"] - per_site["lung"]
    deltas = per_site.reset_index()
    deltas.to_csv(snakemake.output.deltas, sep="\t", index=False, float_format=FMT)
    emit()

    # ------------------------------------------------------------- concordance
    unpaired = models.loc[models["model"] == "primary"].set_index(
        ["signature", "method"]
    )["estimate"]

    rows = []
    for (sig, meth), g in deltas.groupby(["signature", "method"], observed=True):
        delta = g["delta_brain_minus_lung"].to_numpy()
        mean_delta = float(np.mean(delta))
        est = float(unpaired.loc[(sig, meth)])
        # np.sign(0.0) is 0.0, which would silently read as "agrees" against a
        # zero estimate; require both to be non-zero and matching.
        agrees = bool(np.sign(mean_delta) == np.sign(est) and mean_delta != 0 and est != 0)
        rows.append(
            {
                "signature": sig,
                "method": meth,
                "n_patients": int(len(delta)),
                "n_brain_higher": int((delta > 0).sum()),
                "n_lung_higher": int((delta < 0).sum()),
                "mean_delta_brain_minus_lung": mean_delta,
                "median_delta_brain_minus_lung": float(np.median(delta)),
                "paired_direction": "brain_higher" if mean_delta > 0 else "lung_higher",
                "unpaired_estimate": est,
                "unpaired_direction": "brain_higher" if est > 0 else "lung_higher",
                "direction_agrees": agrees,
            }
        )

    concordance = pd.DataFrame(rows).sort_values(["method", "signature"])

    forbidden = {"p", "p_raw", "p_value", "pvalue", "q", "q_bh", "statistic"}
    assert not (forbidden & set(concordance.columns)), (
        "a test statistic reached the paired table; n = 5 does not support one"
    )

    concordance.to_csv(
        snakemake.output.concordance, sep="\t", index=False, float_format=FMT
    )

    emit("within-patient direction vs. the unpaired primary fit (n = 5 patients):")
    for meth in sorted(concordance["method"].unique()):
        emit("")
        emit(f"  method = {meth}")
        for _, r in concordance.loc[concordance["method"] == meth].iterrows():
            emit(
                f"    {r['signature']:22s} "
                f"paired {r['mean_delta_brain_minus_lung']:+.3f} "
                f"({r['n_brain_higher']}/{r['n_patients']} brain-higher)  "
                f"unpaired {r['unpaired_estimate']:+.3f}  "
                f"{'agrees' if r['direction_agrees'] else 'DISAGREES'}"
            )

    n_agree = int(concordance["direction_agrees"].sum())
    emit()
    emit(f"direction agrees in {n_agree} of {len(concordance)} signature x method cells")
    emit()
    emit("NO P-VALUE IS REPORTED HERE AND NONE SHOULD BE ADDED. n = 5 paired")
    emit("patients is a consistency check on direction (PROJECT_PLAN §2.3), not")
    emit("a test. The paired and unpaired results share AOIs and are not")
    emit("independent evidence.")

    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "paired_patients": paired,
                "n_paired_patients": len(paired),
                "n_cells": int(len(concordance)),
                "n_direction_agrees": n_agree,
                "reports_p_values": False,
                "note": (
                    "Direction of effect only. n = 5; paired and unpaired "
                    "results share AOIs and are not independent."
                ),
                "concordance": concordance.to_dict(orient="records"),
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")
