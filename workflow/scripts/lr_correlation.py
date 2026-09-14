"""P4-T3b — Spearman correlation of ligand against paired receptor, across patients.

Owner task: P4-T3. Driven by rule p4t3b_lr_correlation.

For each admitted interaction, in each direction, within each site: the ligand's
collapsed expression in the compartment it is measured in, against the
receptor's in the paired compartment, **across patients**. Ranked by |rho|.

**Spearman, never Pearson.** PROJECT_PLAN §6 P4-T3 is explicit and the reason is
the n: at 13 lung patients and **8 brain** ones, a Pearson coefficient is a coin
flip on a single outlier. `crosstalk.correlation.method` is an enum of one in
the schema so a config edit cannot quietly change it.

**Both directions, labelled.** CellChatDB is directed, and "ligand on the tumour
AOI, receptor on the paired immune AOI" is a different biological claim from its
reverse (ADR 0021 §9). Neither is a control for the other; both are reported.

**The primary keeps every paired patient** — lung n = 13, brain n = **8** — with
`n_detected_both` on every row (ADR 0021 §5). Dropping patients whose values sit
at background is expression-dependent selection: it removes the low values of
*both* partners preferentially, which is a mechanism for inducing a correlation
rather than removing one, and it would make n vary per pair so P4-T4's null
would need recalibrating at every distinct n. The detected-both refit runs as a
**declared sensitivity**, bound by `crosstalk.correlation.min_n_patients`, with
a `delta_vs_primary` column — ADR 0012's shape.

**Two tables, and only one of them will ever carry an FDR** (ADR 0021 §4):

  primary      one-to-one interactions. P4-T4 attaches the empirical FDR here.
  exploratory  complex interactions, every subunit above the floor. It carries
               rho, raw p and n and **NO q of any kind, ever** — the same shape
               `fit_checkpoint_models.R` gives its exploratory table, which
               keeps `p_raw` and sets `q_bh <- NULL`. No Phase 4 sentence may
               rest on it.

**NOTHING HERE IS INTERPRETABLE ON ITS OWN.** With roughly a thousand
interactions at n = 13, |rho| > 0.7 arises by chance, and PROJECT_PLAN §6 says
so in terms: "without this, the ranking is uninterpretable". The ranking this
rule produces is an input to P4-T4, not a result. Gate 4 turns on P4-T4.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    method = snakemake.params.correlation["method"]
    min_n = snakemake.params.correlation["min_n_patients"]
    adjacencies = snakemake.params.adjacencies

    if method != "spearman":
        raise RuntimeError(
            f"crosstalk.correlation.method is '{method}', but this script "
            "computes Spearman and PROJECT_PLAN §6 P4-T3 requires it — at "
            "n = 13 and n = 8 Pearson is a coin flip on one outlier. The schema "
            "restricts this field to an enum of one; changing it is a "
            "stop-and-ask, not a config edit."
        )

    emit("P4-T3b — ligand vs. paired receptor, across patients, within site")
    emit("=" * 70)
    emit()
    emit(f"  method         {method} (never Pearson — PROJECT_PLAN §6 P4-T3)")
    emit("  primary        every paired patient, n_detected_both as a column")
    emit(f"  sensitivity    detected-both only, min n = {min_n}")
    emit("  directions     both, labelled (ADR 0021 §9)")
    emit()
    emit("  NOTHING HERE IS INTERPRETABLE ON ITS OWN. With ~1000 interactions")
    emit("  at n = 13, |rho| > 0.7 arises by chance. This ranking is an INPUT")
    emit("  to P4-T4's empirical FDR, which is what Gate 4 turns on.")
    emit()

    # ---------------------------------------------------------------- inputs
    long = pd.read_csv(
        snakemake.input.expression, sep="\t", float_precision="round_trip"
    )
    filtered = pd.read_csv(
        snakemake.input.filtered, sep="\t", float_precision="round_trip"
    )
    exploratory = pd.read_csv(
        snakemake.input.exploratory, sep="\t", float_precision="round_trip"
    )

    # (site, side, gene) -> Series indexed by patient, in one fixed patient
    # order per site so the ligand and receptor vectors are aligned by
    # construction rather than by a merge nobody checks.
    patients = {
        site: sorted(
            long.loc[long["site"] == site, "patient_id"].unique(),
            key=lambda p: int(str(p).lstrip("P")) if str(p)[1:].isdigit() else 0,
        )
        for site in sorted(adjacencies)
    }
    for site, ps in patients.items():
        emit(f"  {site:5s} n = {len(ps)}  {ps}")
    emit()

    value = {}
    detected = {}
    for (site, side, gene), grp in long.groupby(["site", "side", "gene"]):
        g = grp.set_index("patient_id")
        value[(site, side, gene)] = g["expression"].reindex(patients[site]).to_numpy()
        detected[(site, side, gene)] = (
            g["detected"].reindex(patients[site]).to_numpy().astype(bool)
        )

    side_of = {}
    for site in sorted(adjacencies):
        side_of[(site, adjacencies[site]["tumour"])] = "tumour"
        side_of[(site, adjacencies[site]["immune"])] = "immune"

    def vector(site, aoi_code, subunits):
        """Collapsed value for one partner. A complex is the MEAN of its
        subunits on the log2 scale — arithmetic mean of log2 is the geometric
        mean on the linear scale, which is CellChat's own convention for
        scoring a heteromer."""
        side = side_of[(site, aoi_code)]
        vals = np.vstack([value[(site, side, g)] for g in subunits])
        det = np.vstack([detected[(site, side, g)] for g in subunits])
        return vals.mean(axis=0), det.all(axis=0)

    # ------------------------------------------------------------- correlate
    def correlate(frame, label):
        rows = []
        for r in frame.itertuples(index=False):
            lig = [s for s in str(r.ligand_subunits).split(";") if s]
            rec = [s for s in str(r.receptor_subunits).split(";") if s]
            lv, ld = vector(r.site, r.ligand_aoi_code, lig)
            rv, rd = vector(r.site, r.receptor_aoi_code, rec)
            both = ld & rd
            n = int(len(lv))

            rho, p = stats.spearmanr(lv, rv)

            # Sensitivity: detected-both patients only, at or above min_n.
            n_both = int(both.sum())
            if n_both >= min_n:
                rho_s, p_s = stats.spearmanr(lv[both], rv[both])
                delta = float(rho_s - rho)
                status_s = "fitted"
            else:
                rho_s, p_s, delta = float("nan"), float("nan"), float("nan")
                status_s = f"below_min_n ({n_both} < {min_n})"

            rows.append(
                {
                    "interaction_name": r.interaction_name,
                    "interaction_name_2": r.interaction_name_2,
                    "pathway_name": r.pathway_name,
                    "annotation": r.annotation,
                    "site": r.site,
                    "direction": r.direction,
                    "ligand": r.ligand,
                    "receptor": r.receptor,
                    "ligand_subunits": r.ligand_subunits,
                    "receptor_subunits": r.receptor_subunits,
                    "ligand_aoi_code": r.ligand_aoi_code,
                    "receptor_aoi_code": r.receptor_aoi_code,
                    "table": label,
                    "n": n,
                    "rho": float(rho),
                    "abs_rho": float(abs(rho)),
                    "p_raw": float(p),
                    "n_detected_both": n_both,
                    "frac_detected_both": float(n_both / n) if n else float("nan"),
                    "rho_detected_only": float(rho_s),
                    "p_detected_only": float(p_s),
                    "delta_vs_primary": delta,
                    "sensitivity_status": status_s,
                    # ADR 0014's binding term, carried on every row so a rho can
                    # never be read without it.
                    "ligand_detection": r.ligand_detection,
                    "receptor_detection": r.receptor_detection,
                    "ligand_secreted_type": r.ligand_secreted_type,
                    "ligand_transmembrane": r.ligand_transmembrane,
                    "receptor_transmembrane": r.receptor_transmembrane,
                    "evidence": r.evidence,
                }
            )
        out = pd.DataFrame(rows)
        if len(out):
            out = out.sort_values("abs_rho", ascending=False).reset_index(drop=True)
            out.insert(0, "rank", np.arange(1, len(out) + 1))
        emit(f"{label}: {len(out)} direction-rows correlated")
        return out

    primary = correlate(filtered, "primary")
    explore = correlate(exploratory, "exploratory")
    emit()

    # No FDR reaches either table here — P4-T4 attaches it, and only to the
    # primary. Asserted rather than trusted, the way checkpoint_paired_check.py
    # asserts no p-value can reach its tables.
    forbidden = {"q_bh", "q_value", "empirical_fdr", "fdr"}
    for frame, label in ((primary, "primary"), (explore, "exploratory")):
        leaked = forbidden & set(frame.columns)
        if leaked:
            raise RuntimeError(
                f"{label} carries {sorted(leaked)}. The empirical FDR is "
                "P4-T4's and attaches to the primary alone; the exploratory "
                "table gets none of any kind, ever (ADR 0021 §4)."
            )

    # ------------------------------------------------------------- reporting
    for site in sorted(adjacencies):
        n_site = len(patients[site])
        emit(f"{site} (n = {n_site}) — top 10 by |rho|, PRIMARY, "
             "UNCALIBRATED AND NOT A RESULT")
        sel = primary[primary["site"] == site].head(10)
        emit(f"  {'rho':>7s} {'p':>8s} {'n':>3s} {'det':>6s}  "
             f"{'class':<20s} interaction / direction")
        for r in sel.itertuples(index=False):
            arrow = "T->I" if r.direction.startswith("ligand_tumour") else "I->T"
            emit(
                f"  {r.rho:+7.3f} {r.p_raw:8.4f} {r.n:3d} "
                f"{r.n_detected_both:3d}/{r.n:<2d} {r.annotation:<20s} "
                f"{r.interaction_name_2}  [{arrow}]"
            )
        emit()

    emit("  The `det` column is n_detected_both/n. A rho resting on few")
    emit("  detected patients is a rho computed largely on values at")
    emit("  background — which is why it is a column and not a filter")
    emit("  (ADR 0021 §5).")
    emit()

    n_disagree = int(
        (
            (primary["sensitivity_status"] == "fitted")
            & (np.sign(primary["rho"]) != np.sign(primary["rho_detected_only"]))
        ).sum()
    )
    n_fitted = int((primary["sensitivity_status"] == "fitted").sum())
    emit(f"sensitivity (detected-both only, min n = {min_n}):")
    emit(f"  fitted        {n_fitted} of {len(primary)} primary rows")
    emit(f"  sign flips    {n_disagree}")
    if n_fitted:
        med = float(np.nanmedian(np.abs(primary['delta_vs_primary'])))
        emit(f"  median |delta_vs_primary|  {med:.3f}")
    emit()

    # ------------------------------------------------------------ assertions
    n_checks = 0
    for frame, label in ((primary, "primary"), (explore, "exploratory")):
        if len(frame) and frame.duplicated(
            ["interaction_name", "site", "direction"]
        ).any():
            raise RuntimeError(f"{label}: duplicate interaction/site/direction row.")
        n_checks += 1
        if len(frame) and not np.isfinite(frame["rho"]).all():
            bad = frame.loc[~np.isfinite(frame["rho"]), "interaction_name"].tolist()
            raise RuntimeError(
                f"{label}: non-finite rho for {bad[:10]}. Spearman returns NaN "
                "when a vector is constant; an admitted gene above the floor in "
                "half a compartment's AOIs should not be constant across "
                "patients, so this is a data problem, not a rounding one."
            )
        n_checks += 1

    for site in sorted(adjacencies):
        want = len(patients[site])
        got = set(primary.loc[primary["site"] == site, "n"])
        if got and got != {want}:
            raise RuntimeError(
                f"{site}: primary rows have n in {sorted(got)}, expected all "
                f"{want}. The primary keeps EVERY paired patient so the "
                "permutation null preserves both n's exactly (ADR 0021 §5); a "
                "varying n here means patients were dropped somewhere."
            )
    n_checks += 1
    emit(f"assertions: {n_checks} of {n_checks} passed "
         "(no duplicate rows; rho finite; n fixed per site)")
    emit()

    primary.to_csv(snakemake.output.primary, sep="\t", index=False)
    explore.to_csv(snakemake.output.exploratory, sep="\t", index=False)

    summary = {
        "task": "P4-T3b",
        "aim": "A5",
        "exploratory": True,
        "exploratory_adr": "ADR 0014",
        "pre_registered": "ADR 0021 §§5, 9",
        "method": method,
        "complex_value_rule": (
            "mean of subunit log2 values — the arithmetic mean of log2 is the "
            "geometric mean on the linear scale, CellChat's own convention for "
            "scoring a heteromer"
        ),
        "min_n_patients_sensitivity": min_n,
        "n_primary_rows": len(primary),
        "n_exploratory_rows": len(explore),
        "n_patients": {site: len(patients[site]) for site in sorted(adjacencies)},
        "patients": {site: patients[site] for site in sorted(adjacencies)},
        "n_sensitivity_fitted": n_fitted,
        "n_sensitivity_sign_flips": n_disagree,
        "n_time_b": len(patients.get("brain", [])),
        "power_floor_sd": "1.1-1.3",
        "interpretation_rule": (
            "UNCALIBRATED. With ~1000 interactions at n = 13, |rho| > 0.7 "
            "arises by chance; PROJECT_PLAN §6 P4-T4 says the ranking is "
            "uninterpretable without an empirical null. This table is an INPUT "
            "to P4-T4, not a result, and Gate 4 turns on P4-T4."
        ),
        "reporting_rule": (
            "Every row carries n, n_detected_both and both partners' per-site "
            "detection counts (ADR 0014). TIME-B n = 8 (hard constraint 8). "
            "The exploratory table gets no q of any kind, ever (ADR 0021 §4)."
        ),
        "n_assertions": n_checks,
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit(f"wrote {snakemake.output.primary}, {snakemake.output.exploratory}, "
         f"{snakemake.output.summary}")
