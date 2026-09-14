"""P4-T3a — build the paired expression matrix, and measure the background it rests on.

Owner task: P4-T3. Driven by rule p4t3a_crosstalk_expression.

Two jobs, and the second is the one that could invalidate the first.

**1. The paired matrix.** One row per (site, side, patient, gene), carrying the
collapsed expression value each correlation in P4-T3b will use. Duplicate AOIs
collapse per ADR 0021 §1 — mean of `log2(q3 + 1)` for the primary, the
highest-detection AOI for the sensitivity — and both tables are written so the
claim that the choice is not load-bearing stays checkable.

A gene counts as detected for a patient only if it clears the floor in **every**
constituent AOI. The conservative reading, chosen because `n_detected_both` is
what tells a reader how much of a rho rests on values at background, and an
optimistic definition would understate that.

**2. The background measurement `docs/NEXT_STEPS.md` demands rather than assumes.**

ADR 0018 established that detection and expression move OPPOSITE ways under the
same background gradient. Phases 2 and 3 met that as a *group-contrast* problem.
Phase 4 does not: it correlates expression **across patients**, so the question
that matters here is a different one —

    does AOI-level background vary WITH PATIENT, shared across that patient's
    two paired AOIs?

If it does, then two genes that merely track background will correlate across
patients for no biological reason whatever, and the ranking P4-T3b produces
would be measuring AOI quality. P3-T4 measured its own version of the gradient
and found the *site* component negligible; **that is not this gradient**, and
inheriting the verdict would be assuming exactly what needs checking.

So it is measured directly: across the paired patients, within site, the
Spearman correlation between the tumour-side AOI's `negprobe_log2` and the
immune-side AOI's, and the same for `gene_detection_rate`. A high value means a
patient's two AOIs share their quality, which is the mechanism a spurious
correlation would run through — and it is the empirical case for ADR 0021 §6's
abundance-matched Null B being run at all.

**This rule reports that number. It does not adjust for it and does not gate on
it.** Adjusting would change the pre-registered specification; gating would make
a threshold out of a measurement taken to inform one. P4-T4's Null B is the
pre-registered response, and it works whatever this number turns out to be.
"""

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    adjacencies = snakemake.params.adjacencies
    pairing = snakemake.params.pairing
    gene_floor = snakemake.params.detection["detected_in_aoi_fraction"]
    background_multiple = snakemake.params.background_multiple

    # ADR 0018's threshold vocabulary, reused rather than reinvented so the two
    # phases' verdicts mean the same thing. |rho| below this is `negligible`.
    GRADIENT_NEGLIGIBLE = 0.3

    emit("P4-T3a — paired expression matrix, and the background it rests on")
    emit("=" * 70)
    emit()
    emit(f"  primary collapse      {pairing['duplicate_rule']} (ADR 0021 §1)")
    emit(f"  sensitivity collapse  {pairing['duplicate_sensitivity']}")
    emit(f"  detected for a patient: clears {background_multiple}x NegProbe-WTX")
    emit("                          in EVERY constituent AOI (conservative)")
    emit()

    # ---------------------------------------------------------------- inputs
    pairs = pd.read_csv(snakemake.input.pairs, sep="\t", float_precision="round_trip")
    pairs_s = pd.read_csv(
        snakemake.input.pairs_sensitivity, sep="\t", float_precision="round_trip"
    )
    filtered = pd.read_csv(
        snakemake.input.filtered, sep="\t", float_precision="round_trip"
    )
    exploratory = pd.read_csv(
        snakemake.input.exploratory, sep="\t", float_precision="round_trip"
    )
    qc = pd.read_csv(
        snakemake.input.qc, sep="\t", float_precision="round_trip"
    ).set_index("aoi_label")

    genes = set()
    for frame in (filtered, exploratory):
        for col in ("ligand_subunits", "receptor_subunits"):
            for cell in frame[col].astype(str):
                genes.update(s for s in cell.split(";") if s)
    genes = sorted(genes)
    emit(f"genes reached by an admitted interaction: {len(genes)}")
    emit(f"  (1:1 primary rows {len(filtered)}, complex exploratory rows "
         f"{len(exploratory)})")

    adata = ad.read_h5ad(snakemake.input.h5ad)
    missing = [g for g in genes if g not in set(map(str, adata.var_names))]
    if missing:
        raise RuntimeError(
            f"{len(missing)} admitted genes are absent from the matrix: "
            f"{missing[:10]}. p4t2d admits a gene only if it is on the panel, "
            "so this cannot happen without the two rules disagreeing."
        )

    x_log = pd.DataFrame(
        np.asarray(adata.X), index=adata.obs_names, columns=adata.var_names
    )[genes]
    q3 = pd.DataFrame(
        np.asarray(adata.layers["q3"]), index=adata.obs_names, columns=adata.var_names
    )[genes]
    negprobe = adata.obs["negprobe"]
    detected_aoi = pd.DataFrame(
        q3.to_numpy() > negprobe.to_numpy()[:, None] * background_multiple,
        index=adata.obs_names,
        columns=genes,
    )
    emit()

    # ------------------------------------------------------- collapse to patient
    def build(frame, rule, label):
        rows = []
        for r in frame.itertuples(index=False):
            for side in ("tumour", "immune"):
                aois = str(getattr(r, f"{side}_aoi_labels")).split(";")
                vals = x_log.loc[aois]
                det = detected_aoi.loc[aois]
                if rule == "mean_log2":
                    value = vals.mean(axis=0)
                elif rule == "highest_detection_rate":
                    pick = qc.loc[aois, "gene_detection_rate"].idxmax()
                    value = vals.loc[pick]
                else:
                    raise RuntimeError(
                        f"unknown collapse rule '{rule}'. The schema restricts "
                        "this to mean_log2 or highest_detection_rate (ADR 0021 "
                        "§1); a new rule is a stop-and-ask, not an edit."
                    )
                np_log2 = float(np.log2(negprobe.loc[aois]).mean())
                det_rate = float(qc.loc[aois, "gene_detection_rate"].mean())
                for gene in genes:
                    rows.append(
                        {
                            "site": r.site,
                            "side": side,
                            "patient_id": r.patient_id,
                            "aoi_code": getattr(r, f"{side}_aoi_code"),
                            "aoi_labels": ";".join(aois),
                            "n_aoi": len(aois),
                            "gene": gene,
                            "expression": float(value[gene]),
                            "detected": bool(det[gene].all()),
                            "n_aoi_detected": int(det[gene].sum()),
                            "negprobe_log2": np_log2,
                            "gene_detection_rate": det_rate,
                            "collapse_rule": rule,
                        }
                    )
        out = pd.DataFrame(rows)
        emit(f"{label}: {len(out)} rows "
             f"({out['patient_id'].nunique()} patients x 2 sides x {len(genes)} genes)")
        return out

    long = build(pairs, pairing["duplicate_rule"], "primary")
    long_s = build(pairs_s, pairing["duplicate_sensitivity"], "sensitivity")
    emit()

    # ---------------------------------------------- the background measurement
    emit("BACKGROUND — does AOI quality travel with the patient? (ADR 0018, "
         "Phase 4's form)")
    emit("-" * 70)
    emit()
    emit("  Phase 4 correlates expression ACROSS PATIENTS, so the question is")
    emit("  not whether background differs between groups — it is whether a")
    emit("  patient's two paired AOIs SHARE their background. If they do, two")
    emit("  background-tracking genes correlate for no biological reason.")
    emit("  P3-T4 measured the SITE gradient and found it negligible; that is a")
    emit("  different gradient and is not inherited.")
    emit()

    bg_rows = []
    for site in sorted(adjacencies):
        wide = (
            long[long["site"] == site]
            .drop_duplicates(["side", "patient_id"])
            .pivot(index="patient_id", columns="side")
        )
        for metric in ("negprobe_log2", "gene_detection_rate"):
            t = wide[(metric, "tumour")].to_numpy(dtype=float)
            i = wide[(metric, "immune")].to_numpy(dtype=float)
            n = int(len(t))
            rho, p = stats.spearmanr(t, i)
            verdict = "negligible" if abs(rho) < GRADIENT_NEGLIGIBLE else "shared"
            bg_rows.append(
                {
                    "site": site,
                    "metric": metric,
                    "n_patients": n,
                    "spearman_rho": float(rho),
                    "p_raw": float(p),
                    "tumour_median": float(np.median(t)),
                    "immune_median": float(np.median(i)),
                    "tumour_iqr": float(np.subtract(*np.percentile(t, [75, 25]))),
                    "immune_iqr": float(np.subtract(*np.percentile(i, [75, 25]))),
                    "verdict": verdict,
                    "threshold": GRADIENT_NEGLIGIBLE,
                }
            )
            emit(
                f"  {site:5s} {metric:20s} n={n:2d}  rho={rho:+.3f}  "
                f"p={p:.3f}   {verdict}"
            )
    background = pd.DataFrame(bg_rows)
    emit()
    shared = background.loc[background["verdict"] == "shared"]
    if len(shared):
        emit(f"  {len(shared)} of {len(background)} measurements show SHARED "
             "AOI quality across a patient's paired AOIs.")
        emit("  This is the empirical case for ADR 0021 §6's abundance-matched")
        emit("  Null B, and it is why that control was pre-registered rather")
        emit("  than added after seeing the ranking.")
    else:
        emit("  No measurement reaches the threshold. Null B still runs — it is")
        emit("  pre-registered (ADR 0021 §6) and a control is not contingent on")
        emit("  the confound it guards against turning out to be large.")
    emit()
    emit("  THIS RULE REPORTS THE NUMBER. It does not adjust for it and does")
    emit("  not gate on it: adjusting would change the pre-registered")
    emit("  specification, and gating would make a threshold out of a")
    emit("  measurement taken to inform one.")
    emit()

    # ------------------------------------------------------------ assertions
    n_checks = 0

    for frame, label in ((long, "primary"), (long_s, "sensitivity")):
        if frame.duplicated(["site", "side", "patient_id", "gene"]).any():
            raise RuntimeError(f"{label}: duplicate (site, side, patient, gene) row.")
        n_checks += 1
        if not np.isfinite(frame["expression"]).all():
            raise RuntimeError(f"{label}: non-finite expression value.")
        n_checks += 1

    for site in sorted(adjacencies):
        want = int(snakemake.params.pairing["expected_n_patients"][site])
        got = int(long.loc[long["site"] == site, "patient_id"].nunique())
        if got != want:
            raise RuntimeError(
                f"{site}: {got} patients in the expression matrix, expected "
                f"{want}. p4t1_build_pairs already asserts this; a disagreement "
                "means the two rules read different pairings."
            )
    n_checks += 1

    if (long_s["n_aoi"] != 1).any():
        raise RuntimeError(
            "sensitivity: a row collapses more than one AOI. "
            "highest_detection_rate must resolve to exactly one."
        )
    n_checks += 1

    # A single-AOI patient must give the SAME value under both rules — the mean
    # of one number is that number. Anything else is an indexing bug, and it
    # would silently make the sensitivity look like a real disagreement.
    single = long[long["n_aoi"] == 1].set_index(["site", "side", "patient_id", "gene"])
    other = long_s.set_index(["site", "side", "patient_id", "gene"])
    common = single.index.intersection(other.index)
    delta = float(
        np.abs(single.loc[common, "expression"] - other.loc[common, "expression"]).max()
    )
    if delta > 1e-12:
        raise RuntimeError(
            f"single-AOI patients differ between the two collapse rules by "
            f"{delta:.3g}. With one AOI the two rules are the same operation, "
            "so any difference is an indexing bug — and it would make the "
            "sensitivity look like a disagreement about the science."
        )
    n_checks += 1
    emit(f"single-AOI agreement between the two collapse rules: max |delta| = "
         f"{delta:.3g} over {len(common)} rows")

    emit(f"assertions: {n_checks} of {n_checks} passed")
    emit()

    long.to_csv(snakemake.output.expression, sep="\t", index=False)
    long_s.to_csv(snakemake.output.sensitivity, sep="\t", index=False)
    background.to_csv(snakemake.output.background, sep="\t", index=False)

    summary = {
        "task": "P4-T3a",
        "aim": "A5",
        "exploratory": True,
        "exploratory_adr": "ADR 0014",
        "pre_registered": "ADR 0021 §1",
        "collapse_rule": pairing["duplicate_rule"],
        "collapse_sensitivity": pairing["duplicate_sensitivity"],
        "detected_definition": (
            "clears qc.detection_background_multiple x that AOI's NegProbe-WTX "
            "in EVERY constituent AOI — the conservative reading, because "
            "n_detected_both is what tells a reader how much of a rho rests on "
            "values at background"
        ),
        "n_genes": len(genes),
        "n_rows": len(long),
        "n_patients": {
            site: int(long.loc[long["site"] == site, "patient_id"].nunique())
            for site in sorted(adjacencies)
        },
        "background_check": background.to_dict("records"),
        "background_threshold": GRADIENT_NEGLIGIBLE,
        "n_shared_background_measurements": int(len(shared)),
        "background_rule": (
            "Reported, never adjusted for and never gated on. Phase 4 "
            "correlates across patients, so the relevant gradient is whether a "
            "patient's two paired AOIs share their quality — not the site "
            "gradient P3-T4 measured. ADR 0021 §6's Null B is the "
            "pre-registered response and runs regardless."
        ),
        "single_aoi_max_delta": delta,
        "n_assertions": n_checks,
        "h5ad": str(snakemake.input.h5ad),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit(f"wrote {snakemake.output.expression}, {snakemake.output.sensitivity}, "
         f"{snakemake.output.background}, {snakemake.output.summary}")
