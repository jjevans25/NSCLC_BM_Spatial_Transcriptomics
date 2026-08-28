"""P2-T6 — do deconvolution and signature scoring agree on direction?

Owner task: P2-T6. Driven by rule p2t6_convergence_check.

PROJECT_PLAN §6 asks whether the two methods agree and, where they disagree,
which is trusted. Answering that honestly first requires saying **for which
signatures the question is even askable**, and for most of them it is not. Three
independent things can disqualify a signature, and this script tests all three
rather than letting a disagreement be reported as if it were about biology:

  1. **No counterpart.** `antigen_presentation` is a transcriptional program
     expressed by tumour cells and professional APCs alike, not a lineage
     safeTME resolves. There is no abundance to compare a program against, and
     substituting "B + mDCs + macrophages" would swap the quantity silently.

  2. **Not assessable.** A signature below its detection-coverage floor at
     either site cannot be compared with anything there (ADR 0008). Its score
     is not a measurement of the biology, so agreement or disagreement with
     deconvolution says nothing.

  3. **Degenerate mapping.** `myeloid_m1` and `myeloid_m2` name the *same*
     cell type, because M1/M2 is a polarisation axis *within* a lineage and
     deconvolution resolves the lineage only. Macrophage abundance moving up
     cannot confirm or refute a polarisation score, so a match here would be
     meaningless and a mismatch would be uninformative.

Only the safeTME arm participates. The two-matrix arm uses a different
reference per site — and Brain_Darmanis resolves **no lymphoid type at all** —
so a cross-site direction from it would be the reference rather than the
biology (P2-T5).

Deconvolution direction is descriptive: a difference in mean proportion, with
no test attached. The statistical contrast of record stays the mixed model.
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

    signatures = snakemake.params.signatures
    aoi_codes = list(snakemake.params.aoi_codes)
    lung_code, brain_code = aoi_codes

    models = pd.read_csv(snakemake.input.models, sep="\t")
    coverage = pd.read_csv(snakemake.input.coverage, sep="\t")
    comp = pd.read_csv(snakemake.input.composition, sep="\t")

    comp = comp.loc[comp["arm"] == "safetme"]
    emit(f"deconvolution arm used: safetme ({comp['cell_type'].nunique()} cell types)")
    emit("the tissue arm is excluded: different reference per site, and")
    emit("Brain_Darmanis resolves no lymphoid type (P2-T5).")
    emit()

    # Which mapped sets are shared by more than one signature?
    mapped = {
        name: tuple(sorted(spec["decon_cell_types"]))
        for name, spec in signatures.items()
        if spec.get("decon_cell_types")
    }
    shared = {}
    for name, types in mapped.items():
        shared.setdefault(types, []).append(name)

    known_types = set(comp["cell_type"])
    rows = []
    for name in sorted(signatures):
        types = mapped.get(name)
        has_counterpart = types is not None

        cov = coverage.loc[coverage["signature"] == name]
        assessable = bool(cov["assessable"].all())
        bad_sites = sorted(cov.loc[~cov["assessable"], "aoi_code"])

        peers = [p for p in shared.get(types, []) if p != name] if has_counterpart else []
        unique_mapping = has_counterpart and not peers

        # Signature direction, primary model. Both methods agree in sign
        # throughout, but the disagreement is recorded rather than assumed.
        prim = models.loc[(models["signature"] == name) & (models["model"] == "primary")]
        sig_signs = {r["method"]: float(np.sign(r["estimate"])) for _, r in prim.iterrows()}
        sig_sign = sig_signs.get("zscore", np.nan)
        methods_agree = len(set(sig_signs.values())) == 1

        # Deconvolution direction: summed proportion of the mapped types.
        decon_delta = np.nan
        decon_sign = np.nan
        lung_p = brain_p = np.nan
        if has_counterpart:
            missing = [t for t in types if t not in known_types]
            if missing:
                raise RuntimeError(
                    f"{name}: decon_cell_types {missing} are not produced by the "
                    f"safeTME arm. Valid types: {sorted(known_types)}"
                )
            sel = comp.loc[comp["cell_type"].isin(types)]
            per_aoi = sel.groupby(["aoi_label", "aoi_code"], observed=True)["proportion"].sum()
            per_aoi = per_aoi.reset_index()
            lung_p = float(per_aoi.loc[per_aoi["aoi_code"] == lung_code, "proportion"].mean())
            brain_p = float(per_aoi.loc[per_aoi["aoi_code"] == brain_code, "proportion"].mean())
            decon_delta = brain_p - lung_p
            decon_sign = float(np.sign(decon_delta))

        checkable = bool(has_counterpart and assessable and unique_mapping)
        agrees = bool(checkable and sig_sign == decon_sign and sig_sign != 0)

        if not has_counterpart:
            verdict = "no_counterpart"
        elif not assessable:
            verdict = f"not_assessable_in_{','.join(bad_sites)}"
        elif not unique_mapping:
            verdict = f"degenerate_mapping_with_{','.join(sorted(peers))}"
        else:
            verdict = "checkable"

        rows.append(
            {
                "signature": name,
                "decon_cell_types": ",".join(types) if has_counterpart else "",
                "has_counterpart": has_counterpart,
                "assessable_both_sites": assessable,
                "unique_mapping": unique_mapping,
                "checkable": checkable,
                "verdict": verdict,
                "signature_estimate_zscore": float(
                    prim.loc[prim["method"] == "zscore", "estimate"].iloc[0]
                ),
                "signature_direction": "brain_higher" if sig_sign > 0 else "lung_higher",
                "methods_agree_on_sign": methods_agree,
                "decon_prop_lung": lung_p,
                "decon_prop_brain": brain_p,
                "decon_delta_brain_minus_lung": decon_delta,
                "decon_direction": (
                    "" if not has_counterpart
                    else ("brain_higher" if decon_sign > 0 else "lung_higher")
                ),
                "direction_agrees": agrees if checkable else None,
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(snakemake.output.table, sep="\t", index=False, float_format=FMT)

    emit("per-signature convergence status:")
    for _, r in out.iterrows():
        emit(f"  {r['signature']:22s} {r['verdict']}")
        if r["has_counterpart"]:
            emit(
                f"      signature {r['signature_direction']:12s} "
                f"({r['signature_estimate_zscore']:+.3f} SD)   "
                f"decon {r['decon_direction']:12s} "
                f"({r['decon_prop_lung']:.3f} -> {r['decon_prop_brain']:.3f}, "
                f"delta {r['decon_delta_brain_minus_lung']:+.3f})"
            )
    emit()

    checkables = out.loc[out["checkable"], "signature"].tolist()
    emit(f"genuinely checkable: {len(checkables)} of {len(out)} -> {checkables}")

    if not checkables:
        raise RuntimeError(
            "no signature is checkable, so P2-T6 has nothing to converge. "
            "That is a real possible outcome but it must be looked at, not "
            "written past."
        )

    n_agree = int(out.loc[out["checkable"], "direction_agrees"].sum())
    emit(f"of those, direction agrees in {n_agree}")
    emit()
    emit("Deconvolution proportions are DESCRIPTIVE — no test is attached, and")
    emit("the contrast of record remains the mixed model in signature_models.tsv.")

    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "arm_used": "safetme",
                "n_time_b": 8,
                "n_signatures": int(len(out)),
                "n_checkable": len(checkables),
                "checkable": checkables,
                "n_direction_agrees": n_agree,
                "excluded": {
                    r["signature"]: r["verdict"]
                    for _, r in out.loc[~out["checkable"]].iterrows()
                },
                "trusted_where_they_disagree": (
                    "signatures — they do not depend on a reference matrix, "
                    "whereas deconvolution's answer is conditional on the "
                    "profile matrix chosen"
                ),
                "rows": out.to_dict(orient="records"),
            },
            handle,
            indent=2,
            sort_keys=True,
            default=str,
        )
        handle.write("\n")
