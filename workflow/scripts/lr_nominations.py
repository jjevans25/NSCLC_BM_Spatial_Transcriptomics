"""P4-T5 — nominate 5-10 inferred crosstalk pairs, as hypotheses and nothing more.

Owner task: P4-T5. Driven by rule p4t5_nominations.

Top-ranked, detection-filtered, surviving the empirical FDR, framed explicitly
as hypothesis generation. Each row carries rho, the empirical FDR, the
abundance-matched p, n, `n_detected_both`, and **the per-site detection counts
of both partners** — which is ADR 0014's binding term, not a nicety: a
nomination without them is not reportable.

**IF NOTHING SURVIVES, THAT IS THE DELIVERABLE.** `docs/NEXT_STEPS.md` is
explicit — "write it as such and do not go looking for a softer threshold" — and
Gate 4 licenses it in terms: "no LR pair exceeded chance expectation at n=13" is
an honest, useful result. This script therefore has no fallback path. It reads
the FDR that P4-T4 computed, applies the alpha that was pre-registered, and
writes whatever that leaves, including nothing. There is no second threshold in
it to reach for.

**Gate 4 does not turn on this table.** It turns on P4-T4 — on the empirical FDR
being computed and honoured. A surviving list and an empty list both pass. This
rule must not drift into the gate's subject, and the cap below is the mechanism:
`max_nominations` truncates a long list rather than licensing a longer one.

**The rationale is looked up, never generated.** A one-line biological rationale
per nomination comes from the database's own pathway annotation plus the
interaction class — `pathway_name`, `annotation`, and the evidence string
CellChatDB ships. Writing prose about why a correlation makes biological sense,
after seeing that it correlated, is how a ranked list becomes a story.

**The anticipated null is an assay-sensitivity limit** (ADR 0014), never
evidence that the crosstalk is absent. A pair absent from this table is a pair
that either could not be measured or did not beat chance — and those two are
different, so the summary counts them separately.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    alpha = snakemake.params.null_calibration["fdr_alpha"]
    max_nom = snakemake.params.max_nominations
    adjacencies = snakemake.params.adjacencies

    emit("P4-T5 — nominations, as hypotheses and nothing more")
    emit("=" * 70)
    emit()
    emit(f"  empirical FDR alpha   {alpha}  (pre-registered, ADR 0021 §6)")
    emit(f"  cap per site          {max_nom}")
    emit()
    emit("  IF NOTHING SURVIVES, THAT IS THE DELIVERABLE. There is no second")
    emit("  threshold in this script to reach for, and Gate 4 does not turn on")
    emit("  this table — it turns on P4-T4's FDR being computed and honoured.")
    emit()

    fdr = pd.read_csv(snakemake.input.fdr, sep="\t", float_precision="round_trip")
    membership = pd.read_csv(
        snakemake.input.membership, sep="\t", float_precision="round_trip"
    )
    emit(f"primary direction-rows with an FDR: {len(fdr)}")

    # Detection counts are joined from the membership table rather than trusted
    # from the correlation table's copy: ADR 0014's term is that every
    # nomination STATES the detection of both partners, and a joined value that
    # disagrees with its source is worse than no value at all.
    key = ["interaction_name", "site", "direction"]
    det = membership.set_index(key)[["ligand_detection", "receptor_detection"]]
    joined = fdr.set_index(key)
    mismatch = 0
    for col in ("ligand_detection", "receptor_detection"):
        a = joined[col].astype(str)
        b = det[col].reindex(joined.index).astype(str)
        mismatch += int((a != b).sum())
    if mismatch:
        raise RuntimeError(
            f"{mismatch} detection strings disagree between "
            "crosstalk_fdr.tsv and crosstalk_lr_membership.tsv. ADR 0014 makes "
            "the per-site detection of both partners part of every nomination; "
            "two sources that disagree mean one of them is stale."
        )
    emit("detection strings agree with crosstalk_lr_membership.tsv on all rows")
    emit()

    survivors = fdr[fdr["clears_empirical_fdr"]].copy()
    emit(f"clears the empirical FDR at {alpha}: {len(survivors)} of {len(fdr)}")

    rows = []
    for site in sorted(adjacencies):
        s = (
            survivors[survivors["site"] == site]
            .sort_values("abs_rho", ascending=False)
            .head(max_nom)
        )
        emit(f"  {site:5s} {int((survivors['site'] == site).sum()):4d} survivors "
             f"-> {len(s)} nominated (cap {max_nom})")
        rows.append(s)
    nominations = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    emit()

    if len(nominations):
        # Rationale, LOOKED UP, never generated. The database's own pathway and
        # class annotation plus its evidence string — nothing written after
        # seeing which pairs correlated.
        def rationale(r):
            direction = (
                f"{r.ligand} in the {r.ligand_aoi_code} compartment against "
                f"{r.receptor} in the paired {r.receptor_aoi_code}"
            )
            return (
                f"CellChatDB {r.annotation.lower()} interaction in the "
                f"{r.pathway_name} pathway ({r.evidence}); inferred from "
                f"{direction} across {r.n} patients."
            )

        nominations["rationale"] = [rationale(r) for r in nominations.itertuples()]
        nominations["n_time_b"] = np.where(nominations["site"] == "brain", 8, np.nan)
        nominations["framing"] = "hypothesis generation, not a finding"

        emit("NOMINATIONS — hypothesis generation, NOT findings")
        emit("-" * 70)
        for site in sorted(adjacencies):
            sel = nominations[nominations["site"] == site]
            if not len(sel):
                continue
            n_pat = int(sel["n"].iloc[0])
            tag = "  TIME-B n = 8" if site == "brain" else ""
            emit(f"{site} (n = {n_pat} patients){tag}")
            for r in sel.itertuples(index=False):
                arrow = (
                    "tumour->immune"
                    if r.direction.startswith("ligand_tumour")
                    else "immune->tumour"
                )
                emit(f"  {r.interaction_name_2}  [{arrow}]  {r.annotation}")
                emit(f"    rho {r.rho:+.3f}   empirical FDR {r.empirical_fdr:.4f}   "
                     f"abundance-matched p {r.abundance_matched_p:.4f}")
                emit(f"    n = {r.n}, detected in both {r.n_detected_both}/{r.n}")
                emit(f"    ligand   {r.ligand_detection}")
                emit(f"    receptor {r.receptor_detection}")
                if not r.clears_abundance_matched:
                    emit("    NOTE: clears the permutation null but NOT the")
                    emit("          abundance-matched control — consistent with")
                    emit("          'high-expression genes correlate with")
                    emit("          high-expression genes' (ADR 0021 §6).")
                emit()
    else:
        emit("NO NOMINATION SURVIVES THE EMPIRICAL FDR.")
        emit("-" * 70)
        emit()
        emit("  THIS IS THE DELIVERABLE, not a failure of the analysis.")
        emit("  Gate 4: \"If nothing survives, report that: 'no LR pair exceeded")
        emit("  chance expectation at n=13' is an honest, useful, publishable-")
        emit("  to-blog result, and it's a better outcome than a ranked list")
        emit("  you can't defend.\"")
        emit()
        emit("  And read it the way ADR 0014 would: this is an ASSAY-SENSITIVITY")
        emit("  LIMIT, never evidence that the crosstalk is absent. A")
        emit("  correlation that cannot be resolved at n = 13 and n = 8, on an")
        emit("  assay whose ligand side is largely at background, is not a")
        emit("  correlation of zero.")
        emit()
        nominations = fdr.head(0).copy()
        for col in ("rationale", "framing"):
            nominations[col] = pd.Series(dtype=str)
        nominations["n_time_b"] = pd.Series(dtype=float)

    # ------------------------------------------------------------ assertions
    n_checks = 0
    if len(nominations) and (nominations["empirical_fdr"] > alpha).any():
        raise RuntimeError(
            "a nomination has an empirical FDR above the pre-registered alpha. "
            "There is no softer threshold in this phase to fall back on "
            "(ADR 0021 §6)."
        )
    n_checks += 1

    for site in sorted(adjacencies):
        n_site = int((nominations["site"] == site).sum()) if len(nominations) else 0
        if n_site > max_nom:
            raise RuntimeError(
                f"{site}: {n_site} nominations exceeds the cap of {max_nom}."
            )
    n_checks += 1

    # Every nomination must state both partners' detection. ADR 0014's term,
    # asserted rather than assumed.
    if len(nominations):
        for col in ("ligand_detection", "receptor_detection"):
            if nominations[col].isna().any() or (
                nominations[col].astype(str).str.len() == 0
            ).any():
                raise RuntimeError(
                    f"a nomination has an empty {col}. ADR 0014 makes the "
                    "per-site detection of BOTH partners part of every "
                    "nomination; a row without it is not reportable."
                )
    n_checks += 1
    emit(f"assertions: {n_checks} of {n_checks} passed "
         "(no nomination above alpha; cap honoured; both detections stated)")
    emit()

    nominations.to_csv(snakemake.output.nominations, sep="\t", index=False)

    by_site = {}
    for site in sorted(adjacencies):
        f = fdr[fdr["site"] == site]
        m = membership[(membership["site"] == site) & (membership["table"] == "primary")]
        by_site[site] = {
            "n_interactions_in_database": int(m["interaction_name"].nunique()),
            "n_admitted_rows": int(len(f)),
            "n_not_measurable_rows": int((~m["admitted"]).sum()),
            "n_clears_empirical_fdr": int(f["clears_empirical_fdr"].sum()),
            "n_clears_abundance_matched": int(f["clears_abundance_matched"].sum()),
            "n_nominated": (
                int((nominations["site"] == site).sum()) if len(nominations) else 0
            ),
        }
        emit(f"{site}: {by_site[site]}")

    summary = {
        "task": "P4-T5",
        "aim": "A5",
        "exploratory": True,
        "exploratory_adr": "ADR 0014",
        "pre_registered": "ADR 0021 §6",
        "fdr_alpha": alpha,
        "max_nominations_per_site": max_nom,
        "n_nominations": int(len(nominations)),
        "by_site": by_site,
        "framing": (
            "HYPOTHESIS GENERATION. No A5 result may be a headline claim "
            "(ADR 0014). Every row states the per-site detection of both "
            "partners; brain is TIME-B n = 8 (hard constraint 8)."
        ),
        "empty_is_a_result": (
            "If nothing survives, that IS the deliverable. Gate 4: 'no LR pair "
            "exceeded chance expectation at n=13' is an honest, useful result "
            "and a better outcome than a ranked list you can't defend. This "
            "script has no softer threshold to fall back on."
        ),
        "gate_note": (
            "Gate 4 turns on P4-T4, not on this table. A surviving list and an "
            "empty list both pass; a ranked list without a null does not."
        ),
        "null_reading": (
            "An absence here is an ASSAY-SENSITIVITY LIMIT, never evidence that "
            "the crosstalk is absent. `n_not_measurable_rows` and "
            "`n_clears_empirical_fdr` are counted separately because 'could not "
            "be measured' and 'did not beat chance' are different findings."
        ),
        "rationale_rule": (
            "Looked up from the database's own pathway_name, annotation and "
            "evidence string — never generated. Prose written after seeing "
            "which pairs correlated is how a ranked list becomes a story."
        ),
        "n_assertions": n_checks,
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit()
    emit(f"wrote {snakemake.output.nominations}, {snakemake.output.summary}")
