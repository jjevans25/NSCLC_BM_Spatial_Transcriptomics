"""P3-T2 — the checkpoint detection audit, per gene per compartment.

Owner task: P3-T2. Driven by rule p3t2_checkpoint_detection.

**This table is the Phase 3 deliverable, not a gate to get past.** ADR 0008
demoted A4 to exploratory because most of the checkpoint panel sits at
background in brain, and Gate 3 says in terms that "if everything is
not-assessable, that's still a legitimate finding". What would make the phase
weak is not a null — it is a null reported without the detection table that
licenses it. This script builds that table.

Three things it does NOT do, each for a specific reason:

  * **It does not invent a second detection rule.** A gene is detected in an AOI
    when its `layers['q3']` value exceeds `qc.detection_background_multiple`
    times that AOI's NegProbe-WTX level — ADR 0007's definition, already
    implemented in score_signatures.py. Two detection rules in one project is
    one too many.

  * **It does not key on `compartment`.** It keys on `aoi_code`. `compartment`
    is DEGENERATE across sites — L, LB and mLN are all labelled `tumour` — so
    grouping by it would silently pool lung and brain tumour AOIs and destroy
    the exact contrast Phase 3 exists to measure.

  * **It does not drop QC-flagged AOIs.** The policy is flag-don't-drop
    (CLAUDE.md); the five flagged AOIs are listed in qc_excluded.tsv as
    `flagged_retained`. They are scored, counted per compartment, and a
    `detection_rate_qc_clean` column reports what the rate would be without
    them — so TBME's low detection can be attributed to the compartment or to
    its four flagged AOIs rather than left ambiguous. That column is REPORTED,
    never used to gate `assessable`.

The reporting rule this table exists to enforce (ADR 0008 point 4): a gene below
the floor in a compartment is **"not assessable in <compartment>"**, NEVER
"lower in <compartment>". Brain background is HIGHER (NegProbe-WTX median log2
4.96 in L rising to 5.67 in BC), so a near-background gene reads as depleted in
brain artefactually. Phase 2 observed exactly that: `exhaustion` and `tls`
returned large, nominally significant "reductions in brain" that are not
reportable.

The audit is deliberately wider than the models. All seven compartments are
audited; only the four spanning both DSP runs are modelled (ADR 0016). The
`modelled` column marks which is which.
"""

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

FMT = "%.17g"

# The design table (PROJECT_PLAN §2.1) is authoritative. If these move, the
# subset is wrong or QC dropped something, and either way the power floor this
# phase was planned against no longer holds.
EXPECTED_AOI = {
    "L": 30,
    "LB": 27,
    "mLN": 13,
    "TBME": 20,
    "TIME-L": 15,
    "TIME-B": 8,
    "BC": 7,
}

# ADR 0008's reconnaissance table, transcribed. Eight of the nine panel genes
# appear there (CD276 does not). These counts were produced by a separate
# computation at Gate 0, so exact agreement is real evidence that this script's
# detection code is right — and a disagreement means one of the two is wrong,
# which has to be resolved before any model runs.
#   gene -> (TIME-L detected, TIME-B detected, TIME-L ratio, TIME-B ratio)
ADR0008_RECON = {
    "PDCD1": (8, 3, 2.03, 1.90),
    "CD274": (13, 5, 2.89, 2.32),
    "CTLA4": (10, 2, 2.76, 1.78),
    "LAG3": (6, 3, 1.97, 1.74),
    "HAVCR2": (10, 7, 2.20, 2.34),
    "TIGIT": (8, 1, 2.16, 1.41),
    "IDO1": (8, 2, 2.20, 1.23),
    "VSIR": (14, 7, 3.99, 3.74),
}

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    checkpoints = snakemake.params.checkpoints
    audit_codes = list(snakemake.params.audit_aoi_codes)
    model_codes = list(snakemake.params.model_aoi_codes)
    detection_cfg = snakemake.params.detection
    background_multiple = snakemake.params.background_multiple

    gene_floor = detection_cfg["detected_in_aoi_fraction"]
    panel_floor = detection_cfg["min_assessable_fraction"]

    emit("P3-T2 — checkpoint detection audit")
    emit(f"  detection rule      q3 > {background_multiple} x that AOI's NegProbe-WTX (ADR 0007)")
    emit(f"  gene floor          detected in >= {gene_floor} of a compartment's AOIs")
    emit(f"  panel floor         >= {panel_floor} of genes present (ADR 0016, PRE-REGISTERED)")
    emit(f"  audited             {audit_codes}")
    emit(f"  modelled            {model_codes}")
    emit()

    # ---------------------------------------------------------------- panel
    membership = pd.read_csv(snakemake.input.membership, sep="\t")
    genes = sorted(checkpoints)
    if sorted(membership["gene"]) != genes:
        raise RuntimeError(
            "checkpoint_membership.tsv and the config panel disagree on gene "
            f"membership: {sorted(membership['gene'])} vs {genes}. A panel that "
            "changes between P3-T1 and P3-T2 is a panel chosen after seeing a "
            "result."
        )
    alias_of = dict(zip(membership["gene"], membership["alias"]))
    emit(f"panel: {len(genes)} genes (PRE-REGISTERED, ADR 0015)")

    # --------------------------------------------------------------- subset
    adata = ad.read_h5ad(snakemake.input.h5ad)
    keep = adata.obs["aoi_code"].isin(audit_codes).to_numpy()
    sub = adata[keep].copy()

    counts = sub.obs["aoi_code"].value_counts().to_dict()
    if counts != EXPECTED_AOI:
        raise RuntimeError(
            f"expected {EXPECTED_AOI} AOIs, got {counts}. PROJECT_PLAN §2.1 is "
            "authoritative — do not adjust the expectation to match the data."
        )
    emit(f"AOIs audited: {sub.n_obs} (the full design table)")

    flagged_total = int(sub.obs["qc_flag"].sum())
    emit(
        f"QC-flagged AOIs retained: {flagged_total} "
        "(flag-don't-drop; see results/tables/qc_excluded.tsv)"
    )
    emit()

    # ------------------------------------------------------------- detection
    # Background-relative because the matrix has no zeros at all (Q1). This is
    # the rule from score_signatures.py, lifted rather than rewritten.
    q3 = pd.DataFrame(
        np.asarray(sub.layers["q3"]), index=sub.obs_names, columns=sub.var_names
    )[genes]
    x_log = pd.DataFrame(
        np.asarray(sub.X), index=sub.obs_names, columns=sub.var_names
    )[genes]

    negprobe = sub.obs["negprobe"].to_numpy()
    threshold = negprobe[:, None] * background_multiple
    detected = pd.DataFrame(
        q3.to_numpy() > threshold, index=sub.obs_names, columns=genes
    )
    ratio = pd.DataFrame(
        q3.to_numpy() / negprobe[:, None], index=sub.obs_names, columns=genes
    )

    code_of = sub.obs["aoi_code"].astype(str)
    qc_flag = sub.obs["qc_flag"].astype(bool)

    rows = []
    for code in audit_codes:
        mask = (code_of == code).to_numpy()
        clean = mask & ~qc_flag.to_numpy()
        site = sub.obs.loc[mask, "site"].iloc[0]
        compartment = sub.obs.loc[mask, "compartment"].iloc[0]
        for gene in genes:
            det = detected.loc[mask, gene]
            n_det = int(det.sum())
            rate = float(det.mean())
            det_clean = detected.loc[clean, gene]
            vals = x_log.loc[mask, gene]
            rows.append(
                {
                    "gene": gene,
                    "alias": alias_of[gene],
                    "aoi_code": code,
                    "site": site,
                    "compartment": compartment,
                    "n_aoi": int(mask.sum()),
                    "n_qc_flagged": int(qc_flag.to_numpy()[mask].sum()),
                    "n_detected": n_det,
                    "detection_rate": rate,
                    "n_aoi_qc_clean": int(clean.sum()),
                    "detection_rate_qc_clean": (
                        float(det_clean.mean()) if clean.sum() else float("nan")
                    ),
                    "median_negprobe_ratio": float(ratio.loc[mask, gene].median()),
                    "mean_log2_all": float(vals.mean()),
                    "mean_log2_detected": (
                        float(vals[det.to_numpy()].mean()) if n_det else float("nan")
                    ),
                    "assessable": bool(rate >= gene_floor),
                    "modelled": code in model_codes,
                }
            )

    audit = pd.DataFrame(rows)

    expected_rows = len(genes) * len(audit_codes)
    if len(audit) != expected_rows:
        raise RuntimeError(f"expected {expected_rows} audit rows, got {len(audit)}")

    # ----------------------------------------------------- the audit, printed
    emit("detection per gene x compartment (NEVER pooled):")
    emit()
    header = "  gene     alias    " + "".join(f"{c:>9s}" for c in audit_codes)
    emit(header)
    emit("  " + "-" * (len(header) - 2))
    for gene in genes:
        cells = []
        for code in audit_codes:
            r = audit[(audit["gene"] == gene) & (audit["aoi_code"] == code)].iloc[0]
            cells.append(f"{r['n_detected']:>4d}/{r['n_aoi']:<4d}")
        emit(f"  {gene:8s} {alias_of[gene]:8s} " + "".join(f"{c:>9s}" for c in cells))
    emit()

    # --------------------------------------------------- cross-check 1: total
    # var['detected_in_n_aoi'] was computed by workflow/scripts/export_tsv.py at
    # P0-T7 with the same rule but a SEPARATE implementation, over all 120 AOIs.
    # The audit covers all 120, so the per-gene totals must agree exactly.
    emit("cross-check 1 — audit totals vs var['detected_in_n_aoi'] (P0-T7, independent):")
    totals = audit.groupby("gene")["n_detected"].sum()
    var_counts = adata.var["detected_in_n_aoi"]
    mismatches = []
    for gene in genes:
        got, want = int(totals[gene]), int(var_counts[gene])
        mark = "ok" if got == want else "MISMATCH"
        if got != want:
            mismatches.append((gene, got, want))
        emit(f"  {gene:8s} audit {got:4d}  var {want:4d}   {mark}")
    if mismatches:
        raise RuntimeError(
            "detection disagrees with var['detected_in_n_aoi']: "
            f"{mismatches}. Two implementations of the same rule (this script "
            "and export_tsv.py) do not agree, so one of them is wrong. Resolve "
            "before any model runs."
        )
    emit("  all 9 genes agree exactly")
    emit()

    # ------------------------------------------------ cross-check 2: ADR 0008
    emit("cross-check 2 — TIME-L / TIME-B counts vs ADR 0008's reconnaissance:")
    recon_mismatches = []
    ratio_drift = []
    for gene, (want_l, want_b, ratio_l, ratio_b) in sorted(ADR0008_RECON.items()):
        for code, want, want_ratio in (
            ("TIME-L", want_l, ratio_l),
            ("TIME-B", want_b, ratio_b),
        ):
            r = audit[(audit["gene"] == gene) & (audit["aoi_code"] == code)].iloc[0]
            got = int(r["n_detected"])
            if got != want:
                recon_mismatches.append((gene, code, got, want))
            got_ratio = float(r["median_negprobe_ratio"])
            drift = abs(got_ratio - want_ratio)
            if drift >= 0.005:
                ratio_drift.append((gene, code, round(got_ratio, 3), want_ratio))
            emit(
                f"  {gene:8s} {code:7s} detected {got:2d} vs {want:2d}"
                f"   ratio {got_ratio:5.2f} vs {want_ratio:5.2f}"
                f"   {'ok' if got == want else 'MISMATCH'}"
            )
    if recon_mismatches:
        raise RuntimeError(
            "detection disagrees with ADR 0008's reconnaissance table: "
            f"{recon_mismatches}. That table is not the model of record, but a "
            "divergence means one of the two is wrong. STOP and resolve it "
            "before any Phase 3 model runs."
        )
    emit("  all 16 comparisons agree exactly")
    if ratio_drift:
        # Soft: the ADR published 2 dp, so small drift is transcription, not error.
        emit(f"  NOTE ratio drift beyond 2 dp (logged, not raised): {ratio_drift}")
    emit()

    # ------------------------------------------------- panel-level assessment
    emit("panel assessability per compartment (ADR 0016 floors):")
    panel_rows = {}
    for code in audit_codes:
        sl = audit[audit["aoi_code"] == code]
        n_ass = int(sl["assessable"].sum())
        frac = n_ass / len(genes)
        ok = frac >= panel_floor
        absent = sorted(sl.loc[~sl["assessable"], "gene"])
        present = sorted(sl.loc[sl["assessable"], "gene"])
        panel_rows[code] = {
            "n_assessable": n_ass,
            "n_genes": len(genes),
            "fraction": frac,
            "panel_assessable": ok,
            "assessable_genes": present,
            "not_assessable_genes": absent,
            "modelled": code in model_codes,
        }
        mark = "" if ok else "   <-- PANEL NOT ASSESSABLE"
        emit(f"  {code:7s} {n_ass}/{len(genes)} = {frac:.2f}{mark}")
        if absent:
            emit(f"      not assessable here: {', '.join(absent)}")
    emit()

    emit("REPORTING RULE (ADR 0008 point 4): for every gene x compartment marked")
    emit("not assessable above, the finding is 'not assessable in <compartment>',")
    emit("NEVER 'lower in <compartment>'. Brain background is HIGHER (NegProbe-WTX")
    emit("median log2 4.96 in L rising to 5.67 in BC), so a near-background gene")
    emit("reads as depleted in brain artefactually. Any brain claim additionally")
    emit("states TIME-B n = 8 inline (hard constraint 8).")
    emit()

    # ----------------------------------------------------- long expression TSV
    # Restricted to the modelled compartments, so pooling across the single-run
    # compartments is impossible rather than merely discouraged (ADR 0016).
    model_mask = code_of.isin(model_codes).to_numpy()
    obs_cols = ["patient_id", "aoi_code", "compartment", "site", "dsp_run", "qc_flag"]
    long = (
        x_log.loc[model_mask]
        .join(sub.obs.loc[model_mask, obs_cols])
        .reset_index()
        .melt(
            id_vars=["aoi_label"] + obs_cols,
            value_vars=genes,
            var_name="gene",
            value_name="expression",
        )
    )
    keyed = list(zip(long["aoi_label"], long["gene"]))
    long["q3"] = [q3.at[a, g] for a, g in keyed]
    long["detected"] = [bool(detected.at[a, g]) for a, g in keyed]
    long["negprobe"] = sub.obs.loc[long["aoi_label"], "negprobe"].to_numpy()
    long["negprobe_log2"] = np.log2(long["negprobe"])
    long = long.sort_values(["gene", "aoi_code", "aoi_label"]).reset_index(drop=True)

    expected_long = len(genes) * sum(EXPECTED_AOI[c] for c in model_codes)
    if len(long) != expected_long:
        raise RuntimeError(
            f"expected {expected_long} expression rows, got {len(long)}"
        )

    audit.to_csv(snakemake.output.detection, sep="\t", index=False, float_format=FMT)
    long.to_csv(snakemake.output.expression, sep="\t", index=False, float_format=FMT)
    emit(f"wrote {len(audit)} audit rows and {len(long)} expression rows")

    summary = {
        "n_genes": len(genes),
        "genes": genes,
        "audit_aoi_codes": audit_codes,
        "model_aoi_codes": model_codes,
        "aoi_counts": {k: int(v) for k, v in sorted(counts.items())},
        "n_qc_flagged_retained": flagged_total,
        "detection_background_multiple": background_multiple,
        "detected_in_aoi_fraction": gene_floor,
        "min_assessable_fraction": panel_floor,
        "pre_registered": "ADR 0016",
        "panel_assessability": panel_rows,
        "crosscheck_var_detected_in_n_aoi": {
            "status": "exact agreement",
            "n_genes_compared": len(genes),
            "source": "workflow/scripts/export_tsv.py (P0-T7), independent implementation",
        },
        "crosscheck_adr0008_recon": {
            "status": "exact agreement",
            "n_comparisons": 2 * len(ADR0008_RECON),
            "ratio_drift_beyond_2dp": ratio_drift,
        },
        "figure_colour_basis": (
            "median_negprobe_ratio — background-normalised, so it is immune to "
            "the compartment background gradient by construction (ADR 0016 §5)"
        ),
        "reporting_rule": (
            "A gene below the floor in a compartment is 'not assessable in "
            "<compartment>', never 'lower in <compartment>' (ADR 0008 point 4)."
        ),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
