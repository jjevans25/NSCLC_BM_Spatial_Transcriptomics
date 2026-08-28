"""P2-T2 — score the signatures on the TIME AOIs, two ways, with coverage.

Owner task: P2-T2. Driven by rule p2t2_score_signatures.

Scores the six signatures on the 23 TIME AOIs (TIME-L 15 + TIME-B 8) by two
methods, and reports per-set detection coverage per site.

Choices worth stating, because each could reasonably go the other way:

  * **Two methods, neither of which is the tiebreaker.** ssGSEA is rank-based
    and per-sample; the z-score mean is a linear summary. They fail in
    different directions, so where they agree the direction is not an artefact
    of either. At TIME-B n = 8 that agreement is the robustness evidence, and
    it is the reason the plan asks for both rather than the better one.

  * **`min_size` is passed explicitly.** `gseapy.ssgsea` defaults it to 15,
    and every set here is 4-13 genes. Checked against gseapy 1.3.1: the default
    raises `LookupError("No gene sets passed through filtering condition")`
    rather than returning an empty frame, so this fails loudly rather than
    quietly. The empty-result guard below is kept anyway, because that is a
    behaviour of the installed version rather than a documented contract.

  * **Scored on `X`, and the layer choice does not matter.** ssGSEA ranks genes
    within a sample, and log2(Q3+1) is monotone in Q3, so `X` and `layers['q3']`
    give identical ranks and identical enrichment scores. `X` is used for both
    methods so the two share an input; the z-score mean *is* transform-sensitive
    and log space is the right one for it.

  * **z-scored across the 23 subset AOIs, not all 120.** Scoring, like the
    model, happens within compartment (ADR 0009). Standardising against L, LB,
    TBME and BC would centre these scores on tissue the contrast excludes.

  * **Every set member is scored, including undetected ones.** Dropping
    below-background genes per site would score the two sites on *different*
    gene sets and make the contrast uninterpretable. Detection is reported
    beside the scores instead, as coverage, and the coverage table is what
    licenses or forbids a claim.

  * **Coverage is per site, never pooled.** P0-T5 established that detection
    varies systematically by compartment, and ADR 0008 established that brain
    background is HIGHER — so a near-background gene reads as depleted in brain
    artefactually. A set below the coverage floor at a site is "not assessable"
    there; the phrase "lower in brain" is forbidden for it.

Detection reuses the rule already in the project — a gene is detected in an AOI
when its Q3 value exceeds `qc.detection_background_multiple` times that AOI's
NegProbe-WTX level — rather than inventing a second one.
"""

import json
from pathlib import Path

import anndata as ad
import gseapy as gp
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
    scoring = snakemake.params.scoring
    background_multiple = snakemake.params.background_multiple
    seed = snakemake.params.seed

    # ---------------------------------------------------------------- subset
    adata = ad.read_h5ad(snakemake.input.h5ad)
    keep = adata.obs["aoi_code"].isin(aoi_codes).to_numpy()
    sub = adata[keep].copy()

    counts = sub.obs["aoi_code"].value_counts().to_dict()
    emit(f"AOIs selected: {sub.n_obs} across {aoi_codes} -> {counts}")

    # The design table is authoritative (CLAUDE.md). If these numbers move, the
    # subset is wrong or QC dropped something, and either way the power floor
    # this phase was planned against no longer holds.
    expected = {"TIME-L": 15, "TIME-B": 8}
    if counts != expected:
        raise RuntimeError(
            f"expected {expected} AOIs, got {counts}. PROJECT_PLAN §2.1 is "
            "authoritative — do not adjust the expectation to match the data."
        )

    flagged = sub.obs["qc_flag"].sum() if "qc_flag" in sub.obs else 0
    emit(f"QC-flagged AOIs in the subset: {int(flagged)} (Phase 1 measured 0)")
    if flagged:
        emit(
            "  NOTE: flag-don't-drop policy means these are still scored. "
            "A non-zero count here is a regression to investigate."
        )
    emit()

    # ------------------------------------------------------------- detection
    # Background-relative because the matrix has no zeros at all (Q1).
    q3 = pd.DataFrame(
        np.asarray(sub.layers["q3"]), index=sub.obs_names, columns=sub.var_names
    )
    threshold = sub.obs["negprobe"].to_numpy()[:, None] * background_multiple
    detected = pd.DataFrame(
        q3.to_numpy() > threshold, index=sub.obs_names, columns=sub.var_names
    )

    site_of = sub.obs["aoi_code"].astype(str)
    coverage_rows = []
    for name in sorted(signatures):
        genes = list(signatures[name]["genes"])
        for code in aoi_codes:
            mask = (site_of == code).to_numpy()
            frac = detected.loc[mask, genes].mean(axis=0)
            present = frac >= scoring["detected_in_aoi_fraction"]
            covered = float(present.mean())
            coverage_rows.append(
                {
                    "signature": name,
                    "aoi_code": code,
                    "site": sub.obs.loc[mask, "site"].iloc[0],
                    "n_aoi": int(mask.sum()),
                    "n_genes": len(genes),
                    "n_genes_present": int(present.sum()),
                    "coverage": covered,
                    "assessable": bool(covered >= scoring["min_detected_fraction"]),
                    "genes_absent": ",".join(sorted(frac.index[~present])),
                }
            )

    coverage = pd.DataFrame(coverage_rows)
    emit("detection coverage (per set, per site — never pooled):")
    for _, r in coverage.iterrows():
        mark = "" if r["assessable"] else "   <-- NOT ASSESSABLE"
        emit(
            f"  {r['signature']:22s} {r['aoi_code']:7s} "
            f"{r['n_genes_present']:2d}/{r['n_genes']:2d} = {r['coverage']:.2f}{mark}"
        )
        if r["genes_absent"]:
            emit(f"      below background: {r['genes_absent']}")
    emit()

    not_assessable = coverage.loc[~coverage["assessable"]]
    if len(not_assessable):
        emit("REPORTING RULE (ADR 0008): for the rows above marked NOT ASSESSABLE,")
        emit("the finding is 'not assessable in <site>', never 'lower in <site>'.")
        emit()

    # ---------------------------------------------------------------- scores
    expr = pd.DataFrame(
        np.asarray(sub.X), index=sub.obs_names, columns=sub.var_names
    ).T  # genes x samples, which is what gseapy wants

    gene_sets = {k: list(v["genes"]) for k, v in signatures.items()}

    # ssGSEA. min_size is explicit; see the docstring.
    ss = gp.ssgsea(
        data=expr,
        gene_sets=gene_sets,
        outdir=None,
        min_size=scoring["min_set_size"],
        max_size=max(len(v) for v in gene_sets.values()),
        permutation_num=0,
        no_plot=True,
        seed=seed,
        threads=snakemake.threads,
    )
    res = ss.res2d.copy()
    emit(f"ssgsea returned {len(res)} rows, columns {list(res.columns)}")
    if res.empty:
        raise RuntimeError(
            "gseapy.ssgsea returned no rows. The usual cause is min_size "
            f"({scoring['min_set_size']}) exceeding a set's length."
        )

    ssgsea_scores = (
        res.rename(columns={"Name": "aoi_label", "Term": "signature"})
        .assign(score=lambda d: d["NES"].astype(float), method="ssgsea")
        .loc[:, ["aoi_label", "signature", "method", "score"]]
    )

    # z-score mean, standardised across the 23 subset AOIs only.
    z = (expr.T - expr.T.mean(axis=0)) / expr.T.std(axis=0, ddof=1)
    zrows = []
    for name, genes in gene_sets.items():
        zrows.append(
            pd.DataFrame(
                {
                    "aoi_label": z.index,
                    "signature": name,
                    "method": "zscore",
                    "score": z[genes].mean(axis=1).to_numpy(),
                }
            )
        )
    zscore_scores = pd.concat(zrows, ignore_index=True)

    scores = pd.concat([ssgsea_scores, zscore_scores], ignore_index=True)

    meta = sub.obs.loc[:, ["aoi_code", "site", "patient_id", "dsp_run"]].reset_index()
    meta = meta.rename(columns={meta.columns[0]: "aoi_label"})
    scores = scores.merge(meta, on="aoi_label", how="left", validate="many_to_one")
    if scores["aoi_code"].isna().any():
        raise RuntimeError("a score row failed to join to its AOI metadata")

    scores = scores.sort_values(["method", "signature", "aoi_code", "aoi_label"])
    scores.to_csv(snakemake.output.scores, sep="\t", index=False, float_format=FMT)
    coverage.to_csv(snakemake.output.coverage, sep="\t", index=False, float_format=FMT)

    # ------------------------------------------------- method agreement, per set
    emit("Spearman agreement between the two scoring methods (the robustness")
    emit("evidence at TIME-B n = 8):")
    agreement = {}
    wide = scores.pivot_table(
        index=["aoi_label", "signature"], columns="method", values="score"
    ).reset_index()
    for name in sorted(gene_sets):
        d = wide.loc[wide["signature"] == name]
        rho = float(d["ssgsea"].corr(d["zscore"], method="spearman"))
        agreement[name] = rho
        emit(f"  {name:22s} rho = {rho:+.3f}")
    emit()

    summary = {
        "n_aoi": int(sub.n_obs),
        "aoi_counts": {k: int(v) for k, v in counts.items()},
        "n_qc_flagged": int(flagged),
        "methods": list(scoring["methods"]),
        "min_set_size": scoring["min_set_size"],
        "detection_background_multiple": background_multiple,
        "detected_in_aoi_fraction": scoring["detected_in_aoi_fraction"],
        "min_detected_fraction": scoring["min_detected_fraction"],
        "method_agreement_spearman": agreement,
        "coverage": coverage.to_dict(orient="records"),
        "not_assessable": [
            {"signature": r["signature"], "aoi_code": r["aoi_code"], "coverage": r["coverage"]}
            for _, r in not_assessable.iterrows()
        ],
        "seed": seed,
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit(f"wrote {len(scores)} score rows and {len(coverage)} coverage rows")
