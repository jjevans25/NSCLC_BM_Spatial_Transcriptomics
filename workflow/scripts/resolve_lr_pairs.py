"""P4-T2 — resolve the pinned LR database against the panel, and apply the detection filter.

Owner task: P4-T2. Driven by rule p4t2d_resolve_lr_pairs.

The Phase 4 analogue of resolve_checkpoints.py plus checkpoint_detection.py, in
one rule because the two questions are one question here: an interaction is
usable only if BOTH partners are on the panel AND both clear background, and the
gap between "how many interactions the database has" and "how many survive that"
is the phase's first reportable number.

**That gap is the quantitative form of ADR 0014's claim** — "secreted ligands
and chemokines are systematically undetected" — so it is reported BY INTERACTION
CLASS, which is the whole reason CellChatDB was chosen over CellPhoneDB
(ADR 0021 §2). A count that lumps Secreted Signaling in with Cell-Cell Contact
would hide the finding inside the summary of it.

Three tables, and the split is pre-registered (ADR 0021 §4):

  membership   EVERY usable interaction, with resolution and per-compartment
               detection for both partners. Nothing is dropped from it — an
               interaction that fails is a row saying so, never an absence.
  filtered     ONE-TO-ONE interactions admitted per direction. The only Phase 4
               table that will carry an empirical FDR.
  exploratory  COMPLEX interactions, admitted only when EVERY subunit clears the
               floor. No FDR of any kind, ever. Exploratory-within-exploratory:
               no Phase 4 sentence may rest on it.

What it deliberately does NOT do:

  * **Write a third detection rule.** A gene is detected in an AOI when
    `layers['q3']` exceeds `qc.detection_background_multiple` x that AOI's
    NegProbe-WTX level — ADR 0007's definition, implemented in
    score_signatures.py and lifted by checkpoint_detection.py. This lifts it a
    third time rather than restating it. Two detection rules in one project is
    one too many; three would be farce.

  * **Rank anything.** Ranking is P4-T3's and the empirical FDR is P4-T4's.
    A database filtered by the thing it is about to be ranked on would be a
    panel chosen after seeing the answer.

  * **Compute detection over the paired subset.** Detection is computed over
    ALL AOIs of an `aoi_code` (L 30, TIME-L 15, LB 27, TIME-B 8), so Phase 4's
    coverage numbers are directly comparable with Phase 2's and Phase 3's —
    which are the only coverage numbers a reader has to calibrate against.

  * **Key on `compartment`.** It is degenerate across sites (L, LB and mLN are
    all `tumour`), so keying on it would pool lung and brain.

**Admission is evaluated PER DIRECTION, not per interaction.** CellChatDB is
directed, and "ligand on the tumour AOI, receptor on the paired immune AOI" is a
different biological claim from its reverse (ADR 0021 §9). A pair can be
admissible in one direction and not the other, because the two partners are
being asked about in different compartments.
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
    gene_floor = snakemake.params.detection["detected_in_aoi_fraction"]
    background_multiple = snakemake.params.background_multiple
    expect = snakemake.params.expect
    database = snakemake.params.database
    version = snakemake.params.version

    emit("P4-T2 — resolve the LR database and apply the detection filter")
    emit("=" * 70)
    emit()
    emit(f"  database            {database} {version} (pinned, ADR 0021 §2)")
    emit(f"  detection rule      q3 > {background_multiple} x that AOI's NegProbe-WTX (ADR 0007)")
    emit(f"  gene floor          detected in >= {gene_floor} of a compartment's AOIs")
    emit("  admission           BOTH partners clear the floor, each in the")
    emit("                      compartment that partner is measured in,")
    emit("                      evaluated PER DIRECTION (ADR 0021 §3, §9)")
    emit()
    emit("  A5 is EXPLORATORY (ADR 0014). The before/after gap below is the")
    emit("  quantitative form of ADR 0014's claim, not bookkeeping — a pair")
    emit("  that cannot be measured is NOT a pair that does not interact.")
    emit()

    # ---------------------------------------------------------------- inputs
    ix = pd.read_csv(
        snakemake.input.interactions, sep="\t", float_precision="round_trip"
    )
    cx = pd.read_csv(
        snakemake.input.complexes, sep="\t", float_precision="round_trip"
    )
    emit(f"interactions: {len(ix)} usable (from {snakemake.input.interactions})")
    emit(f"complexes:    {len(cx)} definitions")

    if len(ix) != expect["n_usable"]:
        raise RuntimeError(
            f"got {len(ix)} usable interactions, config/ligand_receptor.yaml "
            f"expects {expect['n_usable']}. The R export already asserts this; "
            "a disagreement here means the TSV was truncated between rules."
        )

    adata = ad.read_h5ad(snakemake.input.h5ad)
    measured = set(map(str, adata.var_names))
    emit(f"genes measured: {len(measured)} (from {snakemake.input.h5ad})")
    emit()

    codes = sorted({c for a in adjacencies.values() for c in a.values()})
    emit(f"compartments: {codes}")

    # ------------------------------------------------------------- detection
    # ADR 0007's rule, lifted from checkpoint_detection.py. Restricted to the
    # genes the database names, because computing it for all 18,694 would be
    # 18,694 numbers nobody reads.
    named = set()
    for col in ("ligand_subunits", "receptor_subunits"):
        for cell in ix[col].astype(str):
            named.update(s for s in cell.split(";") if s)
    panel = sorted(named & measured)
    emit(f"symbols named by the database: {len(named)}")
    emit(f"  of those, on the measured panel: {len(panel)}")
    emit(f"  absent from the panel:           {len(named) - len(panel)}")
    emit()

    sub = adata[adata.obs["aoi_code"].astype(str).isin(codes)].copy()
    q3 = pd.DataFrame(
        np.asarray(sub.layers["q3"]), index=sub.obs_names, columns=sub.var_names
    )[panel]
    negprobe = sub.obs["negprobe"].to_numpy()
    detected = pd.DataFrame(
        q3.to_numpy() > negprobe[:, None] * background_multiple,
        index=sub.obs_names,
        columns=panel,
    )
    code_of = sub.obs["aoi_code"].astype(str)

    det_rate = {}
    det_count = {}
    n_aoi = {}
    for code in codes:
        mask = (code_of == code).to_numpy()
        n_aoi[code] = int(mask.sum())
        sel = detected.loc[mask]
        det_rate[code] = sel.mean(axis=0)
        det_count[code] = sel.sum(axis=0).astype(int)
        emit(f"  {code:7s} {n_aoi[code]:3d} AOIs, "
             f"{int((det_rate[code] >= gene_floor).sum()):4d} of {len(panel)} "
             f"named symbols clear the {gene_floor} floor")
    emit()

    # Cross-check against an INDEPENDENT implementation: var['detected_in_n_aoi']
    # was written by export_tsv.py at P0-T7 over all 120 AOIs. Summing this
    # rule's per-compartment counts over the four compartments cannot exceed it.
    var_counts = adata.var["detected_in_n_aoi"]
    totals = sum(det_count[c] for c in codes)
    over = [g for g in panel if int(totals[g]) > int(var_counts[g])]
    if over:
        raise RuntimeError(
            "detection exceeds var['detected_in_n_aoi'] for "
            f"{over[:10]}. This rule counts a subset of the 120 AOIs that "
            "P0-T7's independent implementation counted over all of them, so "
            "it can never be larger. One of the two is wrong."
        )
    emit(f"cross-check — per-compartment counts <= var['detected_in_n_aoi'] "
         f"(P0-T7, independent) for all {len(panel)} symbols: ok")
    emit()

    # ------------------------------------------------- resolve and admit
    def subunits(cell):
        return [s for s in str(cell).split(";") if s]

    def clears(genes, code):
        """Every subunit clears the floor in this compartment."""
        return all(
            g in measured and float(det_rate[code][g]) >= gene_floor for g in genes
        )

    def describe(genes, code):
        """n_detected/n_aoi per subunit, as a string — never a bare rate."""
        parts = []
        for g in genes:
            if g not in measured:
                parts.append(f"{g}:not_on_panel")
            else:
                parts.append(f"{g}:{int(det_count[code][g])}/{n_aoi[code]}")
        return ";".join(parts)

    rows = []
    for r in ix.itertuples(index=False):
        lig, rec = subunits(r.ligand_subunits), subunits(r.receptor_subunits)
        lig_on = all(g in measured for g in lig)
        rec_on = all(g in measured for g in rec)
        for site in sorted(adjacencies):
            tum = adjacencies[site]["tumour"]
            imm = adjacencies[site]["immune"]
            # Both directions. `ligand_side` names where the LIGAND is measured.
            for direction, lig_code, rec_code in (
                ("ligand_tumour__receptor_immune", tum, imm),
                ("ligand_immune__receptor_tumour", imm, tum),
            ):
                admitted = (
                    lig_on and rec_on and clears(lig, lig_code) and clears(rec, rec_code)
                )
                if not (lig_on and rec_on):
                    status = "not_on_panel"
                elif admitted:
                    status = "admitted"
                else:
                    status = "below_detection_floor"
                rows.append(
                    {
                        "interaction_name": r.interaction_name,
                        "interaction_name_2": r.interaction_name_2,
                        "pathway_name": r.pathway_name,
                        "annotation": r.annotation,
                        "site": site,
                        "direction": direction,
                        "ligand": r.ligand,
                        "receptor": r.receptor,
                        "ligand_subunits": ";".join(lig),
                        "receptor_subunits": ";".join(rec),
                        "ligand_aoi_code": lig_code,
                        "receptor_aoi_code": rec_code,
                        "is_one_to_one": bool(r.is_one_to_one),
                        "table": "primary" if r.is_one_to_one else "exploratory",
                        "ligand_on_panel": lig_on,
                        "receptor_on_panel": rec_on,
                        "ligand_detection": describe(lig, lig_code),
                        "receptor_detection": describe(rec, rec_code),
                        "ligand_min_detection_rate": (
                            min(float(det_rate[lig_code][g]) for g in lig)
                            if lig_on else float("nan")
                        ),
                        "receptor_min_detection_rate": (
                            min(float(det_rate[rec_code][g]) for g in rec)
                            if rec_on else float("nan")
                        ),
                        "admitted": admitted,
                        "status": status,
                        "ligand_secreted_type": r.ligand_secreted_type,
                        "ligand_transmembrane": r.ligand_transmembrane,
                        "receptor_transmembrane": r.receptor_transmembrane,
                        "evidence": r.evidence,
                    }
                )

    membership = pd.DataFrame(rows)

    # ------------------------------------------------------- the reportable gap
    emit("THE BEFORE/AFTER GAP, BY INTERACTION CLASS (ADR 0014, made numeric)")
    emit("-" * 70)
    emit()
    for site in sorted(adjacencies):
        emit(f"{site}  ({adjacencies[site]['tumour']} <-> {adjacencies[site]['immune']})")
        for table in ("primary", "exploratory"):
            m = membership[
                (membership["site"] == site) & (membership["table"] == table)
            ]
            label = "1:1 (primary)" if table == "primary" else "complex (exploratory)"
            emit(f"  {label}")
            emit(f"    {'class':<24s} {'in db':>6s} {'on panel':>9s} "
                 f"{'admitted':>9s} {'%':>6s}")
            for cls in sorted(m["annotation"].unique()):
                c = m[m["annotation"] == cls]
                # Per-interaction, not per-direction, for the "in db" column.
                n_db = c["interaction_name"].nunique()
                n_panel = c.loc[
                    c["ligand_on_panel"] & c["receptor_on_panel"], "interaction_name"
                ].nunique()
                n_adm = c.loc[c["admitted"], "interaction_name"].nunique()
                pct = 100.0 * n_adm / n_db if n_db else 0.0
                emit(f"    {cls:<24s} {n_db:6d} {n_panel:9d} {n_adm:9d} {pct:5.1f}%")
            emit()
    emit("  Read this the way ADR 0014 would: a class that vanishes here is a")
    emit("  class this ASSAY cannot resolve, not a class that is biologically")
    emit("  absent. The anticipated null is an assay-sensitivity limit.")
    emit()

    filtered = membership[
        (membership["table"] == "primary") & membership["admitted"]
    ].copy()
    exploratory = membership[
        (membership["table"] == "exploratory") & membership["admitted"]
    ].copy()

    # No q-value may ever reach the exploratory table (ADR 0021 §4). Asserted
    # here rather than trusted downstream — checkpoint_paired_check.py's
    # no-test-statistic assertion is the pattern.
    forbidden = {"q_bh", "q_value", "empirical_fdr", "p_raw", "rho"}
    leaked = forbidden & set(exploratory.columns)
    if leaked:
        raise RuntimeError(
            f"the exploratory table carries {sorted(leaked)}. ADR 0021 §4 gives "
            "it no FDR of any kind and no test statistic; this rule produces "
            "neither."
        )

    emit(f"admitted, 1:1 primary      {len(filtered):5d} direction-rows "
         f"({filtered['interaction_name'].nunique()} distinct interactions)")
    emit(f"admitted, complex exploratory {len(exploratory):5d} direction-rows "
         f"({exploratory['interaction_name'].nunique()} distinct interactions)")
    emit()

    # ------------------------------------------------------------ assertions
    n_checks = 0

    # Resolution counts are a property of the database crossed with the panel,
    # fixed and declared in advance. What survives DETECTION is deliberately
    # NOT in the manifest — that is the result.
    got_1to1 = int(
        membership.loc[
            membership["is_one_to_one"], "interaction_name"
        ].nunique()
    )
    if got_1to1 != expect["n_one_to_one"]:
        raise RuntimeError(
            f"{got_1to1} one-to-one interactions, manifest expects "
            f"{expect['n_one_to_one']}."
        )
    n_checks += 1

    got_measured = int(
        membership.loc[
            membership["is_one_to_one"]
            & membership["ligand_on_panel"]
            & membership["receptor_on_panel"],
            "interaction_name",
        ].nunique()
    )
    if got_measured != expect["n_one_to_one_measured"]:
        raise RuntimeError(
            f"{got_measured} one-to-one interactions have both partners on the "
            f"panel, manifest expects {expect['n_one_to_one_measured']}. This "
            "is a property of the database crossed with the 18,694 measured "
            "symbols, so it is fixed — a change means the pin or the matrix "
            "moved, and both are stop-and-ask."
        )
    n_checks += 1

    got_cx_measured = int(
        membership.loc[
            ~membership["is_one_to_one"]
            & membership["ligand_on_panel"]
            & membership["receptor_on_panel"],
            "interaction_name",
        ].nunique()
    )
    if got_cx_measured != expect["n_complex_measured"]:
        raise RuntimeError(
            f"{got_cx_measured} complex interactions have every subunit on the "
            f"panel, manifest expects {expect['n_complex_measured']}."
        )
    n_checks += 1

    # Every admitted row must clear the floor on both sides — the filter
    # re-derived from the rates rather than trusted from the boolean.
    for r in membership[membership["admitted"]].itertuples(index=False):
        if (
            r.ligand_min_detection_rate < gene_floor
            or r.receptor_min_detection_rate < gene_floor
        ):
            raise RuntimeError(
                f"{r.interaction_name} ({r.site}, {r.direction}) is admitted but "
                f"has min detection {r.ligand_min_detection_rate:.3f} / "
                f"{r.receptor_min_detection_rate:.3f} against a floor of "
                f"{gene_floor}."
            )
    n_checks += 1

    # Both directions must exist for every interaction at every site.
    per = membership.groupby(["interaction_name", "site"])["direction"].nunique()
    if not (per == 2).all():
        raise RuntimeError(
            "an interaction is missing a direction. CellChatDB is directed and "
            "ADR 0021 §9 tests both; a one-sided row would silently make a "
            "different claim."
        )
    n_checks += 1

    emit(f"assertions: {n_checks} of {n_checks} passed "
         "(1:1 count; 1:1-on-panel count; complex-on-panel count; admitted "
         "rows clear the floor; both directions present)")
    emit()

    membership.to_csv(snakemake.output.membership, sep="\t", index=False)
    filtered.to_csv(snakemake.output.filtered, sep="\t", index=False)
    exploratory.to_csv(snakemake.output.exploratory, sep="\t", index=False)

    def by_class(frame, mask=None):
        f = frame if mask is None else frame[mask]
        return {
            cls: int(f.loc[f["annotation"] == cls, "interaction_name"].nunique())
            for cls in sorted(frame["annotation"].unique())
        }

    summary = {
        "task": "P4-T2",
        "aim": "A5",
        "exploratory": True,
        "exploratory_adr": "ADR 0014",
        "pre_registered": "ADR 0021 §§2-4",
        "database": database,
        "version": version,
        "detection_rule": (
            f"q3 > {background_multiple} x that AOI's NegProbe-WTX (ADR 0007), "
            "lifted from checkpoint_detection.py, not re-derived"
        ),
        "detected_in_aoi_fraction": gene_floor,
        "n_interactions_usable": int(membership["interaction_name"].nunique()),
        "n_symbols_named": len(named),
        "n_symbols_on_panel": len(panel),
        "n_symbols_absent": len(named) - len(panel),
        "n_aoi_per_compartment": {c: n_aoi[c] for c in codes},
        "by_site": {
            site: {
                table: {
                    "in_database": by_class(
                        membership[
                            (membership["site"] == site)
                            & (membership["table"] == table)
                        ]
                    ),
                    "admitted": by_class(
                        membership[
                            (membership["site"] == site)
                            & (membership["table"] == table)
                        ],
                        membership[
                            (membership["site"] == site)
                            & (membership["table"] == table)
                        ]["admitted"],
                    ),
                }
                for table in ("primary", "exploratory")
            }
            for site in sorted(adjacencies)
        },
        "n_admitted_primary_rows": len(filtered),
        "n_admitted_exploratory_rows": len(exploratory),
        "primary_rule": (
            "one-to-one interactions, both partners clearing the floor in the "
            "compartment that partner is measured in, per direction. The only "
            "Phase 4 table that carries an empirical FDR."
        ),
        "exploratory_rule": (
            "complex interactions, every subunit clearing the floor. NO FDR of "
            "any kind, ever. No Phase 4 sentence may rest on it (ADR 0021 §4)."
        ),
        "reporting_rule": (
            "The before/after gap is reported BY INTERACTION CLASS and is the "
            "quantitative form of ADR 0014's claim. A class that vanishes is a "
            "class this assay cannot resolve, NEVER a class that is "
            "biologically absent. TIME-B n = 8 (hard constraint 8)."
        ),
        "n_assertions": n_checks,
        "n_time_b": n_aoi.get("TIME-B"),
        "h5ad": str(snakemake.input.h5ad),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit(
        f"wrote {snakemake.output.membership}, {snakemake.output.filtered}, "
        f"{snakemake.output.exploratory}, {snakemake.output.summary}"
    )
