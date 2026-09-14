"""P4-T6 — informal concordance against an external NSCLC Visium study. NOT validation.

Owner task: P4-T6. Driven by rule p4t6_external_crossref.

The comparator is De Zuani et al. 2024, Nat Commun 15:4388 (PMID 38782901),
pre-registered in `config/ligand_receptor.yaml → external_comparator` with the
paper's own sentence quoted, before any concordance was computed.

**THIS IS NOT VALIDATION, AND THE THREE REASONS ARE STRUCTURAL RATHER THAN
CAUTIONARY.** They are carried on every row of the output because a concordance
table read without them is a replication claim:

  1. **They used CellPhoneDB; this project uses CellChatDB.** A pair missing
     from our panel is a *database* difference, never a biological one.
  2. **They measured co-expression within a Visium spot**, tumour section
     against background section. This project measures a *cross-patient
     correlation between two compartments* and has **no coordinates at all**.
     Those are different quantities. Agreement is encouraging; disagreement is
     evidence against neither result.
  3. **They profiled primary lung (LUAD/LUSC).** The brain arm — the one this
     project exists to ask about — **has no comparator**, and `TIME-B` n = 8.

**Two of the comparators were known to be compromised before this script
existed**, which is why they are reported rather than quietly dropped:

  * **`LGALS9` is not among the 18,694 measured genes.** The panel carries
    `LGALS1`, `LGALS3`, `LGALS8` and `LGALS9C`; `LGALS9` itself is absent. It is
    reported `not_on_panel`. **`LGALS9C` is NOT substituted** — that would be a
    membership change made to keep a comparator alive, and adding a gene to an
    analysis panel is a stop-and-ask (ADR 0021 §8, ADR 0011, ADR 0015).
  * **`TIGIT` sits below the floor across the lung tumour compartment**, so
    `NECTIN2–TIGIT` cannot be evaluated in lung in the direction that matters.

That leaves `VEGFA–NRP1` as the one comparator where both partners are detected
nearly everywhere and the check is genuinely informative. The script does not
decide that in advance — it computes each status and reports it — but a reader
should know it going in.

**A pair we cannot evaluate is not a pair we disagree about.** The output keeps
`not_on_panel`, `not_in_database`, `below_detection_floor` and `evaluable` as
four separate statuses precisely so they cannot collapse into "did not
replicate".
"""

import json
from pathlib import Path

import anndata as ad
import pandas as pd

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    comp = snakemake.params.external_comparator
    alpha = snakemake.params.null_calibration["fdr_alpha"]
    adjacencies = snakemake.params.adjacencies

    emit("P4-T6 — informal concordance against an external study. NOT validation.")
    emit("=" * 70)
    emit()
    emit(f"  comparator   {comp['citation']}  ({comp['source']})")
    emit(f"  assay        {comp['assay']}")
    emit(f"  database     {comp['database_used']}  (this project: CellChatDB v2)")
    emit(f"  tissue       {comp['tissue']}")
    emit()
    emit("  The comparator's own words, quoted rather than paraphrased:")
    for line in (" ".join(comp["source_claim"].split())).split(". "):
        if line.strip():
            emit(f"    {line.strip()}")
    emit()

    membership = pd.read_csv(
        snakemake.input.membership, sep="\t", float_precision="round_trip"
    )
    fdr = pd.read_csv(snakemake.input.fdr, sep="\t", float_precision="round_trip")
    adata = ad.read_h5ad(snakemake.input.h5ad)
    measured = set(map(str, adata.var_names))
    emit(f"genes measured: {len(measured)}")
    emit()

    fdr_key = fdr.set_index(["interaction_name", "site", "direction"])

    rows = []
    for pair in comp["pairs"]:
        lig, rec, reported = pair["ligand"], pair["receptor"], pair["reported"]
        label = f"{lig}-{rec}"
        emit(f"{label}   comparator says: {reported}")

        lig_on, rec_on = lig in measured, rec in measured
        if not (lig_on and rec_on):
            absent = [g for g, on in ((lig, lig_on), (rec, rec_on)) if not on]
            # Name the paralogs that ARE present, so "not on the panel" is a
            # measured fact rather than a bare absence — and so nobody has to
            # go looking for a substitute. Substituting one is a stop-and-ask.
            for g in absent:
                stem = g.rstrip("0123456789")
                near = sorted(x for x in measured if x.startswith(stem[:5]))[:8]
                emit(f"  {g} is NOT on the measured panel.")
                emit(f"    present with a similar symbol: {near}")
            emit("    Reported as not_on_panel. NO PARALOG IS SUBSTITUTED — that")
            emit("    would be a membership change made to keep a comparator")
            emit("    alive (ADR 0021 §8).")

        for site in sorted(adjacencies):
            for direction in (
                "ligand_tumour__receptor_immune",
                "ligand_immune__receptor_tumour",
            ):
                # Either orientation. The comparator measured co-expression
                # within a Visium spot, which is UNDIRECTED — its "NRP1-VEGFA"
                # asserts no ligand/receptor role — so requiring our database's
                # orientation to match a hand-entered one would report a pair
                # the database contains as absent. It did exactly that for
                # CD96-NECTIN1 before this was order-agnostic.
                same = (membership["ligand"] == lig) & (
                    membership["receptor"] == rec
                )
                flipped = (membership["ligand"] == rec) & (
                    membership["receptor"] == lig
                )
                m = membership[
                    (same | flipped)
                    & (membership["site"] == site)
                    & (membership["direction"] == direction)
                ]
                row = {
                    "comparator_pair": label,
                    "comparator_reported": reported,
                    "comparator_citation": comp["citation"],
                    "comparator_source": comp["source"],
                    "comparator_assay": comp["assay"],
                    "comparator_database": comp["database_used"],
                    "comparator_tissue": comp["tissue"],
                    "ligand": lig,
                    "receptor": rec,
                    "site": site,
                    "direction": direction,
                    "ligand_on_panel": lig_on,
                    "receptor_on_panel": rec_on,
                    "in_our_database": bool(len(m)),
                    "db_ligand": "",
                    "db_receptor": "",
                    "orientation_matches_comparator": pd.NA,
                    "rho": float("nan"),
                    "empirical_fdr": float("nan"),
                    "n": pd.NA,
                    "n_detected_both": pd.NA,
                    "ligand_detection": "",
                    "receptor_detection": "",
                    "concordance": "not_assessable",
                    "is_validation": False,
                }
                if not (lig_on and rec_on):
                    row["status"] = "not_on_panel"
                elif not len(m):
                    row["status"] = "not_in_database"
                else:
                    r = m.iloc[0]
                    row["db_ligand"] = r["ligand"]
                    row["db_receptor"] = r["receptor"]
                    row["orientation_matches_comparator"] = bool(
                        r["ligand"] == lig
                    )
                    row["ligand_detection"] = r["ligand_detection"]
                    row["receptor_detection"] = r["receptor_detection"]
                    if not bool(r["admitted"]):
                        row["status"] = "below_detection_floor"
                    else:
                        row["status"] = "evaluable"
                        k = (r["interaction_name"], site, direction)
                        if k in fdr_key.index:
                            f = fdr_key.loc[k]
                            row["rho"] = float(f["rho"])
                            row["empirical_fdr"] = float(f["empirical_fdr"])
                            row["n"] = int(f["n"])
                            row["n_detected_both"] = int(f["n_detected_both"])
                            clears = bool(f["empirical_fdr"] <= alpha)
                            if reported == "enriched_in_tumour":
                                row["concordance"] = (
                                    "agrees" if clears else "not_recovered"
                                )
                            else:
                                row["concordance"] = (
                                    "agrees" if not clears else "diverges"
                                )
                rows.append(row)

        for site in sorted(adjacencies):
            sel = [r for r in rows if r["comparator_pair"] == label
                   and r["site"] == site]
            statuses = sorted({r["status"] for r in sel})
            tag = "  [TIME-B n = 8]" if site == "brain" else ""
            emit(f"  {site:5s} {', '.join(statuses)}{tag}")
            for r in sel:
                if r["status"] == "evaluable":
                    arrow = ("T->I" if r["direction"].startswith("ligand_tumour")
                             else "I->T")
                    emit(f"    [{arrow}] rho {r['rho']:+.3f}  FDR "
                         f"{r['empirical_fdr']:.4f}  n {r['n']}  "
                         f"detected both {r['n_detected_both']}/{r['n']}  "
                         f"-> {r['concordance']}")
                    emit(f"           ligand   {r['ligand_detection']}")
                    emit(f"           receptor {r['receptor_detection']}")
        emit()

    table = pd.DataFrame(rows)

    emit("WHAT AGREEMENT AND DISAGREEMENT EACH MEAN")
    emit("-" * 70)
    emit()
    emit("  AGREEMENT is encouraging and is not confirmation. Two assays that")
    emit("  measure different quantities can agree for a shared reason that is")
    emit("  neither of the biologies proposed — high expression in the same")
    emit("  tissue will do it. That is exactly what ADR 0021 §6's")
    emit("  abundance-matched Null B is for, and a concordant pair that fails")
    emit("  Null B has not been corroborated.")
    emit()
    emit("  DISAGREEMENT IS EVIDENCE AGAINST NEITHER RESULT. Their measurement")
    emit("  is co-expression inside a 55 um Visium spot. Ours is a correlation")
    emit("  across 13 or 8 PATIENTS between two compartments, with no")
    emit("  coordinates at all. Neither design can falsify the other, and the")
    emit("  power floor here is 1.1-1.3 SD (P0-T8) — so 'not_recovered' means")
    emit("  UNINFORMATIVE, never contradicted.")
    emit()
    emit("  `not_on_panel`, `not_in_database` and `below_detection_floor` are")
    emit("  kept as separate statuses from `not_recovered` on purpose. A pair")
    emit("  we could not measure is NOT a pair we disagree about, and")
    emit("  collapsing the four is how an assay limit turns into a biological")
    emit("  claim.")
    emit()
    emit("  THE BRAIN ARM HAS NO COMPARATOR. De Zuani et al. profiled primary")
    emit("  lung; nothing here corroborates or contradicts a brain-metastasis")
    emit("  result, and TIME-B n = 8 regardless.")
    emit()

    counts = table["status"].value_counts().to_dict()
    emit(f"status counts across {len(table)} comparator x site x direction rows:")
    for k, v in sorted(counts.items()):
        emit(f"  {k:24s} {v}")
    emit()

    # ------------------------------------------------------------ assertions
    n_checks = 0
    expected_rows = len(comp["pairs"]) * len(adjacencies) * 2
    if len(table) != expected_rows:
        raise RuntimeError(
            f"{len(table)} rows, expected {expected_rows} "
            f"({len(comp['pairs'])} pairs x {len(adjacencies)} sites x 2 "
            "directions). Every comparator gets a row at every site in both "
            "directions, including the ones that cannot be evaluated — an "
            "absent row reads as a disagreement."
        )
    n_checks += 1

    if bool(table["is_validation"].any()):
        raise RuntimeError(
            "a row claims to be validation. PROJECT_PLAN §6 P4-T6 says informal "
            "concordance check, not validation, and the column exists to make "
            "that structural."
        )
    n_checks += 1

    # LGALS9 must be reported not_on_panel, and LGALS9C must not appear. This is
    # asserted rather than trusted because the substitution is exactly the
    # tempting fix, and it is a stop-and-ask.
    if "LGALS9" in set(table["ligand"]) | set(table["receptor"]):
        got = set(
            table.loc[
                (table["ligand"] == "LGALS9") | (table["receptor"] == "LGALS9"),
                "status",
            ]
        )
        if got != {"not_on_panel"}:
            raise RuntimeError(
                f"LGALS9 rows have status {sorted(got)}, expected only "
                "not_on_panel. LGALS9 is absent from the 18,694 measured "
                "symbols; if that changed, the matrix changed."
            )
    if "LGALS9C" in set(table["ligand"]) | set(table["receptor"]):
        raise RuntimeError(
            "LGALS9C appears in the concordance table. It is a PARALOG of the "
            "absent LGALS9 and substituting it is a stop-and-ask (ADR 0021 §8, "
            "ADR 0011, ADR 0015) — a membership change made to keep a "
            "comparator alive."
        )
    n_checks += 1
    emit(f"assertions: {n_checks} of {n_checks} passed "
         "(every comparator x site x direction present; nothing claims to be "
         "validation; LGALS9 not_on_panel and no paralog substituted)")
    emit()

    table.to_csv(snakemake.output.concordance, sep="\t", index=False)

    summary = {
        "task": "P4-T6",
        "aim": "A5",
        "exploratory": True,
        "exploratory_adr": "ADR 0014",
        "is_validation": False,
        "comparator": {
            "citation": comp["citation"],
            "source": comp["source"],
            "doi": comp["doi"],
            "assay": comp["assay"],
            "database_used": comp["database_used"],
            "tissue": comp["tissue"],
            "source_claim": " ".join(comp["source_claim"].split()),
            "n_pairs": len(comp["pairs"]),
        },
        "three_differences": [
            "They used CellPhoneDB; this project uses CellChatDB. A pair "
            "missing from our panel is a DATABASE difference, never a "
            "biological one.",
            "They measured co-expression within a Visium spot, tumour vs "
            "background section. This project measures a cross-patient "
            "correlation between two compartments with NO coordinates. "
            "Different quantities: agreement is encouraging, disagreement is "
            "evidence against neither.",
            "They profiled primary lung (LUAD/LUSC). The brain arm has NO "
            "comparator, and TIME-B n = 8.",
        ],
        "status_counts": counts,
        "concordance_counts": table["concordance"].value_counts().to_dict(),
        "status_meanings": {
            "not_on_panel": "the symbol is absent from the 18,694 measured genes",
            "not_in_database": "the pair is not an interaction in CellChatDB v2",
            "below_detection_floor": (
                "both partners are on the panel but at least one is below "
                "background in its compartment — not measurable, NOT a "
                "disagreement"
            ),
            "evaluable": "both partners above the floor; a rho and an FDR exist",
        },
        "reading_rule": (
            "`not_recovered` means UNINFORMATIVE, never contradicted: the power "
            "floor is 1.1-1.3 SD (P0-T8) and the two designs measure different "
            "quantities. `not_on_panel`, `not_in_database` and "
            "`below_detection_floor` are kept separate from `not_recovered` so "
            "an assay limit cannot become a biological claim."
        ),
        "lgals9_rule": (
            "LGALS9 is absent from the panel and is reported not_on_panel. "
            "LGALS9C is NOT substituted — adding a gene to an analysis panel is "
            "a stop-and-ask (ADR 0021 §8, ADR 0011, ADR 0015)."
        ),
        "plan_defect": (
            "PROJECT_PLAN §6 P4-T6 names four comparator pairs; the source "
            "reports FIVE — NRP1-VEGFA, NECTIN2-TIGIT, LGALS9-HAVCR2 and "
            "CD96-NECTIN1 enriched, PD1-PDL1 not. CD96-NECTIN1 is included "
            "here because a comparator list that drops one of the source's own "
            "results is a comparison chosen after the fact. Recorded in the "
            "Gate 4 ADR, since the plan is gitignored."
        ),
        "n_assertions": n_checks,
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit(f"wrote {snakemake.output.concordance}, {snakemake.output.summary}")
