"""P4-T4 — empirical null calibration. THE TASK GATE 4 TURNS ON.

Owner task: P4-T4. Driven by rule p4t4_null_calibration.

PROJECT_PLAN §6 P4-T4, verbatim: "With ~1000 LR pairs at n=13, some |rho|>0.7
arises by chance. Permute patient labels to get an empirical null and report an
empirical FDR. **Without this, the ranking is uninterpretable.**"

Gate 4's pass condition is *this number being computed and honoured*. A
surviving list and an empty list both pass; a ranked list without a null does
not. So this script produces the FDR whatever it says, and nothing downstream is
allowed to reach past it for a softer threshold.

TWO NULLS, both pre-registered in ADR 0021 §6, and only the first is the gate.

**Null A — the gate.** Permute *the pairing*, not the expression: shuffle which
patient's immune AOI is matched to which patient's tumour AOI, **within site**,
so the null preserves both n's exactly (13 and **8**). Everything else — each
gene's values, its distribution, its abundance, its detection — is untouched.
What is destroyed is exactly the thing the analysis claims: that the ligand and
the receptor were measured in the same person.

The empirical FDR at a threshold t is the expected number of null rows with
|rho| >= t divided by the observed number, computed **within site across all
primary direction-rows** (ADR 0021 §6) and made monotone in |rho|.

**Null B — the abundance-matched control.** `docs/NEXT_STEPS.md`'s
reconnaissance warns that a list dominated by VEGFA, TGFB1 and CD47 may be
measuring "high-expression genes correlate with high-expression genes" — these
are broadly expressed, and across 13 patients their correlation may track shared
AOI quality rather than crosstalk. P4-T3a measures whether that quality is in
fact shared across a patient's paired AOIs; this answers the question directly
whatever that came out at, because a pre-registered control is not contingent on
the confound turning out to be large.

For each row: random gene pairs drawn from the same mean-expression decile as
the real tumour-side and immune-side partners, **with the true pairing intact**,
and the fraction reaching the observed |rho|. Drawn once per (decile, decile)
cell rather than once per row — the null depends only on the decile pair, so
100 cells serve every row and the arithmetic is identical.

Null B **qualifies** a nomination; it is **not** a second gate, it creates no
second FDR family, and a row that clears A but not B is reported as such rather
than removed.

Both nulls seed from `config["seed"]`, passed explicitly. No implicit RNG
anywhere (PROJECT_PLAN §4.3); `variance_partition.py` is the in-repo precedent
for the `default_rng(seed + i)` convention.

**Why Spearman makes this cheap.** Spearman is Pearson on ranks, and a
permutation only reorders values, so ranks are computed ONCE. With each
partner's rank vector z-scored, rho is a dot product and a whole permutation is
one elementwise multiply — which is what makes 10,000 permutations across
thousands of rows a few seconds rather than an afternoon.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Colour language, stated so it cannot drift. RdBu_r means a signed effect
# (P2-T7) and PRGn centred on the detection multiple means a background ratio
# (P3-T5, ADR 0016 §5). Neither is borrowed here: these are DISTRIBUTIONS, not
# a signed scale, so they get neutral fills and the observed/null contrast
# carries the meaning.
OBS_COLOUR = "#C1666B"
NULL_COLOUR = "#4F6D7A"

def _site_stream(name):
    """A stable integer stream id for a site name.

    NOT `hash(name)`. Python salts string hashing per process
    (PYTHONHASHSEED), so `hash("brain")` differs between runs and the
    permutation null would be reseeded differently every time — which is
    exactly the implicit RNG PROJECT_PLAN §4.3 forbids, wearing the costume of
    an explicit seed. It cost this phase a non-reproducible table: two runs of
    the same commit gave min empirical FDR 0.228 and 0.225.

    sha256 of the name is stable across processes, machines and Python
    versions, which is what "seeded from config['seed'], passed explicitly"
    has to mean if P6-T1's clean-room reproduction is to be byte-identical.
    """
    import hashlib

    return int.from_bytes(
        hashlib.sha256(name.encode("utf-8")).digest()[:4], "big"
    )


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    seed = snakemake.params.seed
    cfg = snakemake.params.null_calibration
    n_perm = cfg["n_permutations"]
    n_draws = cfg["n_matched_draws"]
    n_deciles = cfg["n_expression_deciles"]
    alpha = cfg["fdr_alpha"]
    adjacencies = snakemake.params.adjacencies

    emit("P4-T4 — empirical null calibration. THE TASK GATE 4 TURNS ON.")
    emit("=" * 70)
    emit()
    emit(f"  seed                {seed}  (config['seed'], passed explicitly)")
    emit(f"  Null A permutations {n_perm}  — permutes THE PAIRING, within site")
    emit(f"  Null B draws        {n_draws} per decile cell, {n_deciles} deciles")
    emit(f"  FDR alpha           {alpha}")
    emit("  family              all primary direction-rows, WITHIN SITE")
    emit("                      (ADR 0021 §6 — a null built at n = 13 cannot")
    emit("                       be applied to the brain's n = 8)")
    emit()
    emit("  Gate 4 passes if this FDR is computed and honoured. A surviving")
    emit("  list and an EMPTY list both pass; a ranked list without a null")
    emit("  does not. Nothing downstream may reach past it.")
    emit()

    # ---------------------------------------------------------------- inputs
    long = pd.read_csv(
        snakemake.input.expression, sep="\t", float_precision="round_trip"
    )
    corr = pd.read_csv(
        snakemake.input.correlation, sep="\t", float_precision="round_trip"
    )
    emit(f"primary direction-rows: {len(corr)}")

    side_of = {}
    for site in sorted(adjacencies):
        side_of[(site, adjacencies[site]["tumour"])] = "tumour"
        side_of[(site, adjacencies[site]["immune"])] = "immune"

    patients = {
        site: sorted(
            long.loc[long["site"] == site, "patient_id"].unique(),
            key=lambda p: int(str(p).lstrip("P")) if str(p)[1:].isdigit() else 0,
        )
        for site in sorted(adjacencies)
    }

    value = {}
    for (site, side, gene), grp in long.groupby(["site", "side", "gene"]):
        value[(site, side, gene)] = (
            grp.set_index("patient_id")["expression"]
            .reindex(patients[site])
            .to_numpy(dtype=float)
        )

    def zrank(v):
        """Rank, then z-score. Spearman rho is then a dot product over n.

        Ranks are invariant to permutation, so this is computed once per
        partner and reused across every permutation — which is the whole
        reason 10,000 permutations is cheap.
        """
        r = stats.rankdata(v)
        sd = r.std()
        if sd == 0:
            raise RuntimeError(
                "a partner vector is constant across patients, so its rank "
                "has zero variance and Spearman is undefined. p4t3b already "
                "asserts rho is finite; reaching here means the two rules read "
                "different expression tables."
            )
        return (r - r.mean()) / sd

    def partner_vec(site, aoi_code, subunits):
        side = side_of[(site, aoi_code)]
        vals = np.vstack([value[(site, side, g)] for g in subunits])
        return vals.mean(axis=0)

    # ------------------------------------------- per-site arrays, tumour/immune
    # Every row reduces to the same arithmetic regardless of direction: one
    # partner measured on the TUMOUR side against one measured on the IMMUNE
    # side. `direction` labels which of them is the ligand; it does not change
    # the correlation. So one permutation of the immune side serves both
    # directions, which is also why the null preserves the design exactly.
    results = []
    null_rows = []
    fig, axes = plt.subplots(
        1, len(adjacencies), figsize=(5.2 * len(adjacencies), 4.4)
    )
    axes = np.atleast_1d(axes)

    for ax, site in zip(axes, sorted(adjacencies)):
        sel = corr[corr["site"] == site].reset_index(drop=True)
        n = len(patients[site])
        emit(f"{site} — n = {n} patients, {len(sel)} primary direction-rows")

        tum_code = adjacencies[site]["tumour"]
        imm_code = adjacencies[site]["immune"]

        A = np.empty((len(sel), n))  # tumour-side z-ranks
        B = np.empty((len(sel), n))  # immune-side z-ranks
        mean_t = np.empty(len(sel))
        mean_i = np.empty(len(sel))
        for k, r in enumerate(sel.itertuples(index=False)):
            lig = [s for s in str(r.ligand_subunits).split(";") if s]
            rec = [s for s in str(r.receptor_subunits).split(";") if s]
            if side_of[(site, r.ligand_aoi_code)] == "tumour":
                tv, iv = partner_vec(site, r.ligand_aoi_code, lig), partner_vec(
                    site, r.receptor_aoi_code, rec
                )
            else:
                tv, iv = partner_vec(site, r.receptor_aoi_code, rec), partner_vec(
                    site, r.ligand_aoi_code, lig
                )
            A[k], B[k] = zrank(tv), zrank(iv)
            mean_t[k], mean_i[k] = tv.mean(), iv.mean()

        rho_obs = (A * B).sum(axis=1) / n

        # Identity-permutation check: the machinery must reproduce p4t3b's rho
        # exactly before any permuted value is trusted.
        delta = float(np.abs(rho_obs - sel["rho"].to_numpy()).max())
        if delta > 1e-9:
            raise RuntimeError(
                f"{site}: the null machinery disagrees with p4t3b's rho by "
                f"{delta:.3g} under the identity permutation. The permuted "
                "values are computed the same way, so they cannot be trusted "
                "either."
            )
        emit(f"  identity-permutation check vs p4t3b: max |delta| = {delta:.3g}")

        # ------------------------------------------------ Null A, the gate
        rng = np.random.default_rng([seed, _site_stream(site)])
        # float64 throughout: the null and the observed statistic are separate
        # computations, and at n = 8 rho is heavily tied (42 distinct values
        # across 260 rows), so a float32 null would shift tie counts by
        # rounding alone.
        null_abs = np.empty((n_perm, len(sel)), dtype=np.float64)
        for p_i in range(n_perm):
            perm = rng.permutation(n)
            null_abs[p_i] = np.abs((A * B[:, perm]).sum(axis=1) / n)
        flat_null = null_abs.ravel()
        emit(f"  Null A: {n_perm} permutations x {len(sel)} rows = "
             f"{flat_null.size} null values")

        obs_abs = np.abs(rho_obs)

        # --- empirical FDR, computed at DISTINCT thresholds -----------------
        # Spearman on 8 patients takes few distinct values, so |rho| is heavily
        # tied. Ranking rows 1..N and dividing by the running index would give
        # tied rows DIFFERENT denominators, which is simply the wrong count:
        # at a threshold t, the number observed at or above t is every row that
        # reaches it, not that row's position in an arbitrary tie order.
        #
        # So the FDR is computed once per distinct |rho| and mapped back by
        # value. That makes it a function of the statistic — which is what an
        # FDR is — and makes tied rows share a value by construction rather
        # than by a later pass that happens to equalise them.
        thresholds = np.unique(obs_abs)              # ascending, distinct
        sorted_null = np.sort(flat_null)
        # counts at or above each threshold
        n_obs_ge = len(obs_abs) - np.searchsorted(
            np.sort(obs_abs), thresholds, side="left"
        )
        n_null_ge = flat_null.size - np.searchsorted(
            sorted_null, thresholds, side="left"
        )
        expected = n_null_ge / n_perm
        raw = np.minimum(expected / n_obs_ge, 1.0)

        # q(t) = min over thresholds at or BELOW t, so the FDR cannot fall as
        # the threshold loosens. `thresholds` is ascending, so that is a
        # running minimum from the bottom up.
        fdr_at = np.minimum.accumulate(raw)
        fdr = fdr_at[np.searchsorted(thresholds, obs_abs)]

        # One quantity, not two. `sel["abs_rho"]` came from p4t3b's TSV and
        # this came from the rank arrays here; the identity check above shows
        # they agree to 1e-9, but sorting a table by one column while a
        # monotonicity claim was established on the other is how a spurious
        # assertion failure gets "fixed" by loosening the assertion.
        sel["abs_rho"] = obs_abs

        order = np.argsort(-obs_abs)
        sorted_abs = obs_abs[order]

        n_sig = int((fdr <= alpha).sum())
        emit(f"  empirical FDR <= {alpha}: {n_sig} of {len(sel)} rows")
        emit(f"  null |rho| quantiles: "
             f"50% {np.quantile(flat_null, 0.50):.3f}  "
             f"95% {np.quantile(flat_null, 0.95):.3f}  "
             f"99% {np.quantile(flat_null, 0.99):.3f}  "
             f"max {flat_null.max():.3f}")
        emit(f"  observed |rho|:       "
             f"50% {np.quantile(obs_abs, 0.50):.3f}  "
             f"95% {np.quantile(obs_abs, 0.95):.3f}  "
             f"99% {np.quantile(obs_abs, 0.99):.3f}  "
             f"max {obs_abs.max():.3f}")

        # ----------------------------------- Null B, the abundance control
        # Deciles over the genes reaching an admitted interaction in that
        # compartment — the same admission rule as the real partners, so the
        # comparison is like-for-like. The pool is the LR panel rather than the
        # transcriptome, deliberately: "among genes of this abundance THAT ARE
        # LIGANDS AND RECEPTORS" is the question a nomination has to beat.
        pools = {}
        for side in ("tumour", "immune"):
            genes = sorted(
                long.loc[
                    (long["site"] == site) & (long["side"] == side), "gene"
                ].unique()
            )
            mat = np.vstack([value[(site, side, g)] for g in genes])
            means = mat.mean(axis=1)
            edges = np.quantile(means, np.linspace(0, 1, n_deciles + 1))
            dec = np.clip(np.searchsorted(edges, means, side="right") - 1,
                          0, n_deciles - 1)
            zr = np.vstack([zrank(mat[j]) for j in range(len(genes))])
            pools[side] = {"z": zr, "decile": dec, "edges": edges}
        emit(f"  Null B: pools — tumour {pools['tumour']['z'].shape[0]} genes, "
             f"immune {pools['immune']['z'].shape[0]} genes, "
             f"{n_deciles} deciles")

        def decile_of(side, means):
            """Place a real partner in the pool's own decile bins.

            Reuses the edges computed for the pool above rather than
            recomputing them: two quantile calls that are meant to agree are
            two places for them to stop agreeing.
            """
            edges = pools[side]["edges"]
            return np.clip(
                np.searchsorted(edges, means, side="right") - 1, 0, n_deciles - 1
            )

        dec_t = decile_of("tumour", mean_t)
        dec_i = decile_of("immune", mean_i)

        # One null per (decile, decile) cell, shared by every row in it.
        cell_null = {}
        rng_b = np.random.default_rng([seed, 99, _site_stream(site)])
        for dt in range(n_deciles):
            it = np.flatnonzero(pools["tumour"]["decile"] == dt)
            if not len(it):
                continue
            for di in range(n_deciles):
                ii = np.flatnonzero(pools["immune"]["decile"] == di)
                if not len(ii):
                    continue
                a = pools["tumour"]["z"][rng_b.choice(it, n_draws)]
                b = pools["immune"]["z"][rng_b.choice(ii, n_draws)]
                cell_null[(dt, di)] = np.sort(np.abs((a * b).sum(axis=1) / n))

        p_matched = np.ones(len(sel))
        for k in range(len(sel)):
            nd = cell_null.get((int(dec_t[k]), int(dec_i[k])))
            if nd is None:
                p_matched[k] = float("nan")
                continue
            ge = nd.size - np.searchsorted(nd, obs_abs[k], side="left")
            p_matched[k] = (1 + ge) / (1 + nd.size)
        emit(f"  Null B: {len(cell_null)} decile cells populated; "
             f"abundance-matched p <= {alpha}: "
             f"{int((p_matched <= alpha).sum())} of {len(sel)}")

        both = int(((fdr <= alpha) & (p_matched <= alpha)).sum())
        emit(f"  clears Null A AND Null B: {both} of {len(sel)}")
        emit()

        out = sel.copy()
        out["empirical_fdr"] = fdr
        out["null_a_n_permutations"] = n_perm
        out["abundance_matched_p"] = p_matched
        out["abundance_decile_tumour"] = dec_t
        out["abundance_decile_immune"] = dec_i
        out["clears_empirical_fdr"] = fdr <= alpha
        out["clears_abundance_matched"] = p_matched <= alpha
        results.append(out)

        for q in (0.50, 0.90, 0.95, 0.99, 0.999, 1.0):
            null_rows.append(
                {
                    "site": site,
                    "n_patients": n,
                    "n_rows": len(sel),
                    "n_permutations": n_perm,
                    "quantile": q,
                    "null_abs_rho": float(np.quantile(flat_null, q)),
                    "observed_abs_rho": float(np.quantile(obs_abs, q)),
                }
            )

        # ------------------------------------------------------------ figure
        bins = np.linspace(0, 1, 51)
        ax.hist(flat_null, bins=bins, density=True, color=NULL_COLOUR,
                alpha=0.55, label=f"null ({n_perm} permutations)")
        ax.hist(obs_abs, bins=bins, density=True, histtype="step",
                color=OBS_COLOUR, lw=2.0, label=f"observed ({len(sel)} rows)")
        if n_sig:
            cut = float(sorted_abs[n_sig - 1])
            ax.axvline(cut, color="black", lw=1.1, ls="--")
            ax.text(cut, ax.get_ylim()[1] * 0.92, f"  FDR {alpha}\n  |rho| {cut:.2f}",
                    fontsize=8, va="top")
        ax.set_title(f"{site}   n = {n} patients", fontsize=10)
        ax.set_xlabel("|Spearman rho|", fontsize=9)
        ax.set_ylabel("density", fontsize=9)
        ax.legend(fontsize=7.5, frameon=False)
        ax.tick_params(labelsize=8)

    fdr_table = pd.concat(results, ignore_index=True)
    null_table = pd.DataFrame(null_rows)

    # Keep every title line SHORT — P3-T5's clipped at both canvas edges on
    # first render despite the plan warning about it. Verify by opening the PNG.
    fig.suptitle(
        "P4-T4  inferred crosstalk: observed vs. permuted-pairing null\n"
        "permutation breaks the patient pairing only; brain is TIME-B n = 8\n"
        "exploratory (ADR 0014) — a null here is an assay-sensitivity limit",
        fontsize=9.5, y=0.995, va="top",
    )
    # Explicit geometry, not tight_layout: with a three-line suptitle the
    # rect form left a sixteen-percent band of dead canvas between the
    # title and the panels. Verified by opening the PNG, which is the only
    # way this class of defect is ever caught (P3-T5's lesson).
    fig.subplots_adjust(top=0.76, bottom=0.14, left=0.07, right=0.98,
                        wspace=0.22)
    fig.savefig(snakemake.output.figure, dpi=200)
    plt.close(fig)
    emit(f"wrote {snakemake.output.figure}")

    # ------------------------------------------------------------ assertions
    n_checks = 0
    if len(fdr_table) != len(corr):
        raise RuntimeError(
            f"{len(fdr_table)} rows out, {len(corr)} in. The FDR attaches to "
            "every primary row; a row that fails must carry an FDR of 1, never "
            "be absent."
        )
    n_checks += 1
    if not ((fdr_table["empirical_fdr"] >= 0) & (fdr_table["empirical_fdr"] <= 1)).all():
        raise RuntimeError("empirical_fdr outside [0, 1].")
    n_checks += 1
    for site, grp in fdr_table.groupby("site"):
        g = grp.sort_values("abs_rho", ascending=False)
        if (g["empirical_fdr"].diff().dropna() < -1e-12).any():
            raise RuntimeError(
                f"{site}: empirical_fdr is not monotone in |rho|. An FDR that "
                "falls as the threshold loosens is a bug in the step-down."
            )
    n_checks += 1
    emit(f"assertions: {n_checks} of {n_checks} passed "
         "(every primary row carries an FDR; FDR in [0,1]; FDR monotone in |rho|)")
    emit()

    fdr_table.to_csv(snakemake.output.fdr, sep="\t", index=False)
    null_table.to_csv(snakemake.output.null, sep="\t", index=False)

    by_site = {}
    for site, grp in fdr_table.groupby("site"):
        by_site[site] = {
            "n_patients": len(patients[site]),
            "n_rows": int(len(grp)),
            "n_clears_empirical_fdr": int(grp["clears_empirical_fdr"].sum()),
            "n_clears_abundance_matched": int(grp["clears_abundance_matched"].sum()),
            "n_clears_both": int(
                (grp["clears_empirical_fdr"] & grp["clears_abundance_matched"]).sum()
            ),
            "min_empirical_fdr": float(grp["empirical_fdr"].min()),
            "max_abs_rho": float(grp["abs_rho"].max()),
        }
        emit(f"{site}: {by_site[site]}")

    summary = {
        "task": "P4-T4",
        "aim": "A5",
        "exploratory": True,
        "exploratory_adr": "ADR 0014",
        "pre_registered": "ADR 0021 §6",
        "gate": "This is the task Gate 4 turns on.",
        "gate_rule": (
            "Gate 4 passes if the empirical FDR is computed and honoured. A "
            "surviving list and an EMPTY list both pass; a ranked list without "
            "a null does not. Nothing downstream may reach past it for a "
            "softer threshold."
        ),
        "seed": seed,
        "null_a": {
            "what_is_permuted": (
                "the pairing — which patient's immune AOI is matched to which "
                "patient's tumour AOI — within site, so both n's are preserved "
                "exactly. Each gene's values, distribution, abundance and "
                "detection are untouched."
            ),
            "n_permutations": n_perm,
            "is_the_gate": True,
        },
        "null_b": {
            "what_is_drawn": (
                "random gene pairs from the same mean-expression decile as the "
                "real tumour-side and immune-side partners, with the TRUE "
                "pairing intact. Drawn once per (decile, decile) cell — the "
                "null depends only on the decile pair."
            ),
            "pool": (
                "genes reaching an admitted interaction in that compartment — "
                "the same admission rule as the real partners, so the "
                "comparison is like-for-like. 'Among ligands and receptors of "
                "this abundance' is the question a nomination has to beat."
            ),
            "n_matched_draws": n_draws,
            "n_expression_deciles": n_deciles,
            "is_the_gate": False,
            "role": (
                "qualifies a nomination. NOT a second gate, creates no second "
                "FDR family; a row that clears A but not B is reported as such, "
                "never removed."
            ),
        },
        "fdr_alpha": alpha,
        "family": (
            "all primary direction-rows WITHIN SITE, lung and brain separately "
            "(ADR 0021 §6). A null built at n = 13 cannot be applied to n = 8. "
            "Both directions of an interaction sit in the same family."
        ),
        "by_site": by_site,
        "n_time_b": len(patients.get("brain", [])),
        "power_floor_sd": "1.1-1.3",
        "reporting_rule": (
            "A null is an ASSAY-SENSITIVITY LIMIT, never evidence that the "
            "crosstalk is absent (ADR 0014). TIME-B n = 8 (hard constraint 8). "
            "'no LR pair exceeded chance expectation at n=13' is an honest "
            "result and passes this gate."
        ),
        "n_assertions": n_checks,
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit()
    emit(f"wrote {snakemake.output.fdr}, {snakemake.output.null}, "
         f"{snakemake.output.summary}")
