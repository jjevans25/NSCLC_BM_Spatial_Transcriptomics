"""P0-T8 — power reality-check for the TIME-L vs TIME-B contrast.

Owner task: P0-T8. Driven by rule p0t8_power_check.

§6 asks: given `TIME-L` 15 / `TIME-B` 8 with a patient random effect, what
standardised effect is detectable at 80% power? Simulation-based, not a
closed-form t-test — because the closed form has no way to represent the two
features that actually matter here.

The simulation uses the **real** design read from `samples.tsv`, not an
idealised balanced one:

  * `TIME-L` 15 AOIs across 13 patients (P12 and P24 contribute two each)
  * `TIME-B` 8 AOIs across 8 patients
  * **5 patients appear in both groups** (P5, P12, P15, P19, P35)
  * 16 distinct patients in total

Those five shared patients induce a within-subject correlation across the very
groups being compared, and the two replicate patients mean 23 AOIs are not 23
independent observations. A two-sample t-test on 15 vs 8 would assume both away
and report a more optimistic number than the design can support.

Effects are standardised: patient and residual variances sum to 1, so `effect`
reads as a between-group difference in SD units and ICC is the share of
variance sitting between patients.

**Inference.** statsmodels' MixedLM reports a Wald z. At 16 patients that is
anti-conservative, so significance is judged against a t distribution with
`n_patients - 2` degrees of freedom — a conservative stand-in for the
Satterthwaite approximation `lmerTest` will use in Phase 2. Both are reported,
so the gap between them is visible rather than assumed away.

Phase 2 will fit these models in `lme4` (CLAUDE.md's R/Python boundary; ADR
0006 moved only Phase 0 QC to Python). REML and the df approximation differ
slightly between the two, so treat these as close estimates, not guarantees.
"""

import json
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ALPHA = 0.05
TARGET_POWER = 0.80


def run_cell(args):
    """One (icc, effect) grid cell. Seeded from its own index so the result is
    independent of how the pool schedules it."""
    import statsmodels.formula.api as smf

    cell_index, seed, icc, effect, n_sim, patient_index, is_brain, n_patients = args
    rng = np.random.default_rng([seed, cell_index])
    sd_u = np.sqrt(icc)
    sd_e = np.sqrt(1.0 - icc)
    df_resid = n_patients - 2

    frame = pd.DataFrame({"patient": patient_index, "brain": is_brain.astype(float)})
    hits_t = hits_z = converged = 0

    for _ in range(n_sim):
        u = rng.normal(0.0, sd_u, n_patients) if sd_u > 0 else np.zeros(n_patients)
        frame["y"] = effect * is_brain + u[patient_index] + rng.normal(0.0, sd_e, len(is_brain))
        try:
            fit = smf.mixedlm("y ~ brain", frame, groups=frame["patient"]).fit(reml=True)
            se = float(fit.bse["brain"])
            if not np.isfinite(se) or se <= 0:
                continue
            converged += 1
            t = float(fit.params["brain"]) / se
        except Exception:
            continue
        if 2 * (1 - stats.t.cdf(abs(t), df_resid)) < ALPHA:
            hits_t += 1
        if 2 * (1 - stats.norm.cdf(abs(t))) < ALPHA:
            hits_z += 1

    return {
        "icc": icc,
        "effect": effect,
        "n_sim": n_sim,
        "n_converged": converged,
        "power_t": hits_t / converged if converged else np.nan,
        "power_z": hits_z / converged if converged else np.nan,
    }


def crossing(effects, powers, target=TARGET_POWER):
    """Linear interpolation of the effect size where power first reaches target."""
    for i in range(1, len(effects)):
        lo, hi = powers[i - 1], powers[i]
        if lo < target <= hi:
            if hi == lo:
                return float(effects[i])
            frac = (target - lo) / (hi - lo)
            return float(effects[i - 1] + frac * (effects[i] - effects[i - 1]))
    return float("nan")


if __name__ == "__main__":
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with open(log_path, "w", encoding="utf-8") as log:

        def emit(msg=""):
            print(msg, file=log, flush=True)

        seed = snakemake.params.seed
        n_sim = snakemake.params.n_sim
        emit(f"seed: {seed} (config['seed'], passed explicitly)")

        samples = pd.read_csv(snakemake.input.samples, sep="\t", dtype=str)
        immune = samples[samples["aoi_code"].isin(["TIME-L", "TIME-B"])]

        patients = sorted(immune["patient_id"].unique(), key=lambda p: int(p.lstrip("P")))
        lookup = {p: i for i, p in enumerate(patients)}
        patient_index = immune["patient_id"].map(lookup).to_numpy()
        is_brain = (immune["aoi_code"] == "TIME-B").to_numpy()
        n_patients = len(patients)

        n_l = int((~is_brain).sum())
        n_b = int(is_brain.sum())
        paired = sorted(
            set(immune.loc[immune["aoi_code"] == "TIME-L", "patient_id"])
            & set(immune.loc[immune["aoi_code"] == "TIME-B", "patient_id"]),
            key=lambda p: int(p.lstrip("P")),
        )
        emit(f"design read from samples.tsv: TIME-L n={n_l}, TIME-B n={n_b}, "
             f"{n_patients} patients, {len(paired)} in both ({', '.join(paired)})")
        emit(f"inference: t with df = {n_patients - 2}; Wald z reported alongside\n")

        iccs = [0.0, 0.25, 0.50, 0.75]
        effects = [round(0.2 * i, 2) for i in range(1, 11)]  # 0.2 .. 2.0
        jobs = [
            (i, seed, icc, effect, n_sim, patient_index, is_brain, n_patients)
            for i, (icc, effect) in enumerate(
                (icc, effect) for icc in iccs for effect in effects
            )
        ]
        emit(f"grid: {len(iccs)} ICC x {len(effects)} effects x {n_sim} sims "
             f"= {len(jobs) * n_sim:,} model fits on {snakemake.threads} thread(s)")

        with ProcessPoolExecutor(max_workers=snakemake.threads) as pool:
            rows = list(pool.map(run_cell, jobs))

        grid = pd.DataFrame(rows).sort_values(["icc", "effect"]).reset_index(drop=True)
        grid.to_csv(snakemake.output.power, sep="\t", index=False, float_format="%.6g")
        emit(f"\nwrote {snakemake.output.power}")

        mde = {}
        emit("\nminimum detectable standardised effect at 80% power:")
        for icc in iccs:
            block = grid[grid["icc"] == icc]
            value_t = crossing(block["effect"].to_numpy(), block["power_t"].to_numpy())
            value_z = crossing(block["effect"].to_numpy(), block["power_z"].to_numpy())
            mde[icc] = value_t
            emit(f"  ICC {icc:.2f}:  t-ref {value_t:.2f} SD   (Wald z would say {value_z:.2f})")

        converged_min = float(grid["n_converged"].min() / n_sim)
        emit(f"\nlowest convergence rate across cells: {converged_min:.1%}")

        Path(snakemake.output.summary).write_text(
            json.dumps(
                {
                    "seed": seed,
                    "n_sim_per_cell": n_sim,
                    "alpha": ALPHA,
                    "target_power": TARGET_POWER,
                    "n_time_l": n_l,
                    "n_time_b": n_b,
                    "n_patients": n_patients,
                    "n_paired_patients": len(paired),
                    "paired_patients": paired,
                    "df_used": n_patients - 2,
                    "mde_by_icc": {str(k): (None if np.isnan(v) else round(v, 3))
                                   for k, v in mde.items()},
                    "min_convergence_rate": round(converged_min, 4),
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        emit(f"wrote {snakemake.output.summary}")

        # ------------------------------------------------------------- figure
        fig, ax = plt.subplots(figsize=(8.6, 5.6))
        colours = ["#3E7CB1", "#5B8C5A", "#B0763F", "#8C3F5B"]
        for (icc, colour) in zip(iccs, colours):
            block = grid[grid["icc"] == icc]
            ax.plot(block["effect"], block["power_t"], marker="o", ms=4, lw=1.6,
                    color=colour, label=f"ICC {icc:.2f}  (MDE {mde[icc]:.2f} SD)")
        ax.axhline(TARGET_POWER, ls="--", lw=1.2, color="crimson")
        ax.text(0.22, TARGET_POWER + 0.015, "80% power", color="crimson", fontsize=9)
        ax.set_xlabel("standardised effect (SD units)", fontsize=10)
        ax.set_ylabel("power", fontsize=10)
        ax.set_ylim(0, 1.02)
        ax.set_xlim(0.15, 2.05)
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=9, loc="lower right")
        ax.set_title(
            f"P0-T8 — TIME-L (n={n_l}) vs TIME-B (n={n_b}), patient random intercept\n"
            f"{n_patients} patients, {len(paired)} contributing to both groups; "
            f"{n_sim} simulations per point",
            fontsize=11, fontweight="bold",
        )
        fig.tight_layout()
        fig.savefig(snakemake.output.figure, dpi=200, bbox_inches="tight")
        plt.close(fig)
        emit(f"wrote {snakemake.output.figure}")
