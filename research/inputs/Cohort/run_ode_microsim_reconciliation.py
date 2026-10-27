"""Run clean-room ODE versus stochastic microsimulation reconciliation."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from closed_cohort_ode_microsim import simulate_microsim, simulate_ode, summarize
from calendar_pde_cleanroom import ModelParameters


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260721)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    monthly_dir = output / "source_data" / "monthly"
    monthly_dir.mkdir(parents=True, exist_ok=True)
    p = ModelParameters()
    rows = []
    replicate_rows = []
    for aom_mode in ("local_spn_only", "reference_aligned"):
        for scenario, vaccine in (("no_vaccination", None), ("mature_pcv13", "PCV13")):
            ode_monthly = simulate_ode(vaccine, aom_mode, p)
            ode_monthly.to_csv(monthly_dir / f"ode_{aom_mode}_{scenario}.csv", index=False)
            ode_summary = summarize(ode_monthly, vaccine, p)
            micro_summaries = []
            for replicate in range(args.replicates):
                rng = np.random.default_rng(args.seed + 10_000 * (0 if aom_mode == "local_spn_only" else 1) + 1_000 * (0 if vaccine is None else 1) + replicate)
                micro = summarize(simulate_microsim(vaccine, aom_mode, rng, p), vaccine, p)
                micro["replicate"] = replicate + 1
                micro_summaries.append(micro)
                replicate_rows.append({"aom_mode": aom_mode, "scenario": scenario, **micro})
            frame = pd.DataFrame(micro_summaries)
            for metric, ode_value in ode_summary.items():
                mean = float(frame[metric].mean())
                sd = float(frame[metric].std(ddof=1))
                se = sd / np.sqrt(args.replicates)
                rows.append({
                    "aom_mode": aom_mode, "scenario": scenario, "metric": metric,
                    "ode_value": ode_value, "microsim_mean": mean, "microsim_sd": sd,
                    "microsim_95ci_low": mean - 1.96 * se, "microsim_95ci_high": mean + 1.96 * se,
                    "absolute_difference_ode_minus_mean": ode_value - mean,
                    "z_difference": (ode_value - mean) / se if se else np.nan,
                    "ode_within_microsim_95ci": bool((mean - 1.96 * se) <= ode_value <= (mean + 1.96 * se)) if se else True,
                })
    reconciliation = pd.DataFrame(rows)
    reconciliation.to_csv(output / "ode_microsim_reconciliation.csv", index=False)
    pd.DataFrame(replicate_rows).to_csv(output / "microsim_replicate_summaries.csv", index=False)
    # Build intervention effects from the matched summary rows; pairs are independent
    # stochastic replications, so the uncertainty is descriptive rather than paired.
    effects = []
    for aom_mode in ("local_spn_only", "reference_aligned"):
        subset = reconciliation[reconciliation["aom_mode"] == aom_mode]
        for metric, direction in [("total_cases", -1.0), ("total_deaths", -1.0), ("pv_total_qaly_loss", -1.0), ("pv_total_cost", 1.0)]:
            base = subset[(subset["scenario"] == "no_vaccination") & (subset["metric"] == metric)].iloc[0]
            intervention = subset[(subset["scenario"] == "mature_pcv13") & (subset["metric"] == metric)].iloc[0]
            ode_effect = direction * (float(intervention["ode_value"]) - float(base["ode_value"]))
            mean_effect = direction * (float(intervention["microsim_mean"]) - float(base["microsim_mean"]))
            effects.append({"aom_mode": aom_mode, "outcome": metric, "ode_effect": ode_effect, "microsim_mean_effect": mean_effect, "difference": ode_effect - mean_effect})
    pd.DataFrame(effects).to_csv(output / "ode_microsim_intervention_effects.csv", index=False)
    print(f"Wrote ODE/microsimulation reconciliation for {args.replicates} replicates to {output}")
    print(f"ODE metrics inside microsimulation 95% CI: {int(reconciliation['ode_within_microsim_95ci'].sum())}/{len(reconciliation)}")


if __name__ == "__main__":
    main()
