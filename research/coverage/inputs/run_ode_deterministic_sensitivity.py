"""Deterministic one-way sensitivity analysis for the matched closed-cohort ODE."""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import pandas as pd

from calendar_pde_cleanroom import ModelParameters
from closed_cohort_ode_microsim import simulate_ode, summarize


def scaled(mapping: dict, factor: float) -> dict:
    return {key: float(value) * factor for key, value in mapping.items()}


def utility_scaled(mapping: dict, factor: float) -> dict:
    return {key: 1.0 - (1.0 - float(value)) * factor for key, value in mapping.items()}


def pair(p: ModelParameters) -> dict[str, float]:
    base = summarize(simulate_ode(None, "reference_aligned", p), None, p)
    pcv = summarize(simulate_ode("PCV13", "reference_aligned", p), "PCV13", p)
    qaly = base["pv_total_qaly_loss"] - pcv["pv_total_qaly_loss"]
    cost = pcv["pv_total_cost"] - base["pv_total_cost"]
    return {
        "cases_averted": base["total_cases"] - pcv["total_cases"],
        "deaths_averted": base["total_deaths"] - pcv["total_deaths"],
        "qaly_gained": qaly,
        "incremental_cost": cost,
        "icer": cost / qaly if qaly else float("nan"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = ModelParameters()
    variations = [("base_case", "base", "Base reference-aligned deterministic ODE", base)]
    for level, factor in [("low_75pct", 0.75), ("high_125pct", 1.25)]:
        variations += [
            ("vaccine_price", level, "Documented price-range proxy", replace(base, vaccine_price_2024={**base.vaccine_price_2024, "PCV13": base.vaccine_price_2024["PCV13"] * factor})),
            ("mature_population_effect", level, "IRR/VE uncertainty proxy", replace(base, mature_population_irr={k: scaled(dict(v), factor) for k, v in base.mature_population_irr.items()})),
            ("CFR", level, "Published-range proxy", replace(base, cfr=scaled(dict(base.cfr), factor))),
            ("sequela_cost", level, "Published-range proxy", replace(base, sequela_costs_2024=scaled(dict(base.sequela_costs_2024), factor))),
            ("utility_loss", level, "Utility-range proxy", replace(base, acute_utilities=utility_scaled(dict(base.acute_utilities), factor), sequela_utilities=utility_scaled(dict(base.sequela_utilities), factor))),
            ("disease_burden", level, "Incidence-range proxy", replace(base, annual_incidence_per_100k=scaled(dict(base.annual_incidence_per_100k), factor))),
        ]
    variations += [
        ("coverage", "80pct", "Coverage scenario", replace(base, coverage=0.80)),
        ("coverage", "95pct", "Coverage scenario", replace(base, coverage=0.95)),
        ("discount_rate", "0pct", "Guideline scenario", replace(base, discount_rate=0.0)),
        ("discount_rate", "5pct", "Guideline scenario", replace(base, discount_rate=0.05)),
    ]
    rows = []
    for parameter, level, source, p in variations:
        rows.append({"parameter": parameter, "level": level, "source_or_assumption": source, **pair(p)})
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)
    print(f"Wrote {len(rows)} ODE sensitivity rows to {output}")


if __name__ == "__main__":
    main()
