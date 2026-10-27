"""Validate the PDE transport model in its closed-cohort ODE limiting case.

The scalar ODE characteristic and the age-bin transport implementation are
independently written below. Both use one birth cohort, no continuing births,
no dose transitions and no herd build-up. Equality verifies the algebraic
closed-cohort limit of the PDE transport layer; it is not external validation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from calendar_pde_cleanroom import AGE_MONTHS, ModelParameters, age_weight


def risk(age: int, endpoint: str, p: ModelParameters, multiplier: float) -> float:
    return p.annual_incidence_per_100k[endpoint] * 5.0 * age_weight(age, endpoint) / 100_000.0 * multiplier


def scalar_ode_characteristic(p: ModelParameters, multiplier: float) -> pd.DataFrame:
    alive = p.annual_births * (1.0 - p.neonatal_mortality_per_1000 / 1000.0)
    bg = 1.0 - (1.0 - p.under5_mortality_per_1000 / 1000.0) ** (1.0 / 60.0)
    rows = []
    for age in range(AGE_MONTHS):
        events = {endpoint: alive * risk(age, endpoint, p, multiplier) for endpoint in ("meningitis", "pneumonia", "sepsis", "aom")}
        deaths = events["meningitis"] * p.cfr["meningitis"] + events["pneumonia"] * p.severe_pneumonia_proportion * p.cfr["pneumonia"] + events["sepsis"] * p.cfr["sepsis"]
        alive -= alive * bg + deaths
        rows.append({"month": age + 1, "alive": alive, **events, "deaths": deaths})
    return pd.DataFrame(rows)


def age_time_pde_transport(p: ModelParameters, multiplier: float) -> pd.DataFrame:
    alive = np.zeros(AGE_MONTHS, dtype=float)
    alive[0] = p.annual_births * (1.0 - p.neonatal_mortality_per_1000 / 1000.0)
    bg = 1.0 - (1.0 - p.under5_mortality_per_1000 / 1000.0) ** (1.0 / 60.0)
    rows = []
    for month in range(AGE_MONTHS):
        age = month
        n = alive[age]
        events = {endpoint: n * risk(age, endpoint, p, multiplier) for endpoint in ("meningitis", "pneumonia", "sepsis", "aom")}
        deaths = events["meningitis"] * p.cfr["meningitis"] + events["pneumonia"] * p.severe_pneumonia_proportion * p.cfr["pneumonia"] + events["sepsis"] * p.cfr["sepsis"]
        alive[age] = n - n * bg - deaths
        rows.append({"month": month + 1, "alive": alive[age], **events, "deaths": deaths})
        if age < AGE_MONTHS - 1:
            alive[age + 1] = alive[age]
            alive[age] = 0.0
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    p = ModelParameters()
    rows = []
    for scenario, multiplier in [("baseline", 1.0), ("mature_pcv13_population_irr", p.mature_population_irr["PCV13"]["ipd"])]:
        ode = scalar_ode_characteristic(p, multiplier)
        pde = age_time_pde_transport(p, multiplier)
        for column in ["alive", "meningitis", "pneumonia", "sepsis", "aom", "deaths"]:
            rows.append({"scenario": scenario, "metric": column, "max_absolute_difference": float((ode[column] - pde[column]).abs().max()), "sum_absolute_difference": float((ode[column] - pde[column]).abs().sum())})
    result = pd.DataFrame(rows)
    result["pass_tolerance_1e-8"] = result["max_absolute_difference"] < 1e-8
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"ODE-limit checks passed: {int(result['pass_tolerance_1e-8'].sum())}/{len(result)}")


if __name__ == "__main__":
    main()
