"""Matched closed-cohort ODE and stochastic count microsimulation.

This module is a clean-room phase-1/2 bridge. The stochastic implementation
draws binomial individual events by state; the ODE implementation propagates
their expectation. Both use the same monthly hazards and accounting inputs.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from calendar_pde_cleanroom import ModelParameters, age_weight, discounted_life_years_lost, discount_factor


FLOW_COLUMNS = ("meningitis", "pneumonia_severe", "pneumonia_outpatient", "sepsis", "aom", "deaths")


def _strategy(vaccine: Optional[str], aom_mode: str, p: ModelParameters) -> tuple[float, float, float, float]:
    """Return coverage, IPD VE, AOM baseline-rate factor, and AOM VE."""
    if vaccine is None:
        return 0.0, 0.0, 1.0 if aom_mode == "reference_aligned" else p.spn_aom_proportion, 0.0
    ipd_irr = p.mature_population_irr[vaccine]["ipd"]
    ipd_ve = (1.0 - ipd_irr) / p.coverage
    if aom_mode == "reference_aligned":
        target_aom_irr = 1.0 - p.spn_aom_proportion * (1.0 - ipd_irr)
        return p.coverage, ipd_ve, 1.0, (1.0 - target_aom_irr) / p.coverage
    if aom_mode == "local_spn_only":
        return p.coverage, ipd_ve, p.spn_aom_proportion, ipd_ve
    raise ValueError(f"Unknown AOM mode: {aom_mode}")


def _monthly_probabilities(month: int, aom_rate_factor: float, p: ModelParameters) -> dict[str, float]:
    return {
        "meningitis": p.annual_incidence_per_100k["meningitis"] * 5.0 * age_weight(month, "meningitis") / 100_000.0,
        "pneumonia": p.annual_incidence_per_100k["pneumonia"] * 5.0 * age_weight(month, "pneumonia") / 100_000.0,
        "sepsis": p.annual_incidence_per_100k["sepsis"] * 5.0 * age_weight(month, "pneumonia") / 100_000.0,
        "aom": p.annual_incidence_per_100k["aom"] * aom_rate_factor * 5.0 * age_weight(month, "aom") / 100_000.0,
    }


def _economic_summary(monthly: pd.DataFrame, vaccine: Optional[str], p: ModelParameters) -> dict[str, float]:
    total_cost = total_qaly = 0.0
    for _, row in monthly.iterrows():
        age = float(row["month"] - 0.5) / 12.0
        disc = discount_factor(age, p.discount_rate)
        acute_cost = (
            row["pneumonia_severe"] * p.disease_costs_2024["pneumonia_inpatient"]
            + row["pneumonia_outpatient"] * p.disease_costs_2024["pneumonia_outpatient"]
            + row["meningitis"] * p.disease_costs_2024["meningitis"]
            + row["sepsis"] * p.disease_costs_2024["sepsis"]
            + row["aom"] * p.disease_costs_2024["aom"]
        )
        acute_qaly = (
            row["pneumonia_severe"] * p.acute_los_days["pneumonia_inpatient"] / 365.0 * (1.0 - p.acute_utilities["pneumonia_inpatient"])
            + row["pneumonia_outpatient"] * p.acute_los_days["pneumonia_outpatient"] / 365.0 * (1.0 - p.acute_utilities["pneumonia_outpatient"])
            + row["meningitis"] * p.acute_los_days["meningitis"] / 365.0 * (1.0 - p.acute_utilities["meningitis"])
            + row["sepsis"] * p.acute_los_days["sepsis"] / 365.0 * (1.0 - p.acute_utilities["sepsis"])
            + row["aom"] * p.acute_los_days["aom"] / 365.0 * (1.0 - p.acute_utilities["aom"])
        )
        sequela_cost = (
            row["seq_cognitive"] * p.sequela_costs_2024["cognitive_annual"] / 12.0
            + row["seq_special_education"] * p.sequela_costs_2024["special_education_annual"] / 12.0
            + row["new_meningitis_survivors"] * (
                p.sequela_probabilities["hearing"] * p.sequela_costs_2024["hearing"]
                + p.sequela_probabilities["epilepsy"] * p.sequela_costs_2024["epilepsy"]
                + p.sequela_probabilities["hemiplegia"] * p.sequela_costs_2024["hemiplegia"]
                + p.sequela_probabilities["hearing"] * p.sequela_probabilities["implant"] * p.sequela_costs_2024["implant"]
            )
        )
        sequela_qaly = (
            row["seq_cognitive"] / 12.0 * (1.0 - p.sequela_utilities["cognitive"])
            + row["seq_hearing"] / 12.0 * (1.0 - p.sequela_utilities["hearing"])
            + row["seq_epilepsy"] / 12.0 * (1.0 - p.sequela_utilities["epilepsy"])
            + row["seq_hemiplegia"] / 12.0 * (1.0 - p.sequela_utilities["hemiplegia"])
        )
        death_qaly = row["deaths"] * discounted_life_years_lost(age, p.life_expectancy_years, p.discount_rate)
        total_cost += (acute_cost + sequela_cost + row["deaths"] * p.pv_productivity_loss_per_death) * disc
        total_qaly += (acute_qaly + sequela_qaly) * disc + death_qaly
    if vaccine is not None:
        doses = p.annual_births * p.coverage * 3.0 * (1.0 + p.wastage_rate)
        program_cost = doses * (p.vaccine_price_2024[vaccine] + p.delivery_cost_per_dose)
        program_cost += doses * p.aefi_per_100k_doses / 100_000.0 * p.aefi_cost_2024
        total_cost += program_cost * discount_factor(0.25, p.discount_rate)
    return {"pv_total_cost": total_cost, "pv_total_qaly_loss": total_qaly}


def _make_monthly_row(month: int, alive_v: float, alive_u: float, events: dict[str, float], p: ModelParameters, sequelae: dict[str, float]) -> dict[str, float]:
    pneu_severe = events["pneumonia"] * p.severe_pneumonia_proportion
    deaths = events["meningitis"] * p.cfr["meningitis"] + pneu_severe * p.cfr["pneumonia"] + events["sepsis"] * p.cfr["sepsis"]
    survivor_men = events["meningitis"] * (1.0 - p.cfr["meningitis"])
    return {
        "month": month + 1, "alive": alive_v + alive_u,
        "meningitis": events["meningitis"], "pneumonia_severe": pneu_severe,
        "pneumonia_outpatient": events["pneumonia"] - pneu_severe,
        "sepsis": events["sepsis"], "aom": events["aom"], "deaths": deaths,
        "new_meningitis_survivors": survivor_men,
        "seq_cognitive": sequelae["cognitive"], "seq_hearing": sequelae["hearing"],
        "seq_epilepsy": sequelae["epilepsy"], "seq_hemiplegia": sequelae["hemiplegia"],
        "seq_special_education": sequelae["special_education"],
    }


def simulate_ode(vaccine: Optional[str], aom_mode: str, p: ModelParameters = ModelParameters()) -> pd.DataFrame:
    coverage, ipd_ve, aom_factor, aom_ve = _strategy(vaccine, aom_mode, p)
    alive_v = p.annual_births * coverage * (1.0 - p.neonatal_mortality_per_1000 / 1000.0)
    alive_u = p.annual_births * (1.0 - coverage) * (1.0 - p.neonatal_mortality_per_1000 / 1000.0)
    bg = 1.0 - (1.0 - p.under5_mortality_per_1000 / 1000.0) ** (1.0 / 60.0)
    seq = {key: 0.0 for key in ("cognitive", "hearing", "epilepsy", "hemiplegia", "special_education")}
    rows = []
    for month in range(60):
        alive_v *= 1.0 - bg
        alive_u *= 1.0 - bg
        probs = _monthly_probabilities(month, aom_factor, p)
        events = {
            "meningitis": alive_v * probs["meningitis"] * (1.0 - ipd_ve) + alive_u * probs["meningitis"],
            "pneumonia": alive_v * probs["pneumonia"] * (1.0 - ipd_ve) + alive_u * probs["pneumonia"],
            "sepsis": alive_v * probs["sepsis"] * (1.0 - ipd_ve) + alive_u * probs["sepsis"],
            "aom": alive_v * probs["aom"] * (1.0 - aom_ve) + alive_u * probs["aom"],
        }
        deaths = events["meningitis"] * p.cfr["meningitis"] + events["pneumonia"] * p.severe_pneumonia_proportion * p.cfr["pneumonia"] + events["sepsis"] * p.cfr["sepsis"]
        total_alive = alive_v + alive_u
        if total_alive:
            alive_v -= deaths * alive_v / total_alive
            alive_u -= deaths * alive_u / total_alive
        for key in seq:
            seq[key] *= 1.0 - bg
        men_survivors = events["meningitis"] * (1.0 - p.cfr["meningitis"])
        seq["cognitive"] += men_survivors * p.sequela_probabilities["cognitive"]
        seq["hearing"] += men_survivors * p.sequela_probabilities["hearing"]
        seq["epilepsy"] += men_survivors * p.sequela_probabilities["epilepsy"]
        seq["hemiplegia"] += men_survivors * p.sequela_probabilities["hemiplegia"]
        seq["special_education"] += men_survivors * p.sequela_probabilities["hearing"] * p.sequela_probabilities["special_education"]
        rows.append(_make_monthly_row(month, alive_v, alive_u, events, p, seq))
    return pd.DataFrame(rows)


def simulate_microsim(vaccine: Optional[str], aom_mode: str, rng: np.random.Generator, p: ModelParameters = ModelParameters()) -> pd.DataFrame:
    coverage, ipd_ve, aom_factor, aom_ve = _strategy(vaccine, aom_mode, p)
    cohort = int(round(p.annual_births))
    n_v = int(rng.binomial(cohort, coverage))
    n_u = cohort - n_v
    neonatal = p.neonatal_mortality_per_1000 / 1000.0
    n_v -= int(rng.binomial(n_v, neonatal))
    n_u -= int(rng.binomial(n_u, neonatal))
    bg = 1.0 - (1.0 - p.under5_mortality_per_1000 / 1000.0) ** (1.0 / 60.0)
    seq = {key: 0.0 for key in ("cognitive", "hearing", "epilepsy", "hemiplegia", "special_education")}
    rows = []
    for month in range(60):
        n_v -= int(rng.binomial(n_v, bg))
        n_u -= int(rng.binomial(n_u, bg))
        probs = _monthly_probabilities(month, aom_factor, p)
        events = {}
        for endpoint, ve in [("meningitis", ipd_ve), ("pneumonia", ipd_ve), ("sepsis", ipd_ve), ("aom", aom_ve)]:
            n_event_v = int(rng.binomial(n_v, min(1.0, probs[endpoint] * (1.0 - ve))))
            n_event_u = int(rng.binomial(n_u, min(1.0, probs[endpoint])))
            events[endpoint] = float(n_event_v + n_event_u)
        death_men = int(rng.binomial(int(events["meningitis"]), p.cfr["meningitis"]))
        death_pneu = int(rng.binomial(int(round(events["pneumonia"] * p.severe_pneumonia_proportion)), p.cfr["pneumonia"]))
        death_sep = int(rng.binomial(int(events["sepsis"]), p.cfr["sepsis"]))
        deaths = death_men + death_pneu + death_sep
        total_alive = n_v + n_u
        if total_alive:
            deaths_v = int(rng.binomial(deaths, n_v / total_alive))
            n_v = max(0, n_v - deaths_v)
            n_u = max(0, n_u - (deaths - deaths_v))
        for key in seq:
            seq[key] *= 1.0 - bg
        men_survivors = events["meningitis"] - death_men
        seq["cognitive"] += rng.binomial(int(men_survivors), p.sequela_probabilities["cognitive"])
        seq["hearing"] += rng.binomial(int(men_survivors), p.sequela_probabilities["hearing"])
        seq["epilepsy"] += rng.binomial(int(men_survivors), p.sequela_probabilities["epilepsy"])
        seq["hemiplegia"] += rng.binomial(int(men_survivors), p.sequela_probabilities["hemiplegia"])
        seq["special_education"] += rng.binomial(int(men_survivors), p.sequela_probabilities["hearing"] * p.sequela_probabilities["special_education"])
        rows.append(_make_monthly_row(month, float(n_v), float(n_u), events, p, seq))
    return pd.DataFrame(rows)


def summarize(monthly: pd.DataFrame, vaccine: Optional[str], p: ModelParameters) -> dict[str, float]:
    economics = _economic_summary(monthly, vaccine, p)
    return {
        "total_cases": float(monthly[["meningitis", "pneumonia_severe", "pneumonia_outpatient", "sepsis", "aom"]].sum().sum()),
        "total_deaths": float(monthly["deaths"].sum()),
        **economics,
    }
