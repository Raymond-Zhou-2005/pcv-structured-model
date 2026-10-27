"""Clean-room age--calendar-time PCV model for transparent scenario analysis.

This monthly upwind transport model is a finite-difference representation of an
age--time PDE.  It was written in this workspace from the documented inputs and
assumptions, without modifying or importing the original project source files.
It is a research model, not a clinical or policy recommendation engine.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from itertools import product
from typing import Dict, Iterable, Mapping, MutableMapping, Optional

import numpy as np
import pandas as pd


GROUPS = ("unvaccinated", "enrolled_dose0", "dose1", "dose2", "booster")
G_UNVACCINATED, G_DOSE0, G_DOSE1, G_DOSE2, G_BOOSTER = range(5)
AGE_MONTHS = 60
ENDPOINTS = ("meningitis", "pneumonia", "sepsis", "aom")


@dataclass(frozen=True)
class ModelParameters:
    # Published reference-manuscript demographic and epidemiological inputs.
    annual_births: float = 1_384_610.0
    under5_mortality_per_1000: float = 20.0
    neonatal_mortality_per_1000: float = 10.3
    annual_incidence_per_100k: Mapping[str, float] = field(default_factory=lambda: {
        "pneumonia": 1771.0, "meningitis": 12.0, "sepsis": 124.0, "aom": 9570.0,
    })
    severe_pneumonia_proportion: float = 0.389
    cfr: Mapping[str, float] = field(default_factory=lambda: {
        "pneumonia": 0.014241, "meningitis": 0.032258, "sepsis": 0.08,
    })
    spn_aom_proportion: float = 0.303
    coverage: float = 0.90
    warmup_months: int = 60
    horizon_months: int = 120

    # Base direct VE profiles. PCV10-SII is explicitly a scenario assumption,
    # not an observed Vietnam effectiveness estimate.
    direct_ve: Mapping[str, Mapping[str, float]] = field(default_factory=lambda: {
        "PCV13": {"ipd": 0.602, "pneumonia": 0.19, "aom": 0.46},
        "PCV10-SII": {"ipd": 0.493, "pneumonia": 0.156, "aom": 0.377},
    })
    mature_population_irr: Mapping[str, Mapping[str, float]] = field(default_factory=lambda: {
        "PCV13": {"ipd": 0.2193, "pneumonia": 0.2193, "aom": 0.8660},
        "PCV10-SII": {"ipd": 0.2211, "pneumonia": 0.2211, "aom": 0.8731},
    })
    # Optional fixed H_e* values used only for direct-VE sensitivity analysis.
    # Holding the base indirect multiplier fixed prevents a direct-VE perturbation
    # from automatically rebalancing the inferred herd effect to the same mature
    # population IRR, which would obscure the isolated direct-VE contribution.
    mature_herd_multiplier_override: Optional[Mapping[str, Mapping[str, float]]] = None
    effective_dose_share_reference: float = 0.88
    carriage_reduction: Mapping[str, float] = field(default_factory=lambda: {
        "enrolled_dose0": 0.0, "dose1": 0.15, "dose2": 0.34, "booster": 0.40,
    })
    carriage_age_weights: Mapping[str, float] = field(default_factory=lambda: {
        "0_5": 0.0, "6_17": 1.0, "18_23": 0.8, "24_59": 0.35,
    })
    max_direct_ve: float = 0.95
    waning_severity: float = 1.0  # 0=no waning; 1=base structural waning curve
    reaction_substeps: int = 1  # numerical integration sensitivity control

    # Cost/QALY inputs transcribed from external-manuscript Table 1. Values are
    # used as 2024 US dollars; this removes the local historical 2023 alignment.
    disease_costs_2024: Mapping[str, float] = field(default_factory=lambda: {
        "pneumonia_inpatient": 221.09 + 142.33 + 184.25,
        "pneumonia_outpatient": 16.09 + 21.21 + 139.12,
        "meningitis": 574.50 + 297.95 + 297.82,
        "sepsis": 589.60 + 305.24 + 193.55,
        "aom": 13.31 + 21.21 + 139.12,
    })
    acute_utilities: Mapping[str, float] = field(default_factory=lambda: {
        "pneumonia_inpatient": 0.9941, "pneumonia_outpatient": 0.9963,
        "meningitis": 0.9768, "sepsis": 0.9921, "aom": 0.9950,
    })
    acute_los_days: Mapping[str, float] = field(default_factory=lambda: {
        "pneumonia_inpatient": 8.74, "pneumonia_outpatient": 7.0,
        "meningitis": 18.41, "sepsis": 12.44, "aom": 7.0,
    })
    sequela_probabilities: Mapping[str, float] = field(default_factory=lambda: {
        "cognitive": 0.016, "hearing": 0.046, "epilepsy": 0.022,
        "hemiplegia": 0.032, "implant": 0.025, "special_education": 0.2574,
    })
    sequela_utilities: Mapping[str, float] = field(default_factory=lambda: {
        "cognitive": 0.62, "hearing": 0.91, "epilepsy": 0.83, "hemiplegia": 0.3903,
    })
    sequela_costs_2024: Mapping[str, float] = field(default_factory=lambda: {
        "cognitive_annual": 73.36, "hearing": 1020.0, "epilepsy": 329.06,
        "hemiplegia": 998.475, "implant": 24829.415, "special_education_annual": 2530.115,
    })
    vaccine_price_2024: Mapping[str, float] = field(default_factory=lambda: {"PCV13": 4.22, "PCV10-SII": 2.0})
    delivery_cost_per_dose: float = 2.1
    wastage_rate: float = 0.05
    aefi_per_100k_doses: float = 4.01
    aefi_cost_2024: float = 4800.355
    pv_productivity_loss_per_death: float = 42392.92
    life_expectancy_years: float = 74.0
    discount_rate: float = 0.03


@dataclass(frozen=True)
class Scenario:
    name: str
    vaccine: Optional[str]
    schedule_mode: str = "on_schedule"
    catchup: bool = False
    waning: bool = True
    herd_strength: float = 1.0
    herd_mode: str = "gradual"  # off, immediate_mature, or gradual


def default_scenarios() -> list[Scenario]:
    return [
        Scenario("no_vaccination", None, herd_strength=0.0, herd_mode="off"),
        Scenario("pcv13_on_schedule", "PCV13", "on_schedule", False, False, 1.0),
        Scenario("pcv13_on_schedule_herd_off", "PCV13", "on_schedule", False, False, 0.0, "off"),
        Scenario("pcv13_delayed", "PCV13", "delayed", False, True, 1.0),
        Scenario("pcv13_delayed_herd_off", "PCV13", "delayed", False, True, 0.0, "off"),
        Scenario("pcv13_missed_no_catchup", "PCV13", "missed", False, True, 1.0),
        Scenario("pcv13_delayed_catchup", "PCV13", "delayed", True, True, 1.0),
    ]


def age_weight(age: int, endpoint: str) -> float:
    if endpoint == "meningitis":
        if age < 12:
            return 0.65 / 12.0
        if age < 24:
            return 0.18 / 12.0
        if age < 36:
            return 0.07 / 12.0
        return 0.10 / 24.0
    if endpoint == "aom":
        if age < 12:
            return 0.1869 / 12.0
        if age < 24:
            return 0.1889 / 12.0
        return 0.6242 / 36.0
    if age < 12:
        return 0.50 / 12.0
    if age < 24:
        return 0.25 / 12.0
    if age < 36:
        return 0.15 / 12.0
    return 0.10 / 24.0


def carriage_age_weight(age: int, parameters: ModelParameters) -> float:
    if age <= 5:
        return parameters.carriage_age_weights["0_5"]
    if age <= 17:
        return parameters.carriage_age_weights["6_17"]
    if age <= 23:
        return parameters.carriage_age_weights["18_23"]
    return parameters.carriage_age_weights["24_59"]


def discount_factor(time_years: float, rate: float) -> float:
    return 1.0 if rate == 0 else 1.0 / ((1.0 + rate) ** max(time_years, 0.0))


def discounted_life_years_lost(age_years: float, life_expectancy: float, discount_rate: float) -> float:
    if age_years >= life_expectancy:
        return 0.0
    total = 0.0
    for year in range(int(np.floor(age_years)), int(np.ceil(life_expectancy))):
        left = max(age_years, float(year))
        right = min(life_expectancy, float(year + 1))
        if right > left:
            total += (right - left) * discount_factor((left + right) / 2.0, discount_rate)
    return total


def routine_transition_map(mode: str) -> dict[str, dict[int, float]]:
    if mode == "on_schedule":
        return {"dose1": {2: 1.0}, "dose2": {4: 1.0}, "booster": {12: 1.0}}
    if mode == "delayed":
        return {
            "dose1": {2: 0.65, 3: 0.57, 4: 0.67, 5: 0.60},
            "dose2": {4: 0.70, 5: 0.45, 6: 0.30, 7: 0.20, 8: 0.10},
            "booster": {12: 0.55, 13: 0.35, 14: 0.25, 15: 0.15, 16: 0.10},
        }
    if mode == "missed":
        return {"dose1": {2: 1.0}, "dose2": {}, "booster": {}}
    custom_schedules = {
        "schedule_3_5_12": (3, 5, 12),
        "schedule_2_5_12": (2, 5, 12),
        "schedule_3_6_12": (3, 6, 12),
        "schedule_2_4_15": (2, 4, 15),
        "schedule_3_5_15": (3, 5, 15),
    }
    if mode in custom_schedules:
        dose1, dose2, booster = custom_schedules[mode]
        return {"dose1": {dose1: 1.0}, "dose2": {dose2: 1.0}, "booster": {booster: 1.0}}
    raise ValueError(f"Unknown schedule mode: {mode}")


def dose_relative_ve(group: int, age: int, waning: bool, waning_severity: float = 1.0) -> float:
    def scaled_waning(peak_value: float, age_specific_value: float) -> float:
        """Scale the decline from the dose-specific peak effect.

        A severity of zero must retain the peak effect (rather than forcing
        every dose to 1.0); a severity of one recovers the base age profile.
        """
        return max(0.0, peak_value - waning_severity * (peak_value - age_specific_value))

    if group == G_DOSE1:
        return 0.45
    if group == G_DOSE2:
        if not waning:
            return 1.0
        return scaled_waning(1.0, 0.70 if age >= 36 else 0.85 if age >= 24 else 1.0)
    if group == G_BOOSTER:
        if not waning:
            return 1.10
        return scaled_waning(1.10, 0.70 if age >= 48 else 0.80 if age >= 36 else 0.95 if age >= 24 else 1.10)
    return 0.0


def endpoint_base_ve(vaccine: str, endpoint: str, parameters: ModelParameters) -> float:
    key = "ipd" if endpoint in {"meningitis", "sepsis"} else endpoint
    return float(parameters.direct_ve[vaccine][key])


def apply_transitions(
    alive: np.ndarray,
    scenario: Scenario,
    month: int,
) -> tuple[np.ndarray, dict[str, float]]:
    """Move surviving children through dose states and return administered doses."""
    updated = alive.copy()
    counts = {"dose1": 0.0, "dose2": 0.0, "booster": 0.0}
    if scenario.vaccine is None:
        return updated, counts
    transitions = routine_transition_map(scenario.schedule_mode)
    for label, source, target in [
        ("booster", G_DOSE2, G_BOOSTER),
        ("dose2", G_DOSE1, G_DOSE2),
        ("dose1", G_DOSE0, G_DOSE1),
    ]:
        for age, probability in transitions[label].items():
            moved = updated[source, age] * float(probability)
            updated[source, age] -= moved
            updated[target, age] += moved
            counts[label] += moved
    if scenario.catchup and month < 12:
        for ages, probability in [((7, 12), 0.12), ((12, 24), 0.08)]:
            for age in range(*ages):
                moved = updated[G_UNVACCINATED, age] * probability
                updated[G_UNVACCINATED, age] -= moved
                updated[G_DOSE1, age] += moved
                counts["dose1"] += moved
    return updated, counts


def carriage_index(alive: np.ndarray, parameters: ModelParameters) -> tuple[float, float]:
    weighted_total = 0.0
    weighted_reduction = 0.0
    for age in range(AGE_MONTHS):
        weight = carriage_age_weight(age, parameters)
        if weight == 0:
            continue
        weighted_total += weight * float(alive[:, age].sum())
        for group, name in enumerate(GROUPS):
            weighted_reduction += weight * float(alive[group, age]) * float(parameters.carriage_reduction.get(name, 0.0))
    x = weighted_reduction / weighted_total if weighted_total else 0.0
    reference = parameters.coverage * parameters.effective_dose_share_reference * parameters.carriage_reduction["dose2"]
    return x, min(1.0, max(0.0, x / reference)) if reference else 0.0


def mature_herd_multiplier(vaccine: str, endpoint: str, parameters: ModelParameters) -> float:
    """Return the endpoint-specific mature indirect-effect multiplier H_e*."""
    if parameters.mature_herd_multiplier_override is not None:
        override = parameters.mature_herd_multiplier_override.get(vaccine, {})
        if endpoint in override:
            return min(1.0, max(0.0, float(override[endpoint])))
    target = float(parameters.mature_population_irr[vaccine]["ipd" if endpoint in {"meningitis", "sepsis"} else endpoint])
    direct_reference = parameters.coverage * parameters.effective_dose_share_reference * endpoint_base_ve(vaccine, endpoint, parameters)
    if endpoint == "aom":
        numerator = target - (1.0 - parameters.spn_aom_proportion)
        denominator = parameters.spn_aom_proportion * max(1e-12, 1.0 - direct_reference)
        mature = numerator / denominator
    else:
        mature = target / max(1e-12, 1.0 - direct_reference)
    return min(1.0, max(0.0, mature))


def herd_multiplier(
    scenario: Scenario,
    endpoint: str,
    relative_carriage: float,
    parameters: ModelParameters,
) -> tuple[float, float]:
    """Return H_e(t) and its mature target H_e* for an endpoint.

    H_e(t)=1-s(1-H_e*)r(t), where s is the herd-strength switch and r(t)
    is the normalized carriage-reduction index. For all-cause AOM, only the
    pneumococcal-attributable component is affected by direct/indirect effects.
    """
    if scenario.vaccine is None or scenario.herd_strength <= 0 or scenario.herd_mode == "off":
        return 1.0, 1.0
    vaccine = scenario.vaccine
    mature = mature_herd_multiplier(vaccine, endpoint, parameters)
    if scenario.herd_mode == "immediate_mature":
        build_index = 1.0
    elif scenario.herd_mode == "gradual":
        build_index = relative_carriage
    else:
        raise ValueError(f"Unknown herd mode: {scenario.herd_mode}")
    strength = min(1.0, max(0.0, scenario.herd_strength))
    h = 1.0 - strength * build_index * (1.0 - mature)
    return min(1.0, max(mature, h)), mature


def initialize_population(parameters: ModelParameters) -> np.ndarray:
    """Create a no-program age distribution through a documented 60-month warm-up."""
    alive = np.zeros((len(GROUPS), AGE_MONTHS), dtype=float)
    monthly_births = parameters.annual_births / 12.0
    bg = 1.0 - (1.0 - parameters.under5_mortality_per_1000 / 1000.0) ** (1.0 / 60.0)
    for _ in range(parameters.warmup_months):
        alive[G_UNVACCINATED, 0] += monthly_births * (1.0 - parameters.neonatal_mortality_per_1000 / 1000.0)
        alive *= 1.0 - bg
        alive[:, 1:] = alive[:, :-1]
        alive[:, 0] = 0.0
    return alive


def move_ages(alive: np.ndarray) -> np.ndarray:
    shifted = np.zeros_like(alive)
    shifted[:, 1:] = alive[:, :-1]
    return shifted


def _flow_key(endpoint: str) -> str:
    return {"meningitis": "meningitis", "pneumonia": "pneumonia", "sepsis": "sepsis", "aom": "aom"}[endpoint]


def simulate_calendar_pde(scenario: Scenario, parameters: ModelParameters = ModelParameters()) -> dict[str, pd.DataFrame]:
    """Simulate a calendar under-five population and return auditable source data."""
    alive = initialize_population(parameters)
    sequelae: MutableMapping[str, float] = {key: 0.0 for key in ("cognitive", "hearing", "epilepsy", "hemiplegia", "implant", "special_education")}
    monthly_births = parameters.annual_births / 12.0
    bg = 1.0 - (1.0 - parameters.under5_mortality_per_1000 / 1000.0) ** (1.0 / 60.0)
    month_rows: list[dict[str, float | str]] = []
    age_rows: list[dict[str, float | int | str]] = []
    dose_rows: list[dict[str, float | int | str]] = []
    for month in range(parameters.horizon_months):
        population_before_births = float(alive.sum())
        newborn_survivors = monthly_births * (1.0 - parameters.neonatal_mortality_per_1000 / 1000.0)
        if scenario.vaccine is None:
            alive[G_UNVACCINATED, 0] += newborn_survivors
        else:
            alive[G_DOSE0, 0] += newborn_survivors * parameters.coverage
            alive[G_UNVACCINATED, 0] += newborn_survivors * (1.0 - parameters.coverage)
        alive, doses = apply_transitions(alive, scenario, month)
        population_before_mortality = float(alive.sum())
        x, relative_x = carriage_index(alive, parameters)
        herds = {endpoint: herd_multiplier(scenario, endpoint, relative_x, parameters) for endpoint in ENDPOINTS}
        flows = {endpoint: 0.0 for endpoint in ENDPOINTS}
        death_flows = {endpoint: 0.0 for endpoint in ("meningitis", "pneumonia", "sepsis")}
        age_cases = {endpoint: np.zeros(AGE_MONTHS, dtype=float) for endpoint in ENDPOINTS}
        age_deaths = np.zeros(AGE_MONTHS, dtype=float)
        substeps = max(1, int(parameters.reaction_substeps))
        substep_bg = 1.0 - (1.0 - bg) ** (1.0 / substeps)
        background_deaths = 0.0
        for _ in range(substeps):
            for age in range(AGE_MONTHS):
                population = alive[:, age].copy()
                if population.sum() <= 0:
                    continue
                substep_cases = {}
                for endpoint in ENDPOINTS:
                    baseline = parameters.annual_incidence_per_100k[endpoint] * 5.0 * age_weight(age, endpoint) / 100_000.0 / substeps
                    h, _ = herds[endpoint]
                    group_risk = np.zeros(len(GROUPS), dtype=float)
                    for group in range(len(GROUPS)):
                        if scenario.vaccine is None:
                            ve = 0.0
                        else:
                            ve = min(parameters.max_direct_ve, endpoint_base_ve(scenario.vaccine, endpoint, parameters) * dose_relative_ve(group, age, scenario.waning, parameters.waning_severity))
                        if endpoint == "aom":
                            group_risk[group] = baseline * ((1.0 - parameters.spn_aom_proportion) + parameters.spn_aom_proportion * h * (1.0 - ve))
                        else:
                            group_risk[group] = baseline * h * (1.0 - ve)
                    events = float((population * np.clip(group_risk, 0.0, 1.0)).sum())
                    flows[endpoint] += events
                    age_cases[endpoint][age] += events
                    substep_cases[endpoint] = events
                severe_pneumonia = substep_cases["pneumonia"] * parameters.severe_pneumonia_proportion
                endpoint_deaths = {
                    "meningitis": substep_cases["meningitis"] * parameters.cfr["meningitis"],
                    "pneumonia": severe_pneumonia * parameters.cfr["pneumonia"],
                    "sepsis": substep_cases["sepsis"] * parameters.cfr["sepsis"],
                }
                total_endpoint_deaths = sum(endpoint_deaths.values())
                death_flows = {key: death_flows[key] + value for key, value in endpoint_deaths.items()}
                age_deaths[age] += total_endpoint_deaths
                background_deaths += population.sum() * substep_bg
                mortality_fraction = min(1.0, substep_bg + total_endpoint_deaths / population.sum())
                alive[:, age] = population * (1.0 - mortality_fraction)

        # Current sequelae stock approximation: a transparent extension of the
        # closed-cohort source accounting, with monthly all-cause survival.
        for key in sequelae:
            sequelae[key] *= 1.0 - bg
        meningitis_survivors = flows["meningitis"] * (1.0 - parameters.cfr["meningitis"])
        sequelae["cognitive"] += meningitis_survivors * parameters.sequela_probabilities["cognitive"]
        sequelae["hearing"] += meningitis_survivors * parameters.sequela_probabilities["hearing"]
        sequelae["epilepsy"] += meningitis_survivors * parameters.sequela_probabilities["epilepsy"]
        sequelae["hemiplegia"] += meningitis_survivors * parameters.sequela_probabilities["hemiplegia"]
        sequelae["implant"] += meningitis_survivors * parameters.sequela_probabilities["hearing"] * parameters.sequela_probabilities["implant"]
        sequelae["special_education"] += meningitis_survivors * parameters.sequela_probabilities["hearing"] * parameters.sequela_probabilities["special_education"]

        time_years = month / 12.0
        disc = discount_factor(time_years, parameters.discount_rate)
        pn_inpatient = flows["pneumonia"] * parameters.severe_pneumonia_proportion
        pn_outpatient = flows["pneumonia"] * (1.0 - parameters.severe_pneumonia_proportion)
        acute_cost = (
            pn_inpatient * parameters.disease_costs_2024["pneumonia_inpatient"]
            + pn_outpatient * parameters.disease_costs_2024["pneumonia_outpatient"]
            + flows["meningitis"] * parameters.disease_costs_2024["meningitis"]
            + flows["sepsis"] * parameters.disease_costs_2024["sepsis"]
            + flows["aom"] * parameters.disease_costs_2024["aom"]
        )
        sequela_cost = (
            sequelae["cognitive"] * parameters.sequela_costs_2024["cognitive_annual"] / 12.0
            + sequelae["special_education"] * parameters.sequela_costs_2024["special_education_annual"] / 12.0
            + meningitis_survivors * (
                parameters.sequela_probabilities["hearing"] * parameters.sequela_costs_2024["hearing"]
                + parameters.sequela_probabilities["epilepsy"] * parameters.sequela_costs_2024["epilepsy"]
                + parameters.sequela_probabilities["hemiplegia"] * parameters.sequela_costs_2024["hemiplegia"]
                + parameters.sequela_probabilities["hearing"] * parameters.sequela_probabilities["implant"] * parameters.sequela_costs_2024["implant"]
            )
        )
        administered_doses = sum(doses.values())
        vaccine_cost = 0.0
        if scenario.vaccine is not None:
            procured_doses = administered_doses * (1.0 + parameters.wastage_rate)
            vaccine_cost = procured_doses * (parameters.vaccine_price_2024[scenario.vaccine] + parameters.delivery_cost_per_dose)
            vaccine_cost += procured_doses * parameters.aefi_per_100k_doses / 100_000.0 * parameters.aefi_cost_2024
        death_cost = sum(death_flows.values()) * parameters.pv_productivity_loss_per_death
        acute_qaly = (
            pn_inpatient * parameters.acute_los_days["pneumonia_inpatient"] / 365.0 * (1.0 - parameters.acute_utilities["pneumonia_inpatient"])
            + pn_outpatient * parameters.acute_los_days["pneumonia_outpatient"] / 365.0 * (1.0 - parameters.acute_utilities["pneumonia_outpatient"])
            + flows["meningitis"] * parameters.acute_los_days["meningitis"] / 365.0 * (1.0 - parameters.acute_utilities["meningitis"])
            + flows["sepsis"] * parameters.acute_los_days["sepsis"] / 365.0 * (1.0 - parameters.acute_utilities["sepsis"])
            + flows["aom"] * parameters.acute_los_days["aom"] / 365.0 * (1.0 - parameters.acute_utilities["aom"])
        )
        sequela_qaly = sum(sequelae[key] / 12.0 * (1.0 - parameters.sequela_utilities[key]) for key in parameters.sequela_utilities)
        death_qaly = sum(age_deaths[age] * discounted_life_years_lost(age / 12.0, parameters.life_expectancy_years, parameters.discount_rate) for age in range(AGE_MONTHS))
        total_alive = float(alive.sum())
        disease_deaths = sum(death_flows.values())
        expected_alive = population_before_mortality - background_deaths - disease_deaths
        mass_residual = total_alive - expected_alive
        month_rows.append({
            "scenario": scenario.name,
            "month": month + 1,
            "year": month // 12 + 1,
            "calendar_time_years": time_years,
            "A_integral": total_alive,
            "x_carriage_reduction_index": x,
            "x_relative_to_mature": relative_x,
            "H_ipd": herds["meningitis"][0], "H_pneumonia": herds["pneumonia"][0], "H_aom": herds["aom"][0],
            "H_ipd_mature": herds["meningitis"][1], "H_pneumonia_mature": herds["pneumonia"][1], "H_aom_mature": herds["aom"][1],
            "cases_meningitis": flows["meningitis"], "cases_pneumonia": flows["pneumonia"], "cases_sepsis": flows["sepsis"], "cases_aom_allcause": flows["aom"],
            "cases_total": sum(flows.values()), "deaths_meningitis": death_flows["meningitis"], "deaths_pneumonia": death_flows["pneumonia"], "deaths_sepsis": death_flows["sepsis"], "deaths_total": sum(death_flows.values()),
            "dose1_administered": doses["dose1"], "dose2_administered": doses["dose2"], "booster_administered": doses["booster"], "doses_administered": administered_doses,
            "pv_acute_cost": acute_cost * disc, "pv_sequela_cost": sequela_cost * disc, "pv_vaccine_cost": vaccine_cost * disc, "pv_productivity_death_cost": death_cost * disc,
            "pv_total_cost": (acute_cost + sequela_cost + vaccine_cost + death_cost) * disc,
            "pv_acute_qaly_loss": acute_qaly * disc, "pv_sequela_qaly_loss": sequela_qaly * disc, "pv_death_qaly_loss": death_qaly, "pv_total_qaly_loss": acute_qaly * disc + sequela_qaly * disc + death_qaly,
            "population_before_births": population_before_births,
            "newborn_survivors": newborn_survivors,
            "background_deaths": background_deaths,
            "mass_balance_residual": mass_residual,
        })
        for age in range(AGE_MONTHS):
            for group, group_name in enumerate(GROUPS):
                age_rows.append({"scenario": scenario.name, "month": month + 1, "age_month": age, "group": group_name, "A": float(alive[group, age])})
        dose_rows.append({"scenario": scenario.name, "month": month + 1, **doses})
        alive = move_ages(alive)
    monthly = pd.DataFrame(month_rows)
    age_state = pd.DataFrame(age_rows)
    doses = pd.DataFrame(dose_rows)
    return {"monthly": monthly, "age_state": age_state, "dose_counts": doses}


def summarize_scenario(result: Mapping[str, pd.DataFrame], parameters: ModelParameters) -> pd.DataFrame:
    monthly = result["monthly"]
    totals = monthly[[column for column in monthly.columns if column.startswith(("cases_", "deaths_", "pv_", "doses_"))]].sum(numeric_only=True).to_dict()
    births = parameters.annual_births * parameters.horizon_months / 12.0
    totals.update({"scenario": str(monthly["scenario"].iloc[0]), "normalization_births": births})
    return pd.DataFrame([totals])


def factorial_scenarios(vaccine: str = "PCV13") -> Iterable[Scenario]:
    """Full 2^3 design: delayed dose states, waning, and herd build-up."""
    for delayed, waning, herd in product((False, True), repeat=3):
        label = f"{vaccine.lower()}_delay{int(delayed)}_waning{int(waning)}_gradualherd{int(herd)}"
        herd_mode = "gradual" if herd else "immediate_mature"
        yield Scenario(label, vaccine, "delayed" if delayed else "on_schedule", False, waning, 1.0, herd_mode)
