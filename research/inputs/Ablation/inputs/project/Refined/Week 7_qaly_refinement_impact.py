"""
Economic comparison across PCV model refinements.

The script builds a single output bundle for:

1. the original QALY_new cohort outputs for deterministic ODE and stochastic
   Markov/microsimulation analyses;
2. a calendar under-five PDE model with on-schedule primary-series protection
   and IRR-decomposed herd protection;
3. a calendar under-five PDE model with delayed dose states, waning, seasonal
   births, and carriage-informed herd build-up.

All public-facing model labels describe scientific/model features rather than
development chronology. Historical source file names are used only to locate
the existing project outputs.
"""

from __future__ import annotations

import copy
import importlib.util
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REFINED_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = REFINED_DIR.parent
RESULTS_CSV_DIR = PROJECT_ROOT / "Results" / "Results_csv"
OUTPUT_ROOT = REFINED_DIR / "Outputs"
WEEK7_OUTPUT_DIR = OUTPUT_ROOT / "Week 7_qaly_refinement"
CSV_DIR = WEEK7_OUTPUT_DIR / "csv"
FIGURE_DIR = WEEK7_OUTPUT_DIR / "figures"
SUMMARY_MD = REFINED_DIR / "Week 7_qaly_refinement_methods_and_results.md"


def exactly_one(items: Iterable, label: str):
    matches = list(items)
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one {label}, found {len(matches)}")
    return matches[0]


def newest_path(paths: Iterable[Path], label: str) -> Path:
    matches = list(paths)
    if not matches:
        raise FileNotFoundError(f"Could not find {label}")
    return max(matches, key=lambda path: path.stat().st_mtime)


def calendar_yearly_file(source_dir: Path, scenario: str) -> Path:
    patterns = [
        f"calendar_pde_refined_{scenario}.csv",
        f"calendar_pde_program_realism_{scenario}.csv",
        f"*_pde_calendar_refined_{scenario}.csv",
    ]
    candidates: List[Path] = []
    for pattern in patterns:
        candidates.extend(sorted(source_dir.glob(pattern)))
    return exactly_one(
        candidates,
        f"calendar yearly CSV for {scenario} in {source_dir}",
    )


HISTORICAL_OUTPUT_DIR = OUTPUT_ROOT / "Week 6_calendar_pde_program_realism"
ON_SCHEDULE_SOURCE_DIR = HISTORICAL_OUTPUT_DIR / "on_schedule_reference_csv"
DELAYED_DOSE_SOURCE_DIR = HISTORICAL_OUTPUT_DIR / "csv"

MODEL_ORIGINAL_ODE = "Original cohort ODE + QALY_new"
MODEL_ORIGINAL_MARKOV = "Original cohort Markov + QALY_new"
MODEL_CALENDAR_ON_SCHEDULE = "Calendar PDE: on-schedule primary-series protection + IRR-decomposed herd"
MODEL_CALENDAR_DELAYED = "Calendar PDE: delayed dose states + waning + carriage-informed herd"

MODEL_FEATURES = {
    MODEL_ORIGINAL_ODE: {
        "model_family": "Closed birth-cohort",
        "population_estimand": "One Vietnamese birth cohort followed from birth to 59 months",
        "time_horizon_years": 5,
        "births_denominator": 1_384_610.0,
        "feature_summary": "Deterministic ODE with original total-IRR vaccine effect and QALY_new economics.",
    },
    MODEL_ORIGINAL_MARKOV: {
        "model_family": "Closed birth-cohort",
        "population_estimand": "One Vietnamese birth cohort followed from birth to 59 months",
        "time_horizon_years": 5,
        "births_denominator": 1_384_610.0,
        "feature_summary": "Stochastic Markov/microsimulation counterpart with original QALY_new economics.",
    },
    MODEL_CALENDAR_ON_SCHEDULE: {
        "model_family": "Calendar under-five PDE",
        "population_estimand": "Full 0-59 month calendar population with continuing births",
        "time_horizon_years": 10,
        "births_denominator": 13_846_100.0,
        "feature_summary": "Baseline-only warm-up, continuing births, endpoint-specific VE, protection beginning after primary series, and herd multiplier inferred from annual IRR.",
    },
    MODEL_CALENDAR_DELAYED: {
        "model_family": "Calendar under-five PDE",
        "population_estimand": "Full 0-59 month calendar population with continuing births",
        "time_horizon_years": 10,
        "births_denominator": 13_846_100.0,
        "feature_summary": "Adds mild birth seasonality, explicit dose states, delayed routine vaccination, dose-specific VE, waning, and carriage-informed herd build-up.",
    },
}

VACCINE_ORDER = ["PCV10-GSK", "PCV13", "PCV10-SII"]

CALENDAR_SCENARIOS = [
    ("baseline", "Baseline", None),
    ("pcv10_gsk_rollout_to_stable", "PCV10-GSK", "pcv10"),
    ("pcv13_rollout_to_stable", "PCV13", "pcv13"),
    ("pcv10sii_rollout_to_stable", "PCV10-SII", "pcv10sii"),
]

COHORT_SCENARIO_LABELS = {
    "combined_baseline": ("Baseline", "No vaccination"),
    "combined_pcv10_dynamic": ("PCV10-GSK", "Dynamic rollout IRR"),
    "combined_pcv13_dynamic": ("PCV13", "Dynamic rollout IRR"),
    "combined_pcv10_year5_immediate": ("PCV10-GSK", "Immediate mature IRR"),
    "combined_pcv13_year5_immediate": ("PCV13", "Immediate mature IRR"),
    "combined_pcv10sii_year5_immediate": ("PCV10-SII", "Immediate mature IRR"),
}

CALENDAR_YEARLY_FILES = {
    MODEL_CALENDAR_ON_SCHEDULE: {
        "baseline": calendar_yearly_file(ON_SCHEDULE_SOURCE_DIR, "baseline"),
        "pcv10_gsk_rollout_to_stable": calendar_yearly_file(ON_SCHEDULE_SOURCE_DIR, "pcv10_gsk_rollout_to_stable"),
        "pcv13_rollout_to_stable": calendar_yearly_file(ON_SCHEDULE_SOURCE_DIR, "pcv13_rollout_to_stable"),
        "pcv10sii_rollout_to_stable": calendar_yearly_file(ON_SCHEDULE_SOURCE_DIR, "pcv10sii_rollout_to_stable"),
    },
    MODEL_CALENDAR_DELAYED: {
        "baseline": calendar_yearly_file(DELAYED_DOSE_SOURCE_DIR, "baseline"),
        "pcv10_gsk_rollout_to_stable": calendar_yearly_file(DELAYED_DOSE_SOURCE_DIR, "pcv10_gsk_rollout_to_stable"),
        "pcv13_rollout_to_stable": calendar_yearly_file(DELAYED_DOSE_SOURCE_DIR, "pcv13_rollout_to_stable"),
        "pcv10sii_rollout_to_stable": calendar_yearly_file(DELAYED_DOSE_SOURCE_DIR, "pcv10sii_rollout_to_stable"),
    },
}

EXPECTED_EPI_COLUMNS = [
    "Year",
    "New_Meningitis",
    "New_Pneumo_Sev",
    "New_Pneumo_NonSev",
    "New_Sepsis",
    "New_AOM",
    "New_Deaths_Pneumo",
    "Stock_Seq_Cognitive",
    "Stock_Seq_Hearing",
    "Stock_Seq_Epilepsy",
    "Stock_Seq_Hemiplegia",
    "Stock_Seq_Hearing_Implant",
    "Stock_Seq_Hearing_SpEd",
]

CALENDAR_CONFIG = {
    "BASE_YEAR": 2023,
    "HORIZON_YEARS": 10,
    "AVG_U5_DEATH_AGE_YEARS": 2.5,
    "AVG_VACCINE_EVENT_AGE_YEARS": 0.35,
    "PRICE_YEAR_ALIGNMENT_RATE": 0.03,
    "DISCOUNT_RATE": 0.03,
    "LOS_DAYS": {
        "inpatient_pneumonia": 8.74,
        "outpatient_pneumonia": 7.0,
        "meningitis": 18.41,
        "sepsis": 12.44,
        "aom": 7.0,
    },
}


def load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


qaly_source = load_module(PROJECT_ROOT / "QALY_new.py", "refinement_qaly_source")


def load_calendar_source_modules() -> Tuple[object, object]:
    source_paths = [
        REFINED_DIR / "Week 5_calendar_pde_refined.py",
        REFINED_DIR / "Week 6_calendar_pde_program_realism.py",
    ]
    modules = [load_module(path, f"calendar_feature_source_{idx}") for idx, path in enumerate(source_paths)]
    on_schedule = exactly_one(
        (
            module
            for module in modules
            if not hasattr(module, "GROUP_NONPROGRAM")
        ),
        "on-schedule calendar source module",
    )
    delayed = exactly_one(
        (
            module
            for module in modules
            if hasattr(module, "GROUP_NONPROGRAM") and hasattr(module, "ASSUMED_SEASONAL_BIRTH_SHARES")
        ),
        "delayed-dose calendar source module",
    )
    return on_schedule, delayed


def source_function(module, prefix: str, suffix: str):
    return exactly_one(
        (getattr(module, name) for name in dir(module) if name.startswith(prefix) and name.endswith(suffix)),
        f"source function matching {prefix}...{suffix}",
    )


on_schedule_source, delayed_dose_source = load_calendar_source_modules()
build_delayed_initial_population = source_function(delayed_dose_source, "build_", "_initial_population")
add_delayed_record_births = source_function(delayed_dose_source, "add_record_births_", "")
solve_delayed_month = source_function(delayed_dose_source, "solve_one_month_", "")
age_shift_delayed = source_function(delayed_dose_source, "age_shift_", "")


def ensure_output_dirs() -> None:
    CSV_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def require_inputs() -> None:
    missing = []
    for path in [HISTORICAL_OUTPUT_DIR, ON_SCHEDULE_SOURCE_DIR, DELAYED_DOSE_SOURCE_DIR]:
        if not path.exists():
            missing.append(path)
    for model_files in CALENDAR_YEARLY_FILES.values():
        for path in model_files.values():
            if not path.exists():
                missing.append(path)
    for path in [RESULTS_CSV_DIR / "results_summary_ODE.csv", RESULTS_CSV_DIR / "results_summary_Markov.csv"]:
        if not path.exists():
            missing.append(path)
    if missing:
        missing_text = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing required input files:\n{missing_text}")


def positive_range_midpoint(value: str) -> float:
    text = str(value).replace(",", "")
    nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", text)]
    if not nums:
        raise ValueError(f"Cannot parse range: {value}")
    if len(nums) == 1:
        return nums[0]
    return (nums[0] + nums[1]) / 2.0


def align_money(value: float, source_year: Optional[int]) -> float:
    if source_year is None:
        return float(value)
    delta = int(source_year) - int(CALENDAR_CONFIG["BASE_YEAR"])
    if delta <= 0:
        return float(value)
    return float(value) / ((1.0 + float(CALENDAR_CONFIG["PRICE_YEAR_ALIGNMENT_RATE"])) ** delta)


def build_calendar_params() -> Dict:
    cfg = copy.deepcopy(qaly_source.CONFIG)
    cfg["BASE_COHORT_YEAR"] = CALENDAR_CONFIG["BASE_YEAR"]
    cfg["PRICE_YEAR_ALIGNMENT_RATE"] = CALENDAR_CONFIG["PRICE_YEAR_ALIGNMENT_RATE"]
    cfg["DISCOUNT_RATE"] = CALENDAR_CONFIG["DISCOUNT_RATE"]
    cfg["LOS_DAYS"] = CALENDAR_CONFIG["LOS_DAYS"]
    cfg["HORIZON_YEARS"] = CALENDAR_CONFIG["HORIZON_YEARS"]
    params = qaly_source.build_parameters_from_hardcoded_table(cfg)

    raw = qaly_source.RAW_PARAM_TABLE
    medical = raw["cost_medical_2024"]
    params["costs_sequelae"]["hearing_case"] = align_money(positive_range_midpoint(medical["hearing_case_range"]), 2024)
    params["costs_sequelae"]["hemiplegia_case"] = align_money(positive_range_midpoint(medical["hemiplegia_case_range"]), 2024)
    params["costs_sequelae"]["cochlear_implant_case"] = align_money(positive_range_midpoint(medical["cochlear_implant_case_range"]), 2024)
    params["costs_sequelae"]["special_education_annual"] = align_money(
        positive_range_midpoint(raw["special_education_annual_range_2025"]),
        2025,
    )
    params["ae_cost_per_case"] = align_money(positive_range_midpoint(raw["severe_ae_cost_range_unknown_year"]), None)
    return params


def discount_factor(time_years: float, rate: float) -> float:
    return qaly_source.discount_factor(float(time_years), float(rate))


def future_qaly_from_death_event(age_at_death: float, life_expectancy: float, rate: float) -> float:
    if age_at_death >= life_expectancy:
        return 0.0
    remaining = life_expectancy - age_at_death
    total = 0.0
    full_years = int(math.floor(remaining))
    for year in range(full_years):
        total += discount_factor(year + 0.5, rate)
    tail = remaining - full_years
    if tail > 1e-12:
        total += tail * discount_factor(full_years + tail / 2.0, rate)
    return float(total)


def productivity_loss_from_death_event(age_at_death: float, params: Dict) -> float:
    rate = float(params["discount_rate"])
    birth_anchored = qaly_source.pv_prod_loss_from_death_age(age_at_death, params["prod_discounted_age_contrib"])
    age_df = max(discount_factor(age_at_death, rate), 1e-12)
    return float(birth_anchored / age_df)


def load_calendar_yearly_epi(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = [col for col in EXPECTED_EPI_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing columns: {missing}")
    out = df[EXPECTED_EPI_COLUMNS].copy()
    for col in EXPECTED_EPI_COLUMNS:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    required = ["New_Meningitis", "New_Pneumo_Sev", "New_Pneumo_NonSev", "New_Sepsis", "New_AOM", "New_Deaths_Pneumo"]
    if out[required].isna().any().any():
        raise ValueError(f"{path} has nonnumeric required disease-flow values")
    return out


def vaccine_unit_cost_components(params: Dict, vaccine_key: str) -> Tuple[float, float, float]:
    price = float(params["vaccine_prices"][vaccine_key])
    delivery = float(params["service_delivery_per_dose"])
    wastage = float(params["wastage_rate"])
    ae_rate = float(params["ae_incidence_per_100k_doses"]) / 100000.0
    ae_cost = float(params["ae_cost_per_case"])
    return price * (1.0 + wastage), delivery, ae_rate * ae_cost


def cost_rows_from_dose_counts(
    model_label: str,
    scenario: str,
    vaccine_label: str,
    vaccine_key: Optional[str],
    dose_counts_by_year: Dict[int, Dict[str, float]],
    params: Dict,
) -> List[Dict]:
    rows = []
    rate = float(params["discount_rate"])
    for year in range(1, CALENDAR_CONFIG["HORIZON_YEARS"] + 1):
        counts = dose_counts_by_year.get(year, {})
        dose1 = float(counts.get("dose1", 0.0))
        dose2 = float(counts.get("dose2", 0.0))
        booster = float(counts.get("booster", 0.0))
        total_administered = dose1 + dose2 + booster
        if vaccine_key is None or total_administered <= 0.0:
            raw_procurement = raw_delivery = raw_ae = raw_total = 0.0
        else:
            procurement_unit, delivery_unit, ae_unit = vaccine_unit_cost_components(params, vaccine_key)
            raw_procurement = total_administered * procurement_unit
            raw_delivery = total_administered * delivery_unit
            raw_ae = total_administered * ae_unit
            raw_total = raw_procurement + raw_delivery + raw_ae
        event_time = (year - 1) + float(CALENDAR_CONFIG["AVG_VACCINE_EVENT_AGE_YEARS"])
        rows.append({
            "model_label": model_label,
            "scenario": scenario,
            "vaccine_label": vaccine_label,
            "year": year,
            "vaccine_key": vaccine_key,
            "dose1_administered": dose1,
            "dose2_administered": dose2,
            "booster_administered": booster,
            "total_administered_doses": total_administered,
            "raw_vaccine_procurement_cost": raw_procurement,
            "raw_delivery_cost": raw_delivery,
            "raw_ae_cost": raw_ae,
            "raw_vaccine_program_cost": raw_total,
            "pv_vaccine_program_cost": raw_total * discount_factor(event_time, rate),
        })
    return rows


def on_schedule_dose_counts(vaccine_key: Optional[str]) -> Dict[int, Dict[str, float]]:
    if vaccine_key is None:
        return {}
    p = copy.deepcopy(on_schedule_source.params_vietnam_micro)
    shares = on_schedule_source.base.normalize_birth_shares(p.get("birth_month_shares", [1.0 / 12.0] * 12))
    annual_births = float(p["birth_cohort"])
    coverage = float(p.get("vaccine_coverage", 0.0))
    live_multiplier = float(on_schedule_source.base.neonatal_survival_multiplier(p))
    covered_live_births = []
    for month in range(CALENDAR_CONFIG["HORIZON_YEARS"] * 12):
        gross = annual_births * float(shares[month % 12])
        covered_live_births.append(gross * live_multiplier * coverage)

    dose_counts: Dict[int, Dict[str, float]] = {}
    for month in range(CALENDAR_CONFIG["HORIZON_YEARS"] * 12):
        year = month // 12 + 1
        dose_counts.setdefault(year, {"dose1": 0.0, "dose2": 0.0, "booster": 0.0})
        for label, age_month in [("dose1", 2), ("dose2", 4), ("booster", 12)]:
            source_month = month - age_month
            if source_month >= 0:
                dose_counts[year][label] += covered_live_births[source_month]
    return dose_counts


def counted_transition(A: np.ndarray, source_group: int, target_group: int, probs: Dict[int, float]) -> float:
    total = 0.0
    for age, prob in delayed_dose_source._as_int_float_dict(probs).items():
        if age < 0 or age >= A.shape[1]:
            continue
        move_prob = max(0.0, min(1.0, float(prob)))
        moved = float(A[source_group, age]) * move_prob
        if moved <= 0.0:
            continue
        A[source_group, age] -= moved
        A[target_group, age] += moved
        total += moved
    return total


def apply_delayed_vaccination_events_counted(
    A: np.ndarray,
    p: Dict,
    calendar_month: int,
    intervention: bool,
    vaccine_type: Optional[str],
) -> Dict[str, float]:
    counts = {"dose1": 0.0, "dose2": 0.0, "booster": 0.0}
    if not intervention or vaccine_type is None:
        return counts

    if bool(p.get("catchup_enabled", False)) and int(calendar_month) < int(p.get("catchup_campaign_months", 12)):
        p_7_11 = max(0.0, min(1.0, float(p.get("catchup_7_11_to_dose1_prob", 0.0))))
        p_12_23 = max(0.0, min(1.0, float(p.get("catchup_12_23_to_dose1_prob", 0.0))))
        for age in range(7, min(12, A.shape[1])):
            moved = float(A[delayed_dose_source.GROUP_NONPROGRAM, age]) * p_7_11
            A[delayed_dose_source.GROUP_NONPROGRAM, age] -= moved
            A[delayed_dose_source.GROUP_PROGRAM_D1, age] += moved
            counts["dose1"] += moved
        for age in range(12, min(24, A.shape[1])):
            moved = float(A[delayed_dose_source.GROUP_NONPROGRAM, age]) * p_12_23
            A[delayed_dose_source.GROUP_NONPROGRAM, age] -= moved
            A[delayed_dose_source.GROUP_PROGRAM_D1, age] += moved
            counts["dose1"] += moved

    counts["booster"] += counted_transition(
        A,
        delayed_dose_source.GROUP_PROGRAM_D2,
        delayed_dose_source.GROUP_PROGRAM_D3,
        p.get("routine_booster_transition_probs", {}),
    )
    counts["dose2"] += counted_transition(
        A,
        delayed_dose_source.GROUP_PROGRAM_D1,
        delayed_dose_source.GROUP_PROGRAM_D2,
        p.get("routine_dose2_transition_probs", {}),
    )
    counts["dose1"] += counted_transition(
        A,
        delayed_dose_source.GROUP_PROGRAM_D0,
        delayed_dose_source.GROUP_PROGRAM_D1,
        p.get("routine_dose1_transition_probs", {}),
    )
    return counts


def delayed_dose_counts(
    scenario_name: str,
    intervention: bool,
    vaccine_type: Optional[str],
) -> Dict[int, Dict[str, float]]:
    if not intervention or vaccine_type is None:
        return {}
    p = copy.deepcopy(delayed_dose_source.params_vietnam_micro)
    p["record_months"] = CALENDAR_CONFIG["HORIZON_YEARS"] * 12
    shares = delayed_dose_source.base.normalize_birth_shares(
        p.get("birth_month_shares", delayed_dose_source.ASSUMED_SEASONAL_BIRTH_SHARES)
    )
    annual_births = float(p["birth_cohort"])
    A, S, _ = build_delayed_initial_population(p, shares, scenario_name, retain_preexisting_sequelae=False)
    dose_counts: Dict[int, Dict[str, float]] = {}

    for month in range(CALENDAR_CONFIG["HORIZON_YEARS"] * 12):
        year = month // 12 + 1
        dose_counts.setdefault(year, {"dose1": 0.0, "dose2": 0.0, "booster": 0.0})
        gross_births = annual_births * float(shares[month % 12])
        add_delayed_record_births(A, gross_births, p, intervention, vaccine_type)
        month_counts = apply_delayed_vaccination_events_counted(A, p, month, intervention, vaccine_type)
        for key, value in month_counts.items():
            dose_counts[year][key] += float(value)
        A_after, S, _, _ = solve_delayed_month(A, S, month, p, intervention=intervention, vaccine_type=vaccine_type)
        A, _ = age_shift_delayed(A_after)
    return dose_counts


def build_calendar_vaccine_cost_table(params: Dict) -> pd.DataFrame:
    rows: List[Dict] = []
    for scenario, label, vaccine_key in CALENDAR_SCENARIOS:
        rows.extend(cost_rows_from_dose_counts(
            MODEL_CALENDAR_ON_SCHEDULE,
            scenario,
            label,
            vaccine_key,
            on_schedule_dose_counts(vaccine_key),
            params,
        ))

    delayed_specs = {
        name: (intervention, vaccine_type)
        for name, _, intervention, vaccine_type, _ in delayed_dose_source._scenario_specs(CALENDAR_CONFIG["HORIZON_YEARS"] * 12)
    }
    for scenario, label, vaccine_key in CALENDAR_SCENARIOS:
        intervention, vaccine_type = delayed_specs[scenario]
        rows.extend(cost_rows_from_dose_counts(
            MODEL_CALENDAR_DELAYED,
            scenario,
            label,
            vaccine_key,
            delayed_dose_counts(scenario, intervention, vaccine_type),
            params,
        ))
    return pd.DataFrame(rows)


def positive_increment(current: float, previous: float) -> float:
    return max(float(current) - float(previous), 0.0)


def compute_calendar_economics(
    model_label: str,
    scenario: str,
    vaccine_label: str,
    vaccine_key: Optional[str],
    epi: pd.DataFrame,
    vaccine_costs: pd.DataFrame,
    params: Dict,
) -> pd.DataFrame:
    unit_cost = qaly_source.unit_cost_per_episode(params)
    util = params["utilities"]
    los = CALENDAR_CONFIG["LOS_DAYS"]
    rate = float(params["discount_rate"])
    avg_death_age = float(CALENDAR_CONFIG["AVG_U5_DEATH_AGE_YEARS"])
    life_exp = float(params["life_expectancy"])
    qaly_loss_per_death_event = future_qaly_from_death_event(avg_death_age, life_exp, rate)
    prod_loss_per_death_event = productivity_loss_from_death_event(avg_death_age, params)

    vcost = vaccine_costs[
        (vaccine_costs["model_label"] == model_label)
        & (vaccine_costs["scenario"] == scenario)
    ].set_index("year")

    prev_stocks = {
        "Stock_Seq_Hearing": 0.0,
        "Stock_Seq_Epilepsy": 0.0,
        "Stock_Seq_Hemiplegia": 0.0,
        "Stock_Seq_Hearing_Implant": 0.0,
    }
    rows = []

    for _, row in epi.iterrows():
        year = int(row["Year"])
        event_df = discount_factor(year - 0.5, rate)
        new_mening = float(row["New_Meningitis"])
        new_psev = float(row["New_Pneumo_Sev"])
        new_pnon = float(row["New_Pneumo_NonSev"])
        new_seps = float(row["New_Sepsis"])
        new_aom = float(row["New_AOM"])
        deaths = float(row["New_Deaths_Pneumo"])

        stock_cog = float(row["Stock_Seq_Cognitive"])
        stock_hear = float(row["Stock_Seq_Hearing"])
        stock_epi = float(row["Stock_Seq_Epilepsy"])
        stock_hemi = float(row["Stock_Seq_Hemiplegia"])
        stock_impl = float(row["Stock_Seq_Hearing_Implant"])
        stock_sped = float(row["Stock_Seq_Hearing_SpEd"])

        new_hear = positive_increment(stock_hear, prev_stocks["Stock_Seq_Hearing"])
        new_epi = positive_increment(stock_epi, prev_stocks["Stock_Seq_Epilepsy"])
        new_hemi = positive_increment(stock_hemi, prev_stocks["Stock_Seq_Hemiplegia"])
        new_impl = positive_increment(stock_impl, prev_stocks["Stock_Seq_Hearing_Implant"])
        prev_stocks["Stock_Seq_Hearing"] = stock_hear
        prev_stocks["Stock_Seq_Epilepsy"] = stock_epi
        prev_stocks["Stock_Seq_Hemiplegia"] = stock_hemi
        prev_stocks["Stock_Seq_Hearing_Implant"] = stock_impl

        raw_disease_cost = (
            new_psev * unit_cost["inpatient_pneumonia"]
            + new_pnon * unit_cost["outpatient_pneumonia"]
            + new_mening * unit_cost["inpatient_meningitis"]
            + new_seps * unit_cost["inpatient_sepsis"]
            + new_aom * unit_cost["outpatient_aom"]
        )
        seq_costs = params["costs_sequelae"]
        raw_sequela_cost = (
            stock_cog * seq_costs["cognitive_annual"]
            + stock_sped * seq_costs["special_education_annual"]
            + new_hear * seq_costs["hearing_case"]
            + new_epi * seq_costs["epilepsy_case"]
            + new_hemi * seq_costs["hemiplegia_case"]
            + new_impl * seq_costs["cochlear_implant_case"]
        )
        raw_death_prod_cost = deaths * prod_loss_per_death_event

        def acute_qaly(cases: float, los_days: float, utility: float) -> float:
            return float(cases) * (float(los_days) / 365.0) * (1.0 - float(utility))

        raw_qaly_acute = (
            acute_qaly(new_psev, los["inpatient_pneumonia"], util["inpatient_pneumonia"])
            + acute_qaly(new_pnon, los["outpatient_pneumonia"], util["outpatient_pneumonia"])
            + acute_qaly(new_mening, los["meningitis"], util["meningitis"])
            + acute_qaly(new_seps, los["sepsis"], util["sepsis"])
            + acute_qaly(new_aom, los["aom"], util["aom"])
        )
        raw_qaly_sequela = (
            stock_cog * (1.0 - util["cognitive"])
            + stock_hear * (1.0 - util["hearing"])
            + stock_epi * (1.0 - util["epilepsy"])
            + stock_hemi * (1.0 - util["hemiplegia"])
        )
        raw_qaly_death = deaths * qaly_loss_per_death_event

        if year in vcost.index:
            vrow = vcost.loc[year]
            raw_vaccine_cost = float(vrow["raw_vaccine_program_cost"])
            pv_vaccine_cost = float(vrow["pv_vaccine_program_cost"])
            administered_doses = float(vrow["total_administered_doses"])
        else:
            raw_vaccine_cost = pv_vaccine_cost = administered_doses = 0.0

        pv_disease_cost = raw_disease_cost * event_df
        pv_sequela_cost = raw_sequela_cost * event_df
        pv_death_prod_cost = raw_death_prod_cost * event_df
        pv_total_cost = pv_disease_cost + pv_sequela_cost + pv_death_prod_cost + pv_vaccine_cost

        pv_qaly_acute = raw_qaly_acute * event_df
        pv_qaly_sequela = raw_qaly_sequela * event_df
        pv_qaly_death = raw_qaly_death * event_df
        pv_total_qaly_loss = pv_qaly_acute + pv_qaly_sequela + pv_qaly_death

        rows.append({
            "model_label": model_label,
            "model_family": MODEL_FEATURES[model_label]["model_family"],
            "population_estimand": MODEL_FEATURES[model_label]["population_estimand"],
            "scenario": scenario,
            "vaccine_label": vaccine_label,
            "scenario_variant": "Rollout-to-stable calendar program",
            "vaccine_key": vaccine_key,
            "year": year,
            "New_Meningitis": new_mening,
            "New_Pneumo_Sev": new_psev,
            "New_Pneumo_NonSev": new_pnon,
            "New_Sepsis": new_seps,
            "New_AOM": new_aom,
            "New_Deaths_Pneumo": deaths,
            "cases_total": new_mening + new_psev + new_pnon + new_seps + new_aom,
            "deaths_total": deaths,
            "administered_doses": administered_doses,
            "raw_disease_cost": raw_disease_cost,
            "raw_sequela_cost": raw_sequela_cost,
            "raw_death_prod_cost": raw_death_prod_cost,
            "raw_vaccine_program_cost": raw_vaccine_cost,
            "raw_total_cost": raw_disease_cost + raw_sequela_cost + raw_death_prod_cost + raw_vaccine_cost,
            "pv_disease_cost": pv_disease_cost,
            "pv_sequela_cost": pv_sequela_cost,
            "pv_death_prod_cost": pv_death_prod_cost,
            "pv_vaccine_program_cost": pv_vaccine_cost,
            "pv_total_cost": pv_total_cost,
            "pv_qaly_acute": pv_qaly_acute,
            "pv_qaly_sequela": pv_qaly_sequela,
            "pv_qaly_death": pv_qaly_death,
            "pv_total_qaly_loss": pv_total_qaly_loss,
            "new_seq_hearing_proxy": new_hear,
            "new_seq_epilepsy_proxy": new_epi,
            "new_seq_hemiplegia_proxy": new_hemi,
            "new_seq_implant_proxy": new_impl,
        })
    return pd.DataFrame(rows)


def summarize_calendar_details(details: pd.DataFrame) -> pd.DataFrame:
    grouped = details.groupby(["model_label", "model_family", "population_estimand", "scenario", "vaccine_label", "scenario_variant", "vaccine_key"], dropna=False)
    summary = grouped.agg(
        Total_cases=("cases_total", "sum"),
        Total_deaths=("deaths_total", "sum"),
        Administered_doses=("administered_doses", "sum"),
        Raw_disease_cost=("raw_disease_cost", "sum"),
        Raw_sequela_cost=("raw_sequela_cost", "sum"),
        Raw_death_prod_cost=("raw_death_prod_cost", "sum"),
        Raw_vaccine_program_cost=("raw_vaccine_program_cost", "sum"),
        Raw_total_cost=("raw_total_cost", "sum"),
        PV_disease_cost=("pv_disease_cost", "sum"),
        PV_sequela_cost=("pv_sequela_cost", "sum"),
        PV_death_prod_cost=("pv_death_prod_cost", "sum"),
        PV_vaccine_program_cost=("pv_vaccine_program_cost", "sum"),
        PV_total_cost=("pv_total_cost", "sum"),
        PV_total_qaly_loss=("pv_total_qaly_loss", "sum"),
    ).reset_index()
    return add_baseline_deltas(summary)


def add_baseline_deltas(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model_label, df_model in summary.groupby("model_label"):
        baseline = df_model[df_model["vaccine_label"] == "Baseline"].iloc[0]
        births = MODEL_FEATURES[model_label]["births_denominator"]
        horizon = MODEL_FEATURES[model_label]["time_horizon_years"]
        for _, row in df_model.iterrows():
            out = row.to_dict()
            delta_cost = float(row["PV_total_cost"] - baseline["PV_total_cost"])
            cases_averted = float(baseline["Total_cases"] - row["Total_cases"])
            deaths_averted = float(baseline["Total_deaths"] - row["Total_deaths"])
            qalys_gained = float(baseline["PV_total_qaly_loss"] - row["PV_total_qaly_loss"])
            if row["vaccine_label"] == "Baseline":
                delta_cost = cases_averted = deaths_averted = qalys_gained = 0.0
            out.update({
                "time_horizon_years": horizon,
                "normalization_births": births,
                "Delta_cost_vs_baseline": delta_cost,
                "Cases_averted_vs_baseline": cases_averted,
                "Deaths_averted_vs_baseline": deaths_averted,
                "QALYs_gained_vs_baseline": qalys_gained,
                "USD_per_case_averted": delta_cost / cases_averted if cases_averted else np.nan,
                "USD_per_death_averted": delta_cost / deaths_averted if deaths_averted else np.nan,
                "USD_per_QALY_gained": delta_cost / qalys_gained if qalys_gained else np.nan,
                "Cases_averted_per_100k_births": cases_averted / births * 100000.0,
                "Deaths_averted_per_100k_births": deaths_averted / births * 100000.0,
                "QALYs_gained_per_100k_births": qalys_gained / births * 100000.0,
                "Delta_cost_per_100k_births": delta_cost / births * 100000.0,
            })
            rows.append(out)
    return pd.DataFrame(rows)


def load_original_qaly_summary(path: Path, model_label: str) -> pd.DataFrame:
    raw = pd.read_csv(path)
    rows = []
    for _, row in raw.iterrows():
        scenario_name = str(row["scenario"])
        vaccine_label, scenario_variant = COHORT_SCENARIO_LABELS.get(scenario_name, (scenario_name, "Original cohort output"))
        rows.append({
            "model_label": model_label,
            "model_family": MODEL_FEATURES[model_label]["model_family"],
            "population_estimand": MODEL_FEATURES[model_label]["population_estimand"],
            "scenario": scenario_name,
            "vaccine_label": vaccine_label,
            "scenario_variant": scenario_variant,
            "vaccine_key": row.get("vaccine_used", np.nan),
            "Total_cases": float(row["Total_cases"]),
            "Total_deaths": float(row["Total_deaths"]),
            "Administered_doses": np.nan,
            "Raw_disease_cost": np.nan,
            "Raw_sequela_cost": np.nan,
            "Raw_death_prod_cost": np.nan,
            "Raw_vaccine_program_cost": np.nan,
            "Raw_total_cost": np.nan,
            "PV_disease_cost": np.nan,
            "PV_sequela_cost": np.nan,
            "PV_death_prod_cost": np.nan,
            "PV_vaccine_program_cost": np.nan,
            "PV_total_cost": float(row["PV_total_cost"]),
            "PV_total_qaly_loss": float(row["PV_total_qaly_loss"]),
        })
    return add_baseline_deltas(pd.DataFrame(rows))


def build_calendar_outputs(params: Dict) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    vaccine_costs = build_calendar_vaccine_cost_table(params)
    detail_frames: List[pd.DataFrame] = []
    for model_label, file_map in CALENDAR_YEARLY_FILES.items():
        for scenario, vaccine_label, vaccine_key in CALENDAR_SCENARIOS:
            epi = load_calendar_yearly_epi(file_map[scenario])
            details = compute_calendar_economics(
                model_label=model_label,
                scenario=scenario,
                vaccine_label=vaccine_label,
                vaccine_key=vaccine_key,
                epi=epi,
                vaccine_costs=vaccine_costs,
                params=params,
            )
            detail_frames.append(details)
    calendar_details = pd.concat(detail_frames, ignore_index=True)
    calendar_summary = summarize_calendar_details(calendar_details)
    return vaccine_costs, calendar_details, calendar_summary


def build_cross_refinement_impact(summary: pd.DataFrame) -> pd.DataFrame:
    comparisons = [
        (
            "Cohort stochasticity check",
            MODEL_ORIGINAL_ODE,
            MODEL_ORIGINAL_MARKOV,
            "Original deterministic ODE",
            "Original stochastic Markov",
            ["Dynamic rollout IRR", "Immediate mature IRR"],
        ),
        (
            "Calendar structure and refined vaccine mechanics",
            MODEL_ORIGINAL_ODE,
            MODEL_CALENDAR_ON_SCHEDULE,
            "Original closed birth-cohort ODE",
            "Calendar under-five PDE with on-schedule protection",
            ["Immediate mature IRR"],
        ),
        (
            "Program realism added to calendar PDE",
            MODEL_CALENDAR_ON_SCHEDULE,
            MODEL_CALENDAR_DELAYED,
            "On-schedule protection and IRR-decomposed herd",
            "Delayed dose states, waning, and carriage-informed herd",
            ["Rollout-to-stable calendar program"],
        ),
    ]
    metrics = [
        "Cases_averted_vs_baseline",
        "Deaths_averted_vs_baseline",
        "QALYs_gained_vs_baseline",
        "Delta_cost_vs_baseline",
        "USD_per_QALY_gained",
        "Cases_averted_per_100k_births",
        "Deaths_averted_per_100k_births",
        "QALYs_gained_per_100k_births",
        "Delta_cost_per_100k_births",
    ]
    rows = []
    for comparison_name, reference_model, target_model, reference_label, target_label, variants in comparisons:
        for vaccine_label in VACCINE_ORDER:
            for variant in variants:
                ref = summary[
                    (summary["model_label"] == reference_model)
                    & (summary["vaccine_label"] == vaccine_label)
                    & (summary["scenario_variant"] == variant)
                ]
                tgt = summary[
                    (summary["model_label"] == target_model)
                    & (summary["vaccine_label"] == vaccine_label)
                    & (summary["scenario_variant"] == variant)
                ]
                if ref.empty or tgt.empty:
                    continue
                ref_row = ref.iloc[0]
                tgt_row = tgt.iloc[0]
                out = {
                    "comparison": comparison_name,
                    "vaccine_label": vaccine_label,
                    "scenario_variant": variant,
                    "reference_model": reference_label,
                    "target_model": target_label,
                }
                for metric in metrics:
                    out[f"reference_{metric}"] = ref_row[metric]
                    out[f"target_{metric}"] = tgt_row[metric]
                    out[f"target_minus_reference_{metric}"] = tgt_row[metric] - ref_row[metric]
                rows.append(out)
    return pd.DataFrame(rows)


def write_model_feature_table() -> pd.DataFrame:
    rows = []
    for model_label, meta in MODEL_FEATURES.items():
        row = {"model_label": model_label}
        row.update(meta)
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(CSV_DIR / "model_feature_definitions.csv", index=False)
    return out


def write_assumptions(params: Dict) -> pd.DataFrame:
    rows = [
        ("economic_scope", "mixed comparison with explicit estimand labels", "Original QALY_new rows are closed-cohort 5-year results; calendar PDE rows are 10-year full under-five calendar-program results."),
        ("standardized_reporting", "per 100,000 births or recorded births", "Used to make cohort and calendar estimates easier to compare; raw totals should still be interpreted within each estimand."),
        ("calendar_discount_rate", params["discount_rate"], "Used for newly computed calendar-PDE costs and QALY losses by program year."),
        ("calendar_avg_u5_death_age_years", CALENDAR_CONFIG["AVG_U5_DEATH_AGE_YEARS"], "Neutral under-five average age at pneumococcal death when age-specific death outputs are unavailable."),
        ("on_schedule_calendar_vaccine_costing", "on-time doses at 2, 4, 12 months", "Program-enrolled live births receive scheduled doses if dose month falls inside the 10-year recorded horizon."),
        ("delayed_dose_calendar_vaccine_costing", "model-tracked delayed dose transitions", "Dose counts are re-simulated with explicit dose states, transition probabilities, mortality, and age transport."),
        ("wastage_costing", "procurement only", "Wastage inflates vaccine procurement cost, not service delivery or adverse-event counts."),
        ("range_parser_correction_for_calendar_economics", "positive midpoint parser", "Calendar-PDE economic calculations locally correct hyphenated ranges such as 680-1,360; original imported QALY_new summaries are left as previously generated."),
        ("pcv10sii_price", params["vaccine_prices"]["pcv10sii"], "Retains QALY_new assumption: PCV10-SII uses the PCV10 price row when no separate row is available."),
        ("pcv10_price", params["vaccine_prices"]["pcv10"], "Aligned vaccine price used for PCV10-GSK."),
        ("pcv13_price", params["vaccine_prices"]["pcv13"], "Vaccine price used for PCV13."),
    ]
    out = pd.DataFrame(rows, columns=["parameter", "value", "interpretation"])
    out.to_csv(CSV_DIR / "economic_assumptions.csv", index=False)
    return out


def write_mechanism_mapping() -> pd.DataFrame:
    rows = [
        {
            "mechanism": "Conjugate vaccine immunology",
            "evidence_summary": "PCVs attach pneumococcal capsular polysaccharides to carrier protein, enabling more reliable infant immune response and longer-lasting protection than pure polysaccharide vaccine.",
            "model_implication": "Direct protection should be dose- and product-specific, not an immediate whole-population scalar.",
            "model_refinement_supported": "Explicit program dose states and dose-specific direct VE.",
            "next_data_need": "Vietnam age-at-dose and completion data for 2/4/12 month schedule.",
            "source": "CDC Pink Book, Principles of Vaccination; CDC PCV VIS.",
        },
        {
            "mechanism": "Nasopharyngeal carriage as transmission reservoir",
            "evidence_summary": "S. pneumoniae colonizes the nasopharynx; young children are major reservoirs; PCV reduces vaccine-type carriage.",
            "model_implication": "Herd effect should build through reduced carriage and transmission, not appear at full strength at launch.",
            "model_refinement_supported": "Carriage-informed herd build-up in the calendar PDE.",
            "next_data_need": "Vietnam carriage prevalence and vaccine-type carriage reduction by dose and age.",
            "source": "WHO PCV position paper; PCV carriage systematic reviews; Vietnam PCV10 carriage RCT.",
        },
        {
            "mechanism": "Schedule and booster timing",
            "evidence_summary": "Multiple PCV schedules reduce vaccine-type carriage; Vietnam PCV10 RCT found several schedules reduced vaccine-type carriage, with no clear 2+1 versus 3+0 difference during the first two years.",
            "model_implication": "The model should represent realistic dose timing and allow schedule/catch-up sensitivity rather than assuming instantaneous mature protection.",
            "model_refinement_supported": "Delayed routine vaccination and future catch-up/dropout sensitivity.",
            "next_data_need": "Vietnam national schedule, dropout, delayed vaccination distribution, and catch-up policy.",
            "source": "Smith-Vaughan et al. Vietnam PCV10 carriage RCT; WHO 2025 position paper.",
        },
        {
            "mechanism": "Product serotype coverage and immunogenicity",
            "evidence_summary": "WHO product-choice guidance lists PCV13, PCV10-GSK, and PCV10-SII as WHO-prequalified products; product choice depends on technical and programmatic factors.",
            "model_implication": "PCV10-GSK, PCV13, and PCV10-SII should not share one generic VE/IRR profile.",
            "model_refinement_supported": "Endpoint-specific VE and PCV10-SII GMC-ratio bridge.",
            "next_data_need": "Vietnam serotype distribution and product-specific effectiveness or immunogenicity uncertainty ranges.",
            "source": "WHO PCV product-choice guidance; PCV10-SII immunogenicity evidence bridge.",
        },
        {
            "mechanism": "Serotype replacement",
            "evidence_summary": "After PCV introduction, vaccine-type disease can fall while non-vaccine-type disease rises in some settings.",
            "model_implication": "Long-term model should include vaccine-type and non-vaccine-type compartments or replacement sensitivity.",
            "model_refinement_supported": "Future structural extension beyond the current calendar PDE.",
            "next_data_need": "Vietnam pre/post-PCV serotype-specific IPD, pneumonia, AOM, and carriage data.",
            "source": "CDC Emerging Infectious Diseases serotype replacement analysis; PSERENADE surveillance literature.",
        },
    ]
    out = pd.DataFrame(rows)
    out.to_csv(CSV_DIR / "mechanism_to_model_mapping.csv", index=False)
    return out


def plot_outputs(summary: pd.DataFrame, yearly: pd.DataFrame, refinement: pd.DataFrame) -> List[Path]:
    written: List[Path] = []
    vaccine_summary = summary[summary["vaccine_label"].isin(VACCINE_ORDER)].copy()

    model_order = [
        MODEL_ORIGINAL_ODE,
        MODEL_ORIGINAL_MARKOV,
        MODEL_CALENDAR_ON_SCHEDULE,
        MODEL_CALENDAR_DELAYED,
    ]
    short_labels = {
        MODEL_ORIGINAL_ODE: "Cohort ODE\noriginal economics",
        MODEL_ORIGINAL_MARKOV: "Cohort Markov\noriginal economics",
        MODEL_CALENDAR_ON_SCHEDULE: "Calendar PDE\non-schedule protection",
        MODEL_CALENDAR_DELAYED: "Calendar PDE\ndelayed doses + carriage herd",
    }

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(model_order))
    width = 0.24
    for i, vaccine in enumerate(VACCINE_ORDER):
        values = []
        for model in model_order:
            df = vaccine_summary[(vaccine_summary["model_label"] == model) & (vaccine_summary["vaccine_label"] == vaccine)]
            if model in [MODEL_ORIGINAL_ODE, MODEL_ORIGINAL_MARKOV]:
                df = df[df["scenario_variant"] == "Immediate mature IRR"]
            else:
                df = df[df["scenario_variant"] == "Rollout-to-stable calendar program"]
            values.append(float(df.iloc[0]["QALYs_gained_per_100k_births"]) if not df.empty else np.nan)
        ax.bar(x + (i - 1) * width, values, width=width, label=vaccine)
    ax.set_xticks(x)
    ax.set_xticklabels([short_labels[m] for m in model_order])
    ax.set_ylabel("PV QALYs gained per 100,000 births/recorded births")
    ax.set_title("QALY impact across model refinements, standardized")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    path = FIGURE_DIR / "qalys_gained_standardized_by_model_feature.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)

    fig, ax = plt.subplots(figsize=(12, 6))
    for i, vaccine in enumerate(VACCINE_ORDER):
        values = []
        for model in model_order:
            df = vaccine_summary[(vaccine_summary["model_label"] == model) & (vaccine_summary["vaccine_label"] == vaccine)]
            if model in [MODEL_ORIGINAL_ODE, MODEL_ORIGINAL_MARKOV]:
                df = df[df["scenario_variant"] == "Immediate mature IRR"]
            else:
                df = df[df["scenario_variant"] == "Rollout-to-stable calendar program"]
            values.append(float(df.iloc[0]["USD_per_QALY_gained"]) if not df.empty else np.nan)
        ax.bar(x + (i - 1) * width, values, width=width, label=vaccine)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([short_labels[m] for m in model_order])
    ax.set_ylabel("USD per QALY gained")
    ax.set_title("Cost-effectiveness across model refinements")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    path = FIGURE_DIR / "icer_by_model_feature.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)

    target_refinement = refinement[
        (refinement["comparison"] == "Program realism added to calendar PDE")
        & (refinement["scenario_variant"] == "Rollout-to-stable calendar program")
    ].copy()
    fig, ax = plt.subplots(figsize=(9, 5.4))
    ax.bar(
        target_refinement["vaccine_label"],
        target_refinement["target_minus_reference_QALYs_gained_vs_baseline"],
    )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Change in PV QALYs gained")
    ax.set_title("Effect of delayed-dose and carriage-informed assumptions")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    path = FIGURE_DIR / "delayed_dose_carriage_refinement_qaly_change.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)

    fig, ax = plt.subplots(figsize=(11, 5.8))
    for model, linestyle in [
        (MODEL_CALENDAR_ON_SCHEDULE, "--"),
        (MODEL_CALENDAR_DELAYED, "-"),
    ]:
        for vaccine in VACCINE_ORDER:
            df = yearly[(yearly["model_label"] == model) & (yearly["vaccine_label"] == vaccine)]
            ax.plot(df["year"], df["qalys_gained"], linestyle=linestyle, marker="o", markersize=3, label=f"{vaccine} | {short_labels[model].replace(chr(10), ' ')}")
    ax.set_xlabel("Program year")
    ax.set_ylabel("PV QALYs gained in year")
    ax.set_title("Annual QALYs gained in calendar under-five PDE models")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    path = FIGURE_DIR / "annual_qalys_gained_calendar_pde.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)

    fig, ax = plt.subplots(figsize=(11, 5.8))
    for model, linestyle in [
        (MODEL_CALENDAR_ON_SCHEDULE, "--"),
        (MODEL_CALENDAR_DELAYED, "-"),
    ]:
        for vaccine in VACCINE_ORDER:
            df = yearly[(yearly["model_label"] == model) & (yearly["vaccine_label"] == vaccine)]
            ax.plot(df["year"], df["pv_vaccine_program_cost"], linestyle=linestyle, marker="o", markersize=3, label=f"{vaccine} | {short_labels[model].replace(chr(10), ' ')}")
    ax.set_xlabel("Program year")
    ax.set_ylabel("PV vaccine program cost, USD")
    ax.set_title("Annual vaccine budget impact in calendar under-five PDE models")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    path = FIGURE_DIR / "budget_impact_calendar_pde.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)
    return written


def build_calendar_yearly_impact(details: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model_label, df_model in details.groupby("model_label"):
        baseline = df_model[df_model["vaccine_label"] == "Baseline"].set_index("year")
        for vaccine in VACCINE_ORDER:
            target = df_model[df_model["vaccine_label"] == vaccine].set_index("year")
            for year in sorted(target.index):
                b = baseline.loc[year]
                r = target.loc[year]
                rows.append({
                    "model_label": model_label,
                    "vaccine_label": vaccine,
                    "year": int(year),
                    "cases_averted": float(b["cases_total"] - r["cases_total"]),
                    "deaths_averted": float(b["deaths_total"] - r["deaths_total"]),
                    "qalys_gained": float(b["pv_total_qaly_loss"] - r["pv_total_qaly_loss"]),
                    "delta_pv_cost": float(r["pv_total_cost"] - b["pv_total_cost"]),
                    "pv_vaccine_program_cost": float(r["pv_vaccine_program_cost"]),
                    "pv_disease_cost_saved": float(b["pv_disease_cost"] - r["pv_disease_cost"]),
                    "pv_sequela_cost_saved": float(b["pv_sequela_cost"] - r["pv_sequela_cost"]),
                    "pv_death_prod_cost_saved": float(b["pv_death_prod_cost"] - r["pv_death_prod_cost"]),
                })
    return pd.DataFrame(rows)


def write_markdown_summary(summary: pd.DataFrame, refinement: pd.DataFrame, figures: Iterable[Path]) -> None:
    rows = summary[summary["vaccine_label"].isin(VACCINE_ORDER)].copy()
    compact = rows[[
        "model_label",
        "vaccine_label",
        "scenario_variant",
        "Cases_averted_vs_baseline",
        "Deaths_averted_vs_baseline",
        "QALYs_gained_vs_baseline",
        "QALYs_gained_per_100k_births",
        "Delta_cost_vs_baseline",
        "USD_per_QALY_gained",
    ]]
    ref_cols = [
        "comparison",
        "vaccine_label",
        "scenario_variant",
        "reference_model",
        "target_model",
        "target_minus_reference_QALYs_gained_vs_baseline",
        "target_minus_reference_QALYs_gained_per_100k_births",
        "target_minus_reference_Delta_cost_vs_baseline",
        "target_minus_reference_USD_per_QALY_gained",
    ]
    lines = [
        "# QALY and Economic Impact Across Model Refinements",
        "",
        "This summary uses descriptive model-feature labels rather than development chronology.",
        "",
        "## Main Comparison",
        "",
        compact.to_markdown(index=False, floatfmt=".2f"),
        "",
        "## Refinement Impact",
        "",
        refinement[ref_cols].to_markdown(index=False, floatfmt=".2f"),
        "",
        "## Figure Files",
        "",
    ]
    for fig in figures:
        lines.append(f"- {fig.name}")
    lines.append("")
    SUMMARY_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_output_dirs()
    require_inputs()
    params = build_calendar_params()
    write_model_feature_table()
    write_assumptions(params)
    write_mechanism_mapping()

    vaccine_costs, calendar_details, calendar_summary = build_calendar_outputs(params)
    original_ode = load_original_qaly_summary(RESULTS_CSV_DIR / "results_summary_ODE.csv", MODEL_ORIGINAL_ODE)
    original_markov = load_original_qaly_summary(RESULTS_CSV_DIR / "results_summary_Markov.csv", MODEL_ORIGINAL_MARKOV)

    all_summary = pd.concat([original_ode, original_markov, calendar_summary], ignore_index=True)
    refinement = build_cross_refinement_impact(all_summary)
    yearly_impact = build_calendar_yearly_impact(calendar_details)

    vaccine_costs.to_csv(CSV_DIR / "calendar_pde_vaccine_dose_costs.csv", index=False)
    calendar_details.to_csv(CSV_DIR / "calendar_pde_yearly_economic_details.csv", index=False)
    yearly_impact.to_csv(CSV_DIR / "calendar_pde_yearly_impact_vs_baseline.csv", index=False)
    all_summary.to_csv(CSV_DIR / "qaly_summary_all_model_refinements.csv", index=False)
    refinement.to_csv(CSV_DIR / "refinement_impact_between_model_features.csv", index=False)
    original_ode.to_csv(CSV_DIR / "imported_qaly_new_ode_summary.csv", index=False)
    original_markov.to_csv(CSV_DIR / "imported_qaly_new_markov_summary.csv", index=False)

    figures = plot_outputs(all_summary, yearly_impact, refinement)
    write_markdown_summary(all_summary, refinement, figures)

    print("Economic refinement comparison complete.")
    print(f"Summary: {CSV_DIR / 'qaly_summary_all_model_refinements.csv'}")
    print(f"Refinement impact: {CSV_DIR / 'refinement_impact_between_model_features.csv'}")
    print(f"Methods summary: {SUMMARY_MD}")
    print("Figures:")
    for fig in figures:
        print(f"  {fig}")
    print("\nKey standardized results:")
    cols = [
        "model_label",
        "vaccine_label",
        "scenario_variant",
        "QALYs_gained_per_100k_births",
        "USD_per_QALY_gained",
    ]
    display = all_summary[all_summary["vaccine_label"].isin(VACCINE_ORDER)][cols]
    print(display.to_string(index=False))


if __name__ == "__main__":
    main()
