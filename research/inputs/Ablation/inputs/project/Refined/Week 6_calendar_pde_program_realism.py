"""
Calendar PDE model with programme-realism vaccination timing.

This file builds on Week 5_calendar_pde_refined.py without modifying it. The
on-schedule calendar PDE already used a solver-based age-time framework and
endpoint-specific direct VE. This programme-realism layer adds the next
layer of program realism:

1. seasonal monthly birth shares, used as a transparent assumption until
   verified Vietnam monthly births are available;
2. explicit routine-program dose states: enrolled-no-dose, dose 1, dose 2,
   and booster dose;
3. delayed routine vaccination around the 2-, 4-, and 12-month targets;
4. dose-specific direct VE, with waning after the primary series and partial
   restoration after the booster;
5. herd protection that builds through a carriage index rather than appearing
   immediately at program launch.

The model still simulates ten recorded years by default. It also saves
on-schedule reference outputs and comparison line plots so the effect of the new
assumptions is visible on both short and long time scales.
"""

import copy
import importlib.util
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "Outputs" / "Week 6_calendar_pde_program_realism"
DEFAULT_RECORD_MONTHS = 120


def _load_peer(filename: str, module_name: str):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _resolve_outdir(outdir: Optional[str]) -> Path:
    if outdir is None:
        return DEFAULT_OUTPUT_DIR
    path = Path(outdir)
    if not path.is_absolute():
        path = SCRIPT_DIR / path
    return path


w5 = _load_peer("Week 5_calendar_pde_refined.py", "week5_calendar_pde_refined_reference")
base = w5.base

FLOW_NAMES = list(w5.FLOW_NAMES)
OUTPUT_COLUMNS = list(w5.OUTPUT_COLUMNS)
SEQ_KEYS = list(w5.SEQ_KEYS)


# Alive state groups by age month.
GROUP_NONPROGRAM = 0
GROUP_PROGRAM_D0 = 1
GROUP_PROGRAM_D1 = 2
GROUP_PROGRAM_D2 = 3
GROUP_PROGRAM_D3 = 4
GROUP_NAMES = [
    "NonProgram",
    "Program_Dose0",
    "Program_Dose1",
    "Program_Dose2",
    "Program_Booster",
]
N_GROUPS = len(GROUP_NAMES)


# A mild seasonality assumption. This is deliberately documented as an
# assumption, not observed Vietnam data. Replace this vector when national
# monthly births are available.
ASSUMED_SEASONAL_BIRTH_SHARES = [
    0.080, 0.074, 0.078, 0.082, 0.086, 0.091,
    0.094, 0.092, 0.087, 0.083, 0.078, 0.075,
]


params_vietnam_micro = copy.deepcopy(w5.params_vietnam_micro)
params_vietnam_micro.update({
    "record_months": DEFAULT_RECORD_MONTHS,
    "birth_month_shares": ASSUMED_SEASONAL_BIRTH_SHARES,

    # Routine vaccination around 2, 4, and 12 months.
    # These are monthly transition probabilities applied to children who have
    # reached the relevant dose state but have not yet received the next dose.
    # They intentionally create delayed and incomplete vaccination rather than
    # assuming all enrolled children are vaccinated exactly on schedule.
    "routine_dose1_transition_probs": {2: 0.65, 3: 0.57, 4: 0.67, 5: 0.60},
    "routine_dose2_transition_probs": {4: 0.70, 5: 0.45, 6: 0.30, 7: 0.20, 8: 0.10},
    "routine_booster_transition_probs": {12: 0.55, 13: 0.35, 14: 0.25, 15: 0.15, 16: 0.10},

    # Optional one-year catch-up campaign. It is off by default because national
    # introduction policy has not been specified; keep it here for sensitivity.
    "catchup_enabled": False,
    "catchup_campaign_months": 12,
    "catchup_7_11_to_dose1_prob": 0.12,
    "catchup_12_23_to_dose1_prob": 0.08,

    # Direct VE by dose state. Dose 1 is partial; dose 2 is the on-schedule full
    # primary-series direct VE; the booster restores and slightly boosts VE but
    # is still capped below one.
    "dose1_relative_ve": 0.45,
    "dose2_relative_ve": 1.00,
    "booster_relative_ve": 1.10,
    "max_direct_ve": 0.95,

    # Waning after the primary series and booster. The IPD paper provided by the
    # user notes lower protection with increasing time since vaccination; the
    # values below are transparent structural assumptions for sensitivity.
    "primary_waning_by_age": {24: 0.85, 36: 0.70},
    "booster_waning_by_age": {24: 0.95, 36: 0.80, 48: 0.70},

    # Carriage-informed herd build-up. The carriage RCT supports herd protection
    # through reduced vaccine-type carriage. These values represent relative VT
    # carriage reduction by dose state and are used only to scale how quickly the
    # mature herd multiplier is approached.
    "carriage_reduction_by_group": {
        "Program_Dose0": 0.00,
        "Program_Dose1": 0.15,
        "Program_Dose2": 0.34,
        "Program_Booster": 0.40,
    },
    "carriage_age_weights": {
        "0_5": 0.00,
        "6_17": 1.00,
        "18_23": 0.80,
        "24_59": 0.35,
    },
    "mature_effective_dose_share": 0.88,
    "carriage_to_herd_fraction": 1.00,
})


def _as_int_float_dict(value: Dict) -> Dict[int, float]:
    return {int(k): float(v) for k, v in dict(value).items()}


def _age_weight_for_carriage(age_month: int, p: Dict) -> float:
    weights = p.get("carriage_age_weights", {})
    a = int(age_month)
    if 0 <= a <= 5:
        return float(weights.get("0_5", 0.0))
    if 6 <= a <= 17:
        return float(weights.get("6_17", 1.0))
    if 18 <= a <= 23:
        return float(weights.get("18_23", 0.8))
    return float(weights.get("24_59", 0.35))


def dose_relative_ve_for_group(p: Dict, group: int, age_month: int) -> float:
    if group == GROUP_PROGRAM_D1:
        return float(p.get("dose1_relative_ve", 0.45))
    if group == GROUP_PROGRAM_D2:
        rel = float(p.get("dose2_relative_ve", 1.0))
        for threshold, value in sorted(_as_int_float_dict(p.get("primary_waning_by_age", {})).items()):
            if int(age_month) >= threshold:
                rel = float(value)
        return rel
    if group == GROUP_PROGRAM_D3:
        rel = float(p.get("booster_relative_ve", 1.10))
        for threshold, value in sorted(_as_int_float_dict(p.get("booster_waning_by_age", {})).items()):
            if int(age_month) >= threshold:
                rel = float(value)
        return rel
    return 0.0


def dose_endpoint_ve(p: Dict, vaccine_type: Optional[str], endpoint: str, group: int, age_month: int) -> float:
    if vaccine_type is None:
        return 0.0
    base_ve = w5.endpoint_direct_ve(p, vaccine_type, endpoint)
    rel = dose_relative_ve_for_group(p, group, age_month)
    return max(0.0, min(float(p.get("max_direct_ve", 0.95)), base_ve * rel))


def _transition_by_age(A: np.ndarray, source_group: int, target_group: int, probs: Dict[int, float]) -> None:
    for age, prob in _as_int_float_dict(probs).items():
        if age < 0 or age >= A.shape[1]:
            continue
        p_move = max(0.0, min(1.0, float(prob)))
        moved = A[source_group, age] * p_move
        if moved <= 0.0:
            continue
        A[source_group, age] -= moved
        A[target_group, age] += moved


def apply_vaccination_events(
    A: np.ndarray,
    p: Dict,
    calendar_month: int,
    intervention: bool,
    vaccine_type: Optional[str],
) -> None:
    if not intervention or vaccine_type is None:
        return

    # Optional catch-up is deliberately conservative: it moves a fraction of
    # non-program children into dose 1 during the first program year only.
    if bool(p.get("catchup_enabled", False)) and int(calendar_month) < int(p.get("catchup_campaign_months", 12)):
        p_7_11 = max(0.0, min(1.0, float(p.get("catchup_7_11_to_dose1_prob", 0.0))))
        p_12_23 = max(0.0, min(1.0, float(p.get("catchup_12_23_to_dose1_prob", 0.0))))
        for age in range(7, min(12, A.shape[1])):
            moved = A[GROUP_NONPROGRAM, age] * p_7_11
            A[GROUP_NONPROGRAM, age] -= moved
            A[GROUP_PROGRAM_D1, age] += moved
        for age in range(12, min(24, A.shape[1])):
            moved = A[GROUP_NONPROGRAM, age] * p_12_23
            A[GROUP_NONPROGRAM, age] -= moved
            A[GROUP_PROGRAM_D1, age] += moved

    # Apply later doses first so a child cannot receive two routine doses in the
    # same monthly step.
    _transition_by_age(A, GROUP_PROGRAM_D2, GROUP_PROGRAM_D3, p.get("routine_booster_transition_probs", {}))
    _transition_by_age(A, GROUP_PROGRAM_D1, GROUP_PROGRAM_D2, p.get("routine_dose2_transition_probs", {}))
    _transition_by_age(A, GROUP_PROGRAM_D0, GROUP_PROGRAM_D1, p.get("routine_dose1_transition_probs", {}))


def carriage_index(A: np.ndarray, p: Dict) -> Tuple[float, float]:
    reductions = p.get("carriage_reduction_by_group", {})
    weighted_alive = 0.0
    weighted_reduced = 0.0
    for age in range(A.shape[1]):
        w_age = _age_weight_for_carriage(age, p)
        if w_age <= 0.0:
            continue
        alive_age = float(A[:, age].sum())
        weighted_alive += w_age * alive_age
        for group, name in enumerate(GROUP_NAMES):
            red = float(reductions.get(name, 0.0))
            weighted_reduced += w_age * float(A[group, age]) * red
    index = weighted_reduced / weighted_alive if weighted_alive > 0.0 else 0.0
    reference = (
        float(p.get("vaccine_coverage", 0.0))
        * float(p.get("mature_effective_dose_share", 0.88))
        * float(reductions.get("Program_Dose2", 0.34))
    )
    relative = max(0.0, min(1.0, index / reference)) if reference > 0.0 else 0.0
    return float(index), float(relative)


def herd_multiplier_for_endpoint(
    p: Dict,
    vaccine_type: Optional[str],
    endpoint: str,
    year_idx: int,
    A: np.ndarray,
) -> Tuple[float, float, float]:
    if vaccine_type is None:
        return 1.0, 0.0, 0.0
    series = p.get("vaccine_irrs", {}).get(vaccine_type)
    target_irr = float(base.week4._get_population_irr_for_year(series, int(year_idx)))

    base_ve = w5.endpoint_direct_ve(p, vaccine_type, endpoint)
    mature_direct_ref = (
        float(p.get("vaccine_coverage", 0.0))
        * float(p.get("mature_effective_dose_share", 0.88))
        * float(base_ve)
    )
    mature_direct_ref = max(0.0, min(0.95, mature_direct_ref))
    denom = max(1e-12, 1.0 - mature_direct_ref)
    mature_herd = max(0.0, min(1.0, target_irr / denom))

    _, relative_carriage = carriage_index(A, p)
    herd_fraction = max(0.0, min(1.0, float(p.get("carriage_to_herd_fraction", 1.0))))
    actual_herd = 1.0 - (1.0 - mature_herd) * relative_carriage * herd_fraction
    actual_herd = max(mature_herd, min(1.0, actual_herd))
    return float(actual_herd), float(mature_herd), float(relative_carriage)


def build_week6_month_coefficients(
    p: Dict,
    calendar_month: int,
    A_start: np.ndarray,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    n_age = int(p.get("age_bins", 60))
    year_idx = int(calendar_month) // 12
    prop_sev = float(p["prop_pneumo_severe"])
    cfr_men = float(p["cfr_meningitis"])
    cfr_pneu = float(p["cfr_pneumo_severe"])
    cfr_sep = float(p["cfr_sepsis"])
    lam_bg = base._prob_to_hazard(base.get_background_prob_u5_only(p))

    if intervention and vaccine_type is not None:
        herd_ipd, _, _ = herd_multiplier_for_endpoint(p, vaccine_type, "ipd", year_idx, A_start)
        herd_pneu, _, _ = herd_multiplier_for_endpoint(p, vaccine_type, "pneumonia", year_idx, A_start)
        herd_aom, _, _ = herd_multiplier_for_endpoint(p, vaccine_type, "aom", year_idx, A_start)
    else:
        herd_ipd = herd_pneu = herd_aom = 1.0

    coeffs = {name: np.zeros((N_GROUPS, n_age), dtype=float) for name in [
        "lam_men", "lam_pneu", "lam_sep", "lam_aom", "mu_dis", "mu_all",
    ]}

    for age_month in range(n_age):
        base_prob_men, base_prob_pneu, base_prob_sep, base_prob_aom = base._disease_probabilities_for_age(p, age_month)
        for group in range(N_GROUPS):
            if intervention and vaccine_type is not None:
                ve_ipd = dose_endpoint_ve(p, vaccine_type, "ipd", group, age_month)
                ve_pneu = dose_endpoint_ve(p, vaccine_type, "pneumonia", group, age_month)
                ve_aom = dose_endpoint_ve(p, vaccine_type, "aom", group, age_month)
            else:
                ve_ipd = ve_pneu = ve_aom = 0.0

            p_men = base._clip01(base_prob_men * herd_ipd * (1.0 - ve_ipd))
            p_pneu = base._clip01(base_prob_pneu * herd_pneu * (1.0 - ve_pneu))
            p_sep = base._clip01(base_prob_sep * herd_ipd * (1.0 - ve_ipd))
            p_aom = base._clip01(base_prob_aom * herd_aom * (1.0 - ve_aom))

            coeffs["lam_men"][group, age_month] = base._prob_to_hazard(p_men)
            coeffs["lam_pneu"][group, age_month] = base._prob_to_hazard(p_pneu)
            coeffs["lam_sep"][group, age_month] = base._prob_to_hazard(p_sep)
            coeffs["lam_aom"][group, age_month] = base._prob_to_hazard(p_aom)

            d1 = base._clip01(p_men * cfr_men)
            d2 = base._clip01(p_pneu * prop_sev * cfr_pneu)
            d3 = base._clip01(p_sep * cfr_sep)
            coeffs["mu_dis"][group, age_month] = base._prob_to_hazard(d1) + base._prob_to_hazard(d2) + base._prob_to_hazard(d3)
            coeffs["mu_all"][group, age_month] = lam_bg + coeffs["mu_dis"][group, age_month]

    coeffs["prop_sev"] = np.array([prop_sev], dtype=float)
    return coeffs


def week6_pde_reaction_rhs(_tau: float, y: np.ndarray, coeffs: Dict[str, np.ndarray]) -> np.ndarray:
    n_age = coeffs["mu_all"].shape[1]
    alive_len = N_GROUPS * n_age
    A = np.maximum(y[:alive_len], 0.0).reshape((N_GROUPS, n_age))

    dy = np.zeros_like(y, dtype=float)
    dA = -coeffs["mu_all"] * A
    dy[:alive_len] = dA.reshape(alive_len)

    prop_sev = float(coeffs["prop_sev"][0])
    inc_men = float((coeffs["lam_men"] * A).sum())
    inc_pneu = float((coeffs["lam_pneu"] * A).sum())
    inc_sep = float((coeffs["lam_sep"] * A).sum())
    inc_aom = float((coeffs["lam_aom"] * A).sum())
    inc_death = float((coeffs["mu_dis"] * A).sum())

    dy[alive_len + 0] = inc_men
    dy[alive_len + 1] = prop_sev * inc_pneu
    dy[alive_len + 2] = (1.0 - prop_sev) * inc_pneu
    dy[alive_len + 3] = inc_sep
    dy[alive_len + 4] = inc_aom
    dy[alive_len + 5] = inc_death
    return dy


def solve_one_month_week6(
    A: np.ndarray,
    S: Dict[str, float],
    calendar_month: int,
    p: Dict,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
    rtol: float = 1e-8,
    atol: float = 1e-7,
) -> Tuple[np.ndarray, Dict[str, float], Dict[str, float], float]:
    coeffs = build_week6_month_coefficients(p, calendar_month, A, intervention, vaccine_type)
    n_age = A.shape[1]
    alive_len = N_GROUPS * n_age
    y0 = np.concatenate([A.astype(float).reshape(alive_len), np.zeros(len(FLOW_NAMES), dtype=float)])

    sol = solve_ivp(
        fun=lambda tau, yy: week6_pde_reaction_rhs(tau, yy, coeffs),
        t_span=(0.0, 1.0),
        y0=y0,
        t_eval=[1.0],
        method="RK45",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed at month {calendar_month}: {sol.message}")

    y1 = np.maximum(sol.y[:, -1], 0.0)
    A_after = y1[:alive_len].reshape((N_GROUPS, n_age))
    flow_values = y1[alive_len:alive_len + len(FLOW_NAMES)]
    flows = {name: float(value) for name, value in zip(FLOW_NAMES, flow_values)}

    alive_before = float(A.sum())
    alive_after = float(A_after.sum())
    non_neonatal_all_cause_deaths = max(0.0, alive_before - alive_after)
    surv_seq = max(0.0, min(1.0, alive_after / alive_before)) if alive_before > 0.0 else 1.0
    men_survivor_flow = (1.0 - float(p["cfr_meningitis"])) * flows["New_Meningitis"]

    src_cog = float(p["prob_cog"]) * men_survivor_flow
    src_hear = float(p["prob_hear"]) * men_survivor_flow
    src_epi = float(p["prob_epi"]) * men_survivor_flow
    src_hemi = float(p["prob_hemi"]) * men_survivor_flow
    src_impl = float(p.get("prob_cochlear_implant", 0.0)) * src_hear
    src_sped = float(p.get("prob_special_education", 0.0)) * src_hear

    S_after = base.copy_sequela_stock(S)
    S_after["Stock_Seq_Cognitive"] = S_after["Stock_Seq_Cognitive"] * surv_seq + src_cog
    S_after["Stock_Seq_Hearing"] = S_after["Stock_Seq_Hearing"] * surv_seq + src_hear
    S_after["Stock_Seq_Epilepsy"] = S_after["Stock_Seq_Epilepsy"] * surv_seq + src_epi
    S_after["Stock_Seq_Hemiplegia"] = S_after["Stock_Seq_Hemiplegia"] * surv_seq + src_hemi
    S_after["Stock_Seq_Hearing_Implant"] = S_after["Stock_Seq_Hearing_Implant"] * surv_seq + src_impl
    S_after["Stock_Seq_Hearing_SpEd"] = S_after["Stock_Seq_Hearing_SpEd"] * surv_seq + src_sped
    return A_after, S_after, flows, non_neonatal_all_cause_deaths


def age_shift_week6(A_after: np.ndarray) -> Tuple[np.ndarray, float]:
    A_next = np.zeros_like(A_after, dtype=float)
    A_next[:, 1:] = A_after[:, :-1]
    aged_out = float(A_after[:, -1].sum())
    return A_next, aged_out


def add_record_births_week6(
    A: np.ndarray,
    gross_births: float,
    p: Dict,
    intervention: bool,
    vaccine_type: Optional[str],
) -> Tuple[float, float]:
    alive_births = float(gross_births) * base.neonatal_survival_multiplier(p)
    neonatal_deaths = float(gross_births) - alive_births
    if intervention and vaccine_type is not None:
        cov = max(0.0, min(1.0, float(p.get("vaccine_coverage", 0.0))))
        A[GROUP_PROGRAM_D0, 0] += alive_births * cov
        A[GROUP_NONPROGRAM, 0] += alive_births * (1.0 - cov)
    else:
        A[GROUP_NONPROGRAM, 0] += alive_births
    return alive_births, neonatal_deaths


def build_week6_initial_population(
    p: Dict,
    shares: np.ndarray,
    scenario_name: str,
    retain_preexisting_sequelae: bool = False,
) -> Tuple[np.ndarray, Dict[str, float], Dict[str, float]]:
    n_age = int(p.get("age_bins", 60))
    warmup_months = int(p.get("warmup_months", 60))
    births_annual = float(p["birth_cohort"])
    surv0 = base.neonatal_survival_multiplier(p)

    A = np.zeros((N_GROUPS, n_age), dtype=float)
    S_warm = base.empty_sequela_stock()

    for warm_month in range(warmup_months):
        gross_births = births_annual * float(shares[warm_month % 12])
        A[GROUP_NONPROGRAM, 0] += gross_births * surv0
        A_after, S_warm, _, _ = solve_one_month_week6(A, S_warm, warm_month, p, intervention=False, vaccine_type=None)
        A, _ = age_shift_week6(A_after)

    S_initial = base.copy_sequela_stock(S_warm) if retain_preexisting_sequelae else base.empty_sequela_stock()
    warm_meta = {
        "scenario_name": scenario_name,
        "warmup_months": warmup_months,
        "warmup_mode": "baseline_only_no_vaccination",
        "warmup_alive_after_age_shift": float(A.sum()),
        "warmup_preexisting_sequelae_total": float(sum(S_warm.values())),
        "retain_preexisting_sequelae": bool(retain_preexisting_sequelae),
    }
    return A, S_initial, warm_meta


def _dose_diagnostics(A: np.ndarray, p: Dict) -> Dict[str, float]:
    alive = float(A.sum())
    diag = {
        "NonProgram_Alive": float(A[GROUP_NONPROGRAM].sum()),
        "Program_Dose0_Alive": float(A[GROUP_PROGRAM_D0].sum()),
        "Program_Dose1_Alive": float(A[GROUP_PROGRAM_D1].sum()),
        "Program_Dose2_Alive": float(A[GROUP_PROGRAM_D2].sum()),
        "Program_Booster_Alive": float(A[GROUP_PROGRAM_D3].sum()),
    }
    diag["Program_Alive"] = diag["Program_Dose0_Alive"] + diag["Program_Dose1_Alive"] + diag["Program_Dose2_Alive"] + diag["Program_Booster_Alive"]
    diag["AtLeastDose1_Share"] = (diag["Program_Dose1_Alive"] + diag["Program_Dose2_Alive"] + diag["Program_Booster_Alive"]) / alive if alive else 0.0
    diag["AtLeastDose2_Share"] = (diag["Program_Dose2_Alive"] + diag["Program_Booster_Alive"]) / alive if alive else 0.0
    diag["Booster_Share"] = diag["Program_Booster_Alive"] / alive if alive else 0.0
    ci, cr = carriage_index(A, p)
    diag["Carriage_Index"] = ci
    diag["Relative_Carriage_Index"] = cr
    return diag


def run_pde_calendar_week6_monthly(
    p: Dict,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
    scenario_name: str = "scenario",
    retain_preexisting_sequelae: bool = False,
) -> pd.DataFrame:
    p = copy.deepcopy(p)
    record_months = int(p.get("record_months", DEFAULT_RECORD_MONTHS))
    shares = base.normalize_birth_shares(p.get("birth_month_shares", ASSUMED_SEASONAL_BIRTH_SHARES))
    births_annual = float(p["birth_cohort"])

    A, S, _ = build_week6_initial_population(
        p=p,
        shares=shares,
        scenario_name=scenario_name,
        retain_preexisting_sequelae=retain_preexisting_sequelae,
    )

    cum_men = cum_psev = cum_pnon = cum_sep = cum_aom = cum_ddeath = 0.0
    cum_all_cause_deaths = 0.0
    history = []

    for calendar_month in range(record_months):
        gross_births = births_annual * float(shares[calendar_month % 12])
        _, neonatal_deaths = add_record_births_week6(A, gross_births, p, intervention, vaccine_type)
        cum_all_cause_deaths += neonatal_deaths

        apply_vaccination_events(A, p, calendar_month, intervention, vaccine_type)
        start_diag = _dose_diagnostics(A, p)

        A_after, S, flows, non_neonatal_deaths = solve_one_month_week6(
            A, S, calendar_month, p, intervention=intervention, vaccine_type=vaccine_type
        )
        cum_all_cause_deaths += non_neonatal_deaths

        stock_alive = float(A_after.sum())
        cum_men += flows["New_Meningitis"]
        cum_psev += flows["New_Pneumo_Sev"]
        cum_pnon += flows["New_Pneumo_NonSev"]
        cum_sep += flows["New_Sepsis"]
        cum_aom += flows["New_AOM"]
        cum_ddeath += flows["New_Deaths_Pneumo"]

        row = {
            "Month": calendar_month,
            "Year": calendar_month // 12 + 1,
            "New_Meningitis": flows["New_Meningitis"],
            "New_Pneumo_Sev": flows["New_Pneumo_Sev"],
            "New_Pneumo_NonSev": flows["New_Pneumo_NonSev"],
            "New_Sepsis": flows["New_Sepsis"],
            "New_AOM": flows["New_AOM"],
            "New_Deaths_Pneumo": flows["New_Deaths_Pneumo"],
            "Stock_Alive": stock_alive,
            "Stock_Dead": cum_all_cause_deaths,
            "Stock_New_Meningitis": cum_men,
            "Stock_New_Pneumo_Sev": cum_psev,
            "Stock_New_Pneumo_NonSev": cum_pnon,
            "Stock_New_Sepsis": cum_sep,
            "Stock_New_AOM": cum_aom,
            "Stock_New_Deaths_Pneumo": cum_ddeath,
            "Stock_Seq_Cognitive": S["Stock_Seq_Cognitive"],
            "Stock_Seq_Hearing": S["Stock_Seq_Hearing"],
            "Stock_Seq_Epilepsy": S["Stock_Seq_Epilepsy"],
            "Stock_Seq_Hemiplegia": S["Stock_Seq_Hemiplegia"],
            "Stock_Seq_Hearing_Implant": S["Stock_Seq_Hearing_Implant"],
            "Stock_Seq_Hearing_SpEd": S["Stock_Seq_Hearing_SpEd"],
            "Monthly_Total_Cases": flows["New_Meningitis"] + flows["New_Pneumo_Sev"] + flows["New_Pneumo_NonSev"] + flows["New_Sepsis"] + flows["New_AOM"],
        }
        row.update(start_diag)
        history.append(row)

        A, _ = age_shift_week6(A_after)

    return pd.DataFrame(history)


def _scenario_specs(record_months: int = DEFAULT_RECORD_MONTHS) -> List[Tuple[str, Dict, bool, Optional[str], str]]:
    p = copy.deepcopy(params_vietnam_micro)
    p["record_months"] = int(record_months)
    return [
        ("baseline", p, False, None, "calendar_pde_program_realism_baseline.csv"),
        ("pcv10_gsk_rollout_to_stable", p, True, "PCV10-GSK", "calendar_pde_program_realism_pcv10_gsk_rollout_to_stable.csv"),
        ("pcv13_rollout_to_stable", p, True, "PCV13", "calendar_pde_program_realism_pcv13_rollout_to_stable.csv"),
        ("pcv10sii_rollout_to_stable", p, True, "PCV10-SII", "calendar_pde_program_realism_pcv10sii_rollout_to_stable.csv"),
    ]


def run_all_week6_scenarios(
    outdir: Optional[str] = None,
    record_months: int = DEFAULT_RECORD_MONTHS,
    retain_preexisting_sequelae: bool = False,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    outdir_path = _resolve_outdir(outdir)
    csv_dir = outdir_path / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)

    monthly_outputs: Dict[str, pd.DataFrame] = {}
    yearly_outputs: Dict[str, pd.DataFrame] = {}

    for scenario_name, p, intervention, vaccine_type, filename in _scenario_specs(record_months):
        monthly = run_pde_calendar_week6_monthly(
            p,
            intervention=intervention,
            vaccine_type=vaccine_type,
            scenario_name=scenario_name,
            retain_preexisting_sequelae=retain_preexisting_sequelae,
        )
        yearly = base.aggregate_yearly_like_ode_new(monthly)
        monthly_outputs[scenario_name] = monthly
        yearly_outputs[scenario_name] = yearly
        monthly.to_csv(csv_dir / filename.replace(".csv", "_monthly.csv"), index=False)
        yearly.to_csv(csv_dir / filename, index=False)

    period_summary = base.summarize_transition_and_stable(yearly_outputs)
    stable_annual = base.summarize_stable_annual(yearly_outputs)
    diagnostics = base.stability_diagnostics(yearly_outputs)
    period_summary.to_csv(csv_dir / "program_realism_transition_and_stable_summary.csv", index=False)
    stable_annual.to_csv(csv_dir / "program_realism_stable_annual_summary.csv", index=False)
    diagnostics.to_csv(csv_dir / "program_realism_stability_diagnostics.csv", index=False)
    write_assumption_table(csv_dir / "program_realism_assumptions.csv")
    return monthly_outputs, yearly_outputs


def run_week5_reference_monthly(
    outdir: Optional[str] = None,
    record_months: int = DEFAULT_RECORD_MONTHS,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    outdir_path = _resolve_outdir(outdir) / "on_schedule_reference_csv"
    outdir_path.mkdir(parents=True, exist_ok=True)

    monthly_outputs: Dict[str, pd.DataFrame] = {}
    yearly_outputs: Dict[str, pd.DataFrame] = {}

    for scenario_name, p, intervention, vaccine_type, filename in w5._scenario_specs(record_months=record_months):
        monthly = w5.run_pde_calendar_refined_monthly(
            p,
            intervention=intervention,
            vaccine_type=vaccine_type,
            scenario_name=scenario_name,
            retain_preexisting_sequelae=False,
        )
        yearly = base.aggregate_yearly_like_ode_new(monthly)
        monthly_outputs[scenario_name] = monthly
        yearly_outputs[scenario_name] = yearly
        monthly.to_csv(outdir_path / filename.replace(".csv", "_monthly.csv"), index=False)
        yearly.to_csv(outdir_path / filename, index=False)

    return monthly_outputs, yearly_outputs


def write_assumption_table(path: Path) -> None:
    p = params_vietnam_micro
    rows = [
        ("birth_month_shares", ";".join(f"{x:.3f}" for x in p["birth_month_shares"]), "Assumed mild seasonality; replace with verified Vietnam monthly births."),
        ("routine_dose1_transition_probs", str(p["routine_dose1_transition_probs"]), "Delay around the 2-month target."),
        ("routine_dose2_transition_probs", str(p["routine_dose2_transition_probs"]), "Delay around the 4-month target."),
        ("routine_booster_transition_probs", str(p["routine_booster_transition_probs"]), "Delay around the 12-month booster target."),
        ("catchup_enabled", str(p["catchup_enabled"]), "Off in central run because introduction catch-up policy was not specified."),
        ("dose1_relative_ve", p["dose1_relative_ve"], "Partial direct protection after one dose."),
        ("dose2_relative_ve", p["dose2_relative_ve"], "Primary-series direct protection equals the on-schedule endpoint VE."),
        ("booster_relative_ve", p["booster_relative_ve"], "Booster restores protection, capped by max_direct_ve."),
        ("primary_waning_by_age", str(p["primary_waning_by_age"]), "Structural waning assumption after the primary series."),
        ("booster_waning_by_age", str(p["booster_waning_by_age"]), "Structural waning assumption after booster."),
        ("carriage_reduction_by_group", str(p["carriage_reduction_by_group"]), "Carriage RCT-inspired build-up of herd protection."),
        ("mature_effective_dose_share", p["mature_effective_dose_share"], "Used to translate mature IRR into a mature herd target."),
    ]
    pd.DataFrame(rows, columns=["Parameter", "Value", "Interpretation"]).to_csv(path, index=False)


def _add_total_cases(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["Monthly_Total_Cases"] = out[["New_Meningitis", "New_Pneumo_Sev", "New_Pneumo_NonSev", "New_Sepsis", "New_AOM"]].sum(axis=1)
    return out


def compare_week6_to_week5(
    week6_monthly: Dict[str, pd.DataFrame],
    week6_yearly: Dict[str, pd.DataFrame],
    week5_monthly: Dict[str, pd.DataFrame],
    week5_yearly: Dict[str, pd.DataFrame],
    outdir: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    outdir_path = _resolve_outdir(outdir) / "csv"
    outdir_path.mkdir(parents=True, exist_ok=True)

    flow_cols = ["New_Meningitis", "New_Pneumo_Sev", "New_Pneumo_NonSev", "New_Sepsis", "New_AOM", "New_Deaths_Pneumo"]
    monthly_rows = []
    yearly_rows = []

    for scenario, df6 in week6_monthly.items():
        df5 = week5_monthly[scenario]
        m6 = _add_total_cases(df6)
        m5 = _add_total_cases(df5)
        for month in sorted(set(m6["Month"]).intersection(set(m5["Month"]))):
            r6 = m6.loc[m6["Month"] == month].iloc[0]
            r5 = m5.loc[m5["Month"] == month].iloc[0]
            row = {"Scenario": scenario, "Month": int(month), "Year": int(r6["Year"])}
            for col in flow_cols + ["Monthly_Total_Cases", "Stock_Alive"]:
                row[f"Week5_{col}"] = float(r5[col])
                row[f"Week6_{col}"] = float(r6[col])
                row[f"Delta_{col}"] = float(r6[col] - r5[col])
            monthly_rows.append(row)

        y6 = week6_yearly[scenario]
        y5 = week5_yearly[scenario]
        for year in sorted(set(y6["Year"]).intersection(set(y5["Year"]))):
            r6 = y6.loc[y6["Year"] == year].iloc[0]
            r5 = y5.loc[y5["Year"] == year].iloc[0]
            row = {"Scenario": scenario, "Year": int(year)}
            for col in flow_cols + ["Stock_Alive"]:
                row[f"Week5_{col}"] = float(r5[col])
                row[f"Week6_{col}"] = float(r6[col])
                row[f"Delta_{col}"] = float(r6[col] - r5[col])
            row["Week5_Total_Cases"] = sum(float(r5[col]) for col in flow_cols[:-1])
            row["Week6_Total_Cases"] = sum(float(r6[col]) for col in flow_cols[:-1])
            row["Delta_Total_Cases"] = row["Week6_Total_Cases"] - row["Week5_Total_Cases"]
            yearly_rows.append(row)

    monthly_cmp = pd.DataFrame(monthly_rows)
    yearly_cmp = pd.DataFrame(yearly_rows)
    monthly_cmp.to_csv(outdir_path / "program_realism_vs_on_schedule_monthly_comparison.csv", index=False)
    yearly_cmp.to_csv(outdir_path / "program_realism_vs_on_schedule_yearly_comparison.csv", index=False)
    return monthly_cmp, yearly_cmp


def plot_week6_effects(
    week6_monthly: Dict[str, pd.DataFrame],
    week6_yearly: Dict[str, pd.DataFrame],
    week5_monthly: Dict[str, pd.DataFrame],
    week5_yearly: Dict[str, pd.DataFrame],
    outdir: Optional[str] = None,
) -> List[Path]:
    import matplotlib.pyplot as plt

    outdir_path = _resolve_outdir(outdir) / "figures"
    outdir_path.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    scenario_labels = {
        "baseline": "Baseline",
        "pcv10_gsk_rollout_to_stable": "PCV10-GSK",
        "pcv13_rollout_to_stable": "PCV13",
        "pcv10sii_rollout_to_stable": "PCV10-SII",
    }
    compare_scenarios = ["pcv13_rollout_to_stable", "pcv10sii_rollout_to_stable"]

    # Full ten-year monthly comparison.
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for scenario in compare_scenarios:
        df5 = _add_total_cases(week5_monthly[scenario])
        df6 = _add_total_cases(week6_monthly[scenario])
        ax.plot(df5["Month"] + 1, df5["Monthly_Total_Cases"], linestyle="--", linewidth=1.8, label=f"Previous assumptions {scenario_labels[scenario]}")
        ax.plot(df6["Month"] + 1, df6["Monthly_Total_Cases"], linewidth=2.0, label=f"Timing-refined assumptions {scenario_labels[scenario]}")
    ax.set_title("Monthly total cases, full 10-year horizon")
    ax.set_xlabel("Program month")
    ax.set_ylabel("Monthly cases")
    ax.legend(ncol=2, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = outdir_path / "program_realism_vs_on_schedule_monthly_cases_full_10y.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)

    # Early rollout zoom.
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for scenario in compare_scenarios:
        df5 = _add_total_cases(week5_monthly[scenario]).query("Month < 36")
        df6 = _add_total_cases(week6_monthly[scenario]).query("Month < 36")
        ax.plot(df5["Month"] + 1, df5["Monthly_Total_Cases"], linestyle="--", linewidth=1.8, label=f"Previous assumptions {scenario_labels[scenario]}")
        ax.plot(df6["Month"] + 1, df6["Monthly_Total_Cases"], linewidth=2.0, label=f"Timing-refined assumptions {scenario_labels[scenario]}")
    ax.set_title("Monthly total cases, first 36 months")
    ax.set_xlabel("Program month")
    ax.set_ylabel("Monthly cases")
    ax.legend(ncol=2, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = outdir_path / "program_realism_vs_on_schedule_monthly_cases_first_36m.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)

    # Cases averted compared with each model's own baseline, early zoom.
    fig, ax = plt.subplots(figsize=(11, 5.5))
    base5 = _add_total_cases(week5_monthly["baseline"])
    base6 = _add_total_cases(week6_monthly["baseline"])
    for scenario in ["pcv10_gsk_rollout_to_stable", "pcv13_rollout_to_stable", "pcv10sii_rollout_to_stable"]:
        df5 = _add_total_cases(week5_monthly[scenario])
        df6 = _add_total_cases(week6_monthly[scenario])
        av5 = (base5["Monthly_Total_Cases"] - df5["Monthly_Total_Cases"]).iloc[:36]
        av6 = (base6["Monthly_Total_Cases"] - df6["Monthly_Total_Cases"]).iloc[:36]
        x = np.arange(1, len(av5) + 1)
        ax.plot(x, av5, linestyle="--", linewidth=1.5, label=f"Previous assumptions {scenario_labels[scenario]}")
        ax.plot(x, av6, linewidth=1.9, label=f"Timing-refined assumptions {scenario_labels[scenario]}")
    ax.set_title("Monthly cases averted, first 36 months")
    ax.set_xlabel("Program month")
    ax.set_ylabel("Cases averted vs same-model baseline")
    ax.legend(ncol=2, fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = outdir_path / "program_realism_vs_on_schedule_cases_averted_first_36m.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)

    # Annual total cases over 10 years.
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for scenario in ["pcv10_gsk_rollout_to_stable", "pcv13_rollout_to_stable", "pcv10sii_rollout_to_stable"]:
        y5 = week5_yearly[scenario].copy()
        y6 = week6_yearly[scenario].copy()
        y5["Total_Cases"] = y5[["New_Meningitis", "New_Pneumo_Sev", "New_Pneumo_NonSev", "New_Sepsis", "New_AOM"]].sum(axis=1)
        y6["Total_Cases"] = y6[["New_Meningitis", "New_Pneumo_Sev", "New_Pneumo_NonSev", "New_Sepsis", "New_AOM"]].sum(axis=1)
        ax.plot(y5["Year"], y5["Total_Cases"], linestyle="--", marker="o", markersize=3, label=f"Previous assumptions {scenario_labels[scenario]}")
        ax.plot(y6["Year"], y6["Total_Cases"], marker="o", markersize=3, label=f"Timing-refined assumptions {scenario_labels[scenario]}")
    ax.set_title("Annual total cases, 10-year horizon")
    ax.set_xlabel("Program year")
    ax.set_ylabel("Annual cases")
    ax.legend(ncol=2, fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = outdir_path / "program_realism_vs_on_schedule_annual_cases_10y.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)

    # Dose-state shares during rollout.
    fig, ax = plt.subplots(figsize=(11, 5.5))
    dose_df = week6_monthly["pcv13_rollout_to_stable"].query("Month < 48")
    x = dose_df["Month"] + 1
    ax.plot(x, dose_df["AtLeastDose1_Share"] * 100.0, label="At least dose 1")
    ax.plot(x, dose_df["AtLeastDose2_Share"] * 100.0, label="At least dose 2")
    ax.plot(x, dose_df["Booster_Share"] * 100.0, label="Booster")
    ax.plot(x, dose_df["Relative_Carriage_Index"] * 100.0, label="Relative carriage index")
    ax.set_title("Dose-state and carriage rollout, PCV13 first 48 months")
    ax.set_xlabel("Program month")
    ax.set_ylabel("Percent of reference / population")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = outdir_path / "program_realism_dose_and_carriage_rollout_first_48m.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    written.append(path)

    return written


def main() -> None:
    outdir = DEFAULT_OUTPUT_DIR
    monthly6, yearly6 = run_all_week6_scenarios(outdir=outdir, record_months=DEFAULT_RECORD_MONTHS)
    monthly5, yearly5 = run_week5_reference_monthly(outdir=outdir, record_months=DEFAULT_RECORD_MONTHS)
    monthly_cmp, yearly_cmp = compare_week6_to_week5(monthly6, yearly6, monthly5, yearly5, outdir=outdir)
    figures = plot_week6_effects(monthly6, yearly6, monthly5, yearly5, outdir=outdir)

    period_summary = base.summarize_transition_and_stable(yearly6)
    stable_annual = base.summarize_stable_annual(yearly6)
    diagnostics = base.stability_diagnostics(yearly6)

    print("\nProgramme-realism transition and stable-period summary")
    print(period_summary.to_string(index=False))
    print("\nProgramme-realism stable annual summary")
    print(stable_annual.to_string(index=False))
    print("\nProgramme-realism stability diagnostics, year 10 vs year 9")
    print(diagnostics[["Scenario", "MaxAbsPctChange_FlowsAndAlive"]].to_string(index=False))
    print("\nProgramme-realism minus on-schedule yearly total cases, by year")
    print(yearly_cmp[["Scenario", "Year", "Delta_Total_Cases", "Delta_New_Deaths_Pneumo"]].to_string(index=False))
    print("\nFigures written:")
    for fig in figures:
        print(f"  {fig}")


if __name__ == "__main__":
    main()
