"""
Explicit age-time PDE calendar model, original vaccine logic.

This file keeps the same estimand and vaccine assumptions as
Week 4_calendar_original.py, but writes the monthly update as a solver-based
method-of-lines / finite-volume PDE calculation:

    partial_t A_g(a,t) + partial_a A_g(a,t) = -mu_g(a,t) A_g(a,t).

Within each recorded month, solve_ivp integrates the reaction and event-flow
part for all age bins. Month-end age transport then moves survivors from age
bin a to a+1, while age 59 exits the under-five population. New births enter as
the age-0 boundary condition before each recorded monthly solve.

Warm-up is deliberately baseline-only. It creates the existing under-five age
structure at program start without assuming the vaccination program already
existed during the previous five years.
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
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "Outputs" / "Week 5_calendar_pde_original"
DEFAULT_RECORD_MONTHS = 120
TRANSITION_YEARS = 5
STABLE_START_YEAR = 6
STABLE_END_YEAR = 10


def _load_peer(filename: str, module_name: str):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _resolve_outdir(outdir: Optional[str]) -> str:
    if outdir is None:
        return str(DEFAULT_OUTPUT_DIR)
    path = Path(outdir)
    if not path.is_absolute():
        path = SCRIPT_DIR / path
    return str(path)


week4 = _load_peer("Week 4_calendar_original.py", "week4_calendar_original_reference")

params_vietnam_micro = copy.deepcopy(week4.params_vietnam_micro)
OUTPUT_COLUMNS = list(week4.OUTPUT_COLUMNS)
SEQ_KEYS = list(week4.SEQ_KEYS)

normalize_birth_shares = week4.normalize_birth_shares
get_age_weight = week4.get_age_weight
monthly_prob_from_incidence_per100k = week4.monthly_prob_from_incidence_per100k
infer_ve_from_irr_series = week4.infer_ve_from_irr_series
_clip01 = week4._clip01
_prob_to_hazard = week4._prob_to_hazard
get_background_prob_u5_only = week4.get_background_prob_u5_only
neonatal_survival_multiplier = week4.neonatal_survival_multiplier
empty_sequela_stock = week4.empty_sequela_stock
copy_sequela_stock = week4.copy_sequela_stock
age_shift_for_next_month = week4.age_shift_for_next_month
aggregate_yearly_like_ode_new = week4.aggregate_yearly_like_ode_new
summarize_five_year_totals = week4.summarize_five_year_totals


FLOW_NAMES = [
    "New_Meningitis",
    "New_Pneumo_Sev",
    "New_Pneumo_NonSev",
    "New_Sepsis",
    "New_AOM",
    "New_Deaths_Pneumo",
]


def _with_record_months(p: Dict, record_months: int) -> Dict:
    out = copy.deepcopy(p)
    out["record_months"] = int(record_months)
    return out


def _scenario_specs(record_months: int = DEFAULT_RECORD_MONTHS) -> List[Tuple[str, Dict, bool, Optional[str], str]]:
    p = _with_record_months(params_vietnam_micro, record_months)
    specs: List[Tuple[str, Dict, bool, Optional[str], str]] = [
        ("baseline", p, False, None, "calendar_pde_original_baseline.csv"),
        ("pcv10_gsk_rollout_to_stable", p, True, "PCV10-GSK", "calendar_pde_original_pcv10_gsk_rollout_to_stable.csv"),
        ("pcv13_rollout_to_stable", p, True, "PCV13", "calendar_pde_original_pcv13_rollout_to_stable.csv"),
        ("pcv10sii_rollout_to_stable", p, True, "PCV10-SII", "calendar_pde_original_pcv10sii_rollout_to_stable.csv"),
    ]
    return specs


def _disease_probabilities_for_age(p: Dict, age_month: int) -> Tuple[float, float, float, float]:
    w_men = get_age_weight(age_month, "Meningitis")
    w_pneu = get_age_weight(age_month, "Pneumonia")
    w_aom = get_age_weight(age_month, "AOM")
    base_prob_men = monthly_prob_from_incidence_per100k(p["inc_meningitis"], w_men)
    base_prob_pneu = monthly_prob_from_incidence_per100k(p["inc_pneumonia"], w_pneu)
    base_prob_sep = monthly_prob_from_incidence_per100k(p["inc_sepsis"], w_pneu)
    base_prob_aom = monthly_prob_from_incidence_per100k(p["inc_aom_all"], w_aom) * float(p["prop_aom_spn"])
    return base_prob_men, base_prob_pneu, base_prob_sep, base_prob_aom


def build_original_month_coefficients(
    p: Dict,
    calendar_month: int,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    n_age = int(p.get("age_bins", 60))
    year_idx = int(calendar_month) // 12

    if intervention and vaccine_type is not None:
        cov = float(p.get("vaccine_coverage", 0.0))
        _, ve_year = infer_ve_from_irr_series(p["vaccine_irrs"], vaccine_type, year_idx, cov)
        mult_v = 1.0 - ve_year
    else:
        mult_v = 1.0

    prop_sev = float(p["prop_pneumo_severe"])
    cfr_men = float(p["cfr_meningitis"])
    cfr_pneu = float(p["cfr_pneumo_severe"])
    cfr_sep = float(p["cfr_sepsis"])
    lam_bg = _prob_to_hazard(get_background_prob_u5_only(p))

    coeffs = {name: np.zeros(n_age, dtype=float) for name in [
        "lam_men_v", "lam_pneu_v", "lam_sep_v", "lam_aom_v", "mu_dis_v", "mu_all_v",
        "lam_men_u", "lam_pneu_u", "lam_sep_u", "lam_aom_u", "mu_dis_u", "mu_all_u",
    ]}

    for age_month in range(n_age):
        base_prob_men, base_prob_pneu, base_prob_sep, base_prob_aom = _disease_probabilities_for_age(p, age_month)

        p_men_v = _clip01(base_prob_men * mult_v)
        p_pneu_v = _clip01(base_prob_pneu * mult_v)
        p_sep_v = _clip01(base_prob_sep * mult_v)
        p_aom_v = _clip01(base_prob_aom * mult_v)

        p_men_u = _clip01(base_prob_men)
        p_pneu_u = _clip01(base_prob_pneu)
        p_sep_u = _clip01(base_prob_sep)
        p_aom_u = _clip01(base_prob_aom)

        coeffs["lam_men_v"][age_month] = _prob_to_hazard(p_men_v)
        coeffs["lam_pneu_v"][age_month] = _prob_to_hazard(p_pneu_v)
        coeffs["lam_sep_v"][age_month] = _prob_to_hazard(p_sep_v)
        coeffs["lam_aom_v"][age_month] = _prob_to_hazard(p_aom_v)
        coeffs["lam_men_u"][age_month] = _prob_to_hazard(p_men_u)
        coeffs["lam_pneu_u"][age_month] = _prob_to_hazard(p_pneu_u)
        coeffs["lam_sep_u"][age_month] = _prob_to_hazard(p_sep_u)
        coeffs["lam_aom_u"][age_month] = _prob_to_hazard(p_aom_u)

        d1_v = _clip01(p_men_v * cfr_men)
        d2_v = _clip01(p_pneu_v * prop_sev * cfr_pneu)
        d3_v = _clip01(p_sep_v * cfr_sep)
        d1_u = _clip01(p_men_u * cfr_men)
        d2_u = _clip01(p_pneu_u * prop_sev * cfr_pneu)
        d3_u = _clip01(p_sep_u * cfr_sep)

        coeffs["mu_dis_v"][age_month] = _prob_to_hazard(d1_v) + _prob_to_hazard(d2_v) + _prob_to_hazard(d3_v)
        coeffs["mu_dis_u"][age_month] = _prob_to_hazard(d1_u) + _prob_to_hazard(d2_u) + _prob_to_hazard(d3_u)
        coeffs["mu_all_v"][age_month] = lam_bg + coeffs["mu_dis_v"][age_month]
        coeffs["mu_all_u"][age_month] = lam_bg + coeffs["mu_dis_u"][age_month]

    coeffs["prop_sev"] = np.array([prop_sev], dtype=float)
    return coeffs


def pde_reaction_rhs(_tau: float, y: np.ndarray, coeffs: Dict[str, np.ndarray]) -> np.ndarray:
    n_age = coeffs["mu_all_v"].size
    A_v = np.maximum(y[:n_age], 0.0)
    A_u = np.maximum(y[n_age:2 * n_age], 0.0)

    dy = np.zeros_like(y, dtype=float)
    dy[:n_age] = -coeffs["mu_all_v"] * A_v
    dy[n_age:2 * n_age] = -coeffs["mu_all_u"] * A_u

    prop_sev = float(coeffs["prop_sev"][0])
    inc_men = coeffs["lam_men_v"] * A_v + coeffs["lam_men_u"] * A_u
    inc_pneu = coeffs["lam_pneu_v"] * A_v + coeffs["lam_pneu_u"] * A_u
    inc_sep = coeffs["lam_sep_v"] * A_v + coeffs["lam_sep_u"] * A_u
    inc_aom = coeffs["lam_aom_v"] * A_v + coeffs["lam_aom_u"] * A_u
    inc_death = coeffs["mu_dis_v"] * A_v + coeffs["mu_dis_u"] * A_u

    offset = 2 * n_age
    dy[offset + 0] = float(inc_men.sum())
    dy[offset + 1] = float((prop_sev * inc_pneu).sum())
    dy[offset + 2] = float(((1.0 - prop_sev) * inc_pneu).sum())
    dy[offset + 3] = float(inc_sep.sum())
    dy[offset + 4] = float(inc_aom.sum())
    dy[offset + 5] = float(inc_death.sum())
    return dy


def solve_one_month_pde_original(
    A_v: np.ndarray,
    A_u: np.ndarray,
    S: Dict[str, float],
    calendar_month: int,
    p: Dict,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
    rtol: float = 1e-10,
    atol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], Dict[str, float], float]:
    coeffs = build_original_month_coefficients(p, calendar_month, intervention, vaccine_type)
    n_age = len(A_v)
    y0 = np.concatenate([A_v.astype(float), A_u.astype(float), np.zeros(len(FLOW_NAMES), dtype=float)])

    sol = solve_ivp(
        fun=lambda tau, yy: pde_reaction_rhs(tau, yy, coeffs),
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
    A_v_after = y1[:n_age]
    A_u_after = y1[n_age:2 * n_age]
    flow_values = y1[2 * n_age:2 * n_age + len(FLOW_NAMES)]
    flows = {name: float(value) for name, value in zip(FLOW_NAMES, flow_values)}

    alive_before = float(A_v.sum() + A_u.sum())
    alive_after = float(A_v_after.sum() + A_u_after.sum())
    non_neonatal_all_cause_deaths = max(0.0, alive_before - alive_after)

    surv_seq = max(0.0, min(1.0, alive_after / alive_before)) if alive_before > 0.0 else 1.0
    men_survivor_flow = (1.0 - float(p["cfr_meningitis"])) * flows["New_Meningitis"]

    src_cog = float(p["prob_cog"]) * men_survivor_flow
    src_hear = float(p["prob_hear"]) * men_survivor_flow
    src_epi = float(p["prob_epi"]) * men_survivor_flow
    src_hemi = float(p["prob_hemi"]) * men_survivor_flow
    src_impl = float(p.get("prob_cochlear_implant", 0.0)) * src_hear
    src_sped = float(p.get("prob_special_education", 0.0)) * src_hear

    S_after = copy_sequela_stock(S)
    S_after["Stock_Seq_Cognitive"] = S_after["Stock_Seq_Cognitive"] * surv_seq + src_cog
    S_after["Stock_Seq_Hearing"] = S_after["Stock_Seq_Hearing"] * surv_seq + src_hear
    S_after["Stock_Seq_Epilepsy"] = S_after["Stock_Seq_Epilepsy"] * surv_seq + src_epi
    S_after["Stock_Seq_Hemiplegia"] = S_after["Stock_Seq_Hemiplegia"] * surv_seq + src_hemi
    S_after["Stock_Seq_Hearing_Implant"] = S_after["Stock_Seq_Hearing_Implant"] * surv_seq + src_impl
    S_after["Stock_Seq_Hearing_SpEd"] = S_after["Stock_Seq_Hearing_SpEd"] * surv_seq + src_sped

    return A_v_after, A_u_after, S_after, flows, non_neonatal_all_cause_deaths


def add_record_births_original(
    A_v: np.ndarray,
    A_u: np.ndarray,
    gross_births: float,
    p: Dict,
    intervention: bool,
    vaccine_type: Optional[str],
) -> Tuple[float, float]:
    alive_births = float(gross_births) * neonatal_survival_multiplier(p)
    neonatal_deaths = float(gross_births) - alive_births
    if intervention and vaccine_type is not None:
        cov = float(p.get("vaccine_coverage", 0.0))
        A_v[0] += alive_births * cov
        A_u[0] += alive_births * (1.0 - cov)
    else:
        A_u[0] += alive_births
    return alive_births, neonatal_deaths


def build_pde_initial_population(
    p: Dict,
    shares: np.ndarray,
    scenario_name: str,
    retain_preexisting_sequelae: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], Dict[str, float]]:
    n_age = int(p.get("age_bins", 60))
    warmup_months = int(p.get("warmup_months", 60))
    births_annual = float(p["birth_cohort"])
    surv0 = neonatal_survival_multiplier(p)

    A_v = np.zeros(n_age, dtype=float)
    A_u = np.zeros(n_age, dtype=float)
    S_warm = empty_sequela_stock()

    for warm_month in range(warmup_months):
        gross_births = births_annual * float(shares[warm_month % 12])
        A_u[0] += gross_births * surv0
        A_v_after, A_u_after, S_warm, _, _ = solve_one_month_pde_original(
            A_v, A_u, S_warm, warm_month, p, intervention=False, vaccine_type=None
        )
        A_v, A_u, _ = age_shift_for_next_month(A_v_after, A_u_after)

    S_initial = copy_sequela_stock(S_warm) if retain_preexisting_sequelae else empty_sequela_stock()
    warm_meta = {
        "scenario_name": scenario_name,
        "warmup_months": warmup_months,
        "warmup_mode": "baseline_only_no_vaccination",
        "warmup_alive_after_age_shift": float(A_v.sum() + A_u.sum()),
        "warmup_vaccinated_or_program_alive": float(A_v.sum()),
        "warmup_unvaccinated_alive": float(A_u.sum()),
        "warmup_preexisting_sequelae_total": float(sum(S_warm.values())),
        "retain_preexisting_sequelae": bool(retain_preexisting_sequelae),
    }
    return A_v, A_u, S_initial, warm_meta


def run_pde_calendar_original_monthly(
    p: Dict,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
    scenario_name: str = "scenario",
    retain_preexisting_sequelae: bool = False,
) -> pd.DataFrame:
    p = copy.deepcopy(p)
    record_months = int(p.get("record_months", 60))
    shares = normalize_birth_shares(p.get("birth_month_shares", [1.0 / 12.0] * 12))
    births_annual = float(p["birth_cohort"])

    A_v, A_u, S, _ = build_pde_initial_population(
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
        _, neonatal_deaths = add_record_births_original(A_v, A_u, gross_births, p, intervention, vaccine_type)
        cum_all_cause_deaths += neonatal_deaths

        A_v_after, A_u_after, S, flows, non_neonatal_deaths = solve_one_month_pde_original(
            A_v, A_u, S, calendar_month, p, intervention=intervention, vaccine_type=vaccine_type
        )
        cum_all_cause_deaths += non_neonatal_deaths

        stock_alive = float(A_v_after.sum() + A_u_after.sum())
        cum_men += flows["New_Meningitis"]
        cum_psev += flows["New_Pneumo_Sev"]
        cum_pnon += flows["New_Pneumo_NonSev"]
        cum_sep += flows["New_Sepsis"]
        cum_aom += flows["New_AOM"]
        cum_ddeath += flows["New_Deaths_Pneumo"]

        history.append({
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
        })

        A_v, A_u, _ = age_shift_for_next_month(A_v_after, A_u_after)

    return pd.DataFrame(history)


def summarize_year_window(
    outputs: Dict[str, pd.DataFrame],
    start_year: int,
    end_year: int,
    label: str,
) -> pd.DataFrame:
    rows = []
    for name, df in outputs.items():
        window = df[(df["Year"] >= int(start_year)) & (df["Year"] <= int(end_year))]
        if window.empty:
            raise ValueError(f"No rows found for {name} in years {start_year}-{end_year}")
        years = int(window["Year"].nunique())
        row = {
            "Scenario": name,
            "Window": label,
            "Start_Year": int(start_year),
            "End_Year": int(end_year),
            "Years": years,
            "Meningitis": float(window["New_Meningitis"].sum()),
            "Pneumo_Sev": float(window["New_Pneumo_Sev"].sum()),
            "Pneumo_NonSev": float(window["New_Pneumo_NonSev"].sum()),
            "Sepsis": float(window["New_Sepsis"].sum()),
            "AOM": float(window["New_AOM"].sum()),
            "Disease_Deaths": float(window["New_Deaths_Pneumo"].sum()),
            "Alive_End": float(window.iloc[-1]["Stock_Alive"]),
            "All_Cause_Deaths_End": float(window.iloc[-1]["Stock_Dead"]),
        }
        for col in ["Meningitis", "Pneumo_Sev", "Pneumo_NonSev", "Sepsis", "AOM", "Disease_Deaths"]:
            row[f"Annual_{col}"] = row[col] / float(years)
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_transition_and_stable(
    outputs: Dict[str, pd.DataFrame],
    transition_years: int = TRANSITION_YEARS,
    stable_start_year: int = STABLE_START_YEAR,
    stable_end_year: int = STABLE_END_YEAR,
) -> pd.DataFrame:
    transition = summarize_year_window(outputs, 1, transition_years, f"Y1-Y{transition_years}_rollout_transition")
    stable = summarize_year_window(outputs, stable_start_year, stable_end_year, f"Y{stable_start_year}-Y{stable_end_year}_mature_stable")
    return pd.concat([transition, stable], ignore_index=True)


def summarize_stable_annual(
    outputs: Dict[str, pd.DataFrame],
    stable_year: int = STABLE_END_YEAR,
) -> pd.DataFrame:
    rows = []
    for name, df in outputs.items():
        target = df[df["Year"] == int(stable_year)]
        if target.empty:
            raise ValueError(f"No rows found for {name} in stable year {stable_year}")
        row = target.iloc[-1]
        rows.append({
            "Scenario": name,
            "Stable_Year": int(stable_year),
            "Meningitis_Annual": float(row["New_Meningitis"]),
            "Pneumo_Sev_Annual": float(row["New_Pneumo_Sev"]),
            "Pneumo_NonSev_Annual": float(row["New_Pneumo_NonSev"]),
            "Sepsis_Annual": float(row["New_Sepsis"]),
            "AOM_Annual": float(row["New_AOM"]),
            "Disease_Deaths_Annual": float(row["New_Deaths_Pneumo"]),
            "Alive_End_Stable_Year": float(row["Stock_Alive"]),
            "All_Cause_Deaths_End_Stable_Year": float(row["Stock_Dead"]),
        })
    return pd.DataFrame(rows)


def stability_diagnostics(
    outputs: Dict[str, pd.DataFrame],
    compare_year: int = STABLE_END_YEAR,
    previous_year: int = STABLE_END_YEAR - 1,
) -> pd.DataFrame:
    flow_cols = ["New_Meningitis", "New_Pneumo_Sev", "New_Pneumo_NonSev", "New_Sepsis", "New_AOM", "New_Deaths_Pneumo"]
    rows = []
    for name, df in outputs.items():
        current = df[df["Year"] == int(compare_year)]
        previous = df[df["Year"] == int(previous_year)]
        if current.empty or previous.empty:
            raise ValueError(f"Missing diagnostic years for {name}")
        current_row = current.iloc[-1]
        previous_row = previous.iloc[-1]
        row = {"Scenario": name, "Previous_Year": int(previous_year), "Compare_Year": int(compare_year)}
        pct_changes = []
        for col in flow_cols + ["Stock_Alive"]:
            prev = float(previous_row[col])
            curr = float(current_row[col])
            delta = curr - prev
            pct = delta / prev * 100.0 if prev != 0.0 else 0.0
            row[f"{col}_Delta"] = delta
            row[f"{col}_Pct_Change"] = pct
            pct_changes.append(abs(pct))
        row["MaxAbsPctChange_FlowsAndAlive"] = max(pct_changes) if pct_changes else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def run_all_scenarios_and_save_csv(
    outdir: Optional[str] = None,
    retain_preexisting_sequelae: bool = False,
    record_months: int = DEFAULT_RECORD_MONTHS,
) -> Dict[str, pd.DataFrame]:
    outdir = _resolve_outdir(outdir)
    os.makedirs(outdir, exist_ok=True)
    outputs: Dict[str, pd.DataFrame] = {}

    for scenario_name, p, intervention, vaccine_type, filename in _scenario_specs(record_months):
        monthly = run_pde_calendar_original_monthly(
            p,
            intervention=intervention,
            vaccine_type=vaccine_type,
            scenario_name=scenario_name,
            retain_preexisting_sequelae=retain_preexisting_sequelae,
        )
        yearly = aggregate_yearly_like_ode_new(monthly)
        outputs[scenario_name] = yearly
        monthly.to_csv(os.path.join(outdir, filename.replace(".csv", "_monthly.csv")), index=False)
        yearly.to_csv(os.path.join(outdir, filename), index=False)

    period_summary = summarize_transition_and_stable(outputs)
    stable_annual = summarize_stable_annual(outputs)
    diagnostics = stability_diagnostics(outputs)
    period_summary.to_csv(os.path.join(outdir, "transition_and_stable_summary.csv"), index=False)
    stable_annual.to_csv(os.path.join(outdir, "stable_annual_summary.csv"), index=False)
    diagnostics.to_csv(os.path.join(outdir, "stability_diagnostics.csv"), index=False)
    return outputs


def compare_to_week4_calendar(
    week5_outputs: Optional[Dict[str, pd.DataFrame]] = None,
    outdir: Optional[str] = None,
    retain_preexisting_sequelae: bool = False,
    record_months: int = DEFAULT_RECORD_MONTHS,
) -> pd.DataFrame:
    outdir = _resolve_outdir(outdir)
    if week5_outputs is None:
        week5_outputs = run_all_scenarios_and_save_csv(
            outdir=outdir,
            retain_preexisting_sequelae=retain_preexisting_sequelae,
            record_months=record_months,
        )

    rows = []
    value_cols = [c for c in OUTPUT_COLUMNS if c != "Year"]
    for scenario_name, p, intervention, vaccine_type, _filename in _scenario_specs(record_months):
        ref_monthly = week4.run_age_time_calendar_original_monthly(
            p,
            intervention=intervention,
            vaccine_type=vaccine_type,
            scenario_name=scenario_name,
            retain_preexisting_sequelae=retain_preexisting_sequelae,
            output_initial_state_dir=None,
        )
        ref_yearly = week4.aggregate_yearly_like_ode_new(ref_monthly)
        got_yearly = week5_outputs[scenario_name]
        diff = got_yearly[value_cols].astype(float).to_numpy() - ref_yearly[value_cols].astype(float).to_numpy()
        rows.append({
            "Scenario": scenario_name,
            "MaxAbsRoundedYearlyDiff": float(np.max(np.abs(diff))),
            "RoundedYearlyOutputsMatchWeek4": bool(np.all(diff == 0.0)),
        })

    comparison = pd.DataFrame(rows)
    os.makedirs(outdir, exist_ok=True)
    comparison.to_csv(os.path.join(outdir, "calendar_pde_original_vs_calendar_ode_comparison.csv"), index=False)
    return comparison


if __name__ == "__main__":
    outputs = run_all_scenarios_and_save_csv()
    comparison = compare_to_week4_calendar(outputs)
    period_summary = summarize_transition_and_stable(outputs)
    stable_annual = summarize_stable_annual(outputs)
    diagnostics = stability_diagnostics(outputs)
    print("\nTransition and stable-period summary, Week 5 PDE calendar original model")
    print(period_summary.to_string(index=False))
    print("\nStable annual summary, Week 5 PDE calendar original model")
    print(stable_annual.to_string(index=False))
    print("\nStability diagnostics, year 10 vs year 9")
    print(diagnostics[["Scenario", "MaxAbsPctChange_FlowsAndAlive"]].to_string(index=False))
    print("\nRounded yearly comparison against Week 4 calendar original")
    print(comparison.to_string(index=False))
