import copy
import math
import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import ode_new as ode
except ImportError:
    import ode as ode

try:
    import combine_new as cmb
except ImportError:
    cmb = None

try:
    import QALY_new as qmod
except ImportError:
    qmod = None


SCENARIO_MAP = {
    "baseline": (False, None),
    "pcv10_dynamic": (True, "PCV10-GSK"),
    "pcv13_dynamic": (True, "PCV13"),
    "pcv10_year5_immediate": (True, "PCV10-GSK"),
    "pcv13_year5_immediate": (True, "PCV13"),
    "pcv10sii_year5_immediate": (True, "PCV10-SII"),
}

if cmb is not None:
    SCENARIO_MAP = {name: (intervention, vaccine_type) for name, intervention, vaccine_type in cmb.SCENARIOS}


# ============================================================
# Exact monthly propagator for the current cohort ODE structure
# ============================================================

def clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def prob_to_hazard(prob: float) -> float:
    p = clip01(prob)
    p = min(p, 1.0 - 1e-12)
    return -math.log(1.0 - p)


def phi(mu: float, dt: float = 1.0) -> float:
    if abs(mu) < 1e-14:
        return dt
    return (1.0 - math.exp(-mu * dt)) / mu


def scenario_params(base_params: Dict[str, Any], scenario_name: str) -> Tuple[Dict[str, Any], bool, Optional[str]]:
    params = copy.deepcopy(base_params)
    if scenario_name not in SCENARIO_MAP:
        raise ValueError(f"Unknown scenario: {scenario_name}")
    intervention, vaccine_type = SCENARIO_MAP[scenario_name]
    if ("year5_immediate" in scenario_name) and (vaccine_type is not None):
        patch = ode._make_immediate_irrs_from_year5(vaccine_type)
        merged = dict(params.get("vaccine_irrs", {}))
        merged.update(patch)
        params["vaccine_irrs"] = merged
    return params, intervention, vaccine_type


def month_coefficients(
    p: Dict[str, Any],
    month: int,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
) -> Dict[str, float]:
    month = int(max(0, min(59, month)))
    prob_bg_month = ode.get_background_prob_u5_only(p)
    lam_bg = prob_to_hazard(prob_bg_month)

    w_men = ode.get_age_weight(month, "Meningitis")
    w_pneu = ode.get_age_weight(month, "Pneumonia")
    w_aom = ode.get_age_weight(month, "AOM")

    base_prob_men = ode.monthly_prob_from_incidence_per100k(p["inc_meningitis"], w_men)
    base_prob_pneu = ode.monthly_prob_from_incidence_per100k(p["inc_pneumonia"], w_pneu)
    base_prob_sep = ode.monthly_prob_from_incidence_per100k(p["inc_sepsis"], w_pneu)
    base_prob_aom = ode.monthly_prob_from_incidence_per100k(p["inc_aom_all"], w_aom) * float(p["prop_aom_spn"])

    if intervention and (vaccine_type is not None):
        year_idx = month // 12
        cov = float(p.get("vaccine_coverage", 0.0))
        _, ve_year = ode.infer_ve_from_irr_series(p["vaccine_irrs"], vaccine_type, year_idx, cov)
        mult_v = 1.0 - ve_year
    else:
        mult_v = 1.0

    p_men_v = clip01(base_prob_men * mult_v)
    p_pneu_v = clip01(base_prob_pneu * mult_v)
    p_sep_v = clip01(base_prob_sep * mult_v)
    p_aom_v = clip01(base_prob_aom * mult_v)

    p_men_u = clip01(base_prob_men)
    p_pneu_u = clip01(base_prob_pneu)
    p_sep_u = clip01(base_prob_sep)
    p_aom_u = clip01(base_prob_aom)

    lam_men_v = prob_to_hazard(p_men_v)
    lam_pneu_v = prob_to_hazard(p_pneu_v)
    lam_sep_v = prob_to_hazard(p_sep_v)
    lam_aom_v = prob_to_hazard(p_aom_v)

    lam_men_u = prob_to_hazard(p_men_u)
    lam_pneu_u = prob_to_hazard(p_pneu_u)
    lam_sep_u = prob_to_hazard(p_sep_u)
    lam_aom_u = prob_to_hazard(p_aom_u)

    prop_sev = float(p["prop_pneumo_severe"])
    cfr_men = float(p["cfr_meningitis"])
    cfr_pneu = float(p["cfr_pneumo_severe"])
    cfr_sep = float(p["cfr_sepsis"])

    d1_v = clip01(p_men_v * cfr_men)
    d2_v = clip01(p_pneu_v * prop_sev * cfr_pneu)
    d3_v = clip01(p_sep_v * cfr_sep)

    d1_u = clip01(p_men_u * cfr_men)
    d2_u = clip01(p_pneu_u * prop_sev * cfr_pneu)
    d3_u = clip01(p_sep_u * cfr_sep)

    mu_dis_v = prob_to_hazard(d1_v) + prob_to_hazard(d2_v) + prob_to_hazard(d3_v)
    mu_dis_u = prob_to_hazard(d1_u) + prob_to_hazard(d2_u) + prob_to_hazard(d3_u)

    mu_all_v = lam_bg + mu_dis_v
    mu_all_u = lam_bg + mu_dis_u

    return {
        "lam_bg": lam_bg,
        "lam_men_v": lam_men_v,
        "lam_pneu_v": lam_pneu_v,
        "lam_sep_v": lam_sep_v,
        "lam_aom_v": lam_aom_v,
        "lam_men_u": lam_men_u,
        "lam_pneu_u": lam_pneu_u,
        "lam_sep_u": lam_sep_u,
        "lam_aom_u": lam_aom_u,
        "mu_dis_v": mu_dis_v,
        "mu_dis_u": mu_dis_u,
        "mu_all_v": mu_all_v,
        "mu_all_u": mu_all_u,
        "prop_sev": prop_sev,
        "cfr_men": cfr_men,
    }


def alive_exact(a0: float, mu: float, dt: float = 1.0) -> float:
    return float(a0) * math.exp(-mu * dt)


def cumulative_increment(alpha_v: float, alpha_u: float, a_v0: float, a_u0: float, mu_v: float, mu_u: float, dt: float = 1.0) -> float:
    return alpha_v * a_v0 * phi(mu_v, dt) + alpha_u * a_u0 * phi(mu_u, dt)


def integrated_case_mix(lambda_v: float, lambda_u: float, a_v0: float, a_u0: float, mu_v: float, mu_u: float, dt: float = 1.0) -> float:
    n0 = a_v0 + a_u0
    if n0 <= 0.0:
        return 0.0
    if a_v0 <= 0.0:
        return lambda_u * dt
    if a_u0 <= 0.0:
        return lambda_v * dt

    delta_mu = mu_v - mu_u
    if abs(delta_mu) < 1e-12:
        frac_v = a_v0 / n0
        return (lambda_u + (lambda_v - lambda_u) * frac_v) * dt

    r = a_v0 / a_u0
    term = math.log1p(r) - math.log1p(r * math.exp(-delta_mu * dt))
    return lambda_u * dt + ((lambda_v - lambda_u) / delta_mu) * term


def sequela_exact(
    s0: float,
    q_eff: float,
    cfr_men: float,
    lambda_men_v: float,
    lambda_men_u: float,
    a_v0: float,
    a_u0: float,
    mu_v: float,
    mu_u: float,
    dt: float = 1.0,
) -> float:
    n0 = a_v0 + a_u0
    if n0 <= 0.0:
        return 0.0
    n_dt = a_v0 * math.exp(-mu_v * dt) + a_u0 * math.exp(-mu_u * dt)
    mix_integral = integrated_case_mix(lambda_men_v, lambda_men_u, a_v0, a_u0, mu_v, mu_u, dt)
    return n_dt * ((s0 / n0) + q_eff * (1.0 - cfr_men) * mix_integral)


def run_exact_ode_and_build_monthly(
    p: Dict[str, Any],
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
) -> pd.DataFrame:
    birth = float(p["birth_cohort"])
    y = np.zeros(ode.N_STATE, dtype=float)

    if intervention and (vaccine_type is not None):
        cov = float(p.get("vaccine_coverage", 0.0))
        y[ode.IDX_A_V] = birth * cov
        y[ode.IDX_A_U] = birth * (1.0 - cov)
    else:
        y[ode.IDX_A_V] = 0.0
        y[ode.IDX_A_U] = birth

    p_neonatal = float(p.get("neonatal_mortality_rate", 0.0)) / 1000.0
    y[ode.IDX_A_V] *= (1.0 - p_neonatal)
    y[ode.IDX_A_U] *= (1.0 - p_neonatal)

    history: List[Dict[str, float]] = []

    for month in range(60):
        coeff = month_coefficients(p, month, intervention=intervention, vaccine_type=vaccine_type)
        y0 = np.maximum(y.copy(), 0.0)

        a_v0 = y0[ode.IDX_A_V]
        a_u0 = y0[ode.IDX_A_U]

        mu_v = coeff["mu_all_v"]
        mu_u = coeff["mu_all_u"]

        y1 = y0.copy()
        y1[ode.IDX_A_V] = alive_exact(a_v0, mu_v, 1.0)
        y1[ode.IDX_A_U] = alive_exact(a_u0, mu_u, 1.0)

        inc_men = cumulative_increment(coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        inc_pneu_total = cumulative_increment(coeff["lam_pneu_v"], coeff["lam_pneu_u"], a_v0, a_u0, mu_v, mu_u)
        inc_sep = cumulative_increment(coeff["lam_sep_v"], coeff["lam_sep_u"], a_v0, a_u0, mu_v, mu_u)
        inc_aom = cumulative_increment(coeff["lam_aom_v"], coeff["lam_aom_u"], a_v0, a_u0, mu_v, mu_u)
        inc_death = cumulative_increment(coeff["mu_dis_v"], coeff["mu_dis_u"], a_v0, a_u0, mu_v, mu_u)

        inc_pneu_sev = coeff["prop_sev"] * inc_pneu_total
        inc_pneu_nonsev = (1.0 - coeff["prop_sev"]) * inc_pneu_total

        y1[ode.IDX_C_MEN] += inc_men
        y1[ode.IDX_C_PNEU_SEV] += inc_pneu_sev
        y1[ode.IDX_C_PNEU_NONSEV] += inc_pneu_nonsev
        y1[ode.IDX_C_SEP] += inc_sep
        y1[ode.IDX_C_AOM] += inc_aom
        y1[ode.IDX_C_DEATH] += inc_death

        q_cog = float(p["prob_cog"])
        q_hear = float(p["prob_hear"])
        q_epi = float(p["prob_epi"])
        q_hemi = float(p["prob_hemi"])
        q_impl = float(p.get("prob_cochlear_implant", 0.0)) * q_hear
        q_sped = float(p.get("prob_special_education", 0.0)) * q_hear

        y1[ode.IDX_S_COG] = sequela_exact(y0[ode.IDX_S_COG], q_cog, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_HEAR] = sequela_exact(y0[ode.IDX_S_HEAR], q_hear, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_EPI] = sequela_exact(y0[ode.IDX_S_EPI], q_epi, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_HEMI] = sequela_exact(y0[ode.IDX_S_HEMI], q_hemi, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_IMPLANT] = sequela_exact(y0[ode.IDX_S_IMPLANT], q_impl, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_SPED] = sequela_exact(y0[ode.IDX_S_SPED], q_sped, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)

        y1 = np.maximum(y1, 0.0)

        step_data = {
            "Year": (month // 12) + 1,
            "New_Meningitis": y1[ode.IDX_C_MEN] - y0[ode.IDX_C_MEN],
            "New_Pneumo_Sev": y1[ode.IDX_C_PNEU_SEV] - y0[ode.IDX_C_PNEU_SEV],
            "New_Pneumo_NonSev": y1[ode.IDX_C_PNEU_NONSEV] - y0[ode.IDX_C_PNEU_NONSEV],
            "New_Sepsis": y1[ode.IDX_C_SEP] - y0[ode.IDX_C_SEP],
            "New_AOM": y1[ode.IDX_C_AOM] - y0[ode.IDX_C_AOM],
            "New_Deaths_Pneumo": y1[ode.IDX_C_DEATH] - y0[ode.IDX_C_DEATH],
            "Stock_Alive": y1[ode.IDX_A_V] + y1[ode.IDX_A_U],
            "Stock_Dead": birth - (y1[ode.IDX_A_V] + y1[ode.IDX_A_U]),
            "Stock_New_Meningitis": y1[ode.IDX_C_MEN],
            "Stock_New_Pneumo_Sev": y1[ode.IDX_C_PNEU_SEV],
            "Stock_New_Pneumo_NonSev": y1[ode.IDX_C_PNEU_NONSEV],
            "Stock_New_Sepsis": y1[ode.IDX_C_SEP],
            "Stock_New_AOM": y1[ode.IDX_C_AOM],
            "Stock_New_Deaths_Pneumo": y1[ode.IDX_C_DEATH],
            "Stock_Seq_Cognitive": y1[ode.IDX_S_COG],
            "Stock_Seq_Hearing": y1[ode.IDX_S_HEAR],
            "Stock_Seq_Epilepsy": y1[ode.IDX_S_EPI],
            "Stock_Seq_Hemiplegia": y1[ode.IDX_S_HEMI],
            "Stock_Seq_Hearing_Implant": y1[ode.IDX_S_IMPLANT],
            "Stock_Seq_Hearing_SpEd": y1[ode.IDX_S_SPED],
        }
        history.append(step_data)
        y = y1

    return pd.DataFrame(history)


# ============================================================
# Validation helpers
# ============================================================

def compare_exact_to_rk45(scenario_name: str) -> pd.DataFrame:
    base = copy.deepcopy(ode.params_vietnam_micro)
    params, intervention, vaccine_type = scenario_params(base, scenario_name)

    df_num = run_exact_ode_and_build_monthly(params, intervention=intervention, vaccine_type=vaccine_type)
    df_ref = ode.run_ode_and_build_monthly(params, intervention=intervention, vaccine_type=vaccine_type)

    common_cols = [c for c in df_num.columns if c in df_ref.columns and c != "Year"]
    rows = []
    for col in common_cols:
        diff = np.abs(df_num[col].values - df_ref[col].values)
        scale = np.maximum(np.abs(df_ref[col].values), 1.0)
        rows.append({
            "Scenario": scenario_name,
            "Metric": col,
            "MaxAbsDiff": float(diff.max()),
            "MeanAbsDiff": float(diff.mean()),
            "MaxRelDiff": float((diff / scale).max()),
        })
    return pd.DataFrame(rows)


def validate_all_scenarios(output_csv: str = "exact_vs_rk45_validation.csv") -> pd.DataFrame:
    pieces = [compare_exact_to_rk45(s) for s in SCENARIO_MAP]
    out = pd.concat(pieces, ignore_index=True)
    out.to_csv(output_csv, index=False)
    return out


# ============================================================
# Scenario runners and economic summaries
# ============================================================

def yearly_exact_for_scenario(base_epi_params: Dict[str, Any], scenario_name: str) -> pd.DataFrame:
    params, intervention, vaccine_type = scenario_params(base_epi_params, scenario_name)
    monthly = run_exact_ode_and_build_monthly(params, intervention=intervention, vaccine_type=vaccine_type)
    yearly = ode.aggregate_yearly_like_micro(monthly)
    return yearly


def run_exact_all_scenarios_to_csv(outdir: str = "exact_ode_outputs_csv") -> Dict[str, pd.DataFrame]:
    os.makedirs(outdir, exist_ok=True)
    outputs: Dict[str, pd.DataFrame] = {}
    for scenario_name in SCENARIO_MAP:
        df = yearly_exact_for_scenario(copy.deepcopy(ode.params_vietnam_micro), scenario_name)
        outputs[scenario_name] = df
        df.to_csv(os.path.join(outdir, f"exact_{scenario_name}.csv"), index=False)
    return outputs


def build_economic_inputs() -> Tuple[Dict[str, Any], Dict[str, Any]]:
    if qmod is None:
        raise ImportError("QALY_new.py is required for economic summaries")
    cfg = copy.deepcopy(qmod.CONFIG)
    params = qmod.build_parameters_from_hardcoded_table(cfg)
    return params, cfg


def economic_summary_from_yearly(df_yearly: pd.DataFrame, scenario_name: str, econ_params: Dict[str, Any], cfg: Dict[str, Any]) -> Dict[str, float]:
    subset = df_yearly[qmod.EXPECTED_COLS].copy()
    res = qmod.compute_yearly_pv(subset, econ_params, cfg)
    rows = res.to_dict("records")
    vaccine_key = qmod.scenario_vaccine_key(scenario_name)
    qmod.add_vaccine_program_pv(rows, econ_params, vaccine_key, cfg)
    res2 = pd.DataFrame(rows)
    summ = qmod.summarize_results(res2)
    return {
        "PV_total_cost": float(summ["PV_total_cost"]),
        "PV_total_qaly_loss": float(summ["PV_total_qaly_loss"]),
        "Total_cases": float(summ["Total_cases"]),
        "Total_deaths": float(summ["Total_deaths"]),
    }


def incremental_summary(base_yearly: pd.DataFrame, scenario_yearly: pd.DataFrame, scenario_name: str, econ_params: Dict[str, Any], cfg: Dict[str, Any]) -> Dict[str, float]:
    b = economic_summary_from_yearly(base_yearly, "baseline", econ_params, cfg)
    s = economic_summary_from_yearly(scenario_yearly, scenario_name, econ_params, cfg)
    delta_cost = s["PV_total_cost"] - b["PV_total_cost"]
    qaly_gained = b["PV_total_qaly_loss"] - s["PV_total_qaly_loss"]
    cases_averted = b["Total_cases"] - s["Total_cases"]
    deaths_averted = b["Total_deaths"] - s["Total_deaths"]
    icer = np.nan if abs(qaly_gained) < 1e-14 else delta_cost / qaly_gained
    return {
        "Scenario": scenario_name,
        "Delta_cost": float(delta_cost),
        "QALY_gained": float(qaly_gained),
        "Cases_averted": float(cases_averted),
        "Deaths_averted": float(deaths_averted),
        "ICER": float(icer) if not np.isnan(icer) else np.nan,
        "Baseline_cost": float(b["PV_total_cost"]),
        "Scenario_cost": float(s["PV_total_cost"]),
        "Baseline_qaly_loss": float(b["PV_total_qaly_loss"]),
        "Scenario_qaly_loss": float(s["PV_total_qaly_loss"]),
    }


# ============================================================
# Deterministic one-way sensitivity analysis
# ============================================================

def set_nested(obj: Dict[str, Any], path: Tuple[str, ...], value: Any) -> None:
    cur = obj
    for key in path[:-1]:
        cur = cur[key]
    cur[path[-1]] = value


def get_nested(obj: Dict[str, Any], path: Tuple[str, ...]) -> Any:
    cur = obj
    for key in path:
        cur = cur[key]
    return cur


def default_one_way_specs(primary_scenario: str) -> List[Dict[str, Any]]:
    vaccine_key = qmod.scenario_vaccine_key(primary_scenario) if qmod is not None else None
    specs: List[Dict[str, Any]] = [
        {"label": "inc_pneumonia", "domain": "epi", "path": ("inc_pneumonia",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "inc_meningitis", "domain": "epi", "path": ("inc_meningitis",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "inc_sepsis", "domain": "epi", "path": ("inc_sepsis",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "inc_aom_all", "domain": "epi", "path": ("inc_aom_all",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "cfr_meningitis", "domain": "epi", "path": ("cfr_meningitis",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "cfr_pneumo_severe", "domain": "epi", "path": ("cfr_pneumo_severe",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "cfr_sepsis", "domain": "epi", "path": ("cfr_sepsis",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "prop_pneumo_severe", "domain": "epi", "path": ("prop_pneumo_severe",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "prop_aom_spn", "domain": "epi", "path": ("prop_aom_spn",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "vaccine_coverage", "domain": "epi", "path": ("vaccine_coverage",), "low": 0.85, "high": 0.95, "mode": "abs"},
        {"label": "prob_cog", "domain": "epi", "path": ("prob_cog",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "prob_hear", "domain": "epi", "path": ("prob_hear",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "prob_epi", "domain": "epi", "path": ("prob_epi",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "prob_hemi", "domain": "epi", "path": ("prob_hemi",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "discount_rate", "domain": "econ", "path": ("discount_rate",), "low": 0.0, "high": 0.05, "mode": "abs"},
        {"label": "life_expectancy", "domain": "econ", "path": ("life_expectancy",), "low": 70.0, "high": 78.0, "mode": "abs"},
        {"label": "service_delivery_per_dose", "domain": "econ", "path": ("service_delivery_per_dose",), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "wastage_rate", "domain": "econ", "path": ("wastage_rate",), "low": 0.02, "high": 0.10, "mode": "abs"},
        {"label": "cognitive_annual_cost", "domain": "econ", "path": ("costs_sequelae", "cognitive_annual"), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "cochlear_implant_cost", "domain": "econ", "path": ("costs_sequelae", "cochlear_implant_case"), "low": 0.8, "high": 1.2, "mode": "mult"},
        {"label": "utility_cognitive", "domain": "econ", "path": ("utilities", "cognitive"), "low": 0.55, "high": 0.70, "mode": "abs"},
        {"label": "utility_hearing", "domain": "econ", "path": ("utilities", "hearing"), "low": 0.85, "high": 0.97, "mode": "abs"},
        {"label": "utility_hemiplegia", "domain": "econ", "path": ("utilities", "hemiplegia"), "low": 0.30, "high": 0.50, "mode": "abs"},
    ]
    if vaccine_key is not None:
        specs.append({"label": f"vaccine_price_{vaccine_key}", "domain": "econ", "path": ("vaccine_prices", vaccine_key), "low": 0.8, "high": 1.2, "mode": "mult"})
    return specs


def apply_spec_value(base_value: Any, spec: Dict[str, Any], side: str) -> Any:
    target = spec[side]
    if spec.get("mode", "mult") == "mult":
        return float(base_value) * float(target)
    return float(target)


def one_way_sensitivity_exact(
    primary_scenario: str = "pcv13_dynamic",
    specs: Optional[List[Dict[str, Any]]] = None,
    output_csv: str = "one_way_sensitivity_exact.csv",
) -> pd.DataFrame:
    if qmod is None:
        raise ImportError("QALY_new.py is required for one-way sensitivity analysis")

    if specs is None:
        specs = default_one_way_specs(primary_scenario)

    base_epi = copy.deepcopy(ode.params_vietnam_micro)
    base_econ, cfg = build_economic_inputs()

    rows: List[Dict[str, Any]] = []
    for spec in specs:
        for side in ["low", "high"]:
            epi = copy.deepcopy(base_epi)
            econ = copy.deepcopy(base_econ)
            cfg_local = copy.deepcopy(cfg)

            if spec["domain"] == "epi":
                old = get_nested(epi, spec["path"])
                new = apply_spec_value(old, spec, side)
                set_nested(epi, spec["path"], new)
            elif spec["domain"] == "econ":
                old = get_nested(econ, spec["path"])
                new = apply_spec_value(old, spec, side)
                set_nested(econ, spec["path"], new)
            elif spec["domain"] == "cfg":
                old = get_nested(cfg_local, spec["path"])
                new = apply_spec_value(old, spec, side)
                set_nested(cfg_local, spec["path"], new)
            else:
                raise ValueError(f"Unknown spec domain: {spec['domain']}")

            base_yearly = yearly_exact_for_scenario(epi, "baseline")
            scen_yearly = yearly_exact_for_scenario(epi, primary_scenario)
            summary = incremental_summary(base_yearly, scen_yearly, primary_scenario, econ, cfg_local)

            rows.append({
                "Parameter": spec["label"],
                "Bound": side,
                **summary,
            })

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)
    return out


# ============================================================
# Probabilistic sensitivity analysis
# ============================================================

def beta_sample_from_mean(mean: float, kappa: float, rng: np.random.Generator) -> float:
    mean = min(max(mean, 1e-6), 1.0 - 1e-6)
    alpha = mean * kappa
    beta = (1.0 - mean) * kappa
    return float(rng.beta(alpha, beta))


def lognormal_sample_from_mean_cv(mean: float, cv: float, rng: np.random.Generator) -> float:
    mean = max(float(mean), 1e-9)
    cv = max(float(cv), 1e-6)
    sigma2 = math.log(1.0 + cv * cv)
    sigma = math.sqrt(sigma2)
    mu = math.log(mean) - 0.5 * sigma2
    return float(rng.lognormal(mean=mu, sigma=sigma))


def psa_exact(
    primary_scenario: str = "pcv13_dynamic",
    n_samples: int = 1000,
    seed: int = 20260407,
    output_csv: str = "psa_exact.csv",
) -> pd.DataFrame:
    if qmod is None:
        raise ImportError("QALY_new.py is required for probabilistic sensitivity analysis")

    rng = np.random.default_rng(seed)
    base_epi = copy.deepcopy(ode.params_vietnam_micro)
    base_econ, cfg = build_economic_inputs()

    rows: List[Dict[str, Any]] = []
    for i in range(n_samples):
        epi = copy.deepcopy(base_epi)
        econ = copy.deepcopy(base_econ)
        cfg_local = copy.deepcopy(cfg)

        epi["inc_pneumonia"] = lognormal_sample_from_mean_cv(base_epi["inc_pneumonia"], 0.20, rng)
        epi["inc_meningitis"] = lognormal_sample_from_mean_cv(base_epi["inc_meningitis"], 0.25, rng)
        epi["inc_sepsis"] = lognormal_sample_from_mean_cv(base_epi["inc_sepsis"], 0.20, rng)
        epi["inc_aom_all"] = lognormal_sample_from_mean_cv(base_epi["inc_aom_all"], 0.20, rng)

        epi["prop_pneumo_severe"] = beta_sample_from_mean(base_epi["prop_pneumo_severe"], 200.0, rng)
        epi["prop_aom_spn"] = beta_sample_from_mean(base_epi["prop_aom_spn"], 200.0, rng)
        epi["cfr_meningitis"] = beta_sample_from_mean(base_epi["cfr_meningitis"], 300.0, rng)
        epi["cfr_pneumo_severe"] = beta_sample_from_mean(base_epi["cfr_pneumo_severe"], 600.0, rng)
        epi["cfr_sepsis"] = beta_sample_from_mean(base_epi["cfr_sepsis"], 120.0, rng)
        epi["vaccine_coverage"] = beta_sample_from_mean(base_epi["vaccine_coverage"], 120.0, rng)
        epi["prob_cog"] = beta_sample_from_mean(base_epi["prob_cog"], 120.0, rng)
        epi["prob_hear"] = beta_sample_from_mean(base_epi["prob_hear"], 120.0, rng)
        epi["prob_epi"] = beta_sample_from_mean(base_epi["prob_epi"], 120.0, rng)
        epi["prob_hemi"] = beta_sample_from_mean(base_epi["prob_hemi"], 120.0, rng)

        for k in list(econ["vaccine_prices"].keys()):
            econ["vaccine_prices"][k] = lognormal_sample_from_mean_cv(econ["vaccine_prices"][k], 0.15, rng)
        econ["service_delivery_per_dose"] = lognormal_sample_from_mean_cv(econ["service_delivery_per_dose"], 0.15, rng)
        econ["ae_cost_per_case"] = lognormal_sample_from_mean_cv(econ["ae_cost_per_case"], 0.20, rng)
        econ["wastage_rate"] = beta_sample_from_mean(econ["wastage_rate"], 400.0, rng)
        econ["costs_sequelae"]["cognitive_annual"] = lognormal_sample_from_mean_cv(econ["costs_sequelae"]["cognitive_annual"], 0.20, rng)
        econ["costs_sequelae"]["hearing_case"] = lognormal_sample_from_mean_cv(econ["costs_sequelae"]["hearing_case"], 0.20, rng)
        econ["costs_sequelae"]["epilepsy_case"] = lognormal_sample_from_mean_cv(econ["costs_sequelae"]["epilepsy_case"], 0.20, rng)
        econ["costs_sequelae"]["hemiplegia_case"] = lognormal_sample_from_mean_cv(econ["costs_sequelae"]["hemiplegia_case"], 0.20, rng)
        econ["costs_sequelae"]["cochlear_implant_case"] = lognormal_sample_from_mean_cv(econ["costs_sequelae"]["cochlear_implant_case"], 0.20, rng)
        econ["costs_sequelae"]["special_education_annual"] = lognormal_sample_from_mean_cv(econ["costs_sequelae"]["special_education_annual"], 0.20, rng)
        econ["utilities"]["cognitive"] = beta_sample_from_mean(econ["utilities"]["cognitive"], 200.0, rng)
        econ["utilities"]["hearing"] = beta_sample_from_mean(econ["utilities"]["hearing"], 200.0, rng)
        econ["utilities"]["epilepsy"] = beta_sample_from_mean(econ["utilities"]["epilepsy"], 200.0, rng)
        econ["utilities"]["hemiplegia"] = beta_sample_from_mean(econ["utilities"]["hemiplegia"], 200.0, rng)

        base_yearly = yearly_exact_for_scenario(epi, "baseline")
        scen_yearly = yearly_exact_for_scenario(epi, primary_scenario)
        summary = incremental_summary(base_yearly, scen_yearly, primary_scenario, econ, cfg_local)

        rows.append({
            "Sample": i + 1,
            **summary,
        })

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)
    return out


def ceac_from_psa(psa_df: pd.DataFrame, wtp_values: Optional[List[float]] = None) -> pd.DataFrame:
    if wtp_values is None:
        wtp_values = [0, 500, 1000, 2000, 3000, 5000, 10000]
    rows: List[Dict[str, float]] = []
    for wtp in wtp_values:
        nmb = psa_df["QALY_gained"] * float(wtp) - psa_df["Delta_cost"]
        rows.append({
            "WTP": float(wtp),
            "Pr_CostEffective": float((nmb > 0).mean()),
            "Mean_NMB": float(nmb.mean()),
        })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print("Running exact ODE outputs for all scenarios.")
    run_exact_all_scenarios_to_csv()

    print("Running validation against RK45.")
    validation = validate_all_scenarios()
    print(validation.groupby("Scenario")["MaxAbsDiff"].max())

    if qmod is not None:
        print("Running one-way sensitivity for PCV13 dynamic.")
        ow = one_way_sensitivity_exact(primary_scenario="pcv13_dynamic")
        print(ow[["Parameter", "Bound", "Delta_cost", "QALY_gained", "ICER"]].head())

        print("Running PSA for PCV13 dynamic.")
        psa = psa_exact(primary_scenario="pcv13_dynamic", n_samples=200)
        ceac = ceac_from_psa(psa)
        ceac.to_csv("ceac_exact.csv", index=False)
        print(ceac)
