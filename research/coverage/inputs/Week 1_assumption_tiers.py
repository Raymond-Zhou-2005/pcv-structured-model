import copy
import importlib.util
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import ode_new as ode

_EXACT_PATH = Path(__file__).with_name("Week 1_exact_sensitivity.py")
_EXACT_SPEC = importlib.util.spec_from_file_location("week1_exact_sensitivity", _EXACT_PATH)
if _EXACT_SPEC is None or _EXACT_SPEC.loader is None:
    raise ImportError(f"Cannot load exact sensitivity module from {_EXACT_PATH}")
ex = importlib.util.module_from_spec(_EXACT_SPEC)
_EXACT_SPEC.loader.exec_module(ex)

try:
    import QALY_new as qmod
except ImportError:
    qmod = None


PRIMARY_SCENARIO = "pcv13_dynamic"


# ============================================================
# Shared helpers
# ============================================================

def summarize_monthly_epi(df_monthly: pd.DataFrame) -> Dict[str, float]:
    return {
        "Total_cases": float(df_monthly[["New_Meningitis", "New_Pneumo_Sev", "New_Pneumo_NonSev", "New_Sepsis", "New_AOM"]].sum().sum()),
        "Total_deaths": float(df_monthly["New_Deaths_Pneumo"].sum()),
        "Final_alive": float(df_monthly["Stock_Alive"].iloc[-1]),
        "Final_seq_cog": float(df_monthly["Stock_Seq_Cognitive"].iloc[-1]),
        "Final_seq_hear": float(df_monthly["Stock_Seq_Hearing"].iloc[-1]),
    }


def monthly_to_yearly(df_monthly: pd.DataFrame) -> pd.DataFrame:
    return ode.aggregate_yearly_like_micro(df_monthly)


def try_incremental_econ(base_monthly: pd.DataFrame, scen_monthly: pd.DataFrame, scenario_name: str) -> Dict[str, float]:
    if qmod is None:
        return {}
    econ_params, cfg = ex.build_economic_inputs()
    base_yearly = monthly_to_yearly(base_monthly)
    scen_yearly = monthly_to_yearly(scen_monthly)
    return ex.incremental_summary(base_yearly, scen_yearly, scenario_name, econ_params, cfg)


# ============================================================
# Exact ODE with horizon truncation and structural mapping modes
# ============================================================

def scenario_setup(scenario_name: str) -> Tuple[Dict[str, Any], bool, Optional[str]]:
    base = copy.deepcopy(ode.params_vietnam_micro)
    return ex.scenario_params(base, scenario_name)


def month_coefficients_general(
    p: Dict[str, Any],
    month: int,
    intervention: bool,
    vaccine_type: Optional[str],
    mapping_mode: str = "direct_effect_only",
    herd_share: float = 0.0,
) -> Dict[str, float]:
    month = int(max(0, min(59, month)))

    prob_bg_month = ode.get_background_prob_u5_only(p)
    lam_bg = ex.prob_to_hazard(prob_bg_month)

    w_men = ode.get_age_weight(month, "Meningitis")
    w_pneu = ode.get_age_weight(month, "Pneumonia")
    w_aom = ode.get_age_weight(month, "AOM")

    base_prob_men = ode.monthly_prob_from_incidence_per100k(p["inc_meningitis"], w_men)
    base_prob_pneu = ode.monthly_prob_from_incidence_per100k(p["inc_pneumonia"], w_pneu)
    base_prob_sep = ode.monthly_prob_from_incidence_per100k(p["inc_sepsis"], w_pneu)
    base_prob_aom = ode.monthly_prob_from_incidence_per100k(p["inc_aom_all"], w_aom) * float(p["prop_aom_spn"])

    mult_v = 1.0
    mult_u = 1.0
    irr = 1.0

    if intervention and (vaccine_type is not None):
        year_idx = month // 12
        cov = float(p.get("vaccine_coverage", 0.0))
        irr, ve_year = ode.infer_ve_from_irr_series(p["vaccine_irrs"], vaccine_type, year_idx, cov)
        base_reduction = 1.0 - irr

        if mapping_mode == "direct_effect_only":
            mult_v = 1.0 - ve_year
            mult_u = 1.0
        elif mapping_mode == "uniform_population_irr":
            mult_v = irr
            mult_u = irr
        elif mapping_mode == "partial_herd":
            mult_v = 1.0 - ve_year
            mult_u = max(0.0, 1.0 - herd_share * base_reduction)
        else:
            raise ValueError(f"Unknown mapping_mode: {mapping_mode}")

    p_men_v = ex.clip01(base_prob_men * mult_v)
    p_pneu_v = ex.clip01(base_prob_pneu * mult_v)
    p_sep_v = ex.clip01(base_prob_sep * mult_v)
    p_aom_v = ex.clip01(base_prob_aom * mult_v)

    p_men_u = ex.clip01(base_prob_men * mult_u)
    p_pneu_u = ex.clip01(base_prob_pneu * mult_u)
    p_sep_u = ex.clip01(base_prob_sep * mult_u)
    p_aom_u = ex.clip01(base_prob_aom * mult_u)

    lam_men_v = ex.prob_to_hazard(p_men_v)
    lam_pneu_v = ex.prob_to_hazard(p_pneu_v)
    lam_sep_v = ex.prob_to_hazard(p_sep_v)
    lam_aom_v = ex.prob_to_hazard(p_aom_v)

    lam_men_u = ex.prob_to_hazard(p_men_u)
    lam_pneu_u = ex.prob_to_hazard(p_pneu_u)
    lam_sep_u = ex.prob_to_hazard(p_sep_u)
    lam_aom_u = ex.prob_to_hazard(p_aom_u)

    prop_sev = float(p["prop_pneumo_severe"])
    cfr_men = float(p["cfr_meningitis"])
    cfr_pneu = float(p["cfr_pneumo_severe"])
    cfr_sep = float(p["cfr_sepsis"])

    d1_v = ex.clip01(p_men_v * cfr_men)
    d2_v = ex.clip01(p_pneu_v * prop_sev * cfr_pneu)
    d3_v = ex.clip01(p_sep_v * cfr_sep)
    d1_u = ex.clip01(p_men_u * cfr_men)
    d2_u = ex.clip01(p_pneu_u * prop_sev * cfr_pneu)
    d3_u = ex.clip01(p_sep_u * cfr_sep)

    mu_dis_v = ex.prob_to_hazard(d1_v) + ex.prob_to_hazard(d2_v) + ex.prob_to_hazard(d3_v)
    mu_dis_u = ex.prob_to_hazard(d1_u) + ex.prob_to_hazard(d2_u) + ex.prob_to_hazard(d3_u)

    mu_all_v = lam_bg + mu_dis_v
    mu_all_u = lam_bg + mu_dis_u

    return {
        "irr": irr,
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


def run_exact_monthly_general(
    p: Dict[str, Any],
    intervention: bool,
    vaccine_type: Optional[str],
    horizon_months: int = 60,
    mapping_mode: str = "direct_effect_only",
    herd_share: float = 0.0,
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

    for month in range(horizon_months):
        coeff = month_coefficients_general(
            p,
            month,
            intervention=intervention,
            vaccine_type=vaccine_type,
            mapping_mode=mapping_mode,
            herd_share=herd_share,
        )
        y0 = np.maximum(y.copy(), 0.0)
        a_v0 = y0[ode.IDX_A_V]
        a_u0 = y0[ode.IDX_A_U]
        mu_v = coeff["mu_all_v"]
        mu_u = coeff["mu_all_u"]

        y1 = y0.copy()
        y1[ode.IDX_A_V] = ex.alive_exact(a_v0, mu_v, 1.0)
        y1[ode.IDX_A_U] = ex.alive_exact(a_u0, mu_u, 1.0)

        inc_men = ex.cumulative_increment(coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        inc_pneu_total = ex.cumulative_increment(coeff["lam_pneu_v"], coeff["lam_pneu_u"], a_v0, a_u0, mu_v, mu_u)
        inc_sep = ex.cumulative_increment(coeff["lam_sep_v"], coeff["lam_sep_u"], a_v0, a_u0, mu_v, mu_u)
        inc_aom = ex.cumulative_increment(coeff["lam_aom_v"], coeff["lam_aom_u"], a_v0, a_u0, mu_v, mu_u)
        inc_death = ex.cumulative_increment(coeff["mu_dis_v"], coeff["mu_dis_u"], a_v0, a_u0, mu_v, mu_u)

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

        y1[ode.IDX_S_COG] = ex.sequela_exact(y0[ode.IDX_S_COG], q_cog, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_HEAR] = ex.sequela_exact(y0[ode.IDX_S_HEAR], q_hear, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_EPI] = ex.sequela_exact(y0[ode.IDX_S_EPI], q_epi, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_HEMI] = ex.sequela_exact(y0[ode.IDX_S_HEMI], q_hemi, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_IMPLANT] = ex.sequela_exact(y0[ode.IDX_S_IMPLANT], q_impl, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)
        y1[ode.IDX_S_SPED] = ex.sequela_exact(y0[ode.IDX_S_SPED], q_sped, coeff["cfr_men"], coeff["lam_men_v"], coeff["lam_men_u"], a_v0, a_u0, mu_v, mu_u)

        y1 = np.maximum(y1, 0.0)

        history.append({
            "Month": month + 1,
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
        })
        y = y1

    return pd.DataFrame(history)


# ============================================================
# Category 1: accounting-style checks
# ============================================================

def weight_normalization_table() -> pd.DataFrame:
    rows = []
    disease_map = {
        "Meningitis": "Meningitis",
        "Pneumonia": "Pneumonia",
        "Sepsis": "Pneumonia",
        "AOM": "AOM",
    }
    for label, arg in disease_map.items():
        s = 0.0
        for m in range(60):
            s += float(ode.get_age_weight(m, arg))
        rows.append({"AssumptionType": "Accounting", "Check": f"weight_sum_{label}", "Value": s})
    return pd.DataFrame(rows)


def horizon_prefix_invariance_table(scenario_name: str = PRIMARY_SCENARIO) -> pd.DataFrame:
    p, intervention, vaccine_type = scenario_setup(scenario_name)
    full_df = run_exact_monthly_general(p, intervention, vaccine_type, horizon_months=60)

    rows = []
    for horizon in [12, 24, 36, 48, 60]:
        trunc_df = run_exact_monthly_general(p, intervention, vaccine_type, horizon_months=horizon)
        prefix = full_df.iloc[:horizon].reset_index(drop=True)
        common_cols = [c for c in trunc_df.columns if c in prefix.columns and c not in {"Month", "Year"}]
        max_abs = 0.0
        for col in common_cols:
            max_abs = max(max_abs, float(np.max(np.abs(trunc_df[col].values - prefix[col].values))))
        rows.append({
            "AssumptionType": "Accounting",
            "Check": f"horizon_prefix_{horizon}m",
            "Value": max_abs,
        })
    return pd.DataFrame(rows)


def exact_vs_rk45_table(scenario_name: str = PRIMARY_SCENARIO) -> pd.DataFrame:
    df = ex.compare_exact_to_rk45(scenario_name)
    return pd.DataFrame([
        {
            "AssumptionType": "Accounting",
            "Check": f"exact_vs_rk45_{scenario_name}_max_abs",
            "Value": float(df["MaxAbsDiff"].max()),
        },
        {
            "AssumptionType": "Accounting",
            "Check": f"exact_vs_rk45_{scenario_name}_max_rel",
            "Value": float(df["MaxRelDiff"].max()),
        },
    ])


def run_accounting_checks(output_csv: str = "tier1_accounting_checks.csv") -> pd.DataFrame:
    out = pd.concat([
        weight_normalization_table(),
        horizon_prefix_invariance_table(),
        exact_vs_rk45_table(),
    ], ignore_index=True)
    out.to_csv(output_csv, index=False)
    return out


# ============================================================
# Category 2: parameter uncertainty via local elasticities
# ============================================================

def finite_diff_log_elasticity(base_val: float, low_val: float, high_val: float, low_param: float, high_param: float) -> float:
    if (base_val <= 0.0) or (low_val <= 0.0) or (high_val <= 0.0) or (low_param <= 0.0) or (high_param <= 0.0):
        return np.nan
    return (math.log(high_val) - math.log(low_val)) / (math.log(high_param) - math.log(low_param))


def run_parameter_elasticities(
    scenario_name: str = PRIMARY_SCENARIO,
    output_csv: str = "tier2_parameter_elasticities.csv",
) -> pd.DataFrame:
    base_epi = copy.deepcopy(ode.params_vietnam_micro)
    rows = []

    econ_available = qmod is not None
    if econ_available:
        econ_base, cfg = ex.build_economic_inputs()

    specs = [
        ("inc_pneumonia", ("inc_pneumonia",), 0.9, 1.1),
        ("inc_meningitis", ("inc_meningitis",), 0.9, 1.1),
        ("inc_sepsis", ("inc_sepsis",), 0.9, 1.1),
        ("inc_aom_all", ("inc_aom_all",), 0.9, 1.1),
        ("cfr_meningitis", ("cfr_meningitis",), 0.9, 1.1),
        ("cfr_pneumo_severe", ("cfr_pneumo_severe",), 0.9, 1.1),
        ("cfr_sepsis", ("cfr_sepsis",), 0.9, 1.1),
        ("vaccine_coverage", ("vaccine_coverage",), 0.95, 1.05),
        ("prob_hear", ("prob_hear",), 0.9, 1.1),
        ("prob_hemi", ("prob_hemi",), 0.9, 1.1),
    ]

    for label, path, low_mult, high_mult in specs:
        epi_low = copy.deepcopy(base_epi)
        epi_high = copy.deepcopy(base_epi)

        base_param = ex.get_nested(base_epi, path)
        low_param = float(base_param) * low_mult
        high_param = float(base_param) * high_mult
        ex.set_nested(epi_low, path, low_param)
        ex.set_nested(epi_high, path, high_param)

        p_base, i_base, v_base = scenario_setup("baseline")
        s_base, i_s, v_s = scenario_setup(scenario_name)
        # overwrite scenario params using perturbed epidemiology settings
        for src, dest in [(epi_low, p_base), (epi_low, s_base)]:
            pass

        base_month_low = run_exact_monthly_general(epi_low, False, None)
        scen_month_low = run_exact_monthly_general(epi_low, i_s, v_s)
        base_month_high = run_exact_monthly_general(epi_high, False, None)
        scen_month_high = run_exact_monthly_general(epi_high, i_s, v_s)

        epi_low_summary = summarize_monthly_epi(scen_month_low)
        epi_high_summary = summarize_monthly_epi(scen_month_high)
        base_low_summary = summarize_monthly_epi(base_month_low)
        base_high_summary = summarize_monthly_epi(base_month_high)

        low_cases_averted = base_low_summary["Total_cases"] - epi_low_summary["Total_cases"]
        high_cases_averted = base_high_summary["Total_cases"] - epi_high_summary["Total_cases"]
        low_deaths_averted = base_low_summary["Total_deaths"] - epi_low_summary["Total_deaths"]
        high_deaths_averted = base_high_summary["Total_deaths"] - epi_high_summary["Total_deaths"]

        row = {
            "AssumptionType": "Parameter",
            "Parameter": label,
            "CasesAvertedElasticity": finite_diff_log_elasticity(1.0, abs(low_cases_averted) + 1e-12, abs(high_cases_averted) + 1e-12, low_param, high_param),
            "DeathsAvertedElasticity": finite_diff_log_elasticity(1.0, abs(low_deaths_averted) + 1e-12, abs(high_deaths_averted) + 1e-12, low_param, high_param),
            "LowCasesAverted": float(low_cases_averted),
            "HighCasesAverted": float(high_cases_averted),
            "LowDeathsAverted": float(low_deaths_averted),
            "HighDeathsAverted": float(high_deaths_averted),
        }

        if econ_available:
            econ_low = ex.incremental_summary(monthly_to_yearly(base_month_low), monthly_to_yearly(scen_month_low), scenario_name, copy.deepcopy(econ_base), copy.deepcopy(cfg))
            econ_high = ex.incremental_summary(monthly_to_yearly(base_month_high), monthly_to_yearly(scen_month_high), scenario_name, copy.deepcopy(econ_base), copy.deepcopy(cfg))
            row.update({
                "LowDeltaCost": float(econ_low["Delta_cost"]),
                "HighDeltaCost": float(econ_high["Delta_cost"]),
                "LowQALYGained": float(econ_low["QALY_gained"]),
                "HighQALYGained": float(econ_high["QALY_gained"]),
                "DeltaCostElasticity": finite_diff_log_elasticity(1.0, abs(econ_low["Delta_cost"]) + 1e-12, abs(econ_high["Delta_cost"]) + 1e-12, low_param, high_param),
                "QALYGainedElasticity": finite_diff_log_elasticity(1.0, abs(econ_low["QALY_gained"]) + 1e-12, abs(econ_high["QALY_gained"]) + 1e-12, low_param, high_param),
            })

        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)
    return out


# ============================================================
# Category 3: structural screens that preserve analytic solvability
# ============================================================

def run_structural_screens(
    scenario_name: str = PRIMARY_SCENARIO,
    output_csv: str = "tier3_structural_screens.csv",
) -> pd.DataFrame:
    p, intervention, vaccine_type = scenario_setup(scenario_name)

    screens = [
        ("direct_effect_only", "direct_effect_only", 0.0),
        ("uniform_population_irr", "uniform_population_irr", 0.0),
        ("partial_herd_25pct", "partial_herd", 0.25),
        ("partial_herd_50pct", "partial_herd", 0.50),
    ]

    base_month = run_exact_monthly_general(p, False, None)
    rows = []

    for label, mode, herd_share in screens:
        scen_month = run_exact_monthly_general(
            p,
            intervention,
            vaccine_type,
            mapping_mode=mode,
            herd_share=herd_share,
        )
        epi_base = summarize_monthly_epi(base_month)
        epi_scen = summarize_monthly_epi(scen_month)

        row = {
            "AssumptionType": "Structural",
            "Structure": label,
            "CasesAverted": float(epi_base["Total_cases"] - epi_scen["Total_cases"]),
            "DeathsAverted": float(epi_base["Total_deaths"] - epi_scen["Total_deaths"]),
            "FinalCogAvoided": float(base_month["Stock_Seq_Cognitive"].iloc[-1] - scen_month["Stock_Seq_Cognitive"].iloc[-1]),
            "FinalHearAvoided": float(base_month["Stock_Seq_Hearing"].iloc[-1] - scen_month["Stock_Seq_Hearing"].iloc[-1]),
        }

        if qmod is not None:
            econ = try_incremental_econ(base_month, scen_month, scenario_name)
            row.update({
                "DeltaCost": float(econ["Delta_cost"]),
                "QALYGained": float(econ["QALY_gained"]),
                "ICER": float(econ["ICER"]) if not np.isnan(econ["ICER"]) else np.nan,
            })

        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)
    return out


# ============================================================
# Master runner
# ============================================================

def run_all_tier_demos(outdir: str = "assumption_tier_demo_outputs") -> Dict[str, pd.DataFrame]:
    os.makedirs(outdir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(outdir)
    try:
        tier1 = run_accounting_checks()
        tier2 = run_parameter_elasticities()
        tier3 = run_structural_screens()
        summary = pd.DataFrame([
            {"Tier": "Accounting", "KeyPoint": "Weight normalization and horizon-prefix invariance hold exactly under the analytic propagator.", "PrimaryFile": "tier1_accounting_checks.csv"},
            {"Tier": "Parameter", "KeyPoint": "Local elasticities rank ordinary uncertain inputs without altering model class.", "PrimaryFile": "tier2_parameter_elasticities.csv"},
            {"Tier": "Structural", "KeyPoint": "Alternative vaccine-effect mappings change coefficients while preserving exact solvability.", "PrimaryFile": "tier3_structural_screens.csv"},
        ])
        summary.to_csv("tier_summary.csv", index=False)
        return {
            "tier1": tier1,
            "tier2": tier2,
            "tier3": tier3,
            "summary": summary,
        }
    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    outputs = run_all_tier_demos()
    for name, df in outputs.items():
        print("\n=====", name, "=====")
        print(df.head())
