"""
Explicit age-time PDE calendar model with refined vaccine mechanics.

This file uses the same solver-based calendar PDE framework as
Week 5_calendar_pde_original.py, then replaces the original vaccine logic with
the Week 2 refinements:

1. New recorded-window births are split into program-enrolled and not-enrolled
   children by vaccine_coverage.
2. The epidemiological schedule is represented by direct protection turning on
   at age 4 months, after the 2- and 4-month primary-series doses.
3. Direct VE is endpoint-specific: IPD, pneumonia, and AOM no longer share one
   forced VE.
4. The retained annual total-program IRR is decomposed into a herd multiplier
   and a direct-protection multiplier.

Warm-up is baseline-only and not vaccinated. This means the recorded window
starts at program introduction. Existing children aged 1--59 months are carried
in from the no-program population, while only new recorded-window births enter
the vaccine program. The default run then continues for ten years: years 1--5
are rollout transition years, and years 6--10 approximate the mature program
state after the pre-program under-five cohorts have aged out.
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
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "Outputs" / "Week 5_calendar_pde_refined"


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


base = _load_peer("Week 5_calendar_pde_original.py", "week5_calendar_pde_original_base")
DEFAULT_RECORD_MONTHS = base.DEFAULT_RECORD_MONTHS
STABLE_END_YEAR = base.STABLE_END_YEAR

params_vietnam_micro = copy.deepcopy(base.params_vietnam_micro)


# Evidence distilled in the serotype evidence workbook used for this refinement.
#
# The supplied Vietnam serotype percentages are pooled summary percentages, not a
# strict compositional distribution: the visible categories sum to 109%. We
# therefore keep those percentages as evidence weights, normalize only within the
# rows used for weighting, and do not force them to be a literal 100% serotype
# distribution.
PCV10_SII_EVIDENCE = {
    "raw_visible_total_pct": 109.0,
    "raw_common_vietnam_vaccine_serotype_pct": 76.0,
    "raw_nvt_plus_nontypeable_pct": 33.0,
    "normalized_both_products_cover_visible_pct": 69.72477064220183,
    "weighted_siipl_over_pcv13_ratio_mean": 0.8414473684210527,
    "weighted_siipl_over_pcv13_gmc_ratio": 0.8203164085517026,
}


# Translation choices for PCV10-SII direct VE.
#
# These are relative to the PCV13 endpoint-specific direct VE profile. The
# central value uses the weighted GMC ratio from the Vietnam-relevant shared
# serotypes visible in the supplied materials. GMC is an immunogenicity proxy,
# not a direct efficacy estimate, so the other assumptions are kept for easy
# sensitivity analysis.
PCV10_SII_RELATIVE_VE_ASSUMPTIONS = {
    "conservative_visible_coverage": 0.6972477064220183,
    "central_weighted_gmc": PCV10_SII_EVIDENCE["weighted_siipl_over_pcv13_gmc_ratio"],
    "weighted_ratio_mean": PCV10_SII_EVIDENCE["weighted_siipl_over_pcv13_ratio_mean"],
    "optimistic_near_pcv13": 0.90,
}
DEFAULT_PCV10_SII_ASSUMPTION = "central_weighted_gmc"

BASE_DIRECT_EFFECTIVENESS = {
    "PCV10-GSK": {"ipd": 0.50, "pneumonia": 0.22, "aom": 0.53},
    "PCV13": {"ipd": 0.602, "pneumonia": 0.19, "aom": 0.46},
}


def make_pcv10_sii_direct_effectiveness(assumption: str = DEFAULT_PCV10_SII_ASSUMPTION) -> Dict[str, float]:
    """Build PCV10-SII endpoint VE from Vietnam serotype-weighted evidence."""
    if assumption not in PCV10_SII_RELATIVE_VE_ASSUMPTIONS:
        valid = ", ".join(sorted(PCV10_SII_RELATIVE_VE_ASSUMPTIONS))
        raise ValueError(f"Unknown PCV10-SII assumption: {assumption}. Valid choices: {valid}")
    multiplier = float(PCV10_SII_RELATIVE_VE_ASSUMPTIONS[assumption])
    pcv13 = BASE_DIRECT_EFFECTIVENESS["PCV13"]
    return {endpoint: max(0.0, min(1.0, float(value) * multiplier)) for endpoint, value in pcv13.items()}


def apply_pcv10_sii_assumption(p: Dict, assumption: str = DEFAULT_PCV10_SII_ASSUMPTION) -> Dict:
    """Return a parameter copy with a PCV10-SII-specific direct VE profile."""
    out = copy.deepcopy(p)
    direct = copy.deepcopy(out.get("vaccine_direct_effectiveness", {}))
    direct["PCV10-SII"] = make_pcv10_sii_direct_effectiveness(assumption)
    out["vaccine_direct_effectiveness"] = direct
    out["vaccine_direct_effectiveness_bridge"] = {}
    out["pcv10_sii_relative_ve_assumption"] = assumption
    out["pcv10_sii_relative_ve_multiplier"] = float(PCV10_SII_RELATIVE_VE_ASSUMPTIONS[assumption])
    return out


params_vietnam_micro.update({
    "vaccine_schedule_months": [2, 4, 12],
    "direct_protection_start_month": 4,
    "vaccine_direct_effectiveness": copy.deepcopy(BASE_DIRECT_EFFECTIVENESS),
    "vaccine_direct_effectiveness_bridge": {},
})
params_vietnam_micro = apply_pcv10_sii_assumption(params_vietnam_micro, DEFAULT_PCV10_SII_ASSUMPTION)

OUTPUT_COLUMNS = list(base.OUTPUT_COLUMNS)
SEQ_KEYS = list(base.SEQ_KEYS)
FLOW_NAMES = list(base.FLOW_NAMES)


def _with_record_months(p: Dict, record_months: int) -> Dict:
    out = copy.deepcopy(p)
    out["record_months"] = int(record_months)
    return out


def _scenario_specs(
    pcv10_sii_assumption: str = DEFAULT_PCV10_SII_ASSUMPTION,
    record_months: int = DEFAULT_RECORD_MONTHS,
) -> List[Tuple[str, Dict, bool, Optional[str], str]]:
    p_base = apply_pcv10_sii_assumption(params_vietnam_micro, pcv10_sii_assumption)
    p_base = _with_record_months(p_base, record_months)
    specs: List[Tuple[str, Dict, bool, Optional[str], str]] = [
        ("baseline", p_base, False, None, "calendar_pde_refined_baseline.csv"),
        ("pcv10_gsk_rollout_to_stable", p_base, True, "PCV10-GSK", "calendar_pde_refined_pcv10_gsk_rollout_to_stable.csv"),
        ("pcv13_rollout_to_stable", p_base, True, "PCV13", "calendar_pde_refined_pcv13_rollout_to_stable.csv"),
        ("pcv10sii_rollout_to_stable", p_base, True, "PCV10-SII", "calendar_pde_refined_pcv10sii_rollout_to_stable.csv"),
    ]
    return specs


def endpoint_direct_ve(p: Dict, vaccine_type: Optional[str], endpoint: str) -> float:
    if vaccine_type is None:
        return 0.0
    ve_map = p.get("vaccine_direct_effectiveness", {})
    entry = ve_map.get(vaccine_type, {})
    if endpoint in entry:
        return float(entry[endpoint])
    bridge = p.get("vaccine_direct_effectiveness_bridge", {})
    parent = bridge.get(vaccine_type)
    if parent is not None:
        parent_entry = ve_map.get(parent, {})
        if endpoint in parent_entry:
            return float(parent_entry[endpoint])
    return 0.0


def has_direct_protection(age_month: int, p: Dict) -> bool:
    start_month = int(p.get("direct_protection_start_month", 4))
    return int(age_month) >= start_month


def get_risk_multipliers(
    p: Dict,
    vaccine_type: Optional[str],
    endpoint: str,
    year_idx: int,
    protected_share_alive: float,
) -> Tuple[float, float, float, float, float]:
    if vaccine_type is None:
        return 1.0, 0.0, 1.0, 1.0, 1.0

    series = p.get("vaccine_irrs", {}).get(vaccine_type)
    irr = float(base.week4._get_population_irr_for_year(series, year_idx))
    ve = endpoint_direct_ve(p, vaccine_type, endpoint)
    q = max(0.0, min(1.0, float(protected_share_alive)))
    denom = 1.0 - q * ve

    if denom <= 1e-12:
        herd = 1.0
    else:
        herd = max(0.0, min(1.0, irr / denom))

    mult_protected = herd * (1.0 - ve)
    mult_unprotected = herd
    return irr, ve, herd, mult_protected, mult_unprotected


def protected_share_alive(A_prog: np.ndarray, A_non: np.ndarray, p: Dict) -> float:
    protected_mask = np.array([1.0 if has_direct_protection(age, p) else 0.0 for age in range(len(A_prog))], dtype=float)
    protected = float((A_prog * protected_mask).sum())
    total = float(A_prog.sum() + A_non.sum())
    return protected / total if total > 0.0 else 0.0


def build_refined_month_coefficients(
    p: Dict,
    calendar_month: int,
    A_prog_start: np.ndarray,
    A_non_start: np.ndarray,
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
        q = protected_share_alive(A_prog_start, A_non_start, p)
        _, _, _, mult_ipd_protected, mult_ipd_unprotected = get_risk_multipliers(p, vaccine_type, "ipd", year_idx, q)
        _, _, _, mult_pneu_protected, mult_pneu_unprotected = get_risk_multipliers(p, vaccine_type, "pneumonia", year_idx, q)
        _, _, _, mult_aom_protected, mult_aom_unprotected = get_risk_multipliers(p, vaccine_type, "aom", year_idx, q)
    else:
        mult_ipd_protected = mult_ipd_unprotected = 1.0
        mult_pneu_protected = mult_pneu_unprotected = 1.0
        mult_aom_protected = mult_aom_unprotected = 1.0

    coeffs = {name: np.zeros(n_age, dtype=float) for name in [
        "lam_men_p", "lam_pneu_p", "lam_sep_p", "lam_aom_p", "mu_dis_p", "mu_all_p",
        "lam_men_n", "lam_pneu_n", "lam_sep_n", "lam_aom_n", "mu_dis_n", "mu_all_n",
    ]}

    for age_month in range(n_age):
        base_prob_men, base_prob_pneu, base_prob_sep, base_prob_aom = base._disease_probabilities_for_age(p, age_month)
        direct_on = intervention and vaccine_type is not None and has_direct_protection(age_month, p)

        if direct_on:
            men_p_mult = mult_ipd_protected
            pneu_p_mult = mult_pneu_protected
            sep_p_mult = mult_ipd_protected
            aom_p_mult = mult_aom_protected
        elif intervention and vaccine_type is not None:
            men_p_mult = mult_ipd_unprotected
            pneu_p_mult = mult_pneu_unprotected
            sep_p_mult = mult_ipd_unprotected
            aom_p_mult = mult_aom_unprotected
        else:
            men_p_mult = pneu_p_mult = sep_p_mult = aom_p_mult = 1.0

        men_n_mult = mult_ipd_unprotected
        pneu_n_mult = mult_pneu_unprotected
        sep_n_mult = mult_ipd_unprotected
        aom_n_mult = mult_aom_unprotected

        p_men_p = base._clip01(base_prob_men * men_p_mult)
        p_pneu_p = base._clip01(base_prob_pneu * pneu_p_mult)
        p_sep_p = base._clip01(base_prob_sep * sep_p_mult)
        p_aom_p = base._clip01(base_prob_aom * aom_p_mult)

        p_men_n = base._clip01(base_prob_men * men_n_mult)
        p_pneu_n = base._clip01(base_prob_pneu * pneu_n_mult)
        p_sep_n = base._clip01(base_prob_sep * sep_n_mult)
        p_aom_n = base._clip01(base_prob_aom * aom_n_mult)

        coeffs["lam_men_p"][age_month] = base._prob_to_hazard(p_men_p)
        coeffs["lam_pneu_p"][age_month] = base._prob_to_hazard(p_pneu_p)
        coeffs["lam_sep_p"][age_month] = base._prob_to_hazard(p_sep_p)
        coeffs["lam_aom_p"][age_month] = base._prob_to_hazard(p_aom_p)
        coeffs["lam_men_n"][age_month] = base._prob_to_hazard(p_men_n)
        coeffs["lam_pneu_n"][age_month] = base._prob_to_hazard(p_pneu_n)
        coeffs["lam_sep_n"][age_month] = base._prob_to_hazard(p_sep_n)
        coeffs["lam_aom_n"][age_month] = base._prob_to_hazard(p_aom_n)

        d1_p = base._clip01(p_men_p * cfr_men)
        d2_p = base._clip01(p_pneu_p * prop_sev * cfr_pneu)
        d3_p = base._clip01(p_sep_p * cfr_sep)
        d1_n = base._clip01(p_men_n * cfr_men)
        d2_n = base._clip01(p_pneu_n * prop_sev * cfr_pneu)
        d3_n = base._clip01(p_sep_n * cfr_sep)

        coeffs["mu_dis_p"][age_month] = base._prob_to_hazard(d1_p) + base._prob_to_hazard(d2_p) + base._prob_to_hazard(d3_p)
        coeffs["mu_dis_n"][age_month] = base._prob_to_hazard(d1_n) + base._prob_to_hazard(d2_n) + base._prob_to_hazard(d3_n)
        coeffs["mu_all_p"][age_month] = lam_bg + coeffs["mu_dis_p"][age_month]
        coeffs["mu_all_n"][age_month] = lam_bg + coeffs["mu_dis_n"][age_month]

    coeffs["prop_sev"] = np.array([prop_sev], dtype=float)
    return coeffs


def pde_refined_reaction_rhs(_tau: float, y: np.ndarray, coeffs: Dict[str, np.ndarray]) -> np.ndarray:
    n_age = coeffs["mu_all_p"].size
    A_prog = np.maximum(y[:n_age], 0.0)
    A_non = np.maximum(y[n_age:2 * n_age], 0.0)

    dy = np.zeros_like(y, dtype=float)
    dy[:n_age] = -coeffs["mu_all_p"] * A_prog
    dy[n_age:2 * n_age] = -coeffs["mu_all_n"] * A_non

    prop_sev = float(coeffs["prop_sev"][0])
    inc_men = coeffs["lam_men_p"] * A_prog + coeffs["lam_men_n"] * A_non
    inc_pneu = coeffs["lam_pneu_p"] * A_prog + coeffs["lam_pneu_n"] * A_non
    inc_sep = coeffs["lam_sep_p"] * A_prog + coeffs["lam_sep_n"] * A_non
    inc_aom = coeffs["lam_aom_p"] * A_prog + coeffs["lam_aom_n"] * A_non
    inc_death = coeffs["mu_dis_p"] * A_prog + coeffs["mu_dis_n"] * A_non

    offset = 2 * n_age
    dy[offset + 0] = float(inc_men.sum())
    dy[offset + 1] = float((prop_sev * inc_pneu).sum())
    dy[offset + 2] = float(((1.0 - prop_sev) * inc_pneu).sum())
    dy[offset + 3] = float(inc_sep.sum())
    dy[offset + 4] = float(inc_aom.sum())
    dy[offset + 5] = float(inc_death.sum())
    return dy


def solve_one_month_pde_refined(
    A_prog: np.ndarray,
    A_non: np.ndarray,
    S: Dict[str, float],
    calendar_month: int,
    p: Dict,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
    rtol: float = 1e-10,
    atol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], Dict[str, float], float]:
    coeffs = build_refined_month_coefficients(p, calendar_month, A_prog, A_non, intervention, vaccine_type)
    n_age = len(A_prog)
    y0 = np.concatenate([A_prog.astype(float), A_non.astype(float), np.zeros(len(FLOW_NAMES), dtype=float)])
    sol = solve_ivp(
        fun=lambda tau, yy: pde_refined_reaction_rhs(tau, yy, coeffs),
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
    A_prog_after = y1[:n_age]
    A_non_after = y1[n_age:2 * n_age]
    flow_values = y1[2 * n_age:2 * n_age + len(FLOW_NAMES)]
    flows = {name: float(value) for name, value in zip(FLOW_NAMES, flow_values)}

    alive_before = float(A_prog.sum() + A_non.sum())
    alive_after = float(A_prog_after.sum() + A_non_after.sum())
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
    return A_prog_after, A_non_after, S_after, flows, non_neonatal_all_cause_deaths


def add_record_births_refined(
    A_prog: np.ndarray,
    A_non: np.ndarray,
    gross_births: float,
    p: Dict,
    intervention: bool,
    vaccine_type: Optional[str],
) -> Tuple[float, float]:
    alive_births = float(gross_births) * base.neonatal_survival_multiplier(p)
    neonatal_deaths = float(gross_births) - alive_births
    if intervention and vaccine_type is not None:
        cov = float(p.get("vaccine_coverage", 0.0))
        A_prog[0] += alive_births * cov
        A_non[0] += alive_births * (1.0 - cov)
    else:
        A_non[0] += alive_births
    return alive_births, neonatal_deaths


def build_refined_initial_population(
    p: Dict,
    shares: np.ndarray,
    scenario_name: str,
    retain_preexisting_sequelae: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], Dict[str, float]]:
    n_age = int(p.get("age_bins", 60))
    warmup_months = int(p.get("warmup_months", 60))
    births_annual = float(p["birth_cohort"])
    surv0 = base.neonatal_survival_multiplier(p)

    A_prog = np.zeros(n_age, dtype=float)
    A_non = np.zeros(n_age, dtype=float)
    S_warm = base.empty_sequela_stock()

    for warm_month in range(warmup_months):
        gross_births = births_annual * float(shares[warm_month % 12])
        A_non[0] += gross_births * surv0
        A_prog_after, A_non_after, S_warm, _, _ = solve_one_month_pde_refined(
            A_prog, A_non, S_warm, warm_month, p, intervention=False, vaccine_type=None
        )
        A_prog, A_non, _ = base.age_shift_for_next_month(A_prog_after, A_non_after)

    S_initial = base.copy_sequela_stock(S_warm) if retain_preexisting_sequelae else base.empty_sequela_stock()
    warm_meta = {
        "scenario_name": scenario_name,
        "warmup_months": warmup_months,
        "warmup_mode": "baseline_only_no_vaccination",
        "warmup_program_alive": float(A_prog.sum()),
        "warmup_nonprogram_alive": float(A_non.sum()),
        "warmup_alive_after_age_shift": float(A_prog.sum() + A_non.sum()),
        "warmup_preexisting_sequelae_total": float(sum(S_warm.values())),
        "retain_preexisting_sequelae": bool(retain_preexisting_sequelae),
    }
    return A_prog, A_non, S_initial, warm_meta


def run_pde_calendar_refined_monthly(
    p: Dict,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
    scenario_name: str = "scenario",
    retain_preexisting_sequelae: bool = False,
) -> pd.DataFrame:
    p = copy.deepcopy(p)
    record_months = int(p.get("record_months", 60))
    shares = base.normalize_birth_shares(p.get("birth_month_shares", [1.0 / 12.0] * 12))
    births_annual = float(p["birth_cohort"])

    A_prog, A_non, S, _ = build_refined_initial_population(
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
        _, neonatal_deaths = add_record_births_refined(A_prog, A_non, gross_births, p, intervention, vaccine_type)
        cum_all_cause_deaths += neonatal_deaths

        A_prog_after, A_non_after, S, flows, non_neonatal_deaths = solve_one_month_pde_refined(
            A_prog, A_non, S, calendar_month, p, intervention=intervention, vaccine_type=vaccine_type
        )
        cum_all_cause_deaths += non_neonatal_deaths

        stock_alive = float(A_prog_after.sum() + A_non_after.sum())
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

        A_prog, A_non, _ = base.age_shift_for_next_month(A_prog_after, A_non_after)

    return pd.DataFrame(history)


def run_all_scenarios_and_save_csv(
    outdir: Optional[str] = None,
    retain_preexisting_sequelae: bool = False,
    pcv10_sii_assumption: str = DEFAULT_PCV10_SII_ASSUMPTION,
    record_months: int = DEFAULT_RECORD_MONTHS,
) -> Dict[str, pd.DataFrame]:
    outdir = _resolve_outdir(outdir)
    os.makedirs(outdir, exist_ok=True)
    outputs: Dict[str, pd.DataFrame] = {}

    for scenario_name, p, intervention, vaccine_type, filename in _scenario_specs(pcv10_sii_assumption, record_months):
        monthly = run_pde_calendar_refined_monthly(
            p,
            intervention=intervention,
            vaccine_type=vaccine_type,
            scenario_name=scenario_name,
            retain_preexisting_sequelae=retain_preexisting_sequelae,
        )
        yearly = base.aggregate_yearly_like_ode_new(monthly)
        outputs[scenario_name] = yearly
        monthly.to_csv(os.path.join(outdir, filename.replace(".csv", "_monthly.csv")), index=False)
        yearly.to_csv(os.path.join(outdir, filename), index=False)

    period_summary = base.summarize_transition_and_stable(outputs)
    stable_annual = base.summarize_stable_annual(outputs)
    diagnostics = base.stability_diagnostics(outputs)
    period_summary.to_csv(os.path.join(outdir, "transition_and_stable_summary.csv"), index=False)
    stable_annual.to_csv(os.path.join(outdir, "stable_annual_summary.csv"), index=False)
    diagnostics.to_csv(os.path.join(outdir, "stability_diagnostics.csv"), index=False)
    return outputs


def compare_refined_to_original(
    refined_outputs: Optional[Dict[str, pd.DataFrame]] = None,
    outdir: Optional[str] = None,
    pcv10_sii_assumption: str = DEFAULT_PCV10_SII_ASSUMPTION,
    record_months: int = DEFAULT_RECORD_MONTHS,
) -> pd.DataFrame:
    outdir = _resolve_outdir(outdir)
    if refined_outputs is None:
        refined_outputs = run_all_scenarios_and_save_csv(
            outdir=outdir,
            pcv10_sii_assumption=pcv10_sii_assumption,
            record_months=record_months,
        )
    original_outputs = base.run_all_scenarios_and_save_csv(
        outdir=os.path.join(outdir, "original_reference_outputs"),
        retain_preexisting_sequelae=False,
        record_months=record_months,
    )

    refined_summary = base.summarize_transition_and_stable(refined_outputs).set_index(["Scenario", "Window"])
    original_summary = base.summarize_transition_and_stable(original_outputs).set_index(["Scenario", "Window"])
    rows = []
    for scenario, window in refined_summary.index:
        row = {"Scenario": scenario, "Window": window}
        for col in refined_summary.columns:
            refined_value = refined_summary.loc[(scenario, window), col]
            original_value = original_summary.loc[(scenario, window), col]
            if pd.api.types.is_number(refined_value):
                row[f"Refined_{col}"] = float(refined_value)
                row[f"Original_{col}"] = float(original_value)
                row[f"Delta_{col}"] = float(refined_value - original_value)
        rows.append(row)

    comparison = pd.DataFrame(rows)
    os.makedirs(outdir, exist_ok=True)
    comparison.to_csv(os.path.join(outdir, "calendar_pde_refined_vs_original_summary.csv"), index=False)
    return comparison


def run_pcv10_sii_assumption_sensitivity(
    outdir: Optional[str] = None,
    record_months: int = DEFAULT_RECORD_MONTHS,
) -> pd.DataFrame:
    outdir = _resolve_outdir(outdir)
    rows = []
    for assumption, multiplier in PCV10_SII_RELATIVE_VE_ASSUMPTIONS.items():
        assumption_dir = os.path.join(outdir, f"pcv10sii_{assumption}")
        outputs = run_all_scenarios_and_save_csv(
            outdir=assumption_dir,
            pcv10_sii_assumption=assumption,
            record_months=record_months,
        )
        summary = base.summarize_stable_annual(outputs)
        target = summary.loc[summary["Scenario"] == "pcv10sii_rollout_to_stable"].iloc[0].to_dict()
        rows.append({
            "PCV10_SII_Assumption": assumption,
            "Relative_VE_Multiplier": float(multiplier),
            **target,
        })
    out = pd.DataFrame(rows)
    os.makedirs(outdir, exist_ok=True)
    out.to_csv(os.path.join(outdir, "pcv10sii_assumption_sensitivity.csv"), index=False)
    return out


if __name__ == "__main__":
    outputs = run_all_scenarios_and_save_csv()
    period_summary = base.summarize_transition_and_stable(outputs)
    stable_annual = base.summarize_stable_annual(outputs)
    diagnostics = base.stability_diagnostics(outputs)
    comparison = compare_refined_to_original(outputs)
    sensitivity = run_pcv10_sii_assumption_sensitivity()
    print("\nTransition and stable-period summary, Week 5 PDE calendar refined model")
    print(period_summary.to_string(index=False))
    print("\nStable annual summary, Week 5 PDE calendar refined model")
    print(stable_annual.to_string(index=False))
    print("\nStability diagnostics, year 10 vs year 9")
    print(diagnostics[["Scenario", "MaxAbsPctChange_FlowsAndAlive"]].to_string(index=False))
    print("\nRefined minus original transition/stable summary")
    cols = ["Scenario", "Window"] + [c for c in comparison.columns if c.startswith("Delta_")]
    print(comparison[cols].to_string(index=False))
    print("\nPCV10-SII direct-VE assumption sensitivity")
    print(sensitivity.to_string(index=False))
