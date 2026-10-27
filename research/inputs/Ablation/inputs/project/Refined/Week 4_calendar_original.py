"""
Age-time calendar model, original vaccine logic.

Purpose
-------
This file is the next step after the cohort-comparable age-time PDE validation version.
It turns on the two calendar-time features that are intentionally absent from the
cohort-comparable validation run:

1. baseline warm-up to create a complete 0--59 month under-five age structure at recorded t=0;
2. continuing monthly births during the recorded calendar window.

To keep the comparison clean, this file keeps the same basic epidemiological and vaccine
assumptions as ode_new.py:

- same input parameters;
- same disease age weights;
- same incidence-to-monthly-probability mapping;
- same neonatal and background mortality construction;
- same annual IRR to one direct VE mapping;
- vaccinated/program births receive the direct-effect multiplier immediately, as in the
  original ODE logic;
- unvaccinated children keep baseline disease risk.

Interpretation
--------------
This is NOT a one-birth-cohort model. It represents calendar-time burden in the full
under-five population. The first recorded month contains an existing population generated
by warm-up, plus the first recorded month's newborns. Outputs use the same column names
as ode_new.py so that the files can pass into the same downstream summaries, but the
estimand is different from the old single-cohort ODE.
"""

import copy
import json
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ============================================================
# Parameter construction, aligned to ode_new.py
# ============================================================

def make_pcv10_gsk_series_dynamic() -> List[float]:
    return [0.59, 0.49, 0.39, 0.32, 0.28, 0.26]


def make_pcv13_series_dynamic() -> List[float]:
    a = 0.32
    b = 0.2193
    step = (a - b) / 3.0
    p1 = a - step
    p2 = a - 2.0 * step
    return [0.54, 0.38, 0.32, p1, p2, 0.2193]


def make_constant_series(val: float, desired_length: int = 6) -> List[float]:
    return [float(val)] * int(desired_length)


params_vietnam_micro = {
    'birth_cohort': 1384610,
    'u5_mortality_rate': 20.0,
    'neonatal_mortality_rate': 10.3,
    'inc_pneumonia': 1771,
    'inc_meningitis': 12,
    'inc_sepsis': 124,
    'inc_aom_all': 9570,
    'prop_pneumo_severe': 0.389,
    'prop_aom_spn': 0.303,
    'cfr_pneumo_severe': 0.014241,
    'cfr_meningitis': 0.032258,
    'cfr_sepsis': 0.08,
    'prob_cog': 0.016,
    'prob_epi': 0.022,
    'prob_hemi': 0.032,
    'prob_hear': 0.046,
    'prob_cochlear_implant': 0.025,
    'prob_special_education': 0.2574,
    'vaccine_coverage': 0.90,
    'vaccine_irrs': {
        'PCV10-GSK': make_pcv10_gsk_series_dynamic(),
        'PCV13': make_pcv13_series_dynamic(),
        'PCV10-SII': make_constant_series(0.2211, desired_length=6),
    },
    # New calendar-time controls. Uniform fallback is used when monthly birth shares are unknown.
    'birth_month_shares': [1.0 / 12.0] * 12,
    'warmup_months': 60,
    'record_months': 60,
    'age_bins': 60,
}


OUTPUT_COLUMNS = [
    'New_Meningitis', 'New_Pneumo_Sev', 'New_Pneumo_NonSev', 'New_Sepsis',
    'New_AOM', 'New_Deaths_Pneumo',
    'Stock_Alive', 'Stock_Dead',
    'Stock_New_Meningitis', 'Stock_New_Pneumo_Sev', 'Stock_New_Pneumo_NonSev',
    'Stock_New_Sepsis', 'Stock_New_AOM', 'Stock_New_Deaths_Pneumo',
    'Stock_Seq_Cognitive', 'Stock_Seq_Hearing', 'Stock_Seq_Epilepsy',
    'Stock_Seq_Hemiplegia', 'Stock_Seq_Hearing_Implant', 'Stock_Seq_Hearing_SpEd',
    'Year',
]

SEQ_KEYS = [
    'Stock_Seq_Cognitive',
    'Stock_Seq_Hearing',
    'Stock_Seq_Epilepsy',
    'Stock_Seq_Hemiplegia',
    'Stock_Seq_Hearing_Implant',
    'Stock_Seq_Hearing_SpEd',
]


# ============================================================
# Shared numerical helpers
# ============================================================

def normalize_birth_shares(shares: List[float]) -> np.ndarray:
    arr = np.asarray(list(shares), dtype=float)
    if arr.size != 12:
        raise ValueError('birth_month_shares must contain exactly 12 values')
    arr[arr < 0.0] = 0.0
    total = float(arr.sum())
    if total <= 0.0:
        return np.array([1.0 / 12.0] * 12, dtype=float)
    return arr / total


def get_age_weight(month: int, disease_type: str) -> float:
    if disease_type == 'Meningitis':
        if 0 <= month <= 11:
            return 0.65 / 12.0
        if 12 <= month <= 23:
            return 0.18 / 12.0
        if 24 <= month <= 35:
            return 0.07 / 12.0
        return 0.10 / 24.0

    if disease_type == 'AOM':
        if 0 <= month <= 11:
            return 0.1869 / 12.0
        if 12 <= month <= 23:
            return 0.1889 / 12.0
        if 24 <= month <= 59:
            return 0.6242 / 36.0
        return 0.0

    # Pneumonia and sepsis share the same age-weight pattern.
    if 0 <= month <= 11:
        return 0.50 / 12.0
    if 12 <= month <= 23:
        return 0.25 / 12.0
    if 24 <= month <= 35:
        return 0.15 / 12.0
    return 0.10 / 24.0


def monthly_prob_from_incidence_per100k(inc_per100k: float, w_month: float, years_window: float = 5.0) -> float:
    x = (float(inc_per100k) * float(years_window) * float(w_month)) / 100000.0
    return 1.0 - np.exp(-x)


def _get_population_irr_for_year(series: Optional[List[Optional[float]]], year_idx: int) -> float:
    if series is None:
        return 1.0
    n = len(series)
    if year_idx >= n:
        for j in range(n - 1, -1, -1):
            if series[j] is not None:
                return float(series[j])
        return 1.0
    if series[year_idx] is not None:
        return float(series[year_idx])

    prev_idx = None
    for j in range(year_idx - 1, -1, -1):
        if series[j] is not None:
            prev_idx = j
            break

    next_idx = None
    for j in range(year_idx + 1, n):
        if series[j] is not None:
            next_idx = j
            break

    if prev_idx is not None and next_idx is not None:
        prev_val = float(series[prev_idx])
        next_val = float(series[next_idx])
        span = next_idx - prev_idx
        return prev_val + (next_val - prev_val) * (year_idx - prev_idx) / float(span)
    if prev_idx is not None:
        return float(series[prev_idx])
    if next_idx is not None:
        return float(series[next_idx])
    return 1.0


def infer_ve_from_irr_series(vaccine_irrs: Dict[str, List[float]], vaccine_type: str, year_idx: int, coverage: float) -> Tuple[float, float]:
    series = vaccine_irrs.get(vaccine_type, None)
    irr = float(_get_population_irr_for_year(series, year_idx))
    cov = float(coverage)
    ve = (1.0 - irr) / cov if cov > 0.0 else 0.0
    ve = max(0.0, min(1.0, ve))
    return irr, ve


def _clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def _prob_to_hazard(prob: float) -> float:
    p = _clip01(prob)
    p = min(p, 1.0 - 1e-12)
    return -np.log(1.0 - p)


def get_background_prob_u5_only(p: Dict) -> float:
    p5_total = float(p['u5_mortality_rate']) / 1000.0
    if p5_total > 0.0:
        lambda_month_bg = -np.log(1.0 - p5_total) / 60.0
        return 1.0 - np.exp(-lambda_month_bg)
    return 0.0


def neonatal_survival_multiplier(p: Dict) -> float:
    return 1.0 - float(p.get('neonatal_mortality_rate', 0.0)) / 1000.0


def empty_sequela_stock() -> Dict[str, float]:
    return {k: 0.0 for k in SEQ_KEYS}


def copy_sequela_stock(S: Dict[str, float]) -> Dict[str, float]:
    return {k: float(S.get(k, 0.0)) for k in SEQ_KEYS}


# ============================================================
# One-month exact update over all age bins
# ============================================================

def one_month_exact_age_time_step(
    A_v: np.ndarray,
    A_u: np.ndarray,
    S: Dict[str, float],
    calendar_month: int,
    p: Dict,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], Dict[str, float], float]:
    """
    Advance disease and mortality for one calendar month across age bins.

    The update uses the same original direct-effect vaccine logic as ode_new.py.
    It is exact for a monthwise constant-hazard system.

    Returns
    -------
    A_v_after, A_u_after:
        Alive arrays after within-month mortality but before age shifting.
    S_after:
        Sequela stocks after within-month new sequelae and mortality adjustment.
    flows:
        Monthly incident disease and disease-death counts.
    non_neonatal_all_cause_deaths:
        Background plus disease mortality during this month, excluding neonatal deaths
        among new births. Neonatal deaths are tracked outside this function.
    """
    n_age = len(A_v)
    if len(A_u) != n_age:
        raise ValueError('A_v and A_u must have the same length')

    prob_bg_month = get_background_prob_u5_only(p)
    lam_bg = _prob_to_hazard(prob_bg_month)

    A_v_after = np.zeros_like(A_v, dtype=float)
    A_u_after = np.zeros_like(A_u, dtype=float)

    flows = {
        'New_Meningitis': 0.0,
        'New_Pneumo_Sev': 0.0,
        'New_Pneumo_NonSev': 0.0,
        'New_Sepsis': 0.0,
        'New_AOM': 0.0,
        'New_Deaths_Pneumo': 0.0,
    }

    year_idx = int(calendar_month) // 12
    if intervention and vaccine_type is not None:
        cov = float(p.get('vaccine_coverage', 0.0))
        _, ve_year = infer_ve_from_irr_series(p['vaccine_irrs'], vaccine_type, year_idx, cov)
        mult_v = 1.0 - ve_year
    else:
        mult_v = 1.0

    prop_sev = float(p['prop_pneumo_severe'])
    cfr_men = float(p['cfr_meningitis'])
    cfr_pneu = float(p['cfr_pneumo_severe'])
    cfr_sep = float(p['cfr_sepsis'])

    all_cause_deaths = 0.0
    men_survivor_flow = 0.0

    for age_month in range(n_age):
        w_men = get_age_weight(age_month, 'Meningitis')
        w_pneu = get_age_weight(age_month, 'Pneumonia')
        w_aom = get_age_weight(age_month, 'AOM')

        base_prob_men = monthly_prob_from_incidence_per100k(p['inc_meningitis'], w_men)
        base_prob_pneu = monthly_prob_from_incidence_per100k(p['inc_pneumonia'], w_pneu)
        base_prob_sep = monthly_prob_from_incidence_per100k(p['inc_sepsis'], w_pneu)
        base_prob_aom = monthly_prob_from_incidence_per100k(p['inc_aom_all'], w_aom) * float(p['prop_aom_spn'])

        p_men_v = _clip01(base_prob_men * mult_v)
        p_pneu_v = _clip01(base_prob_pneu * mult_v)
        p_sep_v = _clip01(base_prob_sep * mult_v)
        p_aom_v = _clip01(base_prob_aom * mult_v)

        p_men_u = _clip01(base_prob_men)
        p_pneu_u = _clip01(base_prob_pneu)
        p_sep_u = _clip01(base_prob_sep)
        p_aom_u = _clip01(base_prob_aom)

        lam_men_v = _prob_to_hazard(p_men_v)
        lam_pneu_v = _prob_to_hazard(p_pneu_v)
        lam_sep_v = _prob_to_hazard(p_sep_v)
        lam_aom_v = _prob_to_hazard(p_aom_v)

        lam_men_u = _prob_to_hazard(p_men_u)
        lam_pneu_u = _prob_to_hazard(p_pneu_u)
        lam_sep_u = _prob_to_hazard(p_sep_u)
        lam_aom_u = _prob_to_hazard(p_aom_u)

        d1_v = _clip01(p_men_v * cfr_men)
        d2_v = _clip01(p_pneu_v * prop_sev * cfr_pneu)
        d3_v = _clip01(p_sep_v * cfr_sep)
        d1_u = _clip01(p_men_u * cfr_men)
        d2_u = _clip01(p_pneu_u * prop_sev * cfr_pneu)
        d3_u = _clip01(p_sep_u * cfr_sep)

        mu_dis_v = _prob_to_hazard(d1_v) + _prob_to_hazard(d2_v) + _prob_to_hazard(d3_v)
        mu_dis_u = _prob_to_hazard(d1_u) + _prob_to_hazard(d2_u) + _prob_to_hazard(d3_u)
        mu_all_v = lam_bg + mu_dis_v
        mu_all_u = lam_bg + mu_dis_u

        alive_v = float(A_v[age_month])
        alive_u = float(A_u[age_month])

        surv_v = alive_v * np.exp(-mu_all_v)
        surv_u = alive_u * np.exp(-mu_all_u)

        # Integral of alive count over the month under constant all-cause hazard.
        kern_v = alive_v if mu_all_v <= 0.0 else alive_v * (1.0 - np.exp(-mu_all_v)) / mu_all_v
        kern_u = alive_u if mu_all_u <= 0.0 else alive_u * (1.0 - np.exp(-mu_all_u)) / mu_all_u

        inc_men_v = lam_men_v * kern_v
        inc_men_u = lam_men_u * kern_u
        inc_pneu_v = lam_pneu_v * kern_v
        inc_pneu_u = lam_pneu_u * kern_u
        inc_sep_v = lam_sep_v * kern_v
        inc_sep_u = lam_sep_u * kern_u
        inc_aom_v = lam_aom_v * kern_v
        inc_aom_u = lam_aom_u * kern_u
        inc_death_v = mu_dis_v * kern_v
        inc_death_u = mu_dis_u * kern_u

        flows['New_Meningitis'] += inc_men_v + inc_men_u
        flows['New_Pneumo_Sev'] += prop_sev * (inc_pneu_v + inc_pneu_u)
        flows['New_Pneumo_NonSev'] += (1.0 - prop_sev) * (inc_pneu_v + inc_pneu_u)
        flows['New_Sepsis'] += inc_sep_v + inc_sep_u
        flows['New_AOM'] += inc_aom_v + inc_aom_u
        flows['New_Deaths_Pneumo'] += inc_death_v + inc_death_u

        men_survivor_flow += (1.0 - cfr_men) * (inc_men_v + inc_men_u)
        all_cause_deaths += (alive_v - surv_v) + (alive_u - surv_u)

        A_v_after[age_month] = surv_v
        A_u_after[age_month] = surv_u

    alive_before = float(A_v.sum() + A_u.sum())
    alive_after = float(A_v_after.sum() + A_u_after.sum())
    surv_seq = max(0.0, min(1.0, alive_after / alive_before)) if alive_before > 0.0 else 1.0

    src_cog = float(p['prob_cog']) * men_survivor_flow
    src_hear = float(p['prob_hear']) * men_survivor_flow
    src_epi = float(p['prob_epi']) * men_survivor_flow
    src_hemi = float(p['prob_hemi']) * men_survivor_flow
    src_impl = float(p.get('prob_cochlear_implant', 0.0)) * src_hear
    src_sped = float(p.get('prob_special_education', 0.0)) * src_hear

    S_after = copy_sequela_stock(S)
    S_after['Stock_Seq_Cognitive'] = S_after['Stock_Seq_Cognitive'] * surv_seq + src_cog
    S_after['Stock_Seq_Hearing'] = S_after['Stock_Seq_Hearing'] * surv_seq + src_hear
    S_after['Stock_Seq_Epilepsy'] = S_after['Stock_Seq_Epilepsy'] * surv_seq + src_epi
    S_after['Stock_Seq_Hemiplegia'] = S_after['Stock_Seq_Hemiplegia'] * surv_seq + src_hemi
    S_after['Stock_Seq_Hearing_Implant'] = S_after['Stock_Seq_Hearing_Implant'] * surv_seq + src_impl
    S_after['Stock_Seq_Hearing_SpEd'] = S_after['Stock_Seq_Hearing_SpEd'] * surv_seq + src_sped

    return A_v_after, A_u_after, S_after, flows, all_cause_deaths


# ============================================================
# Calendar-time initialization and running
# ============================================================

@dataclass
class InitialStateSummary:
    scenario_name: str
    mode: str
    warmup_months: int
    retain_preexisting_sequelae: bool
    alive_before_first_record_births: float
    alive_after_first_record_births: float
    first_record_gross_births: float
    first_record_alive_births: float
    first_record_neonatal_deaths: float
    vaccinated_or_program_age0_after_births: float
    unvaccinated_age0_after_births: float
    preexisting_sequelae_total: float

    def to_dict(self) -> Dict[str, float]:
        return {
            'scenario_name': self.scenario_name,
            'mode': self.mode,
            'warmup_months': self.warmup_months,
            'retain_preexisting_sequelae': self.retain_preexisting_sequelae,
            'alive_before_first_record_births': self.alive_before_first_record_births,
            'alive_after_first_record_births': self.alive_after_first_record_births,
            'first_record_gross_births': self.first_record_gross_births,
            'first_record_alive_births': self.first_record_alive_births,
            'first_record_neonatal_deaths': self.first_record_neonatal_deaths,
            'vaccinated_or_program_age0_after_births': self.vaccinated_or_program_age0_after_births,
            'unvaccinated_age0_after_births': self.unvaccinated_age0_after_births,
            'preexisting_sequelae_total': self.preexisting_sequelae_total,
        }


def age_shift_for_next_month(A_v_after: np.ndarray, A_u_after: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Move survivors from age a to age a+1. Age 59 exits the under-five system."""
    n_age = len(A_v_after)
    A_v_next = np.zeros_like(A_v_after, dtype=float)
    A_u_next = np.zeros_like(A_u_after, dtype=float)
    age_out = float(A_v_after[-1] + A_u_after[-1])
    A_v_next[1:] = A_v_after[:-1]
    A_u_next[1:] = A_u_after[:-1]
    return A_v_next, A_u_next, age_out


def add_record_births(
    A_v: np.ndarray,
    A_u: np.ndarray,
    gross_births: float,
    p: Dict,
    intervention: bool,
    vaccine_type: Optional[str],
) -> Tuple[float, float]:
    """Add current month's newborns after neonatal mortality. Return alive births and neonatal deaths."""
    surv0 = neonatal_survival_multiplier(p)
    alive_births = float(gross_births) * surv0
    neonatal_deaths = float(gross_births) * (1.0 - surv0)

    if intervention and vaccine_type is not None:
        cov = float(p.get('vaccine_coverage', 0.0))
        A_v[0] += alive_births * cov
        A_u[0] += alive_births * (1.0 - cov)
    else:
        A_u[0] += alive_births
    return alive_births, neonatal_deaths


def build_calendar_initial_population(
    p: Dict,
    shares: np.ndarray,
    scenario_name: str,
    retain_preexisting_sequelae: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], Dict[str, float]]:
    """
    Create the existing under-five population before the first recorded month.

    Warm-up is always baseline-only in this original-logic calendar rollout version.
    Therefore, t=0 means vaccine introduction starts at the recorded window, and the
    existing children aged 1--59 months have not received a catch-up intervention.
    """
    n_age = int(p.get('age_bins', 60))
    warmup_months = int(p.get('warmup_months', 60))
    births_annual = float(p['birth_cohort'])
    surv0 = neonatal_survival_multiplier(p)

    A_v = np.zeros(n_age, dtype=float)
    A_u = np.zeros(n_age, dtype=float)
    S_warm = empty_sequela_stock()

    for warm_month in range(warmup_months):
        gross_births = births_annual * float(shares[warm_month % 12])
        A_u[0] += gross_births * surv0
        A_v_after, A_u_after, S_warm, _, _ = one_month_exact_age_time_step(
            A_v, A_u, S_warm, warm_month, p, intervention=False, vaccine_type=None
        )
        A_v, A_u, _ = age_shift_for_next_month(A_v_after, A_u_after)

    S_initial = copy_sequela_stock(S_warm) if retain_preexisting_sequelae else empty_sequela_stock()
    warm_meta = {
        'scenario_name': scenario_name,
        'warmup_months': warmup_months,
        'warmup_alive_after_age_shift': float(A_v.sum() + A_u.sum()),
        'warmup_vaccinated_or_program_alive': float(A_v.sum()),
        'warmup_unvaccinated_alive': float(A_u.sum()),
        'warmup_preexisting_sequelae_total': float(sum(S_warm.values())),
        'retain_preexisting_sequelae': bool(retain_preexisting_sequelae),
    }
    return A_v, A_u, S_initial, warm_meta


def run_age_time_calendar_original_monthly(
    p: Dict,
    intervention: bool = False,
    vaccine_type: Optional[str] = None,
    scenario_name: str = 'scenario',
    retain_preexisting_sequelae: bool = False,
    output_initial_state_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Run the full calendar-time under-five age-time model.

    Default interpretation:
    - baseline warm-up creates the existing under-five population at t=0;
    - recorded t=0 starts the observation window;
    - the first recorded month's births are added before the first disease step;
    - intervention scenarios vaccinate/program only new births from recorded t=0 onward;
    - cumulative counters start at zero at recorded t=0;
    - sequela stocks are reset at t=0 by default, so reported sequela stocks are generated
      during the recorded window only. Set retain_preexisting_sequelae=True if total prevalent
      under-five sequela stock is desired.
    """
    p = copy.deepcopy(p)
    n_age = int(p.get('age_bins', 60))
    record_months = int(p.get('record_months', 60))
    shares = normalize_birth_shares(p.get('birth_month_shares', [1.0 / 12.0] * 12))
    births_annual = float(p['birth_cohort'])

    A_v, A_u, S, warm_meta = build_calendar_initial_population(
        p=p,
        shares=shares,
        scenario_name=scenario_name,
        retain_preexisting_sequelae=retain_preexisting_sequelae,
    )

    cum_men = cum_psev = cum_pnon = cum_sep = cum_aom = cum_ddeath = 0.0
    cum_all_cause_deaths = 0.0
    history = []

    # Optional diagnostic tables for the teacher's t=0 question.
    initial_state_written = False

    for calendar_month in range(record_months):
        gross_births = births_annual * float(shares[calendar_month % 12])
        alive_births, neonatal_deaths = add_record_births(
            A_v, A_u, gross_births, p, intervention=intervention, vaccine_type=vaccine_type
        )
        cum_all_cause_deaths += neonatal_deaths

        if calendar_month == 0 and output_initial_state_dir is not None and not initial_state_written:
            os.makedirs(output_initial_state_dir, exist_ok=True)
            init_age_df = pd.DataFrame({
                'AgeMonth': np.arange(n_age, dtype=int),
                'Vaccinated_or_Program_Alive_at_t0_before_first_step': A_v.copy(),
                'Unvaccinated_Alive_at_t0_before_first_step': A_u.copy(),
                'Total_Alive_at_t0_before_first_step': A_v + A_u,
            })
            init_age_df.to_csv(os.path.join(output_initial_state_dir, f'initial_age_distribution_{scenario_name}.csv'), index=False)

            preexisting_total = float(warm_meta.get('warmup_preexisting_sequelae_total', 0.0))
            init_summary = InitialStateSummary(
                scenario_name=scenario_name,
                mode='calendar_rollout_baseline_warmup_new_births_only',
                warmup_months=int(p.get('warmup_months', 60)),
                retain_preexisting_sequelae=retain_preexisting_sequelae,
                alive_before_first_record_births=float(warm_meta['warmup_alive_after_age_shift']),
                alive_after_first_record_births=float(A_v.sum() + A_u.sum()),
                first_record_gross_births=float(gross_births),
                first_record_alive_births=float(alive_births),
                first_record_neonatal_deaths=float(neonatal_deaths),
                vaccinated_or_program_age0_after_births=float(A_v[0]),
                unvaccinated_age0_after_births=float(A_u[0]),
                preexisting_sequelae_total=preexisting_total,
            )
            pd.DataFrame([init_summary.to_dict()]).to_csv(
                os.path.join(output_initial_state_dir, f'initial_summary_{scenario_name}.csv'), index=False
            )
            with open(os.path.join(output_initial_state_dir, f'warmup_meta_{scenario_name}.json'), 'w', encoding='utf-8') as f:
                json.dump(warm_meta, f, indent=2)
            initial_state_written = True

        A_v_after, A_u_after, S, flows, non_neonatal_deaths = one_month_exact_age_time_step(
            A_v, A_u, S, calendar_month, p, intervention=intervention, vaccine_type=vaccine_type
        )
        cum_all_cause_deaths += non_neonatal_deaths

        # Record stocks after within-month mortality but before age shifting.
        # This gives the 0--59 month population that contributed to this calendar month.
        stock_alive = float(A_v_after.sum() + A_u_after.sum())

        cum_men += flows['New_Meningitis']
        cum_psev += flows['New_Pneumo_Sev']
        cum_pnon += flows['New_Pneumo_NonSev']
        cum_sep += flows['New_Sepsis']
        cum_aom += flows['New_AOM']
        cum_ddeath += flows['New_Deaths_Pneumo']

        history.append({
            'Month': calendar_month,
            'Year': calendar_month // 12 + 1,
            'New_Meningitis': flows['New_Meningitis'],
            'New_Pneumo_Sev': flows['New_Pneumo_Sev'],
            'New_Pneumo_NonSev': flows['New_Pneumo_NonSev'],
            'New_Sepsis': flows['New_Sepsis'],
            'New_AOM': flows['New_AOM'],
            'New_Deaths_Pneumo': flows['New_Deaths_Pneumo'],
            'Stock_Alive': stock_alive,
            'Stock_Dead': cum_all_cause_deaths,
            'Stock_New_Meningitis': cum_men,
            'Stock_New_Pneumo_Sev': cum_psev,
            'Stock_New_Pneumo_NonSev': cum_pnon,
            'Stock_New_Sepsis': cum_sep,
            'Stock_New_AOM': cum_aom,
            'Stock_New_Deaths_Pneumo': cum_ddeath,
            'Stock_Seq_Cognitive': S['Stock_Seq_Cognitive'],
            'Stock_Seq_Hearing': S['Stock_Seq_Hearing'],
            'Stock_Seq_Epilepsy': S['Stock_Seq_Epilepsy'],
            'Stock_Seq_Hemiplegia': S['Stock_Seq_Hemiplegia'],
            'Stock_Seq_Hearing_Implant': S['Stock_Seq_Hearing_Implant'],
            'Stock_Seq_Hearing_SpEd': S['Stock_Seq_Hearing_SpEd'],
        })

        A_v, A_u, _ = age_shift_for_next_month(A_v_after, A_u_after)

    return pd.DataFrame(history)


def aggregate_yearly_like_ode_new(df_monthly: pd.DataFrame) -> pd.DataFrame:
    stock_cols = [c for c in df_monthly.columns if 'Stock' in c]
    flow_cols = [c for c in df_monthly.columns if 'New' in c]

    yearly = df_monthly.groupby('Year')[flow_cols].sum()
    yearly[stock_cols] = df_monthly.groupby('Year')[stock_cols].last()
    yearly['Year'] = yearly.index.astype(int)
    yearly = yearly.reset_index(drop=True)
    yearly = yearly[OUTPUT_COLUMNS]

    value_cols = [c for c in yearly.columns if c != 'Year']
    yearly[value_cols] = yearly[value_cols].astype(float)
    return yearly.round(0)


def _make_immediate_irrs_from_year5(vaccine_type: str) -> Dict[str, List[float]]:
    series = params_vietnam_micro['vaccine_irrs'][vaccine_type]
    if len(series) < 6 or series[5] is None:
        raise ValueError(f'Missing year-5 IRR for {vaccine_type}')
    y5_irr = float(series[5])
    return {vaccine_type: [y5_irr] * 6}


def run_all_scenarios_and_save_csv(
    outdir: str = 'ode_age_time_calendar_original_outputs_csv',
    retain_preexisting_sequelae: bool = False,
) -> Dict[str, pd.DataFrame]:
    os.makedirs(outdir, exist_ok=True)
    init_dir = os.path.join(outdir, 'initial_state')
    outputs: Dict[str, pd.DataFrame] = {}

    scenario_specs = []
    scenario_specs.append(('baseline', params_vietnam_micro, False, None, 'ode_age_time_calendar_baseline.csv'))
    scenario_specs.append(('pcv10_dynamic', params_vietnam_micro, True, 'PCV10-GSK', 'ode_age_time_calendar_pcv10_dynamic.csv'))
    scenario_specs.append(('pcv13_dynamic', params_vietnam_micro, True, 'PCV13', 'ode_age_time_calendar_pcv13_dynamic.csv'))

    p_y5_pcv10 = copy.deepcopy(params_vietnam_micro)
    p_y5_pcv10['vaccine_irrs'] = _make_immediate_irrs_from_year5('PCV10-GSK')
    scenario_specs.append(('pcv10_year5_immediate', p_y5_pcv10, True, 'PCV10-GSK', 'ode_age_time_calendar_pcv10_year5_immediate.csv'))

    p_y5_pcv13 = copy.deepcopy(params_vietnam_micro)
    p_y5_pcv13['vaccine_irrs'] = _make_immediate_irrs_from_year5('PCV13')
    scenario_specs.append(('pcv13_year5_immediate', p_y5_pcv13, True, 'PCV13', 'ode_age_time_calendar_pcv13_year5_immediate.csv'))

    p_y5_pcv10sii = copy.deepcopy(params_vietnam_micro)
    p_y5_pcv10sii['vaccine_irrs'] = _make_immediate_irrs_from_year5('PCV10-SII')
    scenario_specs.append(('pcv10sii_year5_immediate', p_y5_pcv10sii, True, 'PCV10-SII', 'ode_age_time_calendar_pcv10sii_year5_immediate.csv'))

    for scenario_name, p, intervention, vaccine_type, filename in scenario_specs:
        monthly = run_age_time_calendar_original_monthly(
            p,
            intervention=intervention,
            vaccine_type=vaccine_type,
            scenario_name=scenario_name,
            retain_preexisting_sequelae=retain_preexisting_sequelae,
            output_initial_state_dir=init_dir,
        )
        yearly = aggregate_yearly_like_ode_new(monthly)
        outputs[scenario_name] = yearly
        monthly.to_csv(os.path.join(outdir, filename.replace('.csv', '_monthly.csv')), index=False)
        yearly.to_csv(os.path.join(outdir, filename), index=False)

    summary = summarize_five_year_totals(outputs)
    summary.to_csv(os.path.join(outdir, 'five_year_summary.csv'), index=False)
    return outputs


def summarize_five_year_totals(outputs: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for name, df in outputs.items():
        rows.append({
            'Scenario': name,
            'Meningitis_5y': float(df['New_Meningitis'].sum()),
            'Pneumo_Sev_5y': float(df['New_Pneumo_Sev'].sum()),
            'Pneumo_NonSev_5y': float(df['New_Pneumo_NonSev'].sum()),
            'Sepsis_5y': float(df['New_Sepsis'].sum()),
            'AOM_5y': float(df['New_AOM'].sum()),
            'Disease_Deaths_5y': float(df['New_Deaths_Pneumo'].sum()),
            'Alive_End_Y5': float(df.loc[df.index[-1], 'Stock_Alive']),
            'All_Cause_Deaths_Record_Window_End_Y5': float(df.loc[df.index[-1], 'Stock_Dead']),
        })
    return pd.DataFrame(rows)


if __name__ == '__main__':
    outputs = run_all_scenarios_and_save_csv(
        outdir='ode_age_time_calendar_original_outputs_csv',
        retain_preexisting_sequelae=False,
    )
    summary = summarize_five_year_totals(outputs)
    print('\nFive-year calendar-window summary, age-time calendar original model')
    print(summary.to_string(index=False))
    print('\nCSV outputs saved to: ode_age_time_calendar_original_outputs_csv')
