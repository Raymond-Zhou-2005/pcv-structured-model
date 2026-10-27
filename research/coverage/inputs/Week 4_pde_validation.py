"""
Cohort-comparable age-time PDE version of the PCV ODE model.

This file is designed for validation against the original ode_new.py.
It keeps the original assumptions, scenarios, output columns, yearly aggregation,
and numerical solver settings. The only structural rewrite is that the cohort is
interpreted as one characteristic of an age-time transport PDE:

    partial_t A(a,t) + partial_a A(a,t) = -mu(a,t) A(a,t).

For direct comparability with ode_new.py, this validation version uses:
    1. no warm-up;
    2. no continuing births after t=0;
    3. one birth cohort initialized at age 0;
    4. age advances by one month after each monthly interval;
    5. outputs recorded before age-out at the end of month 59.

Under these settings, the age-time PDE characteristic is algebraically the same
as the original closed-cohort ODE, so the rounded yearly CSV outputs should match
ode_new.py exactly.
"""

import copy
import os
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


def make_pcv10_gsk_series_dynamic():
    return [0.59, 0.49, 0.39, 0.32, 0.28, 0.26]


def make_pcv13_series_dynamic():
    a = 0.32
    b = 0.2193
    step = (a - b) / 3.0
    p1 = a - step
    p2 = a - 2.0 * step
    return [0.54, 0.38, 0.32, p1, p2, 0.2193]


def make_constant_series(val, desired_length=6):
    return [val] * desired_length


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
    'age_bins': 60,
}


def get_age_weight(month, disease_type):
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

    if 0 <= month <= 11:
        return 0.50 / 12.0
    if 12 <= month <= 23:
        return 0.25 / 12.0
    if 24 <= month <= 35:
        return 0.15 / 12.0
    return 0.10 / 24.0


def monthly_prob_from_incidence_per100k(inc_per100k, w_month, years_window=5.0):
    x = (float(inc_per100k) * float(years_window) * float(w_month)) / 100000.0
    return 1.0 - np.exp(-x)


def _get_population_irr_for_year(series, year_idx):
    if series is None:
        return 1.0

    n = len(series)
    if year_idx >= n:
        for j in range(n - 1, -1, -1):
            if series[j] is not None:
                return series[j]
        return 1.0

    if series[year_idx] is not None:
        return series[year_idx]

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
        prev_val = series[prev_idx]
        next_val = series[next_idx]
        span = next_idx - prev_idx
        step = (next_val - prev_val) / float(span)
        return prev_val + step * (year_idx - prev_idx)

    if prev_idx is not None:
        return series[prev_idx]
    if next_idx is not None:
        return series[next_idx]
    return 1.0


def infer_ve_from_irr_series(vaccine_irrs, vaccine_type, year_idx, coverage):
    series = vaccine_irrs.get(vaccine_type, None)
    irr = float(_get_population_irr_for_year(series, year_idx))
    cov = float(coverage)
    if cov > 0.0:
        ve = (1.0 - irr) / cov
    else:
        ve = 0.0
    ve = max(0.0, min(1.0, ve))
    return irr, ve


def get_background_prob_u5_only(p):
    p5_total = float(p['u5_mortality_rate']) / 1000.0
    if p5_total > 0.0:
        lambda_month_bg = -np.log(1.0 - p5_total) / 60.0
        prob_bg_month = 1.0 - np.exp(-lambda_month_bg)
    else:
        prob_bg_month = 0.0
    return prob_bg_month


def _clip01(x):
    return max(0.0, min(1.0, float(x)))


def _prob_to_hazard(prob):
    p = _clip01(prob)
    p = min(p, 1.0 - 1e-12)
    return -np.log(1.0 - p)


def _safe_nonneg(y):
    y = np.asarray(y, dtype=float).copy()
    y[y < 0.0] = 0.0
    return y


IDX_A_V = 0
IDX_A_U = 1
IDX_C_MEN = 2
IDX_C_PNEU_SEV = 3
IDX_C_PNEU_NONSEV = 4
IDX_C_SEP = 5
IDX_C_AOM = 6
IDX_C_DEATH = 7
IDX_S_COG = 8
IDX_S_HEAR = 9
IDX_S_EPI = 10
IDX_S_HEMI = 11
IDX_S_IMPLANT = 12
IDX_S_SPED = 13
N_STATE = 14


def rhs_age_time_characteristic(t, y, p, intervention=False, vaccine_type=None):
    """
    RHS on one characteristic of the age-time PDE.

    In this validation setting, age equals floor(t), so this RHS is the exact
    characteristic representation of the original closed-cohort ODE.
    """
    y = _safe_nonneg(y)

    A_v = y[IDX_A_V]
    A_u = y[IDX_A_U]
    S_cog = y[IDX_S_COG]
    S_hear = y[IDX_S_HEAR]
    S_epi = y[IDX_S_EPI]
    S_hemi = y[IDX_S_HEMI]
    S_implant = y[IDX_S_IMPLANT]
    S_sped = y[IDX_S_SPED]

    age_month = int(np.floor(float(t)))
    age_month = max(0, min(59, age_month))

    prob_bg_month = get_background_prob_u5_only(p)
    lam_bg = _prob_to_hazard(prob_bg_month)

    w_men = get_age_weight(age_month, 'Meningitis')
    w_pneu = get_age_weight(age_month, 'Pneumonia')
    w_aom = get_age_weight(age_month, 'AOM')

    base_prob_men = monthly_prob_from_incidence_per100k(p['inc_meningitis'], w_men)
    base_prob_pneu = monthly_prob_from_incidence_per100k(p['inc_pneumonia'], w_pneu)
    base_prob_sep = monthly_prob_from_incidence_per100k(p['inc_sepsis'], w_pneu)
    base_prob_aom = monthly_prob_from_incidence_per100k(p['inc_aom_all'], w_aom) * float(p['prop_aom_spn'])

    if intervention and (vaccine_type is not None):
        year_idx = age_month // 12
        cov = float(p.get('vaccine_coverage', 0.0))
        _, ve_year = infer_ve_from_irr_series(p['vaccine_irrs'], vaccine_type, year_idx, cov)
        mult_v = 1.0 - ve_year
    else:
        mult_v = 1.0

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

    prop_sev = float(p['prop_pneumo_severe'])
    cfr_men = float(p['cfr_meningitis'])
    cfr_pneu = float(p['cfr_pneumo_severe'])
    cfr_sep = float(p['cfr_sepsis'])

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

    dA_v = -mu_all_v * A_v
    dA_u = -mu_all_u * A_u

    inc_men = lam_men_v * A_v + lam_men_u * A_u
    inc_pneu_total = lam_pneu_v * A_v + lam_pneu_u * A_u
    inc_sep = lam_sep_v * A_v + lam_sep_u * A_u
    inc_aom = lam_aom_v * A_v + lam_aom_u * A_u

    inc_pneu_sev = prop_sev * inc_pneu_total
    inc_pneu_nonsev = (1.0 - prop_sev) * inc_pneu_total
    inc_death = mu_dis_v * A_v + mu_dis_u * A_u

    dC_men = inc_men
    dC_PSEV = inc_pneu_sev
    dC_PNON = inc_pneu_nonsev
    dC_sep = inc_sep
    dC_aom = inc_aom
    dC_death = inc_death

    A_total = A_v + A_u
    if A_total > 0.0:
        mu_all_avg = (mu_all_v * A_v + mu_all_u * A_u) / A_total
    else:
        mu_all_avg = 0.0

    men_survivor_flow = (1.0 - cfr_men) * inc_men

    src_cog = float(p['prob_cog']) * men_survivor_flow
    src_hear = float(p['prob_hear']) * men_survivor_flow
    src_epi = float(p['prob_epi']) * men_survivor_flow
    src_hemi = float(p['prob_hemi']) * men_survivor_flow
    src_implant = float(p.get('prob_cochlear_implant', 0.0)) * src_hear
    src_sped = float(p.get('prob_special_education', 0.0)) * src_hear

    dS_cog = src_cog - mu_all_avg * S_cog
    dS_hear = src_hear - mu_all_avg * S_hear
    dS_epi = src_epi - mu_all_avg * S_epi
    dS_hemi = src_hemi - mu_all_avg * S_hemi
    dS_implant = src_implant - mu_all_avg * S_implant
    dS_sped = src_sped - mu_all_avg * S_sped

    return np.array([
        dA_v,
        dA_u,
        dC_men,
        dC_PSEV,
        dC_PNON,
        dC_sep,
        dC_aom,
        dC_death,
        dS_cog,
        dS_hear,
        dS_epi,
        dS_hemi,
        dS_implant,
        dS_sped,
    ], dtype=float)


def run_age_time_pde_cohort_comparable(p, intervention=False, vaccine_type=None,
                                       method='RK45', rtol=1e-6, atol=1e-8, max_step=0.1):
    """Run the cohort-comparable age-time PDE characteristic model."""
    p = copy.deepcopy(p)
    birth = float(p['birth_cohort'])
    y = np.zeros(N_STATE, dtype=float)

    if intervention and (vaccine_type is not None):
        cov = float(p.get('vaccine_coverage', 0.0))
        y[IDX_A_V] = birth * cov
        y[IDX_A_U] = birth * (1.0 - cov)
    else:
        y[IDX_A_V] = 0.0
        y[IDX_A_U] = birth

    # Same as ode_new.py: neonatal mortality is applied once before month 0.
    p_neonatal = float(p.get('neonatal_mortality_rate', 0.0)) / 1000.0
    y[IDX_A_V] *= (1.0 - p_neonatal)
    y[IDX_A_U] *= (1.0 - p_neonatal)

    # This age vector is not needed to compute the characteristic, but it makes
    # the age-time interpretation explicit and inspectable.
    age_distribution = np.zeros(int(p.get('age_bins', 60)), dtype=float)
    age_distribution[0] = y[IDX_A_V] + y[IDX_A_U]

    history = []

    for month in range(60):
        y_start = _safe_nonneg(y)

        sol = solve_ivp(
            fun=lambda tt, yy: rhs_age_time_characteristic(tt, yy, p, intervention, vaccine_type),
            t_span=(float(month), float(month + 1)),
            y0=y_start,
            t_eval=[float(month + 1)],
            method=method,
            rtol=rtol,
            atol=atol,
            max_step=max_step,
        )

        if (not sol.success) or (sol.y is None) or (sol.y.shape[1] == 0):
            raise RuntimeError(f'solve_ivp failed at month {month}: {sol.message}')

        y_end = _safe_nonneg(sol.y[:, -1])

        n_new_men = y_end[IDX_C_MEN] - y_start[IDX_C_MEN]
        n_new_pneu_sev = y_end[IDX_C_PNEU_SEV] - y_start[IDX_C_PNEU_SEV]
        n_new_pneu_nonsev = y_end[IDX_C_PNEU_NONSEV] - y_start[IDX_C_PNEU_NONSEV]
        n_new_sepsis = y_end[IDX_C_SEP] - y_start[IDX_C_SEP]
        n_new_aom = y_end[IDX_C_AOM] - y_start[IDX_C_AOM]
        n_new_deaths = y_end[IDX_C_DEATH] - y_start[IDX_C_DEATH]

        alive_total_end = y_end[IDX_A_V] + y_end[IDX_A_U]

        step_data = {
            'Year': (month // 12) + 1,
            'New_Meningitis': n_new_men,
            'New_Pneumo_Sev': n_new_pneu_sev,
            'New_Pneumo_NonSev': n_new_pneu_nonsev,
            'New_Sepsis': n_new_sepsis,
            'New_AOM': n_new_aom,
            'New_Deaths_Pneumo': n_new_deaths,
            'Stock_Alive': alive_total_end,
            'Stock_Dead': birth - alive_total_end,
            'Stock_New_Meningitis': y_end[IDX_C_MEN],
            'Stock_New_Pneumo_Sev': y_end[IDX_C_PNEU_SEV],
            'Stock_New_Pneumo_NonSev': y_end[IDX_C_PNEU_NONSEV],
            'Stock_New_Sepsis': y_end[IDX_C_SEP],
            'Stock_New_AOM': y_end[IDX_C_AOM],
            'Stock_New_Deaths_Pneumo': y_end[IDX_C_DEATH],
            'Stock_Seq_Cognitive': y_end[IDX_S_COG],
            'Stock_Seq_Hearing': y_end[IDX_S_HEAR],
            'Stock_Seq_Epilepsy': y_end[IDX_S_EPI],
            'Stock_Seq_Hemiplegia': y_end[IDX_S_HEMI],
            'Stock_Seq_Hearing_Implant': y_end[IDX_S_IMPLANT],
            'Stock_Seq_Hearing_SpEd': y_end[IDX_S_SPED],
        }
        history.append(step_data)

        # Explicit age transport for transparency. For cohort-comparable output,
        # this transport is not allowed to remove the cohort from the year-5 stock
        # before the recorded output is written.
        if month < 59:
            age_distribution[:] = 0.0
            age_distribution[month + 1] = alive_total_end
        y = y_end

    return pd.DataFrame(history)


def aggregate_yearly_like_micro(df_monthly):
    stock_cols = [c for c in df_monthly.columns if 'Stock' in c]
    flow_cols = [c for c in df_monthly.columns if 'New' in c]

    final_df = df_monthly.groupby('Year')[flow_cols].sum()
    final_df[stock_cols] = df_monthly.groupby('Year')[stock_cols].last()

    col_order = [
        'New_Meningitis', 'New_Pneumo_Sev', 'New_Pneumo_NonSev', 'New_Sepsis', 'New_AOM', 'New_Deaths_Pneumo',
        'Stock_Alive', 'Stock_Dead',
        'Stock_New_Meningitis', 'Stock_New_Pneumo_Sev', 'Stock_New_Pneumo_NonSev', 'Stock_New_Sepsis', 'Stock_New_AOM', 'Stock_New_Deaths_Pneumo',
        'Stock_Seq_Cognitive', 'Stock_Seq_Hearing', 'Stock_Seq_Epilepsy', 'Stock_Seq_Hemiplegia',
        'Stock_Seq_Hearing_Implant', 'Stock_Seq_Hearing_SpEd',
    ]
    final_df = final_df[col_order]
    final_df['Year'] = final_df.index.astype(int)
    final_df = final_df.reset_index(drop=True)
    value_cols = [c for c in final_df.columns if c != 'Year']
    final_df[value_cols] = final_df[value_cols].astype(float)
    return final_df.round(0)


def _make_immediate_irrs_from_year5(vaccine_type):
    series = params_vietnam_micro['vaccine_irrs'][vaccine_type]
    if len(series) < 6 or series[5] is None:
        raise ValueError(f'Missing year5 IRR for {vaccine_type}')
    y5_irr = float(series[5])
    return {vaccine_type: [y5_irr] * 6}


def run_all_scenarios_and_save_csv(outdir='ode_age_time_pde_comparable_outputs_csv'):
    outputs = {}
    os.makedirs(outdir, exist_ok=True)

    df_base = aggregate_yearly_like_micro(
        run_age_time_pde_cohort_comparable(params_vietnam_micro, intervention=False, vaccine_type=None)
    )
    outputs['baseline'] = df_base
    df_base.to_csv(os.path.join(outdir, 'ode_age_time_pde_baseline.csv'), index=False)

    df_pcv10 = aggregate_yearly_like_micro(
        run_age_time_pde_cohort_comparable(params_vietnam_micro, intervention=True, vaccine_type='PCV10-GSK')
    )
    outputs['pcv10_dynamic'] = df_pcv10
    df_pcv10.to_csv(os.path.join(outdir, 'ode_age_time_pde_pcv10_dynamic.csv'), index=False)

    df_pcv13 = aggregate_yearly_like_micro(
        run_age_time_pde_cohort_comparable(params_vietnam_micro, intervention=True, vaccine_type='PCV13')
    )
    outputs['pcv13_dynamic'] = df_pcv13
    df_pcv13.to_csv(os.path.join(outdir, 'ode_age_time_pde_pcv13_dynamic.csv'), index=False)

    params_y5_pcv10 = copy.deepcopy(params_vietnam_micro)
    params_y5_pcv10['vaccine_irrs'] = _make_immediate_irrs_from_year5('PCV10-GSK')
    df_pcv10y5 = aggregate_yearly_like_micro(
        run_age_time_pde_cohort_comparable(params_y5_pcv10, intervention=True, vaccine_type='PCV10-GSK')
    )
    outputs['pcv10_year5_immediate'] = df_pcv10y5
    df_pcv10y5.to_csv(os.path.join(outdir, 'ode_age_time_pde_pcv10_year5_immediate.csv'), index=False)

    params_y5_pcv13 = copy.deepcopy(params_vietnam_micro)
    params_y5_pcv13['vaccine_irrs'] = _make_immediate_irrs_from_year5('PCV13')
    df_pcv13y5 = aggregate_yearly_like_micro(
        run_age_time_pde_cohort_comparable(params_y5_pcv13, intervention=True, vaccine_type='PCV13')
    )
    outputs['pcv13_year5_immediate'] = df_pcv13y5
    df_pcv13y5.to_csv(os.path.join(outdir, 'ode_age_time_pde_pcv13_year5_immediate.csv'), index=False)

    params_y5_pcv10sii = copy.deepcopy(params_vietnam_micro)
    params_y5_pcv10sii['vaccine_irrs'] = _make_immediate_irrs_from_year5('PCV10-SII')
    df_pcv10sii = aggregate_yearly_like_micro(
        run_age_time_pde_cohort_comparable(params_y5_pcv10sii, intervention=True, vaccine_type='PCV10-SII')
    )
    outputs['pcv10sii_year5_immediate'] = df_pcv10sii
    df_pcv10sii.to_csv(os.path.join(outdir, 'ode_age_time_pde_pcv10sii_year5_immediate.csv'), index=False)

    for name, df in outputs.items():
        print('\n=== SCENARIO:', name, '===\n')
        print(df.T)

    print('\nCSV files written to folder:', os.path.abspath(outdir))
    return outputs


def main():
    run_all_scenarios_and_save_csv()


if __name__ == '__main__':
    main()
