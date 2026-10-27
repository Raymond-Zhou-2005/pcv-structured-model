import copy
import os
import numpy as np
import pandas as pd
from tqdm import trange


class VietnamRandSimulation:
    def __init__(self, params):
        """Microsimulation with monthly cycles over 5 years."""
        self.p = params
        self.n_agents = int(params['birth_cohort'])
        self.scaling_factor = 1.0
        self.cycles = 60

    def get_monthly_prob_from_annual_rate(self, annual_rate):
        """Convert a monthly per-person rate into a monthly probability."""
        return 1.0 - np.exp(-(annual_rate / 100000.0))

    def get_age_weight(self, month, disease_type):
        """Return the monthly share of the 5-year burden for each disease."""
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

        # Pneumonia and Sepsis use the same weight pattern.
        if 0 <= month <= 11:
            return 0.50 / 12.0
        if 12 <= month <= 23:
            return 0.25 / 12.0
        if 24 <= month <= 35:
            return 0.15 / 12.0
        return 0.10 / 24.0

    def _get_population_irr_for_year(self, vaccine_type, year_idx):
        """Get the population IRR for a given year index, with interpolation."""
        irr_dict = self.p.get('vaccine_irrs', {})
        series = irr_dict.get(vaccine_type, None)
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

    def infer_ve_from_irr_series(self, vaccine_type, year_idx):
        """Infer VE from IRR and coverage for the given year."""
        irr = float(self._get_population_irr_for_year(vaccine_type, year_idx))
        cov = float(self.p.get('vaccine_coverage', 0.0))
        if cov > 0.0:
            ve = (1.0 - irr) / cov
        else:
            ve = 0.0
        ve = max(0.0, min(1.0, ve))
        return irr, ve

    def run(self, intervention=False, vaccine_type=None):
        """Run one microsimulation and return a monthly output DataFrame."""
        p = self.p

        alive = np.ones(self.n_agents, dtype=bool)

        vaccinated = np.zeros(self.n_agents, dtype=bool)
        if intervention and (vaccine_type is not None):
            cov = float(p.get('vaccine_coverage', 0.0))
            rnd_v = np.random.rand(self.n_agents)
            vaccinated = rnd_v < cov

        has_cognitive = np.zeros(self.n_agents, dtype=bool)
        has_epilepsy = np.zeros(self.n_agents, dtype=bool)
        has_hearing = np.zeros(self.n_agents, dtype=bool)
        has_hemiplegia = np.zeros(self.n_agents, dtype=bool)
        has_cochlear_implant = np.zeros(self.n_agents, dtype=bool)
        has_special_education = np.zeros(self.n_agents, dtype=bool)

        cum_new_men = 0
        cum_new_pneu_sev = 0
        cum_new_pneu_nonsev = 0
        cum_new_sepsis = 0
        cum_new_aom = 0
        cum_new_deaths_pneumo = 0

        # Neonatal mortality is applied once before month 0.
        p_neonatal = float(self.p.get('neonatal_mortality_rate', 0.0)) / 1000.0
        if p_neonatal > 0.0:
            rnd_neonatal = np.random.rand(self.n_agents)
            neonatal_deaths = alive & (rnd_neonatal < p_neonatal)
            alive[neonatal_deaths] = False

        # Background mortality uses the full under-five mortality rate.
        p5_total = float(self.p['u5_mortality_rate']) / 1000.0
        if p5_total > 0.0:
            lambda_month_bg = -np.log(1.0 - p5_total) / 60.0
            prob_bg_month = 1.0 - np.exp(-lambda_month_bg)
        else:
            prob_bg_month = 0.0

        history = []

        for month in trange(self.cycles, desc='Simulating months'):
            alive_before = alive.sum()

            # Background deaths occur at the start of each month.
            rnd_bg = np.random.rand(self.n_agents)
            bg_deaths = alive & (rnd_bg < prob_bg_month)
            alive[bg_deaths] = False

            alive_after_bg = alive.sum()

            # Month-specific disease risks.
            w_men = self.get_age_weight(month, 'Meningitis')
            w_pneu = self.get_age_weight(month, 'Pneumonia')
            w_aom = self.get_age_weight(month, 'AOM')

            p_men = 1.0 - np.exp(-(float(p['inc_meningitis']) * 5.0 * w_men / 100000.0))
            p_pneu = 1.0 - np.exp(-(float(p['inc_pneumonia']) * 5.0 * w_pneu / 100000.0))
            p_sepsis = 1.0 - np.exp(-(float(p['inc_sepsis']) * 5.0 * w_pneu / 100000.0))
            p_aom = (1.0 - np.exp(-(float(p['inc_aom_all']) * 5.0 * w_aom / 100000.0))) * float(p['prop_aom_spn'])

            if intervention and (vaccine_type is not None):
                year_idx = month // 12
                _, ve_year = self.infer_ve_from_irr_series(vaccine_type, year_idx)
                mult_v = 1.0 - ve_year
            else:
                mult_v = 1.0

            p_men_v = min(1.0, p_men * mult_v)
            p_pneu_v = min(1.0, p_pneu * mult_v)
            p_sepsis_v = min(1.0, p_sepsis * mult_v)
            p_aom_v = min(1.0, p_aom * mult_v)

            p_men_u = p_men
            p_pneu_u = p_pneu
            p_sepsis_u = p_sepsis
            p_aom_u = p_aom

            rnd_men = np.random.rand(self.n_agents)
            rnd_pneu = np.random.rand(self.n_agents)
            rnd_sepsis = np.random.rand(self.n_agents)
            rnd_aom = np.random.rand(self.n_agents)

            new_cases_men = alive & (
                (vaccinated & (rnd_men < p_men_v)) | (~vaccinated & (rnd_men < p_men_u))
            )

            new_cases_pneu_total = alive & (
                (vaccinated & (rnd_pneu < p_pneu_v)) | (~vaccinated & (rnd_pneu < p_pneu_u))
            )

            new_cases_sepsis = alive & (
                (vaccinated & (rnd_sepsis < p_sepsis_v)) | (~vaccinated & (rnd_sepsis < p_sepsis_u))
            )

            new_cases_aom = alive & (
                (vaccinated & (rnd_aom < p_aom_v)) | (~vaccinated & (rnd_aom < p_aom_u))
            )

            prop_sev = float(p['prop_pneumo_severe'])
            rnd_sev = np.random.rand(self.n_agents)
            pneu_sev = new_cases_pneu_total & (rnd_sev < prop_sev)
            pneu_nonsev = new_cases_pneu_total & ~pneu_sev

            n_new_men = new_cases_men.sum()
            n_new_pneu_sev = pneu_sev.sum()
            n_new_pneu_nonsev = pneu_nonsev.sum()
            n_new_sepsis = new_cases_sepsis.sum()
            n_new_aom = new_cases_aom.sum()

            cum_new_men += n_new_men
            cum_new_pneu_sev += n_new_pneu_sev
            cum_new_pneu_nonsev += n_new_pneu_nonsev
            cum_new_sepsis += n_new_sepsis
            cum_new_aom += n_new_aom

            cfr_men = float(p['cfr_meningitis'])
            cfr_pneu = float(p['cfr_pneumo_severe'])
            cfr_sepsis = float(p['cfr_sepsis'])

            rnd_death_men = np.random.rand(self.n_agents)
            rnd_death_pneu = np.random.rand(self.n_agents)
            rnd_death_sepsis = np.random.rand(self.n_agents)

            deaths_men = new_cases_men & (rnd_death_men < cfr_men)
            deaths_pneu = pneu_sev & (rnd_death_pneu < cfr_pneu)
            deaths_sepsis = new_cases_sepsis & (rnd_death_sepsis < cfr_sepsis)

            deaths_any = deaths_men | deaths_pneu | deaths_sepsis
            alive[deaths_any] = False

            n_new_deaths = deaths_any.sum()
            cum_new_deaths_pneumo += n_new_deaths

            surv_men = new_cases_men & ~deaths_men
            n_surv_men = surv_men.sum()

            if n_surv_men > 0:
                prob_cog = float(p['prob_cog'])
                prob_epi = float(p['prob_epi'])
                prob_hear = float(p['prob_hear'])
                prob_hemi = float(p['prob_hemi'])

                rnd_cog = np.random.rand(self.n_agents)
                rnd_epi = np.random.rand(self.n_agents)
                rnd_hear = np.random.rand(self.n_agents)
                rnd_hemi = np.random.rand(self.n_agents)

                cog_new = surv_men & (rnd_cog < prob_cog)
                epi_new = surv_men & (rnd_epi < prob_epi)
                hear_new = surv_men & (rnd_hear < prob_hear)
                hemi_new = surv_men & (rnd_hemi < prob_hemi)

                has_cognitive[cog_new] = True
                has_epilepsy[epi_new] = True
                has_hearing[hear_new] = True
                has_hemiplegia[hemi_new] = True

                prob_implant = float(p.get('prob_cochlear_implant', 0.0))
                prob_sped = float(p.get('prob_special_education', 0.0))

                rnd_implant = np.random.rand(self.n_agents)
                rnd_sped = np.random.rand(self.n_agents)

                implant_new = hear_new & (rnd_implant < prob_implant)
                sped_new = hear_new & (rnd_sped < prob_sped)

                has_cochlear_implant[implant_new] = True
                has_special_education[sped_new] = True

            alive_after = alive.sum()

            step_data = {
                'Year': (month // 12) + 1,
                'New_Meningitis': int(n_new_men),
                'New_Pneumo_Sev': int(n_new_pneu_sev),
                'New_Pneumo_NonSev': int(n_new_pneu_nonsev),
                'New_Sepsis': int(n_new_sepsis),
                'New_AOM': int(n_new_aom),
                'New_Deaths_Pneumo': int(n_new_deaths),
                'Stock_Alive': int(alive_after),
                'Stock_Dead': int(self.n_agents - alive_after),
                'Stock_New_Meningitis': int(cum_new_men),
                'Stock_New_Pneumo_Sev': int(cum_new_pneu_sev),
                'Stock_New_Pneumo_NonSev': int(cum_new_pneu_nonsev),
                'Stock_New_Sepsis': int(cum_new_sepsis),
                'Stock_New_AOM': int(cum_new_aom),
                'Stock_New_Deaths_Pneumo': int(cum_new_deaths_pneumo),
                'Stock_Seq_Cognitive': int((has_cognitive & alive).sum()),
                'Stock_Seq_Hearing': int((has_hearing & alive).sum()),
                'Stock_Seq_Epilepsy': int((has_epilepsy & alive).sum()),
                'Stock_Seq_Hemiplegia': int((has_hemiplegia & alive).sum()),
                'Stock_Seq_Hearing_Implant': int((has_cochlear_implant & alive).sum()),
                'Stock_Seq_Hearing_SpEd': int((has_special_education & alive).sum()),
            }
            history.append(step_data)

        df_monthly = pd.DataFrame(history)

        # Apply scaling to value columns only.
        value_cols = [c for c in df_monthly.columns if c != 'Year']
        df_monthly[value_cols] = df_monthly[value_cols] * float(self.scaling_factor)

        return df_monthly


def make_pcv10_gsk_series_dynamic():
    """PCV10-GSK dynamic IRR series."""
    return [0.59, 0.49, 0.39, 0.32, 0.28, 0.26]


def make_pcv13_series_dynamic():
    """PCV13 dynamic IRR series with two interpolated points."""
    a = 0.32
    b = 0.2193
    span = a - b
    t = span / 3.0
    p1 = a - t
    p2 = a - 2.0 * t
    return [0.54, 0.38, 0.32, p1, p2, 0.2193]


def make_constant_series(val, desired_length=6):
    """Return a constant IRR series."""
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
}


base_model = VietnamRandSimulation(params_vietnam_micro)


def aggregate_yearly_like_micro(df_monthly):
    """Aggregate monthly outputs to yearly outputs with ODE-matched columns."""
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


# Backward-compatible name.
summarize_micro = aggregate_yearly_like_micro


def _make_immediate_irrs_from_year5(vaccine_type):
    """Build a constant IRR series from the sixth value of the base series."""
    series = params_vietnam_micro['vaccine_irrs'][vaccine_type]
    if len(series) < 6 or series[5] is None:
        raise ValueError(f'Missing year5 IRR for {vaccine_type}')
    y5_irr = float(series[5])
    return {vaccine_type: [y5_irr] * 6}


def _get_year5_irr(vaccine_type, vaccine_irrs):
    """Return the year-5 IRR (6th element) for the given vaccine type."""
    series = vaccine_irrs.get(vaccine_type, None)
    if (series is not None) and (len(series) >= 6) and (series[5] is not None):
        return float(series[5])
    if series is not None:
        for j in range(len(series) - 1, -1, -1):
            if series[j] is not None:
                return float(series[j])
    return 1.0


def _run_all_scenarios_and_save_csv():
    """Run the six scenarios and save yearly CSV outputs."""
    outputs = {}
    outdir = 'markov_outputs_csv'
    os.makedirs(outdir, exist_ok=True)

    df_base_month = base_model.run(intervention=False, vaccine_type=None)
    df_base = aggregate_yearly_like_micro(df_base_month)
    outputs['baseline'] = df_base
    df_base.to_csv(os.path.join(outdir, 'markov_baseline.csv'), index=False)

    df_pcv10_month = base_model.run(intervention=True, vaccine_type='PCV10-GSK')
    df_pcv10 = aggregate_yearly_like_micro(df_pcv10_month)
    outputs['pcv10_dynamic'] = df_pcv10
    df_pcv10.to_csv(os.path.join(outdir, 'markov_pcv10_dynamic.csv'), index=False)

    df_pcv13_month = base_model.run(intervention=True, vaccine_type='PCV13')
    df_pcv13 = aggregate_yearly_like_micro(df_pcv13_month)
    outputs['pcv13_dynamic'] = df_pcv13
    df_pcv13.to_csv(os.path.join(outdir, 'markov_pcv13_dynamic.csv'), index=False)

    params_y5_pcv10 = copy.deepcopy(params_vietnam_micro)
    params_y5_pcv10['vaccine_irrs'] = _make_immediate_irrs_from_year5('PCV10-GSK')
    model_y5_pcv10 = VietnamRandSimulation(params_y5_pcv10)
    df_pcv10y5_month = model_y5_pcv10.run(intervention=True, vaccine_type='PCV10-GSK')
    df_pcv10y5 = aggregate_yearly_like_micro(df_pcv10y5_month)
    outputs['pcv10_year5_immediate'] = df_pcv10y5
    df_pcv10y5.to_csv(os.path.join(outdir, 'markov_pcv10_year5_immediate.csv'), index=False)

    params_y5_pcv13 = copy.deepcopy(params_vietnam_micro)
    params_y5_pcv13['vaccine_irrs'] = _make_immediate_irrs_from_year5('PCV13')
    model_y5_pcv13 = VietnamRandSimulation(params_y5_pcv13)
    df_pcv13y5_month = model_y5_pcv13.run(intervention=True, vaccine_type='PCV13')
    df_pcv13y5 = aggregate_yearly_like_micro(df_pcv13y5_month)
    outputs['pcv13_year5_immediate'] = df_pcv13y5
    df_pcv13y5.to_csv(os.path.join(outdir, 'markov_pcv13_year5_immediate.csv'), index=False)

    params_y5_pcv10sii = copy.deepcopy(params_vietnam_micro)
    params_y5_pcv10sii['vaccine_irrs'] = _make_immediate_irrs_from_year5('PCV10-SII')
    model_y5_pcv10sii = VietnamRandSimulation(params_y5_pcv10sii)
    df_pcv10sii_month = model_y5_pcv10sii.run(intervention=True, vaccine_type='PCV10-SII')
    df_pcv10sii = aggregate_yearly_like_micro(df_pcv10sii_month)
    outputs['pcv10sii_year5_immediate'] = df_pcv10sii
    df_pcv10sii.to_csv(os.path.join(outdir, 'markov_pcv10sii_year5_immediate.csv'), index=False)

    for name, df in outputs.items():
        print('\n=== SCENARIO:', name, '===\n')
        print(df.T)

    print('\nCSV files written to folder:', os.path.abspath(outdir))


def main():
    """Entry point."""
    _run_all_scenarios_and_save_csv()


if __name__ == '__main__':
    main()