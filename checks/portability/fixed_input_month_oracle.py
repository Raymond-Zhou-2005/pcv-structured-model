"""Diagnose monthly RK45 error on frozen saved inputs without changing the model.

Usage: python fixed_input_month_oracle.py --repo PATH --output PATH
If placed under checks/portability, the repository root is discovered upwards.
The oracle repeats the frozen-coefficient exposure formula in
research/src/experiments.py (simulate_refined); it is not an empirical
validation, a full nonlinear trajectory error bound, or a parameter PSA.
"""
from pathlib import Path
import argparse
import contextlib
from decimal import Decimal, localcontext
import hashlib
import io
import json
import platform
import sys

import numpy as np
import pandas as pd
import scipy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    repo = args.repo.resolve() if args.repo else next(
        p for p in Path(__file__).resolve().parents if (p / 'reproduce.py').exists())
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(repo / 'research/inputs/Ablation/src'))
    import run_ablation as r
    w6 = r.w6
    rows, flow_rows = [], []
    conversions, coefficient_snapshots = {}, {}
    source_paths = ['research/src/experiments.py',
                    'research/inputs/Ablation/src/run_ablation.py',
                    'research/inputs/Ablation/inputs/project/Refined/Week 6_calendar_pde_program_realism.py',
                    'research/inputs/Ablation/inputs/project/Refined/Week 7_qaly_refinement_impact.py']
    sources = {p: hashlib.sha256((repo / p).read_bytes()).hexdigest() for p in source_paths}
    scenarios = ['no_vaccination', 'pcv10_gsk_D1W1G1',
                 'pcv13_D1W1G1', 'pcv10sii_D1W1G1']
    for name in scenarios:
        directory = repo / 'research/inputs/Final/outputs/runs' / name
        p = json.loads((directory / 'parameters.json').read_text())
        e = json.loads((directory / 'economic_parameters.json').read_text())
        e['prod_discounted_age_contrib'] = {
            int(age): value for age, value in e['prod_discounted_age_contrib'].items()}
        with np.load(directory / 'states.npz') as archive:
            states = archive['alive_after_births_and_doses'].copy()
        monthly = pd.read_csv(directory / 'monthly.csv')
        herd = pd.read_csv(directory / 'herd.csv')
        vaccine = next((v for v, (prefix, _) in r.PRODUCTS.items()
                        if name.startswith(prefix + '_')), None)
        death_age = float(r.econ.CALENDAR_CONFIG['AVG_U5_DEATH_AGE_YEARS'])
        prod_weight = r.econ.productivity_loss_from_death_event(death_age, e)
        qaly_weight = r.econ.future_qaly_from_death_event(
            death_age, e['life_expectancy'], e['discount_rate'])
        for disease, incidence, weight_key in [
                ('meningitis', 'inc_meningitis', 'Meningitis'),
                ('pneumonia', 'inc_pneumonia', 'Pneumonia'),
                ('sepsis', 'inc_sepsis', 'Pneumonia'),
                ('aom_all', 'inc_aom_all', 'AOM')]:
            for age in [0, 12, 24, 36, 59]:
                w = w6.base.get_age_weight(age, weight_key)
                x = float(p[incidence]) * 5.0 * float(w) / 100000.0
                legacy = float(w6.base.monthly_prob_from_incidence_per100k(p[incidence], w))
                stable = float(-np.expm1(-x))
                with localcontext() as ctx:
                    ctx.prec = 80
                    reference = Decimal(1) - (-Decimal.from_float(x)).exp()
                    conversions['probability:' + x.hex()] = dict(
                        conversion='1-exp(-x)', example_disease=disease,
                        example_age=age, input=x, input_hex=x.hex(),
                        legacy=legacy, stable=stable, reference_decimal=str(reference),
                        legacy_signed_error=float(Decimal.from_float(legacy) - reference),
                        stable_signed_error=float(Decimal.from_float(stable) - reference))
        for month in [0, 24, 60, 119]:
            A = states[month].copy()
            S = w6.base.empty_sequela_stock()
            if month:
                for key in S:
                    S[key] = float(monthly.loc[month - 1, key])
            original_herd = w6.herd_multiplier_for_endpoint
            original_hazard = w6.base._prob_to_hazard
            def audited_hazard(prob):
                value = float(original_hazard(prob))
                clipped = min(w6.base._clip01(prob), 1 - 1e-12)
                key = 'hazard:' + clipped.hex()
                if key not in conversions:
                    with localcontext() as ctx:
                        ctx.prec = 80
                        reference = -(Decimal(1) - Decimal.from_float(clipped)).ln()
                        stable = float(-np.log1p(-clipped))
                        conversions[key] = dict(
                            conversion='-log(1-p)', example_disease='', example_age='',
                            input=clipped, input_hex=clipped.hex(), legacy=value,
                            stable=stable, reference_decimal=str(reference),
                            legacy_signed_error=float(Decimal.from_float(value) - reference),
                            stable_signed_error=float(Decimal.from_float(stable) - reference))
                return value
            w6.base._prob_to_hazard = audited_hazard
            w6.herd_multiplier_for_endpoint = lambda pp, vv, endpoint, year, aa: (
                float(herd.loc[month, 'H_' + endpoint]),
                float(herd.loc[month, 'H_target_' + endpoint]),
                float(herd.loc[month, 'relative_carriage_index']))
            try:
                c = w6.build_week6_month_coefficients(p, month, A, vaccine is not None, vaccine)
            finally:
                w6.herd_multiplier_for_endpoint = original_herd
                w6.base._prob_to_hazard = original_hazard
            for key, value in c.items():
                coefficient_snapshots[name + '__' + str(month) + '__' + key] = value.copy()
            mu = c['mu_all']
            factor = np.ones_like(mu)
            np.divide(-np.expm1(-mu), mu, out=factor, where=mu != 0)
            exposure = A * factor
            exact_A = A * np.exp(-mu)
            exact_flows = dict(
                New_Meningitis=float(np.sum(exposure * c['lam_men'])),
                New_Pneumo_Sev=float(np.sum(exposure * c['lam_pneu']) * p['prop_pneumo_severe']),
                New_Pneumo_NonSev=float(np.sum(exposure * c['lam_pneu']) * (1 - p['prop_pneumo_severe'])),
                New_Sepsis=float(np.sum(exposure * c['lam_sep'])),
                New_AOM=float(np.sum(exposure * c['lam_aom'])),
                New_Deaths_Pneumo=float(np.sum(exposure * c['mu_dis'])))
            year = month // 12 + 1
            event_df = r.econ.discount_factor(year - .5, e['discount_rate'])
            original_coeff = w6.build_week6_month_coefficients
            original_integrator = w6.solve_ivp
            try:
                # Every solver receives identical saved state and frozen coefficient arrays.
                w6.build_week6_month_coefficients = lambda *a, **k: c
                for method, rtol, atol in [('default_RK45', 1e-8, 1e-7),
                                           ('tighter_RK45', 1e-11, 1e-10)]:
                    solver_info = {}
                    def capture_solver(*a, **k):
                        sol = original_integrator(*a, **k)
                        solver_info.update(nfev=sol.nfev, status=sol.status)
                        return sol
                    w6.solve_ivp = capture_solver
                    after, stock, flows, all_deaths = w6.solve_one_month_week6(
                        A.copy(), dict(S), month, p, vaccine is not None, vaccine,
                        rtol=rtol, atol=atol)
                    delta = after - exact_A
                    death_error = flows['New_Deaths_Pneumo'] - exact_flows['New_Deaths_Pneumo']
                    row = dict(scenario=name, month=month, year=year, method=method,
                               rtol=rtol, atol=atol, **solver_info,
                               initial_alive=float(A.sum()),
                               max_abs_state_error=float(np.max(abs(delta))),
                               L1_state_error=float(np.sum(abs(delta))),
                               total_alive_error=float(after.sum() - exact_A.sum()),
                               max_abs_flow_error=max(abs(flows[k] - exact_flows[k]) for k in exact_flows),
                               exact_deaths=exact_flows['New_Deaths_Pneumo'],
                               computed_deaths=flows['New_Deaths_Pneumo'],
                               signed_death_event_error=death_error,
                               death_cost_weight=prod_weight,
                               calendar_discount=event_df,
                               signed_pv_death_cost_error=death_error * prod_weight * event_df,
                               signed_pv_death_qaly_error=death_error * qaly_weight * event_df,
                               saved_deaths=float(monthly.loc[month, 'New_Deaths_Pneumo']),
                               computed_minus_saved_deaths=flows['New_Deaths_Pneumo'] - float(monthly.loc[month, 'New_Deaths_Pneumo']))
                    rows.append(row)
                    for key in exact_flows:
                        flow_rows.append(dict(scenario=name, month=month, method=method,
                                              endpoint=key, exact=exact_flows[key],
                                              computed=flows[key], signed_error=flows[key] - exact_flows[key]))
            finally:
                w6.build_week6_month_coefficients = original_coeff
                w6.solve_ivp = original_integrator
    frame = pd.DataFrame(rows)
    frame.to_csv(out / 'fixed_input_month_oracle.csv', index=False)
    pd.DataFrame(flow_rows).to_csv(out / 'fixed_input_month_flow_errors.csv', index=False)
    conversion_frame = pd.DataFrame(conversions.values())
    conversion_frame.to_csv(out / 'probability_conversion_audit.csv', index=False)
    np.savez_compressed(out / 'frozen_month_coefficients.npz', **coefficient_snapshots)
    summary = {}
    for method, group in frame.groupby('method'):
        summary[method] = {key: float(group[key].abs().max()) for key in [
            'max_abs_state_error', 'L1_state_error', 'total_alive_error',
            'max_abs_flow_error', 'signed_death_event_error',
            'signed_pv_death_cost_error', 'signed_pv_death_qaly_error',
            'computed_minus_saved_deaths']}
        summary[method]['nfev_range'] = [int(group.nfev.min()), int(group.nfev.max())]
    config = io.StringIO()
    with contextlib.redirect_stdout(config):
        np.show_config()
        scipy.show_config()
    result = dict(scope='Four scenarios x four frozen saved months x two RK45 tolerances',
                  interpretation='Empirical fixed-input numerical diagnostic only; not a full feedback-trajectory bound or empirical calibration.',
                  platform=platform.platform(), python=sys.version,
                  versions=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
                  source_sha256=sources, summary=summary,
                  conversion_scope='80-digit Decimal reference conditional on the exact binary64 input; does not verify epidemiological incidence values.',
                  conversion_max_absolute_error={kind: {
                      key: float(group[key].abs().max()) for key in [
                          'legacy_signed_error', 'stable_signed_error']}
                      for kind, group in conversion_frame.groupby('conversion')},
                  mortality_accounting=r.MORTALITY_ACCOUNTING,
                  build_config=config.getvalue())
    (out / 'fixed_input_month_summary.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(dict(summary=summary, row_count=len(frame)), indent=2))


if __name__ == '__main__':
    main()
