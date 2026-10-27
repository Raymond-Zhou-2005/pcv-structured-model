"""Reconstruct legacy verification and sensitivity evidence from frozen sources.

Writes only within this coverage directory. Legacy cohort economic outputs are
separate estimands and are never used as present calendar-PDE estimates.
"""
from __future__ import annotations
import copy
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'inputs'
OUTPUT = ROOT / 'outputs'
OUTPUT.mkdir(exist_ok=True)
sys.dont_write_bytecode = True
sys.path.insert(0, str(INPUT))

def module(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, INPUT / file)
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj

def dump(frame: pd.DataFrame, name: str):
    frame.to_csv(OUTPUT / name, index=False)

def main():
    checks = []
    for rec in json.loads((ROOT / 'source_manifest.json').read_text()):
        file = INPUT / rec['frozen_file']
        matches = hashlib.sha256(file.read_bytes()).hexdigest() == rec['sha256']
        checks.append(dict(check='frozen_source_hash',item=rec['frozen_file'],passed=matches))
    ode = module('ode_new', 'ode_new.py')
    tier = module('legacy_assumption_tiers', 'Week 1_assumption_tiers.py')
    exact = tier.ex
    accounting = tier.run_accounting_checks(str(OUTPUT / 'week1_accounting.csv'))
    validation = exact.validate_all_scenarios(str(OUTPUT / 'week1_exact_rk45.csv'))
    elasticity = tier.run_parameter_elasticities(output_csv=str(OUTPUT / 'week1_elasticities.csv'))
    structure = tier.run_structural_screens(output_csv=str(OUTPUT / 'week1_structure.csv'))
    for rec in accounting.to_dict('records'):
        label, value = rec['Check'],rec['Value']
        if label.startswith('weight_sum'):
            passed = abs(value-1.0) < 1e-12
        elif label.startswith('horizon_prefix'):
            passed = value == 0.0
        else:
            # Exact/RK45 values are reported, not judged with an invented
            # tolerance detached from the original solver's tolerances.
            passed = np.isfinite(value)
        checks.append(dict(check='week1_accounting',item=label,value=value,passed=bool(passed)))
    # Algebraic identities of local log sensitivities, independently recomputed
    # from saved perturbed cases and deaths.
    for rec in elasticity.to_dict('records'):
        low, high = (0.95,1.05) if rec['Parameter']=='vaccine_coverage' else (0.9,1.1)
        for outcome in ['Cases','Deaths']:
            elasticity_direct=(np.log(abs(rec[f'High{outcome}Averted'])+1e-12)-np.log(abs(rec[f'Low{outcome}Averted'])+1e-12))/(np.log(high)-np.log(low))
            err=abs(elasticity_direct-rec[f'{outcome}AvertedElasticity'])
            checks.append(dict(check='local_log_elasticity',item=rec['Parameter']+'/'+outcome,value=err,passed=bool(err<1e-12)))
    # Independently execute the original PCV six-scenario characteristic code.
    pde = module('legacy_pcv_characteristic', 'Week 4_pde_validation.py')
    rows=[]
    for scenario in exact.SCENARIO_MAP:
        p, intervention, product = exact.scenario_params(copy.deepcopy(ode.params_vietnam_micro),scenario)
        a=ode.run_ode_and_build_monthly(p,intervention=intervention,vaccine_type=product)
        b=pde.run_age_time_pde_cohort_comparable(p,intervention=intervention,vaccine_type=product)
        for col in a.columns.intersection(b.columns):
            if col == 'Year':
                continue
            err=float(np.abs(a[col].to_numpy()-b[col].to_numpy()).max())
            rows.append(dict(scenario=scenario,metric=col,max_absolute_monthly_difference=err))
        ya=ode.aggregate_yearly_like_micro(a)
        yb=pde.aggregate_yearly_like_micro(b)
        for col in ya.columns.intersection(yb.columns):
            if col=='Year':
                continue
            err=float(np.abs(ya[col].to_numpy()-yb[col].to_numpy()).max())
            checks.append(dict(check='week4_rounded_yearly',item=scenario+'/'+col,value=err,passed=bool(err==0.0)))
    dump(pd.DataFrame(rows),'week4_pcv_characteristic.csv')
    # The separate Week8 limit is an explicit monthly scalar-vs-age-bin
    # transport check, not a claim about a fitted epidemic model.
    limit=module('legacy_ode_limit','validate_closed_cohort_ode_limit.py')
    p=limit.ModelParameters();rows=[]
    for scenario,mult in [('baseline',1.0),('mature_pcv13_population_irr',p.mature_population_irr['PCV13']['ipd'])]:
        a=limit.scalar_ode_characteristic(p,mult);b=limit.age_time_pde_transport(p,mult)
        for metric in ['alive','meningitis','pneumonia','sepsis','aom','deaths']:
            err=float((a[metric]-b[metric]).abs().max())
            rows.append(dict(scenario=scenario,metric=metric,max_absolute_difference=err))
            checks.append(dict(check='week8_simplified_ode_limit',item=scenario+'/'+metric,value=err,passed=bool(err<1e-8)))
    dump(pd.DataFrame(rows),'week8_ode_limit.csv')
    # Re-evaluate all 17 historical one-way cohort scenarios. Labels such as
    # "published range proxy" are inherited labels, not independently verified
    # uncertainty intervals.
    dsa=module('legacy_ode_dsa','run_ode_deterministic_sensitivity.py')
    saved_argv=sys.argv
    try:
        sys.argv=['run_ode_deterministic_sensitivity.py','--output',str(OUTPUT/'week8_cohort_dsa.csv')]
        dsa.main()
    finally:
        sys.argv=saved_argv
    before=pd.read_csv(INPUT/'historical_ode_dsa.csv');after=pd.read_csv(OUTPUT/'week8_cohort_dsa.csv')
    merged=before.merge(after,on=['parameter','level'],suffixes=('_historical','_reconstructed'),validate='one_to_one')
    rows=[]
    for metric in ['cases_averted','deaths_averted','qaly_gained','incremental_cost','icer']:
        err=float(np.abs(merged[metric+'_historical']-merged[metric+'_reconstructed']).max())
        rows.append(dict(metric=metric,max_absolute_difference=err))
        checks.append(dict(check='week8_saved_dsa_reconstruction',item=metric,value=err,passed=bool(err<1e-7)))
    dump(pd.DataFrame(rows),'week8_dsa_reconstruction.csv')
    # Retained bridge constants, with a clearly separated reproducible product
    # multiplication and unrecoverable row-level immunogenicity extraction.
    m=0.8203164085517026
    bridge=[]
    for endpoint,ve13 in [('ipd',0.602),('pneumonia',0.19),('aom',0.46)]:
        bridge.append(dict(endpoint=endpoint,pcv13_direct_ve=ve13,sii_relative_multiplier=m,sii_direct_ve=m*ve13,status='retained immunogenicity proxy; original row-level weighting unavailable'))
    dump(pd.DataFrame(bridge),'sii_bridge.csv')
    checks.append(dict(check='visible_coverage_ratio',item='76/109',value=76/109,passed=bool(abs(76/109-0.6972477064220183)<1e-15)))
    result=pd.DataFrame(checks);dump(result,'checks.csv')
    manifest=dict(run_kind='new reconstruction from frozen legacy code, not recovered original raw outputs',
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  python=sys.version,checks=len(result),passed=int(result['passed'].sum()),
                  exact_rk45_max_abs=float(validation['MaxAbsDiff'].max()),
                  exact_rk45_max_scaled=float(validation['MaxRelDiff'].max()),
                  week4_max_monthly_abs=float(pd.read_csv(OUTPUT/'week4_pcv_characteristic.csv')['max_absolute_monthly_difference'].max()),
                  calendar_main_results_changed=False,
                  psa_not_rerun=True)
    (OUTPUT/'provenance.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    if not result['passed'].all():
        raise AssertionError(result.loc[~result['passed']].to_string(index=False))
    print(json.dumps(manifest,indent=2))

if __name__=='__main__':
    main()
