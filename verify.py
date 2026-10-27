"""Check retained/recomputed numerical outputs, seeds and figure completeness."""
from pathlib import Path
import ast,hashlib,json,re
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
CHECKS=[]
def check(name,passed,detail=''):
    CHECKS.append(dict(check=name,passed=bool(passed),detail=str(detail)))

def main():
    for relative in ['research/inputs/Final/outputs/checks.csv',
                     'research/inputs/Final/outputs/independent_checks.csv',
                     'research/outputs/checks.csv','research/coverage/outputs/checks.csv']:
        q=pd.read_csv(ROOT/relative)
        check('numerical_checks_all_pass',q.passed.all(),f'{relative}: {len(q)}')
    info=json.loads((ROOT/'revision/outputs/ctmc_variability_checks.json').read_text())
    check('CTMC_distribution_checks',info['status']=='passed' and all(x['passed'] for x in info['checks']),info['check_count'])
    app=ROOT/'research/inputs/Final/outputs'
    check('application_run_count',len(list((app/'runs').iterdir()))==144)
    check('factorial_core_count',len(pd.read_csv(app/'factorial.csv'))==37)
    check('schedule_count',len(pd.read_csv(app/'schedules.csv'))==27)
    check('vaccinated_sensitivity_count',len(pd.read_csv(app/'sensitivity.csv'))==60)
    state_record=json.loads((ROOT/'checks/reference_state_arrays.json').read_text())
    bitwise_state_matches=0
    for name,record in state_record.items():
        with np.load(app/'runs'/name/'states.npz') as archive:
            a=archive['alive_after_births_and_doses']
        check('run_state_shape',list(a.shape)==record['shape'],name)
        check('run_state_finite',np.isfinite(a).all(),name)
        check('run_integrated_states',np.allclose(a.sum(axis=2),np.asarray(record['integrals']),rtol=1e-12,atol=1e-7),name)
        bitwise_state_matches+=hashlib.sha256(a.tobytes()).hexdigest()==record['array_sha256']
    refs=ROOT/'checks/reference_tables'
    comparisons=[]
    for p in sorted(refs.glob('*.csv')):
        relative=p.name.replace('__','/')
        current=ROOT/relative
        a=pd.read_csv(p);b=pd.read_csv(current)
        if p.name.endswith('checks.csv'):
            sort=[x for x in ['check','detail'] if x in a.columns]
            if sort:a=a.sort_values(sort).reset_index(drop=True);b=b.sort_values(sort).reset_index(drop=True)
        check('reference_shape',a.shape==b.shape,relative)
        check('reference_columns',list(a.columns)==list(b.columns),relative)
        if a.shape!=b.shape or list(a.columns)!=list(b.columns):continue
        maxerr=0.
        for col in a.columns:
            if col in ['seconds','elapsed_seconds']:continue
            if pd.api.types.is_numeric_dtype(a[col]) and pd.api.types.is_numeric_dtype(b[col]):
                av=a[col].to_numpy(dtype=float);bv=b[col].to_numpy(dtype=float)
                okay=np.allclose(av,bv,rtol=1e-12,atol=1e-7,equal_nan=True)
                error=float(np.nanmax(abs(av-bv))) if len(av) and np.isfinite(av-bv).any() else 0.
                maxerr=max(maxerr,error);check('reference_numeric_column',okay,f'{relative}: {col}; max_abs={error:g}')
            else:
                check('reference_identifier_column',a[col].fillna('<NA>').astype(str).equals(b[col].fillna('<NA>').astype(str)),f'{relative}: {col}')
        comparisons.append(dict(path=relative,rows=len(a),max_absolute_numeric_difference=maxerr))
    for p in list((ROOT/'research/src').glob('*.py'))+list((ROOT/'revision').glob('*.py'))+list(ROOT.glob('*.py')):
        tree=ast.parse(p.read_text(encoding='utf-8-sig'))
        check('entry_point_syntax',True,p.relative_to(ROOT))
        literals=[n.value for n in ast.walk(tree) if isinstance(n,ast.Constant) and isinstance(n.value,str)]
        drives=[x for x in literals if re.match(r'^[A-Za-z]:[\\/]',x)]
        check('no_drive_path_dependency',not drives,p.relative_to(ROOT))
    for n in range(1,7):check('main_figure_present',(ROOT/'figures'/f'Fig{n}.pdf').exists(),n)
    for prefix in ['pcv10_gsk','pcv13','pcv10_sii']:
        for kind in ['population_trajectory','dose_integrals','age_calendar_dose_density']:
            check('extended_figure_present',(ROOT/'figures/trajectories'/f'{prefix}_{kind}.pdf').exists(),prefix+'/'+kind)
    failed=[x for x in CHECKS if not x['passed']]
    result=dict(check_count=len(CHECKS),passed=len(CHECKS)-len(failed),failed=failed,
                run_state_arrays=len(state_record),bitwise_equal_state_arrays=bitwise_state_matches,
                reference_table_comparisons=comparisons,
                numeric_tolerance=dict(relative=1e-12,absolute=1e-7),
                exclusions='Execution times; CSV floating point roundtrip differences below tolerance. Reference tests do not prove empirical input validity.')
    (ROOT/'checks/verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    pd.DataFrame(CHECKS).to_csv(ROOT/'checks/verification.csv',index=False)
    print(json.dumps({k:v for k,v in result.items() if k!='reference_table_comparisons'},indent=2))
    if failed:raise AssertionError(failed)

if __name__=='__main__':main()
