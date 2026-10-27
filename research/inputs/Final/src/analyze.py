"""Revised factorial, schedule and perturbation analyses using one corrected solver.

Run from any directory: python analyze.py. No historical files are overwritten.
All perturbations are deterministic scenarios, not probability distributions.
"""
from pathlib import Path
import sys, copy, json, hashlib, itertools, time
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT.parent / 'Ablation'
sys.path.insert(0, str(OLD / 'src'))
import run_ablation as r
OUT = ROOT / 'outputs'
OUT.mkdir(parents=True, exist_ok=True)
r.OUTPUT = OUT
P0, E0 = copy.deepcopy(r.BASE_P), copy.deepcopy(r.ECON_P)
PARAMETERS, HERD = r.parameters, r.w6.herd_multiplier_for_endpoint
checks, rows, results = [], [], {}

def check(name, passed, detail=''):
    checks.append(dict(check=name, passed=bool(passed), detail=str(detail)))

def run(name, vaccine, d=1, w=1, g=1, off=False, p=None, e=None, fixed_herd=False):
    p = copy.deepcopy(p if p is not None else PARAMETERS(d,w))
    r.parameters = lambda *_: copy.deepcopy(p)
    r.ECON_P = copy.deepcopy(e if e is not None else E0)
    if fixed_herd:
        def fixed(pp,vv,endpoint,year,state):
            if vv is None: return 1.,0.,0.
            _, target, _ = HERD(P0,vv,endpoint,year,state)
            _, relative = r.w6.carriage_index(state, pp)
            strength = pp.get('carriage_to_herd_fraction',1.)
            return 1-(1-target)*relative*strength, target, relative
        r.w6.herd_multiplier_for_endpoint = fixed
    try:
        a = r.simulate(vaccine,d,w,g,name,herd_off=off)
        annual,costs,details = r.economics(a,vaccine,name)
        r.save_result(name,a)
        for filename, frame in [('yearly_unrounded',annual),('vaccine_costs',costs),('economic_details',details)]:
            frame.to_csv(OUT/'runs'/name/(filename+'.csv'),index=False)
        (OUT/'runs'/name/'economic_parameters.json').write_text(json.dumps(r.ECON_P,indent=2),encoding='utf-8')
        row = dict(scenario=name,vaccine=vaccine or 'No vaccination',D=d,W=w,G=g,herd_off=off,**r.summarize(a,details))
        checks.extend([dict(check='mass_balance',passed=bool(a['balance'].population_balance_error.abs().max()<1e-7),detail=name),
                       dict(check='dose_conservation',passed=bool(a['balance'].vaccination_balance_error.abs().max()<1e-7),detail=name),
                       dict(check='nonnegative',passed=bool(a['states'].min()>=0),detail=name),
                       dict(check='herd_bounds',passed=bool(a['herd'][['H_ipd','H_pneumonia','H_aom']].ge(0).all().all() and a['herd'][['H_ipd','H_pneumonia','H_aom']].le(1).all().all()),detail=name)])
        results[name]=a
        print(name,flush=True)
        return row
    finally:
        r.parameters=PARAMETERS; r.ECON_P=copy.deepcopy(E0); r.w6.herd_multiplier_for_endpoint=HERD

def contrast(row,nv):
    o=dict(row); scale=100000/row['births_denominator']
    for raw,metric,sgn in [('cases_total','cases_averted',-1),('deaths_total','deaths_averted',-1),('qaly_loss','qalys_gained',-1),('costs_total','incremental_costs',1)]:
        o[metric+'_per100k']=sgn*(row[raw]-nv[raw])*scale
    o['costs_total_per100k']=row['costs_total']*scale
    o['doses_per100k']=row['doses_total']*scale
    o['icer']=o['incremental_costs_per100k']/o['qalys_gained_per100k'] if abs(o['qalys_gained_per100k'])>1e-12 else np.nan
    return o

def main():
    start=time.time()
    nv=run('no_vaccination',None,0,0,0,True); rows.append(dict(nv,design='no_vaccination'))
    for v in r.PRODUCTS:
        for d,w,g in r.BITS:
            rows.append(dict(run(r.scenario_id(v,d,w,g),v,d,w,g),design='factorial'))
        for d,w in itertools.product((0,1),repeat=2):
            rows.append(dict(run(f'{r.PRODUCTS[v][0]}_D{d}W{w}_herd_off',v,d,w,0,True),design='herd_off'))
    f=r.add_contrasts(pd.DataFrame(rows)); f.to_csv(OUT/'factorial.csv',index=False)
    interactions,shapley,indirect=r.decompose(f)
    for name,frame in [('interactions',interactions),('shapley',shapley),('herd_contributions',indirect)]:frame.to_csv(OUT/(name+'.csv'),index=False)
    # Differences from the historical ledger are expected accounting revisions.
    audit=ROOT.parent.parent/'audit';audit.mkdir(parents=True,exist_ok=True)
    original_path=audit/'original_core.csv'
    old=pd.read_csv(original_path if original_path.exists() else OLD/'outputs'/'all_scenarios_summary.csv').set_index('scenario')
    if not original_path.exists():old.reset_index().to_csv(original_path,index=False)
    migration=[]
    metrics=['cases_total','deaths_total','qaly_loss','costs_total','doses_total',
             'cases_averted_per100k','deaths_averted_per100k','qalys_gained_per100k',
             'incremental_costs_per100k','icer']
    for _,row in f.iterrows():
        for metric in metrics:
            before=float(old.loc[row.scenario,metric]);after=float(row[metric])
            migration.append(dict(scenario=row.scenario,vaccine=row.vaccine,metric=metric,
                                  original=before,revised=after,delta=after-before,
                                  relative_change_pct=100*(after-before)/abs(before) if abs(before)>1e-12 else np.nan))
    pd.DataFrame(migration).to_csv(audit/'accounting_migration.csv',index=False)
    core=[]
    for _,row in f.iterrows():
        if row.scenario!='no_vaccination' and (row.D,row.W,row.G,row.herd_off) not in [(0,0,0,False),(1,1,1,False)]:
            continue
        detail=pd.read_csv(OUT/'runs'/row.scenario/'economic_details.csv')
        incident_adjustment=0.0
        for label,cost in [('hearing','hearing_case'),('epilepsy','epilepsy_case'),
                           ('hemiplegia','hemiplegia_case'),('implant','cochlear_implant_case')]:
            difference=detail['new_seq_'+label+'_incident']-detail['legacy_net_increment_seq_'+label]
            incident_adjustment+=float((difference*E0['costs_sequelae'][cost]/
                                        (1+E0['discount_rate'])**(detail.year-.5)).sum())
        before=old.loc[row.scenario]
        entry=dict(scenario=row.scenario,vaccine=row.vaccine,
                   decomposition_order='mortality correction, then incident one-off costs at revised states',
                   incident_flow_cost_adjustment_revised_states_total=incident_adjustment,
                   incident_flow_cost_adjustment_revised_states_per100k=incident_adjustment*1e5/row.births_denominator,
                   mortality_cost_change_before_incident_adjustment_total=row.costs_total-before.costs_total-incident_adjustment,
                   mortality_cost_change_before_incident_adjustment_per100k=(row.costs_total-before.costs_total-incident_adjustment)*1e5/row.births_denominator)
        for metric in metrics:
            entry['original_'+metric]=float(before[metric]);entry['revised_'+metric]=float(row[metric])
            entry['delta_'+metric]=float(row[metric]-before[metric])
        core.append(entry)
    core_frame=pd.DataFrame(core)
    nv_core=core_frame.loc[core_frame.scenario=='no_vaccination'].iloc[0]
    core_frame['incremental_incident_flow_cost_adjustment_per100k']=(
        core_frame.incident_flow_cost_adjustment_revised_states_per100k-
        nv_core.incident_flow_cost_adjustment_revised_states_per100k)
    core_frame['incremental_mortality_cost_change_before_incident_adjustment_per100k']=(
        core_frame.mortality_cost_change_before_incident_adjustment_per100k-
        nv_core.mortality_cost_change_before_incident_adjustment_per100k)
    check('accounting_migration_cost_closure',np.allclose(
        core_frame.delta_incremental_costs_per100k,
        core_frame.incremental_incident_flow_cost_adjustment_per100k+
        core_frame.incremental_mortality_cost_change_before_incident_adjustment_per100k,
        rtol=0.,atol=1e-7))
    core_frame.to_csv(OUT/'accounting_migration.csv',index=False)
    (audit/'mortality_accounting.json').write_text(json.dumps(r.MORTALITY_ACCOUNTING,indent=2),encoding='utf-8')
    check('cause_deleted_background_nonnegative',min(r.MORTALITY_ACCOUNTING['cause_deleted_background_hazard_by_month_of_age'])>=0)
    mu=r.MORTALITY_ACCOUNTING['postneonatal_monthly_all_cause_hazard']
    qN=r.MORTALITY_ACCOUNTING['neonatal_probability'];q5=r.MORTALITY_ACCOUNTING['under_five_probability']
    check('NV_under_five_probability_calibration',abs(1-(1-qN)*np.exp(-60*mu)-q5)<1e-14)
    for v in r.PRODUCTS:
        for metric in shapley.metric.unique():
            s=shapley[(shapley.vaccine==v)&(shapley.metric==metric)]
            check('shapley_closure',abs(s.shapley_delta.sum()-s.joint_delta.iloc[0])<1e-7,f'{v}/{metric}')

    # Schedule experiment holds W=1 and G=1 for every strategy.
    schedules=[]
    timings={'2-4-12':(2,4,12),'3-5-12':(3,5,12),'2-5-12':(2,5,12),'3-6-12':(3,6,12),'2-4-15':(2,4,15),'3-5-15':(3,5,15)}
    keys=['routine_dose1_transition_probs','routine_dose2_transition_probs','routine_booster_transition_probs']
    for v in r.PRODUCTS:
        for label in list(timings)+['delayed','dose1_only','completion_catchup']:
            p=PARAMETERS(1,1)
            if label in timings:
                for key,age in zip(keys,timings[label]):p[key]={age:1.}
            elif label=='dose1_only':
                p[keys[0]]={2:1.};p[keys[1]]={};p[keys[2]]={}
            elif label=='completion_catchup':
                # Hypothetical recovery for enrolled children; no enrolment of older non-programme children.
                # Descending event order prevents same-month receipt of two doses.
                for key,ages in zip(keys,[range(6,12),range(9,18),range(17,24)]):
                    p[key].update({age:.20 for age in ages})
            name=r.PRODUCTS[v][0]+'_schedule_'+label.replace('-','_')
            row=run(name,v,p=p)
            schedules.append(dict(contrast(row,nv),strategy=label))
    s=pd.DataFrame(schedules)
    for v in r.PRODUCTS:
        b=s[(s.vaccine==v)&(s.strategy=='2-4-12')].iloc[0]
        for metric in ['cases_averted_per100k','deaths_averted_per100k','qalys_gained_per100k','incremental_costs_per100k','doses_per100k']:
            s.loc[s.vaccine==v,'delta_'+metric]=s.loc[s.vaccine==v,metric]-b[metric]
    s.to_csv(OUT/'schedules.csv',index=False)

    # Paired NV comparator is recomputed under each changed parameter ledger.
    sensitivity=[]
    specs=[('price',.75),('price',1.25),('direct_VE',.75),('direct_VE',1.25),('CFR',.75),('CFR',1.25),
           ('sequela_cost',.75),('sequela_cost',1.25),('utility_loss',.75),('utility_loss',1.25),
           ('burden',.75),('burden',1.25),('waning',0.),('waning',1.25),('herd_strength',.5),('herd_strength',0.),
           ('delay_probability',.75),('delay_probability',1.25),('birth_seasonality',0.),('healthcare_perspective',0.)]
    nv_cache={}
    for parameter,value in specs:
        p=PARAMETERS(1,1); e=copy.deepcopy(E0)
        if parameter=='price':e['vaccine_prices']={k:x*value for k,x in e['vaccine_prices'].items()}
        elif parameter=='direct_VE':p['vaccine_direct_effectiveness']={v:{k:x*value for k,x in ve.items()} for v,ve in p['vaccine_direct_effectiveness'].items()}
        elif parameter=='CFR':
            for k in ['cfr_meningitis','cfr_pneumo_severe','cfr_sepsis']:p[k]*=value
        elif parameter=='sequela_cost':e['costs_sequelae']={k:x*value for k,x in e['costs_sequelae'].items()}
        elif parameter=='utility_loss':e['utilities']={k:1-(1-x)*value for k,x in e['utilities'].items()}
        elif parameter=='burden':
            for k in ['inc_meningitis','inc_pneumonia','inc_sepsis','inc_aom_all']:p[k]*=value
        elif parameter=='waning':
            for key,peak in [('primary_waning_by_age',1.),('booster_waning_by_age',1.1)]:p[key]={age:peak-(peak-x)*value for age,x in p[key].items()}
        elif parameter=='herd_strength':p['carriage_to_herd_fraction']=value
        elif parameter=='delay_probability':
            for key in keys:p[key]={age:min(1.,prob*value) for age,prob in p[key].items()}
        elif parameter=='birth_seasonality':p['birth_month_shares']=[1/12]*12
        elif parameter=='healthcare_perspective':
            for k in ['costs_nonmedical','costs_caregiver','prod_discounted_age_contrib']:e[k]={a:0. for a in e[k]}
            e['prod_loss_per_death_age0']=0.
        tag=parameter+'_'+str(value).replace('.','p')
        changed_nv=run('nv_'+tag,None,0,0,0,True,p=p,e=e)
        for v in r.PRODUCTS:
            name=r.PRODUCTS[v][0]+'_sens_'+tag
            row=contrast(run(name,v,p=p,e=e,fixed_herd=parameter=='direct_VE'),changed_nv)
            base=f[f.scenario==r.scenario_id(v,1,1,1)].iloc[0]
            for metric in ['cases_averted_per100k','deaths_averted_per100k','qalys_gained_per100k','incremental_costs_per100k','icer']:
                row['delta_'+metric]=row[metric]-base[metric]
            row.update(parameter=parameter,value=value,reference='D1W1G1',fixed_target_herd=parameter=='direct_VE')
            sensitivity.append(row)
    pd.DataFrame(sensitivity).to_csv(OUT/'sensitivity.csv',index=False)

    # Trajectory summaries use saved PRE-reaction states: births and doses applied.
    diagnostics=[];traces=[];stable=[]
    for name,a in results.items():
        if '_sens_' in name or name.startswith('nv_'):continue
        state=a['states'];h=a['herd'];total=state.sum(axis=(1,2));groups=state.sum(axis=2)
        trace=h.copy();trace['scenario']=name
        for i,group in enumerate(r.w6.GROUP_NAMES):trace[group]=groups[:,i]
        trace['alive_pre_reaction']=total
        traces.append(trace)
        row={'scenario':name,'final_C':h.carriage_index.iloc[-1],'final_r':h.relative_carriage_index.iloc[-1],'final_H_ipd':h.H_ipd.iloc[-1]}
        for threshold in [.5,.9]:
            hit=np.flatnonzero(h.relative_carriage_index.to_numpy()>=threshold)
            row[f'first_r_{threshold}_years']=hit[0]/12 if len(hit) else np.nan
        for i,group in enumerate(r.w6.GROUP_NAMES):row['final_share_'+group]=groups[-1,i]/total[-1]
        diagnostics.append(row)
    pd.DataFrame(diagnostics).to_csv(OUT/'trajectory_diagnostics.csv',index=False)
    pd.concat(traces,ignore_index=True).to_csv(OUT/'trajectories.csv',index=False)
    nv_month=results['no_vaccination']['monthly']
    for v in r.PRODUCTS:
        a=results[r.scenario_id(v,1,1,1)]['monthly']
        for endpoint in r.w6.FLOW_NAMES:
            b=nv_month.loc[nv_month.Year==10,endpoint].sum(); n=a.loc[a.Year==10,endpoint].sum()
            stable.append(dict(vaccine=v,endpoint=endpoint,baseline=b,vaccinated=n,reduction_pct=100*(b-n)/b))
    pd.DataFrame(stable).to_csv(OUT/'stable.csv',index=False)
    # State integration and counterfactual identities.
    for name,a in results.items():
        check('state_integral',np.allclose(a['states'].sum(axis=(1,2)),a['states'].sum(axis=2).sum(axis=1)),name)
    pd.DataFrame(checks).to_csv(OUT/'checks.csv',index=False)
    assert all(c['passed'] for c in checks), 'Validation failure; inspect checks.csv'
    source_paths=list((OLD/'inputs'/'project').rglob('*.py'))+[OLD/'src'/'run_ablation.py',Path(__file__)]
    manifest=[dict(path=str(p.relative_to(ROOT.parent.parent)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in source_paths]
    (OUT/'provenance.json').write_text(json.dumps(dict(sources=manifest,elapsed_seconds=time.time()-start,checks=len(checks),python=sys.version,psa=False),indent=2),encoding='utf-8')
    print('DONE',len(checks),'checks passed',flush=True)

if __name__=='__main__':main()
