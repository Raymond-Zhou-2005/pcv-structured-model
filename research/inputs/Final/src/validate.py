"""Independent accounting and frozen-coefficient checks for the final release."""
from pathlib import Path
import sys,json
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'Ablation'/'src'))
import run_ablation as m
OUT=ROOT/'outputs'; checks=[]
def check(name,error,tolerance,detail):
    checks.append(dict(check=name,error=float(error),tolerance=tolerance,passed=bool(abs(error)<=tolerance),detail=detail))
for directory in sorted((OUT/'runs').iterdir()):
    name=directory.name
    p=json.loads((directory/'parameters.json').read_text()); e=json.loads((directory/'economic_parameters.json').read_text())
    a=np.load(directory/'states.npz')['alive_after_births_and_doses']
    monthly=pd.read_csv(directory/'monthly.csv');h=pd.read_csv(directory/'herd.csv');d=pd.read_csv(directory/'doses.csv');details=pd.read_csv(directory/'economic_details.csv')
    v=next((v for v,(prefix,_) in m.PRODUCTS.items() if name.startswith(prefix+'_')),None)
    check('dose_order',max(0,d.dose2.sum()-d.dose1.sum(),d.booster.sum()-d.dose2.sum()),1e-7,name)
    weights=np.array([0.]*6+[1.]*12+[.8]*6+[.35]*36);rho=np.array([0,0,.15,.34,.4])
    index=(a*weights[None,None,:]*rho[None,:,None]).sum((1,2))/(a*weights[None,None,:]).sum((1,2))
    check('independent_C',np.max(abs(index-h.carriage_index)),1e-12,name)
    check('integrated_population',np.max(abs(a.sum((1,2))-a.sum(2).sum(1))),1e-7,name)
    for month in [0,36,119]:
        original=m.w6.herd_multiplier_for_endpoint
        m.w6.herd_multiplier_for_endpoint=lambda pp,vv,endpoint,year,state:(h.loc[month,'H_'+endpoint],h.loc[month,'H_target_'+endpoint],h.loc[month,'relative_carriage_index'])
        try:c=m.w6.build_week6_month_coefficients(p,month,a[month],v is not None,v)
        finally:m.w6.herd_multiplier_for_endpoint=original
        mu=c['mu_all'];exposure=a[month]*(-np.expm1(-mu))/mu
        expected={'New_Meningitis':np.sum(exposure*c['lam_men']),'New_Pneumo_Sev':np.sum(exposure*c['lam_pneu'])*p['prop_pneumo_severe'],
                  'New_Pneumo_NonSev':np.sum(exposure*c['lam_pneu'])*(1-p['prop_pneumo_severe']),'New_Sepsis':np.sum(exposure*c['lam_sep']),
                  'New_AOM':np.sum(exposure*c['lam_aom']),'New_Deaths_Pneumo':np.sum(exposure*c['mu_dis']),'Stock_Alive':np.sum(a[month]*np.exp(-mu))}
        check('closed_form_reaction',max(abs(monthly.loc[month,k]-x) for k,x in expected.items()),1e-5,name+'/'+str(month))
    annual=monthly.groupby('Year')['New_Deaths_Pneumo'].sum()
    survivor=monthly.New_Meningitis*(1-p['cfr_meningitis'])
    for column,fraction in [('Cognitive',p['prob_cog']),('Hearing',p['prob_hear']),
                            ('Epilepsy',p['prob_epi']),('Hemiplegia',p['prob_hemi']),
                            ('Hearing_Implant',p['prob_hear']*p['prob_cochlear_implant']),
                            ('Hearing_SpEd',p['prob_hear']*p['prob_special_education'])]:
        check('actual_incident_sequela_source',np.max(abs(monthly['New_Seq_'+column]-survivor*fraction)),1e-8,name+'/'+column)
    annual_incidence=monthly.groupby('Year')[[c for c in monthly if c.startswith('New_Seq_')]].sum()
    annual_stock=monthly.groupby('Year')[[c for c in monthly if c.startswith('Stock_Seq_')]].last()
    for i,row in details.iterrows():
        year=i+1;rate=e['discount_rate'];weight=sum((1+rate)**(-(k+.5)) for k in range(71))+.5*(1+rate)**(-71.25)
        expected=annual.loc[year]*weight/(1+rate)**(year-.5)
        check('calendar_death_discount',row.pv_qaly_death-expected,1e-7,name+'/'+str(year))
        doses=d.loc[d.Year==year,['dose1','dose2','booster']].to_numpy().sum()
        if v:
            price=e['vaccine_prices'][m.PRODUCTS[v][1]]
            unit=price*(1+e['wastage_rate'])+e['service_delivery_per_dose']+e['ae_incidence_per_100k_doses']/1e5*e['ae_cost_per_case']
            expected=doses*unit/(1+rate)**(year-1+.35)
        else:expected=0.
        check('actual_dose_cost',row.pv_vaccine_program_cost-expected,1e-6,name+'/'+str(year))
        c=e['costs_sequelae'];s=annual_stock.loc[year];inc=annual_incidence.loc[year]
        expected=(s.Stock_Seq_Cognitive*c['cognitive_annual']+s.Stock_Seq_Hearing_SpEd*c['special_education_annual']+
                  inc.New_Seq_Hearing*c['hearing_case']+inc.New_Seq_Epilepsy*c['epilepsy_case']+
                  inc.New_Seq_Hemiplegia*c['hemiplegia_case']+inc.New_Seq_Hearing_Implant*c['cochlear_implant_case'])
        check('incident_based_oneoff_sequela_cost',row.raw_sequela_cost-expected,1e-7,name+'/'+str(year))
        check('economic_cost_component_sum',row.pv_total_cost-(row.pv_disease_cost+row.pv_sequela_cost+
              row.pv_death_prod_cost+row.pv_vaccine_program_cost),1e-6,name+'/'+str(year))
cal=m.MORTALITY_ACCOUNTING
qN=cal['neonatal_probability'];q5=cal['under_five_probability']
mu=np.asarray(cal['cause_deleted_background_hazard_by_month_of_age'])+np.asarray(cal['baseline_NV_disease_hazard_by_month_of_age'])
check('NV_U5_calibration',1-(1-qN)*np.exp(-mu.sum())-q5,1e-14,'baseline all-cause birth-to-60-month probability')
check('background_nonnegative',min(0.,min(cal['cause_deleted_background_hazard_by_month_of_age'])),0.,'all 60 age cells')
f=pd.read_csv(OUT/'factorial.csv');sh=pd.read_csv(OUT/'shapley.csv');j=pd.read_csv(OUT/'interactions.csv');h=pd.read_csv(OUT/'herd_contributions.csv')
for (v,metric),group in j.groupby(['vaccine','metric']):
    b=f[f.scenario==m.scenario_id(v,0,0,0)].iloc[0][metric];end=f[f.scenario==m.scenario_id(v,1,1,1)].iloc[0][metric]
    check('interaction_reconstruction',group.interaction_delta.sum()-(end-b),1e-7,v+'/'+metric)
for metric in ['cases_averted_per100k','deaths_averted_per100k','qalys_gained_per100k']:
    check('indirect_identity',np.max(abs(h['indirect_'+metric]+h['direct_only_'+metric]-h['total_'+metric])),1e-8,metric)
s=pd.read_csv(OUT/'schedules.csv')
for v in m.PRODUCTS:
    for label,bits in [('2-4-12',(0,1,1)),('delayed',(1,1,1))]:
        a=s[(s.vaccine==v)&(s.strategy==label)].iloc[0];b=f[f.scenario==m.scenario_id(v,*bits)].iloc[0]
        for metric in ['cases_averted_per100k','deaths_averted_per100k','qalys_gained_per100k','incremental_costs_per100k']:
            check('schedule_factorial_bridge',a[metric]-b[metric],1e-7,v+'/'+label+'/'+metric)
df=pd.DataFrame(checks);df.to_csv(OUT/'independent_checks.csv',index=False)
print(df.groupby('check').agg(checks=('passed','size'),passed=('passed','sum'),max_error=('error',lambda s:s.abs().max())).to_string())
assert df.passed.all(), df[~df.passed].to_string()
