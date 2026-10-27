"""Mathematical verification and a refinement bridge for the frozen PCV model.

All ages/times are months. Scheduled-dose probabilities are impulses at integer
ages, NOT continuous rates. Refinement subdivides births/transport/feedback;
it preserves every visit probability and every inherited monthly risk mapping.
"""
from pathlib import Path
import sys, json, time, copy, hashlib, ast
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
PCV = ROOT.parent.parent
OUT = ROOT/'outputs'
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT/'inputs'/'Ablation'/'src'))
import run_ablation as r
w6 = r.w6
checks=[]

def check(name, error, tolerance, detail=''):
    checks.append(dict(check=name,error=float(error),tolerance=float(tolerance),
                       passed=bool(abs(error)<=tolerance),detail=str(detail)))

def math_checks():
    # Row generator, explicitly including death; probability is a column vector.
    Q=np.array([[-.21,.20,0,0,.01],[0,-.162,.15,0,.012],
                [0,0,-.098,.08,.018],[0,0,0,-.022,.022],[0,0,0,0,0.]])
    p0=np.array([1.,0,0,0,0]); rows=[]
    for t in [0.25,1.,3.,12.]:
        exact=expm(Q.T*t)@p0
        numerical=solve_ivp(lambda _,x:Q.T@x,[0,t],p0,rtol=1e-11,atol=1e-13).y[:,-1]
        check('CTMC_forward_matrix_exponential',np.max(abs(exact-numerical)),1e-10,t)
        check('CTMC_probability_mass',exact.sum()-1,1e-12,t)
        draws=np.random.default_rng(20261004+int(t*100)).multinomial(100000,exact,size=1000)/100000
        mean=draws.mean(0);se=draws.std(0,ddof=1)/np.sqrt(len(draws))
        for i in range(5):
            rows.append(dict(month=t,state=i,exact=exact[i],ode=numerical[i],
                             simulated_mean=mean[i],mc_se=se[i],z=(mean[i]-exact[i])/se[i] if se[i] else 0.))
    pd.DataFrame(rows).to_csv(OUT/'ctmc.csv',index=False)
    # A visit acts once: a discrete Markov kernel followed by survival.
    J=np.eye(5);J[0,0]=.35;J[1,0]=.65
    x=np.array([1000.,200.,100.,50.,0.]);mu=np.array([.01,.012,.015,.018,0.])
    P=np.diag(np.exp(-mu))@J
    # Add death to the final row, making the transition kernel stochastic.
    P[-1,:]+=1-P.sum(0)
    check('impulse_kernel_mass',np.max(abs(P.sum(0)-1)),1e-12)
    check('impulse_kernel_nonnegative',min(0.,P.min()),1e-12)
    check('impulse_expected_survival',np.max(abs((P@x)[:4]-(J@x)[:4]*np.exp(-mu[:4]))),1e-10)
    # Finite-particle nonlinear feedback: expectation does not commute.
    draws=np.random.default_rng(20261004).binomial(100,.5,size=100000)/100
    (OUT/'nonlinear_expectation.json').write_text(json.dumps(dict(
        E_X_squared=float(np.mean(draws**2)),E_X_squared_exact=.2525,
        square_E_X=.25,explanation='E[X^2] != E[X]^2; this is an illustrative closure counterexample.'),indent=2))

def aggregation():
    p=r.parameters(1,1); rows=[];states=[]
    for label,swap in [('young_primary_old_booster',False),('old_primary_young_booster',True)]:
        A=np.zeros((5,60));A[0,20]=10000;A[1,4]=10000;A[2,8]=10000
        A[3,48 if swap else 16]=10000;A[4,16 if swap else 48]=10000
        states.append(A)
        c=w6.build_week6_month_coefficients(p,60,A,True,'PCV13')
        mu=c['mu_all'];exposure=A*(-np.expm1(-mu))/mu
        C,relative=w6.carriage_index(A,p)
        H=w6.herd_multiplier_for_endpoint(p,'PCV13','ipd',5,A)[0]
        rows.append(dict(configuration=label,alive=A.sum(),C=C,r=relative,H_ipd=H,
                         initial_alive_derivative=-np.sum(mu*A),
                         one_month_cases=sum(np.sum(exposure*c[k]) for k in ['lam_men','lam_pneu','lam_sep','lam_aom']),
                         one_month_deaths=np.sum(exposure*c['mu_dis']),
                         initial_ipd_rate=np.sum(A*(c['lam_men']+c['lam_sep']))))
    check('aggregation_same_dose_totals',np.max(abs(states[0].sum(1)-states[1].sum(1))),1e-12)
    check('aggregation_same_age_totals',np.max(abs(states[0].sum(0)-states[1].sum(0))),1e-12)
    assert abs(rows[0]['one_month_cases']-rows[1]['one_month_cases'])>1
    pd.DataFrame(rows).to_csv(OUT/'aggregation.csv',index=False)
    np.savez_compressed(OUT/'aggregation_states.npz',states=np.asarray(states))

def manufactured():
    # Exact characteristic solution with constant death, time-varying boundary.
    # A0=0, B(t)=100*(1+.2*sin(t)), mu=.03, t=8, a_max=12.
    rows=[];mu=.03;T=8.;L=12.
    from scipy.integrate import quad
    exact_total=quad(lambda a:100*(1+.2*np.sin(T-a))*np.exp(-mu*a),0,T,epsabs=1e-10)[0]
    for k in [1,2,4,8,16,32]:
        h=1/k;A=np.zeros(int(L*k));start=time.perf_counter()
        for n in range(int(T*k)):
            A[0]+=h*100*(1+.2*np.sin(n*h))
            A*=np.exp(-mu*h)
            A[1:]=A[:-1].copy();A[0]=0.
        age=np.arange(len(A))*h
        exact_mass=np.array([quad(lambda a:100*(1+.2*np.sin(T-a))*np.exp(-mu*a),lo,min(lo+h,T),epsabs=1e-10)[0] if lo<T else 0. for lo in age])
        err=np.sum(abs(A-exact_mass))
        rows.append(dict(k=k,h=h,L1_density_error=err,total_error=abs(A.sum()-exact_total),
                         alive=A.sum(),exact_total=exact_total,seconds=time.perf_counter()-start))
    df=pd.DataFrame(rows);df['observed_total_order']=np.log2(df.total_error.shift()/df.total_error)
    df.to_csv(OUT/'manufactured.csv',index=False)
    check('manufactured_refinement_decrease',max(0.,df.total_error.iloc[-1]-df.total_error.iloc[0]),1e-12)

def simulate_refined(vaccine,k,*,d=1,w=1,g=1,mortality_mode='cause_deleted',aom_fraction=None):
    start=time.perf_counter();h=1/k;p=r.parameters(d,w);shares=w6.base.normalize_birth_shares(p['birth_month_shares'])
    A=np.zeros((5,60*k));S=w6.base.empty_sequela_stock();history=[];dosehistory=[];traces=[]
    if aom_fraction is not None:p['prop_aom_spn']=float(aom_fraction)
    p['mortality_accounting']=mortality_mode
    original_herd=w6.herd_multiplier_for_endpoint
    if not g:
        w6.herd_multiplier_for_endpoint=lambda pp,vv,e,y,aa:(original_herd(pp,vv,e,y,aa)[1] if vv else 1.,original_herd(pp,vv,e,y,aa)[1],original_herd(pp,vv,e,y,aa)[2])
    def coefficient(month,intervention):
        collapsed=A.reshape(5,60,k).sum(2)
        cc=w6.build_week6_month_coefficients(p,month,collapsed,intervention,vaccine if intervention else None)
        return {key:np.repeat(val,k,axis=1) if val.ndim==2 else val for key,val in cc.items()}
    def step(n,intervention,record):
        nonlocal A,S
        month=n//k;births=p['birth_cohort']*shares[month%12]*h
        incoming=births*w6.base.neonatal_survival_multiplier(p);before=A.sum()
        if intervention:A[1,0]+=incoming*p['vaccine_coverage'];A[0,0]+=incoming*(1-p['vaccine_coverage'])
        else:A[0,0]+=incoming
        counts=dict(dose1=0.,dose2=0.,booster=0.)
        if intervention:
            for src,target,key,count in [(3,4,'routine_booster_transition_probs','booster'),(2,3,'routine_dose2_transition_probs','dose2'),(1,2,'routine_dose1_transition_probs','dose1')]:
                for age,prob in p[key].items():
                    j=int(age)*k;moved=A[src,j]*prob;A[src,j]-=moved;A[target,j]+=moved;counts[count]+=moved
        if record and n%k==0:
            collapsed=A.reshape(5,60,k).sum(2);C,relative=w6.carriage_index(collapsed,p)
            H=original_herd(p,vaccine,'ipd',month//12,collapsed)[0] if intervention else 1.
            traces.append(dict(Month=month,C=C,r=relative,H_ipd=H,alive=A.sum(),**{f'group{i}':A[i].sum() for i in range(5)}))
        cc=coefficient(month,intervention);mu=cc['mu_all'];exposure=A*(-np.expm1(-mu*h))/mu
        flows=dict(New_Meningitis=np.sum(exposure*cc['lam_men']),
                   New_Pneumo_Sev=np.sum(exposure*cc['lam_pneu'])*p['prop_pneumo_severe'],
                   New_Pneumo_NonSev=np.sum(exposure*cc['lam_pneu'])*(1-p['prop_pneumo_severe']),
                   New_Sepsis=np.sum(exposure*cc['lam_sep']),New_AOM=np.sum(exposure*cc['lam_aom']),
                   New_Deaths_Pneumo=np.sum(exposure*cc['mu_dis']))
        after=A*np.exp(-mu*h);deaths=A.sum()-after.sum();surv=after.sum()/A.sum() if A.sum() else 1.
        m=flows['New_Meningitis']*(1-p['cfr_meningitis'])
        sources=dict(Stock_Seq_Cognitive=p['prob_cog']*m,Stock_Seq_Hearing=p['prob_hear']*m,
                     Stock_Seq_Epilepsy=p['prob_epi']*m,Stock_Seq_Hemiplegia=p['prob_hemi']*m,
                     Stock_Seq_Hearing_Implant=p.get('prob_cochlear_implant',0)*p['prob_hear']*m,
                     Stock_Seq_Hearing_SpEd=p.get('prob_special_education',0)*p['prob_hear']*m)
        S={key:S[key]*surv+sources[key] for key in S}
        agedout=after[:,-1].sum();A=np.zeros_like(after);A[:,1:]=after[:,:-1]
        if record:
            history.append(dict(Month=month,Year=month//12+1,**flows,**r.incident_sequela_flows(p,flows),Stock_Alive=after.sum(),
                                Gross_Births=births,Monthly_Total_Cases=sum(flows[x] for x in w6.FLOW_NAMES if x!='New_Deaths_Pneumo'),**S))
            dosehistory.append(dict(Month=month,Year=month//12+1,**counts))
            check('refined_mass_balance',A.sum()-(before+incoming-deaths-agedout),1e-7,f'{vaccine}/{k}/{n}')
    try:
        for n in range(60*k):step(n,False,False)
        S=w6.base.empty_sequela_stock()
        for n in range(120*k):step(n,vaccine is not None,True)
    finally:
        w6.herd_multiplier_for_endpoint=original_herd
    frame=pd.DataFrame(history);stocks=[x for x in frame if x.startswith('Stock_')]
    flows=[c for c in frame if c.startswith('New_')]+['Gross_Births','Monthly_Total_Cases']
    monthly=frame.groupby(['Month','Year'],as_index=False)[flows].sum()
    last=frame.groupby('Month')[stocks].last().reset_index();monthly=monthly.merge(last,on='Month')
    result=dict(monthly=monthly,doses=pd.DataFrame(dosehistory).groupby(['Month','Year'],as_index=False).sum())
    annual=monthly.groupby('Year')[[c for c in monthly if c.startswith('New_')]].sum().reset_index()
    for key in w6.SEQ_KEYS:annual[key]=monthly.groupby('Year')[key].last().to_numpy()
    counts={int(year):group[['dose1','dose2','booster']].sum().to_dict() for year,group in result['doses'].groupby('Year')}
    key=r.PRODUCTS[vaccine][1] if vaccine else None
    costs=pd.DataFrame(r.econ.cost_rows_from_dose_counts(r.econ.MODEL_CALENDAR_DELAYED,'refinement',vaccine or 'Baseline',key,counts,r.ECON_P))
    details=r.econ.compute_calendar_economics(r.econ.MODEL_CALENDAR_DELAYED,'refinement',vaccine or 'Baseline',key,annual,costs,r.ECON_P)
    summary=dict(vaccine=vaccine or 'No vaccination',k=k,h=h,seconds=time.perf_counter()-start,**r.summarize(result,details))
    return summary,monthly,pd.DataFrame(traces),A

def refinement():
    rows=[]
    for k in [1,2,4,8,16]:
        for v in [None]+list(r.PRODUCTS):
            row,monthly,traces,A=simulate_refined(v,k)
            rows.append(row);tag=(r.PRODUCTS[v][0] if v else 'no_vaccination')+f'_k{k}'
            monthly.to_csv(OUT/(tag+'_monthly.csv'),index=False);traces.to_csv(OUT/(tag+'_trace.csv'),index=False)
            np.savez_compressed(OUT/(tag+'_final_state.npz'),A=A)
            print('refinement',v,k,round(row['seconds'],2),flush=True)
    df=pd.DataFrame(rows)
    for k in df.k.unique():
        base=df[(df.k==k)&(df.vaccine=='No vaccination')].iloc[0]
        for ix,row in df[df.k==k].iterrows():
            for raw,key,sign in [('cases_total','cases_averted',-1),('deaths_total','deaths_averted',-1),('qaly_loss','qalys_gained',-1),('costs_total','incremental_costs',1)]:
                df.loc[ix,key+'_per100k']=sign*(row[raw]-base[raw])*1e5/row.births_denominator
    for v in df.vaccine.unique():
        ref=df[(df.vaccine==v)&(df.k==16)].iloc[0]
        for metric in ['cases_total','deaths_total','qaly_loss','costs_total','qalys_gained_per100k','incremental_costs_per100k']:
            mask=df.vaccine==v
            df.loc[mask,metric+'_relative_error_to_k16']=abs(df.loc[mask,metric]-ref[metric])/max(abs(ref[metric]),1e-12)
    df.to_csv(OUT/'refinement.csv',index=False)
    old=pd.read_csv(ROOT/'inputs'/'Final'/'outputs'/'factorial.csv')
    for row in df[df.k==1].itertuples():
        name='no_vaccination' if row.vaccine=='No vaccination' else r.scenario_id(row.vaccine,1,1,1)
        original=old[old.scenario==name].iloc[0]
        for metric in ['cases_total','deaths_total','qaly_loss','costs_total','doses_total']:
            check('k1_revised_solver_bridge',(getattr(row,metric)-original[metric])/max(1.,abs(original[metric])),1e-9,name+'/'+metric)

def mortality():
    p=r.BASE_P;neon=p['neonatal_mortality_rate']/1000 if 'neonatal_mortality_rate' in p else None
    # Read actual functions rather than guessing units from field names.
    surv=w6.base.neonatal_survival_multiplier(p)
    bg=w6.base.get_background_prob_u5_only(p);all_u5=1-(1-bg)**60
    conditional=1-(1-all_u5)/surv
    # Current accounting convention with historical double-counting comparators.
    rows=[]
    for label in ['cause_deleted','historical_inherited','historical_neonatal_conditioned','zero_background']:
        for v in [None]+list(r.PRODUCTS):
            row,_,_,_=simulate_refined(v,1,mortality_mode=label)
            row.update(background=label,status='current accounting convention' if label=='cause_deleted' else 'structural comparator, not current baseline')
            rows.append(row)
    df=pd.DataFrame(rows)
    for label in df.background.unique():
        nv=df[(df.background==label)&(df.vaccine=='No vaccination')].iloc[0]
        for ix,row in df[df.background==label].iterrows():
            df.loc[ix,'qalys_gained_per100k']=(nv.qaly_loss-row.qaly_loss)*1e5/row.births_denominator
            df.loc[ix,'incremental_costs_per100k']=(row.costs_total-nv.costs_total)*1e5/row.births_denominator
    df.to_csv(OUT/'mortality_brackets.csv',index=False)
    (OUT/'mortality.json').write_text(json.dumps(dict(neonatal_survival=surv,inherited_background_60month_probability=all_u5,
        combined_neonatal_and_background=1-surv*(1-bg)**60,conditioned_residual_probability=conditional,
        revised=r.MORTALITY_ACCOUNTING,
        status='Current: cause-deleted allocation convention. Historical comparators retain double-counting and are not recommended baselines.'),indent=2))

def historical():
    base=ROOT/'inputs/Cohort/historical'
    reps=pd.read_csv(base/'microsim_replicate_summaries.csv');old=pd.read_csv(base/'ode_microsim_reconciliation.csv')
    for x in old.itertuples():
        values=reps[(reps.aom_mode==x.aom_mode)&(reps.scenario==x.scenario)][x.metric]
        check('historical_MC_mean',values.mean()-x.microsim_mean,1e-6,x.metric)
        check('historical_MC_SD',values.std(ddof=1)-x.microsim_sd,1e-6,x.metric)
        check('historical_MC_CI',values.mean()-1.96*values.std(ddof=1)/np.sqrt(len(values))-x.microsim_95ci_low,1e-6,x.metric)
    old.to_csv(OUT/'historical_reconciliation.csv',index=False)
    # Audit retained model inputs without executing historical entry points.
    syntax=[]
    paths=list((ROOT/'inputs/Ablation/inputs/project').rglob('*.py'))+list((ROOT/'inputs/Cohort').rglob('*.py'))
    for path in sorted(paths):
        if any(x in path.parts for x in ['__pycache__','.venv']):continue
        try:ast.parse(path.read_text(encoding='utf-8-sig'));status='PASS';detail=''
        except Exception as exc:status='FAIL';detail=str(exc)
        syntax.append(dict(path=str(path.relative_to(ROOT)),status=status,detail=detail,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    pd.DataFrame(syntax).to_csv(OUT/'syntax.csv',index=False)

def aom_audit():
    rows=[]
    for frac in [.303,.301,.411,0.,1.]:
        for v in [None]+list(r.PRODUCTS):
            row,monthly,_,_=simulate_refined(v,1,aom_fraction=frac)
            row.update(aom_fraction=frac,aom_events=float(monthly.New_AOM.sum()))
            rows.append(row)
    df=pd.DataFrame(rows)
    for frac in df.aom_fraction.unique():
        nv=df[(df.aom_fraction==frac)&(df.vaccine=='No vaccination')].iloc[0]
        for ix,row in df[df.aom_fraction==frac].iterrows():
            df.loc[ix,'qalys_gained_per100k']=(nv.qaly_loss-row.qaly_loss)*1e5/row.births_denominator
            df.loc[ix,'incremental_costs_per100k']=(row.costs_total-nv.costs_total)*1e5/row.births_denominator
    check('zero_AOM_attribution',df[df.aom_fraction==0].aom_events.abs().max(),1e-12)
    df.to_csv(OUT/'aom.csv',index=False)

def main():
    start=time.time();math_checks();aggregation();manufactured();historical();refinement();mortality();aom_audit()
    pd.DataFrame(checks).to_csv(OUT/'checks.csv',index=False)
    assert all(x['passed'] for x in checks),pd.DataFrame(checks).query('not passed').to_string()
    (OUT/'provenance.json').write_text(json.dumps(dict(seed=20261004,elapsed=time.time()-start,checks=len(checks),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=sys.version),indent=2))
    print('PASS',len(checks),'checks',flush=True)

if __name__=='__main__':main()
