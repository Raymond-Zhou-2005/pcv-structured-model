"""Distributional audit of the constructed five-state CTMC (not a PCV risk model).

Reuses the original seed and 1,000 independent terminal multinomial samples at
each time. Samples at different times are independent, not linked trajectories.
Prediction intervals concern terminal counts at fixed rates and population size;
they are neither confidence intervals for the mean nor parameter uncertainty.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.linalg import expm
from scipy.stats import binom

ROOT=Path(__file__).resolve().parent
INPUT=ROOT/'inputs'/'ctmc.csv'
OUT=ROOT/'outputs'
N=100_000
REPLICATES=1_000
SEED=20261004
Q=np.array([[-.21,.20,0,0,.01],
            [0,-.162,.15,0,.012],
            [0,0,-.098,.08,.018],
            [0,0,0,-.022,.022],
            [0,0,0,0,0.]])

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    OUT.mkdir(exist_ok=True)
    original=pd.read_csv(INPUT)
    assert set(original.columns)>={'month','state','exact','ode','simulated_mean','mc_se'}
    assert len(original)==20
    p0=np.array([1.,0,0,0,0])
    rows=[];risks=[];raw=[];checks=[]
    def check(name,passed,detail):
        checks.append(dict(check=name,passed=bool(passed),detail=detail))
    for t in sorted(original.month.unique()):
        exact=expm(Q.T*t)@p0
        seed=SEED+int(t*100)
        draws=np.random.default_rng(seed).multinomial(N,exact,size=REPLICATES)
        check('terminal_probability_mass',abs(exact.sum()-1)<1e-12,float(t))
        check('nonnegative_terminal_probabilities',exact.min()>=0,float(t))
        for state,p in enumerate(exact):
            saved=original[(original.month==t)&(original.state==state)].iloc[0]
            mean=N*p;variance=N*p*(1-p);sd=np.sqrt(variance)
            sample_mean=float(draws[:,state].mean())
            sample_sd=float(draws[:,state].std(ddof=1))
            mean_se=np.sqrt(variance/REPLICATES)
            # Exact fourth central moment of Binomial(N,p), and variance of
            # its unbiased sample variance for REPLICATES independent samples.
            mu4=3*variance**2+variance*(1-6*p*(1-p))
            sample_variance_se=np.sqrt((mu4-(REPLICATES-3)/(REPLICATES-1)*variance**2)/REPLICATES)
            mean_z=(sample_mean-mean)/mean_se if mean_se else 0.
            variance_z=(sample_sd**2-variance)/sample_variance_se if sample_variance_se else 0.
            key=dict(month=float(t),state=state)
            check('matched_forward_ODE_probability',abs(saved.ode-p)<1e-10,key)
            check('matched_frozen_exact_probability',abs(saved.exact-p)<1e-14,key)
            check('reproduced_original_MC_mean',abs(saved.simulated_mean-sample_mean/N)<1e-13,key)
            check('reproduced_original_MC_standard_error',abs(saved.mc_se-sample_sd/N/np.sqrt(REPLICATES))<1e-13,key)
            check('MC_mean_within_five_analytic_SE',abs(mean_z)<5,key)
            check('MC_variance_within_five_analytic_SE',abs(variance_z)<5,key)
            rows.append(dict(**key,state_probability=float(p),expected_count=float(mean),
                             exact_count_SD=float(sd),MC_mean_count=sample_mean,MC_count_SD=sample_sd,
                             mean_z=float(mean_z),sample_variance_z=float(variance_z),seed=seed,
                             population=N,replicates=REPLICATES))
        p=float(exact[4]);counts=draws[:,4]
        low,high=binom.ppf([.025,.975],N,p)
        empirical_low,empirical_high=np.quantile(counts,[.025,.975])
        risks.append(dict(month=float(t),absorbing_event_probability=p,
                          expected_absorbing_count=N*p,exact_count_SD=float(np.sqrt(N*p*(1-p))),
                          MC_count_SD=float(counts.std(ddof=1)),
                          count_prediction_2p5_percent=float(low),count_prediction_97p5_percent=float(high),
                          empirical_count_2p5_percent=float(empirical_low),empirical_count_97p5_percent=float(empirical_high),
                          population=N,replicates=REPLICATES,seed=seed))
        raw.extend(dict(month=float(t),replicate=i+1,absorbing_count=int(count),seed=seed)
                   for i,count in enumerate(counts))
    pd.DataFrame(rows).to_csv(OUT/'ctmc_variability.csv',index=False)
    pd.DataFrame(risks).to_csv(OUT/'ctmc_absorbing_event.csv',index=False)
    pd.DataFrame(raw).to_csv(OUT/'ctmc_absorbing_replicates.csv',index=False)
    summary=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',
                 checks=checks,check_count=len(checks),passed=sum(c['passed'] for c in checks),
                 input_sha256=sha(INPUT),script_sha256=sha(Path(__file__)),
                 python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
                 rates='Constructed fixed five-state CTMC; not calibrated PCV disease hazards',
                 uncertainty='Conditional sampling variation; no parameter PSA',
                 intervals='Central binomial terminal-count prediction intervals at fixed inputs',
                 seed_rule='20261004 + int(month*100); 1000 terminal draws at each time',
                 longitudinal_linkage=False,at_12_months=risks[-1])
    (OUT/'ctmc_variability_checks.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k!='checks'},indent=2))
    assert all(c['passed'] for c in checks),[c for c in checks if not c['passed']]

if __name__=='__main__':main()
