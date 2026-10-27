"""Reproduce revised Figures 1, 2 and 6 from frozen, unrounded plot data.

Run: python reproduction/plot_figures.py
No model coefficients or simulation outputs are changed by this script.
"""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
DATA=Path(__file__).resolve().parent/'data'
FIG=ROOT/'figures'
PRODUCTS=['PCV10-GSK','PCV13','PCV10-SII']
COLORS=['#0072B2','#009E73','#D55E00']
CHECKS=[]
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,
    'axes.titlesize':10,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,
    'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
    'savefig.facecolor':'white'})

def check(name,passed,detail=''):
    CHECKS.append(dict(check=name,passed=bool(passed),detail=detail))
    if not passed:raise AssertionError((name,detail))

def save(fig,number):
    FIG.mkdir(exist_ok=True)
    for suffix,kwargs in [('pdf',{}),('png',{'dpi':200}),('tif',{'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}})]:
        fig.savefig(FIG/f'Fig{number}.{suffix}',bbox_inches='tight',**kwargs)
    with Image.open(FIG/f'Fig{number}.tif') as im:
        im.convert('RGB').save(FIG/f'Fig{number}.tif',compression='tiff_lzw',dpi=(600,600))
    plt.close(fig)

def figure1():
    c=pd.read_csv(DATA/'ctmc.csv');m=pd.read_csv(DATA/'manufactured.csv')
    r=pd.read_csv(DATA/'refinement.csv')
    check('CTMC_ODE_matches_exponential',np.max(abs(c.exact-c.ode))<1e-10)
    fig,axs=plt.subplots(1,3,figsize=(9.0,3.25),layout='constrained')
    state_colors=['#0072B2','#E69F00','#009E73','#CC79A7','#666666']
    for state,g in c.groupby('state'):
        axs[0].plot(g.month,g.exact,'o-',ms=3,lw=1.3,color=state_colors[int(state)],label=f'State {state}')
        axs[0].plot(g.month,g.ode,'x',ms=4,mew=.8,color='black')
    handles,labels=axs[0].get_legend_handles_labels()
    handles.append(Line2D([],[],color='black',marker='x',ls='none',ms=4));labels.append('ODE')
    axs[0].legend(handles,labels,ncol=2,frameon=False,fontsize=7.3,loc='upper right')
    axs[0].set(xlabel='Time (months)',ylabel='State probability',title='A  Matched CTMC and ODE',ylim=(-.025,1.025))
    axs[1].loglog(m.h,m.L1_density_error,'o-',color='#0072B2',ms=4,label='Cell-mass $L^1$ error')
    axs[1].loglog(m.h,m.total_error,'s-',color='#D55E00',ms=4,label='Population error')
    axs[1].set(xlabel='Step $h$ (months)',ylabel='Absolute error (children)',title='B  Manufactured transport')
    axs[1].legend(frameon=False,loc='upper left',fontsize=7.5)
    ticks=[1/16,1/8,1/4,1/2,1]
    for vaccine,color in zip(PRODUCTS,COLORS):
        q=r[r.vaccine==vaccine].sort_values('h')
        check('five_grids_'+vaccine,list(q.k)==[16,8,4,2,1])
        ref=float(q[q.k==16].qalys_gained_per100k.iloc[0])
        expected=abs(q.qalys_gained_per100k-ref)/abs(ref)
        check('refinement_error_identity_'+vaccine,np.max(abs(expected-q.qalys_gained_per100k_relative_error_to_k16))<1e-12)
        axs[2].plot(q.h,q.qalys_gained_per100k,'o-',color=color,lw=1.5,ms=4,label=vaccine)
    axs[2].set_xscale('log',base=2)
    axs[2].set(xticks=ticks,xticklabels=['1/16','1/8','1/4','1/2','1'],xlabel='Step $h$ (months)',
        ylabel='QALYs gained / 100,000 births',title='C  Programme QALY gains',ylim=(1395,1630))
    axs[2].legend(frameon=False,loc='center',bbox_to_anchor=(.48,.4),fontsize=7.7)
    for ax in axs:ax.grid(axis='y',lw=.35,color='#dddddd')
    r[(r.vaccine.isin(PRODUCTS))&(r.k==16)][['vaccine','h','cases_averted_per100k','deaths_averted_per100k','qalys_gained_per100k']].to_csv(DATA/'finest_grid_results.csv',index=False)
    save(fig,1)

def figure2():
    a=np.load(DATA/'aggregation_states.npz')['states'];df=pd.read_csv(DATA/'aggregation.csv')
    p=json.loads((DATA/'parameters.json').read_text())
    check('aggregation_same_dose_marginals',np.array_equal(a[0].sum(1),a[1].sum(1)))
    check('aggregation_same_age_marginals',np.array_equal(a[0].sum(0),a[1].sum(0)))
    rho=np.array([0,0,.15,.34,.4]);w=np.repeat([0,1,.8,.35],[6,12,6,36])
    for i in range(2):
        check('aggregation_total_'+str(i),float(a[i].sum())==50000)
        check('aggregation_each_dose_total_'+str(i),np.all(a[i].sum(1)==10000))
        age=np.argmax(a[i],axis=1)
        expected=[20,4,8,16,48] if i==0 else [20,4,8,48,16]
        check('aggregation_exact_ages_'+str(i),list(age)==expected,expected)
        c=float(np.sum(a[i]*rho[:,None]*w)/np.sum(a[i]*w))
        rel=float(np.clip(c/(p['vaccine_coverage']*p['mature_effective_dose_share']*.34),0,1))
        target=p['vaccine_irrs']['PCV13'][5]
        hstar=float(np.clip(target/(1-p['vaccine_coverage']*p['mature_effective_dose_share']*p['vaccine_direct_effectiveness']['PCV13']['ipd']),0,1))
        h=1-(1-hstar)*rel*p['carriage_to_herd_fraction']
        check('aggregation_C_'+str(i),abs(c-df.C.iloc[i])<1e-12)
        check('aggregation_r_'+str(i),abs(rel-df.r.iloc[i])<1e-12)
        check('aggregation_H_'+str(i),abs(h-df.H_ipd.iloc[i])<1e-12)
    fig,axs=plt.subplots(1,3,figsize=(9.0,3.35),layout='constrained',gridspec_kw={'width_ratios':[1.05,1.05,1.]})
    labels=['Non-programme','No dose','One dose','Two doses','Booster']
    for i,ax in enumerate(axs[:2]):
        ages=np.argmax(a[i],axis=1)
        ax.scatter(ages,np.arange(5),s=65,color=[COLORS[0],COLORS[2]][i],zorder=3)
        for g,age in enumerate(ages):ax.annotate(str(age),(age,g),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
        ax.set(xlim=(-2,61),ylim=(-.7,4.8),xlabel='Age (months)',title=f'{chr(65+i)}  Configuration {chr(65+i)}',
            xticks=[0,12,24,36,48,60],yticks=range(5),yticklabels=labels if i==0 else ['']*5)
        ax.grid(axis='y',lw=.45,color='#dddddd')
    axs[0].text(.02,.98,'Each marker: 10,000 children',transform=axs[0].transAxes,ha='left',va='top',fontsize=7.7)
    axs[2].bar([0,1],df.one_month_cases,color=[COLORS[0],COLORS[2]],width=.58)
    for i,row in df.iterrows():
        axs[2].text(i,row.one_month_cases+2,f'{row.one_month_cases:.2f}\n$H_{{IPD}}$ = {row.H_ipd:.4f}',ha='center',va='bottom',fontsize=8)
    axs[2].set(xticks=[0,1],xticklabels=['Configuration A','Configuration B'],ylim=(0,126),
        ylabel='Expected one-month disease episodes',title='C  Same totals, different outcomes')
    delta=100*(df.one_month_cases.iloc[1]/df.one_month_cases.iloc[0]-1)
    check('aggregation_reported_case_decrement',round(delta,2)==-7.36,delta)
    axs[2].text(.5,.96,f'B relative to A: {delta:.2f}%',transform=axs[2].transAxes,ha='center',va='top',fontsize=8)
    save(fig,2)

def figure6():
    s=pd.read_csv(DATA/'sensitivity.csv');f=pd.read_csv(DATA/'factorial.csv')
    aom=pd.read_csv(DATA/'aom.csv');mort=pd.read_csv(DATA/'mortality_brackets.csv')
    params=list(dict.fromkeys(zip(s.parameter,s.value)))
    check('all_sensitivity_settings',len(params)==20,len(params))
    short={'price':'Price','direct_VE':'Direct VE','CFR':'CFR','sequela_cost':'Sequela cost',
        'utility_loss':'Utility loss','burden':'Disease burden','waning':'Waning severity',
        'herd_strength':'Herd strength','delay_probability':'Dose opportunity',
        'birth_seasonality':'Uniform births','healthcare_perspective':'Restricted healthcare'}
    labels=[short[p]+('' if p in ['birth_seasonality','healthcare_perspective'] else f' ×{v:g}') for p,v in params]
    fig=plt.figure(figsize=(9.2,6.0),layout='constrained')
    gs=fig.add_gridspec(2,2,height_ratios=[3.0,1.35],wspace=.07,hspace=.08)
    axs=[fig.add_subplot(gs[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
    for ax,metric,title in zip(axs[:2],['qalys_gained_per100k','incremental_costs_per100k'],['A  QALY-gain sensitivity','B  Net-cost sensitivity']):
        for k,(v,col,marker) in enumerate(zip(PRODUCTS,COLORS,['o','s','D'])):
            base=float(f[(f.vaccine==v)&(f.D==1)&(f.W==1)&(f.G==1)&(~f.herd_off)][metric].iloc[0])
            ys=[]
            for param,value in params:
                row=s[(s.vaccine==v)&(s.parameter==param)&(s.value==value)].iloc[0]
                delta=float(row[metric]-base)
                check('sensitivity_delta_'+v+'_'+metric+'_'+param+'_'+str(value),abs(delta-float(row['delta_'+metric]))<1e-7)
                ys.append(100*delta/abs(base))
            ax.scatter(ys,np.arange(len(params))+.18*(k-1),s=18,marker=marker,color=col,label=v,zorder=3)
        ax.axvline(0,color='#777777',lw=.7)
        ax.set(yticks=range(len(params)),yticklabels=labels,ylim=(len(params)-.4,-.7),title=title,
            xlabel='Change / |full111 baseline| (%)')
        ax.tick_params(axis='y',labelsize=7.8)
        ax.grid(axis='y',color='#e7e7e7',lw=.4)
    q=aom[aom.vaccine=='PCV13'].sort_values('aom_fraction')
    axs[2].scatter(q.incremental_costs_per100k/1e6,range(len(q)),color=COLORS[1],s=27,zorder=3)
    for i,row in enumerate(q.itertuples()):
        axs[2].annotate(f'{row.incremental_costs_per100k/1e6:.3f}',(row.incremental_costs_per100k/1e6,i),
            xytext=(6,0),textcoords='offset points',va='center',fontsize=7.8)
    axs[2].set(yticks=range(len(q)),yticklabels=[f'θ = {x:g}' for x in q.aom_fraction],
        ylim=(len(q)-.4,-.6),xlim=(-6.9,-2.1),title='C  AOM attribution: PCV13',
        xlabel='Incremental cost (million scenario USD / 100,000 births)')
    names={'cause_deleted':'Cause-deleted','historical_inherited':'Inherited',
        'historical_neonatal_conditioned':'Neonatal-conditioned','zero_background':'Zero background'}
    q=mort[mort.vaccine=='PCV13']
    corrected=float(q[q.background=='cause_deleted'].qalys_gained_per100k.iloc[0])
    axs[3].axvline(corrected,color='#999999',lw=.7,ls='--')
    axs[3].scatter(q.qalys_gained_per100k,range(len(q)),s=27,color=COLORS[1],zorder=3)
    for i,row in enumerate(q.itertuples()):
        axs[3].annotate(f'{row.qalys_gained_per100k:,.2f}',(row.qalys_gained_per100k,i),
            xytext=(6,0),textcoords='offset points',va='center',fontsize=7.8,
            bbox={'facecolor':'white','edgecolor':'none','pad':.3})
    axs[3].set(yticks=range(len(q)),yticklabels=[names[x] for x in q.background],ylim=(len(q)-.5,-.5),
        xlim=(1515,1532),title='D  Mortality accounting: PCV13',xlabel='QALYs gained / 100,000 births')
    for ax in axs[2:]:ax.grid(axis='y',color='#e7e7e7',lw=.4)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False,fontsize=8.5)
    save(fig,6)

def main():
    figure1();figure2();figure6()
    pd.DataFrame(CHECKS).to_csv(DATA.parent/'figure_checks.csv',index=False)
    manifest={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(DATA.glob('*')) if p.is_file()}
    (DATA.parent/'figure_data_hashes.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('PASS',len(CHECKS),'figure/data checks')

if __name__=='__main__':main()
