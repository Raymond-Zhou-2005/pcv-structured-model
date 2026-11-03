"""Plot revised application outputs; never read preceding figures as results."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

PRODUCTS=['PCV10-GSK','PCV13','PCV10-SII']
PREFIXES=['pcv10_gsk','pcv13','pcv10sii']
COLORS=['#0072B2','#009E73','#D55E00']

def application_figures(app, figdir, save):
    f=pd.read_csv(app/'factorial.csv'); sh=pd.read_csv(app/'shapley.csv')
    herd=pd.read_csv(app/'herd_contributions.csv')
    order=[(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]
    fig,axs=plt.subplots(3,2,figsize=(7.4,7.4),layout='constrained')
    metrics=['cases_averted_per100k','deaths_averted_per100k','qalys_gained_per100k','incremental_costs_per100k']
    titles=['A   Cases averted','B   Deaths averted','C   QALYs gained','D   Incremental cost']
    for ax,metric,title in zip(axs.flat,metrics,titles):
        for v,col in zip(PRODUCTS,COLORS):
            q=f[(f.vaccine==v)&(~f.herd_off)].set_index(['D','W','G'])
            baseline=float(q.loc[(0,0,0),metric])
            y=np.array([float(q.loc[b,metric]) for b in order])
            ax.plot(range(8),100*(y-baseline)/abs(baseline),'o-',color=col,lw=1.25,ms=3,label=v)
        ax.axhline(0,color='#888888',lw=.5)
        ax.set(xticks=range(8),xticklabels=[''.join(map(str,b)) for b in order],title=title,xlabel='DWG configuration',ylabel='Change / |000| (%)')
        ax.tick_params(axis='x',labelsize=8.3)
    ax=axs[2,0]
    bottom=np.zeros(3)
    for factor,col,label in [('D','#56B4E9','Delayed / incomplete'),('W','#CC79A7','Attained-age waning'),('G','#E69F00','Gradual herd build-up')]:
        y=np.array([-float(sh[(sh.vaccine==v)&(sh.factor==factor)&(sh.metric=='qalys_gained_per100k')].shapley_delta.iloc[0]) for v in PRODUCTS])
        ax.bar(range(3),y,bottom=bottom,color=col,label=label,width=.6)
        bottom+=y
    ax.set(xticks=range(3),xticklabels=['GSK','PCV13','SII'],ylabel='QALY decrement\nper 100,000 births',title='E   Shapley allocation, 111 − 000')
    ax.legend(fontsize=8.3,frameon=False,loc='upper left')
    ax=axs[2,1]
    off=[]; indirect=[]
    for v in PRODUCTS:
        q=herd[(herd.vaccine==v)&(herd.D==1)&(herd.W==1)&(herd.G==1)].iloc[0]
        off.append(q.direct_only_qalys_gained_per100k);indirect.append(q.indirect_qalys_gained_per100k)
    ax.bar(range(3),off,color='#777777',label='Herd off',width=.6)
    ax.bar(range(3),indirect,bottom=off,color='#009E73',label='Paired indirect contribution',width=.6)
    ax.set(xticks=range(3),xticklabels=['GSK','PCV13','SII'],ylabel='QALY gain\nper 100,000 births',title='F   Full 111: paired herd contrast',ylim=(0,2350))
    ax.legend(fontsize=8.3,frameon=False,loc='upper left')
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False)
    save(fig,4)
    s=pd.read_csv(app/'schedules.csv')
    strategies=list(s.strategy.unique())
    names={'delayed':'Delayed / incomplete','dose1_only':'One dose only','completion_catchup':'Completion recovery'}
    fig,axs=plt.subplots(1,3,figsize=(7.4,4.6),layout='constrained',sharey=True)
    for ax,metric,title,scale in zip(axs,['qalys_gained_per100k','cases_averted_per100k','incremental_costs_per100k'],['A   QALYs gained','B   Cases averted','C   Incremental cost'],[1,1,1e6]):
        for i,(v,col) in enumerate(zip(PRODUCTS,COLORS)):
            y=np.array([s[(s.vaccine==v)&(s.strategy==z)][metric].iloc[0]/scale for z in strategies])
            ax.plot(y,np.arange(len(strategies))+.13*(i-1),'o-',lw=1,ms=3,color=col,label=v)
        ax.set(yticks=np.arange(len(strategies)),yticklabels=[names.get(z,z).replace('Delayed / incomplete','Delayed /\nincomplete').replace('Completion recovery','Completion\nrecovery') for z in strategies],title=title,xlabel='Per 100,000 births'+('\n(million scenario USD)' if scale==1e6 else ''))
        ax.invert_yaxis(); ax.grid(axis='x',color='#dddddd',lw=.4)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False)
    save(fig,5)

def _save_si(fig,figdir,stem):
    fig.savefig(figdir/f'{stem}.pdf',bbox_inches='tight')
    fig.savefig(figdir/f'{stem}.png',dpi=160,bbox_inches='tight')
    fig.savefig(figdir/f'{stem}.tif',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
    with Image.open(figdir/f'{stem}.tif') as im:
        im.convert('RGB').save(figdir/f'{stem}.tif',compression='tiff_lzw',dpi=(600,600))
    plt.close(fig)

def supplementary_trajectories(app,figdir):
    state_labels=['Non-programme','No dose','One dose','Two doses','Booster']
    for j,(v,pre) in enumerate(zip(PRODUCTS,PREFIXES)):
        runs=['no_vaccination',pre+'_D0W1G1',pre+'_D1W1G1']
        labels=['No vaccination','On time (011)','Delayed / incomplete (111)']
        cols=['#555555','#0072B2','#D55E00']; data=[]
        for run in runs:
            a=np.load(app/'runs'/run/'states.npz')['alive_after_births_and_doses']
            h=pd.read_csv(app/'runs'/run/'herd.csv'); m=pd.read_csv(app/'runs'/run/'monthly.csv')
            data.append((a,h,m))
        fig,axs=plt.subplots(2,3,figsize=(7.2,5.5),layout='constrained')
        for (a,h,m),label,col in zip(data,labels,cols):
            t=np.arange(len(a))/12
            vals=[a.sum((1,2))/1e6,h.carriage_index,h.relative_carriage_index,h.H_ipd,h.H_pneumonia,m.Monthly_Total_Cases/1000]
            for ax,y in zip(axs.flat,vals):ax.plot(t,y,color=col,label=label,lw=1.2)
        for ax,title,yl in zip(axs.flat,['A   Alive population','B   C: weighted dose proxy','C   r: target-normalised proxy','D   H: IPD multiplier','E   H: pneumonia multiplier','F   Monthly disease episodes'],['Millions of children','C','r','H','H','Thousands of episodes']):ax.set(title=title,xlabel='Programme year',ylabel=yl)
        fig.suptitle(v,fontsize=10);fig.legend(*axs[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False)
        _save_si(fig,figdir,f'S{j+1}')
        fig,axs=plt.subplots(1,3,figsize=(7.2,3.4),layout='constrained',sharey=True)
        for ax,(a,h,m),label in zip(axs,data,labels):
            for g in range(5):ax.plot(np.arange(len(a))/12,a[:,g].sum(1)/1e6,label=state_labels[g],lw=1.1)
            ax.set(title=label,xlabel='Programme year',ylabel='Millions of children')
        fig.suptitle(v+': age-integrated dose-state populations',fontsize=10)
        fig.legend(*axs[0].get_legend_handles_labels(),loc='outside lower center',ncol=5,frameon=False,fontsize=7)
        _save_si(fig,figdir,f'S{j+4}')
        fig,axs=plt.subplots(3,5,figsize=(7.2,6.2),layout='constrained',sharex=True,sharey=True)
        vmax=max(a.max() for a,_,_ in data)/1000
        for i,((a,h,m),label) in enumerate(zip(data,labels)):
            for g in range(5):
                im=axs[i,g].imshow(a[:,g].T/1000,origin='lower',aspect='auto',extent=[0,10,0,60],cmap='Blues',vmin=0,vmax=vmax,rasterized=True)
                if i==0:axs[i,g].set_title(state_labels[g],fontsize=8)
                if g==0:axs[i,g].set_ylabel(label+'\nAge (months)',fontsize=7)
                if i==2:axs[i,g].set_xlabel('Programme year',fontsize=7)
        fig.suptitle(v+': joint age–calendar density by programme state',fontsize=10)
        fig.colorbar(im,ax=axs.ravel().tolist(),label='Thousand children per month-age bin',fraction=.018,pad=.012)
        _save_si(fig,figdir,f'S{j+7}')
