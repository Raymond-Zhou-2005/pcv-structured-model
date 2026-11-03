"""Draw main Figures 3–5 and nine extended state plots without editing text."""
from pathlib import Path
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parent
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
APP=ROOT/'research/inputs/Final/outputs'
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,
 'legend.fontsize':8.3,'xtick.labelsize':8.3,'ytick.labelsize':8.3,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'ps.fonttype':42})
def save(fig,n):
    fig.savefig(FIG/f'Fig{n}.pdf')
    fig.savefig(FIG/f'Fig{n}.png',dpi=200)
    fig.savefig(FIG/f'Fig{n}.tif',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    with Image.open(FIG/f'Fig{n}.tif') as im:
        im.convert('RGB').save(FIG/f'Fig{n}.tif',compression='tiff_lzw',dpi=(600,600))
    plt.close(fig)
def main(only_main=False):
    fig,axs=plt.subplots(3,3,figsize=(7.4,7.7),layout='constrained')
    names=['no_vaccination','pcv13_D0W1G1','pcv13_D1W1G1'];labels=['No vaccination','On time (011)','Delayed / incomplete (111)'];linecols=['#555555','#0072B2','#D55E00']
    for i,(name,label,col) in enumerate(zip(names,labels,linecols)):
     a=np.load(APP/'runs'/name/'states.npz')['alive_after_births_and_doses'];tr=pd.read_csv(APP/'runs'/name/'herd.csv');mo=pd.read_csv(APP/'runs'/name/'monthly.csv')
     density_image=axs[0,i].imshow(a.sum(1).T/1000,origin='lower',aspect='auto',extent=[0,10,0,60],cmap='Blues',vmin=0,vmax=130); axs[0,i].set(title=f'{chr(65+i)}   '+label.replace('Delayed / incomplete (111)','Delayed / incomplete\n(111)'),xlabel='Programme year',ylabel='Age (months)')
     axs[1,0].plot(np.arange(120)/12,a.sum((1,2))/1e6,label=label,color=col)
     axs[1,1].plot(np.arange(120)/12,tr.relative_carriage_index,color=col,label=label+' r')
     axs[1,1].plot(np.arange(120)/12,tr.carriage_index,color=col,ls='--')
     axs[1,2].plot(np.arange(120)/12,tr.H_ipd,color=col)
     annual=mo.groupby('Year').Monthly_Total_Cases.sum();axs[2,2].plot(annual.index,annual/1000,color=col,label=label)
     if i==2:
      dose_values=[]
      for g in range(5):
       values=a[:,g].sum(1)/1e6
       line,=axs[2,0].plot(np.arange(120)/12,values,label=['Non-program','D0','D1','D2','Booster'][g])
       dose_values.append((values[-1],line.get_color(),line.get_label()))
      for position,(end,col,dose_label) in zip(np.linspace(.15,2.7,5),sorted(dose_values)):
       axs[2,0].annotate(dose_label,xy=(119/12,end),xytext=(10.35,position),fontsize=8,
         va='center',color=col,arrowprops={'arrowstyle':'-','color':col,'lw':.5})
      axs[2,0].set_xlim(0,14.4)
      axs[2,0].set_xticks([0,5,10])
      share=a[:,3:].sum(1)/np.maximum(a.sum(1),1e-12)
      share_image=axs[2,1].imshow(share.T,origin='lower',aspect='auto',extent=[0,10,0,60],cmap='viridis',vmin=0,vmax=1)
    for ax,title,y in [(axs[1,0],'D   Alive population','Children (millions)'),(axs[1,1],'E   C (dashed)\nand r (solid)','Dimensionless index'),(axs[1,2],'F   IPD herd multiplier','H(t)'),(axs[2,0],'G   Delayed dose-state\nintegrals','Children (millions)'),(axs[2,1],'H   Delayed: at least\ntwo doses','Age (months)'),(axs[2,2],'I   Annual disease\nepisodes','Episodes (thousands)')]:ax.set(title=title,xlabel='Programme year',ylabel=y)
    fig.colorbar(density_image,ax=list(axs[0]),label='Children per month-age bin\n(thousands)',fraction=.035,pad=.015,location='right')
    fig.colorbar(share_image,ax=axs[2,1],label='Share',fraction=.045,pad=.02)
    fig.legend(*axs[1,0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False)
    save(fig,3)
    sys.path.insert(0,str(ROOT/'research/src'))
    import revised_plots as plots
    plots.application_figures(APP,FIG,save)
    if only_main:
        print('Reproduced main Figures 3–5; extended figures retained unchanged')
        return
    directory=FIG/'trajectories';directory.mkdir(exist_ok=True)
    original=plots._save_si
    def named_save(fig,figdir,stem):
        i=int(stem[1:])-1
        prefix=['pcv10_gsk','pcv13','pcv10_sii'][i%3]
        kind=['population_trajectory','dose_integrals','age_calendar_dose_density'][i//3]
        original(fig,figdir,prefix+'_'+kind)
    plots._save_si=named_save
    plots.supplementary_trajectories(APP,directory)
    print('Reproduced Figures 3–5 and nine extended state figures')
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only-main',action='store_true',help='Retain the nine extended figures unchanged.')
    main(parser.parse_args().only_main)
