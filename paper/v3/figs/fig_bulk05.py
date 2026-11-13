"""V2: unfiltered late-window densities and explicitly censored decay times."""
import csv
import os
from pathlib import Path
from common import ROOT, PAGE, label, save, project_path
import matplotlib.pyplot as plt
import numpy as np
from campaign.v2.analysis.io import read_timeseries

from common import QSTAR05 as QREF
RESULTS=Path(ROOT)/'campaign/v2/results/bulk'
if not (RESULTS/'summary.csv').exists():
    raise SystemExit('Run python3 -m campaign.v2.analysis.bulk before making the v2 figure.')
rows=list(csv.DictReader((RESULTS/'summary.csv').open()))
for r in rows:
    for k,v in list(r.items()):
        try:r[k]=float(v)
        except ValueError:pass
fig,axs=plt.subplots(1,3,figsize=(PAGE,2.35),gridspec_kw=dict(wspace=.35))
qs=sorted({r['q'] for r in rows});colors={q:plt.get_cmap('coolwarm')(i/(len(qs)-1)) for i,q in enumerate(qs)}
runrows=list(csv.DictReader((RESULTS/'runs.csv').open()))
for q in qs:
    choices=[r for r in runrows if int(r['Lx'])==2048 and float(r['q'])==q]
    if choices:
        _,ts=read_timeseries(project_path(choices[0]['path']));ok=ts[:,0]>=10
        axs[0].plot(ts[ok,0],ts[ok,2],color=colors[q],lw=.8,label=f'${q/QREF:.3f}$')
axs[0].set(xscale='log',xlabel='$t$',ylabel=r'interface density $\rho$',xlim=(10,3e5),ylim=(0,.42))
axs[0].legend(fontsize=5,ncol=2,title=r'$q/q^*$ ($L=2048$)',title_fontsize=5.5,loc='lower left',handlelength=1.2)
marks={256:'v',512:'D',1024:'o',2048:'s',4096:'^'}
for L,mk in marks.items():
    rr=sorted((r for r in rows if r['Lx']==L),key=lambda r:r['q'])
    axs[1].errorbar([r['q']/QREF for r in rr],[r['rho_all_mean'] for r in rr],
                    [r['rho_all_sem'] for r in rr],fmt=mk,ms=3,mfc='none',capsize=1,label=f'$L={L}$')
axs[1].axvline(1,color='.6',lw=.5,ls=':')
axs[1].set(xlabel=r'$q/q^*$',ylabel=r'late-window mean $\rho$',xlim=(.885,1.11),ylim=(-.01,.37))
axs[1].legend(fontsize=5.5,loc='lower right')
axs[1].text(.03,.88,'all completed runs',transform=axs[1].transAxes,fontsize=5.5)
for q in qs:
    if q>=QREF:continue
    rr=sorted((r for r in rows if r['q']==q and r['half_events']==r['n']),key=lambda r:r['Lx'])
    if rr:
        val=np.array([r['half_rmst'] for r in rr]);lo=np.array([r['half_rmst_low'] for r in rr]);hi=np.array([r['half_rmst_high'] for r in rr])
        axs[2].errorbar([r['Lx'] for r in rr],val,yerr=[val-lo,hi-val],fmt='o-',color=colors[q],ms=3,lw=.8,capsize=1)
axs[2].set(xscale='log',yscale='log',xlabel='$L$',ylabel=r'mean half-density time $\tau$',xlim=(45,6000))
axs[2].set_xticks([64,256,1024,4096],['64','256','1024','4096'])
axs[2].text(.96,.05,'95% run bootstrap; upper time brackets',transform=axs[2].transAxes,ha='right',fontsize=5)
for ax,s in zip(axs,'abc'):label(ax,f'({s})',x=.04 if s!='a' else .86)
save(fig,'fig_bulk05')
