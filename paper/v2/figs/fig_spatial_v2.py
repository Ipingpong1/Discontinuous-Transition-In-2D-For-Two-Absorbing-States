"""Signed spatial correlations from complete fixed final-configuration ensembles."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from common import ROOT,PAGE,label,save

path=Path(ROOT)/'campaign/v2/results/spatial_correlations/correlations_summary.csv'
def clear():
    for extension in ('pdf','png'):
        Path(__file__).with_name('fig_spatial_v2.'+extension).unlink(missing_ok=True)
if not path.exists():
    clear();raise SystemExit('Run spatial_correlations first')
rows=[]
for row in csv.DictReader(path.open()):
    rows.append({k:float(v) for k,v in row.items()})
required=[(8192,.04594),(16384,.04592),(16384,.04594),(16384,.04596)]
missing=[(L,q) for L,q in required if not any(r['L']==L and r['q']==q for r in rows)]
if missing:
    clear()
    raise SystemExit('Spatial comparison awaits complete cohorts: '+str(missing))

fig,axs=plt.subplots(1,2,figsize=(PAGE,2.6),gridspec_kw=dict(wspace=.28))
def curve(ax,L,q,color,name):
    rr=sorted((r for r in rows if r['L']==L and r['q']==q and r['r']>0),key=lambda r:r['r'])
    x=np.array([r['r'] for r in rr]);y=np.array([r['C'] for r in rr])
    ax.plot(x,y,color=color,lw=1,label=name)
    ax.fill_between(x,[r['ci_low'] for r in rr],[r['ci_high'] for r in rr],color=color,alpha=.16,lw=0)
curve(axs[0],8192,.04594,'C0','$L=8192$')
curve(axs[0],16384,.04594,'C1','$L=16384$')
for q,color in [(.04592,'C4'),(.04594,'C1'),(.04596,'C2')]:
    curve(axs[1],16384,q,color,f'$q={q:.5f}$')
for ax,s in zip(axs,'ab'):
    ax.set_xscale('log');ax.set_yscale('symlog',linthresh=.001,linscale=.6)
    ax.set(xlabel='separation $r$',ylabel=r'$C(r)=\langle s_i s_{i+r}\rangle$',
           xlim=(1,1024),ylim=(-.003,.7))
    ax.set_xticks([1,4,16,64,256,1024],['1','4','16','64','256','1024'])
    ax.set_yticks([-.001,0,.001,.01,.1],['$-10^{-3}$','0','$10^{-3}$','$10^{-2}$','$10^{-1}$'])
    ax.axhline(0,color='.55',lw=.5)
    ax.legend(fontsize=6,loc='upper right')
    label(ax,f'({s})',x=.04,y=.12)
axs[0].set_title('$q=0.04594$')
axs[1].set_title('$L=16384$')
save(fig,'fig_spatial_v2')
