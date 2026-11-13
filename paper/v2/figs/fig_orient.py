"""V2 orientation figure from true-envelope reductions, with local crossing intervals."""
import csv
import json
from pathlib import Path
from common import ROOT, PAGE, label, save
import matplotlib.pyplot as plt
import numpy as np

RESULTS=Path(ROOT)/'campaign/v2/results/fronts'
if not (RESULTS/'analysis.json').exists():
    raise SystemExit('Run python3 -m campaign.v2.analysis.fronts on the raw profiles first.')
result=json.loads((RESULTS/'analysis.json').read_text())
points=list(csv.DictReader((RESULTS/'points.csv').open()))
for r in points:
    for k in r:
        if k not in ('ab','estimator'):r[k]=float(r[k])
AB={'1-0':('o','C0',r'$(1,0)$'),'2-1':('s','C1',r'$(2,1)$'),'1-1':('^','C3',r'$(1,1)$')}
def selected(r):
    return r['L']==4096 and r['threshold']==.3 and r['t0']==20000 and r['estimator']=='envelope'
fig,axs=plt.subplots(1,3,figsize=(PAGE,2.25),gridspec_kw=dict(wspace=.32))
for ax,p in zip(axs,[.1,.5,.9]):
    for ab,(mk,color,lab) in AB.items():
        rr=sorted((r for r in points if selected(r) and r['p']==p and r['ab']==ab),key=lambda r:r['q'])
        q=np.array([r['q'] for r in rr]);v=np.array([r['v'] for r in rr]);se=np.array([r['sem'] for r in rr])
        ax.errorbar(q,v*1e3,se*1e3,fmt=mk,color=color,ms=3,capsize=1.2,label=lab)
        fs=[f for f in result['fits'] if selected(f) and f['p']==p and f['ab']==ab and f['model']=='bracket']
        if fs:
            f=fs[0];xx=np.linspace(f['q_min'],f['q_max'],25)
            yy=np.polynomial.polynomial.polyval((xx-f['center'])/f['scale'],f['coefficients'])
            ax.plot(xx,yy*1e3,color=color,lw=.8)
            if f['ci_low'] is not None and f['ci_high'] is not None:
                ax.axvspan(f['ci_low'],f['ci_high'],color=color,alpha=.10,lw=0)
            print(f"p={p} {ab}: q*={f['crossing']:.9g}, 95% [{f['ci_low']}, {f['ci_high']}]")
    ax.axhline(0,color='.6',lw=.5)
    ax.set(xlabel='$q$',title=f'$p={p}$')
    ax.ticklabel_format(axis='x',style='sci',scilimits=(-2,2),useOffset=False)
    if p==.5:ax.set_ylim(-4.5,3.6)
axs[0].set_ylabel(r'normal velocity $v\ (10^{-3})$')
axs[0].legend(fontsize=6,loc='upper left')
for ax,s in zip(axs,'abc'):label(ax,f'({s})',x=.86,y=.1)
save(fig,'fig_orient')
