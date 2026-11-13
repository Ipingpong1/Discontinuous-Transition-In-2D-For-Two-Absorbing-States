"""CPU fixed-cohort survival; plots only wholly completed manifest groups."""
import csv,os
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from common import ROOT,PAGE,label,save
path=os.path.join(ROOT,'campaign/v2/results/cpu_lifetimes')
with open(os.path.join(path,'summary_complete.csv')) as f:
    summary=[r for r in csv.DictReader(f) if abs(float(r['q'])-.046359)<1e-10]
if len(summary)!=4:raise SystemExit('Await all four fixed32-seed cohorts at q=.046359; no partial replacementfigure.')
with open(os.path.join(path,'survival_complete.csv')) as f:
    curves=list(csv.DictReader(f))
fig,axes=plt.subplots(1,2,figsize=(PAGE,2.55),gridspec_kw={'wspace':.25})
for ax,endpoint in zip(axes,['half','absorption']):
    for L,color in [(64,'C0'),(128,'C1'),(192,'C2'),(256,'C3')]:
        rows=sorted([r for r in curves if int(r['Lx'])==L and abs(float(r['q'])-.046359)<1e-10 and r['endpoint']==endpoint],key=lambda r:float(r['time']))
        t=np.array([max(1.,float(r['time'])) for r in rows]);s=np.array([float(r['survival']) for r in rows])
        lo=np.array([float(r['low']) for r in rows]);hi=np.array([float(r['high']) for r in rows])
        # The KM table starts at t=0 with certain survival. With no early event,
        # display the exact all-survivor bound throughout the first plotted
        # interval; otherwise an all-censored cohort would show no visible band.
        if len(t)>1 and t[1]>1000:
            t[0]=1000
            lo[0]=.025**(1/int(rows[0]['at_risk']))
        # A cohort whose last event precedes the plannedhorizon has knownzero survival afterward.
        if t[-1]<1e7 and s[-1]==0:
            t=np.r_[t,1e7];s=np.r_[s,0];lo=np.r_[lo,0];hi=np.r_[hi,hi[-1]]
        ax.step(t,s,where='post',color=color,lw=1.2,label=f'$L={L}$')
        ax.fill_between(t,lo,hi,step='post',color=color,alpha=.12,lw=0)
        for r in rows:
            if int(r['censored']):ax.plot(float(r['time']),float(r['survival']),'+',color=color,ms=5,mew=.9)
    ax.set(xscale='log',xlim=(1000,1.2e7),ylim=(-.03,1.04),xlabel='$t$')
    ax.set_ylabel(r'$\Pr(\tau_{1/2}>t)$' if endpoint=='half' else r'$\Pr(\tau_{\rm abs}>t)$')
    ax.legend(loc='lower left',fontsize=6.5)
axes[0].set_title(r'First recorded $\rho<0.145$');axes[1].set_title('Complete consensus')
label(axes[0],'(a)',x=.88,y=.30);label(axes[1],'(b)',x=.88,y=.30)
fig.text(.5,.99,r'$p=\alpha=1/2,\quad q=0.046359;\quad32$ independent runs per size',ha='center',va='bottom',fontsize=8)
save(fig,'fig_lifetimes_v2')
