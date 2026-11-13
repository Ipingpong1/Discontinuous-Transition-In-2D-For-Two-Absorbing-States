"""Two-droplet test at actual q; curvature overlays use corrected front data."""
import csv
import glob
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from common import COL, ROOT, RUNS, label, save
sys.path.insert(0,os.path.join(ROOT,'campaign','round2'))
import drop2

D,R0=.20,1000
NAME=re.compile(r'p0\.5_a0\.5_j2_q([\d.]+)_(disk|hole)R(\d+)_L(\d+)_s(\d+)$')


def corrected_velocities():
    path=Path(ROOT)/'campaign/v2/results/fronts/points.csv'
    rows=list(csv.DictReader(path.open()))
    selected=[r for r in rows if float(r['p'])==.5 and int(float(r['L']))==4096
              and r['ab']=='1-0' and float(r['threshold'])==.3
              and float(r['t0'])==20000 and r['estimator']=='envelope'
              and r.get('t1','inf') in ('inf','nan','')]
    selected.sort(key=lambda r:float(r['q']))
    return np.array([float(r['q']) for r in selected]),np.array([float(r['v']) for r in selected])


def main():
    groups=defaultdict(list)
    for f in sorted(glob.glob(os.path.join(RUNS,'R2_drop2','out','*.txt'))):
        prefix=f[:-4];m=NAME.search(os.path.basename(prefix))
        if not m or not os.path.exists(prefix+'.prof'):continue
        q,kind=float(m[1]),m[2]
        L,ts,prof,done=drop2.load(prefix)
        t,R_ts,R_pr=drop2.radii(L,ts,prof,kind,R0)
        groups[q,kind].append((R_pr,Path(f).stem))
    qs=sorted({q for q,kind in groups});colors=dict(zip(qs,['C0','C9','C1','C3']))
    fq,fv=corrected_velocities()
    fig,axs=plt.subplots(1,2,figsize=(2*COL,2.2),gridspec_kw=dict(wspace=.28),sharey=True)
    measurements=[]
    for ax,kind in zip(axs,['disk','hole']):
        for q in qs:
            rs=groups[q,kind]
            n=min(len(r[0][0]) for r in rs);t=rs[0][0][0][:n]
            mean=np.mean([r[0][1][:n] for r in rs],axis=0)
            for (time,radius),name in rs:
                ax.plot(time/1e3,radius,color=colors[q],lw=.4,alpha=.4)
            ax.plot(t/1e3,mean,color=colors[q],lw=1.2,label=f'${q:.7f}$')
            assert fq.min()<=q<=fq.max(),'Do not extrapolate front velocities'
            v=float(np.interp(q,fq,fv));rate=v-D/R0 if kind=='disk' else -v-D/R0
            sel=t>=1e4
            ax.plot(t[sel]/1e3,mean[sel][0]+rate*(t[sel]-t[sel][0]),'k--',lw=.6)
            slopes=[]
            for (time,radius),name in rs:
                fit=(time>=2e4)&(time<=2e5)&np.isfinite(radius)
                slope=float(np.polyfit(time[fit],radius[fit],1)[0]);slopes.append(slope)
                measurements.append(dict(q=q,kind=kind,run=name,t_start=float(time[fit][0]),
                    t_end=float(time[fit][-1]),radius_slope=slope,front_v=v,curvature_rate=rate,
                    difference=slope-rate,D=D,R0=R0))
            print(f'q={q:.7f} {kind}: measured dR/dt={np.mean(slopes):.8g} +/- {np.std(slopes,ddof=1)/np.sqrt(len(slopes)):.3g}; predicted={rate:.8g}; difference={np.mean(slopes)-rate:.3g}')
        ax.set(xlabel=r'$t$ ($10^3$ steps)',xlim=(0,200));ax.axhline(R0,color='.7',lw=.5,ls=':')
    axs[0].set(ylabel=r'radius $R$',title='active disk in the dead state')
    axs[1].set(title='dead hole in the active state')
    axs[0].legend(fontsize=5.7,title='$q$',title_fontsize=6.5,loc='upper left',ncol=2)
    for ax,s in zip(axs,'ab'):label(ax,f'({s})',x=.9,y=.1)
    save(fig,'fig_drop2')
    out=Path(ROOT)/'campaign/v2/results/droplets';out.mkdir(exist_ok=True)
    with (out/'two_drop_rates.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(measurements[0]));writer.writeheader();writer.writerows(measurements)

if __name__=='__main__':main()
