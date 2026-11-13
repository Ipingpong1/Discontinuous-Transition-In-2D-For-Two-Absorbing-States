"""Unfiltered near-voter magnetization moments; diagnostic active conditioning.

The primary estimate gives every completed independent run equal weight after
trapezoidal time integration over [2e5,T]. A common [2e5,5e5] window tests the
unequal legacy run durations. Conditional estimates reproduce the legacy
activity cutoff as an explicit diagnostic, with retained fractions reported.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path
import numpy as np
from .io import ROOT, read_timeseries, provenance, write_csv, write_json
from .bulk import extend_absorbed
from .stats import window_mean, mean_summary


def conditional_moments(ts, start, end, cutoff):
    """Integrate piecewise-linear moment observations where interpolated rho>cut."""
    t = ts[:,0]
    inside = (t>start)&(t<end)
    tt = np.r_[start,t[inside],end]
    rho = np.interp(tt,t,ts[:,2])
    m = 2*ts[:,1]-1
    yy = np.column_stack([np.interp(tt,t,m**2),np.interp(tt,t,m**4),np.interp(tt,t,abs(m))])
    duration,total = 0.,np.zeros(3)
    for i,dt in enumerate(np.diff(tt)):
        a,b = rho[i]>cutoff,rho[i+1]>cutoff
        if not a and not b: continue
        left,right = 0.,1.
        if a != b:
            cross=(cutoff-rho[i])/(rho[i+1]-rho[i])
            if a: right=cross
            else: left=cross
        ya=yy[i]+left*(yy[i+1]-yy[i]);yb=yy[i]+right*(yy[i+1]-yy[i])
        weight=(right-left)*dt
        total+=weight*(ya+yb)/2;duration+=weight
    return (total/duration if duration else np.full(3,np.nan)),duration/(end-start)


def moments_summary(rows, mode, boot, rng):
    prefix='all_' if mode=='all' else 'conditional_'
    usable=[r for r in rows if np.isfinite(r[prefix+'m2'])]
    n=len(usable)
    if not n:return None
    a=np.array([[r[prefix+'m2'],r[prefix+'m4'],r[prefix+'abs_m']] for r in usable])
    size=rows[0]['Lx']*rows[0]['Ly']
    m2,m4,am=a.mean(0);binder=1-m4/(3*m2*m2) if m2>0 else np.nan
    cs=mean_summary(size*a[:,0],boot=boot,rng=rng)
    result=dict(estimator=mode,n_total=len(rows),n_used=n,n_absorbed=sum(r['absorbed'] for r in rows),
        chi=cs['mean'],chi_sem=cs['sem'],chi_ci_low=cs['ci_low'],chi_ci_high=cs['ci_high'],
        binder=binder,binder_ci_low=np.nan,binder_ci_high=np.nan,mean_abs_m=am,
        retained_time_fraction=float(np.mean([r['active_fraction'] for r in rows])),
        retained_run_fraction=n/len(rows) if mode!='all' else 1.,
        min_T=min(r['T'] for r in rows),max_T=max(r['T'] for r in rows))
    if n>1 and boot:
        means=a[rng.integers(0,n,(boot,n))].mean(1)
        bs=1-means[:,1]/(3*means[:,0]**2)
        result['binder_ci_low'],result['binder_ci_high']=np.quantile(bs,[.025,.975])
    return result


def analyze(paths,boot=2000):
    groups=defaultdict(list);skipped=[];seen=set()
    for path in paths:
        h,ts=read_timeseries(path)
        if not h['done']:
            skipped.append(str(path));continue
        key=(h['p'],h['q'],h['Lx'],h['Ly'],h['seed'])
        if key in seen:raise ValueError(f'duplicate dynamical seed identity: {key}')
        seen.add(key)
        if not h['absorbed'] and ts[-1,0]!=h['T']:
            raise ValueError(f'completed output inconsistent with planned horizon: {path}')
        if h['absorbed']:ts=extend_absorbed(ts,h['T'])
        if h['T']<=200000:continue
        groups[(h['p'],h['q'],h['Lx'],h['Ly'])].append((path,h,ts))
    allrows=[];summaries=[];ratios=[];bygroup={};rng=np.random.default_rng(20261103)
    for key,items in sorted(groups.items()):
        p,q,Lx,Ly=key
        # Deliberately preserve the old, data-dependent pooled-record threshold
        # only in the explicitly labelled diagnostic conditional estimator.
        raw=np.concatenate([ts[ts[:,0]>=200000,2] for _,_,ts in items])
        cutoff=.8*np.median(raw)
        for window in ['planned','common_500k']:
            rows=[]
            for path,h,ts in items:
                end=h['T'] if window=='planned' else min(h['T'],500000)
                t=ts[:,0];m=2*ts[:,1]-1
                cm,fraction=conditional_moments(ts,200000,end,cutoff)
                row=dict(path=str(path.resolve()),p=p,q=q,Lx=Lx,Ly=Ly,seed=h['seed'],
                    T=h['T'],window=window,start=200000,end=end,absorbed=h['absorbed'],
                    all_m2=window_mean(t,m*m,200000,end),all_m4=window_mean(t,m**4,200000,end),
                    all_abs_m=window_mean(t,abs(m),200000,end),
                    conditional_m2=cm[0],conditional_m4=cm[1],conditional_abs_m=cm[2],
                    active_cutoff=cutoff,active_fraction=fraction)
                rows.append(row);allrows.append(row)
            for mode in ['all','active_conditional']:
                sm=moments_summary(rows,mode,boot,rng)
                if sm is None:continue
                result=dict(p=p,q=q,Lx=Lx,Ly=Ly,window=window,**sm)
                summaries.append(result);bygroup[(p,q,Lx,window,mode)]=(rows,result)
    for key,(rows,sm) in sorted(bygroup.items()):
        p,q,L,window,mode=key
        if L!=4096 or (p,q,2048,window,mode) not in bygroup:continue
        lowrows,lowsm=bygroup[(p,q,2048,window,mode)]
        pref='all_' if mode=='all' else 'conditional_'
        hi=np.array([r[pref+'m2']*4096**2 for r in rows]);hi=hi[np.isfinite(hi)]
        lo=np.array([r[pref+'m2']*2048**2 for r in lowrows]);lo=lo[np.isfinite(lo)]
        ci=[np.nan,np.nan]
        if len(hi)>1 and len(lo)>1 and boot:
            samples=hi[rng.integers(0,len(hi),(boot,len(hi)))].mean(1)/lo[rng.integers(0,len(lo),(boot,len(lo)))].mean(1)
            ci=np.quantile(samples,[.025,.975])
        ratios.append(dict(p=p,q=q,window=window,estimator=mode,
            n2048=len(lo),n4096=len(hi),ratio=sm['chi']/lowsm['chi'],ci_low=ci[0],ci_high=ci[1],
            retained2048=lowsm['retained_time_fraction'],retained4096=sm['retained_time_fraction']))
    return allrows,summaries,ratios,skipped


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/'campaign/v2/results/susceptibility')
    ap.add_argument('--boot',type=int,default=2000);args=ap.parse_args()
    paths=sorted(list((ROOT/'data/runs').glob('E6_long*/out/*.txt'))+list((ROOT/'data/runs').glob('E7_binder*/out/*.txt')))
    runs,summary,ratios,skipped=analyze(paths,args.boot)
    write_csv(args.out/'runs.csv',runs);write_csv(args.out/'summary.csv',summary);write_csv(args.out/'ratios.csv',ratios)
    write_json(args.out/'analysis.json',dict(settings=dict(start=200000,windows=['planned_T','common_500000'],
        primary='all completed runs; independent run weighting; trapezoidal time moments',
        conditional='diagnostic only; rho>0.8*group pooled-record median; no minimum-record selection',
        bootstrap=args.boot,ci='95% independent-run percentile bootstrap',
        binder='1-mean_run(m4)/(3*mean_run(m2)^2)'),
        inputs=provenance(paths),skipped_incomplete=skipped,summary=summary,ratios=ratios))
    print(len(paths),'inputs;',len(summary),'summary groups;',len(ratios),'size ratios')
    for r in ratios:
        if r['window']=='planned' and r['estimator']=='all':print(r)
if __name__=='__main__':main()
