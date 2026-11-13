"""Paired window contrasts in the independent, fixed four-seed5m cohorts."""
from collections import defaultdict
import numpy as np
from .stats import mean_summary,group_rng

EARLY=(2_000_000.,3_500_000.)
LATE=(3_500_000.,5_000_000.)
GRID=(.04593,.04594,.04595)


def bracket_zero(q,velocity):
    q=np.asarray(q,float);velocity=np.asarray(velocity,float)
    zeros=np.flatnonzero(velocity==0)
    if len(zeros)>1:return np.nan
    roots=[]
    for k in zeros:
        # An isolated interior zero with same-sign neighbours is a tangent,
        # not a velocity reversal. A boundary zero is explicitly one-sided.
        if k in (0,len(q)-1) or velocity[k-1]*velocity[k+1]<0:roots.append(q[k])
    for k in np.flatnonzero(velocity[:-1]*velocity[1:]<0):
        roots.append(q[k]-velocity[k]*(q[k+1]-q[k])/(velocity[k+1]-velocity[k]))
    return float(roots[0]) if len(roots)==1 else np.nan


def paired_contrasts(valid_front,boot=1000):
    """Input already requires complete declared seed groups and full geometry.

    A contrast additionally requires identical runs valid in both windows.
    Crossing differences use all three declared q groups for each configuration,
    resampling the same run indices across its early/late velocities at each q.
    """
    names=('p','L','Ly','T','ab','threshold','smooth_normal','estimator')
    groups=defaultdict(lambda:defaultdict(lambda:defaultdict(dict)))
    for row in valid_front:
        if row['p']!=.5 or row['T']!=5_000_000 or row['L']!=8192:continue
        window=(row['t0'],row['t1'])
        if window not in (EARLY,LATE):continue
        groups[tuple(row[k] for k in names)][row['q']][window][row['run']]=row['v']
    rng=np.random.default_rng(53829);drifts=[];crossings=[]
    for key,by_q in sorted(groups.items()):
        base=dict(zip(names,key));paired={}
        for q,windows in sorted(by_q.items()):
            first,last=windows.get(EARLY,{}),windows.get(LATE,{})
            if len(first)!=4 or set(first)!=set(last):continue
            runs=sorted(first);values=np.array([[first[i],last[i]] for i in runs])
            paired[q]=values;delta=values[:,1]-values[:,0]
            summary=mean_summary(delta,boot,group_rng(key+(q,),53829))
            drifts.append(dict(**base,q=q,early_start=EARLY[0],early_end=EARLY[1],late_start=LATE[0],late_end=LATE[1],
                v_early_mean=float(values[:,0].mean()),v_late_mean=float(values[:,1].mean()),
                **{'drift_'+k:v for k,v in summary.items()}))
        if tuple(sorted(paired))!=GRID:continue
        means=np.array([paired[q].mean(0) for q in GRID])
        early=bracket_zero(GRID,means[:,0]);late=bracket_zero(GRID,means[:,1])
        if not np.isfinite(early+late):continue
        samples=[];rng=group_rng(key,53829)
        for _ in range(boot):
            values=np.array([paired[q][rng.integers(0,4,4)].mean(0) for q in GRID])
            a,b=bracket_zero(GRID,values[:,0]),bracket_zero(GRID,values[:,1])
            if np.isfinite(a+b):samples.append(b-a)
        low,high=np.quantile(samples,[.025,.975]) if len(samples)>=.8*max(1,boot) else (np.nan,np.nan)
        crossings.append(dict(**base,early_start=EARLY[0],early_end=EARLY[1],late_start=LATE[0],late_end=LATE[1],
            q_early=early,q_late=late,difference=late-early,ci95_low=float(low),ci95_high=float(high),
            bootstrap_requested=boot,bootstrap_valid=len(samples),n_q=3,n_per_q=4,
            interval_scope='paired run bootstrap conditional on tracker, windows, local interpolation, and both resampled curves having a unique in-grid velocity reversal (boundary zeros are one-sided); all declared q groups required; zero-containing intervals do not establish equivalence'))
    return drifts,crossings
