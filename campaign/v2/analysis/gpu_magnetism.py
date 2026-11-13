"""Common-window magnetic size ratios from complete GPU manifest ensembles."""
from collections import defaultdict
from itertools import product
from pathlib import Path
import argparse,csv
import numpy as np
from .io import ROOT,provenance,write_json
from .stats import group_rng


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,default=ROOT/'campaign/v2/results/gpu/bulk_windows.csv')
    ap.add_argument('--runs',type=Path,default=ROOT/'campaign/v2/results/gpu/bulk_runs.csv')
    ap.add_argument('--out',type=Path,default=ROOT/'campaign/v2/results/gpu_magnetism')
    ap.add_argument('--boot',type=int,default=5000);a=ap.parse_args()
    groups=defaultdict(list)
    with a.runs.open() as f:
        runs={r['run']:r for r in csv.DictReader(f)}
    with a.source.open() as f:
        for r in csv.DictReader(f):
            if r['ensemble_complete']!='True' or int(r['Lx'])!=int(r['Ly']):continue
            key=(float(r['p']),float(r['q']),int(r['Lx']),int(float(r['start'])),int(float(r['end'])))
            groups[key].append(r)
    for records in groups.values():
        records.sort(key=lambda row:row['run'])
    rows=[]
    for key,large in sorted(groups.items()):
        p,q,L,start,end=key;small=groups.get((p,q,L//2,start,end))
        if not small:continue
        rng=group_rng(key,20261103)
        hi=np.array([float(r['chi']) for r in large]);lo=np.array([float(r['chi']) for r in small])
        if len(hi)<2 or len(lo)<2:continue
        if len(hi)<=4 and len(lo)<=4:
            # At four runs there are only 256 ordered bootstrap resamples per
            # size. Enumerate their product to remove Monte Carlo variability;
            # this does not make bootstrap confidence coverage exact.
            hm=hi[np.array(list(product(range(len(hi)),repeat=len(hi))))].mean(1)
            lm=lo[np.array(list(product(range(len(lo)),repeat=len(lo))))].mean(1)
            samples=(hm[:,None]/lm[None,:]).ravel()
            bootstrap_method='Exhaustive independent-run percentile bootstrap; approximate confidence coverage'
        else:
            samples=hi[rng.integers(0,len(hi),(a.boot,len(hi)))].mean(1)/lo[rng.integers(0,len(lo),(a.boot,len(lo)))].mean(1)
            bootstrap_method=f'{a.boot} independent-run bootstrap resamples; approximate confidence coverage'
        low,high=np.quantile(samples,[.025,.975])
        def absorbed(rows,time):
            return sum(runs[r['run']]['absorption_event']=='True'
                       and float(runs[r['run']]['absorption_time'])<=time for r in rows)
        ns0,ns1=absorbed(small,start),absorbed(small,end)
        nl0,nl1=absorbed(large,start),absorbed(large,end)
        rows.append(dict(p=p,q=q,L_small=L//2,L_large=L,start=start,end=end,n_small=len(small),n_large=len(large),
            chi_small=lo.mean(),chi_large=hi.mean(),ratio=hi.mean()/lo.mean(),ci_low=low,ci_high=high,
            bootstrap_method=bootstrap_method,bootstrap_ratio_samples=len(samples),
            rho_small=np.mean([float(r['rho']) for r in small]),rho_large=np.mean([float(r['rho']) for r in large]),
            absorbed_by_start_small=ns0,absorbed_by_end_small=ns1,
            absorbed_by_start_large=nl0,absorbed_by_end_large=nl1,
            absorbed_by_horizon_small=sum(runs[r['run']]['absorption_event']=='True' for r in small),
            absorbed_by_horizon_large=sum(runs[r['run']]['absorption_event']=='True' for r in large),
            absorption_free_window=not(ns1 or nl1),
            max_abs_density_drift_small=max(abs(float(r['rho_drift'])) for r in small),
            max_abs_density_drift_large=max(abs(float(r['rho_drift'])) for r in large)))
    a.out.mkdir(parents=True,exist_ok=True)
    if rows:
        with (a.out/'ratios.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    else:
        (a.out/'ratios.csv').unlink(missing_ok=True)
    write_json(a.out/'analysis.json',dict(bootstrap=a.boot,bootstrap_seed=20261103,
        bootstrap_streams='Independent deterministic stream per parameter/size/window key; runs sorted by identity',
        ratios=rows,inputs=provenance([a.source,a.runs]),
        method='Unconditioned independent-run bootstrap ratio; complete fixed ensembles and identical time windows. Exhaustive resampling when both sample sizes are at most four; Monte Carlo otherwise. Bootstrap confidence coverage remains approximate.',
        caveat='A magnetic secondmoment includes consensus m²=1. Interpret a size ratio together with densitystationarity and absorption diagnostics.'))
    print(len(rows),'common-window size ratios')
    for r in rows:
        if r['p']==.995 and r['start']>=5_000_000:print(r)
if __name__=='__main__':main()
