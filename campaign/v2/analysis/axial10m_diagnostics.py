"""Standalone paired-window diagnostic for the independent axial10M declaration.

The original2M and5M figures/contrasts are unchanged. Only complete four-run
cohorts with full geometry in both fixed windows contribute a paired contrast.
"""
import csv,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from .fronts import covers_requested_window
from .front_stationarity import bracket_zero
from .io import ROOT,provenance,write_csv,write_json
from .stats import group_rng,mean_summary

EARLY=(2e6,5e6)
LATE=(5e6,1e7)
GRID=(.04593,.04594,.04595)


def paired_diagnostics(rows,boot=1000):
    groups=defaultdict(lambda:defaultdict(lambda:defaultdict(dict)))
    for row in rows:
        if (row['T'],row['L'],row['ab'],row['estimator'])!=(1e7,8192,'1-0','envelope'):continue
        if not row['ensemble_complete'] or not covers_requested_window(row):continue
        window=(row['t0'],row['t1'])
        if window not in (EARLY,LATE):continue
        groups[row['threshold']][row['q']][window][row['run']]=row['v']
    drifts=[];crossings=[]
    for threshold in (.2,.3,.5):
        paired={}
        for q,windows in sorted(groups[threshold].items()):
            a,b=windows.get(EARLY,{}),windows.get(LATE,{})
            if len(a)!=4 or set(a)!=set(b):continue
            values=np.array([[a[run],b[run]] for run in sorted(a)])
            paired[q]=values
            summary=mean_summary(values[:,1]-values[:,0],boot,group_rng(('axial10m',threshold,q),39171))
            drifts.append(dict(T=10000000,L=8192,ab='1-0',threshold=threshold,q=q,
                early_start=EARLY[0],early_end=EARLY[1],late_start=LATE[0],late_end=LATE[1],
                n=4,v_early=float(values[:,0].mean()),v_late=float(values[:,1].mean()),
                **{'drift_'+k:v for k,v in summary.items()}))
        result=dict(threshold=threshold,n_valid_paired_q=len(paired),required_q=list(GRID),
            available_q=sorted(paired),status='incomplete_or_geometry_invalid_paired_grid')
        if tuple(sorted(paired))==GRID:
            means=np.array([paired[q].mean(0) for q in GRID])
            early=bracket_zero(GRID,means[:,0]);late=bracket_zero(GRID,means[:,1])
            result.update(early_velocities=means[:,0].tolist(),late_velocities=means[:,1].tolist(),
                q_early=early,q_late=late,status='unbracketed_or_nonunique_reversal')
            if np.isfinite(early+late):
                rng=group_rng(('axial10m',threshold,'paired_roots'),39171);samples=[]
                for _ in range(boot):
                    vv=np.array([paired[q][rng.integers(0,4,4)].mean(0) for q in GRID])
                    a,b=bracket_zero(GRID,vv[:,0]),bracket_zero(GRID,vv[:,1])
                    if np.isfinite(a+b):samples.append(b-a)
                sufficient=len(samples)>=.8*max(1,boot)
                low,high=np.quantile(samples,[.025,.975]) if sufficient else (np.nan,np.nan)
                result.update(difference=late-early,ci95_low=float(low),ci95_high=float(high),
                    bootstrap_valid=len(samples),bootstrap_requested=boot,
                    status='paired_crossing_interval' if sufficient else 'paired_crossing_interval_withheld')
        crossings.append(result)
    return drifts,crossings


def main():
    folder=ROOT/'campaign/v2/results/gpu'
    source=folder/'front_runs.csv'
    rows=[]
    with source.open() as stream:
        for row in csv.DictReader(stream):
            if float(row['T'])!=1e7 or row['p']!='0.5':continue
            for key in ('T','L','q','threshold','t0','t1','v','t_start','t_end'):row[key]=float(row[key])
            row['ensemble_complete']=row['ensemble_complete']=='True';rows.append(row)
    drifts,crossings=paired_diagnostics(rows)
    status=json.loads((folder/'status.json').read_text())
    groups=[e for e in status['ensembles'] if e['p']==.5 and e['mode']=='band' and e['T']==1e7]
    if drifts:write_csv(folder/'axial10m_paired_drift.csv',drifts)
    else:(folder/'axial10m_paired_drift.csv').unlink(missing_ok=True)
    write_json(folder/'axial10m_paired_crossings.json',dict(
        n_declared_groups=len(groups),n_complete_groups=sum(e['ensemble_complete'] for e in groups),
        windows=[EARLY,LATE],crossings=crossings,
        per_run_diagnostics=[dict(row,full_window_geometry=covers_requested_window(row))
            for row in rows if row['ensemble_complete'] and row['estimator']=='envelope'
            and (row['t0'],row['t1']) in (EARLY,LATE)],
        interval_scope='Paired run bootstrap conditional on tracker, windows, local interpolation, and both resampled curves having unique in-grid reversals (boundary zeros are one-sided). All three declared four-run groups required; zero-containing intervals do not establish equivalence.',
        mean_interval_method='All256 ordered four-run empirical bootstrap resamples; confidence coverage remains approximate.',
        provenance=provenance([source,folder/'status.json',Path(__file__),
            ROOT/'campaign/v2/analysis/front_stationarity.py',ROOT/'campaign/v2/analysis/stats.py',
            ROOT/'campaign/v2/hpc/gpu/production_manifest_axial10m.json'])))
    print(f'Axial10M: {sum(e["ensemble_complete"] for e in groups)}/3 complete groups; {len(drifts)} paired parameter/threshold means')

if __name__=='__main__':main()
