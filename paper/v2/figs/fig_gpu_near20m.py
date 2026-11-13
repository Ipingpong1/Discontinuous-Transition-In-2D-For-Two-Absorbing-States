"""Largest-system near-voter trajectories and paired late density drift.

Requires all four predeclared four-run20M cohorts. Uses the existing unconditional
curves and time-weighted paired drift; no fit, extrapolation, or survivor filter.
"""
from collections import defaultdict
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from common import ROOT,PAGE,label,save
from fig_gpu_v2 import table,clear_figure
from campaign.v2.analysis.io import provenance,write_json

RESULT=Path(ROOT)/'campaign/v2/results/gpu'
QS=(.000640,.000642,.000644,.000648)
NAME='fig_gpu_near20m'


def main():
    status=json.loads((RESULT/'status.json').read_text())
    cohorts=[e for e in status['ensembles'] if e['p']==.995 and e['mode']=='bulk'
        and e['Lx']==8192 and e['Ly']==8192 and e['T']==20000000 and e['q'] in QS]
    ready=(len(cohorts)==4 and {e['q'] for e in cohorts}==set(QS)
        and all(e['ensemble_complete'] and e['n_complete']==4 for e in cohorts))
    write_json(RESULT/(NAME+'_coverage.json'),dict(n_required_cohorts=4,
        n_complete=sum(e['ensemble_complete'] for e in cohorts),all_complete=ready,
        q=list(QS),L=8192,T=20000000))
    if not ready:
        clear_figure(NAME);raise SystemExit('Largest near-voter figure awaits all four complete fixed cohorts')
    series=[r for r in table('bulk_curves') if r['p']==.995 and r['Lx']==8192 and r['T']==2e7]
    windows=[r for r in table('bulk_window_summary') if r['p']==.995 and r['Lx']==8192
        and r['T']==2e7 and r['start']==1e7 and r['end']==2e7]
    resamples=np.indices((4,4,4,4)).reshape(4,-1).T
    colors=['#542788','#2b8cbe','#238b45','#d95f0e']
    fig,axs=plt.subplots(1,2,figsize=(PAGE,2.65),gridspec_kw=dict(width_ratios=[1.2,1],wspace=.32))
    for q,color in zip(QS,colors):
        runs=defaultdict(list)
        for row in series:
            if row['q']==q:runs[row['run']].append(row)
        if len(runs)!=4:raise ValueError('Curve seed count differs from fixed cohort')
        records=[sorted(runs[key],key=lambda row:row['time']) for key in sorted(runs)]
        times=np.array([row['time'] for row in records[0]])
        if not all(np.array_equal(times,[r['time'] for r in rows]) for rows in records):raise ValueError('Nonmatching curve time grids')
        values=np.array([[r['rho'] for r in rows] for rows in records])
        mean=values.mean(0);lo,hi=np.quantile(values[resamples].mean(1),[.025,.975],axis=0)
        selected=times>=1000
        axs[0].plot(times[selected]/1e6,mean[selected],color=color,label=f'${q:.6f}$')
        axs[0].fill_between(times[selected]/1e6,lo[selected],hi[selected],color=color,alpha=.16,lw=0)
        row=next(r for r in windows if r['q']==q)
        y=row['rho_drift_mean'];low=row['rho_drift_ci_low'];high=row['rho_drift_ci_high']
        axs[1].errorbar(q,y,yerr=[[y-low],[high-y]],fmt='o',color=color,ms=4,capsize=2,lw=.9)
    axs[0].set(xlim=(0,20),ylim=(0,.27),xticks=[0,5,10,15,20],
        xlabel=r'$t\ (10^6)$',ylabel=r'unconditional interface density $\rho$',title=r'$p=0.995,\ L=8192$')
    axs[0].legend(title='$q$',title_fontsize=7,fontsize=6.3,ncol=2,loc='lower left')
    axs[1].axhline(0,color='.5',lw=.7)
    axs[1].set(xlim=(.000639,.000649),xticks=QS,ylim=(-.105,.015),
        xlabel='$q$',ylabel=r'paired late density drift $\Delta\rho$')
    axs[1].xaxis.set_major_formatter(FormatStrFormatter('%.6f'))
    axs[1].tick_params(axis='x',rotation=30,labelsize=6.5)
    for ax,letter in zip(axs,'ab'):label(ax,f'({letter})',x=.88)
    save(fig,NAME)
    write_json(RESULT/(NAME+'_provenance.json'),dict(
        inputs=provenance([RESULT/'status.json',RESULT/'bulk_curves.csv',RESULT/'bulk_window_summary.csv',Path(__file__)]),
        estimator='Unconditional run mean; known absorbing extensions retained; paired time-weighted15–20M minus10–15M means.',
        uncertainty='Pointwise95% empirical run-bootstrap percentile intervals. All256 ordered four-run resamples enumerated; approximate confidence coverage.'))
    print(NAME)

if __name__=='__main__':main()
