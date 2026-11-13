"""Figures from complete fixed GPU ensembles; run gpu_campaign reduction first."""
import csv
import json
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from common import ROOT,PAGE,label,save

RESULT=Path(ROOT)/'campaign/v2/results/gpu'


def clear_figure(name):
    for extension in ('pdf','png'):
        (Path(__file__).parent/(name+'.'+extension)).unlink(missing_ok=True)


def table(name):
    path=RESULT/(name+'.csv')
    if not path.exists():return []
    with path.open() as stream:rows=list(csv.DictReader(stream))
    for row in rows:
        for k,v in row.items():
            if v in ('True','False'):row[k]=v=='True'
            else:
                try:row[k]=float(v)
                except ValueError:pass
    return rows


def grouped_curves(rows,value):
    runs=defaultdict(list)
    for row in rows:runs[row['run']].append(row)
    if not runs:return None
    times=np.unique(np.concatenate([[r['time'] for r in rr] for rr in runs.values()]))
    values=[]
    for rr in runs.values():
        rr.sort(key=lambda r:r['time']);x=np.array([r['time'] for r in rr]);y=np.array([r[value] for r in rr])
        values.append(np.interp(times,x,y))
    values=np.asarray(values);mean=values.mean(0)
    rng=np.random.default_rng(1729)
    bs=values[rng.integers(0,len(values),(1000,len(values)))].mean(1)
    lo,hi=np.quantile(bs,[.025,.975],axis=0)
    return times,mean,lo,hi


def curve(ax,rows,value,color,ls='-',labeltext=None,xscale=1.,yscale=1.,band=True):
    data=grouped_curves(rows,value)
    if data is None:return False
    t,m,lo,hi=data;sel=t>=1000
    ax.plot(t[sel]/xscale,m[sel]/yscale,color=color,ls=ls,lw=.95,label=labeltext)
    if band:ax.fill_between(t[sel]/xscale,lo[sel]/yscale,hi[sel]/yscale,color=color,alpha=.12,lw=0)
    return True


def error_points(ax,rows,x,y,low,high,**kwargs):
    rows=sorted(rows,key=lambda r:r[x]);xx=np.array([r[x] for r in rows]);yy=np.array([r[y] for r in rows])
    lower=np.array([r[low] for r in rows]);upper=np.array([r[high] for r in rows])
    if len(rows):ax.errorbar(xx,yy,yerr=[np.maximum(yy-lower,0),np.maximum(upper-yy,0)],capsize=1.4,lw=.7,**kwargs)


def bulk_figure():
    summaries=[r for r in table('bulk_summary') if r['p']==.5]
    if not summaries:
        print('Flagship bulk: awaiting first complete fixed ensemble');clear_figure('fig_gpu_bulk');return
    series=[r for r in table('bulk_curves') if r['p']==.5]
    if not any(r['q']==.04594 for r in series):
        print('Flagship bulk: awaiting complete q=.04594 group for the time panel');clear_figure('fig_gpu_bulk');return
    fig,axs=plt.subplots(1,3,figsize=(PAGE,2.35),gridspec_kw=dict(wspace=.36))
    marks={2048:'o',4096:'s',8192:'^',16384:'D'}
    colors={2048:'C0',4096:'C1',8192:'C3',16384:'C4'}
    for L in marks:
        curve(axs[0],[r for r in series if r['Lx']==L and r['q']==.04594],'rho',colors[L],labeltext=f'$L={L}$',xscale=1e6)
        rr=[r for r in summaries if r['Lx']==L]
        error_points(axs[1],rr,'q','rho_all_mean','rho_all_ci_low','rho_all_ci_high',fmt=marks[L]+'-',color=colors[L],ms=3,mfc='w')
        error_points(axs[2],rr,'q','rho_drift_mean','rho_drift_ci_low','rho_drift_ci_high',fmt=marks[L]+'-',color=colors[L],ms=3,mfc='w')
    axs[0].set(xlabel=r'$t\ (10^6)$',ylabel=r'interface density $\rho$',xlim=(0,5),title='$q=0.04594$')
    axs[0].legend(fontsize=5.5,loc='best')
    axs[1].set(xlabel='$q$',ylabel=r'mean $\rho$, $t=2.5$–$5\times10^6$',ylim=(-.01,.34))
    axs[2].set(xlabel='$q$',ylabel=r'late-window drift $\Delta\rho$');axs[2].axhline(0,color='.5',lw=.5)
    for ax in axs[1:]:
        ax.set_xlim(.045895,.046005);ax.set_xticks([.04590,.04594,.04600])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'));ax.tick_params(axis='x',rotation=35,labelsize=6)
    for ax,s in zip(axs,'abc'):label(ax,f'({s})',x=.86)
    save(fig,'fig_gpu_bulk')


def orientation_figure():
    points=[r for r in table('front_points') if r['p']==.5 and r['T']==2e6 and r['threshold']==.3 and r['estimator']=='envelope'
            and r['t0']==1e6 and r['t1']==2e6]
    fits=json.loads((RESULT/'front_crossings.json').read_text())
    fits=[r for r in fits if r['p']==.5 and r['T']==2e6 and r['threshold']==.3 and r['estimator']=='envelope' and r['model']=='bracket']
    if not all(any(r['L']==L for r in points) for L in (8192,16384)) or not any(r['t1'] in (5e5,1e6,2e6) for r in fits):
        print('Orientation: awaiting complete valid groups for both sizes and a crossing');clear_figure('fig_gpu_fronts');return
    ab={'1-0':('C0','o',r'$(1,0)$'),'2-1':('C1','s',r'$(2,1)$'),'1-1':('C3','^',r'$(1,1)$')}
    fig,axs=plt.subplots(1,3,figsize=(PAGE,2.35),gridspec_kw=dict(wspace=.55))
    axs[1].sharey(axs[0])
    for ax,L in zip(axs[:2],[8192,16384]):
        for normal,(color,mk,lab) in ab.items():
            rr=sorted((r for r in points if r['L']==L and r['ab']==normal),key=lambda r:r['q'])
            if not rr:continue
            ax.errorbar([r['q'] for r in rr],[1e4*r['v'] for r in rr],[1e4*r['sem'] for r in rr],
                        fmt=mk+'-',color=color,ms=3,lw=.7,capsize=1.3,label=lab)
        ax.axhline(0,color='.5',lw=.5);ax.set(xlabel='$q$',title=f'$L={L}$',xlim=(.045875,.045965),xticks=[.04588,.04592,.04596])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'));ax.tick_params(axis='x',rotation=35,labelsize=6)
    axs[0].set_ylabel(r'$v\ (10^{-4})$');axs[0].legend(fontsize=5.7,loc='best')
    windows=[(1e5,5e5),(5e5,1e6),(1e6,2e6)];lmarks={4096:'o',8192:'s',16384:'^'};plotted_sizes=set()
    for li,(L,mk) in enumerate(lmarks.items()):
        for ai,(normal,(color,_,_)) in enumerate(ab.items()):
            selected=[]
            for wi,(start,end) in enumerate(windows):
                ff=[f for f in fits if f['L']==L and f['ab']==normal and f['t0']==start and f['t1']==end]
                if ff and ff[0]['ci_low'] is not None and ff[0]['ci_high'] is not None:
                    f=ff[0];x=wi+(li-1)*.22+(ai-1)*.055;plotted_sizes.add(L)
                    axs[2].errorbar([x],[f['crossing']],yerr=[[f['crossing']-f['ci_low']],[f['ci_high']-f['crossing']]],
                                    fmt=mk,color=color,ms=3,capsize=1,lw=.6,mfc='w' if L==4096 else color)
    axs[2].set(xticks=range(3),xticklabels=['.1–.5','.5–1','1–2'],xlabel=r'fit window $t\ (10^6)$',ylabel=r'front crossing $q^*$')
    axs[2].yaxis.set_major_formatter(FormatStrFormatter('%.6f'))
    axs[2].tick_params(axis='y',labelsize=6)
    axs[2].legend(handles=[Line2D([],[],marker=mk,ls='',color='.3',mfc='w' if L==4096 else '.3',label=f'$L={L}$') for L,mk in lmarks.items() if L in plotted_sizes],fontsize=5.2,loc='best')
    for ax,s in zip(axs,'abc'):label(ax,f'({s})',x=.86,y=.13)
    save(fig,'fig_gpu_fronts')


def orientation5m_figure():
    points=[r for r in table('front_points') if r['p']==.5 and r['T']==5e6 and r['L']==8192
            and r['threshold']==.3 and r['estimator']=='envelope' and r['t0']==2e6 and r['t1']==5e6]
    fits=json.loads((RESULT/'front_crossings.json').read_text())
    fits=[r for r in fits if r['p']==.5 and r['T']==5e6 and r['L']==8192
          and r['threshold']==.3 and r['estimator']=='envelope' and r['model']=='bracket']
    drifts=[r for r in table('front5m_paired_drift') if r['threshold']==.3 and r['estimator']=='envelope']
    if not points or not drifts or not fits:
        print('Independent5m fronts: awaiting complete groups and a crossing');clear_figure('fig_gpu_fronts5m');return
    normals={'1-0':('C0','o',r'$(1,0)$'),'2-1':('C1','s',r'$(2,1)$'),'1-1':('C3','^',r'$(1,1)$')}
    fig,axs=plt.subplots(1,3,figsize=(PAGE,2.4),gridspec_kw=dict(wspace=.55))
    windows=[(1e6,2e6),(2e6,3.5e6),(3.5e6,5e6)]
    for index,(normal,(color,marker,lab)) in enumerate(normals.items()):
        rows=sorted((r for r in points if r['ab']==normal),key=lambda r:r['q'])
        if rows:axs[0].errorbar([r['q'] for r in rows],[r['v']*1e4 for r in rows],
            [r['sem']*1e4 for r in rows],fmt=marker+'-',ms=3,color=color,lw=.7,capsize=1.4,label=lab)
        for wi,window in enumerate(windows):
            chosen=[r for r in fits if r['ab']==normal and (r['t0'],r['t1'])==window
                    and r['ci_low'] is not None and r['ci_high'] is not None]
            if not chosen:continue
            r=chosen[0];x=wi+(index-1)*.15
            axs[1].errorbar([x],[r['crossing']],yerr=[[max(0,r['crossing']-r['ci_low'])],[max(0,r['ci_high']-r['crossing'])]],
                fmt=marker,color=color,ms=3,capsize=1.4,lw=.7)
        rows=[dict(r,drift_mean=r['drift_mean']*1e4,drift_ci_low=r['drift_ci_low']*1e4,drift_ci_high=r['drift_ci_high']*1e4)
              for r in drifts if r['ab']==normal]
        error_points(axs[2],rows,'q','drift_mean','drift_ci_low','drift_ci_high',fmt=marker+'-',color=color,ms=3,mfc='w')
    axs[0].set(ylabel=r'$v\ (10^{-4})$',title=r'$2$–$5\times10^6$');axs[0].legend(fontsize=5.7)
    axs[1].set(xticks=range(3),xticklabels=['1–2','2–3.5','3.5–5'],xlabel=r'fit window $t\ (10^6)$',ylabel=r'front crossing $q^*$')
    axs[1].yaxis.set_major_formatter(FormatStrFormatter('%.7f'))
    axs[1].tick_params(axis='y',labelsize=6)
    axs[2].set(ylabel=r'paired $\Delta v\ (10^{-4})$',title='late minus early')
    for ax in (axs[0],axs[2]):
        ax.axhline(0,color='.5',lw=.5);ax.set(xlabel='$q$',xlim=(.045928,.045952),xticks=[.04593,.04594,.04595])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'));ax.tick_params(axis='x',rotation=35,labelsize=6)
    for ax,s in zip(axs,'abc'):label(ax,f'({s})',x=.86,y=.13)
    save(fig,'fig_gpu_fronts5m')


def near_figure():
    windows=[r for r in table('bulk_window_summary') if r['p']==.995]
    series=[r for r in table('bulk_curves') if r['p']==.995]
    records=[r for r in table('front_records') if r['p']==.995 and r['ensemble_complete']]
    if (not any(r['Lx']==4096 for r in series) or not records
            or not any(r['start']==5e6 and r['end']==1e7 for r in windows)):
        print('Near voter: awaiting complete groups for both bulk and strip panels');clear_figure('fig_gpu_near');return
    qs=[.000632,.000636,.000640,.000644,.000648]
    colors={q:plt.get_cmap('viridis')(i/4) for i,q in enumerate(qs)}
    fig,axs=plt.subplots(2,2,figsize=(PAGE,4.6),gridspec_kw=dict(wspace=.32,hspace=.52))
    a,b,c,d=axs.flat
    for q in qs:
        curve(a,[r for r in series if r['Lx']==4096 and r['q']==q],'rho',colors[q],labeltext=f'${q:.6f}$',xscale=1e6)
        rr=[r for r in records if r['q']==q]
        curve(c,rr,'envelope',colors[q],labeltext=f'${q:.6f}$',xscale=1e6,yscale=4096.)
        curve(d,rr,'interior_ratio',colors[q],xscale=1e6)
    a.set(xlabel=r'$t\ (10^6)$',ylabel=r'bulk interface density $\rho$',title='$p=0.995$, $L=4096$',xlim=(0,10))
    a.legend(fontsize=5.5,title='$q$',title_fontsize=6,ncol=2)
    for L,mk,color in [(2048,'o','C0'),(4096,'s','C1'),(8192,'^','C3')]:
        rr=[r for r in windows if r['Lx']==L and r['start']==5e6 and r['end']==1e7]
        error_points(b,rr,'q','rho_mean','rho_ci_low','rho_ci_high',fmt=mk+'-',color=color,ms=3,label=rf'$L={L}$, $5$–$10\times10^6$')
    rr=[r for r in windows if r['Lx']==8192 and r['start']==1e7 and r['end']==2e7]
    error_points(b,rr,'q','rho_mean','rho_ci_low','rho_ci_high',fmt='^--',color='C3',ms=3,mfc='w',label=r'$L=8192$, $10$–$20\times10^6$')
    b.set(xlabel='$q$',ylabel=r'unconditional window mean $\rho$',xlim=(.000631,.000649),xticks=[.000632,.000640,.000648],ylim=(-.01,.29));b.xaxis.set_major_formatter(FormatStrFormatter('%.6f'));b.tick_params(axis='x',rotation=30,labelsize=6);b.legend(fontsize=5.3)
    c.set(xlabel=r'$t\ (10^6)$',ylabel=r'band envelope $W/(L_y/2)$',title=r'strips: $L_x=2048$, $L_y=8192$',xlim=(0,10),ylim=(0,2.05))
    d.set(xlabel=r'$t\ (10^6)$',ylabel=r'interior activity $\iota_{\rm in}/\iota_{\rm ref}$',title='interior-density diagnostic',xlim=(0,10),ylim=(0,1.15))
    d.axhline(.3,color='.5',lw=.5,ls=':')
    for ax,s in zip(axs.flat,'abcd'):label(ax,f'({s})',x=.88)
    save(fig,'fig_gpu_near')

def coverage_report():
    status=json.loads((RESULT/'status.json').read_text())
    figures={
        'fig_gpu_bulk':[e for e in status['ensembles'] if e['p']==.5 and e['mode']=='bulk'],
        'fig_gpu_fronts':[e for e in status['ensembles'] if e['p']==.5 and e['mode']=='band' and e['T']==2e6],
        'fig_gpu_fronts5m':[e for e in status['ensembles'] if e['p']==.5 and e['mode']=='band' and e['T']==5e6],
        'fig_gpu_near':[e for e in status['ensembles'] if e['p']==.995],
    }
    coverage={}
    for name,ensembles in figures.items():
        complete=[e for e in ensembles if e['ensemble_complete']]
        missing=[{k:v for k,v in e.items() if k not in ('task_ids','incomplete_ids')} for e in ensembles if not e['ensemble_complete']]
        coverage[name]=dict(n_declared_groups=len(ensembles),n_complete_groups=len(complete),
            full_declared_grid_complete=len(complete)==len(ensembles),missing_groups=missing,
            note='A group requires every declared seed. Front plotting additionally requires full-window geometry in every seed.')
    (RESULT/'figure_coverage.json').write_text(json.dumps(coverage,indent=2)+'\n')


if __name__=='__main__':
    coverage_report();bulk_figure();orientation_figure();orientation5m_figure();near_figure()
