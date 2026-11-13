#!/usr/bin/env python3
"""Extended wall test: fixed ensembles, all-site adsorption, and legacy extrapolation.

Run the gpu_wetting analysis first. This figure refuses stale q* inputs and does
not generate a publication artifact before a complete wall+baseline pair exists.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from common import ROOT, PAGE, QSTAR05, QSTAR05_HALF_RANGE, label, save
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DATA = Path(ROOT)/'campaign/v2/results/gpu_wetting'
QCOLORS = {0.04584:'#246d9c',0.04588:'#bd5b15',0.04590:'#833796'}
STYLES = {4096:'--',16384:'-'}
MARKERS = {4096:'o',16384:'s'}


def read_rows(path):
    with path.open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    for row in rows:
        for key,value in list(row.items()):
            if value in ('True','False'):
                row[key]=value=='True'
            else:
                try:row[key]=float(value)
                except ValueError:pass
    return rows


def draw(data=DATA):
    status=json.loads((data/'status.json').read_text())
    if status['n_adsorption_points']==0:
        raise SystemExit('No complete nonzero-q wall and matching baseline ensemble; figure not generated.')
    analysis=json.loads((data/'analysis.json').read_text())
    if (abs(analysis['reference']['qstar']-QSTAR05)>1e-12 or
        abs(analysis['reference']['systematic_half_range']-QSTAR05_HALF_RANGE)>1e-12):
        raise SystemExit('Shared coexistence reference changed; rerun gpu_wetting before this figure.')
    points=read_rows(data/'points.csv');curves=read_rows(data/'curves.csv')
    if len(points)!=status['n_adsorption_points']:
        raise SystemExit('Wall reduction is being refreshed; retry after the serial reducer completes.')
    model=analysis['legacy_comparison']
    fig,axes=plt.subplots(1,2,figsize=(PAGE,2.65),gridspec_kw={'width_ratios':[1.02,1.]})
    ax=axes[0]
    groups={}
    displayed_sizes={int(row['Lx']) for row in points}
    for row in curves:
        if int(row['Lx']) not in displayed_sizes:continue
        groups.setdefault((row['q'],int(row['Lx'])),[]).append(row)
    for (q,Lx),rows in sorted(groups.items()):
        rows=sorted(rows,key=lambda r:r['time'])
        rows=[r for r in rows if r['time']>0]
        t=np.array([r['time'] for r in rows]);mean=np.array([r['I'] for r in rows])
        lo=np.array([r['ci_low'] for r in rows]);hi=np.array([r['ci_high'] for r in rows])
        color=QCOLORS.get(q,'0.5')
        ax.plot(t/1e6,mean,color=color,ls=STYLES[Lx],lw=1.)
        ax.fill_between(t/1e6,lo,hi,color=color,alpha=.10,lw=0)
    ax.axvspan(5,10,color='0.2',alpha=.04,zorder=-5)
    ax.set(xlabel=r'$t/10^6$',ylabel=r'$I(t)=L_y\iota(t)/2$',xlim=(0,10))
    qhandles=[Line2D([],[],color=QCOLORS[q],label=rf'$q={q:.5f}$')
              for q in sorted({r['q'] for r in points})]
    if any(q==0 for q,_ in groups):qhandles.append(Line2D([],[],color='0.5',label=r'$q=0$ baseline'))
    lg=ax.legend(handles=qhandles,loc='center left',bbox_to_anchor=(.02,.30),handlelength=1.5)
    ax.add_artist(lg)
    ax.legend(handles=[Line2D([],[],color='0.25',ls=STYLES[L],label=rf'$L_x={L}$')
                       for L in sorted({int(r['Lx']) for r in points})],loc='lower right',handlelength=2.)
    label(ax,'(a)')
    ax=axes[1]
    legacy=model.get('fit_points',[]) if model['status']=='ok' else []
    if legacy:
        x=np.array([1-r['q']/QSTAR05 for r in legacy]);y=np.array([r['Gamma'] for r in legacy])
        se=np.array([r['Gamma_sem'] for r in legacy])
        ax.errorbar(x,y,yerr=se,fmt='.',color='0.55',ms=4,lw=.6,capsize=1.5,zorder=2)
        positive=[r['delta_low'] for r in points if r['delta_low']>1e-12]+[r['delta'] for r in points if r['delta']>0]
        xmin=min(positive+[min(x)])*.65;xmax=max(x)*1.15
        delta=np.geomspace(xmin,xmax,400)
        measured=(delta>=model['fit_delta_min'])&(delta<=model['fit_delta_max'])
        for kind,color,style in [('power','0.15',':'),('third_offset','0.5','--')]:
            fit=model['central'][kind]
            if kind=='power':prediction=fit['A']*delta**(-fit['psi'])
            else:prediction=fit['A']*delta**(-1/3)+fit['B']
            ax.plot(delta,np.where(measured,prediction,np.nan),color=color,lw=.9)
            ax.plot(delta,np.where(~measured,prediction,np.nan),color=color,ls=style,lw=.9)
            parameters=np.asarray(model.get('parameter_bootstrap',{}).get(kind,[]),float)
            if len(parameters):
                if kind=='power':draws=parameters[:,0,None]*delta[None,:]**(-parameters[:,1,None])
                else:draws=parameters[:,0,None]*delta[None,:]**(-1/3)+parameters[:,1,None]
                low,high=np.quantile(draws,[.025,.975],axis=0)
                ax.fill_between(delta,low,high,color=color,alpha=.10,lw=0,zorder=-2)
    else:
        positive=[r['delta_low'] for r in points if r['delta_low']>1e-12]+[r['delta'] for r in points if r['delta']>0]
        if not positive:raise SystemExit('All q values exceed the current coexistence reference.')
        xmin=min(positive)*.6;xmax=max(positive)*1.5
    crosses=False;nonpositive=False
    for row in points:
        x,y=row['delta'],row['Gamma']
        if x<=0 or y<=0:
            nonpositive=True
            continue
        Lx=int(row['Lx']);color=QCOLORS[row['q']]
        low=row['delta_low'];high=row['delta_high']
        # q* sensitivity is not a random error: horizontal bars are separate
        # from independent-run sampling intervals. A range reaching delta<=0
        # is explicitly open-ended on this log axis.
        lower=x-low if low>1e-12 else x-xmin
        ax.errorbar(x,y,xerr=[[max(0,lower)],[max(0,high-x)]],
                    yerr=[[max(0,y-row['Gamma_ci_low'])],[max(0,row['Gamma_ci_high']-y)]],
                    fmt=MARKERS[Lx],color=color,ecolor=color,ms=4,lw=.8,capsize=2,
                    mfc=color if row['geometry_valid'] else 'white',mew=.8,zorder=5)
        if low<=1e-12:
            crosses=True
            ax.plot(xmin,y,marker='<',color=color,ms=3,clip_on=False)
    ax.set(xscale='log',yscale='log',xlabel=r'$\delta=1-q/q_\star$',ylabel=r'$\Gamma=\overline{I}(q)-\overline{I}(0)$',xlim=(xmin,xmax))
    handles=[Line2D([],[],color='0.55',marker='.',ls='none',label=r'legacy $L_x=1024$')]
    handles += [Line2D([],[],color='0.3',marker=MARKERS[L],ls='none',label=rf'new $L_x={L}$')
                for L in sorted({int(r['Lx']) for r in points})]
    if legacy:
        handles.append(Line2D([],[],color='0.15',ls=':',label=r'legacy $A\delta^{-\psi}$ extrapolation'))
        handles.append(Line2D([],[],color='0.5',ls='--',label=r'legacy $A\delta^{-1/3}+B$ extrapolation'))
    ax.legend(handles=handles,loc='lower left',fontsize=6.2,handlelength=1.5)
    warnings=[]
    if crosses:warnings.append(r'$q_\star$ range reaches $\delta\leq0$')
    if nonpositive:warnings.append('Nonpositive points omitted from log axes; see table')
    if any(not r['geometry_valid'] for r in points):warnings.append('Open symbols: geometry diagnostic failed')
    if warnings:ax.text(.97,.95,'\n'.join(warnings),transform=ax.transAxes,ha='right',va='top',fontsize=5.8)
    ax.text(.03,.95,'(b)',transform=ax.transAxes,ha='left',va='top',
            bbox=dict(facecolor='white',edgecolor='none',pad=.5))
    fig.subplots_adjust(left=.075,right=.99,bottom=.17,top=.97,wspace=.30)
    return fig


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data',type=Path,default=DATA)
    ap.add_argument('--name',default='fig_gpu_wetting',help='Output PDF/PNG basename')
    args=ap.parse_args()
    fig=draw(args.data)
    save(fig,args.name)
    plt.close(fig)

if __name__=='__main__':main()
