"""V2 wetting: original-space fits and explicit uncertainty of coexistence point."""
import csv
import json
from pathlib import Path
import re
from common import ROOT, FIGDATA, PAGE, QSTAR05, label, save
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

RESULTS=Path(ROOT)/'campaign/v2/results/wetting'
if not (RESULTS/'analysis.json').exists():
    raise SystemExit('Run python3 -m campaign.v2.analysis.wetting before making this figure.')
analysis=json.loads((RESULTS/'analysis.json').read_text())['result']
qstar=analysis['qstar'];dq=analysis['qstar_half_range']
if abs(qstar-QSTAR05)>1e-12:
    raise SystemExit('Re-run the wetting analysis with the shared coexistence reference before plotting.')
points=list(csv.DictReader((RESULTS/'points.csv').open()))
for r in points:
    for k in r:r[k]=float(r[k])
fig,axs=plt.subplots(1,3,figsize=(PAGE,2.3),gridspec_kw=dict(wspace=.32))
groups={}
for path in sorted((Path(FIGDATA)/'E16_wall').glob('p0.5*_Lx1024_Ly1024_s*.csv')):
    if path.stem.endswith('_prof'):continue
    q=float(re.search(r'_q([\d.]+)_',path.stem)[1])
    groups.setdefault(q,[]).append(path)
for i,(q,paths) in enumerate(sorted(groups.items())):
    color='k' if q==0 else plt.get_cmap('plasma')(.85*i/max(len(groups)-1,1))
    series=[np.genfromtxt(p,delimiter=',',names=True) for p in paths]
    t=series[0]['t'];y=np.mean([np.interp(t,s['t'],s['I']) for s in series],axis=0)
    axs[0].plot(t[t>=100],y[t>=100],color=color,lw=.8,label='$0$' if q==0 else f'${q/qstar:.3f}$')
    profiles=[np.genfromtxt(p.with_name(p.stem+'_prof.csv'),delimiter=',',names=True) for p in paths]
    axs[1].plot(profiles[0]['y'],np.mean([s['iface_late_centred'] for s in profiles],axis=0),color=color,lw=.8)
axs[0].set(xscale='log',xlabel='$t$',ylabel=r'interface sites per wall length $I$')
axs[0].legend(fontsize=5.2,title=r'$q/q^*$',title_fontsize=6,loc='upper left',ncol=2)
axs[1].set(xlabel='distance from layer center',ylabel='interface-site fraction',xlim=(-160,160))
mk={256:'v',1024:'o',4096:'s'};colors={1024:'k',2048:'C0',4096:'C3'}
ax=axs[2]
for Lx in mk:
    for Ly in colors:
        rr=sorted((r for r in points if r['Lx']==Lx and r['Ly']==Ly and r['delta']>0
                   and not (Lx==4096 and Ly==1024 and r['delta']<.03)),key=lambda r:r['delta'])
        if not rr:continue
        xx=np.array([r['delta'] for r in rr]);yy=np.array([r['Gamma'] for r in rr]);se=np.array([r['Gamma_sem'] for r in rr])
        ax.errorbar(xx,yy,se,fmt=mk[Lx],color=colors[Ly],mfc='w' if Ly==1024 else colors[Ly],ms=3,capsize=1,lw=.7)
        if Lx==1024 and Ly>=2048:
            q=np.array([r['q'] for r in rr]);lo=1-q/(qstar-dq);hi=1-q/(qstar+dq)
            ax.errorbar(xx,yy,xerr=[xx-lo,hi-xx],fmt='none',ecolor=colors[Ly],lw=.5,alpha=.5)
xx=np.geomspace(min(r['delta'] for r in points if r['delta']>0)*.8,.12,100)
power=analysis['central']['power'];third=analysis['central']['third_offset']
ax.plot(xx,power['A']*xx**(-power['psi']),'k--',lw=.8,label=fr'power fit, $\psi={power["psi"]:.2f}$')
ax.plot(xx,third['A']*xx**(-1/3)+third['B'],'k:',lw=1,label=r'$A\delta^{-1/3}+B$')
ax.set(xscale='log',yscale='log',xlabel=r'$\delta=1-q/q^*$',ylabel=r'excess $\Gamma=I-I_0$')
handles=[Line2D([],[],marker=mk[L],ls='',color='.3',ms=3,label=f'$L_x={L}$') for L in mk]
handles += [Line2D([],[],marker='o',ls='',color=colors[L],mfc='w' if L==1024 else colors[L],ms=3,label=f'$L_y={L}$') for L in colors]
handles += [l for l in ax.get_lines() if l.get_label() in (fr'power fit, $\psi={power["psi"]:.2f}$',r'$A\delta^{-1/3}+B$')]
ax.legend(handles=handles,fontsize=4.8,ncol=2,loc='lower left',handlelength=1.3,columnspacing=.5)
ax.set_title(fr'$q^*={qstar:.5f}$; sensitivity $\pm{dq:.5f}$',fontsize=6)
for a,s in zip(axs,'abc'):label(a,f'({s})',x=.86)
save(fig,'fig_wall05')
