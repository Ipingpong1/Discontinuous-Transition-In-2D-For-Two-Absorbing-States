"""Fig. drop: critical droplets at p = 1/2. (a) boundary-growth probability of a disordered disk of radius R
in the dead state (E14_drop05, 500-1000 replicas per point); (b) R_c (P_grow = 1/2) against the
inverse planar-front velocity 1/v(q), with R_c = D/v + a lattice-scale offset.
Only completed fixed-N replica groups are used. P_grow counts boundary hits
by each run's stated T; a still-alive unescaped seed is not a boundary hit."""
import csv
import glob
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import COL, FIGDATA, RUNS, label, save
from campaign.v2.analysis.seed_tails import read_batch, binomial_interval
from campaign.v2.analysis.droplets import critical_radius, radius_bootstrap
from common import ROOT
import json

from common import QSTAR05 as QS


def surv():
    """P_grow(R) per q from E14_drop05; the R = 0.75 (2x2 block, sites alive with probability 1/2) points are
    replaced by the 4000-replica runs of R2_seed22 (2026-10-31) where available, and the full 2x2 block
    of live sites (R2_seed22 'full') is returned separately."""
    rows, full = {}, {}
    audit = []
    for f in sorted(glob.glob(os.path.join(RUNS, "E14_drop05", "out", "*.txt"))):
        m = re.match(r"p0.5_a0.5_j2_q([\d.]+)_R([\d.]+)_L(\d+)_drop", os.path.basename(f))
        batch=read_batch(Path(f))
        n,s,escaped=map(int,batch['values'][-1])
        eligible=batch['completed'] and batch['reached_fixed_max']
        audit.append(dict(path=f,q=float(m[1]),R=float(m[2]),n=n,surviving_or_escaped=s,
                          escaped=escaped,T=int(batch['times'][-1]),eligible=eligible,
                          seed=batch['seed'],source_sha256=batch['sha256'],n_planned=batch['declared_replicas'],
                          adaptive_target=batch['adaptive_target']))
        if not eligible:
            print('Excluded unfinished replica group:',f,n,batch['declared_replicas'])
            continue
        s=escaped  # Common endpoint: macroscopic boundary growth by T.
        lo,hi,_=binomial_interval(s,n)
        rows.setdefault(float(m[1]), []).append((float(m[2]),s/n,s/n-lo,n,hi-s/n))
    for f in sorted(glob.glob(os.path.join(RUNS, "R2_seed22", "out", "*.txt"))):
        m = re.match(r"p0.5_a0.5_j2_q([\d.]+)_R0.75(full|half)_L1024_drop", os.path.basename(f))
        if not m:
            continue
        batch=read_batch(Path(f))
        n,s,escaped=map(int,batch['values'][-1])
        eligible=batch['completed'] and batch['reached_fixed_max']
        audit.append(dict(path=f,q=float(m[1]),R=.75,n=n,surviving_or_escaped=s,
                          escaped=escaped,T=int(batch['times'][-1]),eligible=eligible,
                          seed=batch['seed'],source_sha256=batch['sha256'],n_planned=batch['declared_replicas'],
                          adaptive_target=batch['adaptive_target']))
        if not eligible:continue
        s=escaped  # Common endpoint: macroscopic boundary growth by T.
        lo,hi,_=binomial_interval(s,n)
        x=float(m[1])
        pt=(.75,s/n,s/n-lo,n,hi-s/n)
        if m[2] == "half" and x in rows:
            rows[x] = [r for r in rows[x] if r[0] != 0.75] + [pt]
        elif m[2] == "full":
            full[x] = pt
    from common import ROOT
    out=Path(ROOT)/'campaign/v3/results/droplets';out.mkdir(exist_ok=True,parents=True)
    with (out/'critical_drop_counts.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(audit[0]));writer.writeheader();writer.writerows(audit)
    return {x: np.array(sorted(v)) for x, v in rows.items()}, full


def rc(d):
    return critical_radius(d[:,0],d[:,1])


def front_v(threshold=.3,t0=2000.):
    """Corrected outer-envelope planar velocities, with independent seed samples."""
    out={}
    path=Path(ROOT)/'campaign/v2/results/drop_fronts/runs.csv'
    with path.open() as stream:
        for r in csv.DictReader(stream):
            if (float(r['threshold'])==threshold and float(r['t0'])==t0
                    and r['estimator']=='envelope' and r['done']=='True'
                    and np.isfinite(float(r['v']))):
                out.setdefault(float(r['q']),[]).append(float(r['v']))
    x=np.array(sorted(out));values=[np.asarray(out[q]) for q in x]
    return x,np.array([v.mean() for v in values]),values


if __name__ == "__main__":
    S, FULL = surv()
    xv, vv, v_runs = front_v()
    rng=np.random.default_rng(3829);bootstrap=5000;radius_results=[]
    fig, (a, b) = plt.subplots(1, 2, figsize=(1.0 * 3.4 * 2, 2.2), gridspec_kw=dict(wspace=0.3))
    cmap = plt.get_cmap("viridis")
    xs = sorted(S)
    res = []
    for i, x in enumerate(xs):
        d = S[x]
        c = cmap(i / (len(xs) - 1))
        a.errorbar(d[:, 0], d[:, 1], np.array([d[:,2],d[:,4]]), fmt="o-", color=c, ms=3, capsize=1.2, lw=0.8, label=f"${x/QS:.3f}\\,q^*$")
        if x in FULL:
            a.errorbar([FULL[x][0]], [FULL[x][1]], [[FULL[x][2]],[FULL[x][4]]], fmt="*", color=c, ms=5, capsize=1.2, mec="k", mew=0.3)
        radius_result=radius_bootstrap(d[:,0],d[:,1],d[:,3],bootstrap,rng)
        radius_result["q"]=x;radius_results.append(radius_result)
        R_c = radius_result["radius"]
        v = np.interp(x, xv, vv) if x <= xv.max() else np.nan
        res.append((x, R_c, v))
        print(f"q={x:.7f} x={x/QS:.3f} R_c={R_c:.2f} v={v:.5f} R_c*v={R_c * v:.3f}  n={d[:, 3].min():.0f}-{d[:, 3].max():.0f}")
    a.axhline(0.5, color="0.6", lw=0.6, ls=":")
    a.set(xscale="log", xlabel=r"initial radius $R$", ylabel=r"$P_{\rm grow}$", ylim=(-0.02, 1.02))
    a.legend(fontsize=5.5, loc="upper left", title=r"$q$", title_fontsize=6)
    label(a, "(a)", x=.9, y=.12)
    res = np.array(res)
    ok = np.isfinite(res[:, 2])
    inv = 1 / res[ok, 2]
    A = np.vstack([inv, np.ones_like(inv)]).T
    (D, dR), *_ = np.linalg.lstsq(A, res[ok, 1], rcond=None)
    print(f"fit R_c = D/v + dR: D = {D:.3f}, dR = {dR:.2f}")
    lower=np.array([r['ci_low'] for r in radius_results])[ok]
    upper=np.array([r['ci_high'] for r in radius_results])[ok]
    b.errorbar(inv,res[ok,1],yerr=[res[ok,1]-lower,upper-res[ok,1]],fmt='ko',ms=4,mfc='w',capsize=2,lw=.7)
    # Each bootstrap draw uses one shared velocity grid for every interpolated q.
    velocity_draws=np.column_stack([r[rng.integers(0,len(r),(bootstrap,len(r)))].mean(1) for r in v_runs])
    fitted=[]
    for draw in range(bootstrap):
        radius=np.array([r['samples'][draw] for r in radius_results])[ok]
        velocity=np.interp(res[ok,0],xv,velocity_draws[draw])
        if np.all(np.isfinite(radius)) and np.all(velocity>0):
            fitted.append(np.linalg.lstsq(np.column_stack([1/velocity,np.ones(len(radius))]),radius,rcond=None)[0])
    ci=np.quantile(fitted,[.025,.975],axis=0)
    sensitivity=[]
    for threshold in (.2,.3,.5):
        for t0 in (1000.,2000.,3000.):
            xx,vel,_=front_v(threshold,t0)
            vvq=np.interp(res[ok,0],xx,vel)
            coef=np.linalg.lstsq(np.column_stack([1/vvq,np.ones(len(vvq))]),res[ok,1],rcond=None)[0]
            sensitivity.append(dict(threshold=threshold,t0=t0,D=float(coef[0]),offset=float(coef[1])))
    analysis=dict(endpoint='boundary growth before extinction by the protocol-specific T',
        radius_interpolation='linear logit P in log R',
        radius_results=[{k:val for k,val in r.items() if k!='samples'} for r in radius_results],
        curvature=dict(D=float(D),offset=float(dR),D_ci95=ci[:,0].tolist(),offset_ci95=ci[:,1].tolist(),
            bootstrap_valid=len(fitted),bootstrap_requested=bootstrap,
            interval_scope='conditional binomial replica and independent planar-run sampling; excludes protocol, interpolation, tracker, and curvature-model systematics'),
        tracker_sensitivity=sensitivity)
    (Path(ROOT)/'campaign/v3/results/droplets/critical_radius.json').write_text(json.dumps(analysis,indent=2)+'\n')
    print('Conditional bootstrap D interval:',ci[:,0],'; tracker/window range:',min(r['D'] for r in sensitivity),max(r['D'] for r in sensitivity))
    for (x, R_c, v) in res[ok]:
        b.annotate(f"{x/QS:.3f}", (1 / v, R_c), textcoords="offset points", xytext=(-8 if 1 / v > 150 else 4, -9), fontsize=6, ha="right" if 1 / v > 150 else "left")
    g = np.linspace(0, inv.max() * 1.05, 50)
    b.plot(g, D * g + dR, "C3-", lw=0.9, label=f"$R_c = D/v + {dR:.1f}$, $D={D:.2f}$")
    b.set(xlabel=r"$1/v(q)$ (steps per site)", ylabel=r"critical radius $R_c$", xlim=(0, None), ylim=(0, None))
    b.legend(fontsize=6.5, loc="upper left")
    label(b, "(b)", x=0.9, y=0.12)
    save(fig, "fig_drop")
