"""Wetting fits with comparable residuals, common-baseline covariance and q* sensitivity.

Both candidate models have two free parameters and are fitted in original Gamma
units. Bootstrap units are independent runs; q* is a separately reported systematic
sensitivity, never absorbed into an invented statistical error distribution.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
import re

import numpy as np

from .io import ROOT, provenance, read_header, write_csv, write_json
from .stats import window_mean

NAME = re.compile(r"p([\d.]+)_a([\d.]+)_j2_q([\d.]+)_Lx(\d+)_Ly(\d+)_s(\d+)")


def collect(figdata=None, experiments=("E16_wall", "R2_wall"), min_time=300000., late_fraction=.5):
    figdata = ROOT / "data/figdata" if figdata is None else Path(figdata)
    # A continued R2 run replaces its E16 predecessor; it is not a new replica.
    paths = {}
    for exp in experiments:
        for path in sorted((figdata / exp).glob("*.csv")):
            if not path.stem.endswith("_prof"):
                paths[path.name] = path
    rows = []
    header_inputs = []
    for path in paths.values():
        match = NAME.match(path.stem)
        if not match:
            continue
        p, alpha, q = map(float, match.groups()[:3])
        Lx, Ly, seed = map(int, match.groups()[3:])
        txt = ROOT / "data/runs" / path.parent.name / "out" / (path.stem + ".txt")
        parameter_source = "filename_fallback"
        source_complete = False
        source_end = None
        if txt.exists():
            hdr = read_header(txt)
            header_inputs.append(txt)
            source_complete = hdr["done"]
            source_end = hdr.get("t_end")
            q, parameter_source = hdr.get("q", q), hdr.get("parameter_source", parameter_source)
        d = np.genfromtxt(path, delimiter=",", names=True)
        t = d["t"]
        T = float(t[-1]); start = max(min_time, T * late_fraction)
        drift = np.flatnonzero(abs(d["alive"] - .5) >= .2)
        end = float(t[drift[0] - 1]) if len(drift) and drift[0] else (0. if len(drift) else T)
        sel = (t >= start) & (t <= end)
        source_consistent = bool(source_complete and source_end == T)
        eligible = source_consistent and T >= 9e5 and end > start and sel.sum() >= 3
        I = window_mean(t, d["I"], start, end) if eligible else np.nan
        W = window_mean(t, d["W"], start, end) if eligible else np.nan
        mid = (start + end) / 2
        change = (window_mean(t, d["I"], mid, end) - window_mean(t, d["I"], start, mid)) if eligible else np.nan
        rows.append(dict(path=str(path.resolve()), run=path.stem, p=p, alpha=alpha, q=q,
                         Lx=Lx, Ly=Ly, seed=seed, parameter_source=parameter_source,
                         T=T, start=start, end=end, truncated=bool(len(drift)),
                         source_complete=source_complete, source_end=source_end, source_consistent=source_consistent,
                         eligible=eligible, n_records=int(sel.sum()), I=I, W=W,
                         I_drift=change, retained_fraction=max(0., end-start) / max(T-start, 1.)))
    return rows, list(paths.values()) + header_inputs


def grouped(rows):
    g = defaultdict(list)
    for r in rows:
        if r["eligible"]:
            g[(r["p"], r["q"], r["Lx"], r["Ly"])].append(r)
    return g


def point_table(rows, qstar=.04596, p=.5):
    g = grouped(rows)
    baselines = {Ly: g[(p, 0., 1024, Ly)] for Ly in (1024, 2048, 4096) if (p, 0., 1024, Ly) in g}
    if 1024 not in baselines:
        raise ValueError("missing q=0 wall reference")
    points = []
    for (pp, q, Lx, Ly), rs in sorted(g.items()):
        if pp != p or q <= 0:
            continue
        base_Ly = Ly if Ly in baselines else 1024
        bv = np.array([r["I"] for r in baselines[base_Ly]])
        iv = np.array([r["I"] for r in rs])
        isem = iv.std(ddof=1) / np.sqrt(len(iv)) if len(iv)>1 else np.nan
        bsem = bv.std(ddof=1) / np.sqrt(len(bv)) if len(bv)>1 else np.nan
        points.append(dict(p=p, q=q, delta=1-q/qstar, Lx=Lx, Ly=Ly, n=len(rs),
                           I=float(iv.mean()), I_sem=float(isem), baseline_Ly=base_Ly,
                           I0=float(bv.mean()), I0_sem=float(bsem), Gamma=float(iv.mean()-bv.mean()),
                           Gamma_sem=float(np.hypot(isem,bsem)), W=float(np.mean([r["W"] for r in rs])),
                           n_truncated=sum(r["truncated"] for r in rs),
                           I_drift=float(np.mean([r["I_drift"] for r in rs]))))
    return points, g, baselines


def golden_minimum(fn, lo=0., hi=2., steps=70):
    ratio = (np.sqrt(5)-1)/2
    a,b = hi-ratio*(hi-lo), lo+ratio*(hi-lo)
    fa,fb = fn(a),fn(b)
    for _ in range(steps):
        if fa < fb:
            hi,b,fb=b,a,fa; a=hi-ratio*(hi-lo); fa=fn(a)
        else:
            lo,a,fa=a,b,fb; b=lo+ratio*(hi-lo); fb=fn(b)
    return (lo+hi)/2


def fit_models(delta, gamma, covariance):
    delta,gamma = np.asarray(delta,float),np.asarray(gamma,float)
    if np.any(delta<=0) or len(delta)<3:
        raise ValueError("positive delta and at least three points required")
    precision = np.linalg.inv(covariance)
    def power(psi):
        z=delta**(-psi); A=float(z@precision@gamma/(z@precision@z))
        res=gamma-A*z
        return float(res@precision@res),A
    psi=golden_minimum(lambda x:power(x)[0]); chi,A=power(psi)
    X=np.column_stack([delta**(-1/3),np.ones(len(delta))])
    coeff=np.linalg.solve(X.T@precision@X,X.T@precision@gamma)
    res=gamma-X@coeff; chi3=float(res@precision@res)
    return dict(power=dict(A=A,psi=psi,chi2=chi,dof=len(delta)-2,
                           max_relative_residual=float(np.max(abs(gamma-A*delta**(-psi))/gamma))),
                third_offset=dict(A=float(coeff[0]),B=float(coeff[1]),chi2=chi3,dof=len(delta)-2,
                                  max_relative_residual=float(np.max(abs(res)/gamma))))


def analyze(rows,qstar=.04596,qstar_half_range=.00005,boot=2000):
    points,g,baselines=point_table(rows,qstar)
    fitpts=[r for r in points if r["Lx"]==1024 and r["Ly"]>=2048 and 0<r["delta"]<=.1001
            and r["n"]>=2 and np.isfinite(r["Gamma_sem"])]
    fitpts.sort(key=lambda r:(r["q"],r["Ly"]))
    if len(fitpts)<4:
        return points,dict(status="insufficient eligible fit points")
    qs=np.array([r["q"] for r in fitpts]); gamma=np.array([r["Gamma"] for r in fitpts])
    C=np.diag([max(r["I_sem"]**2,1e-16) for r in fitpts])
    for i,a in enumerate(fitpts):
        for j,b in enumerate(fitpts):
            if a["baseline_Ly"]==b["baseline_Ly"]:
                C[i,j]+=a["I0_sem"]**2
    central=fit_models(1-qs/qstar,gamma,C)
    rng=np.random.default_rng(1729); values=[]
    for _ in range(boot):
        base={ly:rng.choice([r["I"] for r in rs],len(rs),replace=True).mean() for ly,rs in baselines.items()}
        y=[]
        for r in fitpts:
            rs=g[(r["p"],r["q"],r["Lx"],r["Ly"])]
            y.append(rng.choice([a["I"] for a in rs],len(rs),replace=True).mean()-base[r["baseline_Ly"]])
        f=fit_models(1-qs/qstar,y,C)
        values.append([f["power"]["A"],f["power"]["psi"],f["third_offset"]["A"],f["third_offset"]["B"]])
    if values:
        lo,hi=np.quantile(values,[.025,.975],axis=0)
        for j,(model,param) in enumerate((("power","A"),("power","psi"),("third_offset","A"),("third_offset","B"))):
            central[model][param+"_ci_low"]=float(lo[j]);central[model][param+"_ci_high"]=float(hi[j])
    sensitivity=[]
    qgrid=np.unique(np.r_[np.linspace(qstar-qstar_half_range,qstar+qstar_half_range,11),
                          qstar-qstar_half_range/2,qstar+qstar_half_range/2])
    for q0 in qgrid:
        if np.any(qs>=q0):
            continue
        fit=fit_models(1-qs/q0,gamma,C)
        sensitivity.append(dict(qstar=float(q0),**fit))
    return points,dict(status="ok",qstar=qstar,qstar_half_range=qstar_half_range,
                       central=central,qstar_sensitivity=sensitivity,fit_points=fitpts,
                       gamma_covariance=C.tolist(),bootstrap=boot,
                       notes=["Original-Gamma generalized least squares for both two-parameter models.",
                              "Common q=0 reference covariance retained across points sharing that reference.",
                              "95% independent-run bootstrap intervals conditional on qstar; qstar range is systematic sensitivity.",
                              "Window begins at max(3e5,T/2); truncated before first |alive-.5|>=.2 record."])


def main():
    reference=ROOT/'campaign/v2/results/coexistence_reference.json'
    shared=json.loads(reference.read_text())
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figdata",type=Path,default=ROOT/"data/figdata")
    ap.add_argument("--out",type=Path,default=ROOT/"campaign/v2/results/wetting")
    ap.add_argument("--qstar",type=float,default=shared["qstar"])
    ap.add_argument("--qstar-half-range",type=float,default=shared["systematic_half_range"])
    ap.add_argument("--boot",type=int,default=2000)
    args=ap.parse_args()
    rows,paths=collect(args.figdata)
    paths += [reference,Path(__file__),Path(__file__).with_name("io.py"),Path(__file__).with_name("stats.py")]
    points,result=analyze(rows,args.qstar,args.qstar_half_range,args.boot)
    write_csv(args.out/"runs.csv",rows);write_csv(args.out/"points.csv",points)
    write_json(args.out/"analysis.json",dict(result=result,inputs=provenance(paths)))
    print(f"{len(rows)} runs, {len(points)} points; {args.out}")
    if result["status"]=="ok":
        print(result["central"])


if __name__=="__main__":
    main()
