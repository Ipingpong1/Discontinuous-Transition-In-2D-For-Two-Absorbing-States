#!/usr/bin/env python3
"""Orientation test with direct front tracking (R2_orient, R2_orient2k; dnsim v1.5 --pdir a,b).

    python3 orient2.py EXP_DIR [--t0 10000] [--csv OUT.csv] [--boot 2000]

The profile of a run is binned along the band normal, s = (a y - b x) mod L, L bins of L cells each,
1/sqrt(a^2+b^2) apart in normal distance.  The band starts at s in [0, L/2); the profile is rolled by
L/4 so that it sits in the middle.  Per record the band width is measured in two ways:
  outer   : outermost bins (walking out from the band's interface-weighted centre) whose 9-bin running
            mean of the interface density exceeds 0.3 x iota_ref, iota_ref the interior interface
            fraction (central 20% of the band) at t ~ 1e3 (the estimator of the strip campaign);
  integ   : (sum over bins of the interface count) / (L iota_in(t)), iota_in(t) the central-20% value
            at that time (the band's interior relaxes, but the width uses its instantaneous density).
Normal width W = width in bins / sqrt(a^2 + b^2); the normal front velocity is v = (dW/dt)/2 by least
squares over t >= t0 while 0.15 L < width < 0.85 L (bins).  Per (p, a, b): seed mean and s.e.m. of v
at each q, the zero crossing x* of the weighted straight line through the points (its error from the
fit and from a seed bootstrap), chi^2 of the line.
"""
import argparse
import glob
import math
import os
import re
import sys
from collections import defaultdict

import numpy as np

NAME = re.compile(r"p([\d.]+)_a0\.5_q([\d.]+)_ab(\d)-(\d)_L(\d+)_s(\d+)$")
QS = {0.1: 0.03680, 0.5: 0.04596, 0.9: 0.01194}


def read(prefix):
    L = None
    for line in open(prefix + ".txt"):
        if line.startswith("# Lx"):
            L = int(line.split()[2])
            break
    raw = np.fromfile(prefix + ".prof", dtype=np.uint8)
    rec = 8 + 8 * L
    n = raw.size // rec
    raw = raw[: n * rec].reshape(n, rec)
    t = raw[:, :8].copy().view(np.int64).ravel()
    alive = raw[:, 8: 8 + 4 * L].copy().view(np.uint32).astype(float)
    iface = raw[:, 8 + 4 * L:].copy().view(np.uint32).astype(float)
    return L, t, alive, iface


def widths(L, t, iface, frac=0.5):
    """Per record: (outer width, integrated width) in bins, and the interior interface fraction."""
    d = np.roll(iface, L // 4, axis=1) / L            # interface fraction per bin, band centred at L/2
    n = len(t)
    c0 = L // 4 + int(frac * L) // 2                  # initial centre bin
    half = int(0.1 * frac * L)
    iota_in = d[:, c0 - half: c0 + half].mean(1)
    i_ref = np.interp(1000.0, t, iota_in)
    s = np.arange(L)
    W_out = np.zeros(n)
    sm = np.apply_along_axis(lambda r: np.convolve(r, np.ones(9) / 9, "same"), 1, d)
    for i in range(n):
        tot = d[i].sum()
        if tot <= 0:
            continue
        c = int(round((d[i] * s).sum() / tot))
        c = min(max(c, 1), L - 2)
        thr = 0.3 * i_ref
        lo = c
        while lo > 0 and sm[i, lo - 1] > thr:
            lo -= 1
        hi = c
        while hi < L - 1 and sm[i, hi + 1] > thr:
            hi += 1
        W_out[i] = hi - lo + 1 if sm[i, c] > thr else 0
    with np.errstate(invalid="ignore", divide="ignore"):
        W_int = np.where(iota_in > 0, d.sum(1) / iota_in, 0.0)
    return W_out, W_int, iota_in


def velocity(t, W, L, ab_norm, t0):
    sel = (t >= t0) & (W > 0.15 * L) & (W < 0.85 * L)
    bad = np.flatnonzero((t >= t0) & ~sel)
    if bad.size:
        sel &= t < t[bad[0]]
    if sel.sum() < 5:
        return np.nan, 0, 0
    slope = np.polyfit(t[sel], W[sel] / ab_norm, 1)[0]
    return slope / 2, t[sel][-1], int(sel.sum())


def crossing(x, v, se):
    sel = np.isfinite(v)
    se = np.where(np.isfinite(se) & (se > 0), se, np.nan)
    w = 1 / se[sel] ** 2 if np.isfinite(se[sel]).all() else np.ones(sel.sum())   # single seed: equal weights
    A = np.vstack([np.ones(sel.sum()), x[sel] - 1]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    a, b = cov @ (A.T @ (w * v[sel]))
    xs = 1 - a / b
    J = np.array([-1 / b, a / b ** 2])
    chi2 = float((w * (v[sel] - a - b * (x[sel] - 1)) ** 2).sum())
    return xs, float(np.sqrt(J @ cov @ J)), b, chi2, int(sel.sum()) - 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exp")
    ap.add_argument("--t0", type=float, default=10000)
    ap.add_argument("--csv")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--est", default="outer", choices=["outer", "integ"])
    a = ap.parse_args()
    rows = []
    for f in sorted(glob.glob(os.path.join(a.exp, "out", "*.txt"))):
        if not any(ln.startswith("# done") for ln in open(f)):
            continue
        prefix = f[:-4]
        m = NAME.search(os.path.basename(prefix))
        if not m or not os.path.exists(prefix + ".prof"):
            continue
        p, q, aa, bb, L, s = float(m[1]), float(m[2]), int(m[3]), int(m[4]), int(m[5]), int(m[6])
        L_, t, alive, iface = read(prefix)
        W_out, W_int, iota = widths(L, t, iface)
        nrm = math.hypot(aa, bb)
        vo, te_o, n_o = velocity(t, W_out, L, nrm, a.t0)
        vi, te_i, n_i = velocity(t, W_int, L, nrm, a.t0)
        rows.append(dict(p=p, q=q, x=q / QS[p], ab=f"{aa}-{bb}", L=L, s=s, v_out=vo, v_int=vi, t_end_fit=te_o,
                         npts=n_o, iota_in=float(iota[np.argmin(np.abs(t - 1e4))]),
                         W0=W_out[0] / nrm, W_end=W_out[-1] / nrm))
    if a.csv:
        import csv
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    key = "v_" + ("out" if a.est == "outer" else "int")
    g = defaultdict(list)
    for r in rows:
        g[(r["p"], r["L"], r["ab"], round(r["x"], 4))].append(r[key])
    print(f"estimator {a.est}, fit from t0 = {a.t0:g}")
    print("p     L     ab    x        v_mean      sem       n   seeds")
    rng = np.random.default_rng(1)
    for (p, L, ab), _ in sorted({(k[0], k[1], k[2]): 0 for k in g}.items()):
        xs = sorted(k[3] for k in g if k[:3] == (p, L, ab))
        vm, se, vlists = [], [], []
        for x in xs:
            vs = np.array(g[(p, L, ab, x)])
            fin = vs[np.isfinite(vs)]
            vm.append(fin.mean() if fin.size else np.nan)
            se.append(fin.std(ddof=1) / math.sqrt(fin.size) if fin.size > 1 else np.nan)
            vlists.append(fin)
            print(f"{p:<5} {L:<5} {ab:<5} {x:<8.4f} {vm[-1]:+.6f}  {se[-1]:.6f}  {fin.size:<3d} "
                  + " ".join(f"{v:+.5f}" for v in vs))
        xs, vm, se = np.array(xs), np.array(vm), np.array(se)
        if np.isfinite(vm).sum() >= 3:
            xstar, dx, slope, chi2, dof = crossing(xs, vm, se)
            boots = []
            for _ in range(a.boot):
                bv = np.array([rng.choice(v, size=v.size).mean() if v.size else np.nan for v in vlists])
                try:
                    boots.append(crossing(xs, bv, se)[0])
                except np.linalg.LinAlgError:
                    pass
            boots = np.array(boots)
            print(f"  -> p={p} L={L} ab={ab}: x* = {xstar:.5f} +- {dx:.5f} (fit), +- {boots.std():.5f} (seed bootstrap, "
                  f"{len(boots)} resamples); q* dv/dq = {slope:.4f}; chi2/dof = {chi2:.1f}/{dof}")


if __name__ == "__main__":
    sys.exit(main())
