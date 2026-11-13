#!/usr/bin/env python3
"""Seed-averaged tables per parameter set for strip (front) and bulk experiments.

    python3 fronts.py $SCRATCH/dn/runs/E2_bulk $SCRATCH/dn/runs/E1_front ...

Experiment dirs whose runs have profiles (.prof) get front tracking; the others get
bulk rho(t) tables with the half-life of the disordered state.

Front tracking: the strip occupies rows [0, Ly/2) at t=0.  Roll so the strip centre
sits at Ly/2, smooth the per-row interface density over 32 rows, and walk outward from
the centre until it drops below half the bulk value rho_b (taken from the strip
interior at the first record with t >= 1000).  W(t) = width between the two edges.
v = dW/dt / 2 fitted over records with W in [0.15, 0.85] Ly, interior still disordered
(centre density > 0.7 rho_b) and t >= 2000.  t_order = first time the strip centre
orders (density < 0.5 rho_b).
"""
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from profile import read  # noqa: E402

KEY = re.compile(r"(?:^|_)(?P<k>p|a|r|j|q|Lx|Ly|L|R|s)(?P<v>[0-9.]+)")


def params(name):
    return {m["k"]: m["v"] for m in KEY.finditer(name)}


def group_key(name):
    p = params(name)
    return tuple((k, float(v)) for k, v in sorted(p.items()) if k != "s")


def fmt_key(key):
    return " ".join(f"{k}={v:g}" for k, v in key)


def load_ts(path):
    return np.array([[float(x) for x in l.split()] for l in open(path) if not l.startswith("#")])


def bulk(exp):
    print(f"=== {exp}: seed-mean rho(t) (bond density) ===")
    g = defaultdict(list)
    for p in glob.glob(f"{exp}/out/*.txt"):
        g[group_key(os.path.basename(p)[:-4])].append(load_ts(p))
    tmax = max(ts[-1, 0] for runs in g.values() for ts in runs)
    times = [t for t in (1e3, 3e3, 1e4, 2e4, 3e4, 5e4, 1e5, 2e5, 5e5, 1e6) if t <= tmax]
    if tmax > 1e5:
        times = [t for t in times if t not in (2e4, 5e4)]
    print("params".ljust(37) + " ".join(f"{t:>7.0e}" for t in times) + "   tau_half")
    for key in sorted(g):
        rows, taus = [], []
        for ts in g[key]:
            rows.append([np.interp(t, ts[:, 0], ts[:, 2]) for t in times])
            r0 = np.interp(1000, ts[:, 0], ts[:, 2])
            below = ts[(ts[:, 0] > 1000) & (ts[:, 2] < 0.5 * r0)]
            taus.append(below[0, 0] if len(below) else np.inf)
        r = np.mean(rows, 0)
        tt = " ".join(f"{x:.0f}" if np.isfinite(x) else f">{g[key][0][-1, 0]:.0e}" for x in taus)
        print(fmt_key(key).ljust(36) + " " + " ".join(f"{x:7.4f}" for x in r) + f"   {tt}")


def track(txt, prof):
    Lx, Ly, t, alive, iface = read(txt, prof)
    d = np.roll(iface / Lx, Ly // 2 - Ly // 4, axis=1)  # strip centre -> Ly/2
    k = np.ones(32) / 32
    ds = np.array([np.convolve(np.concatenate([r[-16:], r, r[:15]]), k, "valid") for r in d])
    i0 = min(np.searchsorted(t, 1000), len(t) - 1)
    c = Ly // 2
    rho_b = d[i0, c - Ly // 8: c + Ly // 8].mean()
    W, cen = [], []
    for row in ds:
        lo = row[:c][::-1] < 0.5 * rho_b
        hi = row[c:] < 0.5 * rho_b
        yl = np.argmax(lo) if lo.any() else c
        yr = np.argmax(hi) if hi.any() else Ly - c
        W.append(yl + yr)
        cen.append(row[c - 128: c + 128].mean() / rho_b)
    return t, np.array(W, float), np.array(cen), Ly, rho_b


def front_v(txt):
    t, W, cen, Ly, rho_b = track(txt, txt[:-4] + ".prof")
    sel = (t >= 2000) & (W > 0.15 * Ly) & (W < 0.85 * Ly) & (cen > 0.7)
    v = np.polyfit(t[sel], W[sel], 1)[0] / 2 if sel.sum() >= 4 else np.nan
    nuc = t[(cen < 0.5) & (t > 1000)]
    return (v, sel.sum(), t[sel][-1] if sel.any() else 0, nuc[0] if len(nuc) else np.inf,
            W[-1] / Ly, rho_b)


def fronts(exp):
    print(f"=== {exp}: front tracking ===")
    g = defaultdict(list)
    for p in sorted(glob.glob(f"{exp}/out/*.txt")):
        g[group_key(os.path.basename(p)[:-4])].append(front_v(p))
    print("params".ljust(37) + "v_mean    v_seeds                     npts  t_fitmax  t_order  W_end/Ly  rho_b")
    for key in sorted(g):
        rs = g[key]
        vs = np.array([r[0] for r in rs])
        fin = vs[np.isfinite(vs)]
        vm = f"{fin.mean():+.5f}" if len(fin) else "   nan  "
        print(fmt_key(key).ljust(36) + " " + f"{vm}  "
              + " ".join(f"{x:+.5f}" for x in vs).ljust(28)
              + f" {'/'.join(str(r[1]) for r in rs):8s} {'/'.join(f'{r[2]:.0f}' for r in rs):12s}"
              + f" {'/'.join(f'{r[3]:.0f}' if np.isfinite(r[3]) else '-' for r in rs):14s}"
              + f" {'/'.join(f'{r[4]:.2f}' for r in rs):14s} {rs[0][5]:.3f}")


if __name__ == "__main__":
    for exp in sys.argv[1:]:
        exp = exp.rstrip("/")
        (fronts if glob.glob(f"{exp}/out/*.prof") else bulk)(exp)
