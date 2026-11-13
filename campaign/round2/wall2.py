#!/usr/bin/env python3
"""Wetting layer at p = 1/2, round 2 (R2_wall; wall.py series): stationarity, the excess adsorption
Gamma = I - I_0 against delta = 1 - q/q*, and fits.

    python3 wall2.py EXP_DIR [--e16 E16_DIR] [--qs 0.04596] [--csv OUT.csv]

Per run: I(t) (interface cells per unit wall length, mean of the two walls), W(t), rc(t) from wall.py;
window means of I over [3e5, 1e6], [1e6, 2e6], [2e6, 3e6] (stationarity), and the late window
[T/2, T] used for Gamma.  I_0 from the bare-wall runs (q = 0) of the same Ly (E16 Ly 1024: 8.7).
Fits of Gamma(delta) for delta <= 0.1: A delta^-psi; A delta^-psi + B; A delta^-1/3 + B."""
import argparse
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, "/home/mo2273/Desktop/dn/campaign")
sys.path.insert(0, "/scratch/mo2273/dn/campaign2")
import wall  # noqa: E402

NAME = re.compile(r"p0\.5_a0\.5_j2_q([\d.]+)_Lx(\d+)_Ly(\d+)_s(\d+)$")


def window(t, y, a, b):
    sel = (t >= a) & (t <= b)
    return y[sel].mean() if sel.sum() >= 3 else np.nan


def collect(dirs):
    runs = {}
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "out", "*.prof"))):
            prefix = f[:-5]
            m = NAME.search(os.path.basename(prefix))
            if not m or not os.path.exists(prefix + ".txt"):
                continue
            done = any(l.startswith("# done") for l in open(prefix + ".txt"))
            t, I, sig, rows, rin, al, W, rc = wall.run_series(prefix)
            runs[(float(m[1]), int(m[2]), int(m[3]), int(m[4]), os.path.basename(d))] = dict(
                t=t, I=I, W=W, rc=rc, T=t[-1], done=done)
    return runs


def fit_power(d, G, e):
    w = 1 / np.maximum(e, 0.3) ** 2
    A = np.vstack([np.ones(len(d)), -np.log(d)]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    lnA, psi = cov @ (A.T @ (w * np.log(G)))
    chi2 = (w * (np.log(G) - lnA + psi * np.log(d)) ** 2).sum()
    return psi, np.sqrt(cov[1, 1]), chi2


def fit_offset(d, G, e, psi=None):
    """Gamma = A d^-psi + B; psi free (grid) or fixed."""
    w = 1 / np.maximum(e, 0.3) ** 2
    def lin(ps):
        X = np.vstack([d ** -ps, np.ones(len(d))]).T
        cov = np.linalg.inv(X.T @ (X * w[:, None]))
        A, B = cov @ (X.T @ (w * G))
        return (w * (G - A * d ** -ps - B) ** 2).sum(), A, B
    if psi is not None:
        c, A, B = lin(psi)
        return psi, A, B, c
    grid = np.linspace(0.15, 1.2, 211)
    res = [lin(ps) for ps in grid]
    k = int(np.argmin([r[0] for r in res]))
    return grid[k], res[k][1], res[k][2], res[k][0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exp")
    ap.add_argument("--e16")
    ap.add_argument("--qs", type=float, default=0.04596)
    ap.add_argument("--csv")
    a = ap.parse_args()
    runs = collect([a.exp] + ([a.e16] if a.e16 else []))
    print("run                                             T      done  I[3e5,1e6] I[1e6,2e6] I[2e6,3e6]  I_late    W_late  rc_late")
    g = defaultdict(list)
    for (q, Lx, Ly, s, exp), r in sorted(runs.items()):
        t, I = r["t"], r["I"]
        late = (t >= r["T"] / 2)
        w1, w2, w3 = window(t, I, 3e5, 1e6), window(t, I, 1e6, 2e6), window(t, I, 2e6, 3e6)
        Il, Wl, rl = I[late].mean(), r["W"][late].mean(), r["rc"][late].mean()
        print(f"{exp[:8]:<8s} q{q:<10} Lx{Lx:<5d} Ly{Ly:<5d} s{s} {r['T']:8.0f}  {int(r['done'])}   {w1:8.2f}  {w2:8.2f}  {w3:8.2f}  "
              f"{Il:8.2f}  {Wl:7.1f}  {rl:6.3f}")
        if r["done"] or r["T"] >= 9e5:
            g[(q, Lx, Ly)].append((Il, Wl, r["T"]))
    I0 = {}
    for (q, Lx, Ly), v in g.items():
        if q == 0:
            I0[Ly] = np.mean([x[0] for x in v])
    I0.setdefault(1024, 8.7)
    print("\nbare wall I_0:", {k: round(v, 2) for k, v in I0.items()})
    print("\nq          delta     Lx    Ly    n  Gamma      err     W       T_min")
    pts = defaultdict(list)
    for (q, Lx, Ly), v in sorted(g.items()):
        if q == 0:
            continue
        Il = np.array([x[0] for x in v])
        i0 = I0.get(Ly, I0[1024])
        G, e = Il.mean() - i0, Il.std(ddof=1) / np.sqrt(len(Il)) if len(Il) > 1 else np.nan
        d = 1 - q / a.qs
        print(f"{q:<10} {d:<9.4f} {Lx:<5d} {Ly:<5d} {len(v):<2d} {G:8.2f}  {e:6.2f}  {np.mean([x[1] for x in v]):6.1f}  {min(x[2] for x in v):.0e}")
        pts[(Lx, Ly)].append((d, G, e))
    for (Lx, Ly), p in sorted(pts.items()):
        p = np.array(sorted(p))
        sel = p[:, 0] <= 0.1001
        if sel.sum() < 4:
            continue
        d, G, e = p[sel].T
        e = np.where(np.isfinite(e), e, 1.0)
        psi, dpsi, c1 = fit_power(d, G, e)
        ps2, A2, B2, c2 = fit_offset(d, G, e)
        ps3, A3, B3, c3 = fit_offset(d, G, e, psi=1 / 3)
        print(f"\nLx {Lx} Ly {Ly}, delta <= 0.1 ({sel.sum()} points): Gamma ~ delta^-psi: psi = {psi:.3f} +- {dpsi:.3f} (chi2 {c1:.1f}); "
              f"A delta^-psi + B: psi = {ps2:.2f}, B = {B2:.1f} (chi2 {c2:.1f}); psi = 1/3 fixed: A = {A3:.1f}, B = {B3:.1f} (chi2 {c3:.1f})")
    if a.csv:
        with open(a.csv, "w") as f:
            f.write("exp,q,Lx,Ly,s,T,done,I_late,W_late,rc_late\n")
            for (q, Lx, Ly, s, exp), r in sorted(runs.items()):
                late = r["t"] >= r["T"] / 2
                f.write(f"{exp},{q},{Lx},{Ly},{s},{r['T']:.0f},{int(r['done'])},{r['I'][late].mean():.4f},"
                        f"{r['W'][late].mean():.3f},{r['rc'][late].mean():.4f}\n")


if __name__ == "__main__":
    main()
