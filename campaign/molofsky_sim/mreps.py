#!/usr/bin/env python3
"""Seed spreading (dnsim --reps --init disk:p0:R): <A>(t) = mean live cells per replica,
P_surv(t), <I>(t), and the late-time exponential rate of <A>: r = d ln<A>/dt over the last
decade in t.  At the voter point <A> is a martingale (= initial 4 for the 2x2 block); r < 0
means the absorbing state is linearly stable (the mean population of a small seed decays),
r > 0 means it is unstable.

    python3 mreps.py reps/*.txt
"""
import glob
import os
import sys

import numpy as np


def load(path):
    rows = [l.split() for l in open(path) if not l.startswith("#")]
    a = np.array(rows, float)
    done = any(l.startswith("# done") for l in open(path))
    ts = np.array(sorted(set(a[:, 1])))
    n = np.array([a[a[:, 1] == t, 2].sum() for t in ts])
    surv = np.array([a[a[:, 1] == t, 3].sum() for t in ts])
    A = np.array([a[a[:, 1] == t, 7].sum() for t in ts])
    I = np.array([a[a[:, 1] == t, 5].sum() for t in ts])
    # jackknife over groups for <A>(T)
    G = sorted(set(a[:, 0]))
    tT = ts[-1]
    rows = a[a[:, 1] == tT]
    Ag, ng = rows[:, 7], rows[:, 2]
    jk = [(Ag.sum() - Ag[i]) / (ng.sum() - ng[i]) for i in range(len(G))]
    errA = np.sqrt((len(G) - 1) / len(G) * ((np.array(jk) - np.mean(jk)) ** 2).sum())
    return dict(t=ts, n=n, surv=surv, A=A / n, I=I / n, P=surv / n, done=done, errA=errA)


def main(files):
    print("run                         n      <A>(100) <A>(1e3)  <A>(T)   err    P_s(1e3) P_s(T)   r(last decade)  done")
    for f in files:
        r = load(f)
        t, A = r["t"], r["A"]
        T = t[-1]
        sel = (t >= T / 10) & (A > 0)
        rate = np.polyfit(t[sel], np.log(A[sel]), 1)[0] if sel.sum() >= 3 else np.nan
        g = lambda tt, arr: float(np.interp(tt, t, arr))
        print(f"{os.path.basename(f)[:-4]:26s} {int(r['n'][0]):6d}  {g(100, A):8.3f} {g(1000, A):8.3f} {A[-1]:8.3f} "
              f"{r['errA']:6.3f}  {g(1000, r['P']):.4f}  {r['P'][-1]:.4f}   {rate:+.2e}       {int(r['done'])}")


if __name__ == "__main__":
    files = [f for pat in sys.argv[1:] for f in sorted(glob.glob(pat))] or sys.argv[1:]
    main(files)
