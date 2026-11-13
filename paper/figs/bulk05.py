"""E13_bulk05 reduction: per run the plateau, the half-life of the active phase and chi_L; per (x, L)
the pooled values.  Used by fig_bulk05.py and printed as a table.

    python3 bulk05.py            table of every (x, L)
"""
import glob
import os
import re
from collections import defaultdict

import numpy as np

from common import RUNS

QS = 0.04597         # strip-front threshold (E12_phase, 4 seeds); run names use 0.0459 x grid
RHO_REF = 0.29          # bulk interface density just above q* (E4_voter_bulk: 0.290 at 1.002 q*)


def load(path):
    T = None
    done = absorbed = None
    rows = []
    for ln in open(path):
        if ln.startswith("# Lx"):
            T = int(ln.split()[6])
        elif ln.startswith("# done"):
            w = ln.split()
            done, absorbed = int(w[3]), int(w[5])
        elif ln.strip() and not ln.startswith("#"):
            rows.append([float(v) for v in ln.split()])
    return np.array(rows), T, done, absorbed


def per_run(path):
    ts, T, t_end, absorbed = load(path)
    if t_end is None:
        return None
    t, n, rho = ts[:, 0], ts[:, 1], ts[:, 2]
    m = 2 * n - 1
    below = np.flatnonzero((t > 1000) & (rho < RHO_REF / 2))
    tau = t[below[0]] if below.size else np.inf
    late = (t >= t_end / 2) & (t >= 1000)
    alive_late = late & (rho > 0.8 * RHO_REF)
    L = int(re.search(r"_L(\d+)_", os.path.basename(path))[1])
    plateau = rho[alive_late].mean() if alive_late.sum() >= 10 and np.isinf(tau) else np.nan
    chi = L * L * np.mean(m[alive_late] ** 2) if alive_late.sum() >= 10 and np.isinf(tau) else np.nan
    return dict(t=t, rho=rho, m=m, tau=tau, t_end=t_end, absorbed=absorbed, plateau=plateau, chi=chi, T=T)


def collect():
    g = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(RUNS, "E13_bulk05", "out", "*.txt"))):
        mt = re.match(r"p0.5_a0.5_j2_q([\d.]+)_L(\d+)_s(\d+)", os.path.basename(f))
        r = per_run(f)
        if r is None:
            continue
        r["seed"] = int(mt[3])
        g[(round(float(mt[1]) / QS, 3), int(mt[2]))].append(r)
    return g


def summary(g):
    out = {}
    for k, rs in g.items():
        taus = np.array([r["tau"] for r in rs])
        fin = np.isfinite(taus)
        pl = np.array([r["plateau"] for r in rs])
        chis = np.array([r["chi"] for r in rs])
        out[k] = dict(n=len(rs), n_decayed=int(fin.sum()), tau_med=np.median(taus) if fin.all() else np.nan,
                      tau_mean=taus[fin].mean() if fin.any() else np.nan,
                      tau_sem=taus[fin].std(ddof=1) / np.sqrt(fin.sum()) if fin.sum() > 1 else np.nan,
                      plateau=np.nanmean(pl) if np.isfinite(pl).any() else np.nan,
                      plateau_sem=np.nanstd(pl, ddof=1) / np.sqrt(np.isfinite(pl).sum()) if np.isfinite(pl).sum() > 1 else np.nan,
                      chi=np.nanmean(chis) if np.isfinite(chis).any() else np.nan,
                      absorbed=sum(r["absorbed"] for r in rs), T=rs[0]["T"],
                      t_abs=np.array([r["t_end"] for r in rs if r["absorbed"]]))
    return out


if __name__ == "__main__":
    s = summary(collect())
    print("   x      L   runs decayed  tau_mean(sem)       plateau(sem)     chi    absorbed  T")
    for (x, L) in sorted(s):
        d = s[(x, L)]
        print(f"{x:5.3f} {L:6d} {d['n']:5d} {d['n_decayed']:5d}   {d['tau_mean']:9.0f} ({d['tau_sem']:6.0f})  "
              f"{d['plateau']:.4f} ({d['plateau_sem']:.4f})  {d['chi']:7.1f}  {d['absorbed']:4d}  {d['T']:.0e}")
