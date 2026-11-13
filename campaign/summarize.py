#!/usr/bin/env python3
"""Summarise one experiment's outputs: per run status and late-time averages, then
seed-averaged rows per parameter set.  For strip runs with profiles (E1) also the
front velocity v = (dA/dt)/2 from a linear fit of the active width A(t) (in rows of
bulk-disordered-equivalent, see profile.py) over the second half of the run.

    python3 summarize.py ~/dn/runs/E1_front [--rho-if 0.71]
"""
import argparse
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from profile import read as read_profile  # noqa: E402

KEY = re.compile(r"(?:^|_)(?P<k>p|a|j|q|Lx|Ly|L|R|s)(?P<v>[0-9.]+)")


def params(name):
    return {m["k"]: m["v"] for m in KEY.finditer(name)}


def load_ts(path):
    rows, done, T = [], None, None
    for line in open(path):
        if line.startswith("# done"):
            w = line.split()
            done = dict(t_end=int(w[3]), absorbed=int(w[5]), wall=float(w[7]))
        elif line.startswith("# Lx"):
            T = int(line.split()[6])
        elif not line.startswith("#"):
            rows.append([float(x) for x in line.split()])
    return np.array(rows), done, T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("expdir")
    ap.add_argument("--rho-if", type=float, default=None,
                    help="bulk iface density of the disordered phase for front widths")
    a = ap.parse_args()
    runs = sorted(glob.glob(os.path.join(a.expdir, "out", "*.txt")))
    groups = defaultdict(list)
    n_done = 0
    for path in runs:
        name = os.path.basename(path)[:-4]
        ts, done, T = load_ts(path)
        if done is None or len(ts) == 0:
            continue
        n_done += 1
        late = ts[ts[:, 0] >= ts[-1, 0] / 2]
        rec = dict(rho=late[:, 2].mean(), act=late[:, 3].mean(),
                   absm=np.abs(2 * late[:, 1] - 1).mean(), absorbed=done["absorbed"],
                   t_end=done["t_end"], v=np.nan)
        prof = path[:-4] + ".prof"
        if os.path.exists(prof):
            Lx, Ly, t, alive, iface = read_profile(path, prof)
            dens = iface / Lx
            rho_if = a.rho_if or float(np.median(dens[1][dens[1] > 0]))
            A = dens.sum(1) / rho_if
            sel = t >= t[-1] / 2
            if sel.sum() >= 3 and not done["absorbed"]:
                rec["v"] = np.polyfit(t[sel], A[sel], 1)[0] / 2
            elif done["absorbed"]:
                rec["v"] = -np.inf  # strip died out completely
        p = params(name)
        key = tuple((k, p[k]) for k in sorted(p) if k != "s")
        groups[key].append(rec)
    print(f"# {a.expdir}: {n_done}/{len(runs)} runs done")
    print("# params  nseeds  rho  act  |m|  absorbed  v_front(+-se)")
    for key in sorted(groups, key=lambda k: [(kk, float(vv)) for kk, vv in k]):
        rs = groups[key]
        f = lambda x: np.array([r[x] for r in rs], float)  # noqa: E731
        v = f("v")
        vtxt = ""
        if not np.all(np.isnan(v)):
            fin = v[np.isfinite(v)]
            vtxt = f"{fin.mean():+.5f}" if len(fin) else ""
            if len(fin) > 1:
                vtxt += f"({fin.std(ddof=1) / np.sqrt(len(fin)):.5f})"
            if np.isneginf(v).any():
                vtxt += f" died:{int(np.isneginf(v).sum())}"
        ktxt = " ".join(f"{k}={v_}" for k, v_ in key)
        print(f"{ktxt:42s} {len(rs):2d}  {f('rho').mean():.4f} {f('act').mean():.4f} "
              f"{f('absm').mean():.3f}  {int(f('absorbed').sum())}/{len(rs)}  {vtxt}")


if __name__ == "__main__":
    main()
