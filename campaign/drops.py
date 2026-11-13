#!/usr/bin/env python3
"""Droplet test (disk of the disordered phase in the all-dead absorbing state).

    python3 drops.py $SCRATCH/dn/runs/E5_drop [--rho-b 0.70]

Active area A(t) = iface cells / rho_b, with rho_b the bulk interface fraction of the
disordered phase (default: interpolated from the alpha=0.5 E1 strip interiors), and
R(t) = sqrt(A/pi).  Per (q, R0): fraction of seeds that died (absorbed), grew past
2 R0, or neither by T; mean R(T); time to die.  First order predicts a critical
radius R_c ~ 1/(q - q*): P(survive) steps from 0 to 1 around R_c.
"""
import argparse
import glob
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fronts import group_key, fmt_key  # noqa: E402

# E1 (alpha=0.5, j=2) strip-interior interface fraction vs q
Q_E1 = [0.016, 0.017, 0.0175, 0.018, 0.0185, 0.019, 0.0195, 0.020, 0.021, 0.023, 0.026]
RB_E1 = [0.592, 0.622, 0.640, 0.648, 0.660, 0.680, 0.687, 0.701, 0.722, 0.758, 0.800]


def run(path, rho_b):
    L, T, done = None, None, None
    rows = []
    for line in open(path):
        if line.startswith("# Lx"):
            w = line.split()
            L, T = int(w[2]) * int(w[4]), int(w[6])
        elif line.startswith("# done"):
            done = line.split()
        elif not line.startswith("#"):
            rows.append([float(x) for x in line.split()])
    if done is None:
        return None
    ts = np.array(rows)
    R = np.sqrt(ts[:, 4] * L / rho_b / np.pi)
    return dict(t=ts[:, 0], R=R, absorbed=int(done[5]), t_end=int(done[3]), T=T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("expdir")
    ap.add_argument("--rho-b", type=float)
    a = ap.parse_args()
    g = defaultdict(list)
    for p in sorted(glob.glob(os.path.join(a.expdir, "out", "*.txt"))):
        key = group_key(os.path.basename(p)[:-4])
        kd = dict(key)
        rb = a.rho_b or float(np.interp(kd["q"], Q_E1, RB_E1))
        r = run(p, rb)
        if r:
            g[key].append(r)
    print("params".ljust(37) + "n   died grew neither  P_surv  R(T)_mean   t_die_median")
    for key in sorted(g, key=lambda k: (dict(k)["q"], dict(k)["R"])):
        rs, R0 = g[key], dict(key)["R"]
        died = sum(r["absorbed"] for r in rs)
        grew = sum((not r["absorbed"]) and r["R"][-1] > 2 * R0 for r in rs)
        surv = [r for r in rs if not r["absorbed"]]
        Rm = np.mean([r["R"][-1] for r in surv]) if surv else 0.0
        td = [r["t_end"] for r in rs if r["absorbed"]]
        tdm = f"{np.median(td):.0f}" if td else "-"
        print(fmt_key(key).ljust(36) + f" {len(rs):<3d} {died:4d} {grew:4d} {len(rs) - died - grew:6d}"
              f"   {1 - died / len(rs):5.2f}  {Rm:9.1f}   {tdm}")


if __name__ == "__main__":
    main()
