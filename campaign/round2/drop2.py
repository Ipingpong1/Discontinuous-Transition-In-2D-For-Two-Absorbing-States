#!/usr/bin/env python3
"""Two-droplet coexistence test at p = 1/2 (R2_drop2): an active disk of radius R0 in the dead state
(init disk:0.5:R0) and a dead hole of radius R0 in the active state (init hole:0.5:R0), at
q/q* = 1 -+ 0.003, 0.006, L 4096.

    python3 drop2.py EXP_DIR [--t0 10000] [--bulk BULK_DIR]

The active phase has alive density 1/2 (Z2 symmetry), the absorbing phases 0 or 1, so from the time
series alone: disk area A = 2 alive_frac L^2, hole area A = L^2 - 2 alive_frac L^2, R = sqrt(A/pi).
From the row profile, the hole's chord in row y is c(y) = Lx - 2 alive(y) (the disk's 2 alive(y)), and
A_prof = sum over the rows within 1.3 R0 of the centre of max(c, 0): the same with less noise from the
active phase's magnetisation.  dR/dt by least squares over t >= t0 (while 0.3 R0 < R < 1.8 R0).
Growth law with the Table I slope q* dv/dq = 0.22 and D = 0.20 (Fig. 4): disk dR/dt = v - D/R,
hole dR/dt = -(v + D/R), v = 0.22 (x - 1).  Generic coexistence of width > 0.6% would make both
droplets shrink at both signs of x - 1."""
import argparse
import glob
import math
import os
import re
from collections import defaultdict

import numpy as np

QS = 0.04596
NAME = re.compile(r"p0\.5_a0\.5_j2_q([\d.]+)_(disk|hole)R(\d+)_L(\d+)_s(\d+)$")


def load(prefix):
    L = None
    rows = []
    done = False
    for line in open(prefix + ".txt"):
        if line.startswith("# Lx"):
            L = int(line.split()[2])
        elif line.startswith("# done"):
            done = True
        elif not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    ts = np.array(rows)
    prof = None
    if os.path.exists(prefix + ".prof"):
        raw = np.fromfile(prefix + ".prof", dtype=np.uint8)
        rec = 8 + 8 * L
        n = raw.size // rec
        raw = raw[: n * rec].reshape(n, rec)
        tp = raw[:, :8].copy().view(np.int64).ravel()
        alive = raw[:, 8: 8 + 4 * L].copy().view(np.uint32).astype(float)
        prof = (tp, alive)
    return L, ts, prof, done


def radii(L, ts, prof, kind, R0):
    t = ts[:, 0]
    A_ts = 2 * ts[:, 1] * L * L if kind == "disk" else L * L - 2 * ts[:, 1] * L * L
    R_ts = np.sqrt(np.maximum(A_ts, 0) / math.pi)
    R_pr = None
    if prof is not None:
        tp, alive = prof
        y = np.arange(L)
        sel = np.abs(y - L / 2) < 1.3 * R0 + 50
        chord = 2 * alive[:, sel] if kind == "disk" else L - 2 * alive[:, sel]
        A_pr = np.maximum(chord, 0).sum(1)
        R_pr = (tp, np.sqrt(A_pr / math.pi))
    return t, R_ts, R_pr


def fit(t, R, R0, t0):
    sel = (t >= t0) & (R > 0.3 * R0) & (R < 1.8 * R0)
    bad = np.flatnonzero((t >= t0) & ~sel)
    if bad.size:
        sel &= t < t[bad[0]]
    if sel.sum() < 5:
        return np.nan, np.nan, 0
    b, a = np.polyfit(t[sel], R[sel], 1)
    res = R[sel] - (a + b * t[sel])
    err = np.sqrt(res.var(ddof=2) / ((t[sel] - t[sel].mean()) ** 2).sum())
    return b, err, t[sel][-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exp")
    ap.add_argument("--t0", type=float, default=10000)
    ap.add_argument("--slope", type=float, default=0.22)
    ap.add_argument("--D", type=float, default=0.20)
    a = ap.parse_args()
    g = defaultdict(list)
    print("run                                               R(t0)   R(end)  dR/dt(ts)     dR/dt(prof)   t_fit_end  done")
    for f in sorted(glob.glob(os.path.join(a.exp, "out", "*.txt"))):
        prefix = f[:-4]
        m = NAME.search(os.path.basename(prefix))
        if not m:
            continue
        q, kind, R0, L, s = float(m[1]), m[2], int(m[3]), int(m[4]), int(m[5])
        L_, ts, prof, done = load(prefix)
        t, R_ts, R_pr = radii(L, ts, prof, kind, R0)
        b1, e1, te = fit(t, R_ts, R0, a.t0)
        b2, e2, _ = fit(R_pr[0], R_pr[1], R0, a.t0) if R_pr else (np.nan, np.nan, 0)
        i0 = np.argmin(np.abs(t - a.t0))
        print(f"{os.path.basename(prefix):<50s} {R_ts[i0]:6.0f}  {R_ts[-1]:6.0f}  {b1:+.2e}+-{e1:.0e}  "
              f"{b2:+.2e}+-{e2:.0e}  {te:9.0f}  {int(done)}")
        g[(q, kind, R0)].append((b1, b2))
    print("\nq         x        kind  n   dR/dt (time series)      dR/dt (profile)         predicted v -+ D/R")
    for (q, kind, R0), vs in sorted(g.items()):
        x = q / QS
        v = a.slope * (x - 1)
        pred = v - a.D / R0 if kind == "disk" else -(v + a.D / R0)
        b1 = np.array([v_[0] for v_ in vs])
        b2 = np.array([v_[1] for v_ in vs])
        def ms(b):
            b = b[np.isfinite(b)]
            return (f"{b.mean():+.2e} +- {b.std(ddof=1) / math.sqrt(len(b)):.1e}" if len(b) > 1
                    else (f"{b.mean():+.2e}" if len(b) else "  nan"))
        print(f"{q:<9} {x:<8.4f} {kind:<5} {len(vs):<3d} {ms(b1):<24s} {ms(b2):<24s} {pred:+.2e}")


if __name__ == "__main__":
    main()
