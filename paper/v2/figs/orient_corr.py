"""Orientation test with the interior drift removed.  A band's interface count is
I(t) = iota(t) A(t) + (front terms), and the interior interface fraction iota(t) of a band started from
a random half relaxes like a bulk system at the same q.  tilt.vel takes iota fixed at its t = 1000 value,
so the slow relaxation of iota enters dI/dt as A d(iota)/dt, a velocity bias proportional to
1/sqrt(a^2 + b^2) (same area, longer fronts).  Here A(t) = I(t) / iota_bulk(t) with iota_bulk(t) the
seed-mean interface fraction of bulk runs (L 2048) at the same (p, q), and v = dA/dt / (2 L sqrt(a^2+b^2))."""
import glob
import math
import os
import re
from collections import defaultdict

import numpy as np

from common import RUNS
import tilt

ORIENT_BULK = "/tmp/claude-1000/-home-matias-Desktop-silly-dn/7f0d9a3b-c595-4ffd-a22c-cdbe8fd1d247/scratchpad/p8/orientbulk/out"
QREF = {0.1: 0.0368, 0.5: 0.0459, 0.9: 0.0119}


def bulk_iota():
    """{(p, q): (t, iota)} from E13_bulk05 (p 0.5, L 2048) and the local orientation references."""
    runs = defaultdict(list)
    for f in glob.glob(os.path.join(RUNS, "E13_bulk05", "out", "p0.5_a0.5_j2_q*_L2048_s*.txt")):
        q = float(re.search(r"_q([\d.]+)_", os.path.basename(f))[1])
        runs[(0.5, q)].append(f)
    for f in glob.glob(os.path.join(ORIENT_BULK, "*.txt")):
        m = re.match(r"p([\d.]+)_q([\d.]+)_L2048_s\d", os.path.basename(f))
        if any(ln.startswith("# done") for ln in open(f)):
            runs[(float(m[1]), float(m[2]))].append(f)
    out = {}
    for k, fs in runs.items():
        ts = [tilt.load_ts(f) for f in fs]
        t = ts[0][:, 0]
        n = min(len(x) for x in ts)
        out[k] = (t[:n], np.mean([x[:n, 4] for x in ts], axis=0))
    return out


def vel_corr(path, iota, t0=2000):
    m = tilt.NAME.search(os.path.basename(path))
    a, b = int(m[1]), int(m[2])
    L = int(re.search(r"_L(\d+)_", os.path.basename(path))[1])
    ts = tilt.load_ts(path)
    t, I = ts[:, 0], ts[:, 4] * L * L
    ib = np.interp(t, *iota)
    A = I / ib / (L * L)
    sel = (t >= t0) & (A > 0.15) & (A < 0.85)
    bad = np.flatnonzero((t >= 1000) & ~((A > 0.15) & (A < 0.85)))
    if bad.size:
        sel &= t < t[bad[0]]
    if sel.sum() < 4:
        return np.nan
    return np.polyfit(t[sel], A[sel] * L * L, 1)[0] / (2 * L * math.hypot(a, b))


def collect(t0=2000):
    io = bulk_iota()
    g = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(RUNS, "E15_orient", "out", "*.txt"))):
        m = re.match(r"p([\d.]+)_a0.5_q([\d.]+)_ab(\d-\d)_L\d+_s(\d)", os.path.basename(f))
        p, q = float(m[1]), float(m[2])
        if (p, q) not in io:
            continue
        g[(p, m[3], q)].append(vel_corr(f, io[(p, q)], t0))
    return g


if __name__ == "__main__":
    from fig_orient import crossing
    for t0 in (2000, 10000):
        g = collect(t0)
        print(f"--- corrected, fit from t >= {t0}")
        for p in sorted(QREF):
            row = []
            for ab in ("1-0", "2-1", "1-1"):
                ks = sorted(k for k in g if k[0] == p and k[1] == ab)
                if len(ks) < 3:
                    continue
                x = np.array([k[2] / QREF[p] for k in ks])
                vs = [np.array(g[k]) for k in ks]
                vm = np.array([np.nanmean(v) for v in vs])
                se = np.array([np.nanstd(v, ddof=1) / np.sqrt(np.isfinite(v).sum()) for v in vs])
                xs, dxs, bb = crossing(x, vm, se)
                row.append(f"{ab}: {xs:.4f}({dxs * 1e4:.0f}) [{len(ks)} q]")
            print(f"p={p}: " + "  ".join(row))
