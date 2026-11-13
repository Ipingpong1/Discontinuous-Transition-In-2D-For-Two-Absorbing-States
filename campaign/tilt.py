#!/usr/bin/env python3
"""Tilted active bands for the orientation test (E15_orient).

    python3 tilt.py make OUT.ckpt --L 4096 --a 2 --b 1 [--frac 0.5] [--p0 0.5] [--seed 1]
    python3 tilt.py vel RUN.txt [...]          front velocity per run (time series only)
    python3 tilt.py table EXP_DIR              seed-mean v(q) per (p, a, b)

A band of direction (a, b) (gcd 1) on the L x L torus is the set of cells with
u = (a y - b x) mod L < frac L.  Its two edges are closed lines of length L sqrt(a^2 + b^2)
each, its normal width is frac L / sqrt(a^2 + b^2) and its normal period L / sqrt(a^2 + b^2).
Cells inside are alive with probability p0 (numpy RNG keyed by --seed), cells outside dead:
(a, b) = (1, 0) is dnsim's --init strip:p0:frac --bg 0 up to the random draw.

The checkpoint (dnsim format: "DNSIM1\\0\\0", u64 Lx Ly t, then Lx/64 * Ly u64 words, bit b of
word x in row y = cell (64x + b, y)) is started with `dnsim --load`, at t = 0.

Velocity: the number of interface cells I(t) = iface_frac * L^2 changes as the band grows or
shrinks, dI/dt = rho_in * 2 L sqrt(a^2 + b^2) * v, with rho_in the interface fraction of the
active interior, taken as I(t0) / (frac L^2) at t0 = T0 (the band has not moved yet; |v| T0 is
a few rows).  v is the slope of I(t) over [T0, t_end] by least squares, per run; the fit stops
early if the band fills more than 85% or less than 15% of the torus.  The zero crossing of v(q)
does not depend on rho_in.
"""
import argparse
import glob
import math
import os
import re
import sys
from collections import defaultdict

import numpy as np

MAGIC = b"DNSIM1\0\0"
T0 = 1000


def make(path, L, a, b, frac=0.5, p0=0.5, seed=1):
    if math.gcd(a, b) != 1:
        raise SystemExit("a, b must be coprime")
    if L % 64:
        raise SystemExit("L must be a multiple of 64")
    y, x = np.mgrid[0:L, 0:L]
    u = (a * y - b * x) % L
    band = u < int(round(frac * L))
    rng = np.random.default_rng(seed)
    g = (band & (rng.random((L, L)) < p0)).astype(np.uint8)
    words = np.packbits(g.reshape(L, L // 64, 64), axis=2, bitorder="little").view(np.uint64)
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        f.write(MAGIC)
        f.write(np.array([L, L, 0], dtype=np.uint64).tobytes())
        f.write(words.tobytes())
    os.replace(tmp, path)


def load_ts(path):
    return np.array([[float(v) for v in ln.split()] for ln in open(path) if ln.strip() and not ln.startswith("#")])


NAME = re.compile(r"_ab(\d+)-(\d+)_")


def vel(path, frac=0.5):
    """(v, rho_in, t_fit_end, npts) for one run; a, b, L from the name and header."""
    m = NAME.search(os.path.basename(path))
    a, b = int(m[1]), int(m[2])
    L = None
    for ln in open(path):
        if ln.startswith("# Lx"):
            L = int(ln.split()[2])
            break
    ts = load_ts(path)
    t, iface = ts[:, 0], ts[:, 4] * L * L
    ok = t >= T0
    if ok.sum() < 4:
        return np.nan, np.nan, 0, 0
    rho_in = np.interp(T0, t, iface) / (frac * L * L)
    area = iface / rho_in / (L * L)
    sel = ok & (area > 0.15) & (area < 0.85)
    # stop at the first excursion: a band that filled or vanished does not come back
    bad = np.flatnonzero(ok & ~sel)
    if bad.size:
        sel &= t < t[bad[0]]
    if sel.sum() < 4:
        return np.nan, rho_in, 0, int(sel.sum())
    slope = np.polyfit(t[sel], iface[sel], 1)[0]
    v = slope / (rho_in * 2 * L * math.hypot(a, b))
    return v, rho_in, t[sel][-1], int(sel.sum())


KEY = re.compile(r"(?:^|_)(?P<k>p|a|eta|q|ab|L|s)(?P<v>[0-9.\-]+)")


def table(exp):
    g = defaultdict(list)
    for p in sorted(glob.glob(f"{exp}/out/*.txt")):
        n = os.path.basename(p)[:-4]
        kv = {mm["k"]: mm["v"] for mm in KEY.finditer(n)}
        g[(kv.get("p", "-"), kv["ab"], float(kv.get("q", kv.get("eta", "nan"))))].append(vel(p))
    print("p      ab    q          v_mean      sem       seeds                                   rho_in")
    for k in sorted(g, key=lambda k: (k[0], k[1], k[2])):
        vs = np.array([r[0] for r in g[k]])
        fin = vs[np.isfinite(vs)]
        sem = fin.std(ddof=1) / math.sqrt(len(fin)) if len(fin) > 1 else np.nan
        print(f"{k[0]:6s} {k[1]:5s} {k[2]:.7f}  {fin.mean() if len(fin) else np.nan:+.6f}  {sem:.6f}  "
              + " ".join(f"{v:+.5f}" for v in vs).ljust(40) + f" {np.nanmean([r[1] for r in g[k]]):.3f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("make")
    m.add_argument("out")
    m.add_argument("--L", type=int, required=True)
    m.add_argument("--a", type=int, required=True)
    m.add_argument("--b", type=int, required=True)
    m.add_argument("--frac", type=float, default=0.5)
    m.add_argument("--p0", type=float, default=0.5)
    m.add_argument("--seed", type=int, default=1)
    v = sub.add_parser("vel")
    v.add_argument("runs", nargs="+")
    t = sub.add_parser("table")
    t.add_argument("exp")
    a = ap.parse_args()
    if a.cmd == "make":
        make(a.out, a.L, a.a, a.b, a.frac, a.p0, a.seed)
    elif a.cmd == "vel":
        for r in a.runs:
            v_, rho, te, n = vel(r)
            print(f"{os.path.basename(r)}  v {v_:+.6f}  rho_in {rho:.3f}  t_fit_end {te:.0f}  npts {n}")
    else:
        table(a.exp.rstrip("/"))


if __name__ == "__main__":
    sys.exit(main())
