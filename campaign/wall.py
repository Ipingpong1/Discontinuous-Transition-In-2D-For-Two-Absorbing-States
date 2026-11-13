#!/usr/bin/env python3
"""E8b domain-wall analysis (dnsim --init strip:1:0.5 --bg 0 --profile).

    python3 wall.py RUNS_DIR [--every 1]     (RUNS_DIR/E8b_wall/out/*.txt + .prof)
    python3 wall.py RUNS_DIR --exp E11_lazy_wall   lazy runs (_r<r> in the name): the summary uses
                                             physical time t r, and q/q* is relative to q*(r)
    python3 wall.py --run PREFIX             one run in detail

The lattice starts all alive in rows [0, Ly/2) and all dead above, so there are two flat walls,
at y = Ly/2 and at y = 0 (across the seam).  Each wall is analysed in its own half of the torus
(rows within Ly/4 of it).  From the per-row interface counts iface_y (cells with a disagreeing
neighbour), per wall and output time:
  I    = sum_y iface_y / Lx, interface cells per unit wall length;
  sig  = rms distance of the interface cells from their mean row, the layer half-width;
  rows = rows whose interface density exceeds RHO_HALF (half the bulk active value, ~0.3);
  rin  = I / rows, the interface density inside that layer (0 when rows = 0);
  W    = distance between the outermost rows whose 9-row running mean exceeds RHO_HALF, the
         width of the active slab, and rc = mean interface density over those W rows, holes
         included (first order: rc stays at the active-phase value; GV: rc ~ 1/ln t).
Averaged over the two walls and over seeds.  Local exponents over a factor 4 in t:
a_sig = dln sig / dln t (ballistic 1, diffusive 1/2, rough wall ~1/3 or less),
a_I = dln I / dln t.  First order: above q* the layer is active (rin ~ 0.6) and grows
ballistically, below q* it stays a thin wall.  GV at q_c: sig ~ t^1/2, rin ~ 1/ln t.
"""
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

RHO_HALF = 0.3
QSTAR = {0.9: 0.0119, 0.98: 0.00252, 0.99: 0.00127, 0.995: 0.00062}
QSTAR_LAZY = {(0.98, 0.5): 0.0021773}   # make_tasks.QSTAR_LAZY (E10 fronts)


def read(prefix):
    Lx = Ly = None
    for line in open(prefix + ".txt"):
        if line.startswith("# Lx"):
            w = line.split()
            Lx, Ly = int(w[2]), int(w[4])
            break
    raw = np.fromfile(prefix + ".prof", dtype=np.uint8)
    rec = 8 + 8 * Ly
    n = raw.size // rec
    raw = raw[: n * rec].reshape(n, rec)
    t = raw[:, :8].copy().view(np.int64).ravel()
    alive = raw[:, 8: 8 + 4 * Ly].copy().view(np.uint32).astype(float)
    iface = raw[:, 8 + 4 * Ly:].copy().view(np.uint32).astype(float)
    return Lx, Ly, t, alive, iface


def walls(Lx, Ly, iface):
    """Per record: arrays (I, sig, rows) for the walls at y = Ly/2 and y = 0."""
    d = iface / Lx
    out = []
    for c in (Ly // 2, 0):
        y = (np.arange(Ly) - c + Ly // 2) % Ly - Ly // 2          # signed distance from the wall row
        sel = np.flatnonzero(np.abs(y) < Ly // 4)
        sel = sel[np.argsort(y[sel])]                              # rows in order of y
        dd, yy = d[:, sel], y[sel].astype(float)
        I = dd.sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            m = (dd * yy).sum(1) / I
            sig = np.sqrt(np.maximum((dd * yy ** 2).sum(1) / I - m ** 2, 0))
        rows = (dd > RHO_HALF).sum(1)
        W, rc = np.zeros(len(dd)), np.zeros(len(dd))
        sm = np.apply_along_axis(lambda r: np.convolve(r, np.ones(9) / 9, "same"), 1, dd)
        for i in range(len(dd)):
            k = np.flatnonzero(sm[i] > RHO_HALF)
            if k.size:
                W[i] = k[-1] - k[0] + 1
                rc[i] = dd[i, k[0]: k[-1] + 1].mean()
        out.append((I, sig, rows, W, rc))
    return out


def slope(t, y, w=4.0):
    ok = (t > 0) & (y > 0) & np.isfinite(y)
    lt, ly = np.log(t[ok]), np.log(y[ok])
    res = np.full(len(t), np.nan)
    h = np.log(w) / 2
    for i, ti in enumerate(t):
        if ti > 0 and len(lt) > 1 and lt[0] <= np.log(ti) - h and np.log(ti) + h <= lt[-1]:
            res[i] = (np.interp(np.log(ti) + h, lt, ly) - np.interp(np.log(ti) - h, lt, ly)) / (2 * h)
    return res


def run_series(prefix):
    Lx, Ly, t, alive, iface = read(prefix)
    ws = walls(Lx, Ly, iface)
    I = np.mean([w[0] for w in ws], axis=0)
    sig = np.mean([w[1] for w in ws], axis=0)
    rows = np.mean([w[2] for w in ws], axis=0)
    W = np.mean([w[3] for w in ws], axis=0)
    rc = np.mean([w[4] for w in ws], axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        rin = np.where(rows > 0, I / np.maximum(rows, 1), 0.0)
    return t, I, sig, rows, rin, alive.sum(1) / (Lx * Ly), W, rc


def key(prefix):
    n = os.path.basename(prefix)
    mt = re.match(r"p([\d.]+)_a[\d.]+(?:_r([\d.]+))?_j2_q([\d.]+)_Lx\d+_Ly\d+_s(\d+)$", n)
    return float(mt[1]), float(mt[3]), int(mt[4]), float(mt[2] or 1)


def detail(prefix, every=1):
    t, I, sig, rows, rin, al, W, rc = run_series(prefix)
    a_s, a_I, a_W = slope(t, sig), slope(t, I), slope(t, W)
    print(f"== {os.path.basename(prefix)}")
    print(f"{'t':>8} {'I':>8} {'sig':>8} {'rows':>7} {'rin':>6} {'W':>7} {'rc':>6} {'a_sig':>6} {'a_I':>6} "
          f"{'a_W':>6} {'alive':>7}")
    for i in range(0, len(t), every):
        print(f"{t[i]:8d} {I[i]:8.2f} {sig[i]:8.2f} {rows[i]:7.1f} {rin[i]:6.3f} {W[i]:7.1f} {rc[i]:6.3f} "
              f"{a_s[i]:6.2f} {a_I[i]:6.2f} {a_W[i]:6.2f} {al[i]:7.4f}")


def summary(root, exp="E8b_wall"):
    g = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(root, exp, "out", "*.prof"))):
        prefix = f[:-5]
        p, q, s, rr = key(prefix)
        rs = run_series(prefix)
        g[(p, q, rr)].append((rs[0] * rr,) + rs[1:])
    marks = [1e3, 1e4, 1e5, 3e5, 1e6]
    print("per (p, q/q*): seed means at t = " + ", ".join(f"{m:.0e}" for m in marks) + "  (W | rc), then over "
          "t 1e5-1e6 the exponent a_W and the front speed v = (dW/dt)/2 (fit while W < 0.8 of the window)")
    for (p, q, rr), rs in sorted(g.items()):
        def at(idx, m):
            vals = [r[idx][np.argmin(np.abs(r[0] - m))] for r in rs if r[0][-1] >= m]
            return np.mean(vals) if vals else np.nan
        cells = [f"{at(6, m):6.0f} {at(7, m):5.2f}" for m in marks]
        aW, v = [], []
        for r in rs:
            t, W = r[0], r[6]
            sel = (t >= 1e5) & (t <= 1e6) & (W > 0) & (W < 0.8 * 2048)
            if sel.sum() >= 3:
                aW.append(np.polyfit(np.log(t[sel]), np.log(W[sel]), 1)[0])
                v.append(np.polyfit(t[sel], W[sel], 1)[0] / 2)
        tend = min(r[0][-1] for r in rs)
        qs = QSTAR.get(p) if rr == 1 else QSTAR_LAZY.get((p, rr))
        x = f"x {q / qs:5.3f}" if qs else f"q {q:<7g}"
        print(f"p {p:<5}{f' r {rr:g}' if rr != 1 else ''} {x} n {len(rs)} | " + " | ".join(cells) +
              f" | a_W {np.mean(aW) if aW else np.nan:5.2f} v {np.mean(v) if v else np.nan:9.2e}"
              + (f" ±{np.std(v) / np.sqrt(len(v) - 1):.1e}" if len(v) > 1 else "")
              + ("" if tend >= 1e6 else f"  (t_end {tend:.0e})"))


def main():
    a = sys.argv[1:]
    every = int(a[a.index("--every") + 1]) if "--every" in a else 1
    if a[0] == "--run":
        detail(a[1], every)
    else:
        summary(a[0], a[a.index("--exp") + 1] if "--exp" in a else "E8b_wall")


if __name__ == "__main__":
    main()
