#!/usr/bin/env python3
"""v(p1) at p2 = 0.45 and 0.50 (seed-pooled, adaptive-threshold velocity; open symbols where the
strip interior had decayed below 60 % by the end of the run), with the closure's rays
theta_a, theta_front, theta_half at T = 16 and 64 mapped to p1, and the asymptotic line."""
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mfront import track, fit_v, KEY  # noqa: E402

RAYS = {16: (104.39, 107.33, 108.82), 64: (105.40, 106.69, 107.33)}
ASY = 106.12
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, p2 in zip(axes, (0.45, 0.50)):
    pts = {}
    for f in sorted(glob.glob(f"scan/p2_{p2:.2f}_p1_*_s*.txt")):
        m = KEY.search(os.path.basename(f))
        r = track(f)
        v, va, vA, n, tm = fit_v(r)
        pts.setdefault(float(m["p1"]), []).append((va, r["cen"][-1]))
    p1 = np.array(sorted(pts))
    vm = np.array([np.nanmean([x[0] for x in pts[p]]) for p in p1])
    ve = np.array([np.nanstd([x[0] for x in pts[p]], ddof=1) / np.sqrt(sum(np.isfinite(x[0]) for x in pts[p]))
                   if sum(np.isfinite(x[0]) for x in pts[p]) > 1 else 0 for p in p1])
    stable = np.array([np.mean([x[1] for x in pts[p]]) > 0.6 for p in p1])
    ok = np.isfinite(vm)
    ax.errorbar(p1[ok & stable], vm[ok & stable], ve[ok & stable], fmt="o", color="C0", ms=5, label="interior stationary")
    ax.errorbar(p1[ok & ~stable], vm[ok & ~stable], ve[ok & ~stable], fmt="o", mfc="none", color="C0", ms=5,
                label="interior decaying (v biased)")
    d2 = p2 - 0.4
    for T, (a, fr, h), c in zip(RAYS, RAYS.values(), ("C3", "C1")):
        for th, ls, lab in ((a, ":", "theta_a"), (fr, "-", "theta_front"), (h, "--", "theta_half")):
            x = 0.2 + d2 / np.tan(np.radians(th))
            ax.axvline(x, color=c, ls=ls, lw=1, alpha=0.8, label=f"closure T={T} {lab}" if lab == "theta_front" else None)
    ax.axvline(0.2 + d2 / np.tan(np.radians(ASY)), color="k", ls="-", lw=0.8, alpha=0.5, label="closure T->inf")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("p1"); ax.set_ylabel("front velocity v (rows/step, + = active invades)")
    ax.set_title(f"p2 = {p2}: v(p1), {len(p1)} values")
    lo, hi = p1.min() - 0.0005, p1.max() + 0.0005
    ax.set_xlim(lo, hi)
    ax.legend(fontsize=7, loc="upper left")
    ins = ax.inset_axes([0.58, 0.08, 0.4, 0.42])
    sel = (p1 > (0.1835 if p2 == 0.45 else 0.1685)) & (p1 < (0.1875 if p2 == 0.45 else 0.1735)) & ok
    ins.errorbar(p1[sel & stable], vm[sel & stable], ve[sel & stable], fmt="o", color="C0", ms=4)
    ins.errorbar(p1[sel & ~stable], vm[sel & ~stable], ve[sel & ~stable], fmt="o", mfc="none", color="C0", ms=4)
    for T, (a, fr, h), c in zip(RAYS, RAYS.values(), ("C3", "C1")):
        for th, ls in ((a, ":"), (fr, "-"), (h, "--")):
            ins.axvline(0.2 + d2 / np.tan(np.radians(th)), color=c, ls=ls, lw=1, alpha=0.8)
    ins.axvline(0.2 + d2 / np.tan(np.radians(ASY)), color="k", lw=0.8, alpha=0.5)
    ins.axhline(0, color="k", lw=0.5); ins.set_ylim(-0.012, 0.025); ins.tick_params(labelsize=6)
fig.tight_layout()
fig.savefig("fronts.png", dpi=130)
print("fronts.png written")
