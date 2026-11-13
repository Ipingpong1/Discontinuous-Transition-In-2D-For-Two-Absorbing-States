#!/usr/bin/env python3
"""Compact v(p1) tables for NOTES.md from scan/ and scan_wide/: adaptive-threshold velocity
pooled over seeds (mean +- error of the mean), number of seeds, interior status at the end of
the run (stationary: interior interface density > 0.8 of its t = 1000 value; decaying: 0.6-0.8;
gone: < 0.6 or no fit possible)."""
import glob
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mfront import track, fit_v, KEY  # noqa: E402


def table(pattern):
    g = defaultdict(list)
    for f in sorted(glob.glob(pattern)):
        m = KEY.search(os.path.basename(f))
        r = track(f)
        v, va, vA, n, tm = fit_v(r)
        g[(float(m["p2"]), float(m["p1"]))].append((va, r["cen"][-1], r["Lx"]))
    rows = []
    for key in sorted(g):
        va = np.array([x[0] for x in g[key]])
        fin = va[np.isfinite(va)]
        cen = np.mean([x[1] for x in g[key]])
        st = "stationary" if cen > 0.8 else "decaying" if cen > 0.6 else "gone"
        if len(fin) == 0:
            vs = "-"
        elif len(fin) == 1:
            vs = f"{fin[0]:+.4f}"
        else:
            vs = f"{fin.mean():+.4f} +- {fin.std(ddof=1) / np.sqrt(len(fin)):.4f}"
        rows.append((key, vs, len(g[key]), st, cen, g[key][0][2]))
    return rows


if __name__ == "__main__":
    for p2 in ("0.45", "0.50"):
        print(f"\n**p_2 = {p2}** (Lx 1024 x Ly 4096, T 6e4; v_adapt in rows/step per front, + = active invades)\n")
        print("| p_1 | v | seeds | strip interior at the end (fraction of its t = 1000 interface density) |")
        print("|---|---|---|---|")
        for key, vs, n, st, cen, Lx in table(f"scan/p2_{p2}_p1_*_s*.txt"):
            print(f"| {key[1]:.5f} | {vs} | {n} | {st} ({cen:.2f}) |")
        wide = table(f"scan_wide/p2_{p2}_p1_*_s*.txt")
        if wide:
            print("\nLx 2048 x Ly 4096: " + "; ".join(f"p_1 = {k[1]:.5f}: v = {vs} ({n} seeds, {st} {cen:.2f})" for k, vs, n, st, cen, Lx in wide))
