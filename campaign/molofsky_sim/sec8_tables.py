#!/usr/bin/env python3
"""Tables for NOTES.md Sec. 8 in the column layout of Secs. 3 and 4 (python3 sec8_tables.py > sec8_tables.md)."""
import glob
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mfront import track, fit_v, KEY  # noqa: E402
from mfine import header, at, pooled, CEN_T  # noqa: E402

PSTAR = {0.45: 0.18535, 0.50: 0.17057}


def strip_rows(pattern):
    g = defaultdict(list)
    for f in sorted(glob.glob(pattern)):
        m = KEY.search(os.path.basename(f))
        if not m or "t_end" not in header(f):
            continue
        r = track(f)
        v, va, vA, n, tm = fit_v(r)
        v2, va2, vA2, n2, tm2 = fit_v(r, tmin=30000)
        g[(float(m["p2"]), float(m["p1"]))].append(dict(va=va, va2=va2, vA=vA, Lx=r["Lx"],
                                                      cen=[at(r["t"], r["cen"], tt) for tt in CEN_T]))
    return g


def strip_table(p2, g, gw):
    print(f"**p_2 = {p2:.2f}** (Lx 2048 x Ly 4096, T 1e5, 4 seeds; v_adapt in rows/step per front over records with t >= 2e3 and, "
          "separately, t >= 3e4, and v_int from the integrated width; + = active invades; interior = central interface density "
          "as a fraction of its t = 1000 value at t = 1e4 / 3e4 / 6e4 / 1e5; stationary = no decline beyond 0.03 between 3e4 and 1e5, "
          "gone = below 0.6 at 1e5, decaying otherwise)\n")
    print("| p_1 | p_1 - p_1* | v (t >= 2e3) | v (t >= 3e4) | v_int | seeds | strip interior at 1e4 / 3e4 / 6e4 / 1e5 (fraction of its t = 1000 interface density) |")
    print("|---|---|---|---|---|---|---|")
    for key in sorted(k for k in g if k[0] == p2):
        rs = g[key]
        s1, _, _ = pooled([x["va"] for x in rs]); s2, _, _ = pooled([x["va2"] for x in rs]); s3, _, _ = pooled([x["vA"] for x in rs])
        cen = np.nanmean([x["cen"] for x in rs], axis=0)
        st = "gone" if cen[-1] <= 0.6 else "stationary" if cen[1] - cen[-1] <= 0.03 else "decaying"
        print(f"| {key[1]:.5f} | {(key[1] - PSTAR[p2]) * 1e4:+.1f}e-4 | {s1} | {s2} | {s3} | {len(rs)} | {st} ({' / '.join(f'{c:.2f}' for c in cen)}) |")
    rows = []
    for key in sorted(k for k in gw if k[0] == p2):
        rs = gw[key]
        s1, _, _ = pooled([x["va"] for x in rs]); s2, _, _ = pooled([x["va2"] for x in rs]); s3, _, _ = pooled([x["vA"] for x in rs])
        cen = np.nanmean([x["cen"] for x in rs], axis=0)
        rows.append(f"p_1 = {key[1]:.5f}: v = {s1} / {s2} / {s3} ({len(rs)} seeds, interior {' / '.join(f'{c:.2f}' for c in cen)})")
    print("\nLx 4096 x Ly 4096, T 1e5: " + "; ".join(rows) + "\n")


if __name__ == "__main__":
    g = strip_rows("scan_fine/p2_*_s*.txt")
    gw = strip_rows("scan_fine/wide/p2_*_s*.txt")
    strip_table(0.45, g, gw)
    strip_table(0.50, g, gw)
