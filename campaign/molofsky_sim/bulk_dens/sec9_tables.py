#!/usr/bin/env python3
"""Condensed seed-averaged density table for NOTES.md Sec. 9.3 (one row per p_1 and L): the mean
over the last half of each run (the protocol) with its error, the seed-to-seed sd, Q2 - Q1 and the
drift over that half, the integrated autocorrelation time, the means in the fixed windows, the
status and the sorted seed means.  Run from campaign/molofsky_sim: python3 bulk_dens/sec9_tables.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np  # noqa: E402
import mdens  # noqa: E402

pts = mdens.points(mdens.load_runs(window=None), require_done=True)
W = mdens.WINDOWS
print("| p_1 | dp_1 (1e-4) | L | seeds | T | rho (last half) | err | seed sd | Q2 - Q1 | drift /1e5 | tau | 1.5e5-3e5 | 3e5-6e5 | 6e5-1.2e6 | status | seed means (last half, sorted) |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for p in pts:
    if "rho" not in p:
        continue
    wins = " | ".join(("%.4f (%d)" % (p["win"][w], p["win_n"][w])) if not np.isnan(p["win"][w]) else "-" for w in W)
    seeds = " ".join("%.3f" % m for m, s in p["seed_means"])
    print(f"| {p['p1']:.5f} | {p['dp1'] * 1e4:.1f} | {p['L']} | {p['n']} | {p['T']:.3g} | {p['rho']:.4f} | {p['err']:.4f} | {p['sd']:.4f} | "
          f"{p['dq']:+.4f} +- {p['dq_se']:.4f} | {p['drift']:+.4f} +- {p['drift_se']:.4f} | {p['tau']:.0f} | {wins} | {p['status']} | {seeds} |")
