#!/usr/bin/env python3
"""Figure for NOTES.md Sec. 8 (fine.png): (a) v(p1) at 1e-4 steps on the two cuts (seed-pooled
adaptive-threshold velocity; open symbols where the strip interior decayed below 80 % of its
t = 1000 value by 1e5 steps; triangles Lx 4096); (b) interior interface density of the strips
(fraction of the t = 1000 value, seed-averaged) against t at p2 = 0.45; (c) stationary rho_bond
of the bulk active phase against dp1 = p1 - p1* with the log and power fits."""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mfront import track, fit_v, KEY  # noqa: E402
from mfine import collect, fit_forms, header, P1STAR, DP_MIN  # noqa: E402

PSTAR = {0.45: 0.18535, 0.50: 0.17057}
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

# (a) velocities, (b) interior histories
hist = {}
for p2, mk, col in ((0.45, "o", "C0"), (0.50, "s", "C3")):
    for sub, mk2 in (("scan_fine", mk), ("scan_fine/wide", "^")):
        pts = {}
        for f in sorted(glob.glob(f"{sub}/p2_{p2:.2f}_p1_*_s*.txt")):
            m = KEY.search(os.path.basename(f))
            if not m or "t_end" not in header(f):
                continue
            r = track(f)
            v, va, vA, n, tm = fit_v(r)
            pts.setdefault(float(m["p1"]), []).append((va, r["cen"][-1], r["t"], r["cen"]))
        if not pts:
            continue
        p1 = np.array(sorted(pts))
        vm = np.array([np.nanmean([x[0] for x in pts[p]]) for p in p1])
        ve = np.array([np.nanstd([x[0] for x in pts[p]], ddof=1) / np.sqrt(max(1, sum(np.isfinite(x[0]) for x in pts[p])))
                       if sum(np.isfinite(x[0]) for x in pts[p]) > 1 else 0 for p in p1])
        cen = np.array([np.mean([x[1] for x in pts[p]]) for p in p1])
        dp = (p1 - PSTAR[p2]) * 1e4
        ok = cen > 0.8
        axes[0].errorbar(dp[ok], vm[ok], ve[ok], fmt=mk2, color=col, ms=6, label=f"$p_2 = {p2}$" + (" ($L_x$ 4096)" if sub.endswith("wide") else ""))
        axes[0].errorbar(dp[~ok], vm[~ok], ve[~ok], fmt=mk2, color=col, ms=6, mfc="none")
        if sub == "scan_fine" and p2 == 0.45:
            for p in p1:
                t = pts[p][0][2]
                c = np.mean([x[3] for x in pts[p]], axis=0)
                axes[1].plot(t, c, label=f"{(p - PSTAR[p2]) * 1e4:+.1f}")
axes[0].axhline(0, color="k", lw=0.5); axes[0].axvline(0, color="k", lw=0.5)
axes[0].set_xlabel(r"$(p_1 - p_1^*)\times 10^4$"); axes[0].set_ylabel("v (rows/step per front)"); axes[0].legend(fontsize=8, loc="lower right")
axes[1].set_xscale("log"); axes[1].set_xlim(1e3, 1.2e5); axes[1].set_ylim(0, 1.05)
axes[1].set_xlabel("t"); axes[1].set_ylabel("strip interior interface density / value at t = 1000")
axes[1].legend(title=r"$(p_1 - p_1^*)\times10^4$", fontsize=8, ncol=2); axes[1].set_title("$p_2 = 0.45$, $L_x$ 2048 x $L_y$ 4096", fontsize=9)

# (c) density against dp1
pts = collect(["bulk_fine/p2_0.45_p1_*_L2048_s*.txt", "bulk_fine/p2_0.45_p1_*_L1024_T1e6_s*.txt",
               "bulk_long/p2_0.45_p1_*_L1024.txt", "bulk/p2_0.45_p1_*_s1.txt"])
best = {}
for (p1, L), rs in pts.items():
    T = max(r[2] for r in rs)
    cand = (L == 2048, T)
    if p1 not in best or cand > best[p1][0]:
        best[p1] = (cand, np.mean([r[0] for r in rs if r[2] == T]))
p1s = np.array(sorted(p for p in best if p >= P1STAR + DP_MIN - 1e-9)); rho = np.array([best[p][1] for p in p1s]); dp = p1s - P1STAR
exc = np.array(sorted(p for p in best if P1STAR < p < P1STAR + DP_MIN - 1e-9))
axes[2].plot(dp, rho, "ko", label="bulk, stationary (L 2048, 3e5 steps)")
if len(exc):
    axes[2].plot(exc - P1STAR, [best[p][1] for p in exc], "ko", mfc="none", label="still decaying at 3e5 (L 2048)")
sel = dp <= 1.0e-3 + 1e-9
if sel.sum() >= 4:
    r = fit_forms(dp[sel], rho[sel])
    x = np.logspace(-5, -2, 200)
    axes[2].plot(x, r["A"] / np.log(r["B"] / x), "C0-", label=f"A/ln(B/dp), p1* = 0.18535: A = {r['A']:.2f}, B = {r['B']:.1f}")
    axes[2].plot(x, r["rho0"] + r["C"] * x ** r["beta"], "C3--", label=f"rho0 + C dp^b, p1* = 0.18535: rho0 = {r['rho0']:.2f}, b = {r['beta']:.2f}")
    # p1* free: best of each form over a grid
    bl = bp = None
    for ps in np.arange(0.18515, 0.185399, 0.0000025):
        d2 = p1s[sel] - ps
        if (d2 <= 0).any():
            continue
        rr = fit_forms(d2, rho[sel])
        if bl is None or rr["chi_log"] < bl[0]:
            bl = (rr["chi_log"], ps, rr)
        if bp is None or rr["chi_pow"] < bp[0]:
            bp = (rr["chi_pow"], ps, rr)
    for b, st, lab in ((bl, "C0:", "log"), (bp, "C3:", "power")):
        chi, ps, rr = b
        xx = np.logspace(-5, -2, 200)
        d2 = xx + (P1STAR - ps)  # plotted against dp = p1 - 0.18535
        ok = d2 > 0
        y = rr["A"] / np.log(rr["B"] / d2[ok]) if lab == "log" else rr["rho0"] + rr["C"] * d2[ok] ** rr["beta"]
        axes[2].plot(xx[ok], y, st, label=f"{lab}, p1* free = {ps:.5f}" + (f", rho0 = {rr['rho0']:.2f}, b = {rr['beta']:.2f}" if lab == "power" else ""))
axes[2].set_xscale("log"); axes[2].set_xlim(2e-5, 1e-2); axes[2].set_ylim(0.1, 0.33)
axes[2].set_xlabel(r"$\Delta p_1 = p_1 - 0.18535$"); axes[2].set_ylabel(r"$\rho_{\rm bond}$ of the active phase"); axes[2].legend(fontsize=7)
axes[2].set_title("$p_2 = 0.45$ (fits over $\\Delta p_1 \\leq 10^{-3}$)", fontsize=9)
for ax, lab in zip(axes, "abc"):
    ax.text(0.02, 0.95, f"({lab})", transform=ax.transAxes, va="top")
fig.tight_layout(); fig.savefig("fine.png", dpi=130)
print("fine.png")
