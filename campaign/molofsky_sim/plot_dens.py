#!/usr/bin/env python3
"""Figure for NOTES.md Sec. 9 (dens.png): (a) seed-averaged rho_bond of the active phase vs
dp_1 = p_1 - p_1* (p_1* = 0.18535) on the cut p_2 = 0.45, in fixed time windows, with the fitted
forms; (b) residuals of the fits with p_1* fixed; (c) the run-by-run history of every run within
0.5e-4 of the reversal (L 2048 thin, L 4096 thick), where the seeds split; (d) the seed-averaged
history at every p_1.  The Sec. 4 single-run points (grey) and the Sec. 8.4 single L 2048 runs
(open) are shown for comparison.  Everything is recomputed from the raw runs in bulk_dens/
through mdens.py.  Usage: python3 plot_dens.py [out.png]
"""
import glob
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mdens  # noqa: E402

out = sys.argv[1] if len(sys.argv) > 1 else "dens.png"
P1STAR = mdens.P1STAR
W1, W2, W3 = mdens.WINDOWS           # (1.5e5, 3e5), (3e5, 6e5), (6e5, 1.2e6)


def window_points(w, L, require_done=False):
    pts = mdens.points(mdens.load_runs(window=w), require_done=require_done)
    return [p for p in pts if p["L"] == L and "rho" in p]


plh = window_points(None, 2048, require_done=True)   # last half of each finished run: the protocol of Sec. 9 (= bulk_dens/points.txt)
pw2 = window_points(W2, 2048)
pw3 = window_points(W3, 2048)
pw1_big = window_points(W1, 4096)
pw2_big = window_points(W2, 4096)
pw3_big = window_points(W3, 4096)
old = mdens.old_points()
sec8 = []                               # Sec. 8.4 single L 2048 runs to 3e5 (mean over the last third)
for f in sorted(glob.glob("bulk_fine/p2_0.45_p1_0.185[4-6]*_L2048_s1.txt")):
    t, r = mdens.series(f)
    sec8.append((float(f.split("_p1_")[1].split("_")[0]), r[t >= 2e5].mean()))

fig, ax = plt.subplots(2, 2, figsize=(13, 9.6))
ax = ax.ravel()

# ---- (a) rho(dp1) and fits (fits: stationary L 2048 points, last half of each run, + Sec. 4 up to 1e-3)
stat = [(p["p1"], p["rho"], p["err"]) for p in plh if p["stationary"]]
pts = sorted(stat + [(p[0], p[1], p[2]) for p in old if p[0] - P1STAR <= 1e-3 + 1e-9])
p1 = np.array([p[0] for p in pts]); rho = np.array([p[1] for p in pts]); err = np.array([p[2] for p in pts])
dp = p1 - P1STAR; w = 1 / err ** 2
fits = [mdens.fit_log(dp, rho, w), mdens.fit_pow(dp, rho, w, rho0=0), mdens.fit_pow(dp, rho, w), mdens.fit_pow(dp, rho, w, beta=0.3)]
labels = ["A/ln(B/dp)", "C dp^beta", "rho0 + C dp^beta (rho0 >= 0)", "rho0 + C dp^0.3"]
styles = ["-", "--", "-.", ":"]
xx = np.logspace(np.log10(1e-5), np.log10(5e-3), 400)
for F, lab, st in zip(fits, labels, styles):
    ax[0].plot(xx, F["f"](xx), st, lw=1.3, label=f"{lab}: {F['label']}, rms {mdens.rms(dp, rho, F['f']):.4f}")
# the same log form with p_1* free (best value), to show where the continuous forms put the end of the branch
best = None
for ps in np.arange(0.18500, p1.min() - 1e-6, 2e-6):
    F = mdens.fit_log(p1 - ps, rho, w)
    if best is None or F["chi"] < best[1]["chi"]:
        best = (ps, F)
ps, F = best
xs = np.logspace(np.log10(1.1e-6), np.log10(5e-3), 400)
ax[0].plot(xs + ps - P1STAR, F["f"](xs), "-", color="0.3", lw=0.8, alpha=0.8,
           label=f"A/ln(B/(p1 - p1*)) with p1* free: p1* = {ps:.6f}, chi2 {F['chi']:.1f}")
if sec8:
    ax[0].plot([p - P1STAR for p, r in sec8], [r for p, r in sec8], "o", mfc="none", mec="0.4", ms=6, label="Sec. 8.4 single L 2048 runs (2e5-3e5)")
ax[0].errorbar([p[0] - P1STAR for p in old], [p[1] for p in old], [p[2] for p in old], fmt="s", color="0.5", ms=5, label="Sec. 4 single runs (L 2048, 6e4), +-0.003")
for sub, mfc, lab in ((True, "C3", "L 2048, 3e5-6e5, stationary"), (False, "none", "L 2048, 3e5-6e5, not stationary (drifting or seeds split)")):
    q = [p for p in pw2 if p["stationary"] == sub]
    if q:
        ax[0].errorbar([p["dp1"] for p in q], [p["rho"] for p in q], [p["err"] for p in q], fmt="o", color="C3", mfc=mfc, ms=6, capsize=2, label=lab)
if pw3:
    ax[0].errorbar([p["dp1"] * 1.04 for p in pw3], [p["rho"] for p in pw3], [p["err"] for p in pw3], fmt="D", color="C1", mfc="none", ms=6, capsize=2,
                   label="L 2048, 6e5-1.2e6 (continued runs)")
for q, mk, lab in ((pw1_big, "^", "L 4096, 1.5e5-3e5"), (pw2_big, "v", "L 4096, 3e5-6e5")):
    if q:
        ax[0].errorbar([p["dp1"] * 0.96 for p in q], [p["rho"] for p in q], [p["err"] for p in q], fmt=mk, color="C2", ms=7, capsize=2, label=lab)
ax[0].plot([p[0] - P1STAR for p in stat], [p[1] for p in stat], "o", mfc="none", mec="k", ms=10, mew=0.8,
           label="fitted: L 2048, last half of each run, stationary (Sec. 9 protocol)")
# individual seeds of the 3e5-6e5 window near the edge
for p in pw2:
    if p["dp1"] <= 0.5e-4 + 1e-9:
        ax[0].plot([p["dp1"]] * len(p["seed_means"]), [m for m, s in p["seed_means"]], ".", color="C3", ms=3, alpha=0.5)
ax[0].axvspan(2.5e-5, 5e-5, color="0.92", zorder=0)
ax[0].text(3.5e-5, 0.305, "stability edge\n(Sec. 8)", fontsize=7, ha="center", color="0.4")
ax[0].set_xscale("log"); ax[0].set_xlim(1.5e-5, 6e-3); ax[0].set_ylim(0.12, 0.32)
ax[0].set_xlabel(r"$\delta p_1 = p_1 - 0.18535$"); ax[0].set_ylabel(r"$\rho_{\rm bond}$ of the active phase")
ax[0].set_title(r"(a) density of the active phase above the line, $p_2 = 0.45$ (dots: single seeds)", fontsize=10)
ax[0].legend(fontsize=6.3, loc="lower right")

# ---- (b) residuals
mks = ["o", "s", "D", "^"]
for F, lab, mk, c in zip(fits, labels, mks, ["C0", "C1", "C2", "C3"]):
    ax[1].errorbar(dp * (1 + 0.03 * (mks.index(mk) - 1.5)), (rho - F["f"](dp)) / err, 1, fmt=mk, color=c, ms=4, capsize=0, alpha=0.8, label=lab)
ax[1].axhline(0, color="k", lw=0.5)
ax[1].set_xscale("log"); ax[1].set_xlabel(r"$\delta p_1$"); ax[1].set_ylabel("(data - fit) / err")
ax[1].set_title("(b) residuals of the fits with $p_1^*$ fixed at 0.18535 (stationary points, last half of each run, up to $10^{-3}$)", fontsize=10)
ax[1].legend(fontsize=7)


def smooth(t, r, n=10):
    tt, rr = mdens.uniform(t, r, 0, t[-1])
    if len(rr) < n:
        return tt, rr
    k = np.ones(n) / n
    return tt[n - 1:], np.convolve(rr, k, mode="valid")


# ---- (c) every run within 0.5e-4 of the reversal
cols = {0.18538: "C0", 0.18539: "C1", 0.18540: "C3"}
tmax_c = 0
seen = set()
for f in sorted(glob.glob("bulk_dens/p2_0.45_p1_*_s*.txt")):
    h = mdens.header(f)
    p = float(f.split("_p1_")[1].split("_")[0])
    if p - P1STAR > 0.5e-4 + 1e-9:
        continue
    t, r = mdens.series(f)
    tt, rr = smooth(t, r)
    tmax_c = max(tmax_c, t[-1])
    big = h["Lx"] == 4096
    lab = f"{p:.5f}, L {h['Lx']}" if (p, big) not in seen else None
    seen.add((p, big))
    ax[2].plot(tt, rr, lw=(2.2 if big else 0.7), ls=("--" if big else "-"), color=cols[round(p, 5)], alpha=(1 if big else 0.75), label=lab)
    if h.get("absorbed"):
        ax[2].plot([h["t_end"]], [0.12], "x", color=cols[round(p, 5)], ms=7)
ax[2].axvspan(W2[0], W2[1], color="0.92", zorder=0)
if tmax_c > 6e5:
    ax[2].axvspan(W3[0], min(W3[1], tmax_c), color="0.85", zorder=0)
ax[2].set_xlim(0, max(6.2e5, tmax_c * 1.02)); ax[2].set_ylim(0.11, 0.22)
ax[2].set_xlabel("t (steps)"); ax[2].set_ylabel(r"$\rho_{\rm bond}(t)$, running mean over $2\times10^4$ steps")
ax[2].set_title(r"(c) every run within $0.5\times10^{-4}$ of the reversal (thin: L 2048 seeds; thick dashed: L 4096)", fontsize=10)
ax[2].legend(fontsize=7, ncol=2, loc="lower left")

# ---- (d) seed-averaged history at every p_1 (L 2048, runs that reached the same length)
g = defaultdict(list)
for f in sorted(glob.glob("bulk_dens/p2_0.45_p1_*_L2048_s*.txt")):
    t, r = mdens.series(f)
    g[float(f.split("_p1_")[1].split("_")[0])].append((t, r))
tmax_all = 0
for i, (p, runs) in enumerate(sorted(g.items())):
    tmax = max(t[-1] for t, r in runs)
    tmax_all = max(tmax_all, tmax)
    tt = np.arange(2000, tmax + 1, 2000)
    rows = []
    for t, r in runs:
        rr = np.interp(tt, t, r); rr[tt > t[-1]] = np.nan
        rows.append(rr)
    rr = np.nanmean(rows, axis=0)
    cnt = np.sum(~np.isnan(rows), axis=0)
    # average only where every run of the group (to 6e5), or every continued run (beyond), contributes
    n6 = cnt[tt <= 6e5].min()
    nx = cnt[tt > 6e5].max() if tmax > 6e5 else 0
    keep = ((tt <= 6e5) & (cnt == n6)) | ((tt > 6e5) & (cnt == nx))
    rr[~keep] = np.nan
    ax[3].plot(tt, rr, lw=1, color=plt.cm.viridis(i / max(1, len(g) - 1)),
               label=f"{p:.5f} ({n6} seeds to 6e5" + (f", {nx} to {tt[keep].max():.1e})" if nx else ")"))
ax[3].axvspan(W2[0], W2[1], color="0.92", zorder=0)
if tmax_all > 6e5:
    ax[3].axvspan(W3[0], W3[1], color="0.85", zorder=0)
ax[3].set_xscale("log"); ax[3].set_xlim(1e3, max(6.5e5, tmax_all * 1.1)); ax[3].set_ylim(0.14, 0.26)
ax[3].set_xlabel("t (steps)"); ax[3].set_ylabel(r"$\rho_{\rm bond}(t)$, seed-averaged, L 2048")
ax[3].set_title("(d) seed-averaged relaxation at each $p_1$; shaded = averaging windows", fontsize=10)
ax[3].legend(fontsize=6.5, ncol=2)
fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote", out)
