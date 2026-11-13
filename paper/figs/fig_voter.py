"""Fig. voter: near the voter point. (a) chi_L = L^2 <m^2> against L (E6/E7); (b) seed survival
t P_s(t) (E8a, E8a_long); (c) width of the active layer on a wall between the absorbing phases (E8b)."""
import csv
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np

from common import FIGDATA, PAGE, RUNS, label, save
import spread

fig, axs = plt.subplots(1, 3, figsize=(PAGE, 2.25), gridspec_kw=dict(wspace=0.32))

# (a) chi_L: campaign/binder.py on the E6 + E7 runs (t >= 2e5, active window), seed-scatter errors;
# m = 2n - 1 as in the paper (recomputed 2026-10-31: until then binder.py used m = n - 1/2 and this
# table was 4x too small)
CHI = {  # (p, q/q*): [(L, chi, err)]
    (0.9, 1.03): [(512, 1293, 102), (1024, 1413, 171), (2048, 1253, 16), (4096, 1649, 210)],
    (0.98, 1.02): [(512, 10450, 1388), (1024, 14669, 1503), (2048, 10147, 1186), (4096, 11553, 1409)],
    (0.98, 1.2): [(512, 1604, 353), (1024, 1998, 174), (2048, 1790, 135), (4096, 1763, 165)],
    (0.99, 1.01): [(512, 16702, 7572), (1024, 33878, 17398), (2048, 59188, 35734), (4096, 48851, 10701)],
    (0.99, 1.02): [(512, 12512, 1922), (1024, 21164, 2253), (2048, 33359, 6154), (4096, 25428, 14376)],
    (0.99, 1.05): [(512, 10638, 703), (1024, 15091, 1778), (2048, 17479, 1473), (4096, 14185, 3209)],
    (0.99, 1.2): [(512, 3188, 289), (1024, 4366, 485), (2048, 3695, 558), (4096, 3385, 672)],
}
ax = axs[0]
cols = {0.9: "C2", 0.98: "C0", 0.99: "C3"}
mk = {1.01: "v", 1.02: "o", 1.03: "o", 1.05: "s", 1.2: "^"}
for (p, x), rows in CHI.items():
    L, c, e = np.array(rows).T
    ax.errorbar(L, c, e, fmt=mk[x] + "-", color=cols[p], ms=3, capsize=1.2, lw=0.8)
Lg = np.array([1500, 4500])
ax.plot(Lg, 8e4 * (Lg / 1500) ** 1.75, "k--", lw=0.8)
ax.text(1400, 8e4, r"$L^{7/4}$", fontsize=7, ha="right", va="bottom")
from matplotlib.lines import Line2D
h1 = [Line2D([], [], color=c, lw=1.2, label=f"$p={p}$") for p, c in cols.items()]
h2 = [Line2D([], [], color="0.3", marker=m, ls="", ms=3, label=f"${x:g}\\,q^*$")
      for x, m in ((1.01, "v"), (1.02, "o"), (1.05, "s"), (1.2, "^"))]
ax.legend(handles=h1 + h2, ncol=4, fontsize=5.3, loc="upper left", handlelength=1.0, columnspacing=0.35,
          handletextpad=0.3)
ax.set(xscale="log", yscale="log", xlabel=r"$L$", ylabel=r"$\chi_L = L^2\langle m^2\rangle$",
       ylim=(800, 3e6), xlim=(400, 5500))
ax.set_xticks([512, 1024, 2048, 4096], ["512", "1024", "2048", "4096"])
label(ax, "(a)", x=0.88, y=0.08)

# (b) seeds: t P_s(t) for a 2x2 live block in the dead state
ax = axs[1]
qstar = {0.98: 0.00252, 0.99: 0.00127, 0.995: 0.00062}   # thresholds the runs were placed with
qtrue = {0.98: 0.00252, 0.99: 0.001282, 0.995: 0.00064}  # Table I (2026-10-31); legends use these


def pooled(pattern):
    files = sorted(glob.glob(os.path.join(RUNS, pattern)))
    t, data, meta = spread.load(files if len(files) > 1 else files[0])
    S = data.sum(0)
    return t, S[:, 1] / S[:, 0]


t, ps = pooled("E8a_spread/out/p1_a0.5_j2_q0_L2048_seed.txt")
ok = t >= 10
ax.plot(t[ok], (t * ps)[ok], "k-", lw=1.4, label="voter")
pts = [(0.98, 1.02, "C0", "-"), (0.98, 1.05, "C0", "--"), (0.99, 1.02, "C3", "-"), (0.99, 1.05, "C3", "--"),
       (0.995, 1.02, "C1", "-"), (0.995, 1.05, "C1", "--")]
for p, x, c, ls in pts:
    q = round(qstar[p] * x, 7)
    pat = f"E8a_long/out/p{p}_a0.5_j2_q{q:.7f}".rstrip("0") + "_L8192_T*_k*.txt"
    files = sorted(glob.glob(os.path.join(RUNS, pat)))
    if not files:
        print("missing", pat)
        continue
    t, data, meta = spread.load(files)
    S = data.sum(0)
    P = S[:, 1] / S[:, 0]
    ok = (t >= 10) & (P > 0)
    ax.plot(t[ok], (t * P)[ok], ls, color=c, label=f"$p={p}$, ${q / qtrue[p]:.2f}\\,q^*$")
ax.set(xscale="log", yscale="log", xlabel=r"$t$", ylabel=r"$t\,P_s(t)$", xlim=(10, 1.2e6), ylim=(0.5, 300))
ax.legend(fontsize=5.2, loc="upper left", handlelength=1.8)
label(ax, "(b)", x=0.88, y=0.95)

# (c) walls: width W(t) of the active layer (E8b; wall.py W, mean over the two walls and 2 seeds)
ax = axs[2]
qs = {0.9: 0.0119, 0.98: 0.00252, 0.99: 0.00127, 0.995: 0.00062}
pc = {0.9: "C2", 0.98: "C0", 0.99: "C3", 0.995: "C1"}
for p in qs:
    for x, ls in ((0.98, "-"), (1.05, "--")):
        q = round(qs[p] * x, 7)
        files = sorted(glob.glob(os.path.join(FIGDATA, "E8b_wall", f"p{p}_a0.5_j2_q{q:.7f}".rstrip("0")
                                              + "_Lx1024_Ly4096_s?.csv")))
        if not files:
            print("missing wall", p, x)
            continue
        runs = [np.array([[float(v) for v in r] for r in list(csv.reader(open(f)))[1:]]) for f in files]
        tt = runs[0][:, 0]
        W = np.mean([r[:, 6] for r in runs], axis=0)
        ok = tt >= 100
        ax.plot(tt[ok], W[ok], ls, color=pc[p], label=f"$p={p}$" if x == 0.98 else None)
ref = sorted(glob.glob(os.path.join(FIGDATA, "E16_wall_ref", "p1_a0.5_j2_q0_Lx1024_Ly4096_s?.csv")))
if ref:
    runs = [np.array([[float(v) for v in r] for r in list(csv.reader(open(f)))[1:]]) for f in ref]
    tt = runs[0][:, 0]
    W = np.mean([r[:, 6] for r in runs], axis=0)
    ok = tt >= 100
    ax.plot(tt[ok], W[ok], "k-", lw=1.4, label="voter")
ax.set(xscale="log", yscale="log", xlabel=r"$t$", ylabel=r"layer width $W$", xlim=(100, 1.2e6), ylim=(8, 2500))
ax.legend(fontsize=5.5, loc="upper left", bbox_to_anchor=(0.0, 0.9))
ax.text(0.97, 0.05, r"solid: $0.95$--$0.98\,q^*$; dashed: $1.02$--$1.05\,q^*$", transform=ax.transAxes, ha="right", fontsize=5.5)
label(ax, "(c)", x=0.04, y=0.97)
save(fig, "fig_voter")
