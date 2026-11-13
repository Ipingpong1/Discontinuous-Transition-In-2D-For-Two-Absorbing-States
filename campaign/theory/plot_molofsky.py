#!/usr/bin/env python3
"""Figure for MOLOFSKY.md: the lineage closure applied to the five-cell rule of Molofsky et al.
Reads molofsky.csv (from `python3 vpt.py molofsky --load patterns5.txt --boot 300 --csv molofsky.csv`)
and patterns.txt (Moore run, for the >= 5-lineage fraction).  Writes molofsky.png.
    python3 plot_molofsky.py [molofsky.csv]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import vpt

HERE = Path(__file__).resolve().parent
csv = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "molofsky.csv"
d = np.genfromtxt(csv, delimiter=",", names=True)
d = d[d["T"] > 0]
T0 = 0.45                                             # where the T = 0 (mean-field) row is drawn
C_A, C_F, C_H = "#2a78d6", "#eb6834", "#1baf7a"       # validated categorical slots 1-3
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#dcdbd6"
th_a0, th_f0, th_h0 = vpt.MF_RAYS
dyad = vpt.DYAD_RAY
chord = np.degrees(np.arctan2(1 - 0.4, 0.024 - 0.2))  # voter point -> Molofsky's (0.024, 1)

# Moore (D&N) >= 5-lineage fraction among >= 3-class patterns, for comparison
moore = vpt.read_patterns(HERE / "patterns.txt")
mT, m5 = [], []
for T in sorted(moore):
    trials, pats = moore[T]
    n = sum(c for c, _ in pats)
    if T > 0 and n:
        mT.append(T); m5.append(sum(c for c, ns in pats if len(ns) >= 4) / n)

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 9.5,
                     "axes.titleweight": "bold", "axes.titlelocation": "left"})
fig, ax = plt.subplots(2, 2, figsize=(8.4, 6.6), dpi=160)
for a in ax.flat:
    a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
    for s in ("top", "right"): a.spines[s].set_visible(False)

# (a) the three rays
a = ax[0, 0]
a.set_xscale("log")
a.axhline(dyad, color=MUTED, lw=1, ls="--", label="their dyad rule 3p₁ + p₂ = 1 (108.4°)")
a.axhline(chord, color=MUTED, lw=1, ls=":", label="chord to their point (0.024, 1) (106.4°)")
for key, col, lab, mk in (("th_a", C_A, "u = 0 marginal (absorbing loses stability)", "o"),
                          ("th_f", C_F, "stationary front", "s"),
                          ("th_h", C_H, "u = ½ marginal (active loses stability)", "^")):
    a.plot(d["T"], d[key], color=col, lw=2, marker=mk, ms=3.5, label=lab)
a.plot(d["T"], d["th_3"], color=INK, lw=1, ls="-.", label="three-lineage patterns only")
a.plot([T0, T0, T0], [th_a0, th_f0, th_h0], ls="none", marker="D", ms=4,
       markerfacecolor="white", markeredgecolor=MUTED)
a.text(T0 * 1.25, th_h0 + 0.5, "mean field (T = 0): 90°, 111.8°, 123.7°", color=MUTED, fontsize=7.5)
a.set_xlim(0.3, 1.3e5); a.set_ylim(88, 126)
a.set_xlabel("walk time T (steps)"); a.set_ylabel("direction in the (p₁ − 1/5, p₂ − 2/5) plane (°)")
a.set_title("(a) The three rays of the closure vs walk time")
a.legend(fontsize=7, loc="lower right", frameon=False)

# (b) sector width
a = ax[0, 1]
a.set_xscale("log"); a.set_yscale("log")
a.errorbar(d["T"], d["width"], yerr=np.nan_to_num(d["sd_w"]), color=C_F, lw=2, marker="s", ms=3.5, capsize=0)
a.plot([T0], [th_h0 - th_a0], ls="none", marker="D", ms=4, markerfacecolor="white", markeredgecolor=MUTED)
a.text(T0 * 1.6, 36, "mean-field wedge 33.7°", color=MUTED, fontsize=7.5)
a.set_xlim(0.3, 1.3e5); a.set_ylim(0.05, 60)
a.set_xlabel("walk time T (steps)"); a.set_ylabel("angular width of the bistable sector (°)")
a.set_title("(b) Angular width of the bistable sector")

# (c) five-lineage fractions, both rules
a = ax[1, 0]
a.set_xscale("log"); a.set_yscale("log")
a.plot(d["T"], d["f5"], color=C_F, lw=2, marker="s", ms=3.5, label="five-cell rule: P(5 classes)/P(≥3)")
a.plot(mT, m5, color=C_A, lw=2, marker="o", ms=3.5, label="D&N, Moore: P(≥5 classes)/P(≥3)")
a.set_xlim(0.7, 1.3e5)
a.set_xlabel("walk time T (steps)"); a.set_ylabel("fraction of ≥3-class patterns")
a.set_title("(c) Fraction of patterns with five or more lineages")
a.legend(fontsize=7.5, loc="lower left", frameon=False)

# (d) the (p1, p2) plane
a = ax[1, 1]
v1, v2 = 0.2, 0.4
L = 0.7
def wedge(t1, t2, **kw):
    ts = np.radians(np.linspace(t1, t2, 50))
    a.fill(np.r_[v1, v1 + L * np.cos(ts)], np.r_[v2, v2 + L * np.sin(ts)], **kw)
wedge(th_a0, th_h0, color="#e9e8e3", lw=0, label="mean-field multiple outcome")
for T, al in ((16, 0.35), (1024, 0.9)):
    r = d[d["T"] == T][0]
    wedge(r["th_a"], r["th_h"], color=C_F, alpha=al, lw=0, label=f"closure, T = {T} ({r['width']:.1f}°)")
ts = np.radians(np.linspace(0, 360, 400))
for ang, ls in ((dyad, "--"), (chord, ":")):
    a.plot([v1 + L * np.cos(np.radians(ang)), v1 - 0.3 * np.cos(np.radians(ang))],
           [v2 + L * np.sin(np.radians(ang)), v2 - 0.3 * np.sin(np.radians(ang))], color=MUTED, lw=1, ls=ls)
a.plot([0, 0.2], [0.4, 0.4], color=INK, lw=2.5, solid_capstyle="butt")
a.plot([0.2, 0.2], [0.5, 1.0], color=INK, lw=2.5, solid_capstyle="butt")
a.text(0.02, 0.36, "transect p₂ = 0.4: clustering", fontsize=7, color=INK, va="top")
a.text(0.215, 0.52, "transect p₁ = 0.2:\nergodic", fontsize=7, color=INK)
a.plot([0.024, 0.27, 0.31, 0.35], [1, 0, 0, 0], ls="none", marker="o", ms=4, color=INK)
a.annotate("their boundary\npoint (0.024, 1)", (0.024, 1.0), (0.075, 0.84), fontsize=7, color=INK,
           arrowprops=dict(arrowstyle="-", color=INK, lw=0.7))
a.text(0.25, 0.04, "their p₂ = 0 points:\n0.27 clustering, 0.31 phase\nseparation, 0.35 ergodic", fontsize=6.3, color=INK)
a.plot([v1], [v2], marker="*", ms=9, color=INK, ls="none")
a.text(0.21, 0.405, "voter point", fontsize=7, color=INK)
a.set_xlim(0, 0.5); a.set_ylim(0, 1.02)
a.set_xlabel("p₁"); a.set_ylabel("p₂")
a.set_title("(d) Molofsky's (p₁, p₂) plane")
a.legend(fontsize=6.8, loc="upper right", frameon=False)

fig.tight_layout(h_pad=1.5, w_pad=1.0)
out = HERE / "molofsky.png"
fig.savefig(out)
print("wrote", out)
