#!/usr/bin/env python3
"""kappa*/alpha from voter-model perturbation theory vs walk time T (patterns.txt), against
the measured q*/((1-p) alpha) at alpha = 0.5, p = 0.9-0.995."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pa
import vpt

d = vpt.read_patterns("patterns.txt")
Ts = [T for T in sorted(d) if T >= 1]
th = [vpt.thresholds(vpt.drift_basis(*d[T][::-1])) for T in Ts]
meas = [q / ((1 - p) * al) for p, al, q in pa.SIM if p >= 0.9]

INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=150)
ax.grid(True, color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.axhspan(min(meas), max(meas), color="#2a78d6", alpha=0.15, lw=0,
           label="simulation, p = 0.9–0.995 (range)")
ax.plot(Ts, [t["kappa_a"] for t in th], color="#eb6834", lw=2, label="f′(0) = 0  (absorbing loses stability)")
ax.plot(Ts, [t["kappa_front"] for t in th], color="#1baf7a", lw=2, label="∫ f du over [0, ½] = 0  (stationary front)")
ax.set_xscale("log")
ax.set_ylim(0.1, 0.6)
ax.set_xlabel("dual walk time T  (~ 1/(1−p))")
ax.set_ylabel("κ*/α = q* / ((1−p) α)")
ax.set_title("Voter-point perturbation theory", color=INK, loc="left")
ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="upper center",
          bbox_to_anchor=(0.5, -0.16), ncol=1)
fig.tight_layout()
fig.savefig("vpt_kappa.png")
