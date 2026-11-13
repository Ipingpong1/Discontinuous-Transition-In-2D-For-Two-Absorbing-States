#!/usr/bin/env python3
"""q*(p) at alpha = 0.5: simulation vs the PA front criterion (curve_a05.csv)."""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pa

rows = list(csv.DictReader(open("curve_a05.csv")))
p = [float(r["p"]) for r in rows]
get = lambda k: [float(r[k]) if r[k] else float("nan") for r in rows]
sim = [(pp, q) for pp, al, q in pa.SIM if al == 0.5]

INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED})
fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=150)
ax.grid(True, color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.plot(p, get("disc_q_front"), color="#eb6834", lw=2, label="PA front criterion, synchronous map")
ax.plot(p, get("cont_q_front"), color="#1baf7a", lw=2, label="PA front criterion, continuous time")
ax.plot(*zip(*sim), "o", ms=7, color="#2a78d6", mec="white", mew=1.5, label="simulation q* (fronts)")
ax.set_xlabel("voter mixing p")
ax.set_ylabel("threshold q")
ax.set_title("Noisy D&N threshold, α = 0.5", color=INK, loc="left")
ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
fig.tight_layout()
fig.savefig("pa_vs_sim.png")
