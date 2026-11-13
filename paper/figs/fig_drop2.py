"""Fig. drop2: the two-droplet coexistence test at p = 1/2 (R2_drop2, L 4096, 4 seeds).  Radius of an active
disk in the dead state (left) and of a dead hole in the active state (right), both of initial radius
1000, at q/q* = 0.994, 0.997, 1.003, 1.006; R from the row profile (drop2.py).  Dashed: the growth law
dR/dt = +-v - D/R with v = 0.278 (q/q* - 0.999), the local slope and crossing of the (1,0) band of the
orientation test (R2_orient, L 4096), and D = 0.20 (Fig. 4)."""
import glob
import os
import re
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from common import COL, RUNS, label, save
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "campaign", "round2"))
import drop2  # noqa: E402

QS, SLOPE, X0, D, R0 = 0.04596, 0.278, 0.999, 0.20, 1000   # slope and crossing from the (1,0) bands, L 4096
NAME = re.compile(r"p0\.5_a0\.5_j2_q([\d.]+)_(disk|hole)R(\d+)_L(\d+)_s(\d+)$")
COLS = {0.994: "C0", 0.997: "C9", 1.003: "C1", 1.006: "C3"}


def main():
    g = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(RUNS, "R2_drop2", "out", "*.txt"))):
        prefix = f[:-4]
        m = NAME.search(os.path.basename(prefix))
        if not m or not os.path.exists(prefix + ".prof"):
            continue
        q, kind = float(m[1]), m[2]
        L, ts, prof, done = drop2.load(prefix)
        t, R_ts, R_pr = drop2.radii(L, ts, prof, kind, R0)
        g[(round(q / QS, 3), kind)].append(R_pr)
    fig, axs = plt.subplots(1, 2, figsize=(2 * COL, 2.1), gridspec_kw=dict(wspace=0.28), sharey=True)
    for ax, kind in zip(axs, ("disk", "hole")):
        for x in sorted(COLS):
            rs = g.get((x, kind), [])
            if not rs:
                continue
            n = min(len(r[0]) for r in rs)
            t = rs[0][0][:n]
            Rm = np.mean([r[1][:n] for r in rs], axis=0)
            for r in rs:
                ax.plot(r[0] / 1e3, r[1], color=COLS[x], lw=0.4, alpha=0.4)
            ax.plot(t / 1e3, Rm, color=COLS[x], lw=1.2, label=f"${x:.3f}$")
            v = SLOPE * (x - X0)
            rate = v - D / R0 if kind == "disk" else -(v + D / R0)
            sel = t >= 1e4
            ax.plot(t[sel] / 1e3, Rm[sel][0] + rate * (t[sel] - t[sel][0]), "--", color="k", lw=0.6)
        ax.set(xlabel=r"$t$ ($10^3$ steps)", xlim=(0, 200))
        ax.axhline(R0, color="0.7", lw=0.5, ls=":")
    axs[0].set(ylabel=r"radius $R$", title="active disk in the dead state")
    axs[1].set(title="dead hole in the active state")
    axs[0].legend(fontsize=6, title=r"$q/q^*$", title_fontsize=6.5, loc="upper left", ncol=2)
    label(axs[0], "(a)", x=0.9, y=0.1)
    label(axs[1], "(b)", x=0.9, y=0.1)
    save(fig, "fig_drop2")


if __name__ == "__main__":
    main()
