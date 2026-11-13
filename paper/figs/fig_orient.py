"""Fig. orient: normal front velocity of bands with normal direction (1,0), (2,1), (1,1) from
direct front tracking (R2_orient, L 4096, T 2e5, 8 seeds; orient2.py --csv) at p = 0.1, 0.5, 0.9, with the
zero crossings.  Reads data/figdata/R2/orient_runs.csv (or the file given as argv[1])."""
import csv
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from common import FIGDATA, PAGE, label, save

AB = {"1-0": ("o", "C0", r"$(1,0)$, $0^\circ$"), "2-1": ("s", "C1", r"$(2,1)$, $26.6^\circ$"),
      "1-1": ("^", "C3", r"$(1,1)$, $45^\circ$")}
KEY = "v_out"


def crossing(x, v, se):
    w = 1 / se ** 2
    A = np.vstack([np.ones(len(x)), x - 1]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    a, b = cov @ (A.T @ (w * v))
    xs = 1 - a / b
    J = np.array([-1 / b, a / b ** 2])
    return xs, float(np.sqrt(J @ cov @ J)), b


def main(path, L):
    g = defaultdict(list)
    for r in csv.DictReader(open(path)):
        if int(r["L"]) == L and r[KEY] not in ("nan", ""):
            g[(float(r["p"]), r["ab"], round(float(r["x"]), 4))].append(float(r[KEY]))
    ps = sorted({k[0] for k in g})
    fig, axs = plt.subplots(1, len(ps), figsize=(PAGE, 2.1), gridspec_kw=dict(wspace=0.3), squeeze=False)
    for ax, p in zip(axs[0], ps):
        for ab, (mk, c, lab) in AB.items():
            ks = sorted(k for k in g if k[0] == p and k[1] == ab)
            if not ks:
                continue
            x = np.array([k[2] for k in ks])
            vm = np.array([np.mean(g[k]) for k in ks])
            se = np.array([np.std(g[k], ddof=1) / np.sqrt(len(g[k])) for k in ks])
            ax.errorbar(x, vm * 1e3, se * 1e3, fmt=mk, color=c, ms=3, capsize=1.2, label=lab)
            # the line is fitted to x >= 0.995: at 0.99 q* the band interior nucleates the absorbing
            # phase within the fit window (lifetime ~ 7e4 steps at p = 1/2), which speeds the retreat
            fit = x >= 0.9949
            xs, dxs, b = crossing(x[fit], vm[fit], se[fit])
            xx = np.linspace(0.987, 1.013, 10)
            ax.plot(xx, b * (xx - xs) * 1e3, "-", color=c, lw=0.6)
            print(f"p={p} L={L} ab={ab}: x*={xs:.4f} +- {dxs:.4f}  slope={b:.3f}  n={[len(g[k]) for k in ks]}")
        ax.axhline(0, color="0.6", lw=0.5)
        ax.axvline(1, color="0.6", lw=0.5, ls=":")
        ax.set(xlabel=r"$q/q^*$", title=f"$p={p}$", xlim=(0.987, 1.013))
        if p == 0.5:
            ax.set_ylim(-4.5, 3.6)
        ax.set_xticks([0.99, 1.0, 1.01])
    axs[0][0].set_ylabel(r"normal front velocity $v\ (10^{-3})$")
    axs[0][0].legend(fontsize=6, loc="upper left")
    for ax, s in zip(axs[0], "abc"):
        label(ax, f"({s})", x=0.86, y=0.1)
    save(fig, "fig_orient" if L == 4096 else f"fig_orient_L{L}")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(FIGDATA, "R2", "orient_runs.csv")
    L = int(sys.argv[2]) if len(sys.argv) > 2 else 4096
    main(path, L)
