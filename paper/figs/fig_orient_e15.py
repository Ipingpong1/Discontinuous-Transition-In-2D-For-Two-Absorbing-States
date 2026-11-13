"""(E15, superseded 2026-10-31 by the direct front tracking of fig_orient.py / campaign/round2) Fig. orient: front velocity of bands with normal direction (1,0), (2,1), (1,1) (E15_orient, L 2048,
T 1e5, 4 seeds; campaign/tilt.py) at p = 0.1, 0.5, 0.9, and the zero crossings."""
import glob
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from common import PAGE, RUNS, label, save
import tilt

QREF = {0.1: 0.0368, 0.5: 0.0459, 0.9: 0.0119}
AB = {"1-0": ("o", "C0", r"$(1,0)$, $0^\circ$"), "2-1": ("s", "C1", r"$(2,1)$, $26.6^\circ$"),
      "1-1": ("^", "C3", r"$(1,1)$, $45^\circ$")}


def collect():
    g = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(RUNS, "E15_orient", "out", "*.txt"))):
        if not any(ln.startswith("# done") for ln in open(f)):
            continue
        m = re.match(r"p([\d.]+)_a0.5_q([\d.]+)_ab(\d-\d)_L\d+_s(\d)", os.path.basename(f))
        v, rho, te, n = tilt.vel(f)
        g[(float(m[1]), m[3], float(m[2]))].append(v)
    return g


def crossing(x, v, se):
    """Zero of the weighted straight line through the points with |x - 1| <= 0.02 (x*, its error)."""
    sel = np.isfinite(v) & (np.abs(x - 1) <= 0.0201)
    w = 1 / np.maximum(se[sel], 1e-5) ** 2
    A = np.vstack([np.ones(sel.sum()), x[sel] - 1]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    a, b = cov @ (A.T @ (w * v[sel]))
    xs = 1 - a / b
    J = np.array([-1 / b, a / b ** 2])
    return xs, float(np.sqrt(J @ cov @ J)), b


if __name__ == "__main__":
    g = collect()
    fig, axs = plt.subplots(1, 3, figsize=(PAGE, 2.1), gridspec_kw=dict(wspace=0.3))
    for ax, p in zip(axs, sorted(QREF)):
        for ab, (mk, c, lab) in AB.items():
            ks = sorted(k for k in g if k[0] == p and k[1] == ab)
            if not ks:
                continue
            x = np.array([k[2] / QREF[p] for k in ks])
            vs = [np.array(g[k]) for k in ks]
            vm = np.array([np.nanmean(v) for v in vs])
            se = np.array([np.nanstd(v, ddof=1) / np.sqrt(np.isfinite(v).sum()) if len(v) > 1 else np.nan for v in vs])
            n = [len(v) for v in vs]
            ax.errorbar(x, vm * 1e3, se * 1e3, fmt=mk, color=c, ms=3, capsize=1.2, label=lab)
            if len(x) >= 3 and np.isfinite(se).all():
                xs, dxs, b = crossing(x, vm, se)
                xx = np.linspace(0.965, 1.035, 10)
                ax.plot(xx, b * (xx - xs) * 1e3, "-", color=c, lw=0.6)
                print(f"p={p} ab={ab}: x*={xs:.4f} +- {dxs:.4f}  slope dv/dx={b:.4f}  n={n}")
        ax.axhline(0, color="0.6", lw=0.5)
        ax.axvline(1, color="0.6", lw=0.5, ls=":")
        ax.set(xlabel=r"$q/q^*_{\rm ref}$", title=f"$p={p}$")
        if p == 0.1:
            ax.set_ylabel(r"normal front velocity $v\ (10^{-3})$")
            ax.legend(fontsize=6, loc="upper left")
    save(fig, "fig_orient")
