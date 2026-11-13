"""Fig. rough: the edges of the wetting layer at p = 1/2 (R2_wall snapshots, rough.py reductions in
data/figdata/R2).  (a) Height-difference correlation C(r) = <[h(x+r) - h(x)]^2> of the layer edges (mean of
the lower and upper edge, seeds averaged) at delta = 0.011 for Lx = 256, 1024, 4096 (Ly 2048), and of a free
front (the edges of a (1,0) band at q*, L 4096); guide: C ~ r (zeta = 1/2).  (b) Distribution of the local
layer width w (per column) at Lx 1024, Ly 2048 for delta = 0.051, 0.021, 0.011, 0.0065."""
import csv
import glob
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from common import COL, FIGDATA, label, save

from common import QSTAR05 as QS
R2 = os.path.join(FIGDATA, "R2")


def load_rough():
    g = defaultdict(list)
    for f in glob.glob(os.path.join(R2, "rough_*.csv")):
        rows = list(csv.DictReader(open(f)))
        if not rows:
            continue
        name = rows[0]["prefix"]
        m = re.match(r"p0\.5_a0\.5_j2_q([\d.]+)_Lx(\d+)_Ly(\d+)_s(\d+)$", name)
        mb = re.match(r"p([\d.]+)_a0\.5_q([\d.]+)_ab1-0_L(\d+)_s(\d+)$", name)
        if m:
            key = ("wall", float(m[1]), int(m[2]), int(m[3]))
        elif mb:
            key = ("band", float(mb[1]), float(mb[2]), int(mb[3]))
        else:
            continue
        r = np.array([int(x["r"]) for x in rows])
        C = 0.5 * (np.array([float(x["C_lo"]) for x in rows]) + np.array([float(x["C_hi"]) for x in rows]))
        g[key].append((r, C, float(rows[0]["sig_lo"]), float(rows[0]["sig_hi"]), float(rows[0]["w_mean"]),
                       float(rows[0]["w_sd"])))
    return g


def load_hist():
    g = defaultdict(lambda: defaultdict(float))
    n = defaultdict(int)
    for f in glob.glob(os.path.join(R2, "whist_*.csv")):
        rows = list(csv.DictReader(open(f)))
        if not rows:
            continue
        m = re.match(r"p0\.5_a0\.5_j2_q([\d.]+)_Lx(\d+)_Ly(\d+)_s(\d+)$", rows[0]["prefix"])
        if not m:
            continue
        key = (float(m[1]), int(m[2]), int(m[3]))
        n[key] += 1
        for x in rows:
            g[key][int(x["w"])] += float(x["P"])
    return {k: {w: p / n[k] for w, p in v.items()} for k, v in g.items()}


if __name__ == "__main__":
    g = load_rough()
    fig, (a, b) = plt.subplots(1, 2, figsize=(2 * COL, 2.15), gridspec_kw=dict(wspace=0.3))
    q011 = 0.045441
    for Lx, mk, c in ((256, "v", "C2"), (1024, "o", "C0"), (4096, "s", "C3")):
        rs = g.get(("wall", q011, Lx, 2048), [])
        if not rs:
            continue
        r = rs[0][0]
        C = np.mean([x[1] for x in rs], axis=0)
        a.plot(r, C, mk + "-", color=c, ms=3, lw=0.8, label=f"layer edges, $L_x = {Lx}$")
        print(f"delta 0.011 Lx {Lx} Ly 2048: n={len(rs)} sigma_lo={np.mean([x[2] for x in rs]):.1f} sigma_hi={np.mean([x[3] for x in rs]):.1f} "
              f"w={np.mean([x[4] for x in rs]):.1f} sd={np.mean([x[5] for x in rs]):.1f}; zeta_loc(r=2-8)={np.polyfit(np.log(r[(r >= 2) & (r <= 8)]), np.log(C[(r >= 2) & (r <= 8)]), 1)[0] / 2:.3f}")
    bands = [k for k in g if k[0] == "band" and k[1] == 0.5]
    for k in bands:
        rs = g[k]
        r = rs[0][0]
        C = np.mean([x[1] for x in rs], axis=0)
        a.plot(r, C, "k^-", ms=3, lw=0.8, mfc="w", label=f"free front, $L = {k[3]}$")
        print(f"free front p=0.5 L={k[3]}: n={len(rs)} sigma={np.mean([x[2] for x in rs]):.1f}/{np.mean([x[3] for x in rs]):.1f}; "
              f"zeta_loc(r=2-8)={np.polyfit(np.log(r[(r >= 2) & (r <= 8)]), np.log(C[(r >= 2) & (r <= 8)]), 1)[0] / 2:.3f}")
    rr = np.array([1.5, 30])
    a.plot(rr, 25 * rr, "k--", lw=0.6)
    a.text(8, 120, r"$r^{1}$", fontsize=7)
    a.set(xscale="log", yscale="log", xlabel=r"$r$", ylabel=r"$\langle[h(x+r)-h(x)]^2\rangle$")
    a.legend(fontsize=5.5, loc="lower right")
    label(a, "(a)")
    h = load_hist()
    cmap = plt.get_cmap("plasma")
    qs = [0.043605, 0.044982, 0.045441, 0.0456705]
    for i, q in enumerate(qs):
        d = h.get((q, 1024, 2048))
        if not d:
            continue
        w = np.array(sorted(d))
        P = np.array([d[x] for x in w])
        sm = np.convolve(P, np.ones(5) / 5, "same")
        b.plot(w, sm, color=cmap(0.85 * i / (len(qs) - 1)), lw=0.9, label=f"${1 - q / QS:.3f}$")
    b.set(xlabel=r"local layer width $w$", ylabel=r"$P(w)$", xlim=(0, 600), yscale="log", ylim=(1e-5, 0.05))
    b.legend(fontsize=6, title=r"$\delta$", title_fontsize=6.5)
    label(b, "(b)", x=0.9)
    save(fig, "fig_rough")
