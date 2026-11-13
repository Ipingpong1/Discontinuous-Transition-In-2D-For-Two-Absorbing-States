"""Fig. wall: wetting of the wall between the two absorbing phases at p = 1/2 (E16_wall and, since
2026-10-31, R2_wall: Ly 2048 and 4096 series, longer runs; wall.py series and profiles reduced by
extract_prof.py).  (a) I(t), interface sites per unit wall length; (b) the mean interface profile across
the wall; (c) the excess adsorption Gamma = I - I(q=0) against delta = 1 - q/q*, for Lx = 256, 1024,
4096 and Ly = 1024, 2048, 4096, with the fits Gamma = A delta^-psi and A delta^-1/3 + B."""
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from common import FIGDATA, PAGE, label, save

QS = {0.5: float(os.environ.get("QSTAR05", 0.04597)), 0.9: float(os.environ.get("QSTAR09", 0.01203))}
TMIN = 3e5      # stationary window t >= max(TMIN, T/2) (E16: T = 1e6; R2_wall: T up to 3e6)


def read(f):
    rows = list(csv.reader(open(f)))
    return np.array([[float(v) for v in r] for r in rows[1:]])


def collect():
    """{(p, q, Lx, Ly): [runs]}; a R2_wall run replaces the E16 run of the same name (the Lx 4096 seeds 1, 2
    were continued from their E16 checkpoints to 2e6 steps)."""
    runs = {}
    for exp in ("E16_wall", "R2_wall"):
        for f in sorted(glob.glob(os.path.join(FIGDATA, exp, "*.csv"))):
            if f.endswith("_prof.csv"):
                continue
            m = re.match(r"p([\d.]+)_a0.5_j2_q([\d.]+)_Lx(\d+)_Ly(\d+)_s(\d)", os.path.basename(f))
            d = read(f)
            pr = read(f[:-4] + "_prof.csv")
            runs[os.path.basename(f)] = ((float(m[1]), float(m[2]), int(m[3]), int(m[4])),
                                         dict(t=d[:, 0], I=d[:, 1], W=d[:, 6], rc=d[:, 7], alive=d[:, 5], y=pr[:, 0],
                                              prof=pr[:, 3], T=d[-1, 0]))
    g = defaultdict(list)
    for k, r in runs.values():
        g[k].append(r)
    return g


def stationary(rs):
    """Late-window means of I and W per run.  A wetted wall wanders; a run in which it has drifted so far that
    an absorbing slab is being consumed (alive fraction outside 0.5 +- 0.2) is cut at the first such record
    and dropped if fewer than three records of its window remain."""
    vals = []
    for r in rs:
        sel = r["t"] >= (r["T"] / 2 if r["T"] > 1.1e6 else TMIN)   # E16 (T 1e6): [3e5, 1e6] as published
        bad = np.flatnonzero(np.abs(r["alive"] - 0.5) >= 0.2)    # a wall more than Ly/5 from its row leaves the analysis window
        if bad.size:
            sel &= r["t"] < r["t"][bad[0]]
        if sel.sum() >= 3 and r["T"] >= 9e5:
            vals.append((r["I"][sel].mean(), r["W"][sel].mean()))
    if not vals:
        return np.nan, np.nan, np.nan, 0
    v = np.array(vals)
    return v[:, 0].mean(), v[:, 0].std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan, v[:, 1].mean(), len(v)


if __name__ == "__main__":
    g = collect()
    I0 = stationary(g.get((0.5, 0.0, 1024, 1024), []))[0]
    I0b = {Ly: stationary(g.get((0.5, 0.0, 1024, Ly), []))[0] for Ly in (1024, 2048) if (0.5, 0.0, 1024, Ly) in g}
    print(f"bare wall I0 = {I0:.2f}  (per Ly: {I0b})")
    fig, axs = plt.subplots(1, 3, figsize=(PAGE, 2.2), gridspec_kw=dict(wspace=0.32))
    xs = sorted({k[1] / QS[0.5] for k in g if k[0] == 0.5 and k[1] > 0 and k[3] == 1024})
    cmap = plt.get_cmap("plasma")
    for i, x in enumerate(xs):
        c = cmap(i / max(1, len(xs) - 1) * 0.85)
        k = min((kk for kk in g if kk[0] == 0.5 and kk[2] == 1024 and kk[3] == 1024 and kk[1] > 0),
                key=lambda kk: abs(kk[1] / QS[0.5] - x))
        rs = g.get(k, [])
        if not rs:
            continue
        t = rs[0]["t"]
        I = np.mean([r["I"] for r in rs], axis=0)
        ok = t >= 100
        axs[0].plot(t[ok], I[ok], color=c, label=f"${x:.3f}$")
        prof = np.mean([r["prof"] for r in rs], axis=0)
        axs[1].plot(rs[0]["y"], prof, color=c)
    if (0.5, 0.0, 1024, 1024) in g:
        rs = g[(0.5, 0.0, 1024, 1024)]
        t = rs[0]["t"]
        axs[0].plot(t[t >= 100], np.mean([r["I"] for r in rs], axis=0)[t >= 100], "k-", label="$q=0$")
        axs[1].plot(rs[0]["y"], np.mean([r["prof"] for r in rs], axis=0), "k-")
    axs[0].set(xscale="log", xlabel=r"$t$", ylabel=r"interface sites per unit length $I$")
    axs[0].legend(fontsize=5.5, title=r"$q/q^*$", title_fontsize=6, loc="upper left", ncol=2)
    label(axs[0], "(a)", x=0.88, y=0.1)
    axs[1].set(xlabel=r"distance from the layer centre $y$", ylabel=r"interface-site fraction", xlim=(-160, 160))
    label(axs[1], "(b)")
    ax = axs[2]
    mk = {256: "v", 1024: "o", 4096: "s"}
    lyc = {1024: "k", 2048: "C0", 4096: "C3"}
    fitx, fity, fite = [], [], []
    for Ly in (1024, 2048, 4096):
        for Lx in (256, 1024, 4096):
            pts = []
            for k, rs in g.items():
                if k[0] == 0.5 and k[2] == Lx and k[3] == Ly and k[1] > 0:
                    if Lx == 4096 and Ly == 1024 and 1 - k[1] / QS[0.5] < 0.03:
                        continue   # box too small: the wandering wall's excursions reach the other wall at delta <= 0.02
                    I, e, W, n = stationary(rs)
                    if n:
                        pts.append((1 - k[1] / QS[0.5], I - I0b.get(Ly, I0), e, W, n))
            if not pts:
                continue
            pts = np.array(sorted(pts))
            ax.errorbar(pts[:, 0], pts[:, 1], pts[:, 2], fmt=mk[Lx], color=lyc[Ly], mfc="w" if Ly == 1024 else lyc[Ly],
                        ms=3.2, capsize=1.2, lw=0.7, ls="-" if (Lx, Ly) in ((1024, 1024), (1024, 2048)) else "")
            for row in pts:
                print(f"Lx={Lx} Ly={Ly} delta={row[0]:.4f} Gamma={row[1]:.2f}+-{row[2]:.2f} W={row[3]:.1f} n={row[4]:.0f}")
            if Lx == 1024 and Ly >= 2048:
                fitx += list(pts[:, 0]); fity += list(pts[:, 1]); fite += list(np.nan_to_num(pts[:, 2], nan=1.0))
    fitx, fity, fite = np.array(fitx), np.array(fity), np.maximum(np.array(fite), 0.3)
    if len(fitx) >= 4:
        sel = fitx <= 0.1001
        w = 1 / fite[sel] ** 2
        A = np.vstack([np.ones(sel.sum()), -np.log(fitx[sel])]).T
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        lnA, psi = cov @ (A.T @ (w * np.log(fity[sel])))
        xx = np.logspace(np.log10(fitx.min() * 0.8), np.log10(0.12), 30)
        ax.plot(xx, np.exp(lnA) * xx ** -psi, "k--", lw=0.7, label=fr"$\delta^{{-{psi:.2f}}}$")
        fit_labels = []
        X = np.vstack([fitx[sel] ** (-1 / 3), np.ones(sel.sum())]).T
        cov3 = np.linalg.inv(X.T @ (X * w[:, None]))
        A3, B3 = cov3 @ (X.T @ (w * fity[sel]))
        ax.plot(xx, A3 * xx ** (-1 / 3) + B3, "k:", lw=0.8, label=fr"${A3:.1f}\,\delta^{{-1/3}} {B3:+.0f}$")
        chi_p = (w * (np.log(fity[sel]) - lnA + psi * np.log(fitx[sel])) ** 2).sum()
        chi_3 = (w * (fity[sel] - A3 * fitx[sel] ** (-1 / 3) - B3) ** 2).sum()
        print(f"fits over delta <= 0.1, Lx 1024, Ly >= 2048 ({sel.sum()} points): Gamma ~ delta^-{psi:.3f} +- {np.sqrt(cov[1, 1]):.3f} "
              f"(chi2 {chi_p:.1f}); {A3:.2f} delta^-1/3 + {B3:.2f} (chi2 {chi_3:.1f})")
    ax.set(xscale="log", yscale="log", xlabel=r"$\delta = 1 - q/q^*$", ylabel=r"excess $\Gamma = I - I_{q=0}$")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color="0.3", marker=mk[Lx], ls="", ms=3.2, mfc="0.3", label=(r"$L_x=$" if Lx == 256 else "") + f"${Lx}$")
         for Lx in (256, 1024, 4096)]
    h += [Line2D([], [], color=lyc[Ly], marker="o", ls="", ms=3.2, mfc="w" if Ly == 1024 else lyc[Ly],
                 label=(r"$L_y=$" if Ly == 1024 else "") + f"${Ly}$") for Ly in (1024, 2048, 4096)]
    h += [l for l in ax.get_lines() if l.get_label().startswith("$")]
    ax.legend(handles=h, fontsize=5, ncol=2, loc="lower left", handlelength=1.3, columnspacing=0.5, labelspacing=0.25,
              handletextpad=0.4, borderpad=0.3)
    label(ax, "(c)", x=0.85, y=0.95)
    save(fig, "fig_wall05")
