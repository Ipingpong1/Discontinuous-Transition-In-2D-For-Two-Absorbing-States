#!/usr/bin/env python3
"""strips.png: (top) strip width W(t)/Ly of every run at p = 0.995 and p = 0.99, coloured by q,
with the crossover time t_x marked; (middle) the interior interface density (central 256 rows,
relative to its t ~ 1e3 value); (bottom) seed-mean front velocity against q in the early window
[2e3, 1e5] (the Table I estimator's range, open symbols) and in the late window [t_x, T] (filled),
with the old q* (dashed) and the late-window line fit (solid).  Widths from the outer-envelope
estimator (analyze.py).

    python3 plot_strips.py          # writes strips.png next to this script
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, D)
import analyze  # noqa: E402

PS = (0.995, 0.99)
fig, axes = plt.subplots(3, 2, figsize=(12, 11))
rows = analyze.strips([f"{D}/strips", f"{D}/wide", f"{D}/tall"], out=open(os.devnull, "w"))
res = analyze.cross(rows, out=open(os.devnull, "w"), nboot=200)
for j, p in enumerate(PS):
    runs = sorted(glob.glob(f"{D}/strips/out/p{p}_*.txt") + glob.glob(f"{D}/wide/out/p{p}_*.txt")
                  + glob.glob(f"{D}/tall/out/p{p}_*.txt"))
    qs = sorted(set(float(analyze.fronts.params(os.path.basename(r)[:-4])["q"]) for r in runs))
    cmap = plt.get_cmap("viridis")
    col = {q: cmap(i / max(len(qs) - 1, 1)) for i, q in enumerate(qs)}
    for txt in runs:
        if not os.path.exists(txt[:-4] + ".prof") or os.path.getsize(txt[:-4] + ".prof") == 0:
            continue
        P = analyze.fronts.params(os.path.basename(txt)[:-4])
        q, Lx, Ly = float(P["q"]), int(P["Lx"]), int(P["Ly"])
        try:
            t, Wd, cen, Lx, Ly, rho_b = analyze.track(txt); W = Wd["o"]
        except Exception:
            continue
        ls = "-" if (Lx, Ly) == (2048, 4096) else ("--" if Lx == 4096 else ":")
        axes[0, j].plot(t, W / Ly, ls, color=col[q], lw=0.9, alpha=0.85)
        axes[1, j].plot(t, cen, ls, color=col[q], lw=0.9, alpha=0.85)
    for q in qs:
        axes[0, j].plot([], [], "-", color=col[q], label=f"q = {q:g}")
    axes[0, j].plot([], [], "k-", label="Lx 2048 x Ly 4096"); axes[0, j].plot([], [], "k--", label="Lx 4096 x Ly 4096")
    axes[0, j].plot([], [], "k:", label="Lx 1024 x Ly 8192, strip Ly/4")
    for ax in axes[:2, j]:
        ax.axvline(analyze.TX[p], color="k", lw=0.8, alpha=0.5)
        ax.set_xscale("log"); ax.set_xlim(1e3, 2e6)
    axes[0, j].set_ylim(0, 1); axes[0, j].set_ylabel("W / Ly"); axes[0, j].set_title(f"p = {p}: strip width (t_x marked)")
    axes[0, j].legend(fontsize=7, loc="upper left", ncol=2)
    axes[1, j].set_ylim(0, 1.1); axes[1, j].set_ylabel("interior / value at t ~ 1e3"); axes[1, j].set_xlabel("t (steps)")
    ax = axes[2, j]
    for Lx, mk in ((2048, "o"), (4096, "s"), (1024, "^")):
        rs = [r for r in rows if r["p"] == p and r["Lx"] == Lx]
        if not rs:
            continue
        qq = sorted(set(r["q"] for r in rs))
        Ly = rs[0]["Ly"]
        for w, fill, lab in (("vE", "none", "early [2e3, 1e5]"), ("vL", None, "late [t_x, T]")):
            m = [analyze.seed_mean([r["v"][("o", w[1:])] for r in rs if r["q"] == q]) for q in qq]
            x = np.array(qq); y = np.array([a[0] for a in m]); e = np.array([a[1] if np.isfinite(a[1]) else 0 for a in m])
            ok = np.isfinite(y)
            ax.errorbar(x[ok], y[ok], e[ok], fmt=mk, mfc=fill, color="C0" if w == "vE" else "C3", ms=6,
                        label=f"Lx {Lx}: {lab}")
        key = (p, Lx, Ly, ("o", "L"))
        if Lx == 2048 and key in res and "qstar" in res[key]:
            r = res[key]
            xx = np.linspace(min(qq), max(qq), 50)
            ax.plot(xx, r["slope"] * (xx - r["qstar"]), "-", color="C3", lw=1,
                    label=f"late line: q* = {r['qstar']:.6f} +- {r['sd']:.6f}")
    ax.axvline(analyze.QS_OLD[p], color="k", ls="--", lw=0.8, label=f"Table I q* = {analyze.QS_OLD[p]:g}")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("q"); ax.set_ylabel("v (rows/step per front)"); ax.set_title(f"p = {p}: v(q)")
    ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(f"{D}/strips.png", dpi=120)
print("strips.png written")
