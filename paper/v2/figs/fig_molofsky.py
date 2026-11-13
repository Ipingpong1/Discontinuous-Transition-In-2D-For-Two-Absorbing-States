"""Fig. molofsky: the five-cell rule of Molofsky et al. (Appendix D): the lineage closure and
the direct simulation of the rule.

(a), (b) read campaign/theory/molofsky_1e7.csv, written by
    python3 vpt.py molofsky --load patterns5_1e7.txt --boot 300 --csv molofsky_1e7.csv
(see campaign/theory/MOLOFSKY.md).  (c) reads the seed-pooled planar-front tables
campaign/molofsky_sim/scan/table_p2_*.txt (Lx 1024) and scan_wide/table_p2_*.txt (Lx 2048),
written by `python3 mfront.py table scan/p2_X_p1_*_s*.txt` (see campaign/molofsky_sim/NOTES.md);
its inset reads the fine scan at 1e-4 steps, scan_fine/table.txt (Lx 2048) and
scan_fine/wide/table.txt (Lx 4096), written by `python3 sec8_tables.py` (NOTES.md Sec. 8)."""
import os

import matplotlib.pyplot as plt
import numpy as np

from common import COL, ROOT, THEORY, label, save

SIM = os.path.join(ROOT, "campaign", "molofsky_sim")

# simulation (NOTES.md Secs. 3-6): front reversal on the two cuts, and the band of p_1 above it in
# which the mean population of a 2x2 seed decreases (absorbing state linearly stable)
TH_FRONT = (106.33, 106.40)        # p_2 = 0.45, 0.50
TH_A = (104.5, 105.0)              # seed decay stops at p_1 = 0.187-0.188 and 0.173-0.174
P1STAR = {0.45: 0.18535, 0.50: 0.17057}   # front reversal on the two cuts (NOTES.md Sec. 8.1)


def theta(p1, p2):
    return np.degrees(np.arctan2(p2 - 0.4, p1 - 0.2))


def read_table(path):
    rows = []
    for line in open(path):
        w = line.split()
        if not w or w[0] == "p2":
            continue
        rows.append((float(w[0]), float(w[1]), float(w[4]), float(w[5]), float(w[-1])))
    return np.array(rows)  # p2, p1, v_adapt, err, interior density at the end / its t = 1000 value


d = np.genfromtxt(os.path.join(THEORY, "molofsky_1e7.csv"), delimiter=",", names=True)
d = d[d["T"] >= 1]
T = d["T"]

fig = plt.figure(figsize=(COL, 5.6))
gs = fig.add_gridspec(2, 1, height_ratios=[2.5, 1.55], hspace=0.32)
gs0 = gs[0].subgridspec(2, 1, height_ratios=[1.5, 1], hspace=0.08)
a = fig.add_subplot(gs0[0])
b = fig.add_subplot(gs0[1], sharex=a)
c = fig.add_subplot(gs[1])

# (a) the closure's rays against T, with the simulation's line and seed-decay edge
h_band = a.axhspan(min(TH_A), max(TH_A), color="C0", alpha=0.15, lw=0, label="seeds stop decaying (sim.)")
h_rev = a.axhline(np.mean(TH_FRONT), color="k", lw=1.0, label="front reversal (sim.)")
h_dyad = a.axhline(108.43, color="0.5", ls="--", lw=0.8, label=r"dyad rule $3p_1+p_2=1$")
h_a, = a.plot(T, d["th_a"], "o-", color="C0", ms=2.5, label=r"$f'(0)=0$")
h_f, = a.plot(T, d["th_f"], "s-", color="C3", ms=2.5, label=r"$\int_0^{1/2} f = 0$")
h_h, = a.plot(T, d["th_h"], "^-", color="C2", ms=2.5, label=r"$f'(1/2)=0$")
h_3, = a.plot(T, d["th_3"], "-.", color="k", lw=0.8, label="three lineages only")
a.set_ylabel(r"direction ($^\circ$)")
a.set_ylim(95, 125)
leg1 = a.legend(handles=[h_a, h_f, h_h, h_3], loc="upper right", fontsize=6, ncol=2, columnspacing=0.8,
                handlelength=1.5, labelspacing=0.3)
a.add_artist(leg1)
a.legend(handles=[h_rev, h_band, h_dyad], loc="lower right", fontsize=6, ncol=1, handlelength=1.5,
         labelspacing=0.3)
a.tick_params(labelbottom=False)
label(a, "(a)")

# (b) sector width
b.plot(T, d["width"], "o-", color="k", ms=2.5)
b.set_yscale("log")
b.set_xscale("log", base=2)
b.set_xlabel(r"walk time $T$ (steps)")
b.set_ylabel(r"sector width ($^\circ$)")
b.set_xlim(0.7, 1e5)
label(b, "(b)")

# (c) simulated front velocity against direction, on the two cuts
c.axhline(0, color="k", lw=0.5)
c.axvspan(np.mean(TH_A), np.mean(TH_FRONT), color="0.9", lw=0, zorder=0)
for Tc, ls in ((16, ":"), (1024, "--")):
    r = d[d["T"] == Tc][0]
    for key, col in (("th_a", "C0"), ("th_f", "C3"), ("th_h", "C2")):
        c.axvline(r[key], color=col, ls=ls, lw=0.8, zorder=1)
for p2, mk, col, lab in ((0.45, "o", "k", r"$p_2 = 0.45$"), (0.50, "s", "0.45", r"$p_2 = 0.50$")):
    for sub, m, lw in (("scan", mk, 0), ("scan_wide", "^", 0)):
        t = read_table(os.path.join(SIM, sub, f"table_p2_{p2:.2f}.txt"))
        th, v, e, cen = theta(t[:, 1], t[:, 0]), t[:, 2], t[:, 3], t[:, 4]
        ok = np.isfinite(v)
        e = np.where(np.isfinite(e), e, 0)
        for stat, mfc in ((cen > 0.8, col), (cen <= 0.8, "w")):
            s = ok & stat
            c.errorbar(th[s], v[s], yerr=e[s], fmt=m, color=col, mfc=mfc, ms=3.2, mew=0.8, capsize=0,
                       elinewidth=0.6, ls="none", label=lab if (sub == "scan" and mfc == col) else None)
c.plot([], [], "^", color="0.3", mfc="0.3", ms=3.2, label=r"$L_x = 2048$")
c.plot([], [], "o", color="0.3", mfc="w", ms=3.2, label="interior decays")
c.set_xlim(103.3, 109.0)
c.set_ylim(-0.012, 0.047)
c.set_xlabel(r"direction $\theta$ ($^\circ$)")
c.set_ylabel(r"front velocity $v$")
c.legend(loc="lower left", fontsize=6, handlelength=1.2, labelspacing=0.3)
c.text(105.4, 0.036, "seeds decay,\nactive phase\ninvades", ha="center", va="top", fontsize=6)
c.text(107.9, 0.012, "active phase\ncoarsens", ha="center", va="top", fontsize=6)
label(c, "(c)")


# inset of (c): the fine scan across the reversal (NOTES.md Sec. 8.1), v against p_1 - p_1^*
def read_fine(path):
    """Rows of scan_fine/table.txt: p_1, v (records with t >= 2e3) and its error, L_x, and the
    drop of the strip's interior interface density between 3e4 and 1e5 steps (fractions of its
    t = 1000 value; the interior counts as stationary when the drop is at most 0.03)."""
    rows = []
    for line in open(path):
        w = [s.strip() for s in line.strip().strip("|").split("|")]
        if len(w) < 8 or not w[0][:1].isdigit():
            continue
        val, err = w[1].split("(")
        e = int(err.rstrip(")")) * 10.0 ** -len(val.split(".")[1])
        frac = [float(s) for s in w[6].split("(")[1].rstrip(")").split("/")]
        rows.append((float(w[0]), float(val), e, int(w[5]), frac[1] - frac[3]))
    return rows


ci = c.inset_axes([0.6, 0.57, 0.39, 0.41])
ci.axhline(0, color="k", lw=0.5)
ci.axvline(0, color="k", lw=0.5)
for sub, m in (("", None), ("wide", "^")):
    for p1, v, e, lx, drop in read_fine(os.path.join(SIM, "scan_fine", sub, "table.txt")):
        p2 = 0.45 if p1 > 0.18 else 0.50
        col, dx = ("k", -0.05) if p2 == 0.45 else ("0.45", 0.05)
        mk = m or ("o" if p2 == 0.45 else "s")
        ci.errorbar((p1 - P1STAR[p2]) * 1e4 + dx, v * 1e3, yerr=e * 1e3, fmt=mk, color=col,
                    mfc=col if drop <= 0.03 else "w", ms=2.6, mew=0.7, capsize=0, elinewidth=0.6,
                    ls="none")
ci.set_xlim(-3.5, 2.5)
ci.set_ylim(-2.6, 3.6)
ci.set_xticks([-3, -2, -1, 0, 1, 2])
ci.set_yticks([-2, 0, 2])
ci.tick_params(labelsize=5.5, pad=1.5, length=2)
ci.set_xlabel(r"$\Delta p_1\times10^{4}$", fontsize=6, labelpad=1)
ci.set_ylabel(r"$v\times10^{3}$", fontsize=6, labelpad=1)
save(fig, "fig_molofsky")
