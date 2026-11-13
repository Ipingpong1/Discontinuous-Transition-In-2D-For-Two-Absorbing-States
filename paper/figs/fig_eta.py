"""Fig. eta: the single-noise rule w(j) = alpha R(j) + eta [j in {1,2,4,5}] at alpha = 1/2 (E17).
(a) strip-front velocity against eta (3 seeds); (b) bulk interface density against time, L = 2048."""
import glob
import os

import matplotlib.pyplot as plt
import numpy as np

from common import COL, RUNS, label, save
import fronts_util as fu

import csv
from collections import defaultdict
from common import FIGDATA
g = defaultdict(list)
for r in csv.DictReader(open(os.path.join(FIGDATA, "E17_eta_front_fronts.csv"))):
    if int(r["npts"]) >= 20:   # fronts fitted over >= 20 records; below eta* the interior orders first
        g[float(r["run"].split("_q")[1].split("_")[0])].append(float(r["v"]))
q = np.array(sorted(g))
v = np.array([np.mean(g[k]) for k in q])
se = np.array([np.std(g[k], ddof=1) / np.sqrt(len(g[k])) if len(g[k]) > 1 else np.nan for k in q])
fig, (a, b) = plt.subplots(2, 1, figsize=(COL, 3.6), gridspec_kw=dict(hspace=0.35))
ok = np.isfinite(v)
a.errorbar(q[ok] * 1e3, v[ok] * 1e3, se[ok] * 1e3, fmt="o", color="C0", ms=3.5, capsize=1.2)
a.axhline(0, color="0.6", lw=0.5)
qs, sl = fu.crossing(q, v)
print(f"eta* = {qs:.5f}, slope {sl:.2f}, eta* dv/deta = {qs * sl:.3f}")
a.axvline(qs * 1e3, color="0.6", lw=0.5, ls=":")
a.set(xlabel=r"$\eta\ (10^{-3})$", ylabel=r"front velocity $v\ (10^{-3})$")
label(a, "(a)")
cmap = plt.get_cmap("coolwarm")
files = sorted(glob.glob(os.path.join(RUNS, "E17_eta_bulk", "out", "eta_a0.5_q*_L2048_s1.txt")),
               key=lambda f: float(os.path.basename(f).split("_q")[1].split("_")[0]))
for i, f in enumerate(files):
    eta = float(os.path.basename(f).split("_q")[1].split("_")[0])
    r = np.array([[float(x) for x in ln.split()] for ln in open(f) if ln.strip() and not ln.startswith("#")])
    sel = r[:, 0] >= 10
    b.plot(r[sel, 0], r[sel, 2], color=cmap(i / (len(files) - 1)), lw=0.8, label=f"{eta * 1e3:.4g}")
b.set(xscale="log", xlabel=r"$t$", ylabel=r"interface density $\rho$", ylim=(0, 0.42))
b.legend(fontsize=5.5, title=r"$\eta\ (10^{-3})$", title_fontsize=6, ncol=2, loc="lower left")
label(b, "(b)", x=0.88, y=0.95)
save(fig, "fig_eta")
