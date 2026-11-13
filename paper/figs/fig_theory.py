"""Fig. theory: voter-point perturbation theory thresholds against walk time T, with simulations.

Panel (b) revised 2026-10-30: shows the >= 5-lineage fraction (the terms that can make f bistable)."""
import os

import matplotlib.pyplot as plt
import numpy as np

from common import COL, THEORY, label, save


def table(path):
    rows = []
    for ln in open(path):
        a = ln.split()
        if len(a) == 8 and a[0].isdigit():
            rows.append([float(x) for x in a])
    return np.array(rows)


d = table(os.path.join(THEORY, "vpt_run_th3.txt"))
d = d[d[:, 0] >= 4]
T, p3, ka, kf, kth = d[:, 0], d[:, 2], d[:, 3], d[:, 4], d[:, 5]

fig, (a, b) = plt.subplots(2, 1, figsize=(COL, 3.9), gridspec_kw=dict(height_ratios=[1.6, 1], hspace=0.08),
                           sharex=True)
a.axhspan(0.245, 0.275, color="0.88", lw=0)
a.plot(T, ka, "o-", color="C3", ms=2.5, label=r"$\kappa_a/\alpha$: $f'(0)=0$")
a.plot(T, kf, "s-", color="C0", ms=2.5, label=r"$\kappa_{\rm front}/\alpha$: $\int_0^{1/2} f = 0$")
a.plot(T, kth, "^-", color="C2", ms=2.5, label=r"three-lineage terms only")
# simulations: c = q*/(alpha (1-p)) at alpha = 1/2 (strip fronts) and the alpha scan at p 0.98, 0.99
sims = {0.9: 0.239, 0.95: 0.248, 0.98: 0.252, 0.99: 0.256, 0.995: 0.256}  # Table I (2026-10-31: p >= 0.99 from strips run to 4 t_x)
for p, c in sims.items():   # plotted at T = 1/(1-p), the time scale on which the perturbation acts
    a.plot(1 / (1 - p), c, "kD", ms=3, mfc="w", mew=0.8)
a.plot([], [], "kD", ms=3, mfc="w", mew=0.8, label=r"simulation at $T=1/(1-p)$")
a.text(3.3e4, 0.243, r"$c=0.26\pm0.015$", fontsize=6.5, ha="center", va="top", color="0.35")
a.set_ylim(0.15, 0.5)
a.set_ylabel(r"$\kappa^*/\alpha$")
a.legend(loc="upper right", fontsize=6.3)
label(a, "(a)")
# (b) Only patterns with five or more lineage classes contribute terms of degree five and
# higher to f(u); four-class patterns give the same cubic as three-class ones (vpt.pattern_poly).
# The bistable window therefore lives in the >= 5-class fraction, read from the pattern
# histogram itself; the >= 4-class fraction (1 - p3, from vpt_run_th3.txt) is shown for reference.
import sys
sys.path.insert(0, THEORY)
from vpt import read_patterns
pats = read_patterns(os.path.join(THEORY, "patterns.txt"))
f5 = []
for t in T:
    trials, pl = pats[int(t)]
    n3 = sum(c for c, ns in pl if len(ns) >= 2)
    n5 = sum(c for c, ns in pl if len(ns) >= 4)
    f5.append(n5 / n3)
b.plot(T, 1 - p3, "o--", color="0.6", ms=2.5, label=r"$\geq 4$ lineages / $\geq 3$")
b.plot(T, f5, "o-", color="k", ms=2.5, label=r"$\geq 5$ lineages / $\geq 3$")
b.set_xscale("log", base=2)
b.set_xlabel(r"walk time $T$ (steps)")
b.set_ylabel(r"fraction of patterns")
b.set_ylim(0, 1)
b.set_xlim(3, 1e5)
b.legend(loc="upper right", fontsize=6.3)
label(b, "(b)", y=0.16)
save(fig, "fig_theory")
