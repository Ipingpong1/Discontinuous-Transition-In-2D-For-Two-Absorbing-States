"""Fig. bulk: p = 1/2 bulk runs (E13_bulk05). (a) rho(t) at L = 2048; (b) stationary interface density
against q/q* for L = 256-4096; (c) lifetime of the active state against L: half-life below q* (filled),
mean time to absorption of small systems above q* (open)."""
import matplotlib.pyplot as plt
import numpy as np

from bulk05 import collect, summary
from common import PAGE, label, save

g = collect()
s = summary(g)
fig, axs = plt.subplots(1, 3, figsize=(PAGE, 2.25), gridspec_kw=dict(wspace=0.34))
xs = sorted({k[0] for k in g})
cmap = plt.get_cmap("coolwarm")
col = {x: cmap(i / (len(xs) - 1)) for i, x in enumerate(xs)}
ax = axs[0]
for x in xs:
    rs = g.get((x, 2048), [])
    if not rs:
        continue
    r = rs[0]
    ok = r["t"] >= 10
    ax.plot(r["t"][ok], r["rho"][ok], color=col[x], lw=0.8, label=f"${x:.3f}$")
ax.set(xscale="log", xlabel=r"$t$", ylabel=r"interface density $\rho$", ylim=(0, 0.42), xlim=(10, 3e5))
ax.legend(fontsize=5, ncol=2, title=r"$q/q^*$ ($L=2048$)", title_fontsize=5.5, loc="lower left",
          handlelength=1.2, columnspacing=0.8)
label(ax, "(a)", x=0.86, y=0.95)

ax = axs[1]
Ls = sorted({k[1] for k in g})
mk = {256: "v", 512: "D", 1024: "o", 2048: "s", 4096: "^"}
for L in Ls:
    if L < 256:
        continue
    pts = [(x, s[(x, L)]["plateau"], s[(x, L)]["plateau_sem"]) for x in xs if (x, L) in s]
    pts = np.array([p for p in pts if np.isfinite(p[1])])
    if len(pts):
        ax.errorbar(pts[:, 0], pts[:, 1], pts[:, 2], fmt=mk[L], ms=3, mfc="none", capsize=1, label=f"$L={L}$")
dec = [x for x in xs if x < 1]
ax.plot(dec, np.zeros(len(dec)), "k_", ms=6)
ax.axvline(1, color="0.6", lw=0.5, ls=":")
ax.set(xlabel=r"$q/q^*$", ylabel=r"stationary $\rho$", ylim=(-0.01, 0.37), xlim=(0.885, 1.11))
ax.legend(fontsize=5.5, loc="lower right")
label(ax, "(b)")

ax = axs[2]
for x in xs:
    pts, pts_abs = [], []
    for L in Ls:
        d = s.get((x, L))
        if d is None:
            continue
        if x < 1 and d["n_decayed"] == d["n"]:
            pts.append((L, d["tau_mean"], d["tau_sem"]))

    if pts:
        pts = np.array(pts)
        ax.errorbar(pts[:, 0], pts[:, 1], pts[:, 2], fmt="o-", color=col[x], ms=3, capsize=1, lw=0.8,
                    label=f"${x:.3f}$")
    if pts_abs:
        pts_abs = np.array(pts_abs)
        ax.errorbar(pts_abs[:, 0], pts_abs[:, 1], pts_abs[:, 2], fmt="s--", color=col[x], mfc="w", ms=3,
                    capsize=1, lw=0.8, label=f"${x:.3f}$")
ax.set(xscale="log", yscale="log", xlabel=r"$L$", ylabel=r"half-life $\tau$ of the active state", xlim=(45, 6000))
ax.set_xticks([64, 256, 1024, 4096], ["64", "256", "1024", "4096"])
ax.text(0.97, 0.05, "colors as in (a)", transform=ax.transAxes, ha="right", fontsize=6)
label(ax, "(c)", x=0.04, y=0.95)
save(fig, "fig_bulk05")
