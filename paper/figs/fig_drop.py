"""Fig. drop: critical droplets at p = 1/2. (a) survival probability of a disordered disk of radius R
in the dead state (E14_drop05, 500-1000 replicas per point); (b) R_c (P_surv = 1/2) against the
inverse planar-front velocity 1/v(q), with dR/dt = v - D/R + offset."""
import csv
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np

from common import COL, FIGDATA, RUNS, label, save
import spread

QS = 0.04597    # strip-front threshold at p = 1/2 (E12_phase)


def surv():
    """P_surv(R) per x from E14_drop05; the R = 0.75 (2x2 block, sites alive with probability 1/2) points are
    replaced by the 4000-replica runs of R2_seed22 (2026-10-31) where available, and the full 2x2 block
    of live sites (R2_seed22 'full') is returned separately."""
    rows, full = {}, {}
    for f in sorted(glob.glob(os.path.join(RUNS, "E14_drop05", "out", "*.txt"))):
        m = re.match(r"p0.5_a0.5_j2_q([\d.]+)_R([\d.]+)_L(\d+)_drop", os.path.basename(f))
        t, data, meta = spread.load1(f)
        S = data.sum(0)
        n, s = S[-1, 0], S[-1, 1]
        rows.setdefault(round(float(m[1]) / QS, 4), []).append((float(m[2]), s / n, np.sqrt(s * (n - s) / n) / n, n))
    for f in sorted(glob.glob(os.path.join(RUNS, "R2_seed22", "out", "*.txt"))):
        m = re.match(r"p0.5_a0.5_j2_q([\d.]+)_R0.75(full|half)_L1024_drop", os.path.basename(f))
        if not m:
            continue
        t, data, meta = spread.load1(f)
        S = data.sum(0)
        n, s = S[-1, 0], S[-1, 1]
        x = round(float(m[1]) / QS, 4)
        pt = (0.75, s / n, np.sqrt(s * (n - s) / n) / n, n)
        if m[2] == "half" and x in rows:
            rows[x] = [r for r in rows[x] if r[0] != 0.75] + [pt]
        elif m[2] == "full":
            full[x] = pt
    return {x: np.array(sorted(v)) for x, v in rows.items()}, full


def rc(d):
    """R where P_surv = 1/2, by linear interpolation of logit P in ln R between the bracketing points."""
    R, P = d[:, 0], np.clip(d[:, 1], 1e-3, 1 - 1e-3)
    k = np.flatnonzero((P[:-1] < 0.5) & (P[1:] >= 0.5))[0]
    lg = np.log(P / (1 - P))
    f = -lg[k] / (lg[k + 1] - lg[k])
    return float(np.exp(np.log(R[k]) + f * (np.log(R[k + 1]) - np.log(R[k]))))


def front_v():
    """Seed-mean strip-front velocity against x = q/q* at p = 1/2 (E4_voter, E12_phase if present)."""
    out = {}
    for fn in ("E4_voter_fronts.csv", "E12_phase_fronts.csv"):
        path = os.path.join(FIGDATA, fn)
        if not os.path.exists(path):
            continue
        for r in csv.DictReader(open(path)):
            m = re.match(r"p0\.5_a0\.5_j2_q([\d.]+)_Lx2048_Ly4096_s\d", r["run"])
            if m and r["v"] not in ("nan", ""):
                out.setdefault(float(m[1]) / QS, []).append(float(r["v"]))
    x = np.array(sorted(out))
    return x, np.array([np.mean(out[k]) for k in x])


if __name__ == "__main__":
    S, FULL = surv()
    xv, vv = front_v()
    fig, (a, b) = plt.subplots(1, 2, figsize=(1.0 * 3.4 * 2, 2.2), gridspec_kw=dict(wspace=0.3))
    cmap = plt.get_cmap("viridis")
    xs = sorted(S)
    res = []
    for i, x in enumerate(xs):
        d = S[x]
        c = cmap(i / (len(xs) - 1))
        a.errorbar(d[:, 0], d[:, 1], d[:, 2], fmt="o-", color=c, ms=3, capsize=1.2, lw=0.8, label=f"${x:.3f}\\,q^*$")
        if x in FULL:
            a.errorbar([FULL[x][0]], [FULL[x][1]], [FULL[x][2]], fmt="*", color=c, ms=5, capsize=1.2, mec="k", mew=0.3)
        R_c = rc(d)
        v = np.interp(x, xv, vv) if x <= xv.max() else np.nan
        res.append((x, R_c, v))
        print(f"x={x:.3f} R_c={R_c:.2f} v={v:.5f} R_c*v={R_c * v:.3f}  n={d[:, 3].min():.0f}-{d[:, 3].max():.0f}")
    a.axhline(0.5, color="0.6", lw=0.6, ls=":")
    a.set(xscale="log", xlabel=r"initial radius $R$", ylabel=r"$P_{\rm surv}$", ylim=(-0.02, 1.02))
    a.legend(fontsize=5.5, loc="upper left", title=r"$q$", title_fontsize=6)
    label(a, "(a)")
    res = np.array(res)
    ok = np.isfinite(res[:, 2])
    inv = 1 / res[ok, 2]
    A = np.vstack([inv, np.ones_like(inv)]).T
    (D, dR), *_ = np.linalg.lstsq(A, res[ok, 1], rcond=None)
    print(f"fit R_c = D/v + dR: D = {D:.3f}, dR = {dR:.2f}")
    b.plot(inv, res[ok, 1], "ko", ms=4, mfc="w")
    for (x, R_c, v) in res[ok]:
        b.annotate(f"{x:.3f}", (1 / v, R_c), textcoords="offset points", xytext=(-8 if 1 / v > 150 else 4, -9), fontsize=6, ha="right" if 1 / v > 150 else "left")
    g = np.linspace(0, inv.max() * 1.05, 50)
    b.plot(g, D * g + dR, "C3-", lw=0.9, label=f"$R_c = D/v + {dR:.1f}$, $D={D:.2f}$")
    b.set(xlabel=r"$1/v(q)$ (steps per site)", ylabel=r"critical radius $R_c$", xlim=(0, None), ylim=(0, None))
    b.legend(fontsize=6.5, loc="upper left")
    label(b, "(b)", x=0.9, y=0.12)
    save(fig, "fig_drop")
