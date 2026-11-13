#!/usr/bin/env python3
"""E7 magnetisation finite-size scaling, pooling E6 and E7 bulk runs near q*.

    python3 binder.py RUNS_DIR [--hist OUT.json]   (RUNS_DIR/{E6_long*,E7_binder*}/out/*.txt)

Per (p, q, L), over t >= 2e5 and only while the run is active (rho > 0.8 x the group's median
rho over the same window, so an absorbing or collapsing run stops counting when it drops):
  chi   = L^2 <m^2>, m = 2 alive - 1 (the paper's convention; until 2026-10-31 this script used
          alive - 1/2, so its chi was 4x smaller); error from the scatter of per-run values / sqrt(runs)
          (independent of tau_m, so it stays honest when a run holds only a few samples);
  U4    = 1 - <m^4> / (3 <m^2>^2) pooled over runs; error by jackknife over runs;
  |m|   = <|m|>; a magnetised active phase has |m| that does not shrink with L;
  tau   = integrated autocorrelation time of m, averaged over runs; n_eff = summed over runs;
  x2    = chi_L / chi_{L/2} and gamma/nu = log2(x2).
Finite xi (first order): x2 -> 1, U4 -> 0.  2D Ising critical: x2 -> 2^1.75 = 3.36 and U4 ->
0.61.  Ordered (Ising below T_c): x2 -> 4 and U4 -> 2/3.
--hist writes P(m L^{1/8}) histograms per group (for the plots page).
"""
import glob
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

QSTAR = {0.9: 0.0119, 0.98: 0.00252, 0.99: 0.00127, 0.995: 0.00062}
T_MIN = 2e5


def load(f):
    return np.array([[float(x) for x in l.split()] for l in open(f) if not l.startswith("#")])


def tau_int(m, dt):
    mc = m - m.mean()
    if len(mc) < 10 or not mc.any():
        return float("nan")
    ac = np.correlate(mc, mc, "full")[len(mc) - 1:]
    ac /= ac[0]
    k0 = np.argmax(ac < 0) or len(ac)
    return (0.5 + ac[1:k0].sum()) * dt


def collect(root, hist=False):
    """Rows (one per (p, q, L), sorted) with per-run chi in row["chi_runs"]; histograms if asked."""
    groups = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(root, "E6_long*", "out", "*.txt")) +
                    glob.glob(os.path.join(root, "E7_binder*", "out", "*.txt"))):
        n = os.path.basename(f)
        mt = re.match(r"^p([\d.]+)_a[\d.]+_j2_q([\d.]+)_L(\d+)_s(\d+)\.txt$", n)
        if not mt:
            continue
        p, q, L = float(mt[1]), float(mt[2]), int(mt[3])
        if p not in QSTAR or q / QSTAR[p] < 1.005:
            continue
        ts = load(f)
        done = any(l.startswith("# done") for l in open(f))
        groups[(p, q, L)].append((ts[ts[:, 0] >= T_MIN], done, os.path.dirname(os.path.dirname(f))))

    rows, hists = [], {}
    for (p, q, L), runs in sorted(groups.items()):
        win = [r for r, _, _ in runs if len(r)]
        if not win:
            continue
        cut = 0.8 * np.median(np.concatenate([r[:, 2] for r in win]))
        per = []
        for r in win:
            act = r[r[:, 2] > cut]
            if len(act) < 20:
                continue
            m = 2 * act[:, 1] - 1  # paper convention m = 2n - 1 (before 2026-10-31: n - 1/2, chi 4x smaller)
            dt = np.median(np.diff(act[:, 0]))
            ta = tau_int(m, dt)
            per.append({"m2": (m ** 2).mean(), "m4": (m ** 4).mean(), "am": np.abs(m).mean(), "m": m,
                        "tau": ta, "neff": (act[-1, 0] - act[0, 0]) / (2 * ta), "w": len(m)})
        if not per:
            continue
        k = len(per)
        m2 = np.array([x["m2"] for x in per])
        m4 = np.array([x["m4"] for x in per])
        w = np.array([x["w"] for x in per], float)
        chi = L * L * np.average(m2, weights=w)
        chi_err = L * L * m2.std(ddof=1) / np.sqrt(k) if k > 1 else float("nan")

        def u4(sel):
            a2, a4 = np.average(m2[sel], weights=w[sel]), np.average(m4[sel], weights=w[sel])
            return 1 - a4 / (3 * a2 ** 2)
        allsel = np.ones(k, bool)
        U = u4(allsel)
        if k > 1:
            jk = np.array([u4(np.arange(k) != i) for i in range(k)])
            U_err = np.sqrt((k - 1) / k * ((jk - jk.mean()) ** 2).sum())
        else:
            U_err = float("nan")
        absorbed = sum(1 for r, d, _ in runs if len(r) and r[-1, 2] < 0.05)
        running = sum(1 for _, d, _ in runs if not d)
        rows.append({"p": p, "x": q / QSTAR[p], "L": L, "runs": len(runs), "used": k, "abs": absorbed,
                     "run": running, "chi": chi, "chi_err": chi_err, "U4": U, "U4_err": U_err,
                     "am": np.average([x["am"] for x in per], weights=w),
                     "tau": np.nanmean([x["tau"] for x in per]), "neff": np.nansum([x["neff"] for x in per]),
                     "chi_runs": (L * L * m2).tolist()})
        if hist:
            y = np.concatenate([x["m"] for x in per]) * L ** 0.125
            h, e = np.histogram(y, bins=41, range=(-1.2, 1.2), density=True)
            hists[f"p{p}_x{q / QSTAR[p]:.3f}_L{L}"] = {"edges": e.round(4).tolist(), "P": h.round(4).tolist()}
    return rows, hists


def main():
    hist_out = sys.argv[sys.argv.index("--hist") + 1] if "--hist" in sys.argv else None
    rows, hists = collect(sys.argv[1], hist=bool(hist_out))
    print(f"{'p':6}{'q/q*':7}{'L':>5}{'runs':>9}{'chi':>9}{'±':>7}{'x2':>6}{'g/nu':>6}{'U4':>7}{'±':>6}"
          f"{'<|m|>':>7}{'tau':>8}{'n_eff':>7}")
    prev = {}
    for r in rows:
        key = (r["p"], round(r["x"], 3))
        x2 = r["chi"] / prev[key][1] if key in prev and prev[key][0] * 2 == r["L"] else float("nan")
        prev[key] = (r["L"], r["chi"])
        note = (f"{r['used']}/{r['runs']}" + (f" a{r['abs']}" if r["abs"] else "") +
                (f" r{r['run']}" if r["run"] else ""))
        print(f"{r['p']:<6}{r['x']:<7.3f}{r['L']:>5}{note:>9}{r['chi']:9.0f}{r['chi_err']:7.0f}{x2:6.2f}"
              f"{np.log2(x2) if x2 == x2 else float('nan'):6.2f}{r['U4']:7.2f}{r['U4_err']:6.2f}{r['am']:7.3f}"
              f"{r['tau']:8.0f}{r['neff']:7.0f}")
    print("runs = used/total (a = absorbed, r = still running); x2 = chi_L / chi_{L/2}; "
          "Ising: x2 3.36, g/nu 1.75, U4 0.61")
    if hist_out:
        json.dump(hists, open(hist_out, "w"), separators=(",", ":"))


if __name__ == "__main__":
    main()
