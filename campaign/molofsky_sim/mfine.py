#!/usr/bin/env python3
"""Analysis for the fine scan (NOTES.md Sec. 8).

    python3 mfine.py strips scan_fine/*.txt      seed-pooled v(p1) and interior history per (p2, p1)
    python3 mfine.py runs   scan_fine/*.txt      per run: v over two windows, interior density vs t
    python3 mfine.py bulk   bulk_fine/p2_*.txt   rho_bond(t) at fixed times, late drift, absorption
    python3 mfine.py life   bulk_fine/life_*.txt lifetimes per (p2, p1, L) and plateau diagnostics
    python3 mfine.py fit                         rho(delta p1) fits (log vs power) from bulk data

Velocities as in mfront.py (adaptive threshold = half the current interior Moore-interface
density, records with t >= tmin, W in [0.1, 0.9] Ly, interior above 60 % of its t = 1000 value);
'late' uses tmin = 3e4.  cen(t) = interior interface density / its t = 1000 value.
"""
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mfront import track, fit_v, KEY  # noqa: E402

BKEY = re.compile(r"p2_(?P<p2>[0-9.]+)_p1_(?P<p1>[0-9.]+)_L(?P<L>[0-9]+)(?P<tag>_T1e6)?_s(?P<s>[0-9]+)")
CEN_T = (10000, 30000, 60000, 100000)


def header(path):
    h = {}
    for line in open(path):
        if line.startswith("# done"):
            w = line.split()
            h["t_end"], h["absorbed"] = int(w[3]), int(w[5])
        elif line.startswith("# Lx"):
            w = line.split()
            h["Lx"], h["Ly"] = int(w[2]), int(w[4])
    return h


def series(path):
    a = np.array([[float(x) for x in l.split()] for l in open(path) if not l.startswith("#")])
    return a[:, 0], a[:, 2], a[:, 1]  # t, rho_bond, alive


def at(t, y, tt):
    return float(np.interp(tt, t, y)) if tt <= t[-1] else np.nan


def cmd_runs(files):
    print("run                              v_adapt   v_late    v_int    npts  W_end/Ly rho_b   cen@1e4 3e4   6e4   1e5   abs t_end")
    for f in files:
        r = track(f)
        v, va, vA, n, tm = fit_v(r)
        v2, va2, vA2, n2, tm2 = fit_v(r, tmin=30000)
        h = header(f)
        cens = " ".join(f"{at(r['t'], r['cen'], tt):5.2f}" for tt in CEN_T)
        print(f"{os.path.basename(f)[:-4]:32s} {va:+.5f}  {va2:+.5f}  {vA:+.5f}  {n:3d}  {r['W'][-1] / r['Ly']:.3f}   {r['rho_b']:.3f}  {cens}  "
              f"{h.get('absorbed', -1)} {h.get('t_end', 0)}")


def pooled(vals):
    a = np.array(vals, float); a = a[np.isfinite(a)]
    if len(a) == 0:
        return "   -    ", np.nan, np.nan
    m = a.mean(); e = a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else np.nan
    return f"{m:+.4f}" + (f"({e * 1e4:.0f})" if np.isfinite(e) else "    "), m, e


def cmd_strips(files):
    g = defaultdict(list)
    for f in files:
        m = KEY.search(os.path.basename(f))
        if not m:
            continue
        r = track(f)
        v, va, vA, n, tm = fit_v(r)
        v2, va2, vA2, n2, tm2 = fit_v(r, tmin=30000)
        h = header(f)
        g[(float(m["p2"]), float(m["p1"]))].append(dict(va=va, va2=va2, vA=vA, rho_b=r["rho_b"], Lx=r["Lx"],
                                                      cen=[at(r["t"], r["cen"], tt) for tt in CEN_T],
                                                      Wend=r["W"][-1] / r["Ly"], absorbed=h.get("absorbed", -1), t_end=h.get("t_end", 0)))
    print("| p_1 | v (t >= 2e3) | v (t >= 3e4) | v_int | seeds | Lx | interior at 1e4 / 3e4 / 6e4 / 1e5 (fraction of t = 1000) | W_end/Ly |")
    print("|---|---|---|---|---|---|---|---|")
    for key in sorted(g):
        rs = g[key]
        s1, m1, e1 = pooled([x["va"] for x in rs]); s2, m2, e2 = pooled([x["va2"] for x in rs]); s3, m3, e3 = pooled([x["vA"] for x in rs])
        cen = np.nanmean([x["cen"] for x in rs], axis=0)
        last = cen[np.isfinite(cen)][-1]
        st = "stationary" if last > 0.8 else "decaying" if last > 0.6 else "gone"
        print(f"| {key[1]:.5f} | {s1} | {s2} | {s3} | {len(rs)} | {rs[0]['Lx']} | {st} ({' / '.join(f'{c:.2f}' for c in cen)}) | "
              f"{np.mean([x['Wend'] for x in rs]):.3f} |")


def cmd_bulk(files):
    times = (100, 1000, 10000, 30000, 100000, 200000, 300000, 500000, 1000000)
    print("| run | " + " | ".join(f"{t:.0e}" for t in times) + " | drift (last third, per 1e5) | absorbed at |")
    print("|---|" + "---|" * (len(times) + 2))
    for f in files:
        t, rho, alive = series(f)
        h = header(f)
        T = t[-1]
        sel = t >= 2 * T / 3
        if sel.sum() < 3:
            sel = np.arange(len(t)) >= len(t) - 4
        drift = np.polyfit(t[sel], rho[sel], 1)[0] * 1e5 if sel.sum() >= 3 else np.nan
        cells = " | ".join(f"{at(t, rho, tt):.4f}" if tt <= T else "" for tt in times)
        print(f"| {os.path.basename(f)[:-4]} | {cells} | {drift:+.4f} | {h.get('t_end', 0) if h.get('absorbed') else '-'} |")


def cmd_bulk2(files):
    """Seed-pooled bulk table (one row per p2, p1, L): rho_bond at fixed times, ratio rho(3e5)/rho(1e5),
    absorption times."""
    times = (100, 1000, 3000, 10000, 30000, 100000, 200000, 300000)
    g = defaultdict(list)
    for f in files:
        m = BKEY.search(os.path.basename(f))
        if not m or m["tag"]:
            continue
        t, rho, alive = series(f)
        h = header(f)
        g[(float(m["p2"]), float(m["p1"]), int(m["L"]))].append((t, rho, h))
    print("| p_2 | p_1 | L | seeds | " + " | ".join(f"{t:.0e}" for t in times) + " | ratio 3e5/1e5 | absorbed at |")
    print("|---|---|---|---|" + "---|" * (len(times) + 2))
    for key in sorted(g):
        rs = g[key]
        cells = []
        for tt in times:
            vals = [at(t, rho, tt) for t, rho, h in rs]
            cells.append(f"{np.mean(vals):.4f}" if np.all(np.isfinite(vals)) else "-")
        r1, r3 = [at(t, rho, 1e5) for t, rho, h in rs], [at(t, rho, 3e5) for t, rho, h in rs]
        ratio = f"{np.mean(r3) / np.mean(r1):.2f}" if np.all(np.isfinite(r1 + r3)) and np.mean(r1) > 0 else "-"
        absd = [f"{h['t_end']:.1e}" for t, rho, h in rs if h.get("absorbed")]
        print(f"| {key[0]:.2f} | {key[1]:.5f} | {key[2]} | {len(rs)} | " + " | ".join(cells) + f" | {ratio} | {', '.join(absd) if absd else '-'} |")


def cmd_life(files):
    g = defaultdict(list)
    for f in files:
        m = BKEY.search(os.path.basename(f))
        t, rho, alive = series(f)
        h = header(f)
        te = h.get("t_end", t[-1])
        r1e4 = at(t, rho, 1e4)
        # plateau diagnostics: rho at 1e4, at te/2, and the fraction of the lifetime with rho > 0.8 rho(1e4)
        above = t[rho > 0.8 * r1e4]
        g[(float(m["p2"]), float(m["p1"]), int(m["L"]))].append(dict(te=te, absorbed=h.get("absorbed", 0), r1e4=r1e4, rhalf=at(t, rho, te / 2),
                                                                      frac=(above.max() / te) if len(above) else 0.0, done="t_end" in h))
    print("| p_2 | p_1 | L | lifetime per seed | mean | rho at 1e4 | rho at t_end/2 | fraction of lifetime with rho > 0.8 rho(1e4) |")
    print("|---|---|---|---|---|---|---|---|")
    for key in sorted(g):
        rs = g[key]
        tes = [f"{x['te']:.2e}" + ("" if x["absorbed"] else "+") for x in rs]
        fin = [x["te"] for x in rs if x["absorbed"]]
        mean = f"{np.mean(fin):.2e}" + (f" ({len(fin)}/{len(rs)} absorbed)" if len(fin) < len(rs) else "") if fin else f"> {min(x['te'] for x in rs):.0e}"
        print(f"| {key[0]:.2f} | {key[1]:.5f} | {key[2]} | {' '.join(tes)} | {mean} | {np.mean([x['r1e4'] for x in rs]):.3f} | "
              f"{np.mean([x['rhalf'] for x in rs]):.3f} | {np.mean([x['frac'] for x in rs]):.2f} |")


def late_mean(path, frac=1 / 3):
    t, rho, alive = series(path)
    T = t[-1]
    sel = t >= (1 - frac) * T
    if sel.sum() < 3:
        sel = np.arange(len(t)) >= len(t) - 4
    drift = np.polyfit(t[sel], rho[sel], 1)[0] * 1e5 if sel.sum() >= 3 else np.nan
    return rho[sel].mean(), drift, T, header(path).get("absorbed", 0)


def fit_forms(dp, rho, w=None):
    """1/rho = a + b ln dp  <=>  rho = A / ln(B / dp), A = -1/b, B = exp(-a / b); and
    rho = rho0 + C dp^beta by a grid over beta with linear least squares.  Returns dict."""
    w = np.ones_like(dp) if w is None else w
    X = np.vstack([np.ones_like(dp), np.log(dp)]).T
    coef, *_ = np.linalg.lstsq(X * w[:, None], (1 / rho) * w, rcond=None)
    a, b = coef
    pred = 1 / (a + b * np.log(dp))
    chi_log = float((((rho - pred) * w) ** 2).sum())
    best = None
    for beta in np.linspace(0.05, 1.5, 291):
        Xb = np.vstack([np.ones_like(dp), dp ** beta]).T
        c, *_ = np.linalg.lstsq(Xb * w[:, None], rho * w, rcond=None)
        chi = float((((rho - Xb @ c) * w) ** 2).sum())
        if best is None or chi < best[0]:
            best = (chi, beta, c[0], c[1])
    return dict(A=-1 / b, B=np.exp(-a / b), chi_log=chi_log, chi_pow=best[0], beta=best[1], rho0=best[2], C=best[3])


P1STAR = 0.18535
DP_MIN = 0.5e-4   # points closer to the line than this are not stationary within 3e5 (L 2048) / 1e6 (L 1024) steps


def collect(pattern_sets):
    pts = {}
    for pat in pattern_sets:
        for f in sorted(glob.glob(pat)):
            m = re.search(r"p1_(?P<p1>[0-9.]+)_L(?P<L>[0-9]+)|p1_(?P<p1b>[0-9.]+)_s", os.path.basename(f))
            if "t_end" not in header(f):
                continue  # unfinished run
            p1 = float(m["p1"] or m["p1b"]); L = int(m["L"]) if m["L"] else 2048
            mean, drift, T, ab = late_mean(f)
            pts.setdefault((p1, L), []).append((mean, drift, T, ab, os.path.basename(f)))
    return pts


def fit_report(dp, rho, label):
    r = fit_forms(dp, rho)
    n = len(dp)
    print(f"  {label}: n = {n}")
    print(f"    log   rho = A/ln(B/dp):        A = {r['A']:.4f}  B = {r['B']:.3e}   rms = {np.sqrt(r['chi_log'] / n):.5f}")
    print(f"    power rho = rho0 + C dp^beta:  rho0 = {r['rho0']:.4f}  C = {r['C']:.3f}  beta = {r['beta']:.2f}   rms = {np.sqrt(r['chi_pow'] / n):.5f}")
    return r


def cmd_fit(_):
    print("# stationary rho_bond above the line at p2 = 0.45 (mean over the last third of each run, drift per 1e5 steps)")
    pts = collect(["bulk_fine/p2_0.45_p1_*_L2048_s*.txt", "bulk_fine/p2_0.45_p1_*_L1024_T1e6_s*.txt",
                   "bulk_long/p2_0.45_p1_*_L1024.txt", "bulk/p2_0.45_p1_*_s1.txt"])
    print("| p_1 | dp_1 (p_1* = 0.18535) | L | T | rho (last third) | drift per 1e5 | file |")
    print("|---|---|---|---|---|---|---|")
    for key in sorted(pts):
        for mean, drift, T, ab, f in pts[key]:
            print(f"| {key[0]:.5f} | {key[0] - P1STAR:+.5f} | {key[1]} | {T:.0e} | {mean:.4f} | {drift:+.4f} | {f} |")
    # one value per p1: the longest L 2048 run (fallback: L 1024, longest), averaged over seeds
    best = {}
    for (p1, L), rs in pts.items():
        T = max(r[2] for r in rs)
        cand = (L == 2048, T)
        if p1 not in best or cand > best[p1][0]:
            best[p1] = (cand, np.mean([r[0] for r in rs if r[2] == T]), L, T)
    excl = sorted(p for p in best if P1STAR < p < P1STAR + DP_MIN - 1e-9)
    if excl:
        print("\n# excluded as not stationary (dp_1 < 0.5e-4): " + ", ".join(f"({p - P1STAR:.5f}, {best[p][1]:.4f} still decaying)" for p in excl))
    p1s = np.array(sorted(p for p in best if p >= P1STAR + DP_MIN - 1e-9))
    rho = np.array([best[p][1] for p in p1s])
    print("# points used (one per p_1: longest L 2048 run, seed-averaged; L 1024 if no L 2048 run): "
          + ", ".join(f"({p - P1STAR:.5f}, {r:.4f}, L {best[p][2]})" for p, r in zip(p1s, rho)))
    print("# fits of rho(dp_1), dp_1 = p_1 - p_1*, p_1* = 0.18535 fixed")
    for dmax in (1.0e-3, 5e-3):
        sel = (p1s - P1STAR) <= dmax + 1e-9
        if sel.sum() >= 4:
            fit_report(p1s[sel] - P1STAR, rho[sel], f"dp_1 <= {dmax:.0e}")
    print("\n# p_1* free (grid 0.18515 .. 0.185397), points with p_1 - 0.18535 <= 1e-3")
    sel = (p1s - P1STAR) <= 1.0e-3 + 1e-9
    for form in ("log", "power"):
        bestc = None
        for ps in np.arange(0.18515, 0.185399, 0.0000025):
            dp = p1s[sel] - ps
            if (dp <= 0).any():
                continue
            r = fit_forms(dp, rho[sel])
            chi = r["chi_log"] if form == "log" else r["chi_pow"]
            if bestc is None or chi < bestc[0]:
                bestc = (chi, ps, r)
        chi, ps, r = bestc
        extra = f"A = {r['A']:.4f} B = {r['B']:.3e}" if form == "log" else f"rho0 = {r['rho0']:.4f} C = {r['C']:.3f} beta = {r['beta']:.2f}"
        print(f"  {form:5s}: best p_1* = {ps:.5f}  rms = {np.sqrt(chi / sel.sum()):.5f}  {extra}")


if __name__ == "__main__":
    cmd, files = sys.argv[1], sys.argv[2:]
    files = [f for pat in files for f in sorted(glob.glob(pat))] or files
    files = [f for f in files if re.match(r"(life_)?p2_[0-9.]+_p1_", os.path.basename(f))]
    dict(strips=cmd_strips, runs=cmd_runs, bulk=cmd_bulk, bulk2=cmd_bulk2, life=cmd_life, fit=cmd_fit)[cmd](files)
