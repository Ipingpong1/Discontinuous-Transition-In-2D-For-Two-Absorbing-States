#!/usr/bin/env python3
"""Analysis for NOTES.md (strip crossings at p = 0.99, 0.995 beyond the crossover time).

    python3 analyze.py strips     # per-run table: v in windows, interior history, status
    python3 analyze.py cross      # seed-mean v(q) per window, crossings with bootstrap errors
    python3 analyze.py bulk       # bulk rho(t) tables, bulk/out together with E6_long*
    python3 analyze.py all        # all three, written to analysis.md and printed

Front tracking from the profiles (per-row interface-site density d_y, strip centre rolled to Ly/2,
32-row running mean), three width estimators per record:
  c  "centre": walk outward from the centre until the smoothed density drops below rho_b/2, with
     rho_b the interior value at t ~ 1e3 (campaign/fronts.py, the Table I estimator; it jumps when a
     hole opens in a thinning interior);
  o  "outer": the outermost row on each side whose smoothed density exceeds 0.3 rho_b (the E8b wall
     width definition of campaign/wall.py; insensitive to interior holes) -- the primary estimator;
  a  "adaptive": as c, but the threshold is half of the current central density (holes matter less);
  A  "integrated": sum_y d_y / rho_b (profile.py; biased low when the interior thins);
  An "integrated, normalized": sum_y d_y / (rho_b cen), i.e. A divided by the current interior
     density (unbiased by a uniformly thinning interior).
W = width between the two edges; v = slope of a linear fit of W(t)/2 per run in fixed windows:
E = [2e3, 1e5] (the Table I estimator's range; for c with its condition cen > 0.7), M = [1e5, t_x],
L = [t_x, T] and its halves L1, L2, with t_x = 3e5 (p 0.99), 4e5 (p 0.995), over records with W in
[0.15, 0.85] Ly (for every estimator, so a strip that has filled the box or dissolved drops out); a
run's v in the M, L, L1, L2 windows counts only if the window holds at least NMIN = 20 records.
cen = central 256 rows / rho_b.
Status from the interior and the width: absorbed (run ended absorbed), gone (last-quarter cen <= 0.5
or W_o < 0.15 Ly), filled (W_o >= 0.85 Ly), decaying (cen fell by more than 0.03 from t_x to the
last quarter), else stationary; "*" = run not finished.
"""
import glob
import os
import sys
from collections import defaultdict

import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(D)
sys.path.insert(0, CAMP)
sys.path.insert(0, os.path.dirname(CAMP))
import fronts  # noqa: E402
from profile import read  # noqa: E402

TX = {0.995: 4e5, 0.99: 3e5}
QS_OLD = {0.995: 0.00062, 0.99: 0.00128}
INTERIOR_T = (1e4, 1e5, 3e5, 6e5, 1e6, 1.3e6, 1.6e6)
EST = ("c", "a", "o", "A", "An")
NMIN = 20  # records (4e4 steps) a run needs in a late window for its v to count
WIN = ("E", "M", "L", "L1", "L2")
ROOT = os.path.dirname(CAMP)
E6 = {0.995: [f"{ROOT}/data/runs/E6_long_p995", f"{ROOT}/data/runs/E6_long_L4096_p995"],
      0.99: [f"{ROOT}/data/runs/E6_long", f"{ROOT}/data/runs/E6_long_L4096"]}


def done_line(txt):
    for line in open(txt):
        if line.startswith("# done"):
            return line[2:].strip()
    return None


def strip_frac(txt):
    """Initial strip fraction of Ly from the header line '# init strip:p0:frac ...'."""
    for line in open(txt):
        if line.startswith("# init"):
            spec = line.split()[2]
            return float(spec.split(":")[2]) if spec.startswith("strip") else 0.5
        if not line.startswith("#"):
            break
    return 0.5


def track(txt):
    """t, widths {c, o, A}, cen, Lx, Ly, rho_b."""
    Lx, Ly, t, alive, iface = read(txt, txt[:-4] + ".prof")
    fr = strip_frac(txt)
    d = np.roll(iface / Lx, Ly // 2 - int(fr * Ly / 2), axis=1)  # strip centre -> Ly/2
    k = np.ones(32) / 32
    ds = np.array([np.convolve(np.concatenate([r[-16:], r, r[:15]]), k, "valid") for r in d])
    i0 = min(np.searchsorted(t, 1000), len(t) - 1)
    c = Ly // 2
    h = int(fr * Ly / 4)  # interior window: central frac*Ly/2 rows (fronts.py: Ly/4)
    rho_b = d[i0, c - h: c + h].mean()
    Wc, Wa, Wo, WA, WAn, cen = [], [], [], [], [], []
    for row, raw in zip(ds, d):
        cc = row[c - 128: c + 128].mean() / rho_b
        lo = row[:c][::-1] < 0.5 * rho_b
        hi = row[c:] < 0.5 * rho_b
        Wc.append((np.argmax(lo) if lo.any() else c) + (np.argmax(hi) if hi.any() else Ly - c))
        lo = row[:c][::-1] < 0.5 * rho_b * cc
        hi = row[c:] < 0.5 * rho_b * cc
        Wa.append((np.argmax(lo) if lo.any() else c) + (np.argmax(hi) if hi.any() else Ly - c))
        lo = np.nonzero(row[:c][::-1] > 0.3 * rho_b)[0]
        hi = np.nonzero(row[c:] > 0.3 * rho_b)[0]
        Wo.append((lo[-1] + 1 if len(lo) else 0) + (hi[-1] + 1 if len(hi) else 0))
        WA.append(raw.sum() / rho_b)
        WAn.append(raw.sum() / rho_b / max(cc, 0.05))
        cen.append(cc)
    return (t.astype(float), {"c": np.array(Wc, float), "a": np.array(Wa, float), "o": np.array(Wo, float),
                              "A": np.array(WA, float), "An": np.array(WAn, float)}, np.array(cen), Lx, Ly, rho_b)


def vfit(t, W, Ly, a, b, cen=None, cmin=None, box=True):
    sel = (t >= a) & (t <= b)
    if box:
        sel &= (W > 0.15 * Ly) & (W < 0.85 * Ly)
    if cmin is not None:
        sel &= cen > cmin
    if sel.sum() < 4:
        return np.nan, int(sel.sum()), 0.0
    return np.polyfit(t[sel], W[sel], 1)[0] / 2, int(sel.sum()), t[sel][-1]


def run_row(txt):
    P = fronts.params(os.path.basename(txt)[:-4])
    p, q = float(P["p"]), float(P["q"])
    t, W, cen, Lx, Ly, rho_b = track(txt)
    T = t[-1]
    tx = TX[p]
    dl = done_line(txt)
    absorbed = dl is not None and "absorbed 1" in dl
    r = dict(p=p, q=q, Lx=Lx, Ly=Ly, s=int(P["s"]), T=T, rho_b=rho_b, done=dl is not None, absorbed=absorbed,
             frac=strip_frac(txt), v={}, n={}, tL={})
    mid = (tx + max(T, tx)) / 2
    wins = {"E": (2e3, 1e5), "M": (1e5, tx), "L": (tx, 2e7), "L1": (tx, mid), "L2": (mid, 2e7)}
    for e in EST:
        for w, (a, b) in wins.items():
            v, n, tl = vfit(t, W[e], Ly, a, b, cen, 0.7 if (e == "c" and w == "E") else None)
            if w != "E" and n < NMIN:
                v = np.nan
            r["v"][(e, w)], r["n"][(e, w)], r["tL"][(e, w)] = v, n, tl
    r["cen"] = [np.interp(x, t, cen) if x <= T else np.nan for x in INTERIOR_T]
    q4 = t >= T - (T - tx) / 4 if T > tx else t >= 0.75 * T
    cq = cen[q4].mean() if q4.any() else cen[-1]
    cx = np.interp(tx, t, cen) if T >= tx else cen[-1]
    r["cen_x"], r["cen_q"] = cx, cq
    r["Wend"] = W["o"][-1] / Ly
    if absorbed or cq <= 0.5 or W["o"][-1] < 0.15 * Ly:
        st = "absorbed" if absorbed else "gone"
    elif W["o"][-1] >= 0.85 * Ly:
        st = "filled"
    elif cq < cx - 0.03:
        st = "decaying"
    else:
        st = "stationary"
    if not r["done"]:
        st += "*"
    r["status"] = st
    return r


def f(v):
    return "   nan  " if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:+.5f}"


def strips(dirs, out=sys.stdout):
    rows = []
    for d in dirs:
        for txt in sorted(glob.glob(os.path.join(d, "out", "*.txt"))):
            if os.path.exists(txt[:-4] + ".prof") and os.path.getsize(txt[:-4] + ".prof") > 0:
                try:
                    rows.append(run_row(txt))
                except Exception as e:  # partial file
                    print(f"skip {txt}: {e}", file=sys.stderr)
    rows.sort(key=lambda r: (r["p"], r["Lx"], r["q"], r["s"]))
    print("## per run: v per front in rows/step; outer-envelope estimator in windows E [2e3,1e5], M [1e5,t_x], "
          "L [t_x,T], L1/L2 = halves of L; v_L of the centre (Table I), adaptive, integrated and normalized-integrated estimators; n_L, t_L = "
          "records and last time in L (outer); interior = cen at t = "
          + ", ".join(f"{x:.2g}" for x in INTERIOR_T) + " (central 256 rows / rho_b at t ~ 1e3)", file=out)
    print("| p | q | Lx x Ly | s | T | v_E | v_M | v_L | v_L1 | v_L2 | v_L c | v_L a | v_L A | v_L An | n_L | t_L | W(T)/Ly | interior | status |", file=out)
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|", file=out)
    for r in rows:
        v = r["v"]
        print(f"| {r['p']} | {r['q']:.6g} | {r['Lx']}x{r['Ly']} | {r['s']} | {r['T']:.2g} | {f(v[('o', 'E')])} | {f(v[('o', 'M')])} "
              f"| {f(v[('o', 'L')])} | {f(v[('o', 'L1')])} | {f(v[('o', 'L2')])} | {f(v[('c', 'L')])} | {f(v[('a', 'L')])} | {f(v[('A', 'L')])} | {f(v[('An', 'L')])} "
              f"| {r['n'][('o', 'L')]} | {r['tL'][('o', 'L')]:.2g} | {r['Wend']:.2f} | "
              + " ".join("-" if np.isnan(x) else f"{x:.2f}" for x in r["cen"]) + f" | {r['status']} |", file=out)
    return rows


def seed_mean(vals):
    v = np.array([x for x in vals if np.isfinite(x)])
    if len(v) == 0:
        return np.nan, np.nan, 0
    return v.mean(), (v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan), len(v)


def line_fit(rs, ok, key):
    """Weighted line through the seed means at the q in ok; returns (q*, slope) or None."""
    x, y, wt = [], [], []
    for q in ok:
        vv = np.array([r["v"][key] for r in rs if r["q"] == q and np.isfinite(r["v"][key])])
        if len(vv) < 2:
            continue
        x.append(q); y.append(vv.mean()); wt.append(len(vv) / max(vv.var(ddof=1), 1e-10))
    if len(x) < 2:
        return None
    x, y, wt = map(np.array, (x, y, wt))
    A = np.vstack([np.ones_like(x), x - x.mean()]).T
    coef = np.linalg.lstsq(A * np.sqrt(wt)[:, None], y * np.sqrt(wt), rcond=None)[0]
    return x.mean() - coef[0] / coef[1], coef[1]


COMBOS = [("o", "E"), ("o", "M"), ("o", "L"), ("o", "L1"), ("o", "L2"), ("c", "E"), ("c", "L"), ("a", "L"), ("A", "L"), ("An", "L"), ("a", "E")]


def cross(rows, out=sys.stdout, nboot=2000, seed=0, usable=("stationary", "filled")):
    """Seed-mean tables and crossings.  Line fits use the q whose seeds are all in `usable`
    (status without the '*'), i.e. strips whose interior held; the two-point crossing uses every q."""
    rng = np.random.default_rng(seed)
    by = defaultdict(list)
    for r in rows:
        by[(r["p"], r["Lx"], r["Ly"])].append(r)
    res = {}
    for (p, Lx, Ly), rs in sorted(by.items()):
        qs = sorted(set(r["q"] for r in rs))
        print(f"\n## p = {p}, Lx x Ly = {Lx}x{Ly}: seed mean +- sem (n) of v (outer estimator) per window; "
              f"then v_L of the c, a, A, An estimators; status per seed", file=out)
        print("| q | x = q/q*_old | v_E | v_M | v_L | v_L1 | v_L2 | v_L c | v_L a | v_L A | v_L An | status |", file=out)
        print("|---|---|---|---|---|---|---|---|---|---|---|---|", file=out)
        tab = {}
        for q in qs:
            rq = [r for r in rs if r["q"] == q]
            cells = []
            for key in (("o", "E"), ("o", "M"), ("o", "L"), ("o", "L1"), ("o", "L2"), ("c", "L"), ("a", "L"), ("A", "L"), ("An", "L"), ("c", "E"), ("a", "E")):
                m, e, n = seed_mean([r["v"][key] for r in rq])
                if key[1] != "E" or key[0] == "o":
                    cells.append("-" if n == 0 else f"{m:+.5f} +- {e:.5f} ({n})" if n > 1 else f"{m:+.5f} ({n})")
                tab[(q, key)] = (m, e, n)
            st = "/".join(r["status"] for r in sorted(rq, key=lambda r: r["s"]))
            print(f"| {q:.6g} | {q / QS_OLD[p]:.3f} | " + " | ".join(cells) + f" | {st} |", file=out)
        for key in COMBOS:
            allq = [q for q in qs if tab[(q, key)][2] >= 2]
            ok = [q for q in allq if all(r["status"].rstrip("*") in usable for r in rs if r["q"] == q)]
            tp = None
            for q1, q2 in zip(allq[:-1], allq[1:]):
                v1, v2 = tab[(q1, key)][0], tab[(q2, key)][0]
                if v1 < 0 <= v2 and tp is None:
                    tp = q1 - v1 * (q2 - q1) / (v2 - v1)
            fit = line_fit(rs, ok, key) if len(ok) >= 2 else None
            entry = dict(twopt=tp, pts=ok)
            if fit is not None and np.isfinite(fit[0]):
                seeds = sorted(set(r["s"] for r in rs))
                bs = []
                for _ in range(nboot):
                    pick = rng.choice(seeds, len(seeds), replace=True)
                    sel = [r for s in pick for r in rs if r["s"] == s]
                    l = line_fit(sel, ok, key)
                    if l is not None and np.isfinite(l[0]):
                        bs.append(l[0])
                bs = np.array(bs)
                entry.update(qstar=fit[0], slope=fit[1], sd=bs.std() if len(bs) > 10 else np.nan)
            res[(p, Lx, Ly, key)] = entry
            name = f"{key[0]}:{key[1]}"
            print(f"- {name:5s} two-point crossing (all q): {'-' if tp is None else f'{tp:.6f}'}; "
                  + (f"line through {[f'{q:.6g}' for q in ok]}: q* = {entry['qstar']:.6f} +- {entry['sd']:.6f} "
                     f"(seed bootstrap), dv/dq = {entry['slope']:.1f}, q* dv/dq = {entry['qstar'] * entry['slope']:.3f}, "
                     f"c = {entry['qstar'] / (0.5 * (1 - p)):.4f}" if "qstar" in entry
                     else f"line: fewer than 2 usable points ({[f'{q:.6g}' for q in ok]})"), file=out)
    return res


def bulk(dirs, out=sys.stdout):
    times = (1e3, 1e4, 1e5, 2e5, 5e5, 1e6, 1.2e6, 1.5e6)
    print("## bulk: rho_bond per run at fixed t; last = mean over the last quarter of the run, drift = slope over "
          "that quarter per 1e5 steps; stationary = |drift| x (quarter length) < 0.01 and last > 0.15", file=out)
    print("| run (dir) | T | " + " | ".join(f"{x:.1e}" for x in times) + " | last | drift | status |", file=out)
    print("|---|---|" + "---|" * len(times) + "---|---|---|", file=out)
    rows = []
    for d in dirs:
        for txt in sorted(glob.glob(os.path.join(d, "out", "*.txt"))):
            name = os.path.basename(txt)[:-4]
            if "_L" not in name or "Lx" in name:
                continue
            ts = fronts.load_ts(txt)
            if len(ts) < 5:
                continue
            P = fronts.params(name)
            p, q = float(P["p"]), float(P["q"])
            if p not in (0.99, 0.995):
                continue
            T = ts[-1, 0]
            dl = done_line(txt)
            absorbed = dl is not None and "absorbed 1" in dl
            r = [np.interp(x, ts[:, 0], ts[:, 2]) if x <= T else np.nan for x in times]
            q4 = ts[:, 0] >= 0.75 * T
            last = ts[q4, 2].mean()
            drift = np.polyfit(ts[q4, 0], ts[q4, 2], 1)[0] * 1e5 if q4.sum() > 3 else np.nan
            st = "absorbed" if absorbed else ("stationary" if abs(drift) * (T * 0.25 / 1e5) < 0.01 and last > 0.15 else "decaying")
            if dl is None:
                st += "*"
            rows.append(dict(p=p, q=q, L=int(P["L"]), s=int(P["s"]), T=T, last=last, drift=drift, status=st, src=os.path.basename(d)))
            print(f"| {name} ({os.path.basename(d)}) | {T:.2g} | " + " | ".join("-" if np.isnan(x) else f"{x:.3f}" for x in r)
                  + f" | {last:.3f} | {drift:+.4f} | {st} |", file=out)
    return rows


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    sdirs = [f"{D}/strips", f"{D}/wide", f"{D}/tall"]
    bdirs = [f"{D}/bulk"] + E6[0.995] + E6[0.99]
    if what == "all":
        with open(f"{D}/analysis.md", "w") as out:
            print("# analysis.md (regenerated by `python3 analyze.py all`)\n", file=out)
            rows = strips(sdirs, out=out)
            cross(rows, out=out)
            print("\n", file=out)
            bulk(bdirs, out=out)
        print(open(f"{D}/analysis.md").read())
    elif what == "strips":
        strips(sdirs)
    elif what == "cross":
        cross(strips(sdirs, out=open(os.devnull, "w")))
    elif what == "bulk":
        bulk(bdirs)
