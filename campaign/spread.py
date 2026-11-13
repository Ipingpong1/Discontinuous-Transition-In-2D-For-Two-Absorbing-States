#!/usr/bin/env python3
"""E8a seed-spreading analysis (dnsim --reps output).

    python3 spread.py FILE.txt [...]          one table per file
    python3 spread.py --summary RUNS_DIR      full tables for every run in RUNS_DIR/E8a_spread/out
    python3 spread.py --brief RUNS_DIR        one line per run
    python3 spread.py --exp E8a_long --brief RUNS_DIR    (any experiment dir; default E8a_spread)
Runs whose names differ only in a trailing _k<n> (E8a_long sub-runs, separate seeds) are pooled:
their groups are stacked, so the jackknife runs over all of them.
Lazy runs (--ftab, _r<r> in the name, E11_lazy_seed) are read in physical time t r, and their
q/q* is relative to q*(r).

Per output time t, pooled over the G groups (errors: jackknife over groups):
  P_s(t) = surviving replicas / replicas (escaped replicas count as surviving);
  N(t)   = <I>, mean interface cells over all replicas (dead ones count 0);
  A(t)   = <A>, mean live cells (pure voter: exactly 4, a martingale);
  R2(t)  = sum r^2 over interface cells / sum I, the spread of the active region;
  local slopes over a factor `--w` in t (default 4): delta_eff = -dln P_s/dln t,
  eta_eff = dln N/dln t, z_eff = dln R2/dln t.
Escaped replicas stop contributing to N, A and R2, so those are flagged once escapes exceed
1% of the survivors.  Reference: DP delta 0.451, eta 0.230, z 1.13 (R2 ~ t^z);
GV / voter: P_s ~ ln t / t (delta_eff -> 1 slowly), eta 0, z 1.
"""
import glob
import os
import re
import sys

import numpy as np

COLS = ["n", "surv", "esc", "I", "I2", "A", "A2", "IR2", "area"]


def load1(path):
    rows = np.array([[float(x) for x in l.split()] for l in open(path) if l.strip() and not l.startswith("#")])
    G = int(rows[:, 0].max()) + 1
    t = np.unique(rows[:, 1])
    data = rows[:, 2:].reshape(G, len(t), len(COLS))
    meta = {}
    nm = re.match(r"p([\d.]+)_a[\d.]+_r([\d.]+)_j2_q([\d.]+)_", os.path.basename(path))
    for l in open(path):
        if l.startswith("# ftab") and nm:
            meta["p"], meta["jflip"], meta["r"] = float(nm[1]), f"2:{nm[3]}", float(nm[2])
        elif l.startswith("# rule"):
            mt = re.search(r"voter (\S+) jflip (\S+)", l)
            meta["p"], meta["jflip"] = float(mt[1]), mt[2]
        elif l.startswith("# replicas"):
            meta["status"] = l[2:].strip()
        elif l.startswith("# done"):
            meta["done"] = True
    return t * meta.get("r", 1), data, meta


def load(paths):
    """One file, or a list of sub-runs of one point (same t grid): groups stacked."""
    if isinstance(paths, str):
        return load1(paths)
    parts = [load1(p) for p in paths]
    t = parts[0][0]
    if any(len(pt[0]) != len(t) or (pt[0] != t).any() for pt in parts):
        raise ValueError(f"t grids differ: {paths}")
    meta = dict(parts[0][2])
    meta["done"] = all(pt[2].get("done") for pt in parts)
    meta["status"] = f"{len(parts)} sub-runs, {sum(pt[2].get('done', False) for pt in parts)} done"
    return t, np.concatenate([pt[1] for pt in parts]), meta


def points(files):
    """Group files by name minus a trailing _k<n>; returns [(label, file or [files])]."""
    by = {}
    for f in files:
        by.setdefault(re.sub(r"_k\d+\.txt$", ".txt", os.path.basename(f)), []).append(f)
    return [(k, v[0] if len(v) == 1 else sorted(v)) for k, v in by.items()]


def stats(data, sel=None):
    """Pooled observables from per-group sums (all groups, or those with sel True)."""
    d = data if sel is None else data[sel]
    S = d.sum(axis=0)
    n = S[:, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        return {"Ps": S[:, 1] / n, "N": S[:, 3] / n, "A": S[:, 5] / n, "R2": S[:, 7] / S[:, 3],
                "esc": S[:, 2], "surv": S[:, 1]}


def jack(data, key):
    G = data.shape[0]
    full = stats(data)[key]
    if G < 2:
        return full, np.full_like(full, np.nan)
    js = np.array([stats(data, np.arange(G) != g)[key] for g in range(G)])
    return full, np.sqrt((G - 1) / G * ((js - js.mean(axis=0)) ** 2).sum(axis=0))


def slope(t, y, w):
    """d ln y / d ln t between t/sqrt(w) and t*sqrt(w), by interpolation in ln t."""
    ok = (t > 0) & (y > 0) & np.isfinite(y)
    lt, ly = np.log(t[ok]), np.log(y[ok])
    out = np.full(len(t), np.nan)
    h = np.log(w) / 2
    for i, ti in enumerate(t):
        if ti <= 0 or np.log(ti) - h < lt[0] or np.log(ti) + h > lt[-1]:
            continue
        a, b = np.interp(np.log(ti) - h, lt, ly), np.interp(np.log(ti) + h, lt, ly)
        out[i] = (b - a) / (2 * h)
    return out


def table(path, w, label=None):
    t, data, meta = load(path)
    st = stats(data)
    Ps, ePs = jack(data, "Ps")
    N, eN = jack(data, "N")
    A, eA = jack(data, "A")
    R2, eR2 = jack(data, "R2")
    dl, et, zz = -slope(t, Ps, w), slope(t, N, w), slope(t, R2, w)
    print(f"== {label or os.path.basename(path)}  p={meta.get('p')} jflip={meta.get('jflip')}  {meta.get('status', '')}"
          f"{'' if meta.get('done') else '  (running)'}")
    print(f"{'t':>8} {'P_s':>10} {'±':>8} {'N':>9} {'±':>7} {'A':>8} {'±':>6} {'R2':>9} {'d_eff':>6} {'eta':>6} "
          f"{'z_eff':>6} {'tPs/lnt':>8} esc")
    for i in range(len(t)):
        if t[i] == 0:
            continue
        flag = "*" if st["esc"][i] > 0.01 * max(st["surv"][i], 1) else ""
        bg = t[i] * Ps[i] / np.log(t[i]) if t[i] > 1 else np.nan
        print(f"{t[i]:8.0f} {Ps[i]:10.4g} {ePs[i]:8.2g} {N[i]:9.4g} {eN[i]:7.2g} {A[i]:8.4g} {eA[i]:6.2g} "
              f"{R2[i]:9.4g} {dl[i]:6.3f} {et[i]:6.3f} {zz[i]:6.3f} {bg:8.4g} {st['esc'][i]:.0f}{flag}")


QSTAR = {0.9: 0.0119, 0.98: 0.00252, 0.99: 0.00127, 0.995: 0.00062}
QSTAR_LAZY = {(0.98, 0.5): 0.0021773}   # make_tasks.QSTAR_LAZY (E10 fronts)
BRIEF_T = [1e3, 1e4, 3e4, 1e5]
BRIEF_T_LONG = [1e4, 1e5, 3e5, 1e6]


def brief(files, w):
    """One line per point.  Runs to T >= 3e5 get P_s at 1e4-1e6, d_eff at 3e4, 1e5, 2.5e5 and
    the ratio P_s(T)/P_s(1e5) (voter, from ln t / t: 0.12 for T = 1e6, 0.37 for T = 3e5)."""
    for label, f in points(files):
        t, data, meta = load(f)
        long_ = t[-1] >= 3e5
        BT = BRIEF_T_LONG if long_ else BRIEF_T
        DT = (3e4, 1e5, 2.5e5) if long_ else (3e3, 1e4, 3e4)
        if label == points(files)[0][0]:
            print(f"{'p':>5} {'q/q*':>6} {'reps':>8} {'alive':>5} | P_s at t = " + ", ".join(f"{x:.0e}" for x in BT) +
                  " | d_eff at t = " + ", ".join(f"{x:.1e}" for x in DT) + " | eta_eff at the same t | " +
                  ("P_s(T)/P_s(1e5) ± jack" if long_ else "A(1e4)"))
        Ps, ePs = jack(data, "Ps")
        st = stats(data)
        dl, et = -slope(t, Ps, w), slope(t, st["N"], w)
        at = lambda y, T: y[np.argmin(np.abs(t - T))]
        p, jf = meta["p"], meta["jflip"]
        q = float(jf.split(":")[1]) if ":" in jf else 0.0
        rr = meta.get("r", 1)
        qs = QSTAR.get(p) if rr == 1 else QSTAR_LAZY.get((p, rr))
        x = f"{q / qs:6.3f}" if qs else f"q{q:<5g}"
        n = int(data[:, -1, 0].sum())
        if long_:
            G = data.shape[0]
            i5, i6 = np.argmin(np.abs(t - 1e5)), len(t) - 1
            def ratio(d):
                S = d.sum(axis=0)
                return S[i6, 1] / S[i5, 1] if S[i5, 1] else np.nan
            r = ratio(data)
            js = np.array([ratio(data[np.arange(G) != g]) for g in range(G)])
            er = np.sqrt((G - 1) / G * np.nansum((js - np.nanmean(js)) ** 2))
            tail = f" | {r:5.3f} ± {er:5.3f}" + ("" if t[-1] >= 1e6 else f" (T {t[-1]:.0e})")
        else:
            tail = f" | {at(st['A'], 1e4):5.2f}"
        print(f"{p:5g} {x:>6} {n:8d} {int(st['surv'][-1]):5d} | " +
              " ".join(f"{at(Ps, T):8.2e}" if T <= t[-1] else f"{'-':>8}" for T in BT) + " | " +
              " ".join(f"{at(dl, T):5.2f}" for T in DT) + " | " +
              " ".join(f"{at(et, T):5.2f}" for T in DT) + tail +
              ("" if rr == 1 else f"  [r {rr:g}]") +
              ("" if meta.get("done") else f"  (running: {meta.get('status', '')})"))


def main():
    w = 4.0
    args = sys.argv[1:]
    if "--w" in args:
        i = args.index("--w")
        w = float(args[i + 1])
        del args[i:i + 2]
    exp = "E8a_spread"
    if "--exp" in args:
        i = args.index("--exp")
        exp = args[i + 1]
        del args[i:i + 2]
    key = lambda f: [float(v) for v in re.findall(r"_[pqk]([\d.]+)", os.path.basename(f))]
    if args and args[0] == "--brief":
        files = sorted(glob.glob(os.path.join(args[1], exp, "out", "*.txt")), key=key)
        return brief(files, w)
    if args and args[0] == "--summary":
        files = sorted(glob.glob(os.path.join(args[1], exp, "out", "*.txt")), key=key)
    else:
        files = args
    for label, f in points(files):
        table(f, w, label)


if __name__ == "__main__":
    main()
