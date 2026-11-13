#!/usr/bin/env python3
"""Analysis for the seed-averaged density campaign (NOTES.md Sec. 9; raw data in bulk_dens/).

    python3 mdens.py runs            per run: mean over the last half of the run and over its two
                                     quarters, drift, autocorrelation time, absorption
    python3 mdens.py points [lo,hi]  seed-averaged rho(dp_1) per (p_1, L) with errors, the quarter
                                     check and the drift; writes bulk_dens/points.txt, or, with a
                                     fixed window lo,hi (steps), bulk_dens/points_<lo>-<hi>.txt from
                                     every run that covers the window (finished or not)
    python3 mdens.py fit [file] [--all] [--L n]
                                     fits of rho(dp_1) (log, power with rho0 >= 0, jump with beta
                                     fixed), p_1* fixed and free, from bulk_dens/points.txt (or file)
                                     plus the Sec. 4 points at dp_1 >= 4e-4; only points flagged
                                     stationary unless --all; L 2048 unless --L
    python3 mdens.py slopes [file] [--all] [--L n]
                                     finite-difference slopes between adjacent points vs the fits
    python3 mdens.py blocks [dmax]   means over consecutive 1e5-step blocks of every run with
                                     dp_1 <= dmax (default 0.5e-4): the seed-by-seed history
    python3 mdens.py all

Window: the last half [T/2, T] of each run (T = length reached by a finished run, nominal length
otherwise), or a fixed window [lo, hi], sampled on the uniform output grid (every 2000 steps); its
quarters are the two halves of the window.  Per-run error: std x sqrt(2 tau / n)
with tau the integrated autocorrelation time (Sokal windowing, c = 5).  Per-point error: the larger
of the seed-to-seed standard error and the combined per-run error, with a floor of 0.001.
Fits are weighted least squares, chi^2 = sum ((rho - f) / err)^2; rms is the plain rms residual.
"""
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

P1STAR = 0.18535           # front reversal on the cut p_2 = 0.45 (NOTES.md Sec. 8.1), +-3e-5
P1STAR_ERR = 0.00003
ERR_FLOOR = 0.001
OLD_ERR = 0.003            # Sec. 4 single runs (L 2048, 6e4 steps, mean over the last third)
OLD_P1 = ["0.18575", "0.186", "0.187", "0.188", "0.190"]
KEY = re.compile(r"p2_(?P<p2>[0-9.]+)_p1_(?P<p1>[0-9.]+)_L(?P<L>[0-9]+)_s(?P<s>[0-9]+)")
WINDOWS = ((1.5e5, 3e5), (3e5, 6e5), (6e5, 1.2e6))
WLAB = ["1.5e5-3e5", "3e5-6e5", "6e5-1.2e6"]


def header(path):
    h = {}
    for line in open(path):
        if not line.startswith("#"):
            continue
        w = line.split()
        if line.startswith("# done"):
            h["t_end"], h["absorbed"] = int(w[3]), int(w[5])
            h["wall_s"] = float(w[7])
        elif line.startswith("# Lx"):
            h["Lx"], h["Ly"], h["T"], h["seed"] = int(w[2]), int(w[4]), int(w[6]), int(w[8])
    return h


def series(path):
    a = np.array([[float(x) for x in l.split()] for l in open(path) if not l.startswith("#")])
    return a[:, 0], a[:, 2]


def uniform(t, rho, lo, hi):
    """Points of the uniform output grid (the most common spacing) within [lo, hi]."""
    sel = (t >= lo) & (t <= hi)
    if sel.sum() < 4:
        return t[sel], rho[sel]
    d = np.diff(t[sel])
    every = int(np.bincount(d.astype(int)).argmax())
    sel &= (t % every == 0)
    return t[sel], rho[sel]


def tau_int(x, c=5.0):
    """Integrated autocorrelation time in samples (0.5 for white noise) by Sokal's windowing."""
    x = x - x.mean()
    n = len(x)
    v = (x * x).mean()
    if n < 8 or v == 0:
        return 0.5
    tau = 0.5
    for M in range(1, n // 2):
        tau += np.dot(x[:n - M], x[M:]) / (n - M) / v
        if M >= c * tau:
            break
    return max(tau, 0.5)


def analyse(path, window=None):
    """Statistics of one run over its last half, or over the fixed window (lo, hi) if given (NaN if the
    run has not reached hi, or absorbed before it)."""
    h = header(path)
    t, rho = series(path)
    T = h["t_end"] if "t_end" in h and not h.get("absorbed") else h["T"]
    r = dict(file=os.path.basename(path), L=h["Lx"], T=T, seed=h["seed"], t_last=t[-1],
             done="t_end" in h, absorbed=h.get("absorbed", 0), t_end=h.get("t_end", t[-1]), wall=h.get("wall_s", np.nan))
    m = KEY.search(r["file"])
    r["p1"], r["p2"] = float(m["p1"]), float(m["p2"])
    lo, hi = window if window else (T / 2, T)
    r["window"] = (lo, hi)
    tw, rw = uniform(t, rho, lo, hi)
    r["n"] = len(tw)
    r["win"] = {}
    for wlo, whi in WINDOWS:                      # fixed windows for the comparisons across L and time
        tm, rm = uniform(t, rho, wlo, whi)
        r["win"][(wlo, whi)] = rm.mean() if len(rm) >= 8 and t[-1] >= whi and not (r["absorbed"] and r["t_end"] < whi) else np.nan
    r["mean_mid"] = r["win"][WINDOWS[0]]
    r["partial"] = (not r["done"]) and t[-1] < T
    if (r["absorbed"] and r["t_end"] < hi) or len(tw) < 8 or t[-1] < hi:
        r.update(mean=np.nan, q1=np.nan, q2=np.nan, drift=np.nan, std=np.nan, tau=np.nan, se=np.nan)
        return r
    r["mean"], r["std"] = rw.mean(), rw.std(ddof=1)
    q1, q2 = rw[tw < (lo + hi) / 2], rw[tw >= (lo + hi) / 2]
    r["q1"], r["q2"] = (q1.mean() if len(q1) else np.nan), (q2.mean() if len(q2) else np.nan)
    r["drift"] = np.polyfit(tw, rw, 1)[0] * 1e5
    r["tau"] = tau_int(rw)                        # in samples
    r["dt"] = tw[1] - tw[0]
    r["se"] = r["std"] * np.sqrt(2 * r["tau"] / len(rw))
    return r


def load_runs(pattern="bulk_dens/p2_*.txt", window=None):
    return [analyse(f, window) for f in sorted(glob.glob(pattern)) if KEY.search(os.path.basename(f))]


def parse_window(args):
    """'3e5,6e5' -> (300000.0, 600000.0); None if the first argument is not a window."""
    if args and "," in args[0] and not args[0].startswith("-"):
        lo, hi = (float(x) for x in args[0].split(","))
        return (lo, hi)
    return None


def cmd_runs(_):
    rs = load_runs()
    print("| run | L | T | t reached | status | mean (last half) | Q1 | Q2 | Q2-Q1 | drift /1e5 | std | tau (steps) | se | " + " | ".join("mean " + w for w in WLAB) + " | wall s |")
    print("|" + "---|" * 17)
    for r in sorted(rs, key=lambda r: (r["p1"], r["L"], r["seed"])):
        st = "absorbed at %.3g" % r["t_end"] if r["absorbed"] else ("done" if r["done"] else "running")
        if np.isnan(r["mean"]):
            print(f"| {r['file'][:-4]} | {r['L']} | {r['T']:.3g} | {r['t_last']:.3g} | {st} | - | - | - | - | - | - | - | - | - | - | - | {r['wall']:.0f} |")
            continue
        wins = " | ".join(f"{r['win'][w]:.4f}" if not np.isnan(r['win'][w]) else "-" for w in WINDOWS)
        print(f"| {r['file'][:-4]} | {r['L']} | {r['T']:.3g} | {r['t_last']:.3g} | {st} | {r['mean']:.4f} | {r['q1']:.4f} | {r['q2']:.4f} | "
              f"{r['q2'] - r['q1']:+.4f} | {r['drift']:+.4f} | {r['std']:.4f} | {r['tau'] * r['dt']:.0f} | {r['se']:.4f} | "
              f"{wins} | {r['wall']:.0f} |")
    done = [r for r in rs if r["done"]]
    print(f"\n{len(done)}/{len(rs)} runs finished; summed wall time of finished runs {sum(r['wall'] for r in done) / 3600:.1f} job-hours")


def points(rs, require_done=True):
    """Seed averages per (p_1, L).  With require_done only finished runs count; otherwise every run
    with a complete window (absorbed runs count as absorbed either way)."""
    g = defaultdict(list)
    for r in rs:
        if require_done and not r["done"]:
            continue
        g[(r["p1"], r["L"])].append(r)
    out = []
    for (p1, L), runs in sorted(g.items()):
        ok = [r for r in runs if not np.isnan(r["mean"])]
        absd = [r for r in runs if r["absorbed"]]
        p = dict(p1=p1, dp1=p1 - P1STAR, L=L, n=len(ok), n_abs=len(absd), T=runs[0]["T"],
                 t_abs=[r["t_end"] for r in absd], window=runs[0]["window"],
                 seed_means=sorted((r["mean"], r["seed"]) for r in ok))
        if len(ok) == 0:
            out.append(p)
            continue
        m = np.array([r["mean"] for r in ok])
        d = np.array([r["q2"] - r["q1"] for r in ok])
        dr = np.array([r["drift"] for r in ok])
        p["rho"] = m.mean()
        p["sd"] = m.std(ddof=1) if len(ok) > 1 else np.nan
        p["se_seed"] = p["sd"] / np.sqrt(len(ok)) if len(ok) > 1 else np.nan
        p["se_tau"] = np.sqrt(sum(r["se"] ** 2 for r in ok)) / len(ok)
        p["err"] = max(np.nan_to_num(p["se_seed"]), p["se_tau"], ERR_FLOOR)
        p["dq"] = d.mean()
        p["dq_se"] = d.std(ddof=1) / np.sqrt(len(ok)) if len(ok) > 1 else np.sqrt(2) * p["se_tau"]
        p["drift"] = dr.mean()
        p["drift_se"] = dr.std(ddof=1) / np.sqrt(len(ok)) if len(ok) > 1 else np.nan
        p["mid"] = np.nanmean([r["mean_mid"] for r in ok])
        p["win"] = {w: (np.nanmean([r["win"][w] for r in ok]) if not all(np.isnan(r["win"][w]) for r in ok) else np.nan) for w in WINDOWS}
        p["win_n"] = {w: sum(not np.isnan(r["win"][w]) for r in ok) for w in WINDOWS}
        p["tau"] = np.mean([r["tau"] * r["dt"] for r in ok])
        # seed consistency: spread of the seed means in units of the within-run standard deviation
        # (seeds sampling one stationary state spread by about std x sqrt(2 tau / n) << std)
        p["std_run"] = np.mean([r["std"] for r in ok])
        p["spread"] = (m.max() - m.min()) / p["std_run"]
        # status: 'absorbing' if a seed reached an absorbing state; 'split' if the seed means spread by
        # more than 3 within-run standard deviations (the seeds are not sampling one state); 'drifting'
        # if the quarters differ by more than 2 standard errors or the drift exceeds 2.5 se; else
        # 'stationary'
        if len(absd):
            p["status"] = "absorbing"
        elif p["spread"] > 3:
            p["status"] = "split"
        elif not (abs(p["dq"]) <= 2 * max(p["dq_se"], np.sqrt(2) * p["se_tau"]) and (np.isnan(p["drift_se"]) or abs(p["drift"]) <= 2.5 * p["drift_se"])):
            p["status"] = "drifting"
        else:
            p["status"] = "stationary"
        p["stationary"] = p["status"] == "stationary"
        out.append(p)
    return out


def cmd_points(args, write=True):
    window = parse_window(args)
    rs = load_runs(window=window)
    pts = points(rs, require_done=window is None)
    wlab = "last half" if window is None else f"{int(window[0])}-{int(window[1])}"
    out = "bulk_dens/points.txt" if window is None else f"bulk_dens/points_{int(window[0])}-{int(window[1])}.txt"
    print(f"| p_1 | dp_1 | L | T | seeds (absorbed) | rho ({wlab}) | seed sd | se (seeds) | se (tau) | err used | Q2 - Q1 | drift /1e5 | tau (steps) | "
          + " | ".join("mean " + w for w in WLAB) + " | seed spread / std | status | seed means (sorted) |")
    print("|" + "---|" * 19)
    lines = []
    for p in pts:
        if "rho" not in p:
            print(f"| {p['p1']:.5f} | {p['dp1']:+.5f} | {p['L']} | {p['T']:.3g} | {p['n']} ({p['n_abs']}) | "
                  + (("all absorbed at " + ", ".join(f"{x:.2g}" for x in p["t_abs"])) if p["t_abs"] else "no run covers the window")
                  + " | | | | | | | | | | | | - | |")
            continue
        sd = f"{p['sd']:.4f}" if not np.isnan(p["sd"]) else "-"
        ses = f"{p['se_seed']:.4f}" if not np.isnan(p["se_seed"]) else "-"
        dse = f"{p['drift_se']:.4f}" if not np.isnan(p["drift_se"]) else "-"
        print(f"| {p['p1']:.5f} | {p['dp1']:+.5f} | {p['L']} | {p['T']:.3g} | {p['n']} ({p['n_abs']}) | {p['rho']:.4f} | {sd} | {ses} | "
              f"{p['se_tau']:.4f} | {p['err']:.4f} | {p['dq']:+.4f} +- {p['dq_se']:.4f} | {p['drift']:+.4f} +- {dse} | {p['tau']:.0f} | "
              + " | ".join(f"{p['win'][w]:.4f} ({p['win_n'][w]})" if not np.isnan(p['win'][w]) else "-" for w in WINDOWS)
              + f" | {p['spread']:.1f} | {p['status']} | " + " ".join(f"{m:.3f}" for m, s in p["seed_means"]) + " |")
        lines.append(f"{p['p1']:.5f} {p['dp1']:+.5f} {p['L']} {p['n']} {p['rho']:.5f} {p['err']:.5f} {p['dq']:+.5f} {p['drift']:+.5f} "
                     f"{p['status']} {p['n_abs']}\n")
    if write:
        with open(out, "w") as f:
            f.write(f"# p1 dp1 L n rho err q2-q1 drift status n_absorbed   (mdens.py points, window {wlab}; p1* = {P1STAR:.5f})\n")
            f.writelines(lines)
        print(f"\n(written to {out})")
    return pts


# ---------------------------------------------------------------- fits

def old_points():
    """Sec. 4 single runs at dp_1 >= 4e-4 (bulk/, L 2048, 6e4 steps): mean over the last third."""
    out = []
    for p1 in OLD_P1:
        f = f"bulk/p2_0.45_p1_{p1}_s1.txt"
        t, rho = series(f)
        sel = t >= 2 * t[-1] / 3
        out.append((float(p1), rho[sel].mean(), OLD_ERR, "Sec. 4"))
    return out


def read_points(path="bulk_dens/points.txt", L=2048, all_status=False):
    out = []
    for line in open(path):
        if line.startswith("#"):
            continue
        w = line.split()
        if int(w[2]) == L and (all_status or w[8] == "stationary"):
            out.append((float(w[0]), float(w[4]), float(w[5]), f"L {L}" + ("" if w[8] == "stationary" else f" ({w[8]})")))
    return out


def fit_args(args):
    """[file] [--all] [--L n] -> (file, all_status, L)"""
    path, all_status, L = "bulk_dens/points.txt", False, 2048
    a = list(args)
    while a:
        x = a.pop(0)
        if x == "--all":
            all_status = True
        elif x == "--L":
            L = int(a.pop(0))
        else:
            path = x
    return path, all_status, L


def wlsq(X, y, w):
    """Weighted least squares; returns coefficients, covariance, chi^2."""
    Xw = X * np.sqrt(w)[:, None]
    yw = y * np.sqrt(w)
    c, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    cov = np.linalg.inv(Xw.T @ Xw)
    return c, cov, float(((yw - Xw @ c) ** 2).sum())


BETAS = np.concatenate([np.linspace(0.01, 0.4, 157), np.linspace(0.41, 1.5, 110)])


def fit_log(dp, rho, w):
    """rho = A / ln(B / dp): grid over ln B, A linear."""
    lmax = np.log(dp).max()
    best = None
    for grid in (np.linspace(lmax + 0.05, lmax + 30, 600), None):
        if grid is None:
            grid = np.linspace(best[1] - 0.06, best[1] + 0.06, 241)
        for lb in grid:
            g = 1 / (lb - np.log(dp))
            A = (w * g * rho).sum() / (w * g * g).sum()
            chi = float((w * (rho - A * g) ** 2).sum())
            if best is None or chi < best[0]:
                best = (chi, lb, A)
    chi, lb, A = best
    f = lambda x: A / (lb - np.log(x))
    return dict(form="log", chi=chi, k=2, A=A, B=np.exp(lb), f=f, label=f"A = {A:.3f}, B = {np.exp(lb):.3g}")


def fit_pow(dp, rho, w, rho0=None, beta=None):
    """rho = rho0 + C dp^beta with rho0 >= 0 (rho0 = None: free; a number: fixed), beta free or fixed."""
    best = None
    for b in (BETAS if beta is None else [beta]):
        xb = dp ** b
        if rho0 is None:
            c, cov, chi = wlsq(np.vstack([np.ones_like(dp), xb]).T, rho, w)
            r0, C, r0e = c[0], c[1], np.sqrt(cov[0, 0])
            if r0 < 0:
                C = (w * xb * rho).sum() / (w * xb * xb).sum()
                r0, r0e = 0.0, np.nan
                chi = float((w * (rho - C * xb) ** 2).sum())
        else:
            C = (w * xb * (rho - rho0)).sum() / (w * xb * xb).sum()
            r0, r0e = rho0, 0.0
            chi = float((w * (rho - rho0 - C * xb) ** 2).sum())
        if best is None or chi < best[0]:
            best = (chi, b, r0, C, r0e)
    chi, b, r0, C, r0e = best
    k = 1 + (rho0 is None) + (beta is None)
    f = lambda x: r0 + C * x ** b
    name = "pow0" if rho0 == 0 else ("jump%.1f" % beta if beta is not None else "pow")
    return dict(form=name, chi=chi, k=k, rho0=r0, rho0_err=r0e, C=C, beta=b, f=f,
                label=f"rho0 = {r0:.4f}" + (f" +- {r0e:.4f}" if r0e and not np.isnan(r0e) else "") + f", C = {C:.3f}, beta = {b:.3f}")


def rho0_profile(dp, rho, w, chimin):
    """Profile chi^2 over rho0 (beta, C free): intervals at dchi^2 = 1 and 4 (rho0 >= 0)."""
    grid = np.linspace(0, 0.25, 126)
    prof = np.array([fit_pow(dp, rho, w, rho0=r0)["chi"] for r0 in grid])
    out = {}
    for d in (1, 4):
        ok = grid[prof <= chimin + d]
        out[d] = (ok.min(), ok.max()) if len(ok) else (np.nan, np.nan)
    return out


def rms(dp, rho, f):
    return np.sqrt(np.mean((rho - f(dp)) ** 2))


def fit_set(pts, p1s, label, show_profile=True):
    """All forms at fixed p_1* for the points pts = [(p1, rho, err, src)]; prints a table block."""
    pts = [p for p in pts if p[0] - p1s > 1e-9]
    p1 = np.array([p[0] for p in pts]); rho = np.array([p[1] for p in pts]); err = np.array([p[2] for p in pts])
    dp = p1 - p1s; w = 1 / err ** 2
    n = len(dp)
    fits = [fit_log(dp, rho, w), fit_pow(dp, rho, w, rho0=0), fit_pow(dp, rho, w),
            fit_pow(dp, rho, w, beta=0.3), fit_pow(dp, rho, w, beta=0.5)]
    print(f"\n{label}: p_1* = {p1s:.5f}, n = {n} points, dp_1 = {dp.min():.1e} .. {dp.max():.1e}, median err {np.median(err):.4f}")
    print("| form | parameters | chi^2 | dof | chi^2/dof | rms | largest |resid|/err |")
    print("|---|---|---|---|---|---|---|")
    names = dict(log="A/ln(B/dp)", pow0="C dp^beta (rho0 = 0)", pow="rho0 + C dp^beta, rho0 >= 0",
                 **{"jump0.3": "rho0 + C dp^0.3 (jump, beta fixed)", "jump0.5": "rho0 + C dp^0.5 (jump, beta fixed)"})
    for F in fits:
        dof = n - F["k"]
        res = (rho - F["f"](dp)) / err
        print(f"| {names[F['form']]} | {F['label']} | {F['chi']:.1f} | {dof} | {F['chi'] / dof:.2f} | {rms(dp, rho, F['f']):.4f} | {np.abs(res).max():.1f} |")
    if show_profile:
        pr = rho0_profile(dp, rho, w, fits[2]["chi"])
        print(f"  rho0 profile (beta, C free): dchi^2 <= 1: [{pr[1][0]:.3f}, {pr[1][1]:.3f}]; dchi^2 <= 4: [{pr[4][0]:.3f}, {pr[4][1]:.3f}]")
    return fits, dp, rho, err


def free_p1star(pts, dmax, label):
    """Profile of chi^2 over p_1* for each form; points with p_1 - P1STAR <= dmax."""
    pts = [p for p in pts if p[0] - P1STAR <= dmax + 1e-9]
    p1 = np.array([p[0] for p in pts]); rho = np.array([p[1] for p in pts]); w = 1 / np.array([p[2] for p in pts]) ** 2
    grid = np.arange(0.18500, p1.min() - 1e-6, 2e-6)
    forms = {"log": lambda d: fit_log(d, rho, w), "pow0": lambda d: fit_pow(d, rho, w, rho0=0), "pow": lambda d: fit_pow(d, rho, w),
             "jump0.3": lambda d: fit_pow(d, rho, w, beta=0.3)}
    print(f"\n{label}: p_1* free on a 2e-6 grid from 0.18500 to {grid[-1]:.6f}, n = {len(p1)} points with dp_1 <= {dmax:.0e} (dp_1 from p_1* = {P1STAR})")
    print("| form | best p_1* | dchi^2 <= 1 range | chi^2 at best | chi^2 at p_1* = 0.18535 | parameters at best |")
    print("|---|---|---|---|---|---|")
    for name, fn in forms.items():
        prof = np.array([fn(p1 - ps)["chi"] for ps in grid])
        i = prof.argmin()
        ok = grid[prof <= prof[i] + 1]
        F = fn(p1 - grid[i])
        at = fn(p1 - P1STAR)["chi"]
        print(f"| {name} | {grid[i]:.6f} | [{ok.min():.6f}, {ok.max():.6f}] | {prof[i]:.1f} | {at:.1f} | {F['label']} |")


def cmd_fit(args):
    path, all_status, L = fit_args(args)
    new = read_points(path, L, all_status)
    old = old_points()
    print(f"# new points from {path}" + (" (all rows, including non-stationary ones)" if all_status else " (stationary rows only)"))
    allp = sorted(new + old)
    print("# points: " + ", ".join(f"({p[0]:.5f}: {p[1]:.4f} +- {p[2]:.4f}, {p[3]})" for p in allp))
    for p1s in (P1STAR, P1STAR - P1STAR_ERR, P1STAR + P1STAR_ERR):
        ranges = ((3.0e-4, "new points only"), (1.0e-3, "new + Sec. 4 up to 1e-3"), (5e-3, "all points")) if p1s == P1STAR \
            else ((1.0e-3, "new + Sec. 4 up to 1e-3 (p_1* at the edge of its error)"),)
        for dmax, lab in ranges:
            pts = [p for p in allp if p[0] - P1STAR <= dmax + 1e-9]
            if len(pts) >= 5:
                fit_set(pts, p1s, f"## {lab}", show_profile=(p1s == P1STAR))
    for dmax, lab in ((3.0e-4, "new points only"), (1.0e-3, "new + Sec. 4 up to 1e-3")):
        free_p1star(allp, dmax, f"## {lab}")


def cmd_slopes(args):
    path, all_status, L = fit_args(args)
    new = sorted(read_points(path, L, all_status))
    old = old_points()
    allp = sorted(new + old)
    p1 = np.array([p[0] for p in new]); rho = np.array([p[1] for p in new]); err = np.array([p[2] for p in new])
    dp = p1 - P1STAR
    pts = [p for p in allp if p[0] - P1STAR <= 1e-3 + 1e-9]
    fits, *_ = fit_set(pts, P1STAR, "## fits used for the slope comparison (new + Sec. 4 up to 1e-3)", show_profile=False)
    print("\n| segment dp_1 | drho | err | drho/d(dp_1) per 1e-4 | drho/d(ln dp_1) | log fit | pow0 fit | pow fit | jump0.3 fit |")
    print("|---|---|---|---|---|---|---|---|---|")
    for i in range(len(dp) - 1):
        a, b = dp[i], dp[i + 1]
        d = rho[i + 1] - rho[i]; e = np.hypot(err[i], err[i + 1])
        pred = " | ".join(f"{(F['f'](b) - F['f'](a)) / (b - a) * 1e-4:+.4f}" for F in (fits[0], fits[1], fits[2], fits[3]))
        print(f"| {a:.1e} -> {b:.1e} | {d:+.4f} | {e:.4f} | {d / (b - a) * 1e-4:+.4f} | {d / np.log(b / a):+.4f} | {pred} |")
    print("(the fit columns are the same finite differences of each fitted curve, per 1e-4 in dp_1)")


def cmd_blocks(args):
    """Means over consecutive 1e5-step blocks for every run with dp_1 <= dmax (default 0.5e-4): the
    seed-by-seed history near the stability edge."""
    dmax = float(args[0]) if args else 0.5e-4
    rs = [r for r in load_runs() if r["p1"] - P1STAR <= dmax + 1e-9]
    nb = int(max(r["t_last"] for r in rs) // 1e5)
    print("| run | status | " + " | ".join(f"{i:g}-{i + 1:g}e5" for i in range(nb)) + " |")
    print("|" + "---|" * (nb + 2))
    for r in sorted(rs, key=lambda r: (r["p1"], r["L"], r["seed"])):
        t, rho = series(os.path.join("bulk_dens", r["file"]))
        st = "absorbed at %.3g" % r["t_end"] if r["absorbed"] else ("done" if r["done"] else "running")
        cells = []
        for i in range(nb):
            tb, rb = uniform(t, rho, i * 1e5, (i + 1) * 1e5)
            cells.append(f"{rb.mean():.3f}" if len(rb) >= 8 and t[-1] >= (i + 1) * 1e5 else "-")
        print(f"| {r['file'][:-4]} | {st} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "all":
        cmd_runs([]); print(); cmd_points([]); print(); cmd_fit([]); print(); cmd_slopes([])
    else:
        dict(runs=cmd_runs, points=cmd_points, fit=cmd_fit, slopes=cmd_slopes, blocks=cmd_blocks)[cmd](args)
