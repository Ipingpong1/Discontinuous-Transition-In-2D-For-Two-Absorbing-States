#!/usr/bin/env python3
"""Front tracking for Molofsky-rule strip runs (dnsim --molofsky p1,p2 --init strip:0.5:0.5 --bg 0
--profile f.prof), and bulk summaries.

    python3 mfront.py front  out/*.txt      per run: v (rows/step per front), fit range, rho_b
    python3 mfront.py table  out/*.txt      seed-pooled v(p1) per p2 (from file names p2_X_p1_Y_sZ.txt)
    python3 mfront.py wall   out/x.txt      active width A(t) of a wall run (strip:1:0.5): growth law
    python3 mfront.py bulk   out/*.txt      rho_bond(t), alive(t) at log-spaced times

The strip occupies rows [0, Ly/2) at t = 0, so its centre is at Ly/4; roll the per-row Moore
interface density so the centre sits at Ly/2, smooth over SMOOTH rows, and walk outward from
the centre until the density drops below half the bulk value rho_b (the strip interior at the
first record with t >= T_B).  W(t) = distance between the two edges; v = dW/dt / 2 fitted over
records with t >= T_FIT, W in [0.1, 0.9] Ly and the interior still active.  A second estimate
uses the integrated width A(t) = sum_y iface_y / (Lx rho_b(t)) with rho_b(t) the running
interior density (immune to a slow drift of the bulk density).
"""
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
from profile import read  # noqa: E402

SMOOTH = 32
T_B = 1000
T_FIT = 2000
KEY = re.compile(r"p2_(?P<p2>[0-9.]+)_p1_(?P<p1>[0-9.]+)_s(?P<s>[0-9]+)")


def hdr_done(txt):
    for line in open(txt):
        if line.startswith("# done"):
            return True
    return False


def track(txt, prof=None):
    prof = prof or txt[:-4] + ".prof"
    Lx, Ly, t, alive, iface = read(txt, prof)
    d = np.roll(iface / Lx, Ly // 2 - Ly // 4, axis=1)  # strip centre -> Ly/2
    k = np.ones(SMOOTH) / SMOOTH
    h = SMOOTH // 2
    ds = np.array([np.convolve(np.concatenate([r[-h:], r, r[:h - 1]]), k, "valid") for r in d])
    c = Ly // 2
    i0 = min(np.searchsorted(t, T_B), len(t) - 1)
    rho_b = d[i0, c - Ly // 8: c + Ly // 8].mean()
    W, Wa, cen, A = [], [], [], []
    for i, row in enumerate(ds):
        rb_t = d[i, c - Ly // 16: c + Ly // 16].mean()
        for thr, out in ((0.5 * rho_b, W), (0.5 * rb_t, Wa)):
            lo = row[:c][::-1] < thr
            hi = row[c:] < thr
            yl = np.argmax(lo) if lo.any() else c
            yr = np.argmax(hi) if hi.any() else Ly - c
            out.append(yl + yr)
        cen.append(rb_t / rho_b)
        A.append(d[i].sum() / rb_t if rb_t > 0 else np.nan)
    return dict(t=t, W=np.array(W, float), Wa=np.array(Wa, float), A=np.array(A), cen=np.array(cen), Ly=Ly,
                Lx=Lx, rho_b=rho_b, alive=alive.sum(1) / (Lx * Ly))


def fit_v(r, tmin=T_FIT, cen_min=0.6):
    """v per front from the fixed-threshold width W, the adaptive-threshold width Wa (half the
    current interior density) and the integrated width A, over records with t >= tmin, the
    strip neither collapsed nor filling the box, and the interior above cen_min of its t = T_B
    density."""
    t, W, Wa, A, cen, Ly = r["t"], r["W"], r["Wa"], r["A"], r["cen"], r["Ly"]
    sel = (t >= tmin) & (W > 0.1 * Ly) & (W < 0.9 * Ly) & (cen > cen_min)
    if sel.sum() < 4:
        return np.nan, np.nan, np.nan, 0, 0
    pW = np.polyfit(t[sel], W[sel], 1)
    pWa = np.polyfit(t[sel], Wa[sel], 1)
    pA = np.polyfit(t[sel], A[sel], 1)
    return pW[0] / 2, pWa[0] / 2, pA[0] / 2, int(sel.sum()), int(t[sel][-1])


def cmd_front(files):
    print("run                                   v_edge    v_adapt   v_int      npts t_max   W_end/Ly rho_b  cen_end done")
    for f in files:
        r = track(f)
        v, va, vA, n, tm = fit_v(r)
        print(f"{os.path.basename(f)[:-4]:36s} {v:+.5f}  {va:+.5f}  {vA:+.5f}  {n:4d} {tm:6d}  {r['W'][-1] / r['Ly']:.3f}  "
              f"{r['rho_b']:.3f}  {r['cen'][-1]:.2f}  {int(hdr_done(f))}")


def cmd_table(files):
    g = defaultdict(list)
    for f in files:
        m = KEY.search(os.path.basename(f))
        if not m:
            continue
        r = track(f)
        v, va, vA, n, tm = fit_v(r)
        g[(float(m["p2"]), float(m["p1"]))].append((int(m["s"]), v, va, vA, n, tm, r["rho_b"], r["cen"][-1]))
    print("p2     p1        v_edge    err       v_adapt   err       v_int     n  v_edge seeds                         rho_b  cen_end")
    for key in sorted(g):
        rs = g[key]
        def ms(i):
            a = np.array([x[i] for x in rs]); a = a[np.isfinite(a)]
            return (a.mean() if len(a) else np.nan, a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else np.nan)
        (vm, ve), (vam, vae), (vAm, _) = ms(1), ms(2), ms(3)
        print(f"{key[0]:.3f}  {key[1]:.5f}  {vm:+.5f}  {ve:.5f}   {vam:+.5f}  {vae:.5f}   {vAm:+.5f}  {len(rs)}  "
              + " ".join(f"{x[1]:+.5f}" for x in rs).ljust(36)
              + f" {np.mean([x[6] for x in rs]):.3f}  {np.mean([x[7] for x in rs]):.2f}")


def cmd_wall(files):
    for f in files:
        Lx, Ly, t, alive, iface = read(f, f[:-4] + ".prof")
        d = iface / Lx
        print(f"# {os.path.basename(f)}: Lx {Lx} Ly {Ly}; active rows = rows with Moore-interface fraction > 0.25")
        print("# t  active_rows  sum_iface_frac  alive_frac")
        for i in range(len(t)):
            print(f"{t[i]:7d} {(d[i] > 0.25).sum():6d} {d[i].sum():9.2f} {alive[i].sum() / (Lx * Ly):.5f}")


def cmd_bulk(files):
    times = [10, 100, 1000, 3000, 10000, 30000, 100000]
    print("run".ljust(30) + " ".join(f"{t:>8d}" for t in times) + "   (rho_bond; second line alive_frac)")
    for f in files:
        ts = np.array([[float(x) for x in l.split()] for l in open(f) if not l.startswith("#")])
        tt = [t for t in times if t <= ts[-1, 0]]
        print(os.path.basename(f)[:-4].ljust(30) + " ".join(f"{np.interp(t, ts[:, 0], ts[:, 2]):8.4f}" for t in tt))
        print("".ljust(30) + " ".join(f"{np.interp(t, ts[:, 0], ts[:, 1]):8.4f}" for t in tt))


if __name__ == "__main__":
    cmd, files = sys.argv[1], sys.argv[2:]
    files = [f for pat in files for f in sorted(glob.glob(pat))] or files
    dict(front=cmd_front, table=cmd_table, wall=cmd_wall, bulk=cmd_bulk)[cmd](files)
