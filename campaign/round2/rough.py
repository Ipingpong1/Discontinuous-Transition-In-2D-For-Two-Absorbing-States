#!/usr/bin/env python3
"""Roughness of the edges of the wetting layer (R2_wall snapshots) or of a (1,0) band's fronts
(R2_orient snapshots).

    python3 rough.py RUN_PREFIX [...] [--tmin 1e5] [--csv OUT.csv]

Walls (init strip:1:0.5): the alive phase is rows [0, Ly/2), the dead phase above; the two walls are
at y = Ly/2 and y = 0.  For each wall and snapshot, per column x: the lowest and highest interface
rows within +-Ly/4 of the wall (the two edges of the active layer), their difference w(x) (the local
layer width).  Per edge: sigma^2 = Var_x(h), the height-difference correlation C(r) = <(h(x+r)-h(x))^2>
on a log grid of r, and the local exponent 2 zeta = dln C/dln r.  Averaged over snapshots (t >= tmin)
and over the two walls.  Bands (init band, (1,0)): the same with the band's two edges at y = 0 and
y = L/2 (half-torus windows of L/4), i.e. two free active-absorbing fronts.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import snap  # noqa: E402


def analyse(prefix, tmin, rs):
    Lx, Ly = snap.lattice_size(prefix)
    acc = dict(n=0, w=0.0, w2=0.0, var_lo=0.0, var_hi=0.0, C_lo=np.zeros(len(rs)), C_hi=np.zeros(len(rs)),
               whist=np.zeros(4096), I=0.0, nan_cols=0.0)
    for t, g in snap.read_snaps(prefix, tmin):
        m = snap.iface_mask(g)
        for yc in (Ly // 2, 0):
            lo, hi = snap.column_edges(m, yc, Ly // 4)
            w = hi - lo + 1
            ok = np.isfinite(w)
            acc["nan_cols"] += (~ok).mean()
            acc["n"] += 1
            acc["w"] += np.nanmean(w)
            acc["w2"] += np.nanmean(w * w)
            acc["var_lo"] += np.nanvar(lo)
            acc["var_hi"] += np.nanvar(hi)
            acc["C_lo"] += snap.height_diff(lo, rs)
            acc["C_hi"] += snap.height_diff(hi, rs)
            wi = np.clip(w[ok].astype(int), 0, 4095)
            acc["whist"] += np.bincount(wi, minlength=4096)
            rows = (np.arange(yc - Ly // 4, yc + Ly // 4)) % Ly
            acc["I"] += m[rows].sum() / Lx
    n = max(acc["n"], 1)
    return dict(prefix=prefix, Lx=Lx, Ly=Ly, nsnap=acc["n"] // 2, w=acc["w"] / n,
                w_sd=np.sqrt(max(acc["w2"] / n - (acc["w"] / n) ** 2, 0)), sig_lo=np.sqrt(acc["var_lo"] / n),
                sig_hi=np.sqrt(acc["var_hi"] / n), C_lo=acc["C_lo"] / n, C_hi=acc["C_hi"] / n,
                whist=acc["whist"] / n, I=acc["I"] / n, nan_cols=acc["nan_cols"] / n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefixes", nargs="+")
    ap.add_argument("--tmin", type=float, default=1e5)
    ap.add_argument("--csv")
    ap.add_argument("--hist", help="CSV of the per-column layer-width distribution P(w) (one file, all runs)")
    a = ap.parse_args()
    rows = []
    for pr in a.prefixes:
        pr = pr[:-5] if pr.endswith(".snap") else pr
        Lx, Ly = snap.lattice_size(pr)
        rs = np.unique(np.round(np.logspace(0, np.log10(Lx / 2), 25)).astype(int))
        r = analyse(pr, a.tmin, rs)
        C = 0.5 * (r["C_lo"] + r["C_hi"])
        lc, lr = np.log(C), np.log(rs)
        zeta = np.gradient(lc, lr) / 2
        print(f"== {os.path.basename(pr)}  Lx {Lx} Ly {Ly}  snapshots {r['nsnap']} (t >= {a.tmin:g})")
        print(f"   layer width per column: mean {r['w']:.1f}  sd {r['w_sd']:.1f}   edges sigma: lower {r['sig_lo']:.2f} "
              f"upper {r['sig_hi']:.2f}   I {r['I']:.2f}   columns without interface {100 * r['nan_cols']:.2f}%")
        print("   r      C(r)      zeta_loc")
        for i in range(len(rs)):
            print(f"   {rs[i]:<6d} {C[i]:9.2f}  {zeta[i]:6.3f}")
        rows.append((pr, Lx, Ly, r, rs, C))
    if a.csv:
        with open(a.csv, "w") as f:
            f.write("prefix,Lx,Ly,nsnap,w_mean,w_sd,sig_lo,sig_hi,I,r,C_lo,C_hi\n")
            for pr, Lx, Ly, r, rs, C in rows:
                for i, ri in enumerate(rs):
                    f.write(f"{os.path.basename(pr)},{Lx},{Ly},{r['nsnap']},{r['w']:.3f},{r['w_sd']:.3f},"
                            f"{r['sig_lo']:.4f},{r['sig_hi']:.4f},{r['I']:.4f},{ri},{r['C_lo'][i]:.5f},{r['C_hi'][i]:.5f}\n")
    if a.hist:
        with open(a.hist, "w") as f:
            f.write("prefix,w,P\n")
            for pr, Lx, Ly, r, rs, C in rows:
                h = r["whist"] / max(r["whist"].sum(), 1)
                nz = np.flatnonzero(h > 0)
                for wi in range(nz.min(), nz.max() + 1):
                    f.write(f"{os.path.basename(pr)},{wi},{h[wi]:.6f}\n")


if __name__ == "__main__":
    main()
