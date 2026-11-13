#!/usr/bin/env python3
"""Two-point correlation function of the active phase at p = 1/2 from bulk snapshots (R2_corr, L 2048,
t >= 1e5): m = 2 n - 1, C(r) = <m_i m_{i+r}> - <m>^2 by FFT, averaged over the two axes, snapshots and
seeds.  Reports chi = sum_r C(r) = L^2 <m^2>_c, the exponential length from ln C(r) over 4 <= r <= 40,
and the paper's estimate (chi / 2 pi)^{1/2}.

    python3 corr.py DIR [--tmin 1e5]"""
import argparse
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import snap  # noqa: E402

QS = 0.04596


def corr_axes(g):
    m = 2.0 * g.astype(np.float32) - 1.0
    mbar = m.mean()
    F = np.fft.rfft2(m)
    P = np.fft.irfft2(F * np.conj(F), s=m.shape) / m.size   # <m_i m_{i+r}> (periodic)
    C = P - mbar ** 2
    L = m.shape[0]
    axis = 0.5 * (C[0, :L // 2] + C[:L // 2, 0])
    return axis, C, mbar


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--tmin", type=float, default=1e5)
    a = ap.parse_args()
    # chi_L = L^2 <m^2> from the time series (alive_frac every 1000 steps, t >= tmin): the periodic sum of
    # the connected function over one snapshot is identically zero, so chi is an ensemble average
    acc = defaultdict(lambda: [0, None, [], 0])
    for f in sorted(glob.glob(os.path.join(a.dir, "*.txt"))):
        pr = f[:-4]
        m = re.search(r"_q([\d.]+)_L(\d+)_s(\d+)$", os.path.basename(pr))
        if not m or not os.path.exists(pr + ".snap"):
            continue
        q, L = float(m[1]), int(m[2])
        ts = np.array([[float(v) for v in l.split()] for l in open(f) if l.strip() and not l.startswith("#")])
        sel = ts[:, 0] >= a.tmin
        e = acc[q]
        e[2].extend(((2 * ts[sel, 1] - 1) ** 2 * L * L).tolist())
        e[3] += 1
        for t, g in snap.read_snaps(pr, a.tmin):
            axis, C, mbar = corr_axes(g)
            e[0] += 1
            e[1] = axis if e[1] is None else e[1] + axis
    print("q         x       runs snaps  chi_L=L^2<m^2> (n_t)   xi_exp(4..40)  xi_exp(10..60)  (chi/2pi)^1/2   C(1)    C(5)    C(10)   C(20)   C(40)   C(60)")
    for q, (n, axis, chis, nr) in sorted(acc.items()):
        axis = axis / n
        r = np.arange(len(axis))
        def xi_fit(r1, r2):
            s = (r >= r1) & (r <= r2) & (axis > 0)
            return -1 / np.polyfit(r[s], np.log(axis[s]), 1)[0]
        chi = np.mean(chis)
        chi_e = np.std(chis, ddof=1) / np.sqrt(len(chis) / 10)    # ~10 samples per correlation time (guess)
        print(f"{q:<9} {q / QS:<7.3f} {nr:<4d} {n:<5d} {chi:8.1f} +- {chi_e:5.1f} ({len(chis)})  {xi_fit(4, 40):8.2f}       {xi_fit(10, 60):8.2f}      "
              f"{np.sqrt(chi / (2 * np.pi)):8.2f}     " + " ".join(f"{axis[k]:.4f}" for k in (1, 5, 10, 20, 40, 60)))


if __name__ == "__main__":
    main()
