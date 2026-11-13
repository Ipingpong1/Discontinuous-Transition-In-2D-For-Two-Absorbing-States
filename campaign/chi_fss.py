#!/usr/bin/env python3
"""Magnetisation susceptibility near q* from long bulk runs (E6; E7 uses the same estimator).

    python3 chi_fss.py RUNS_DIR      (RUNS_DIR/E6_long*/out/*.txt)

Per run with 1.005 < q/q* < 1.06 that is still active at the end, over t >= 2e5 (records every 1000):
  chi   = L^2 <m^2>, m = alive - 1/2;
  U4    = 1 - <m^4> / (3 <m^2>^2);
  tau   = integrated autocorrelation time of m (window to the first zero of the ACF);
  n_eff = n / (2 tau / 1000); the chi error printed is chi / sqrt(n_eff).
"""
import glob
import os
import re
import sys

import numpy as np

QSTAR = {0.98: 0.00252, 0.99: 0.00127, 0.995: 0.00062}


def load(f):
    return np.array([[float(x) for x in l.split()] for l in open(f) if not l.startswith("#")])


def main():
    print(f"{'p':6}{'q/q*':7}{'L':6}{'chi':>8}{'±':>7}{'U4':>7}{'<m>':>8}{'tau_int':>9}{'n_eff':>7}")
    for f in sorted(glob.glob(os.path.join(sys.argv[1], "E6_long*", "out", "*.txt"))):
        n = os.path.basename(f)
        p = float(re.search(r"^p([\d.]+)_", n)[1])
        q = float(re.search(r"_q([\d.]+)_", n)[1])
        L = int(re.search(r"_L(\d+)_", n)[1])
        if p not in QSTAR or not 1.005 < q / QSTAR[p] < 1.06:
            continue
        ts = load(f)
        late = ts[(ts[:, 0] >= 2e5) & (ts[:, 0] % 1000 == 0)]
        if late[-1, 2] < 0.05:
            continue
        m = late[:, 1] - 0.5
        mc = m - m.mean()
        ac = np.correlate(mc, mc, "full")[len(mc) - 1:]
        ac /= ac[0]
        k0 = np.argmax(ac < 0) or len(ac)
        tint = 0.5 + ac[1:k0].sum()
        neff = len(m) / (2 * tint)
        chi = L * L * (m ** 2).mean()
        u4 = 1 - (m ** 4).mean() / (3 * (m ** 2).mean() ** 2)
        print(f"{p:<6}{q / QSTAR[p]:<7.3f}{L:<6}{chi:8.0f}{chi / np.sqrt(neff):7.0f}{u4:7.2f}{m.mean():+8.3f}"
              f"{tint * 1000:9.0f}{neff:7.1f}")


if __name__ == "__main__":
    main()
