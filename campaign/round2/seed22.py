#!/usr/bin/env python3
"""Survival of a 2x2 block (R2_seed22, --reps): the full block (disk:1:0.75, four live sites) against
the half-filled block (disk:0.5:0.75, the E14_drop05 R = 0.75 initial condition).

    python3 seed22.py DIR
P_surv = survivors (alive or escaped at T) / replicas; error: jackknife over the 16 groups."""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, "/home/mo2273/Desktop/dn/campaign")
import spread  # noqa: E402

QS, QPL = 0.04596, 0.0459


def main():
    print("kind   q         q/q*(TabI) q/0.0459  replicas survivors  P_surv    err      P_s(t) at t = 10, 100, 1000, 10000")
    for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.txt"))):
        m = re.match(r"p0.5_a0.5_j2_q([\d.]+)_R0.75(full|half)_L1024_drop", os.path.basename(f))
        if not m or not any(l.startswith("# done") for l in open(f)):
            continue
        t, data, meta = spread.load1(f)
        n, s = data[:, -1, 0], data[:, -1, 1]
        P = s.sum() / n.sum()
        G = len(n)
        jk = np.array([(s.sum() - s[g]) / (n.sum() - n[g]) for g in range(G)])
        err = np.sqrt((G - 1) / G * ((jk - jk.mean()) ** 2).sum())
        S = data.sum(0)
        Ps = S[:, 1] / S[:, 0]
        at = [Ps[np.argmin(np.abs(t - tt))] for tt in (10, 100, 1000, 10000)]
        q = float(m[1])
        print(f"{m[2]:<6} {q:<9} {q / QS:<10.3f} {q / QPL:<9.3f} {int(n.sum()):<8d} {int(s.sum()):<9d} {P:.4f}   {err:.4f}   "
              + "  ".join(f"{v:.3f}" for v in at))


if __name__ == "__main__":
    main()
