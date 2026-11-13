#!/usr/bin/env python3
"""Phase 3 (E6) long-run analysis: plateau drift and magnetisation fluctuations near q*.

    python3 longruns.py RUNS_DIR      (RUNS_DIR/E6_long*/out/*.txt)

Per run (bulk, T = 1e6):
  drift  = <rho>[7e5, 1e6] - <rho>[1e5, 2e5], error from 5 blocks per window.  The pure voter
           model (p = 1, q = 0, L 4096) follows 1/rho = 0.184 ln t + 2.13 and drifts -0.013
           between these windows; a stationary plateau drifts 0.
  rms m  = sqrt(<(alive - 1/2)^2>) over t >= 3e5: ~ xi / L for a finite correlation length
           xi, so it halves from L 2048 to 4096 if xi is finite.
q/q* is relative to the q* the grid was built on (E6 QSTAR_LONG / QSTAR_P995).
"abs." = absorbed or collapsing (rho < 0.05 at the end); "run." = still running.
"""
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

QSTAR = {0.98: 0.00252, 0.99: 0.00127, 0.995: 0.00062}


def load(f):
    return np.array([[float(x) for x in l.split()] for l in open(f) if not l.startswith("#")])


def key(path):
    n = os.path.basename(path)
    return (float(re.search(r"^p([\d.]+)_", n).group(1)), float(re.search(r"_q([\d.]+)_", n).group(1)),
            int(re.search(r"_L(\d+)_", n).group(1)))


def stats(ts):
    out = {"t_end": ts[-1, 0], "rho_end": ts[-1, 2]}
    a, b = ts[(ts[:, 0] >= 1e5) & (ts[:, 0] <= 2e5), 2], ts[ts[:, 0] >= 7e5, 2]
    if ts[-1, 0] >= 1e6 and len(a) >= 5 and len(b) >= 5:
        err = [np.std([x.mean() for x in np.array_split(w, 5)]) / np.sqrt(5) for w in (a, b)]
        out.update(ra=a.mean(), rb=b.mean(), drift=b.mean() - a.mean(), derr=np.hypot(*err))
    late = ts[ts[:, 0] >= 3e5]
    if len(late) >= 10 and late[-1, 2] > 0.05:
        m = late[:, 1] - 0.5
        out.update(rho_late=late[:, 2].mean(), rms_m=np.sqrt((m ** 2).mean()))
    return out


def main():
    root = sys.argv[1]
    g = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(root, "E6_long*", "out", "*.txt"))):
        g[key(f)].append(stats(load(f)))
    print("p      q          q/q*   L     drift[1-2e5 -> 7-10e5] per seed       <rho>[3e5-1e6]  rms m      t_end")
    for k in sorted(g):
        p, q, L = k
        rel = f"{q / QSTAR[p]:.3f}" if p in QSTAR else "  -  "
        rs = g[k]
        dr = "  ".join(f"{r['drift']:+.4f}±{r['derr']:.4f}" if "drift" in r else "     -       " for r in rs)
        rl = " ".join(f"{r['rho_late']:.3f}" if "rho_late" in r else "abs." if r["rho_end"] < 0.05 else "run." for r in rs)
        mm = " ".join(f"{r['rms_m']:.3f}" if "rms_m" in r else "  -  " for r in rs)
        te = " ".join(f"{r['t_end']:.0e}" for r in rs)
        print(f"{p:<6} {q:<10} {rel}  {L:<5} {dr:30s}  {rl:14s}  {mm:10s} {te}")


if __name__ == "__main__":
    main()
