"""Front velocities from the extract_prof.py CSVs: seed means per (p, q) and zero crossings q*(p)."""
import csv
import os
import re
from collections import defaultdict

import numpy as np

from common import FIGDATA

NAME = re.compile(r"^(?:p(?P<p>[\d.]+)_)?(?P<model>a[\d.]+|eta_a[\d.]+)(?:_r(?P<r>[\d.]+))?_j?2?_?q(?P<q>[\d.]+)_Lx(?P<Lx>\d+)_Ly(?P<Ly>\d+)_s(?P<s>\d+)$")


def load(exps):
    """{(p, a, r, q): [v per seed]} from <EXP>_fronts.csv (p = 0 for the pure D&N runs, -1 for eta)."""
    out = defaultdict(list)
    for exp in exps:
        path = os.path.join(FIGDATA, f"{exp}_fronts.csv")
        if not os.path.exists(path):
            continue
        for r in csv.DictReader(open(path)):
            n = r["run"]
            m = re.match(r"^(?:p(?P<p>[\d.]+)_)?a(?P<a>[\d.]+)(?:_r(?P<r>[\d.]+))?_j2_q(?P<q>[\d.]+)_Lx\d+_Ly\d+_s\d+$", n)
            if m:
                key = (float(m["p"] or 0), float(m["a"]), float(m["r"] or 1), float(m["q"]))
            else:
                m = re.match(r"^eta_a(?P<a>[\d.]+)_q(?P<q>[\d.]+)_Lx\d+_Ly\d+_s\d+$", n)
                if not m:
                    continue
                key = (-1.0, float(m["a"]), 1.0, float(m["q"]))
            v = float(r["v"]) if r["v"] not in ("", "nan") else np.nan
            out[key].append(v)
    return out


def curve(data, p, a=0.5, r=1.0):
    """Sorted arrays q, seed-mean v, sem, n for one (p, alpha, r)."""
    ks = sorted(k for k in data if k[0] == p and k[1] == a and k[2] == r)
    q = np.array([k[3] for k in ks])
    vs = [np.array(data[k], float) for k in ks]
    vm = np.array([np.nanmean(v) if np.isfinite(v).any() else np.nan for v in vs])
    se = np.array([np.nanstd(v, ddof=1) / np.sqrt(np.isfinite(v).sum()) if np.isfinite(v).sum() > 1 else np.nan
                   for v in vs])
    return q, vm, se, np.array([np.isfinite(v).sum() for v in vs])


def crossing(q, v, npts=2):
    """Zero of a straight line through the npts finite points nearest the sign change (q*, slope)."""
    ok = np.isfinite(v)
    q, v = q[ok], v[ok]
    if len(q) < 2:
        return np.nan, np.nan
    k = np.flatnonzero((v[:-1] < 0) & (v[1:] >= 0))
    c = k[0] + 0.5 if k.size else np.argmin(np.abs(v))
    idx = np.argsort(np.abs(np.arange(len(q)) - c))[:npts]
    a, b = np.polyfit(q[idx], v[idx], 1)
    return -b / a, a
