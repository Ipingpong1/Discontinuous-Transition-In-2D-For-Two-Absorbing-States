#!/usr/bin/env python3
"""Cost of the campaign from the '# done' lines of every output (wall_s per run, single thread).

    python3 cost.py
"""
import glob
import os
import re

D = os.path.dirname(os.path.abspath(__file__))
tot, n, per = 0.0, 0, {}
for txt in sorted(glob.glob(f"{D}/*/out/*.txt")):
    sub = txt.split("/")[-3]
    for line in open(txt):
        if line.startswith("# done"):
            w = float(re.search(r"wall_s ([0-9.]+)", line)[1])
            tot += w; n += 1
            per[sub] = per.get(sub, 0.0) + w
for sub, w in sorted(per.items()):
    print(f"{sub:8s} {w / 3600:7.1f} job-h")
print(f"total {n} finished runs, {tot / 3600:.1f} job-hours (sum of single-thread wall times under 12-15-job load)")
