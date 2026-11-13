#!/usr/bin/env python3
"""Cost of the Sec. 9 runs in job-hours: the sum of the wall times of every run segment.  A run that
was continued in place (--resume with a larger --T) carries in its output only the done line of the
last segment; the done lines of the first segments of the 48 runs the main chain finished are kept
in bulk_dens/walltimes_main.txt, so the total is  sum(walltimes_main) + sum(done lines in the
outputs that are not the same first segment again).  Also prints the wall-clock span of the chains
from bulk_dens/run_dens.log.  Usage: python3 bulk_dens/cost.py (from campaign/molofsky_sim)."""
import glob
import os
import re
import subprocess

WALL = re.compile(r"wall_s ([0-9.]+)")
first = {}                                   # file -> wall_s of the first segment (main chain)
for line in open("bulk_dens/walltimes_main.txt"):
    f, rest = line.split(":", 1)
    first[os.path.basename(f)] = float(WALL.search(rest).group(1))
tot_first = sum(first.values())
tot_other = 0.0; n_other = 0; n_dup = 0; rows = []
for f in sorted(glob.glob("bulk_dens/p2_*.txt")):
    done = [l for l in open(f) if l.startswith("# done")]
    if not done:
        rows.append((os.path.basename(f), "running", 0.0)); continue
    w = float(WALL.search(done[-1]).group(1))
    b = os.path.basename(f)
    if b in first and abs(first[b] - w) < 1e-6:
        n_dup += 1; rows.append((b, "main first segment", w))
    else:
        tot_other += w; n_other += 1
        rows.append((b, "continued segment" if b in first else "later chain", w))
print(f"main chain first segments (walltimes_main.txt): {len(first)} runs, {tot_first / 3600:.1f} job-hours")
print(f"other segments (missing runs, redo, L 4096 side runs, continuations): {n_other} done lines, {tot_other / 3600:.1f} job-hours")
print(f"({n_dup} outputs repeat a main first segment; {sum(1 for r in rows if r[1] == 'running')} runs still running)")
print(f"TOTAL so far: {(tot_first + tot_other) / 3600:.1f} job-hours")
log = open("bulk_dens/run_dens.log").read().splitlines()
marks = [l for l in log if re.match(r"^(START|BIG_START|DONE|EXT_START|EXT_DONE|MISS_START|MISS_DONE|ONE_START|ONE_DONE|BIG_DONE|EXT2_START|EXT2_DONE)", l)]
for m in marks:
    print("  ", m[:60])
