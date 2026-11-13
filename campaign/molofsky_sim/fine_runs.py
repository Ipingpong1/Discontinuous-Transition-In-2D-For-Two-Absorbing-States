#!/usr/bin/env python3
"""Run lists for the fine scan (NOTES.md Sec. 8): strips at 1e-4 steps across the front
reversal (scan_fine/), bulk and lifetime runs at the same p1 (bulk_fine/), and the density
points above the line at p2 = 0.45.  Writes scan_fine/runs.txt, bulk_fine/runs.txt and the
ordered combined list runs_fine.txt (phase order: strips s1-2, lifetimes, small bulks, strips
s3-4, L 2048 bulks, L 1024 1e6 density runs) plus runs_fine_wide.txt (Lx 4096 strips at the
two points nearest each crossing).  Finished runs (a '# done' line) are left out.
"""
import os

GRID = {"0.45": ["0.18505", "0.18515", "0.18525", "0.18535", "0.18545", "0.18555"],
        "0.50": ["0.17027", "0.17037", "0.17047", "0.17057", "0.17067", "0.17077"]}
DENS = ["0.18540", "0.18545", "0.18550", "0.18565"]          # p1* + 0.5e-4, 1e-4, 1.5e-4, 3e-4 at p2 = 0.45


def done(path):
    return os.path.exists(path) and any(l.startswith("# done") for l in open(path))


def strip(p2, p1, s, Lx, out):
    f = f"{out}/p2_{p2}_p1_{p1}_s{s}"
    return (f"--molofsky {p1},{p2} --Lx {Lx} --Ly 4096 --T 100000 --init strip:0.5:0.5 --bg 0 "
            f"--stop-absorbed --seed {s} --ppd 10 --every 2000 --out {f}.txt --profile {f}.prof"), f + ".txt"


def bulk(p2, p1, L, s, T, every, tag=""):
    f = f"bulk_fine/p2_{p2}_p1_{p1}_L{L}{tag}_s{s}.txt"
    return (f"--molofsky {p1},{p2} --L {L} --T {T} --init random:0.5 --seed {s} --ppd 10 --every {every} "
            f"--stop-absorbed --out {f}"), f


def life(p2, p1, L, s):
    f = f"bulk_fine/life_p2_{p2}_p1_{p1}_L{L}_s{s}.txt"
    return f"--molofsky {p1},{p2} --L {L} --T 1000000 --init random:0.5 --seed {s} --ppd 10 --stop-absorbed --out {f}", f


phases = {k: [] for k in "ABCDEFG"}
for s in (1, 2):
    for p2, p1s in GRID.items():
        for p1 in p1s:
            phases["A"].append(strip(p2, p1, s, 2048, "scan_fine"))
for p2, p1s in GRID.items():
    for p1 in p1s:
        for L in (256, 512):
            for s in (1, 2, 3, 4):
                phases["B"].append(life(p2, p1, L, s))
for p2, p1s in GRID.items():
    for p1 in p1s:
        for L in (512, 1024):
            for s in (1, 2):
                phases["C"].append(bulk(p2, p1, L, s, 300000, 10000))
for s in (3, 4):
    for p2, p1s in GRID.items():
        for p1 in p1s:
            phases["D"].append(strip(p2, p1, s, 2048, "scan_fine"))
for p2, p1s in GRID.items():
    for p1 in p1s:
        phases["E"].append(bulk(p2, p1, 2048, 1, 300000, 10000))
for p1 in DENS:
    if p1 not in GRID["0.45"]:
        phases["E"].append(bulk("0.45", p1, 2048, 1, 300000, 10000))
for p1 in DENS:
    for s in (1, 2):
        phases["F"].append(bulk("0.45", p1, 1024, s, 1000000, 20000, tag="_T1e6"))
for s in (1, 2):
    for p2, p1s in {"0.45": ["0.18525", "0.18535"], "0.50": ["0.17047", "0.17057"]}.items():
        for p1 in p1s:
            phases["G"].append(strip(p2, p1, s, 4096, "scan_fine/wide"))

os.makedirs("scan_fine/wide", exist_ok=True); os.makedirs("bulk_fine", exist_ok=True)
with open("runs_fine.txt", "w") as main, open("runs_fine_wide.txt", "w") as wide:
    for k, runs in phases.items():
        todo = [c for c, f in runs if not done(f)]
        (wide if k == "G" else main).write("".join(c + "\n" for c in todo))
        print(f"phase {k}: {len(todo)}/{len(runs)} to run")
with open("scan_fine/runs.txt", "w") as f:
    f.write("".join(c + "\n" for k in "ADG" for c, _ in phases[k]))
with open("bulk_fine/runs.txt", "w") as f:
    f.write("".join(c + "\n" for k in "BCEF" for c, _ in phases[k]))
