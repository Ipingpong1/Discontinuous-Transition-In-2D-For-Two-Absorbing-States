#!/usr/bin/env python3
"""Run lists for NOTES.md (strip crossings at p = 0.99, 0.995 beyond the crossover time).

    python3 make_runs.py            # dry run: prints the counts and the cost estimate
    python3 make_runs.py write      # writes runs_wide.txt, runs_inner.txt, runs_outer.txt, runs_bulk.txt, runs_extra.txt
    python3 make_runs.py write runs_extra.txt   # writes only that list (the others are consumed by a running chain)

One dnsim argument line per run (absolute paths); run_all.sh feeds them to xargs.  Names follow
E4_voter / E6_long so that campaign/fronts.py and vwin.py read the folders unchanged.
"""
import os
import sys

D = os.path.dirname(os.path.abspath(__file__))
BIN_RATE = 2.0e9  # cell updates / s, single thread, idle machine (test/rate_p995_T3000.txt)

STRIP = ("--rule B3678/S34678 --alpha 0.5 --voter {p} --jflip 2:{q} --Lx {Lx} --Ly {Ly} --T {T} "
         "--init strip:0.5:{frac} --bg 0 --every 2000 --ppd 10 --stop-absorbed --seed {s} "
         "--out {P}.txt --profile {P}.prof --ckpt {P}.ckpt --ckpt-every 50000 --resume")
BULK = ("--rule B3678/S34678 --alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --T {T} "
        "--init random:0.5 --every 1000 --ppd 20 --stop-absorbed --seed {s} "
        "--out {P}.txt --ckpt {P}.ckpt --ckpt-every 50000 --resume")

T_STRIP = {0.995: 1600000, 0.99: 1200000}
T_BULK = {0.995: 1500000, 0.99: 1200000}
INNER = {0.995: ["0.00060", "0.00062", "0.00064"], 0.99: ["0.00125", "0.00128", "0.00131"]}
OUTER = {0.995: ["0.00058", "0.00066"], 0.99: ["0.00122", "0.00134"]}
BULKQ = {0.995: ["0.00063", "0.000636", "0.00064"], 0.99: ["0.00128", "0.00130"]}


def strip(p, q, Lx, Ly, s, sub, frac=0.5, T=None):
    P = f"{D}/{sub}/out/p{p}_a0.5_j2_q{q}_Lx{Lx}_Ly{Ly}_s{s}"
    T = T or T_STRIP[p]
    return STRIP.format(p=p, q=q, Lx=Lx, Ly=Ly, T=T, s=s, P=P, frac=frac), Lx * Ly * T


def bulk(p, q, L, s):
    P = f"{D}/bulk/out/p{p}_a0.5_j2_q{q}_L{L}_s{s}"
    return BULK.format(p=p, q=q, L=L, T=T_BULK[p], s=s, P=P), L * L * T_BULK[p]


def lists():
    wide = [strip(0.995, q, 4096, 4096, s, "wide") for q in ("0.00062", "0.00064") for s in (1, 2)]
    inner = [strip(p, q, 2048, 4096, s, "strips") for p in (0.995, 0.99) for q in INNER[p] for s in (1, 2, 3, 4)]
    outer = [strip(p, q, 2048, 4096, s, "strips") for p in (0.995, 0.99) for q in OUTER[p] for s in (1, 2, 3, 4)]
    blk = ([bulk(p, q, 2048, s) for p in (0.995, 0.99) for q in BULKQ[p] for s in (3, 4)]
           + [bulk(0.995, "0.00064", 4096, 3)])
    # extra (added 10-31 00:40 after the first wave showed the p 0.995 crossing moving above 0.00064):
    # one more q above each grid, 4 seeds, and bulk seeds at 0.00066
    extra = ([strip(0.995, "0.00068", 2048, 4096, s, "strips") for s in (1, 2, 3, 4)]
             + [strip(0.99, "0.00137", 2048, 4096, s, "strips") for s in (1, 2, 3, 4)]
             + [bulk(0.995, "0.00066", 2048, s) for s in (3, 4)])
    # fill (added 10-31 02:10, run by run_fill.sh in the slots the other chains free): a taller box for
    # the slope above q*, Lx 1024 x Ly 8192 with the strip = Ly/4 (room for 2vT ~ 4900 rows of growth
    # against 1434 in the Ly 4096 box, whose strips at 0.00068 filled it before t_x); T 1.6e6 at the two
    # q nearest the crossing, 1.2e6 above
    fill = [strip(0.995, q, 1024, 8192, s, "tall", 0.25, T) for q, T in
            (("0.00064", 1600000), ("0.00066", 1600000), ("0.00068", 1200000), ("0.00070", 1200000)) for s in (1, 2)]
    # fill2 (added 10-31 08:45): the same tall box at p 0.99 for the slope above q*(0.99) (the Ly 4096
    # strips at 0.00131-0.00134 fill the box by 5-8e5; the queued Ly 4096 runs at 0.00137 were skipped
    # by lock files, see run_all.log)
    fill2 = [strip(0.99, q, 1024, 8192, s, "tall", 0.25, 1200000) for q in ("0.00128", "0.00131", "0.00134", "0.00137")
             for s in (1, 2)]
    return {"runs_wide.txt": wide, "runs_inner.txt": inner, "runs_outer.txt": outer, "runs_bulk.txt": blk,
            "runs_extra.txt": extra, "runs_fill.txt": fill, "runs_fill2.txt": fill2}


if __name__ == "__main__":
    L = lists()
    tot = 0
    for name, runs in L.items():
        ch = sum(u for _, u in runs) / BIN_RATE / 3600
        tot += ch
        print(f"{name:16s} {len(runs):3d} runs  {ch:6.1f} core-h at {BIN_RATE:.1e}/s")
        if len(sys.argv) > 1 and sys.argv[1] == "write" and (len(sys.argv) < 3 or name in sys.argv[2:]):
            for sub in ("strips", "wide", "bulk", "tall"):
                os.makedirs(os.path.join(D, sub, "out"), exist_ok=True)
            with open(os.path.join(D, name), "w") as f:
                f.write("\n".join(a for a, _ in runs) + "\n")
    print(f"total {sum(len(r) for r in L.values())} runs, {tot:.1f} core-h")
