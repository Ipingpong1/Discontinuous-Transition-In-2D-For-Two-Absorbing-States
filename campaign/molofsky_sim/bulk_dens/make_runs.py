#!/usr/bin/env python3
"""Run list for the seed-averaged density campaign (NOTES.md Sec. 9): the stationary interface
density of the active phase of Molofsky's rule above the transition line on the cut p2 = 0.45
(front reversal p1* = 0.18535(3)).  Writes bulk_dens/runs.txt in the order the chain
(run_dens.sh, 12 niced single-thread jobs) should start them:

  (ii)  L 4096, seeds 11-12, T 3e5 at dp1 = 0.5e-4 and 1e-4 (longest jobs first)
  (i)   L 2048, seeds 11-18, T 6e5 at dp1 = 0.4, 0.5, 0.7, 1, 1.5, 2, 3e-4
  (iii) L 2048, seeds 11-14, T 6e5 at dp1 = 0.3e-4 (0.18538), interleaved with (i)

in seed pairs (11-12 for every p1, then 13-14, ...) so that partial results cover every p1.
Seeds 11-18 are new (Secs. 4 and 8 used seeds 1-4), so no run duplicates an earlier one.
Every run checkpoints every 5e4 steps and resumes from its checkpoint if restarted; runs that
already carry a '# done' line are left out, so re-running the chain only finishes what is missing.
Invocation as in Sec. 8 (random initial condition at density 1/2, --stop-absorbed, outputs every
2000 steps plus 10 per decade).

NEVER rewrite a run list while a chain is reading it: xargs reads its input lazily (8 kB at a
time), so a shorter file written underneath it drops runs from the stream and can hand xargs a
truncated line.  That happened at 17:00 on 30 Oct 2026 (make_runs.py re-run without arguments to
check the budget, when it still defaulted to --main): the main chain skipped 16 runs (seeds 17-18
at every p1, seed 16 at 0.18555 and 0.18565), and one command arrived cut off before its --out
(0.18555 seed 16): that run went to completion with its series on /dev/null and only its checkpoint
saved, which --resume then took for a finished run.  The 15 skipped runs were run afterwards from
bulk_dens/runs_missing.txt (run_missing.sh); the orphan checkpoint was deleted and that run redone
(run_one.sh).  Every other output's header was checked against its file name (p1, p2, L, seed,
init, nbhd).  Without arguments the script is now a dry run that writes nothing.

Extensions (decided after the first two seeds, see NOTES.md Sec. 9):
  python3 make_runs.py --missing  bulk_dens/runs_missing.txt: main runs that were never started
                                  (run_missing.sh)
  python3 make_runs.py --big   bulk_dens/runs_big.txt: L 4096, seeds 11-12, T 6e5 at dp1 = 0.4e-4
                               (started while the main chain was running, run_big.sh)
  python3 make_runs.py --ext   bulk_dens/runs_ext.txt: the finished runs at dp1 = 0.3, 0.4, 0.5e-4
                               (L 2048) continued in place to 1.2e6 steps, and the L 4096 runs at
                               0.5e-4 and 1e-4 to 6e5, all through --resume (run_ext.sh, after the
                               main chain)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)

P2 = "0.45"
P1_MAIN = ["0.18539", "0.18540", "0.18542", "0.18545", "0.18550", "0.18555", "0.18565"]  # dp1 = 0.4 .. 3e-4
P1_EDGE = "0.18538"                                                                      # dp1 = 0.3e-4
SEEDS_MAIN = range(11, 19)
SEEDS_EDGE = range(11, 15)
SEEDS_BIG = range(11, 13)


def done(path):
    return os.path.exists(path) and any(l.startswith("# done") for l in open(path))


def bulk(p1, L, s, T, every=2000):
    name = f"p2_{P2}_p1_{p1}_L{L}_s{s}"
    out = f"bulk_dens/{name}.txt"
    cmd = (f"--molofsky {p1},{P2} --nbhd vn --L {L} --T {T} --init random:0.5 --seed {s} --ppd 10 "
           f"--every {every} --stop-absorbed --ckpt bulk_dens/ckpt/{name}.ck --ckpt-every 50000 --resume --out {out}")
    return cmd, out


def done_to(path, T):
    """Finished with t_end >= T (a run continued in place carries one '# done' line, at its end)."""
    if not os.path.exists(path):
        return False
    for l in open(path):
        if l.startswith("# done") and int(l.split()[3]) >= T:
            return True
    return False


mode = sys.argv[1] if len(sys.argv) > 1 else "--dry"
if mode in ("--main", "--dry", "--missing"):
    runs = []
    for s in SEEDS_BIG:
        for p1 in ("0.18540", "0.18545"):
            runs.append(bulk(p1, 4096, s, 300000))
    for s0 in range(11, 19, 2):
        for s in (s0, s0 + 1):
            for p1 in P1_MAIN:
                runs.append(bulk(p1, 2048, s, 600000))
            if s in SEEDS_EDGE:
                runs.append(bulk(P1_EDGE, 2048, s, 600000))
    todo = [c for c, f in runs if not done(f)]
    if mode == "--main":
        out = "bulk_dens/runs.txt"
        with open("bulk_dens/runs_all.txt", "w") as f:
            f.write("".join(c + "\n" for c, _ in runs))
    elif mode == "--missing":
        # main runs the chain never started (not finished and not being continued by the extension chain)
        ext = open("bulk_dens/runs_ext.txt").read() if os.path.exists("bulk_dens/runs_ext.txt") else ""
        todo = [c for c, f in runs if not done(f) and os.path.basename(f) not in ext]
        out = "bulk_dens/runs_missing.txt"
    else:
        out = None
elif mode == "--big":
    runs = [bulk("0.18539", 4096, s, 600000) for s in SEEDS_BIG]
    todo = [c for c, f in runs if not done_to(f, 600000)]
    out = "bulk_dens/runs_big.txt"
elif mode == "--ext":
    runs = []
    for s in SEEDS_BIG:
        for p1 in ("0.18540", "0.18545"):
            runs.append(bulk(p1, 4096, s, 600000))
    for s in SEEDS_MAIN:
        for p1 in ("0.18539", "0.18540"):
            runs.append(bulk(p1, 2048, s, 1200000))
        if s in SEEDS_EDGE:
            runs.append(bulk(P1_EDGE, 2048, s, 1200000))
    # only runs whose first segment is finished are continued; the rest are listed in the log
    todo = [c for c, f in runs if done(f) and not done_to(f, int(c.split("--T ")[1].split()[0]))]
    skipped = [f for c, f in runs if not done(f)]
    if skipped:
        print("not finished, not extended:", " ".join(skipped))
    out = "bulk_dens/runs_ext.txt"
else:
    sys.exit("usage: make_runs.py [--dry | --main | --missing | --big | --ext]")

if out is not None:
    with open(out, "w") as f:
        f.write("".join(c + "\n" for c in todo))
else:
    out = "(dry run, nothing written)"


def work(c):
    L = int(c.split("--L ")[1].split()[0]); T = int(c.split("--T ")[1].split()[0])
    f = c.split("--out ")[1].strip()
    t0 = 0
    if os.path.exists(f):   # continued runs: only the remaining steps cost
        for l in open(f):
            if l.startswith("# done"):
                t0 = max(t0, int(l.split()[3]))
    return L * L * max(T - t0, 0)


cells = sum(work(c) for c in todo)
print(f"{out}: {len(todo)}/{len(runs)} runs to do, {cells:.3g} cell updates "
      f"(~{cells / 1.1e9 / 3600:.0f} job-hours at 1.1e9/s per job under full load)")
