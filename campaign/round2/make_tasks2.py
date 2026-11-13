#!/usr/bin/env python3
"""Round-2 run lists (2026-10-31): the experiments left open by the referee round, for dnsim v1.5.

    python3 make_tasks2.py --root /scratch/mo2273/dn --stage STAGE_DIR [--only R2_orient ...]

Writes STAGE_DIR/runs/<EXP>/tasks/shard_NNN.txt (lines "<out prefix>\t<dnsim args>", absolute paths
under --root, i.e. as seen on Torch) and manifest.tsv, in the format of campaign/run_shard.sh.
Experiments (q* from Table I unless noted):
  R2_orient     orientation test with direct front tracking: tilted bands (a,b) = (1,0), (2,1), (1,1)
                on L 4096 (tilt.py checkpoints, 8 seeds), q/q* = 0.99 ... 1.01, T 2e5, profile binned
                along the band normal (--pdir); snapshots for the x = 1 bands.
  R2_orient2k   the same at L 2048, p = 0.5 only, T 1e5 (size check).
  R2_drop2      two-droplet coexistence test at p = 1/2: an active disk (R 1000) in the dead state and
                a dead hole (R 1000) in the active state, q/q* = 1 -+ 0.003, 0.006, L 4096, T 2e5.
  R2_dropbulk   bulk reference (L 2048) at the R2_drop2 q values: iota(t) of the active phase.
  R2_wall       wetting at p = 1/2: Ly 2048 series (Lx 1024, 4 seeds, delta 0.05 ... 0.0065), Ly 4096 at
                delta 0.0065 and 0.004, extra seeds at Lx 4096, Lx 256 at Ly 2048, bare wall at Ly 2048;
                snapshots every 1e5 steps for the roughness of the layer edges.  Placed at the E16 q
                values (multiples of 0.0459) so that they extend the existing points.
Local (this machine, run_local.sh): R2_seed22 (2x2 live block, --reps) and R2_corr (bulk snapshots for
the correlation length at p = 1/2).
"""
import argparse
import math
import os
import zlib

MODEL = "--rule B3678/S34678"
QS = {0.1: 0.03680, 0.5: 0.04596, 0.9: 0.01194}      # Table I
QS_E16 = 0.0459                                        # the E16_wall / E14_drop05 placement value
ROOT = "/scratch/mo2273/dn"


def q_str(q, digits=7):
    return f"{q:.{digits}f}".rstrip("0").rstrip(".")


def band_ckpt(L, a, b, s):
    return f"{ROOT}/runs/R2_orient/init/band_L{L}_ab{a}-{b}_s{s}.ckpt"


ORIENT_X = (0.99, 0.995, 1.0, 1.005, 1.01)
ORIENT_AB = ((1, 0), (2, 1), (1, 1))


def r2_orient(L=4096, T=200000, seeds=8, ps=(0.1, 0.5, 0.9), snap=True):
    for p in ps:
        for a, b in ORIENT_AB:
            for x in ORIENT_X:
                q = round(QS[p] * x, 7)
                for s in range(1, seeds + 1):
                    sn = (f" --snap {{P}}.snap --snap-every {T // 5} --snap-from {T // 5}"
                          if snap and x == 1.0 else "")
                    yield (f"p{p}_a0.5_q{q_str(q)}_ab{a}-{b}_L{L}_s{s}",
                           f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --load {band_ckpt(L, a, b, s)} "
                           f"--T {T} --every 1000 --ppd 10 --stop-absorbed --profile {{P}}.prof --pdir {a},{b}{sn}",
                           L * L * T, 50000)


def r2_orient2k():
    yield from r2_orient(L=2048, T=100000, seeds=8, ps=(0.5,), snap=False)


DROP_X = (0.994, 0.997, 1.003, 1.006)


def r2_drop2(L=4096, T=200000, R=1000, seeds=4):
    for x in DROP_X:
        q = round(QS[0.5] * x, 7)
        for kind in ("disk", "hole"):
            for s in range(1, seeds + 1):
                yield (f"p0.5_a0.5_j2_q{q_str(q)}_{kind}R{R}_L{L}_s{s}",
                       f"--alpha 0.5 --voter 0.5 --jflip 2:{q} --L {L} --T {T} --init {kind}:0.5:{R} --bg 0 "
                       f"--every 1000 --ppd 10 --stop-absorbed --profile {{P}}.prof",
                       L * L * T, 50000)


def r2_dropbulk(L=2048, T=200000, seeds=2):
    for x in DROP_X:
        q = round(QS[0.5] * x, 7)
        for s in range(1, seeds + 1):
            yield (f"p0.5_a0.5_j2_q{q_str(q)}_L{L}_s{s}",
                   f"--alpha 0.5 --voter 0.5 --jflip 2:{q} --L {L} --T {T} --init random:0.5 "
                   f"--every 1000 --ppd 10 --stop-absorbed",
                   L * L * T, 50000)


def wall_args(q, Lx, Ly, T, snap_every=100000, snap_from=100000):
    jf = f"--jflip 2:{q} " if q else ""
    return (f"--alpha 0.5 --voter 0.5 {jf}--Lx {Lx} --Ly {Ly} --T {T} --init strip:1:0.5 --bg 0 "
            f"--every 10000 --ppd 20 --stop-absorbed --profile {{P}}.prof "
            f"--snap {{P}}.snap --snap-every {snap_every} --snap-from {snap_from}")


def r2_wall():
    def q_of(x):
        return round(QS_E16 * x, 7)
    runs = []
    # (a) Ly 2048 series, Lx 1024, 4 seeds
    for x in (0.95, 0.97, 0.98, 0.99, 0.995):
        for s in (1, 2, 3, 4):
            runs.append((q_of(x), 1024, 2048, s, 3000000 if x >= 0.995 else 2000000))
    # (b) Ly 4096 at the two smallest delta, 4 seeds
    for x in (0.995, 0.9975):
        for s in (1, 2, 3, 4):
            runs.append((q_of(x), 1024, 4096, s, 3000000))
    # (c) Lx 4096 x Ly 1024: new seeds 3, 4 to 2e6; seeds 1, 2 are extended from their E16
    #     checkpoints (copied into this experiment's out/ by setup_torch.sh) to 2e6
    for x in (0.98, 0.99, 0.995):
        for s in (1, 2, 3, 4):
            runs.append((q_of(x), 4096, 1024, s, 2000000))
    # (d) Lx 4096 x Ly 2048 at x 0.99 (roughness against Lx at one Ly)
    for s in (1, 2):
        runs.append((q_of(0.99), 4096, 2048, s, 2000000))
    # (e) Lx 256 x Ly 2048 at x 0.99, 0.995
    for x in (0.99, 0.995):
        for s in (1, 2, 3, 4):
            runs.append((q_of(x), 256, 2048, s, 2000000))
    # (f) bare wall at Ly 2048
    for s in (1, 2):
        runs.append((0, 1024, 2048, s, 1000000))
    for q, Lx, Ly, s, T in runs:
        yield (f"p0.5_a0.5_j2_q{q_str(q)}_Lx{Lx}_Ly{Ly}_s{s}", wall_args(q, Lx, Ly, T),
               Lx * Ly * T / 2, 200000)   # absorbing domains step 2-4x faster


EXPERIMENTS = {"R2_orient": r2_orient, "R2_orient2k": r2_orient2k, "R2_drop2": r2_drop2,
               "R2_dropbulk": r2_dropbulk, "R2_wall": r2_wall}


def main():
    global ROOT
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=ROOT)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--rate", type=float, default=1.9e9)
    ap.add_argument("--cpus", type=int, default=32)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    ROOT = a.root
    grand = 0
    for ei, (exp, gen) in enumerate(EXPERIMENTS.items()):
        if a.only and exp not in a.only:
            continue
        out_dir = f"{ROOT}/runs/{exp}/out"
        runs = []
        for name, args, cost, ckpt_every in gen():
            prefix = f"{out_dir}/{name}"
            seed = ((41 + ei) << 32) | zlib.crc32(f"{exp}/{name}".encode())
            full = (f"{MODEL} {args.replace('{P}', prefix)} --seed {seed} --out {prefix}.txt "
                    f"--ckpt {prefix}.ckpt --ckpt-every {ckpt_every} --resume")
            runs.append((cost / a.rate, prefix, full))
        core_h = sum(r[0] for r in runs) / 3600
        longest = max(r[0] for r in runs) / 3600
        nshard = max(1, math.ceil(len(runs) / a.cpus))
        shards = [[0.0, []] for _ in range(nshard)]
        for cost_s, prefix, full in sorted(runs, reverse=True):
            sh = min(shards, key=lambda s: (len(s[1]) >= a.cpus, s[0]))
            sh[0] += cost_s
            sh[1].append((prefix, full))
        print(f"{exp:12s} runs {len(runs):4d}  est {core_h:7.1f} core-h  longest {longest:4.2f} h  "
              f"shards {nshard} x {a.cpus}")
        grand += core_h
        tdir = f"{a.stage}/runs/{exp}/tasks"
        os.makedirs(tdir, exist_ok=True)
        for f in os.listdir(tdir):
            os.remove(os.path.join(tdir, f))
        for i, (_, lines) in enumerate(shards):
            with open(f"{tdir}/shard_{i:03d}.txt", "w") as f:
                for prefix, full in lines:
                    f.write(f"{prefix}\t{full}\n")
        with open(f"{a.stage}/runs/{exp}/manifest.tsv", "w") as f:
            f.write("prefix\test_core_s\targs\n")
            for cost_s, prefix, full in runs:
                f.write(f"{prefix}\t{cost_s:.0f}\t{full}\n")
    print(f"total est {grand:.1f} core-h at {a.rate:.2g} updates/s/core")


if __name__ == "__main__":
    main()
