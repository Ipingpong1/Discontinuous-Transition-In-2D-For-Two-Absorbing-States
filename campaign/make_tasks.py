#!/usr/bin/env python3
"""Generate dnsim run lists for the Day & Night symmetric-noise campaign.

Each experiment becomes runs/<EXP>/tasks/shard_NNN.txt (one run per line:
"<out prefix>\t<dnsim args>") plus runs/<EXP>/manifest.tsv.  Runs are spread over
shards by longest-processing-time packing, so each shard (one Slurm array task with
--cpus slots running single-threaded runs) has a similar estimated wall time.
Every run checkpoints and resumes, so a shard that hits its time limit is simply
resubmitted.

    python3 make_tasks.py --rate 7e8 [--root $SCRATCH/dn] [--only E1_front] [--dry]
"""
import argparse
import math
import os
import zlib

MODEL = "--rule B3678/S34678"
ROOT = None   # set by main(); E15_orient points its --load at checkpoints under it


def q_str(q, digits=5):
    """q as written in run names.  NB: E4 names round to 5 decimals while --jflip gets the
    full value (e.g. name q0.00136 ran q=0.001365); E6 names use 7 digits."""
    return f"{q:.{digits}f}".rstrip("0").rstrip(".")


def e0_map():
    """Coarse map: where does coarsening give way to the disordered phase, for each
    synchronicity alpha and each position j of the minority flip?"""
    qs = [0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018, 0.020,
          0.022, 0.025, 0.03, 0.04, 0.06, 0.1]
    for a in (0.25, 0.5, 0.75, 1.0):
        yield f"a{a}_q0_L1024", f"--alpha {a} --L 1024 --T 20000 --ppd 20 --stop-absorbed", 1024 * 1024 * 20000
        for j in (1, 2):
            for q in qs:
                yield (f"a{a}_j{j}_q{q_str(q)}_L1024",
                       f"--alpha {a} --jflip {j}:{q} --L 1024 --T 20000 --ppd 20 --stop-absorbed",
                       1024 * 1024 * 20000)


def e1_front():
    """Front velocity v(q): disordered strip against the all-dead absorbing state.
    First order: v changes sign linearly at q*, with the absorbing phase invading below.
    Continuous: v -> 0 as a power law at q_c and no retreating front below."""
    qs = [0.016, 0.017, 0.0175, 0.018, 0.0185, 0.019, 0.0195, 0.020, 0.021, 0.023, 0.026]
    for Lx, Ly, seeds in ((2048, 4096, 3), (8192, 4096, 2)):
        for q in qs:
            for s in range(1, seeds + 1):
                yield (f"a0.5_j2_q{q_str(q)}_Lx{Lx}_Ly{Ly}_s{s}",
                       f"--alpha 0.5 --jflip 2:{q} --Lx {Lx} --Ly {Ly} --T 100000 "
                       f"--init strip:0.5:0.5 --bg 0 --every 1000 --ppd 10 --stop-absorbed --profile {{P}}.prof",
                       Lx * Ly * 100000)


def e2_bulk():
    """Homogeneous random start vs system size: steady interface density on the
    disordered side (does it vanish at the transition or jump?), lifetime of the
    disordered state below q*, magnetisation statistics."""
    qs = [0.017, 0.0175, 0.018, 0.0185, 0.019, 0.0195, 0.020, 0.021, 0.023]
    for L, seeds in ((512, 4), (1024, 4), (2048, 4), (4096, 2), (8192, 2)):
        for q in qs:
            for s in range(1, seeds + 1):
                yield (f"a0.5_j2_q{q_str(q)}_L{L}_s{s}",
                       f"--alpha 0.5 --jflip 2:{q} --L {L} --T 100000 --every 1000 --ppd 20 "
                       f"--stop-absorbed",
                       L * L * 100000)


def e0_fine():
    """Finer map along the second knob (alpha), j=2 only: locates q*(alpha) and how
    the size of the jump changes, to aim E3 at a possible tricritical point."""
    qs = [round(0.003 * 1.15 ** i, 5) for i in range(25)]  # 0.003 .. 0.087
    for a in (0.125, 0.15, 0.2, 0.25, 0.3, 0.375, 0.625, 0.75):
        for q in qs:
            yield (f"a{a}_j2_q{q_str(q)}_L1024",
                   f"--alpha {a} --jflip 2:{q} --L 1024 --T 20000 --ppd 20 --stop-absorbed",
                   1024 * 1024 * 20000)


def e3_alpha():
    """Front velocity along the second knob.  q*(alpha) from the local E0_fine map
    (L=1024, 1 seed); grids span ~0.85-1.4 q*.  T scales as 1/alpha so every alpha
    gets the same number of update attempts per cell as E1 (alpha=0.5, T=1e5): a smaller
    alpha only slows the clock, and the coarse map's apparently softer jump at small
    alpha could be that alone."""
    qstar = {0.125: 0.0058, 0.25: 0.0105, 0.75: 0.037}
    grids = {a: [round(qs * f, 5) for f in (0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4)]
             for a, qs in qstar.items()}
    grids[1.0] = [0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.1]
    Lx, Ly = 2048, 4096
    for a, qs in grids.items():
        T = int(round(1e5 * 0.5 / a))
        for q in qs:
            for s in (1, 2):
                yield (f"a{a}_j2_q{q_str(q)}_Lx{Lx}_Ly{Ly}_s{s}",
                       f"--alpha {a} --jflip 2:{q} --Lx {Lx} --Ly {Ly} --T {T} "
                       f"--init strip:0.5:0.5 --bg 0 --every {T // 100} --ppd 10 --stop-absorbed "
                       f"--profile {{P}}.prof",
                       Lx * Ly * T)


# q*(p) along the D&N -> voter path (alpha=0.5, j=2), from the local strip pilot (10-14,
# Lx 512, T 1e4, one seed; +-10%).  Voter steps flip j=5 cells, which D&N never does, so
# q* first rises (0.018 -> ~0.047 at p~0.25-0.5) and then falls as ~0.13 (1-p).
QSTAR_VOTER = {0.5: 0.046, 0.75: 0.026, 0.9: 0.0122, 0.95: 0.0065, 0.98: 0.0027, 0.99: 0.0013}
VOTER_GRID = (0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.35, 1.5)


def e4_voter():
    """Where does the first-order line end?  Mix voter steps (prob p) into alpha=0.5 D&N
    with the j=2 minority flip.  p=0 is first order (E1/E2), p=1 is the voter point (GV).
    For each p: strip fronts (v(q) sign change, Lx 2048 x Ly 4096) and bulk runs (does the
    active density stay finite at q* or vanish?) on a grid around q*(p)."""
    Lx, Ly, T = 2048, 4096, 100000
    for p, qs in sorted(QSTAR_VOTER.items()):
        for q in [round(qs * f, 6) for f in VOTER_GRID]:
            for s in (1, 2):
                yield (f"p{p}_a0.5_j2_q{q_str(q)}_Lx{Lx}_Ly{Ly}_s{s}",
                       f"--alpha 0.5 --voter {p} --jflip 2:{q} --Lx {Lx} --Ly {Ly} --T {T} "
                       f"--init strip:0.5:0.5 --bg 0 --every 1000 --ppd 10 --stop-absorbed "
                       f"--profile {{P}}.prof",
                       Lx * Ly * T)


def e4_voter_bulk():
    """Bulk partner of E4_voter: random start, L 2048 (2 seeds) and 4096 (1 seed)."""
    T = 100000
    for p, qs in sorted(QSTAR_VOTER.items()):
        for q in [round(qs * f, 6) for f in VOTER_GRID]:
            for L, seeds in ((2048, 2), (4096, 1)):
                for s in range(1, seeds + 1):
                    yield (f"p{p}_a0.5_j2_q{q_str(q)}_L{L}_s{s}",
                           f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --T {T} --every 1000 "
                           f"--ppd 20 --stop-absorbed",
                           L * L * T)


# q -> R0 grid (factors of sqrt 2) spanning ~0.35-3 x the critical radius seen in the local
# droplet pilot (10-14, L=1024, 4 seeds): R_c ~ 45 (q=0.019), 23 (0.020-0.021), 11 (0.023)
DROP_GRID = {0.019: [16, 23, 32, 45, 64, 90, 128], 0.0195: [12, 16, 23, 32, 45, 64, 90],
             0.020: [8, 12, 16, 23, 32, 45, 64], 0.021: [6, 8, 12, 16, 23, 32, 45],
             0.023: [4, 6, 8, 12, 16, 23, 32]}


def e5_drop():
    """Critical droplet: disk of radius R0 (random 0.5 inside) in the all-dead state,
    alpha=0.5, q > q*.  First order: survival steps from 0 to 1 at R_c ~ 1/(q - q*)."""
    T = 50000
    for q, radii in sorted(DROP_GRID.items()):
        for R in radii:
            L = max(1024, 1 << (int(8 * R) - 1).bit_length())
            for s in range(1, 33):
                yield (f"a0.5_j2_q{q_str(q)}_R{R}_L{L}_s{s}",
                       f"--alpha 0.5 --jflip 2:{q} --L {L} --T {T} --init disk:0.5:{R} --bg 0 "
                       f"--every 500 --ppd 10 --stop-absorbed",
                       L * L * T)


# q*(p) near the voter point: 0.98, 0.99 from the E4_voter fronts (v = 0 by interpolation,
# +-2%); 0.995 from a local strip scan (10-14, Lx=Ly=2048, T 1e5, v1.2, 2 seeds: v = -0.0016,
# -0.0001, +0.0007 at q = 0.00058, 0.00061, 0.00064; +-2%), which finished after 0.98/0.99
# were submitted, hence the separate E6_long_p995 (+ its L 4096 run alone, so a 16-core
# shard does not idle while one 8 h run finishes).
QSTAR_LONG = {0.98: 0.00252, 0.99: 0.00127}
QSTAR_P995 = {0.995: 0.00062}
LONG_GRID = (0.97, 1.0, 1.01, 1.02, 1.035, 1.05, 1.1, 1.2)


def e6_long(qstar=QSTAR_LONG, refs=True):
    """Does the jump survive near the voter point, or does a GV segment open up?  Bulk runs
    to T = 1e6 (10x Phase 2), densely just above q*(p).  First order: the active density
    settles at a finite plateau for every q > q*.  GV: it decays like the voter model
    (~1/ln t) at q* and the plateau -> 0 continuously as q -> q*+.  p = 1 references: the
    pure voter model (q = 0) and the voter model with the j=2 flip (q = 0.0003)."""
    L, T = 2048, 1000000
    for p, qs in sorted(qstar.items()):
        for q in [round(qs * f, 7) for f in LONG_GRID]:
            for s in (1, 2):
                yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_L{L}_s{s}",
                       f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --T {T} --every 1000 "
                       f"--ppd 20 --stop-absorbed",
                       L * L * T)
    for q in (0, 0.0003) if refs else ():
        for s in (1, 2):
            jf = f"--jflip 2:{q} " if q else ""
            yield (f"p1_a0.5_j2_q{q_str(q, 7)}_L{L}_s{s}",
                   f"--alpha 0.5 --voter 1 {jf}--L {L} --T {T} --every 1000 --ppd 20 --stop-absorbed",
                   L * L * T)


def e6_long_p995():
    """E6_long at p = 0.995."""
    yield from e6_long(QSTAR_P995, refs=False)


def e6_long_l4096_p995():
    """E6_long_L4096 at p = 0.995 (1.02 q*)."""
    yield from e6_long_l4096(QSTAR_P995, extra=False)


def e6_long_l4096(qstar=QSTAR_LONG, extra=True):
    """Size check for E6_long at the points that matter most (just above q*), and the pure
    voter reference, where coarsening domains would reach L = 2048 first."""
    L, T = 4096, 1000000
    pts = [(p, qstar[p] * 1.02) for p in sorted(qstar)] + ([(0.99, qstar[0.99] * 1.05)] if extra else [])
    for p, q in pts:
        q = round(q, 7)
        yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_L{L}_s1",
               f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --T {T} --every 1000 --ppd 20 "
               f"--stop-absorbed",
               L * L * T)
    if extra:
        yield (f"p1_a0.5_j2_q0_L{L}_s1",
               f"--alpha 0.5 --voter 1 --L {L} --T {T} --every 1000 --ppd 20 --stop-absorbed", L * L * T)


# E7 (README "Phase 4 design"): magnetisation FSS where E6 left the size check open.  Tier 1 is
# p = 0.99 at 1.02 and 1.05 q*.  Seed numbers continue E6's (L 2048 s1-s2, L 4096 s1) so the
# two experiments pool without name clashes; the RNG seeds differ anyway (crc32 of EXP/name).
E7_POINTS = [(0.99, 1.02), (0.99, 1.05)]
E7_SEEDS = {512: range(1, 9), 1024: range(1, 9), 2048: range(3, 7), 4096: range(2, 6)}


def e7_binder(sizes=(512, 1024, 2048)):
    """Is there an Ising transition inside the active phase near the voter point?  Bulk random
    start to T = 1e6, m sampled every 500 steps; chi_L = L^2 <m^2> and U4 against L over
    t >= 2e5.  First order (finite xi): chi flat in L, U4 -> 0.  Ising: chi ~ L^1.75 (x3.4 per
    doubling) and U4 curves crossing near 0.61."""
    T = 1000000
    for p, x in E7_POINTS:
        q = round(QSTAR_LONG[p] * x, 7)
        for L in sizes:
            for s in E7_SEEDS[L]:
                yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_L{L}_s{s}",
                       f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --T {T} --every 500 --ppd 20 "
                       f"--stop-absorbed",
                       L * L * T)


def e7_binder_l4096():
    """E7_binder at L 4096, separate so it can run 2 threads per run (THREADS=2 in run_shard.sh)."""
    yield from e7_binder((4096,))


# E8b (README "Phase 4 design"): q* from fronts (E4_voter at p 0.9-0.99; local strip scan at 0.995)
QSTAR_FRONT = {0.9: 0.0119, 0.98: 0.00252, 0.99: 0.00127, 0.995: 0.00062}


def e8b_wall():
    """Domain wall between the two absorbing states (all alive below Ly/2, all dead above;
    two flat walls on the torus), no active phase at t = 0.  Does an active layer form, and
    how does it grow?  First order: width W ~ 2 v(q) t, bounded at q*, interior density at
    the plateau.  GV: W ~ t^1/2 and interior density ~ 1/ln t at q_c.  Per-row profiles."""
    Lx, Ly, T = 1024, 4096, 1000000
    for p, qs in sorted(QSTAR_FRONT.items()):
        for x in (0.98, 1.0, 1.02, 1.05):
            q = round(qs * x, 7)
            for s in (1, 2):
                yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_Lx{Lx}_Ly{Ly}_s{s}",
                       f"--alpha 0.5 --voter {p} --jflip 2:{q} --Lx {Lx} --Ly {Ly} --T {T} "
                       f"--init strip:1:0.5 --bg 0 --every 10000 --ppd 20 --stop-absorbed "
                       f"--profile {{P}}.prof",
                       Lx * Ly * T)


def e8a_spread():
    """Seed spreading: a 2x2 live block on an all-dead L 2048 torus, many replicas (dnsim v1.3
    --reps: bounding-box stepping, exact).  Survival P_s(t), interface count N(t) and spread
    R2(t) to T = 1e5.  First order: no power-law window; a seed must reach R_c ~ 1/(q - q*)
    to survive.  GV: P_s ~ ln t / t, R2 ~ t at q_c.  DP: delta 0.451, eta 0.230, z 1.13.
    Replicas stop at `target` survivors to T (300; 100 at x >= 1.1, where survivors grow
    ballistically and cost ~1 min each) or at 3e6 replicas.  Calibration (local, 10-15): a
    dying replica costs 0.1-0.4 ms, a survivor 20-50 s; P_s(1e5) ~ 1e-3 at 1.2 q*."""
    L, T, cap = 2048, 100000, 3000000
    pts = [(p, round(qs * x, 7), 100 if x >= 1.1 else 300) for p, qs in sorted(QSTAR_FRONT.items())
           for x in (0.95, 1.0, 1.02, 1.05, 1.1, 1.2)]
    pts += [(1, 0, 300), (1, 0.0003, 100)]
    for p, q, target in pts:
        jf = f"--jflip 2:{q} " if q else ""
        yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_L{L}_seed",
               f"--alpha 0.5 --voter {p} {jf}--L {L} --T {T} --init disk:1:0.75 --reps {cap} "
               f"--target {target} --ppd 20",
               3 * 3600 * 8e8)   # cost unknown up front: pack as ~3 core-h each


# E7 rest (README "E7_binder" design, run 10-17): the near points tier 1 left out and the far
# points plus the p 0.9 control.  Seeds continue E6's where E6 has runs at that (p, q, L): L 2048
# s1-s2 at every LONG_GRID point; L 4096 s1 at 1.02 q* (p 0.98, 0.995) and 1.05 q* (p 0.99).
E7_REST_NEAR = [(0.98, 1.02), (0.99, 1.01), (0.995, 1.05)]
E7_REST_FAR = [(0.98, 1.2), (0.99, 1.2), (0.995, 1.2), (0.9, 1.03)]
QSTAR_E7 = {**QSTAR_LONG, **QSTAR_P995, 0.9: QSTAR_FRONT[0.9]}
E6_L4096 = {(0.98, 1.02), (0.99, 1.02), (0.99, 1.05), (0.995, 1.02)}


def e7_binder_rest(sizes=(512, 1024, 2048)):
    """E7_binder at the remaining points: near (8/8/4/4-5 seeds at L 512/1024/2048/4096, T 1e6)
    and far + control (4/4/2/2 seeds, T 5e5).  (0.99, 1.01) is the delta = 0.01 point: is xi
    still finite at the bulk threshold?"""
    for p, x in E7_REST_NEAR + E7_REST_FAR:
        near = (p, x) in E7_REST_NEAR
        T = 1000000 if near else 500000
        q = round(QSTAR_E7[p] * x, 7)
        for L in sizes:
            if L == 2048:
                seeds = range(3, 7) if near else (range(1, 3) if p == 0.9 else range(3, 5))
            elif L == 4096:
                seeds = (range(2, 6) if (p, x) in E6_L4096 else range(1, 6)) if near else range(1, 3)
            else:
                seeds = range(1, 9) if near else range(1, 5)
            for s in seeds:
                yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_L{L}_s{s}",
                       f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --T {T} --every 500 --ppd 20 "
                       f"--stop-absorbed",
                       L * L * T)


def e7_binder_rest_l4096():
    """E7_binder_rest at L 4096 (THREADS=2)."""
    yield from e7_binder_rest((4096,))


E8A_LONG_POINTS = [(0.98, 1.02), (0.98, 1.05), (0.99, 1.02), (0.99, 1.05), (0.995, 1.02), (0.995, 1.05)]
E8A_LONG_SHORT = {(0.98, 1.05): (300000, 10)}   # (T, target): ballistic survivors


def e8a_long():
    """E8a to T = 1e6: does P_s plateau (P_inf > 0) just above the bulk threshold, where E8a's
    delta_eff was still rising at 3e4?  Bulk thresholds (E6): 1.02, 1.01, 1.035 q* at p 0.98,
    0.99, 0.995, so (0.995, 1.02) is a below-threshold control that should keep falling faster
    than the voter curve.  L 8192 so that diffusive (rms ~1600) or ballistic (v ~ 2e-3) clusters
    do not reach the edge (escape) before 1e6; box stepping makes L free.  8 sub-runs per point
    (different seeds, pooled by spread.py), each up to 375k replicas or 5 survivors to T, so 40
    survivors or 3e6 replicas per point.  Survivors that spread ballistically cost ~ (2v)^2 T^3 / 3
    updates: at (0.98, 1.05), v = 0.0027, one survivor to 1e6 ran > 2.6 h, so that point (already
    turning over at 1e5 in E8a) stops at T = 3e5 with 10 survivors per sub-run (E8A_LONG_SHORT)."""
    L, K, cap = 8192, 8, 375000
    for p, x in E8A_LONG_POINTS:
        T, target = E8A_LONG_SHORT.get((p, x), (1000000, 5))
        q = round(QSTAR_FRONT[p] * x, 7)
        for k in range(1, K + 1):
            yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_L{L}_T{T:.0e}".replace("e+0", "e") + f"_k{k}",
                   f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --T {T} --init disk:1:0.75 "
                   f"--reps {cap} --target {target} --ppd 20",
                   5 * 3600 * 8e8)   # ~5 core-h each (unknown up front)


ALPHA_GRID = (0.9, 0.95, 1.0, 1.05, 1.1)


def e9_alpha():
    """q* proportional to alpha near the voter point (theory/vpt.py: q* = c alpha (1-p)): E4 strip
    fronts at alpha 0.25, 0.5, 1 and p 0.98, 0.99, q = x * 2 alpha q*(p, alpha 0.5), 4 seeds.
    alpha 0.5 is rerun so that all three use one protocol."""
    Lx, Ly, T = 2048, 4096, 100000
    for p in (0.98, 0.99):
        for a in (0.25, 0.5, 1):
            for x in ALPHA_GRID:
                q = round(x * 2 * a * QSTAR_FRONT[p], 7)
                for s in range(1, 5):
                    yield (f"p{p}_a{a}_j2_q{q_str(q, 7)}_Lx{Lx}_Ly{Ly}_s{s}",
                           f"--alpha {a} --voter {p} --jflip 2:{q} --Lx {Lx} --Ly {Ly} --T {T} "
                           f"--init strip:0.5:0.5 --bg 0 --every 1000 --ppd 10 --stop-absorbed "
                           f"--profile {{P}}.prof",
                           Lx * Ly * T)


# E10 (10-17): asynchrony check.  Each cell updates with probability r per step: the flip table
# is scaled by r and passed verbatim (--ftab), so the kernel is unchanged and r = 1 is
# bit-identical to --rule/--jflip.  r -> 0 is continuous time; t r is the physical time.
LAZY_B, LAZY_S = {3, 6, 7, 8}, {3, 4, 6, 7, 8}
# Local pilot (Lx 256 x Ly 2048, 1 seed, p 0.98): v changes sign near 1.11 (r 0.5) and 1.15-1.25
# (r 0.25, 0.125) x q*(r = 1), so the grids below are in units of QSTAR_FRONT and span 1.0-1.36.
LAZY_FRONT_POINTS = [(0.98, 0.5), (0.98, 0.25), (0.98, 0.125), (0.99, 0.25)]
LAZY_FRONT_GRID = (1.0, 1.06, 1.12, 1.18, 1.24, 1.30, 1.36)


def lazy_ftab(p, alpha, q, r):
    """f[s][k] = r [(1-p) alpha det(s, k) + p j/8 + q [j == 2]], k alive neighbours, j disagreeing."""
    rows = []
    for st in (0, 1):
        row = []
        for k in range(9):
            j = 8 - k if st else k
            det = (k not in LAZY_S) if st else (k in LAZY_B)
            row.append(r * ((1 - p) * alpha * det + p * j / 8 + (q if j == 2 else 0)))
        rows.append(",".join(f"{f:.12g}" for f in row))
    return ";".join(rows)


def e10_lazy_front(grid=LAZY_FRONT_GRID, seeds=range(1, 5)):
    """q*(r): E4 strip fronts (Lx 2048 x Ly 4096) to physical time 1e5 (T = 1e5 / r) at
    q = x q*(p, r = 1), x in LAZY_FRONT_GRID, 4 seeds.  v from fronts.py is per step (physical v / r)."""
    Lx, Ly = 2048, 4096
    for p, r in LAZY_FRONT_POINTS:
        T = int(round(1e5 / r))
        for x in grid:
            q = round(x * QSTAR_FRONT[p], 7)
            for sd in seeds:
                yield (f"p{p}_a0.5_r{r}_j2_q{q_str(q, 7)}_Lx{Lx}_Ly{Ly}_s{sd}",
                       f"--ftab {lazy_ftab(p, 0.5, q, r)} --Lx {Lx} --Ly {Ly} --T {T} "
                       f"--init strip:0.5:0.5 --bg 0 --every {int(1000 / r)} --ppd 10 --stop-absorbed "
                       f"--profile {{P}}.prof",
                       Lx * Ly * T)


LAZY_BULK_GRID = (1.06, 1.12, 1.18, 1.24, 1.30)   # x q*(p, r = 1), around the pilot's q*(0.25)


def e10_lazy_bulk(grid=LAZY_BULK_GRID):
    """Bulk under asynchrony, as E6 (random init, L 2048, physical T 1e6) at r = 0.25, p = 0.98:
    is there a stationary plateau with a finite interface density just above q*(r) (first order),
    or does rho -> 0 continuously as q -> q*+ (GV)?  Grid in units of q*(p, r = 1)."""
    L, p, r = 2048, 0.98, 0.25
    T = int(round(1e6 / r))
    for x in grid:
        q = round(x * QSTAR_FRONT[p], 7)
        for sd in (1, 2):
            yield (f"p{p}_a0.5_r{r}_j2_q{q_str(q, 7)}_L{L}_s{sd}",
                   f"--ftab {lazy_ftab(p, 0.5, q, r)} --L {L} --T {T} --every {int(1000 / r)} "
                   f"--ppd 20 --stop-absorbed",
                   L * L * T)


# Lazy theory (CWALK_R, vpt.py) predicts c(r) = 0.208 (r 0.25), 0.201 (r 0.125) vs 0.26 at r = 1,
# i.e. q*(r) ~ 0.8 q*(1); pilot 2 (Lx 512, 2 seeds) instead puts q*(0.25) at ~1.0-1.05 q*(1).
# The _lo experiments extend the grids down to the theory's value so E10 decides between them.
LAZY_FRONT_GRID_LO = (0.76, 0.82, 0.88, 0.94)
LAZY_BULK_GRID_LO = (0.94, 1.0)


def e10_lazy_front_lo():
    """E10_lazy_front below x = 1 (3 seeds)."""
    yield from e10_lazy_front(LAZY_FRONT_GRID_LO, range(1, 4))


def e10_lazy_bulk_lo():
    """E10_lazy_bulk at x = 0.94 and 1.0."""
    yield from e10_lazy_bulk(LAZY_BULK_GRID_LO)


# (p, x, T): P_inf(delta) between 1.02 (still barrier-limited at 1e6) and 1.05 (plateau) at p 0.98.
E8A_PINF_POINTS = [(0.98, 1.03, 500000), (0.98, 1.04, 500000), (0.98, 1.07, 300000)]


def e8a_pinf():
    """E8a_long's protocol (L 8192, 8 sub-runs of <= 375k replicas) at intermediate x, 10 survivors
    per sub-run: is the onset of P_inf as steep as a critical droplet requires?"""
    L, K, cap, target = 8192, 8, 375000, 10
    for p, x, T in E8A_PINF_POINTS:
        q = round(QSTAR_FRONT[p] * x, 7)
        for k in range(1, K + 1):
            yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_L{L}_T{T:.0e}".replace("e+0", "e") + f"_k{k}",
                   f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --T {T} --init disk:1:0.75 "
                   f"--reps {cap} --target {target} --ppd 20",
                   2 * 3600 * 1.5e9)   # ~2 core-h each (E8a_long: ~23 ms per replica at T 1e6)


# E11 (10-18): the full first-order signature under asynchrony.  E10 put q*(r) at 0.864 q*(1)
# for p 0.98, r 0.5 (the shift is saturated by r 0.5: 0.858 at 0.25, 0.870 at 0.125), with a
# fit-window systematic of up to +0.04, so the seed and wall grids are wider than E8's.
QSTAR_LAZY = {(0.98, 0.5): 0.0021773}
LAZY_SEED_GRID = (0.95, 1.0, 1.05, 1.1, 1.2)
LAZY_WALL_GRID = (0.95, 1.0, 1.05, 1.1)


def e11_lazy_seed():
    """E8a at r 0.5, p 0.98: seeds on an all-dead L 2048 torus to physical t 1e5 (T 2e5),
    x q*(0.5), plus the lazy voter (p 1, q 0) as the reference curve.  4 sub-runs per point
    (pooled by spread.py), each up to 750k replicas or 75 survivors (25 at x >= 1.1)."""
    L, p, r, K, cap = 2048, 0.98, 0.5, 4, 750000
    T = int(round(1e5 / r))
    pts = [(p, round(QSTAR_LAZY[(p, r)] * x, 7), 25 if x >= 1.1 else 75) for x in LAZY_SEED_GRID]
    for pp, q, target in pts + [(1, 0, 75)]:
        for k in range(1, K + 1):
            yield (f"p{pp}_a0.5_r{r}_j2_q{q_str(q, 7)}_L{L}_k{k}",
                   f"--ftab {lazy_ftab(pp, 0.5, q, r)} --L {L} --T {T} --init disk:1:0.75 "
                   f"--reps {cap} --target {target} --ppd 20",
                   3 * 3600 * 9e8)   # ~3 core-h each (unknown up front)


def e11_lazy_wall():
    """E8b at r 0.5, p 0.98: wall between the two absorbing states (Lx 1024 x Ly 4096) to physical
    t 1e6 (T 2e6), x q*(0.5), 2 seeds, plus the lazy voter wall (p 1, q 0)."""
    Lx, Ly, p, r = 1024, 4096, 0.98, 0.5
    T = int(round(1e6 / r))
    pts = [(p, round(QSTAR_LAZY[(p, r)] * x, 7)) for x in LAZY_WALL_GRID] + [(1, 0)]
    for pp, q in pts:
        for s in (1, 2):
            yield (f"p{pp}_a0.5_r{r}_j2_q{q_str(q, 7)}_Lx{Lx}_Ly{Ly}_s{s}",
                   f"--ftab {lazy_ftab(pp, 0.5, q, r)} --Lx {Lx} --Ly {Ly} --T {T} "
                   f"--init strip:1:0.5 --bg 0 --every {int(1e4 / r)} --ppd 20 --stop-absorbed "
                   f"--profile {{P}}.prof",
                   Lx * Ly * T)


LAZY_FRONT_P99_GRID = (0.84, 0.88, 0.92, 0.96, 1.0)


def e11_lazy_front_p99():
    """Is the p 0.99, r 0.25 crossing (x* 0.93-0.95 at physical t 1e5, against 0.86-0.87 at p 0.98)
    converged?  E10 strips to physical t 3e5 (T 1.2e6), 4 seeds."""
    Lx, Ly, p, r = 2048, 4096, 0.99, 0.25
    T = int(round(3e5 / r))
    for x in LAZY_FRONT_P99_GRID:
        q = round(x * QSTAR_FRONT[p], 7)
        for sd in range(1, 5):
            yield (f"p{p}_a0.5_r{r}_j2_q{q_str(q, 7)}_Lx{Lx}_Ly{Ly}_s{sd}",
                   f"--ftab {lazy_ftab(p, 0.5, q, r)} --Lx {Lx} --Ly {Ly} --T {T} "
                   f"--init strip:0.5:0.5 --bg 0 --every {int(1000 / r)} --ppd 10 --stop-absorbed "
                   f"--profile {{P}}.prof",
                   Lx * Ly * T)


# Phase 8 (10-27): the clean model.  p > 0 has exactly two absorbing states (every j >= 1 flips
# with prob >= p j/8), p = 0 has exponentially many frozen ones, so the paper's flagship moves to
# p = 0.5 (strongly first order; q* = 0.0459 from E4) and the phase diagram gets the small-p
# points that connect it to the p = 0 limit.  q* below: E4 (0.5); local pilots 10-27 (scratchpad
# p8/: Lx 512 strips, 1 seed, for 0.01-0.25; Lx 1024, 2 seeds for 0.1, where q* = 0.0368).
QSTAR_P8 = {0.01: 0.021, 0.02: 0.0212, 0.05: 0.029, 0.1: 0.0368, 0.25: 0.049, 0.5: 0.0459}
PHASE_GRID = (0.9, 0.94, 0.97, 1.0, 1.03, 1.06, 1.1)
PHASE_GRID_05 = (0.96, 0.98, 0.99, 1.0, 1.01, 1.02, 1.04)


def e12_phase():
    """q*(p) between the D&N limit and p = 0.5: E4 strip fronts (Lx 2048 x Ly 4096, T 1e5), 3 seeds
    (p = 0.5: 4 seeds on a fine grid, the flagship crossing)."""
    Lx, Ly, T = 2048, 4096, 100000
    for p, qs in sorted(QSTAR_P8.items()):
        grid, seeds = (PHASE_GRID_05, 4) if p == 0.5 else (PHASE_GRID, 3)
        for x in grid:
            q = round(qs * x, 7)
            for s in range(1, seeds + 1):
                yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_Lx{Lx}_Ly{Ly}_s{s}",
                       f"--alpha 0.5 --voter {p} --jflip 2:{q} --Lx {Lx} --Ly {Ly} --T {T} "
                       f"--init strip:0.5:0.5 --bg 0 --every 1000 --ppd 10 --stop-absorbed "
                       f"--profile {{P}}.prof",
                       Lx * Ly * T)


def e12_phase_bulk():
    """Bulk partner of E12_phase at the new p: the active density just above q* (the jump) and the
    collapse below.  L 2048, T 1e5, 2 seeds."""
    L, T = 2048, 100000
    for p, qs in sorted(QSTAR_P8.items()):
        if p == 0.5:
            continue   # E13_bulk05
        for x in (0.97, 1.0, 1.03, 1.06, 1.1):
            q = round(qs * x, 7)
            for s in (1, 2):
                yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_L{L}_s{s}",
                       f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --T {T} --every 1000 --ppd 20 "
                       f"--stop-absorbed",
                       L * L * T)


BULK05_GRID = (0.9, 0.95, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.05, 1.1)
BULK05_SIZES = {64: (32, 1000000), 128: (32, 1000000), 256: (16, 1000000), 512: (8, 1000000),
                1024: (8, 1000000), 2048: (4, 300000), 4096: (2, 200000)}   # L: (seeds, T)


def e13_bulk05():
    """Flagship bulk FSS at p = 0.5: (i) the active density above q* against L (a jump that does not
    shrink with L); (ii) the lifetime of the metastable active state below q* against L: single-droplet
    (tau ~ 1 / (I L^2)) at small L, L-independent multi-droplet (Avrami) decay at large L; (iii) the
    lifetime above q* growing with L; (iv) chi_L = L^2 <m^2> flat (paramagnetic active phase)."""
    for L, (seeds, T) in sorted(BULK05_SIZES.items()):
        for x in BULK05_GRID:
            q = round(QSTAR_P8[0.5] * x, 7)
            for s in range(1, seeds + 1):
                yield (f"p0.5_a0.5_j2_q{q_str(q, 7)}_L{L}_s{s}",
                       f"--alpha 0.5 --voter 0.5 --jflip 2:{q} --L {L} --T {T} --every 500 --ppd 20 "
                       f"--stop-absorbed",
                       L * L * T)


# x -> (radii, L, T, replicas).  Local pilot (10-27, L 2048, 200 replicas, T 2e4): P_surv = 1/2 near
# R ~ 20 (x 1.05), 11 (1.1), 8 (1.2), 4 (1.4); 2x2 seeds: 0 of 200 survive at every x <= 1.4.  L is
# small so survivors escape (box within 1 word of the edge: radius ~ L/2 - 64 ~ 8-10 R_c) soon after
# they are clearly supercritical; escaped replicas count as survivors.  disk:0.5:R fills the cells
# whose centres lie within R of the lattice centre (R 0.75 = the 2x2 block, R 2 = 12 cells).
DROP05 = {1.02: ([8, 11, 16, 23, 32, 45, 64, 90], 1024, 300000, 500),
          1.05: ([4, 6, 8, 11, 16, 23, 32, 45], 512, 100000, 1000),
          1.1: ([2, 3, 4, 6, 8, 11, 16, 23, 32], 512, 50000, 1000),
          1.2: ([2, 3, 4, 6, 8, 11, 16, 23], 512, 50000, 1000),
          1.4: ([0.75, 2, 3, 4, 6, 8, 11, 16], 512, 30000, 1000),
          2.0: ([0.75, 2, 3, 4, 6, 8, 11], 512, 20000, 1000)}


def e14_drop05():
    """Critical droplets at p = 0.5 (E5 for the clean model, with dnsim v1.3 --reps): a disordered
    disk of radius R (alive with prob 1/2 inside) on the all-dead torus.  First order: P_surv(R)
    steps from 0 to 1 at R_c(q), with R_c v(q) ~ const (curvature-limited growth, dR/dt = v - D/R)."""
    for x, (radii, L, T, reps) in sorted(DROP05.items()):
        q = round(QSTAR_P8[0.5] * x, 7)
        for R in radii:
            yield (f"p0.5_a0.5_j2_q{q_str(q, 7)}_R{R}_L{L}_drop",
                   f"--alpha 0.5 --voter 0.5 --jflip 2:{q} --L {L} --T {T} --init disk:0.5:{R} "
                   f"--reps {reps} --ppd 10",
                   reps * 0.5 * (L / 2) ** 2 * (L / 2) / max(0.21 * (x - 1), 1e-3) / 3)


# Orientation (generic coexistence) test: the same band at three normal directions (a, b) =
# (1, 0), (2, 1), (1, 1) on an L x L torus, from checkpoints made by tilt.py (one per (a, b, seed),
# shared by every q, so the q-dependence is measured on a common initial condition).
ORIENT_POINTS = {0.1: 0.0368, 0.5: 0.0459, 0.9: QSTAR_FRONT[0.9]}
ORIENT_GRID = (0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03)
ORIENT_AB = ((1, 0), (2, 1), (1, 1))
ORIENT_L, ORIENT_T, ORIENT_SEEDS = 2048, 100000, 4


def orient_ckpt(a, b, s):
    return os.path.join(ROOT, "runs", "E15_orient", "init", f"band_L{ORIENT_L}_ab{a}-{b}_s{s}.ckpt")


def e15_orient():
    """Front velocity v(q) of a tilted band (tilt.py: dI/dt of the interface count) at p 0.1, 0.5, 0.9.
    Generic coexistence would show up as different zero crossings q*(a, b)."""
    L, T = ORIENT_L, ORIENT_T
    for p, qs in sorted(ORIENT_POINTS.items()):
        for a, b in ORIENT_AB:
            for x in ORIENT_GRID:
                q = round(qs * x, 7)
                for s in range(1, ORIENT_SEEDS + 1):
                    yield (f"p{p}_a0.5_q{q_str(q, 7)}_ab{a}-{b}_L{L}_s{s}",
                           f"--alpha 0.5 --voter {p} --jflip 2:{q} --L {L} --load {orient_ckpt(a, b, s)} "
                           f"--T {T} --every 500 --ppd 10 --stop-absorbed",
                           L * L * T)


WALL05_GRID = (0.8, 0.9, 0.95, 0.97, 0.98, 0.99, 0.995)
WALL05_LX = {256: 4, 1024: 4, 4096: 2}   # Lx: seeds


def e16_wall():
    """Wetting of the wall between the two absorbing states by the active phase below q*: E8b walls
    at p = 0.5 (Lx 256 / 1024 / 4096 x Ly 1024, T 1e6), the stationary layer thickness against
    delta = 1 - q/q* (and against Lx: wall wandering), plus the bare wall (q = 0); and p = 0.9
    (Lx 1024 x Ly 2048, 2 seeds) to compare with E8b."""
    Ly, T = 1024, 1000000
    pts = [(0.5, x, Lx, s) for x in WALL05_GRID for Lx, ns in WALL05_LX.items() for s in range(1, ns + 1)]
    pts += [(0.5, 0, 1024, s) for s in (1, 2)]
    for p, x, Lx, s in pts:
        q = round(QSTAR_P8[0.5] * x, 7)
        jf = f"--jflip 2:{q} " if q else ""
        yield (f"p{p}_a0.5_j2_q{q_str(q, 7)}_Lx{Lx}_Ly{Ly}_s{s}",
               f"--alpha 0.5 --voter {p} {jf}--Lx {Lx} --Ly {Ly} --T {T} --init strip:1:0.5 --bg 0 "
               f"--every 10000 --ppd 20 --stop-absorbed --profile {{P}}.prof",
               Lx * Ly * T / 2)   # absorbing domains run ~2-4x faster
    Lx, Ly = 1024, 2048
    for x in (0.9, 0.95, 0.98, 0.99):
        q = round(QSTAR_FRONT[0.9] * x, 7)
        for s in (1, 2):
            yield (f"p0.9_a0.5_j2_q{q_str(q, 7)}_Lx{Lx}_Ly{Ly}_s{s}",
                   f"--alpha 0.5 --voter 0.9 --jflip 2:{q} --Lx {Lx} --Ly {Ly} --T {T} --init strip:1:0.5 "
                   f"--bg 0 --every 10000 --ppd 20 --stop-absorbed --profile {{P}}.prof",
                   Lx * Ly * T / 2)


def e16_wall_ref():
    """The synchronous pure voter wall (p = 1, q = 0) in the E8b geometry (Lx 1024 x Ly 4096, T 1e6),
    the reference that E8b compared against only through a local one-seed run."""
    Lx, Ly, T = 1024, 4096, 1000000
    for s in (1, 2):
        yield (f"p1_a0.5_j2_q0_Lx{Lx}_Ly{Ly}_s{s}",
               f"--alpha 0.5 --voter 1 --Lx {Lx} --Ly {Ly} --T {T} --init strip:1:0.5 --bg 0 "
               f"--every 10000 --ppd 20 --stop-absorbed --profile {{P}}.prof",
               Lx * Ly * T)


# Single-noise variant (in-house check 10-27): w(j) = alpha [j in {3,6,7,8}] + eta [j in {1,2,4,5}],
# exactly two absorbing states, no voter admixture; eta* ~ 0.0125 at alpha 0.5 (Lx 1024 strips, 1 seed).
ETA_STAR = 0.0125
ETA_GRID = (0.9, 0.95, 0.98, 1.0, 1.02, 1.05, 1.1)


def eta_ftab(alpha, eta):
    rows = []
    for st in (0, 1):
        row = []
        for k in range(9):
            j = 8 - k if st else k
            row.append(alpha * (j in (3, 6, 7, 8)) + eta * (j in (1, 2, 4, 5)))
        rows.append(",".join(f"{f:.12g}" for f in row))
    return ";".join(rows)


def e17_eta_front():
    """The single-noise variant: strip fronts (Lx 2048 x Ly 4096, T 1e5, 3 seeds).  Run names carry
    eta as q."""
    Lx, Ly, T = 2048, 4096, 100000
    for x in ETA_GRID:
        eta = round(ETA_STAR * x, 7)
        for s in (1, 2, 3):
            yield (f"eta_a0.5_q{q_str(eta, 7)}_Lx{Lx}_Ly{Ly}_s{s}",
                   f"--ftab {eta_ftab(0.5, eta)} --Lx {Lx} --Ly {Ly} --T {T} --init strip:0.5:0.5 --bg 0 "
                   f"--every 1000 --ppd 10 --stop-absorbed --profile {{P}}.prof",
                   Lx * Ly * T)


def e17_eta_bulk():
    """The single-noise variant in bulk: L 1024 and 2048, T 2e5, 2 seeds, on the E17_eta_front grid."""
    T = 200000
    for x in ETA_GRID:
        eta = round(ETA_STAR * x, 7)
        for L in (1024, 2048):
            for s in (1, 2):
                yield (f"eta_a0.5_q{q_str(eta, 7)}_L{L}_s{s}",
                       f"--ftab {eta_ftab(0.5, eta)} --L {L} --T {T} --every 1000 --ppd 20 --stop-absorbed",
                       L * L * T)


# E18 (10-27): asynchrony at the flagship.  Local pilot (Lx 512 strips, 2 seeds, physical t 1.5e4):
# at p 0.5, r 0.25 the per-step v is -0.0034, +0.0018, +0.0039 at x = 0.85, 0.90, 0.95 x q*(1), so
# q*(0.25) ~ 0.88 q*(1) (p 0.98: 0.86).  Grids in units of q*(p 0.5, r 1) = 0.0459.
LAZY05_FRONT_GRID = (0.84, 0.86, 0.87, 0.88, 0.89, 0.9, 0.92, 0.96)
LAZY05_BULK_GRID = (0.86, 0.88, 0.9, 0.92, 0.96)


def e18_lazy05_front():
    """E4 strips at p 0.5, r 0.25 (Lx 2048 x Ly 4096) to physical t 5e4 (T 2e5), 3 seeds."""
    Lx, Ly, p, r = 2048, 4096, 0.5, 0.25
    T = int(round(5e4 / r))
    for x in LAZY05_FRONT_GRID:
        q = round(x * QSTAR_P8[0.5], 7)
        for sd in (1, 2, 3):
            yield (f"p{p}_a0.5_r{r}_j2_q{q_str(q, 7)}_Lx{Lx}_Ly{Ly}_s{sd}",
                   f"--ftab {lazy_ftab(p, 0.5, q, r)} --Lx {Lx} --Ly {Ly} --T {T} "
                   f"--init strip:0.5:0.5 --bg 0 --every {int(1000 / r)} --ppd 10 --stop-absorbed "
                   f"--profile {{P}}.prof",
                   Lx * Ly * T)


def e18_lazy05_bulk():
    """Bulk at p 0.5, r 0.25: L 2048 to physical t 1e5 (T 4e5), 2 seeds."""
    L, p, r = 2048, 0.5, 0.25
    T = int(round(1e5 / r))
    for x in LAZY05_BULK_GRID:
        q = round(x * QSTAR_P8[0.5], 7)
        for sd in (1, 2):
            yield (f"p{p}_a0.5_r{r}_j2_q{q_str(q, 7)}_L{L}_s{sd}",
                   f"--ftab {lazy_ftab(p, 0.5, q, r)} --L {L} --T {T} --every {int(1000 / r)} --ppd 20 "
                   f"--stop-absorbed",
                   L * L * T)


def steps_left(txt):
    """None if the run never started, 0 if finished, else T minus its last output time."""
    if not os.path.exists(txt):
        return None
    T, last = None, 0
    for line in open(txt):
        if line.startswith("# done"):
            return 0
        if line.startswith("# Lx"):
            T = int(line.split()[6])
        elif line.strip() and not line.startswith("#"):
            last = int(line.split()[0])
    return None if T is None else T - last


EXPERIMENTS = {"E0_map": e0_map, "E0_fine": e0_fine, "E1_front": e1_front, "E2_bulk": e2_bulk,
               "E3_alpha": e3_alpha, "E4_voter": e4_voter, "E4_voter_bulk": e4_voter_bulk,
               "E5_drop": e5_drop, "E6_long": e6_long, "E6_long_L4096": e6_long_l4096,
               "E6_long_p995": e6_long_p995, "E6_long_L4096_p995": e6_long_l4096_p995,
               "E7_binder": e7_binder, "E7_binder_L4096": e7_binder_l4096,
               "E8b_wall": e8b_wall, "E8a_spread": e8a_spread,
               "E7_binder_rest": e7_binder_rest, "E7_binder_rest_L4096": e7_binder_rest_l4096,
               "E8a_long": e8a_long, "E9_alpha": e9_alpha,
               "E10_lazy_front": e10_lazy_front, "E10_lazy_bulk": e10_lazy_bulk, "E8a_pinf": e8a_pinf,
               "E10_lazy_front_lo": e10_lazy_front_lo, "E10_lazy_bulk_lo": e10_lazy_bulk_lo,
               "E11_lazy_seed": e11_lazy_seed, "E11_lazy_wall": e11_lazy_wall,
               "E11_lazy_front_p99": e11_lazy_front_p99,
               "E12_phase": e12_phase, "E12_phase_bulk": e12_phase_bulk, "E13_bulk05": e13_bulk05,
               "E14_drop05": e14_drop05, "E15_orient": e15_orient, "E16_wall": e16_wall,
               "E17_eta_front": e17_eta_front, "E17_eta_bulk": e17_eta_bulk,
               "E16_wall_ref": e16_wall_ref,
               "E18_lazy05_front": e18_lazy05_front, "E18_lazy05_bulk": e18_lazy05_bulk}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(os.environ.get("SCRATCH", os.path.expanduser("~")), "dn"))
    ap.add_argument("--rate", type=float, default=7e8, help="cell updates/s per core")
    ap.add_argument("--cpus", type=int, default=16, help="runs in parallel per shard")
    ap.add_argument("--hours", type=float, default=3.0, help="target wall hours per shard")
    ap.add_argument("--ckpt-min", type=float, default=10.0, help="minutes between checkpoints")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--remaining", action="store_true",
                    help="resubmit: skip finished runs, cost partial ones by the steps left "
                         "(keeps manifest.tsv)")
    a = ap.parse_args()
    global ROOT
    ROOT = a.root

    grand = 0.0
    for ei, (exp, gen) in enumerate(EXPERIMENTS.items()):
        if a.only and exp not in a.only:
            continue
        runs = []
        out_dir = os.path.join(a.root, "runs", exp, "out")
        for name, args, cost in gen():
            prefix = os.path.join(out_dir, name)
            seed = ((ei + 1) << 32) | zlib.crc32(f"{exp}/{name}".encode())
            cells = cost / int(args.split("--T ")[1].split()[0])
            ckpt_every = max(1000, int(round(a.rate * 60 * a.ckpt_min / cells, -3)))
            if "--reps" in args:
                ckpt_every = 1000   # replicas between checkpoints (plus every 10 min)
            full = (f"{MODEL} {args.replace('{P}', prefix)} --seed {seed} --out {prefix}.txt "
                    f"--ckpt {prefix}.ckpt --ckpt-every {ckpt_every} --resume")
            if a.remaining:
                left = steps_left(prefix + ".txt")
                if left == 0:
                    continue
                if left is not None:
                    T = int(args.split("--T ")[1].split()[0])
                    cost *= left / T
            runs.append((cost / a.rate, prefix, full))
        if not runs:
            print(f"{exp:9s} nothing left to run")
            continue
        core_h = sum(r[0] for r in runs) / 3600
        longest = max(r[0] for r in runs) / 3600
        nshard = max(1, math.ceil(core_h / (a.cpus * a.hours * 0.85)))
        shards = [[0.0, []] for _ in range(nshard)]
        for cost_s, prefix, full in sorted(runs, reverse=True):
            sh = min(shards, key=lambda s: s[0])
            sh[0] += cost_s
            sh[1].append((prefix, full))
        span = max(s[0] for s in shards) / a.cpus / 3600
        print(f"{exp:9s} runs {len(runs):4d}  est {core_h:8.1f} core-h  longest run {longest:5.2f} h  "
              f"shards {nshard:3d} x {a.cpus} cpus  est wall/shard ~{max(span, longest):.2f} h")
        grand += core_h
        if a.dry:
            continue
        tdir = os.path.join(a.root, "runs", exp, "tasks")
        os.makedirs(tdir, exist_ok=True)
        os.makedirs(out_dir, exist_ok=True)
        for f in os.listdir(tdir):
            if f.startswith("shard_"):
                os.remove(os.path.join(tdir, f))
        for i, (_, lines) in enumerate(shards):
            with open(os.path.join(tdir, f"shard_{i:03d}.txt"), "w") as f:
                for prefix, full in lines:
                    f.write(f"{prefix}\t{full}\n")
        if a.remaining:
            continue
        with open(os.path.join(a.root, "runs", exp, "manifest.tsv"), "w") as f:
            f.write("prefix\test_core_s\targs\n")
            for cost_s, prefix, full in runs:
                f.write(f"{prefix}\t{cost_s:.0f}\t{full}\n")
    print(f"total est {grand:.1f} core-h at {a.rate:.2g} updates/s/core")


if __name__ == "__main__":
    main()
