# Direct simulation of Molofsky's five-cell rule near its voter point (10-30)

Test of the lineage-closure prediction in `campaign/theory/MOLOFSKY.md` (Secs. 3-7): a thin
bistable window around the line p_2 - 2/5 = -3.46 (p_1 - 1/5), in which the absorbing state is
linearly stable, the active state u = 1/2 is stable, and the planar front between them
reverses. Everything here is in `campaign/molofsky_sim/`; `paper/main.tex` was not touched.

## 1. Code: dnsim v1.4 (`dnsim.c`, `dnsim_v14.patch` against the root `dnsim.c`)

    cd campaign/molofsky_sim && make && ./dnsim --selftest          # SELFTEST PASSED
    ./dnsim --molofsky p1,p2 [--Lx --Ly --T --init --bg --reps ...]   # five-cell rule
    ./dnsim --nbhd vn --ftab ...                                      # any von Neumann table

- `--molofsky p1,p2`: square lattice, synchronous; the next state of a site is 1 with
  probability p_n, n = number of 1s among the site and its four nearest neighbours, p_0 = 0,
  p_5 = 1, p_3 = 1 - p_2, p_4 = 1 - p_1 (Molofsky, Durrett, Dushoff, Griffeath & Levin, Theor.
  Popul. Biol. 55, 270 (1999)). Over k = live von Neumann neighbours the flip table is
  f[0][k] = p_k and f[1][k] = 1 - p_{k+1} = p_{4-k}, so f[1][k] = f[0][4-k] (exact Z2 symmetry)
  and f(j = 0) = 0 for both states: both uniform states are absorbing. Voter point (1/5, 2/5).
- Implementation: a `vn` flag in `Model`, a bit-sliced von Neumann count `count_eq_vn` (k in
  0..4, eq[5..8] = 0) selected in `step_word`; everything else (Philox streams, shared
  bit-serial sampler, strips, disks, profiles, checkpoints, `--reps` box stepping) is reused
  unchanged. The Moore path is untouched: a D&N voter-mixed strip run (Lx 256 x Ly 256, T 500,
  seed 5) gives a byte-identical time series and profile with the root `dnsim` (v1.3).
- Selftest additions: bit-sliced vs naive von Neumann count for 5 deterministic rules x 4 lattice
  shapes (0 mismatches in 1.8e6 cells); the Molofsky table at the voter point (f = j/5) and its
  Z2 symmetry at (0.17, 0.53); per-class flip frequencies of the 4-value von Neumann table
  against the naive count (max |z| = 2.0); box stepping == full stepping for the Molofsky rule
  over 120 seeds (0 mismatches). Speed: 1.1e9 cell updates/s per thread (same as D&N), 7.8e9 on
  8 threads; 16 threads is slower (8 physical cores), so all scans run 8 single-thread jobs.
- Observables as before: `rho_bond` = fraction of disagreeing nearest-neighbour bonds (the
  natural order parameter for this rule), `iface_frac` = cells with a disagreeing Moore
  neighbour (a superset of the cells that can flip, so `--stop-absorbed` and the `--reps` box
  invariant stay exact). Profiles hold per-row alive and Moore-interface counts.

Scripts: `mfront.py` (front tracking: `front`, seed-pooled `table`, `wall`, `bulk`),
`mreps.py` (seed spreading / droplets: <A>(t), P_surv, late-time rate), `scan.sh` (strip
scans, skips finished runs), `run_all*.sh` (the chains actually run).

## 2. Verification (`verify/`)

- **Voter point (0.2, 0.4).** Density is conserved on average: from `random:0.3` on L 1024,
  the alive fraction over 8 seeds is 0.3000 (t 0), 0.2978 +- 0.0024 (t 1e3), 0.3047 +- 0.0076
  (t 1e4) (seed values 0.283-0.345: the single-run drift is the voter martingale noise). A
  2 x 2 seed (`--reps 20000`, L 1024, T 1e4): <A>(t) = 4.00, 3.98, 3.93, 4.11, 3.95, 3.20 +- 0.98
  at t = 0, 1, 10, 100, 1e3, 1e4 (martingale, = 4); P_surv = 0.119 (1e2), 0.0170 (1e3), 0.0015 (1e4),
  i.e. t P_s / ln t = 2.6, 2.5, 1.6, the Bramson-Griffeath form seen for the Moore voter model.
  Coarsening from random 1/2 on L 2048: 1/rho_bond = 2.71, 3.42, 4.19, 4.96, 5.68 at t = 10, 1e2,
  1e3, 1e4, 6.3e4, i.e. linear in ln t with slope 0.31-0.39 per unit ln t (1/rho ~ ln t).
- **Their bottom-edge points, wall geometry** (all-1 strip against all-0, Lx 1024 x Ly 4096,
  their Figs. 4-6): at (0.35, 0) the mixed zone grows linearly, 121, 348, 1158, 2321, 3476 active
  rows at t = 100, 300, 1e3, 2e3, 3e3 (0.58 rows/step per front), and the bulk from random 1/2
  sits at rho_bond = 0.344 from t = 100 to 1e4 (ergodic). At (0.31, 0) the mixed zone grows at
  490, 1036, 1583 rows at t = 1e3, 2e3, 3e3 (0.27 rows/step per front, "about half the speed of
  the ergodic instance" as they write), with the bulk density drifting (phase separation). At
  (0.27, 0) the wall stays tight (15-33 active rows up to t = 2e4) and the bulk coarsens
  (rho_bond 0.065, 0.022, 0.010 at t = 1e3, 1e4, 1e5): clustering.

## 3. Planar fronts across the predicted line (`pilot/`, `scan/`, `scan_wide/`)

Protocol (as E4/E12 for D&N): a strip of the active phase (random density 1/2, rows [0, Ly/2))
in the all-0 background, Lx 1024 x Ly 4096, T 6e4, `--stop-absorbed`, profiles every 2000 steps
plus 10 per decade; `mfront.py` tracks the two strip edges (half of the interior Moore-interface
density, fixed at its t = 1000 value, `v_edge`; or half of the *current* interior density,
`v_adapt`, which is unbiased when the interior slowly decays) and fits dW/dt / 2 over records
with t >= 2000, W in [0.1, 0.9] Ly and the interior above 60 % of its t = 1000 density. The
integrated width A = sum_y iface_y / (Lx rho_interior(t)) gives `v_int`. Pilot: one seed, T 2e4,
steps of 2e-3 over p_1 in [0.176, 0.196] (p_2 = 0.45) and [0.160, 0.182] (p_2 = 0.50)
(`pilot/table_*.txt`). Then 1e-3 steps with 3 seeds over +-0.01 (5 of the outer values with 2
seeds), 2.5e-4 steps with 3-8 seeds within 1e-3 of the crossing, and Lx 2048 strips at the two
points nearest the crossing. Full seed-pooled tables: `scan/table_p2_0.45.txt`,
`scan/table_p2_0.50.txt` (columns: v_edge, error over seeds, v_adapt, v_int, per-seed v, the
interior interface density at t = 1000 and its final value relative to that).

**p_2 = 0.45** (Lx 1024 x Ly 4096, T 6e4; v_adapt in rows/step per front, + = active invades)

| p_1 | v | seeds | strip interior at the end (fraction of its t = 1000 interface density) |
|---|---|---|---|
| 0.17600 | - | 2 | gone (0.13) |
| 0.17700 | - | 2 | gone (0.19) |
| 0.17800 | - | 2 | gone (0.16) |
| 0.17900 | - | 2 | gone (0.12) |
| 0.18000 | - | 2 | gone (0.15) |
| 0.18100 | - | 3 | gone (0.15) |
| 0.18200 | - | 3 | gone (0.15) |
| 0.18300 | -0.0016 | 3 | gone (0.15) |
| 0.18400 | -0.0038 +- 0.0006 | 3 | gone (0.13) |
| 0.18500 | -0.0063 +- 0.0016 | 8 | gone (0.33) |
| 0.18525 | -0.0028 +- 0.0010 | 8 | gone (0.60) |
| 0.18550 | +0.0019 +- 0.0003 | 8 | stationary (0.85) |
| 0.18575 | +0.0052 +- 0.0007 | 3 | stationary (0.93) |
| 0.18600 | +0.0090 +- 0.0004 | 3 | stationary (0.95) |
| 0.18700 | +0.0180 +- 0.0004 | 3 | stationary (0.99) |
| 0.18800 | +0.0262 +- 0.0005 | 3 | stationary (0.99) |
| 0.18900 | +0.0328 +- 0.0004 | 3 | stationary (0.99) |
| 0.19000 | +0.0384 +- 0.0003 | 3 | stationary (1.00) |
| 0.19100 | +0.0437 +- 0.0006 | 2 | stationary (1.00) |
| 0.19200 | +0.0516 +- 0.0004 | 2 | stationary (1.00) |
| 0.19300 | +0.0544 +- 0.0015 | 2 | stationary (1.00) |
| 0.19400 | +0.0606 +- 0.0010 | 2 | stationary (1.00) |
| 0.19500 | +0.0645 +- 0.0004 | 2 | stationary (1.00) |
| 0.19600 | +0.0678 +- 0.0003 | 2 | stationary (1.00) |

Lx 2048 x Ly 4096: p_1 = 0.18525: v = +0.0001 +- 0.0003 (3 seeds, decaying 0.64); p_1 = 0.18550: v = +0.0021 +- 0.0003 (3 seeds, stationary 0.84)

**p_2 = 0.50** (Lx 1024 x Ly 4096, T 6e4; v_adapt in rows/step per front, + = active invades)

| p_1 | v | seeds | strip interior at the end (fraction of its t = 1000 interface density) |
|---|---|---|---|
| 0.16000 | - | 2 | gone (0.12) |
| 0.16100 | - | 2 | gone (0.16) |
| 0.16200 | - | 2 | gone (0.15) |
| 0.16300 | - | 2 | gone (0.15) |
| 0.16400 | - | 2 | gone (0.16) |
| 0.16500 | - | 2 | gone (0.12) |
| 0.16600 | - | 3 | gone (0.14) |
| 0.16700 | - | 3 | gone (0.14) |
| 0.16800 | - | 3 | gone (0.14) |
| 0.16900 | -0.0063 | 3 | gone (0.18) |
| 0.17000 | -0.0039 +- 0.0019 | 3 | gone (0.24) |
| 0.17025 | -0.0065 +- 0.0017 | 8 | gone (0.36) |
| 0.17050 | -0.0013 +- 0.0003 | 8 | decaying (0.69) |
| 0.17075 | +0.0023 +- 0.0002 | 8 | stationary (0.90) |
| 0.17100 | +0.0059 +- 0.0001 | 3 | stationary (0.93) |
| 0.17200 | +0.0155 +- 0.0001 | 3 | stationary (0.99) |
| 0.17300 | +0.0217 +- 0.0005 | 3 | stationary (1.00) |
| 0.17400 | +0.0294 +- 0.0002 | 3 | stationary (1.00) |
| 0.17500 | +0.0352 +- 0.0003 | 3 | stationary (1.00) |
| 0.17600 | +0.0402 +- 0.0002 | 2 | stationary (1.00) |
| 0.17700 | +0.0452 +- 0.0002 | 2 | stationary (1.00) |
| 0.17800 | +0.0507 +- 0.0010 | 2 | stationary (0.99) |
| 0.17900 | +0.0549 +- 0.0009 | 2 | stationary (1.00) |
| 0.18000 | +0.0603 +- 0.0009 | 2 | stationary (1.00) |
| 0.18100 | +0.0645 +- 0.0008 | 2 | stationary (1.00) |
| 0.18200 | +0.0687 +- 0.0008 | 2 | stationary (1.00) |

Lx 2048 x Ly 4096: p_1 = 0.17050: v = -0.0006 +- 0.0003 (3 seeds, decaying 0.68); p_1 = 0.17075: v = +0.0022 +- 0.0001 (3 seeds, stationary 0.89)

Reading of the tables:
- For p_1 above the crossing the front has a definite velocity: seeds agree to +-5e-4 or better,
  v is linear in p_1 (slope dv/dp_1 ~ 9-16 per unit p_1 within 1e-3 of the crossing, flattening to
  ~6 further out), and the strip interior is stationary (`cen_end` ~ 0.93-1.00).
- v crosses zero at **p_1* = 0.18535 +- 0.0001 (p_2 = 0.45)** and **p_1* = 0.17057 +- 0.0001 (p_2 = 0.50)**
  (linear interpolation of v_adapt between the 2.5e-4 points on either side, 8 seeds at Lx 1024 and
  3 at Lx 2048: v_adapt = +0.0019(3) / +0.0021(3) at 0.1855 and -0.0028(10) / +0.0001(3) at 0.18525;
  +0.0023(2) / +0.0022(1) at 0.17075 and -0.0013(3) / -0.0006(3) at 0.1705).
- Below the crossing the strip does not retreat with a stationary interior: 2.5e-4 below p_1*
  the interior interface density decays slowly (to 0.6-0.7 of its t = 1000 value by 6e4) while the
  edge velocity is small and not robust (Lx 1024, 8 seeds: -0.0028 +- 0.0010 at 0.18525 and
  -0.0013 +- 0.0003 at 0.1705; Lx 2048, 3 seeds: +0.0001 +- 0.0003 and -0.0006 +- 0.0003); 5e-4 below,
  the interior decays to 0.2-0.35 within 6e4 steps and the strip dissolves rather than retreats;
  1e-3 and more below, the interior coarsens away within ~1e4 steps (`cen` 0.13-0.15), exactly as
  the bulk does (Sec. 4). The fixed-threshold `v_edge` is strongly biased there (the interior falls
  through the threshold) and is listed only for completeness.

## 4. Are both phases stable near the crossing?

**Absorbing state** (`reps/`): 2 x 2 seeds, 20000 replicas, L 2048, T 1e4 (`mreps.py`). <A>(t)/4
at t = 10, 30, 100, 300, 1e3, 3e3, 1e4 and P_surv(1e3), P_surv(1e4):

| p_2 | p_1 | <A>/4 at t = 10 | 30 | 100 | 300 | 1e3 | 3e3 | 1e4 | P_surv(1e3) | P_surv(1e4) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.45 | 0.183 | 0.860 | 0.746 | 0.577 | 0.346 | 0.119 | 0.000 | 0.000 | 0.0030 | 0.0000 |
| 0.45 | 0.184 | 0.878 | 0.793 | 0.668 | 0.464 | 0.211 | 0.022 | 0.000 | 0.0050 | 0.0000 |
| 0.45 | 0.185 | 0.896 | 0.842 | 0.768 | 0.632 | 0.486 | 0.345 | 0.060 | 0.0092 | 0.0003 |
| 0.45 | 0.18525 | 0.901 | 0.854 | 0.813 | 0.734 | 0.624 | 0.345 | 0.083 | 0.0106 | 0.0003 |
| 0.45 | 0.1855 | 0.907 | 0.868 | 0.818 | 0.798 | 0.852 | 0.787 | 0.653 | 0.0129 | 0.0013 |
| 0.45 | 0.18575 | 0.911 | 0.877 | 0.848 | 0.882 | 0.845 | 1.100 | 2.042 | 0.0126 | 0.0029 |
| 0.45 | 0.186 | 0.916 | 0.888 | 0.872 | 0.911 | 1.090 | 1.708 | 5.964 | 0.0148 | 0.0046 |
| 0.45 | 0.187 | 0.939 | 0.933 | 0.983 | 1.199 | 2.238 | 8.964 | 98.223 | 0.0238 | 0.0163 |
| 0.45 | 0.188 | 0.965 | 0.987 | 1.142 | 1.729 | 4.360 | 31.091 | 464 | 0.0348 | 0.0301 |
| 0.45 | 0.190 | 0.997 | 1.120 | 1.510 | 3.260 | 14.038 | 190 | 2853 | 0.0698 | 0.0686 |
| 0.50 | 0.167 | 0.719 | 0.557 | 0.370 | 0.156 | 0.020 | 0.000 | 0.000 | 0.0007 | 0.0000 |
| 0.50 | 0.168 | 0.736 | 0.587 | 0.431 | 0.230 | 0.052 | 0.005 | 0.000 | 0.0017 | 0.0000 |
| 0.50 | 0.169 | 0.754 | 0.631 | 0.506 | 0.309 | 0.133 | 0.028 | 0.000 | 0.0030 | 0.0000 |
| 0.50 | 0.170 | 0.772 | 0.674 | 0.602 | 0.446 | 0.263 | 0.083 | 0.003 | 0.0049 | 0.0001 |
| 0.50 | 0.171 | 0.789 | 0.718 | 0.696 | 0.580 | 0.536 | 0.599 | 1.119 | 0.0084 | 0.0013 |
| 0.50 | 0.172 | 0.808 | 0.755 | 0.768 | 0.884 | 1.260 | 3.674 | 31.566 | 0.0151 | 0.0081 |
| 0.50 | 0.173 | 0.831 | 0.803 | 0.895 | 1.156 | 2.475 | 14.056 | 194 | 0.0222 | 0.0170 |
| 0.50 | 0.174 | 0.860 | 0.872 | 1.055 | 1.664 | 5.265 | 44.237 | 677 | 0.0356 | 0.0323 |
| 0.50 | 0.176 | 0.901 | 0.989 | 1.382 | 2.960 | 14.810 | 197 | 2823 | 0.0598 | 0.0593 |

At the voter point <A> = 4 at all t. Here <A> decays from the first steps for every p_1 up to
0.187 (p_2 = 0.45) and 0.173 (p_2 = 0.50) and grows from the first steps at 0.188 / 0.174: the
mean population of a small seed decays, i.e. the absorbing state is linearly stable, up to
~1.5-2e-3 *above* the front reversal. In that band the late-time mean turns up again (the few
seeds that reach a supercritical droplet grow at the front speed; P_surv(1e4) = 0.003-0.016) and
the droplet runs (Sec. 5) show a critical radius. The time at which <A> turns from decay to
growth moves down toward p_1* as t grows (crossing of <A>(t) = 4: 0.1871 at t = 100, 0.1863 at
300, 0.1859 at 1e3, 0.1856 at 1e4), because the growing droplets dominate the mean once they
exist; the early-time crossing (t <= 300) is the linear-stability edge, theta_a ~ 104-105 deg.

**Active phase** (`bulk/`, `bulk_long/`, `bulk_L/`): random density 1/2 on L 2048, T 6e4; rho_bond
at t = 1e2, 1e3, 3e3, 1e4, 2e4, 3e4, 4e4, 5e4, 6e4 and the ratio rho(6e4)/rho(2e4):

| run | 1e2 | 1e3 | 3e3 | 1e4 | 2e4 | 3e4 | 4e4 | 5e4 | 6e4 | ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| p2_0.45_p1_0.183 | 0.2869 | 0.2032 | 0.1399 | 0.0719 | 0.0502 | 0.0415 | 0.0372 | 0.0352 | 0.0340 | 0.677 |
| p2_0.45_p1_0.184 | 0.2912 | 0.2215 | 0.1773 | 0.1018 | 0.0672 | 0.0531 | 0.0435 | 0.0372 | 0.0306 | 0.455 |
| p2_0.45_p1_0.18525 | 0.2965 | 0.2444 | 0.2234 | 0.2031 | 0.1830 | 0.1731 | 0.1632 | 0.1665 | 0.1533 | 0.838 |
| p2_0.45_p1_0.1855 | 0.2982 | 0.2466 | 0.2307 | 0.2179 | 0.2169 | 0.2126 | 0.2131 | 0.2101 | 0.2045 | 0.943 |
| p2_0.45_p1_0.18575 | 0.2987 | 0.2511 | 0.2381 | 0.2321 | 0.2304 | 0.2306 | 0.2272 | 0.2327 | 0.2302 | 0.999 |
| p2_0.45_p1_0.185 | 0.2951 | 0.2412 | 0.2121 | 0.1778 | 0.1549 | 0.1254 | 0.1067 | 0.0922 | 0.0841 | 0.543 |
| p2_0.45_p1_0.186 | 0.2995 | 0.2557 | 0.2464 | 0.2418 | 0.2436 | 0.2418 | 0.2430 | 0.2406 | 0.2453 | 1.007 |
| p2_0.45_p1_0.187 | 0.3049 | 0.2689 | 0.2671 | 0.2668 | 0.2672 | 0.2674 | 0.2667 | 0.2675 | 0.2677 | 1.002 |
| p2_0.45_p1_0.188 | 0.3098 | 0.2824 | 0.2817 | 0.2812 | 0.2815 | 0.2816 | 0.2812 | 0.2815 | 0.2815 | 1.000 |
| p2_0.45_p1_0.190 | 0.3161 | 0.3006 | 0.3011 | 0.3015 | 0.3011 | 0.3012 | 0.3011 | 0.3022 | 0.3011 | 1.000 |
| p2_0.50_p1_0.167 | 0.2872 | 0.1828 | 0.1085 | 0.0559 | 0.0387 | 0.0325 | 0.0300 | 0.0290 | 0.0249 | 0.643 |
| p2_0.50_p1_0.168 | 0.2930 | 0.2019 | 0.1303 | 0.0671 | 0.0503 | 0.0364 | 0.0284 | 0.0266 | 0.0252 | 0.500 |
| p2_0.50_p1_0.169 | 0.2975 | 0.2227 | 0.1695 | 0.0911 | 0.0621 | 0.0522 | 0.0470 | 0.0430 | 0.0388 | 0.625 |
| p2_0.50_p1_0.170 | 0.3015 | 0.2412 | 0.2132 | 0.1657 | 0.1229 | 0.0938 | 0.0787 | 0.0725 | 0.0631 | 0.514 |
| p2_0.50_p1_0.171 | 0.3058 | 0.2611 | 0.2483 | 0.2406 | 0.2443 | 0.2453 | 0.2416 | 0.2434 | 0.2426 | 0.993 |
| p2_0.50_p1_0.172 | 0.3105 | 0.2773 | 0.2722 | 0.2738 | 0.2730 | 0.2718 | 0.2739 | 0.2724 | 0.2730 | 1.000 |
| p2_0.50_p1_0.173 | 0.3152 | 0.2885 | 0.2885 | 0.2859 | 0.2888 | 0.2879 | 0.2860 | 0.2886 | 0.2878 | 0.997 |
| p2_0.50_p1_0.174 | 0.3194 | 0.2998 | 0.2991 | 0.2992 | 0.2994 | 0.2990 | 0.2972 | 0.2988 | 0.2987 | 0.998 |
| p2_0.50_p1_0.176 | 0.3278 | 0.3144 | 0.3149 | 0.3147 | 0.3143 | 0.3158 | 0.3162 | 0.3164 | 0.3147 | 1.002 |

Longer and larger/smaller lattices: L 1024 to T 3e5 (`bulk_long/`): rho_bond at t = 1e4, 1e5, 2e5,
3e5 = 0.2212, 0.2112, 0.2124, 0.2101 (p_1 = 0.1855, p_2 = 0.45: stationary at 0.210), 0.2452, 0.2465,
0.2398, 0.2447 (0.186), 0.2644, 0.2668, 0.2675, 0.2690 (0.187); at p_2 = 0.50: 0.2127, 0.1866, 0.1620,
0.1390 (0.1705: still decaying), 0.2415, 0.2426, 0.2477, 0.2453 (0.171: stationary), 0.2739, 0.2689,
0.2690, 0.2724 (0.172).

The active phase is stationary (ratio 1.00) for p_1 >= 0.18575 (p_2 = 0.45) and >= 0.171
(p_2 = 0.50), drifts slowly at 0.1855 (ratio 0.94) and decays at 0.18525 (0.84), 0.185 (0.54) and
below, where it coarsens like the clustering phase (rho ~ 0.03 by 6e4 at 0.183-0.184). Its
stationary density falls steeply toward the line (0.301, 0.282, 0.268, 0.245, 0.230, ~0.205 at
p_1 = 0.190, 0.188, 0.187, 0.186, 0.18575, 0.1855), with no plateau value at the crossing and no
metastable plateau below it: the decay 2.5e-4 below the crossing is a steady slow decline, not a
plateau followed by a collapse.

**Lifetime vs L** (`bulk_L/`, L 256 and 512, T 1e6, `--stop-absorbed`): just below the crossing the
active phase is absorbed sooner on smaller lattices: at p_1 = 0.18525 (p_2 = 0.45) at t = 4.8e4
(L 256), 1.3e5 (L 512), and it is still at 75 % of its density at 6e4 on L 2048; at 0.185: 4.3e4,
1.5e5, still 40 % at 6e4 on L 2048; at p_2 = 0.50: 0.17025: 6.1e4, 1.0e5; 0.1705: 1.6e5, 3.1e5.
Its density decays steadily from the start in every case (e.g. 0.18525, L 2048: 0.223, 0.203,
0.183, 0.173, 0.163, 0.153 at 3e3, 1e4, 2e4, 3e4, 4e4, 6e4, only slightly faster than the voter
point's own logarithmic coarsening, 0.202 -> 0.176 over 1e4 -> 6e4). A lifetime growing with L
and a monotone decay without a plateau are coarsening of an *unstable* active phase whose domains
have to reach the system size, not nucleation of the absorbing phase inside a metastable one (that
would be L-independent, or shorter at larger L, after a plateau). Just above the crossing the
phase is stable: at 0.1855 the L 256 lattice keeps rho = 0.22 flat until a fluctuation absorbs it
at 6e5, and at 0.17075 (p_2 = 0.50) L 256 is still active at 1e6 (rho 0.225).

## 5. Critical droplets on the invading side (`drop/`, `drop/table.txt`)

Active disks (random density 1/2) of radius R in the all-0 state, L 2048, T 2e4, up to 400
replicas per R (stopped once 100 are alive at T), `--reps` box stepping:

| p_2 | p_1 | v | P_surv(2e4) at R = 1.5, 2.5, 3.5, 5, 7, 10, 14, 20, 28 | R_c (P = 1/2) | R_c v |
|---|---|---|---|---|---|
| 0.45 | 0.1855 | +0.0028 | 0.000 0.003 0.000 0.010 0.030 0.028 0.055 0.158 0.280 (still falling) | > 28 (D/v = 70) | - |
| 0.45 | 0.186 | +0.0090 | 0.003 0.007 0.010 0.055 0.110 0.160 0.323 0.498 0.820 | 20 | 0.18 |
| 0.45 | 0.187 | +0.0180 | 0.007 0.030 0.068 0.175 0.309 0.493 0.685 0.943 1.000 | 10 | 0.18 |
| 0.50 | 0.171 | +0.0059 | 0.000 0.003 0.013 0.025 0.033 0.087 0.128 0.296 0.524 (still falling) | ~27-30 | 0.16-0.18 |
| 0.50 | 0.172 | +0.0155 | 0.007 0.020 0.045 0.120 0.225 0.405 0.588 0.862 0.980 | 12 | 0.19 |

For R below R_c the survival keeps falling between 1e4 and 2e4 (the disks shrink and die); above
it, it is flat (the disks grow at the front speed). R_c v = 0.16-0.19 at every point, i.e.
R_c = D/v with D ~ 0.19 (the lazy von Neumann walk has D = 1/5), the curvature-limited nucleus
dR/dt = v - D/R found for D&N (R_c v = 0.21-0.26 at p = 0.5, D = 0.20). Together with the decay of
the 2 x 2 seed's mean population this is Case 4A of the closure: on the invading side, within
~2e-3 of the line, the absorbing state is stable against small seeds and invaded only by droplets
larger than D/v.

## 6. Comparison with the closure (`campaign/theory/MOLOFSKY.md`)

Measured, in the angle theta from the +delta_1 axis at the voter point (p_1 = 0.2 + (p_2 - 0.4) cot theta):

| | p_2 = 0.45 (rho = 0.052) | p_2 = 0.50 (rho = 0.104) |
|---|---|---|
| front reversal p_1* | 0.18535 +- 0.0001 -> theta = 106.33 deg | 0.17057 +- 0.0001 -> 106.40 deg |
| absorbing state linearly stable below (early-time seed decay) | 0.187-0.1875 -> 104.0-104.6 deg | 0.1735-0.174 -> 104.6-104.8 deg |
| active phase stable above | 0.1855 (rho = 0.210 stationary 1e5-3e5 on L 1024) | 0.17075 (L 256 active at 1e6), 0.171 (0.245 stationary to 3e5) |
| active phase unstable (coarsens, lifetime grows with L) below | 0.18525 | 0.1705 |
| 4A width (theta_a - theta_front) | ~1.5-2e-3 in p_1 (1.7-2.4 deg) | ~3e-3 (1.7-1.9 deg) |
| 4B width (theta_front - theta_half) | < 2.5e-4 (< 0.3 deg) | < 2.5e-4 (< 0.15 deg) |

Closure (table in MOLOFSKY.md Sec. 4, mapped to p_1 at the two p_2; T = 1/rho is 19 and 10):

| T | p_1(theta_a) | p_1(theta_front) | p_1(theta_half) at p_2 = 0.45 | 4A | 4B | same at p_2 = 0.50: a / front / half | 4A | 4B |
|---|---|---|---|---|---|---|---|---|
| 8 | 0.18813 | 0.18388 | 0.18163 | 4.3e-3 | 2.3e-3 | 0.17627 / 0.16776 / 0.16327 | 8.5e-3 | 4.5e-3 |
| 16 | 0.18717 | 0.18440 | 0.18296 | 2.8e-3 | 1.4e-3 | 0.17434 / 0.16880 / 0.16592 | 5.5e-3 | 2.9e-3 |
| 32 | 0.18658 | 0.18477 | 0.18383 | 1.8e-3 | 0.9e-3 | 0.17317 / 0.16954 / 0.16766 | 3.6e-3 | 1.9e-3 |
| 64 | 0.18623 | 0.18501 | 0.18440 | 1.2e-3 | 0.6e-3 | 0.17246 / 0.17002 / 0.16880 | 2.4e-3 | 1.2e-3 |
| 256 | 0.18583 | 0.18526 | 0.18497 | 0.6e-3 | 0.3e-3 | 0.17167 / 0.17051 / 0.16994 | 1.2e-3 | 0.6e-3 |
| 1024 | 0.18563 | 0.18535 | 0.18521 | 0.3e-3 | 0.14e-3 | 0.17127 / 0.17070 / 0.17042 | 0.6e-3 | 0.3e-3 |
| inf | 0.18560 | 0.18555 | 0.18552 | 0 | 0 | 0.17119 / 0.17110 / 0.17104 | 0 | 0 |

- **Direction of the line: confirmed.** The two crossings and the voter point are collinear
  (106.33 and 106.40 deg, slope dp_2/dp_1 = -3.42 and -3.40, i.e. -3.41 +- 0.02), 0.26 deg from the closure's
  T -> inf line (106.12 deg, -3.46) and within 0.1 deg of the chord from the voter point to
  Molofsky's empirical top-edge point (0.024, 1) (106.4 deg, -3.41). Extended linearly the measured
  line meets p_2 = 1 at p_1 = 0.0235 (theirs: 0.024) and p_2 = 0 at p_1 = 0.318 (their phase
  separation at 0.31, ergodic at 0.35). The closure's finite-T theta_front at T ~ 1/rho (107.3-107.9
  deg) is 1-1.5 deg too high; the measured line matches theta_front at T ~ 300-500.
- **The window is one-sided and narrower.** The absorbing state is linearly stable up to ~1-2 deg
  above the reversal (theta_a ~ 104-105 deg against the closure's 104.4 deg at T = 16, 105.0 deg
  at T = 32: the 4A sector has about the width the closure gives at T ~ 1/rho), with critical
  droplets R_c = D/v. But the active phase loses stability at the reversal itself: no point with a
  stationary active slab and a retreating front was found at 2.5e-4 resolution (0.3 deg at
  p_2 = 0.45, 0.15 deg at 0.50), where the closure at T = 16-32 predicts a 4B sector of 1.4-2.9e-3
  (1.0-1.5 deg) and even at T = 256-1024 of 1.5-6e-4. So theta_half coincides with theta_front to
  within the resolution, and the closure's form f = h(u)[a + c u(1-u)], in which theta_front sits
  two thirds of the way from theta_a to theta_half, is not what the rule does: the measured
  reaction term is "more concave", the active state going unstable as soon as the front stops.
- **Window width versus the predicted 5e-3 at distance 1/16.** At rho = 0.052 (1/19) the whole
  bistable band (theta_a to theta_half) is ~2e-3 wide in p_1 (~1.7e-3 perpendicular to the line),
  all of it on the ergodic side of the reversal, against the closure's 4.2e-3 (perpendicular
  3.9e-3) at T = 16, with 1.4e-3 of it on the clustering side. At rho = 0.104: ~3e-3 (perpendicular
  2.8e-3) against 8.5e-3 (T = 8) or 5.5e-3 (T = 16).

## 7. Verdict

Molofsky's five-cell rule has a sharp transition line through its voter point in the direction
the lineage closure predicts, 106.3 deg from the +delta_1 axis (slope dp_2/dp_1 = -3.41 +- 0.02,
against the closure's asymptotic -3.46 and the chord to their empirical point (0.024, 1), -3.41),
located at p_1* = 0.18535 +- 0.0001 (p_2 = 0.45) and 0.17057 +- 0.0001 (p_2 = 0.50). On its ergodic
side the transition looks first order in the D&N sense: the front velocity is linear in p_1
through the line (dv/dp_1 ~ 13-14 per unit p_1 at the crossing), the active phase has a finite
stationary density there (rho_bond ~ 0.2), small seeds die on average while droplets larger than
R_c = D/v (D ~ 0.19) grow, and this metastable-absorbing band is ~2e-3 wide in p_1 at distance
0.05 from the voter point, about what the closure gives at T ~ 1/rho. But the other half of the
closure's window is missing: nowhere within 2.5e-4 of the line (0.15-0.3 deg) is there a stable
active slab with a retreating front. Below the line the active phase decays by coarsening (its
lifetime grows with L and its density declines from the start), so the active state's stability
edge coincides with the front reversal to within the resolution, where the closure puts it
1.4e-3 (p_2 = 0.45) or 2.9e-3 (p_2 = 0.50) further, i.e. a 4B sector of 1.0-1.5 deg at T = 16; even
its T = 1024 value (0.14-0.3e-3) is at the edge of what was resolved. Compared with the predicted
perpendicular width of ~5e-3 at distance 1/16, the measured bistable band at distance 1/19 is
~1.7e-3 wide, entirely on the ergodic side of the reversal. So the closure gets the line (to
0.3 deg) and the order of magnitude and side of the absorbing-state metastability right, but
overestimates the stability of the active state beyond the reversal: for this rule the "bistable
window" is the metastability of the absorbing state alone, with no metastable active branch at
the 2.5e-4 level. That is also consistent with what Molofsky et al. saw at O(1) distance (patches
nucleated by dyads, no multiple-outcome region): the critical nucleus is D/v, a few cells once
v ~ 0.05, and nothing happens on the clustering side.

Resolution reached: 2.5e-4 in p_1 (0.3 deg at p_2 = 0.45, 0.15 deg at 0.50) with 8 seeds, Lx 1024
(and 2048) x Ly 4096, T 6e4; velocities +-3e-4 at the crossing; active-phase stability from L 256
to 2048 over 6e4-1e6 steps; seed survival over 20000 replicas to 1e4; droplets 400 replicas to 2e4.

**Draft manuscript sentences** (for the Discussion, replacing the Molofsky paragraph of
MOLOFSKY.md Sec. 8; not applied):

> We have simulated the five-cell rule of Molofsky \emph{et al.}~\cite{molofsky1999} near its
> voter point with the protocol of Secs.~\ref{...}. It has a sharp transition line through the
> voter point in the direction the closure predicts, $dp_2/dp_1 = -3.41(2)$ against $-3.46$, which
> also extrapolates to their empirical boundary point $(0.024, 1)$. On the ergodic side the
> transition is discontinuous in the sense of this paper: the front velocity crosses zero linearly,
> the active phase keeps a finite density, and the absorbing state is metastable, invaded only by
> droplets larger than $R_c = D/v$, over a band about $2\times10^{-3}$ wide in $p_1$ at distance
> $0.05$ from the voter point, as the closure gives for $T\sim1/\rho$. The other half of the
> closure's window is absent: the active phase coarsens away as soon as the front stops, so no
> retreating front with a stable active phase exists within $2.5\times10^{-4}$ of the line, where
> the closure predicts a sector of $1.4\times10^{-3}$. The closure thus locates the transition line
> and the metastability of the absorbing state but overstates the stability of the active state
> for a five-cell neighbourhood; the difference from D\&N, whose active phase is metastable down to
> $0.9\,q^*$, lies in the weight of the five-lineage patterns, not in their sign.


## 8. Fine scan across the reversal and the density at the line (`scan_fine/`, `bulk_fine/`)

Tightening of the two open numbers of Secs. 3, 4 and 7: the bound on an active-state branch
(a stable active phase with a retreating front) beyond the front reversal, and the stationary
density of the active phase at the line. Same code (`dnsim` v1.4, selftest passed before the
runs) and the same observables; new scripts `fine_runs.py` (run lists `runs_fine.txt`,
`runs_fine_wide.txt`, `runs_fine_extra.txt`), `run_fine.sh` / `run_fine_extra.sh` (the chains:
12 niced single-thread jobs, which on the 8 cores with hyperthreading gave about 25 % more
throughput than 8), `mfine.py` (analysis: `strips`, `runs`, `bulk`, `bulk2`, `life`, `fit`),
`sec8_tables.py` (the strip tables below) and `plot_fine.py` (`fine.png`: v(p_1), interior
histories, rho(delta p_1) with the fits). Raw outputs: `scan_fine/` (per-run table `runs_table.txt`,
pooled `table.txt`, `wide/` for Lx 4096), `bulk_fine/` (all tables in `tables.txt`). Cost 18
core-hours (1 h 45 min wall). Throughout, delta p_1 = p_1 - p_1* with p_1* = 0.18535 (p_2 = 0.45) and
0.17057 (p_2 = 0.50), and theta = atan2(p_2 - 0.4, p_1 - 0.2): 1e-4 in p_1 is 0.106 deg on the first
cut and 0.053 deg on the second.

Protocol. (i) Strips as in Sec. 3 (`--init strip:0.5:0.5 --bg 0 --stop-absorbed --profile`, profiles
every 2000 steps plus 10 per decade) but Lx 2048 x Ly 4096, T 1e5, 4 seeds, p_1 in steps of 1e-4
from 0.18505 to 0.18555 (p_2 = 0.45) and 0.17027 to 0.17077 (p_2 = 0.50), plus the half-step points
0.18540 and 0.17062 (4 seeds) and Lx 4096 x Ly 4096 at the two points nearest each crossing
(2 seeds). Velocities as in Sec. 3 (adaptive threshold, records with t >= 2e3, W in [0.1, 0.9] Ly,
interior above 60 % of its t = 1000 value), plus the same fit restricted to t >= 3e4 and the
integrated-width velocity v_int; the strip interior (central Ly/8 rows) is followed to 1e5.
(ii) Bulk runs from random 1/2 at the same p_1: L 512 and 1024 (2 seeds each) and L 2048 (1 seed) to
3e5 steps, output every 1e4, `--stop-absorbed`; lifetimes on L 256 and 512, 4 seeds, T 1e6,
`--stop-absorbed`. (iii) Density runs at p_2 = 0.45: L 2048 to 3e5 at delta p_1 = 0.5, 1, 1.5, 2, 3e-4
(0.18540, 0.18545, 0.18550, 0.18555, 0.18565) and L 1024 to 1e6 (2 seeds) at 0.5, 1, 1.5, 3e-4, plus
one run each at 2.5e-5 (0.185375); the stationary density is the mean over the last third of the
run and its drift the slope over that window.

### 8.1 Planar fronts at 1e-4 steps (`scan_fine/`)

**p_2 = 0.45** (Lx 2048 x Ly 4096, T 1e5, 4 seeds; v_adapt in rows/step per front over records with t >= 2e3 and, separately, t >= 3e4, and v_int from the integrated width; + = active invades; interior = central interface density as a fraction of its t = 1000 value at t = 1e4 / 3e4 / 6e4 / 1e5; stationary = no decline beyond 0.03 between 3e4 and 1e5, gone = below 0.6 at 1e5, decaying otherwise)

| p_1 | p_1 - p_1* | v (t >= 2e3) | v (t >= 3e4) | v_int | seeds | strip interior at 1e4 / 3e4 / 6e4 / 1e5 (fraction of its t = 1000 interface density) |
|---|---|---|---|---|---|---|
| 0.18505 | -3.0e-4 | -0.0016(3) |    -     | -0.0019(12) | 4 | gone (0.78 / 0.57 / 0.39 / 0.26) |
| 0.18515 | -2.0e-4 | -0.0008(1) | +0.0003     | +0.0002(9) | 4 | gone (0.81 / 0.62 / 0.44 / 0.35) |
| 0.18525 | -1.0e-4 | -0.0006(4) | -0.0015(8) | +0.0001(1) | 4 | gone (0.83 / 0.73 / 0.64 / 0.52) |
| 0.18535 | +0.0e-4 | +0.0002(1) | -0.0008(1) | +0.0005(3) | 4 | decaying (0.84 / 0.78 / 0.75 / 0.70) |
| 0.18540 | +0.5e-4 | +0.0008(3) | +0.0006(4) | +0.0005(5) | 4 | stationary (0.87 / 0.81 / 0.81 / 0.79) |
| 0.18545 | +1.0e-4 | +0.0016(1) | +0.0014(2) | +0.0013(2) | 4 | stationary (0.89 / 0.83 / 0.83 / 0.83) |
| 0.18555 | +2.0e-4 | +0.0031(1) | +0.0030(1) | +0.0032(1) | 4 | stationary (0.89 / 0.88 / 0.88 / 0.86) |

Lx 4096 x Ly 4096, T 1e5: p_1 = 0.18525: v = -0.0004(4) / -0.0012(5) / +0.0000(6) (2 seeds, interior 0.84 / 0.72 / 0.62 / 0.48); p_1 = 0.18535: v = -0.0001(1) / -0.0006(4) / +0.0000(2) (2 seeds, interior 0.86 / 0.81 / 0.76 / 0.74)

**p_2 = 0.50** (Lx 2048 x Ly 4096, T 1e5, 4 seeds; v_adapt in rows/step per front over records with t >= 2e3 and, separately, t >= 3e4, and v_int from the integrated width; + = active invades; interior = central interface density as a fraction of its t = 1000 value at t = 1e4 / 3e4 / 6e4 / 1e5; stationary = no decline beyond 0.03 between 3e4 and 1e5, gone = below 0.6 at 1e5, decaying otherwise)

| p_1 | p_1 - p_1* | v (t >= 2e3) | v (t >= 3e4) | v_int | seeds | strip interior at 1e4 / 3e4 / 6e4 / 1e5 (fraction of its t = 1000 interface density) |
|---|---|---|---|---|---|---|
| 0.17027 | -3.0e-4 | -0.0008(4) |    -     | +0.0014(8) | 4 | gone (0.78 / 0.60 / 0.38 / 0.27) |
| 0.17037 | -2.0e-4 | -0.0012(4) | -0.0020(21) | -0.0016(7) | 4 | gone (0.83 / 0.67 / 0.46 / 0.35) |
| 0.17047 | -1.0e-4 | -0.0003(3) | -0.0005(7) | +0.0005(10) | 4 | gone (0.84 / 0.69 / 0.58 / 0.53) |
| 0.17057 | +0.0e-4 | +0.0000(3) | -0.0001(5) | +0.0004(2) | 4 | decaying (0.87 / 0.82 / 0.76 / 0.75) |
| 0.17062 | +0.5e-4 | +0.0009(1) | +0.0006(3) | +0.0007(2) | 4 | stationary (0.89 / 0.83 / 0.82 / 0.83) |
| 0.17067 | +1.0e-4 | +0.0012(2) | +0.0013(2) | +0.0012(1) | 4 | stationary (0.88 / 0.87 / 0.86 / 0.86) |
| 0.17077 | +2.0e-4 | +0.0030(2) | +0.0031(3) | +0.0032(2) | 4 | stationary (0.90 / 0.88 / 0.90 / 0.89) |

Lx 4096 x Ly 4096, T 1e5: p_1 = 0.17047: v = -0.0008(6) / -0.0012(0) / -0.0012(6) (2 seeds, interior 0.82 / 0.78 / 0.67 / 0.53); p_1 = 0.17057: v = -0.0002(2) / -0.0006(0) / -0.0003(1) (2 seeds, interior 0.88 / 0.83 / 0.80 / 0.79)

Reading:
- **Velocity crossing.** Where the interior is stationary (delta p_1 >= +0.5e-4) the three estimators
  agree and v is linear in p_1: a weighted line through the three stationary points gives slope
  dv/dp_1 = 15.1 and v = 0 at 0.18534(1) (p_2 = 0.45), slope 13.6 and 0.17056(2) (p_2 = 0.50). At the
  grid point itself, 0.18535 / 0.17057, v = +0.0002(1) / +0.0000(3) (Lx 2048) and -0.0001(1) /
  -0.0002(2) (Lx 4096); the late-window values there, -0.0008(1) / -0.0001(5) and -0.0006(4) /
  -0.0006(0), carry the slow decay of the interior (below), so the spread between estimators,
  +-0.0005, is the systematic error: **p_1* = 0.18535(3) and 0.17057(3)**, the Sec. 3 values with a
  third of the error. theta_front = 106.33(3) and 106.40(2) deg; slopes from the voter point
  -3.413(7) and -3.398(3) (mean -3.41(1)).
- **Below the reversal the strip dissolves, it does not retreat.** At delta p_1 = -1e-4 the interior
  is at 0.83, 0.73, 0.64, 0.52 of its t = 1000 value at 1e4, 3e4, 6e4, 1e5 (Lx 4096: 0.84, 0.72,
  0.62, 0.48), the same decline as the bulk at that point (L 2048: rho_bond 0.203, 0.175, 0.132,
  0.097, 0.074 at 1e4, 3e4, 1e5, 2e5, 3e5, Sec. 8.2), while the edge velocity is small and
  estimator-dependent (-0.0006(4), -0.0015(8), +0.0001(1); Lx 4096: -0.0004(4), -0.0012(5),
  +0.0000(6)): the edges do not move at the -0.0015 that a linear v(p_1) would give there, because
  the whole slab coarsens at once. At -2e-4 and -3e-4 the interior is at 0.35 and 0.26 by 1e5 and the
  fixed-threshold width has collapsed to zero.
- **At the reversal the interior still declines**, slowly (0.84, 0.78, 0.75, 0.70; Lx 4096: 0.86,
  0.81, 0.76, 0.74; p_2 = 0.50: 0.87, 0.82, 0.76, 0.75 and 0.88, 0.83, 0.80, 0.79), at the relative
  rate of the bulk there (L 2048: 0.209 -> 0.137 from 1e4 to 3e5). **From +0.5e-4 on it is
  stationary** (0.81, 0.81, 0.79 at 3e4, 6e4, 1e5 at 0.18540; 0.83, 0.82, 0.83 at 0.17062) and v is
  positive in every estimator and every seed but one (+0.0008(3), +0.0006(4), +0.0005(5); +0.0009(1),
  +0.0006(3), +0.0007(2); per-seed values in `scan_fine/runs_table.txt`).
- So, to the question of Sec. 7: **no p_1 below the reversal has a stationary strip interior with a
  retreating front**, now at 1e-4 steps (0.5e-4 on the upper side). Wherever v < 0 the interior
  decays (to <= 0.53 of its early value by 1e5 steps), and wherever the interior is stationary v > 0.

### 8.2 Bulk at the same p_1: no plateau, lifetimes grow with L (`bulk_fine/`)

Random 1/2 on L 512 (2 seeds), 1024 (2) and 2048 (1), T 3e5, `--stop-absorbed`; rho_bond at fixed
times, seed-averaged, the ratio rho(3e5)/rho(1e5), and the absorption times (p_1 = 0.18538 is
0.185375):

| p_2 | p_1 | L | seeds | 1e+02 | 1e+03 | 3e+03 | 1e+04 | 3e+04 | 1e+05 | 2e+05 | 3e+05 | ratio 3e5/1e5 | absorbed at |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.45 | 0.18505 | 512 | 2 | 0.2946 | 0.2402 | 0.2196 | 0.1973 | 0.1479 | - | - | - | - | 1.6e+05, 7.9e+04 |
| 0.45 | 0.18505 | 1024 | 2 | 0.2954 | 0.2387 | 0.2138 | 0.1878 | 0.1363 | 0.0781 | 0.0654 | 0.0434 | 0.56 | - |
| 0.45 | 0.18505 | 2048 | 1 | 0.2950 | 0.2422 | 0.2160 | 0.1821 | 0.1344 | 0.0619 | 0.0424 | 0.0366 | 0.59 | - |
| 0.45 | 0.18515 | 512 | 2 | 0.2957 | 0.2404 | 0.2218 | 0.2053 | 0.1799 | 0.0953 | - | - | - | 1.1e+05, 2.9e+05 |
| 0.45 | 0.18515 | 1024 | 2 | 0.2958 | 0.2423 | 0.2224 | 0.1915 | 0.1568 | 0.0912 | 0.0596 | 0.0258 | 0.28 | - |
| 0.45 | 0.18515 | 2048 | 1 | 0.2955 | 0.2422 | 0.2205 | 0.1899 | 0.1491 | 0.0838 | 0.0464 | 0.0397 | 0.47 | - |
| 0.45 | 0.18525 | 512 | 2 | 0.2967 | 0.2445 | 0.2319 | 0.2030 | 0.1662 | 0.0324 | - | - | - | 1.3e+05, 1.0e+05 |
| 0.45 | 0.18525 | 1024 | 2 | 0.2965 | 0.2437 | 0.2245 | 0.2052 | 0.1699 | 0.0780 | 0.0432 | - | - | 2.8e+05 |
| 0.45 | 0.18525 | 2048 | 1 | 0.2965 | 0.2444 | 0.2234 | 0.2031 | 0.1746 | 0.1324 | 0.0974 | 0.0737 | 0.56 | - |
| 0.45 | 0.18535 | 512 | 2 | 0.2958 | 0.2429 | 0.2216 | 0.2023 | 0.1855 | 0.1924 | 0.1417 | 0.1181 | 0.61 | - |
| 0.45 | 0.18535 | 1024 | 2 | 0.2978 | 0.2455 | 0.2246 | 0.2050 | 0.1967 | 0.1784 | 0.1738 | 0.1653 | 0.93 | - |
| 0.45 | 0.18535 | 2048 | 1 | 0.2969 | 0.2440 | 0.2278 | 0.2086 | 0.1970 | 0.1734 | 0.1524 | 0.1369 | 0.79 | - |
| 0.45 | 0.18538 | 2048 | 1 | 0.2969 | 0.2469 | 0.2277 | 0.2096 | 0.2054 | 0.1864 | 0.1712 | 0.1731 | 0.93 | - |
| 0.45 | 0.18540 | 2048 | 1 | 0.2970 | 0.2479 | 0.2258 | 0.2113 | 0.2041 | 0.1900 | 0.1906 | 0.1887 | 0.99 | - |
| 0.45 | 0.18545 | 512 | 2 | 0.2967 | 0.2443 | 0.2239 | 0.2101 | 0.2220 | 0.1971 | 0.2037 | 0.1938 | 0.98 | - |
| 0.45 | 0.18545 | 1024 | 2 | 0.2989 | 0.2476 | 0.2335 | 0.2182 | 0.2000 | 0.2058 | 0.2045 | 0.2092 | 1.02 | - |
| 0.45 | 0.18545 | 2048 | 1 | 0.2973 | 0.2475 | 0.2307 | 0.2176 | 0.2139 | 0.2064 | 0.2035 | 0.2072 | 1.00 | - |
| 0.45 | 0.18550 | 2048 | 1 | 0.2982 | 0.2466 | 0.2307 | 0.2179 | 0.2167 | 0.2107 | 0.2092 | 0.2155 | 1.02 | - |
| 0.45 | 0.18555 | 512 | 2 | 0.2979 | 0.2471 | 0.2312 | 0.2203 | 0.2113 | 0.2290 | 0.2194 | 0.2149 | 0.94 | - |
| 0.45 | 0.18555 | 1024 | 2 | 0.2998 | 0.2492 | 0.2347 | 0.2242 | 0.2163 | 0.2092 | 0.2121 | 0.2191 | 1.05 | - |
| 0.45 | 0.18555 | 2048 | 1 | 0.2986 | 0.2486 | 0.2338 | 0.2237 | 0.2187 | 0.2165 | 0.2164 | 0.2136 | 0.99 | - |
| 0.45 | 0.18565 | 2048 | 1 | 0.2983 | 0.2510 | 0.2371 | 0.2269 | 0.2257 | 0.2240 | 0.2232 | 0.2253 | 1.01 | - |
| 0.50 | 0.17027 | 512 | 2 | 0.3024 | 0.2431 | 0.2282 | 0.1593 | 0.1296 | - | - | - | - | 7.6e+04, 1.2e+05 |
| 0.50 | 0.17027 | 1024 | 2 | 0.3030 | 0.2474 | 0.2262 | 0.1979 | 0.1433 | 0.0787 | 0.0642 | 0.0636 | 0.81 | - |
| 0.50 | 0.17027 | 2048 | 1 | 0.3030 | 0.2492 | 0.2232 | 0.1983 | 0.1443 | 0.0589 | 0.0357 | 0.0151 | 0.26 | - |
| 0.50 | 0.17037 | 512 | 2 | 0.3027 | 0.2480 | 0.2305 | 0.1955 | 0.1915 | 0.0543 | - | - | - | 1.2e+05, 1.3e+05 |
| 0.50 | 0.17037 | 1024 | 2 | 0.3031 | 0.2526 | 0.2280 | 0.2025 | 0.1720 | 0.0718 | - | - | - | 1.5e+05 |
| 0.50 | 0.17037 | 2048 | 1 | 0.3033 | 0.2507 | 0.2284 | 0.2010 | 0.1566 | 0.0942 | 0.0646 | 0.0627 | 0.67 | - |
| 0.50 | 0.17047 | 512 | 2 | 0.3029 | 0.2511 | 0.2215 | 0.2009 | 0.2034 | 0.1709 | - | - | - | 1.4e+05, 2.5e+05 |
| 0.50 | 0.17047 | 1024 | 2 | 0.3045 | 0.2521 | 0.2329 | 0.2058 | 0.1888 | 0.1297 | 0.1028 | 0.0918 | 0.71 | - |
| 0.50 | 0.17047 | 2048 | 1 | 0.3045 | 0.2511 | 0.2309 | 0.2112 | 0.1841 | 0.1025 | 0.0591 | 0.0525 | 0.51 | - |
| 0.50 | 0.17057 | 512 | 2 | 0.3022 | 0.2573 | 0.2320 | 0.2289 | 0.2167 | 0.1882 | 0.1771 | 0.1511 | 0.80 | - |
| 0.50 | 0.17057 | 1024 | 2 | 0.3046 | 0.2514 | 0.2334 | 0.2161 | 0.1965 | 0.1825 | 0.1450 | 0.1004 | 0.55 | - |
| 0.50 | 0.17057 | 2048 | 1 | 0.3044 | 0.2524 | 0.2347 | 0.2154 | 0.2060 | 0.1955 | 0.1841 | 0.1667 | 0.85 | - |
| 0.50 | 0.17062 | 2048 | 1 | 0.3046 | 0.2539 | 0.2377 | 0.2224 | 0.2127 | 0.1986 | 0.2052 | 0.1997 | 1.01 | - |
| 0.50 | 0.17067 | 512 | 2 | 0.3051 | 0.2571 | 0.2428 | 0.2398 | 0.2170 | 0.2174 | 0.2023 | 0.2341 | 1.08 | - |
| 0.50 | 0.17067 | 1024 | 2 | 0.3056 | 0.2523 | 0.2373 | 0.2205 | 0.2119 | 0.2213 | 0.2179 | 0.2132 | 0.96 | - |
| 0.50 | 0.17067 | 2048 | 1 | 0.3053 | 0.2546 | 0.2417 | 0.2276 | 0.2189 | 0.2167 | 0.2159 | 0.2105 | 0.97 | - |
| 0.50 | 0.17077 | 512 | 2 | 0.3059 | 0.2521 | 0.2426 | 0.2377 | 0.2243 | 0.2394 | 0.2328 | 0.2404 | 1.00 | - |
| 0.50 | 0.17077 | 1024 | 2 | 0.3059 | 0.2549 | 0.2436 | 0.2370 | 0.2323 | 0.2243 | 0.2307 | 0.2234 | 1.00 | - |
| 0.50 | 0.17077 | 2048 | 1 | 0.3057 | 0.2568 | 0.2420 | 0.2306 | 0.2290 | 0.2304 | 0.2290 | 0.2282 | 0.99 | - |

Lifetimes (L 256 and 512, T 1e6, `--stop-absorbed`; "+" = still active at 1e6; the last column is
the fraction of the lifetime over which rho stayed above 0.8 of its t = 1e4 value, ~1 for a plateau
followed by nucleation):

| p_2 | p_1 | L | lifetime per seed | mean | rho at 1e4 | rho at t_end/2 | fraction of lifetime with rho > 0.8 rho(1e4) |
|---|---|---|---|---|---|---|---|
| 0.45 | 0.18505 | 256 | 4.85e+04 5.76e+04 4.69e+04 3.92e+04 | 4.80e+04 | 0.197 | 0.120 | 0.36 |
| 0.45 | 0.18505 | 512 | 1.62e+05 7.86e+04 2.10e+05 8.72e+04 | 1.34e+05 | 0.187 | 0.122 | 0.31 |
| 0.45 | 0.18515 | 256 | 7.94e+04 6.07e+04 1.56e+05 5.11e+04 | 8.69e+04 | 0.208 | 0.130 | 0.37 |
| 0.45 | 0.18515 | 512 | 1.15e+05 2.91e+05 1.01e+05 1.06e+05 | 1.53e+05 | 0.201 | 0.111 | 0.25 |
| 0.45 | 0.18525 | 256 | 4.78e+04 6.84e+04 7.52e+04 1.41e+05 | 8.30e+04 | 0.199 | 0.156 | 0.48 |
| 0.45 | 0.18525 | 512 | 1.28e+05 1.02e+05 2.01e+05 1.63e+05 | 1.49e+05 | 0.205 | 0.136 | 0.36 |
| 0.45 | 0.18535 | 256 | 3.88e+05 2.98e+05 1.41e+05 1.39e+05 | 2.41e+05 | 0.218 | 0.185 | 0.61 |
| 0.45 | 0.18535 | 512 | 5.28e+05 3.42e+05 5.07e+05 1.00e+06+ | 4.59e+05 (3/4 absorbed) | 0.208 | 0.170 | 0.51 |
| 0.45 | 0.18540 | 256 | 1.79e+05 6.32e+05 2.57e+05 3.09e+05 | 3.44e+05 | 0.222 | 0.186 | 0.59 |
| 0.45 | 0.18540 | 512 | 5.41e+05 1.00e+06+ 1.00e+06+ 1.00e+06+ | 5.41e+05 (1/4 absorbed) | 0.209 | 0.190 | 0.77 |
| 0.45 | 0.18545 | 256 | 2.31e+05 4.58e+05 7.44e+05 1.00e+06+ | 4.78e+05 (3/4 absorbed) | 0.198 | 0.194 | 0.77 |
| 0.45 | 0.18545 | 512 | 1.00e+06+ 1.00e+06+ 1.00e+06+ 1.00e+06+ | > 1e+06 | 0.219 | 0.209 | 1.00 |
| 0.45 | 0.18555 | 256 | 1.00e+06+ 1.00e+06+ 1.00e+06+ 1.00e+06+ | > 1e+06 | 0.222 | 0.216 | 1.00 |
| 0.45 | 0.18555 | 512 | 1.00e+06+ 1.00e+06+ 1.00e+06+ 1.00e+06+ | > 1e+06 | 0.219 | 0.216 | 1.00 |
| 0.50 | 0.17027 | 256 | 5.58e+04 4.80e+04 6.08e+04 3.61e+04 | 5.02e+04 | 0.177 | 0.145 | 0.46 |
| 0.50 | 0.17027 | 512 | 7.60e+04 1.16e+05 6.92e+04 1.60e+05 | 1.05e+05 | 0.177 | 0.114 | 0.31 |
| 0.50 | 0.17037 | 256 | 8.62e+04 6.56e+04 4.94e+04 3.60e+04 | 5.93e+04 | 0.188 | 0.146 | 0.39 |
| 0.50 | 0.17037 | 512 | 1.24e+05 1.30e+05 7.71e+04 1.14e+05 | 1.11e+05 | 0.201 | 0.127 | 0.35 |
| 0.50 | 0.17047 | 256 | 2.41e+05 8.29e+04 1.03e+05 1.26e+05 | 1.38e+05 | 0.223 | 0.164 | 0.43 |
| 0.50 | 0.17047 | 512 | 1.38e+05 2.46e+05 1.28e+05 2.96e+05 | 2.02e+05 | 0.207 | 0.154 | 0.41 |
| 0.50 | 0.17057 | 256 | 2.24e+05 4.04e+05 1.19e+05 1.81e+05 | 2.32e+05 | 0.220 | 0.188 | 0.61 |
| 0.50 | 0.17057 | 512 | 7.16e+05 9.40e+05 1.00e+06+ 7.60e+05 | 8.05e+05 (3/4 absorbed) | 0.225 | 0.199 | 0.69 |
| 0.50 | 0.17062 | 256 | 4.66e+05 1.00e+06+ 1.08e+05 1.00e+06+ | 2.87e+05 (2/4 absorbed) | 0.225 | 0.193 | 0.81 |
| 0.50 | 0.17062 | 512 | 1.00e+06+ 1.00e+06+ 1.00e+06+ 1.00e+06+ | > 1e+06 | 0.224 | 0.207 | 0.95 |
| 0.50 | 0.17067 | 256 | 2.07e+05 8.17e+05 1.00e+06+ 2.38e+05 | 4.21e+05 (3/4 absorbed) | 0.227 | 0.224 | 0.67 |
| 0.50 | 0.17067 | 512 | 1.00e+06+ 1.00e+06+ 1.00e+06+ 1.00e+06+ | > 1e+06 | 0.235 | 0.214 | 0.95 |
| 0.50 | 0.17077 | 256 | 1.00e+06+ 1.00e+06+ 1.00e+06+ 1.00e+06+ | > 1e+06 | 0.238 | 0.233 | 1.00 |
| 0.50 | 0.17077 | 512 | 1.00e+06+ 1.00e+06+ 1.00e+06+ 1.00e+06+ | > 1e+06 | 0.234 | 0.231 | 1.00 |

Reading:
- **No plateau at or below the reversal.** On L 2048 the density declines from the start at every
  delta p_1 <= 0 (rho(3e5)/rho(1e5) = 0.59, 0.47, 0.56, 0.79 at -3, -2, -1, 0e-4 on the first cut;
  0.26, 0.67, 0.51, 0.85 on the second) and is stationary from +0.5e-4 on (0.99, 1.00, 1.02, 0.99,
  1.01 at +0.5, 1, 1.5, 2, 3e-4; 1.01, 0.97, 0.99 at +0.5, 1, 2e-4). At +2.5e-5 (0.185375) the L 2048
  run is still declining at 3e5 (0.210, 0.205, 0.186, 0.171, 0.173 at 1e4, 3e4, 1e5, 2e5, 3e5) and
  the L 1024 run coarsens away (0.177 at 1e5, 0.127 at 3e5, 0.079 at 5e5, 0.037 at 1e6).
- **Lifetimes grow with L** everywhere at or below the reversal: L 256 -> 512 gives x2.8, x1.8,
  x1.8, x1.9 at -3, -2, -1, 0e-4 (p_2 = 0.45) and x2.1, x1.9, x1.5, x3.5 (p_2 = 0.50), and the density
  spends only 0.25-0.48 of the lifetime above 0.8 of its t = 1e4 value below the line (0.5-0.7 at it):
  coarsening of an unstable phase, as in Sec. 4, not nucleation inside a metastable one (which would
  give an L-independent or decreasing lifetime after a plateau, as for D&N). From +0.5e-4 on, L 512
  survives to 1e6 in 7 of 8 runs (one absorbed at 5.4e5 at 0.18540) and the L 256 absorptions
  (1e5-7e5) come out of a stationary density (fraction 0.6-0.8): fluctuations of a small stable
  system.
- So the active phase is unstable at the reversal itself, on both cuts, and stable from 0.5e-4
  above it.

### 8.3 Bound on the active-state branch

The stability edge of the active phase, p_1(theta_half), lies between the reversal and the first
point above it on both cuts: 0 < p_1(theta_half) - p_1* <= 5e-5, with p_1* itself known to 3e-5. So
the edge coincides with the front reversal to within 5e-5 in p_1, and if anything lies on the
*invading* side of it (the active phase coarsens away while the front is still stationary), not
beyond it as the closure's 4B sector requires. In degrees: **|theta_half - theta_front| <= 0.05 deg
at p_2 = 0.45 and <= 0.03 deg at p_2 = 0.50**, with theta_half - theta_front <= 0 (previously
< 0.3 and < 0.15 deg from the 2.5e-4 grid). The closure's 4B width (MOLOFSKY.md Sec. 4, Sec. 6
above) is 1.49 deg at T = 16, 0.65 at T = 64 and 0.15 at T = 1024 (1.4e-3, 0.6e-3, 0.14e-3 in p_1 at
p_2 = 0.45; 2.9e-3, 1.2e-3, 0.3e-3 at p_2 = 0.50): the measured bound is 30 (50) times below the
T = 16 value and 3 (5) times below the T = 1024 value; with the width falling as ~T^-0.55 between
T = 64 and 1024, a sector narrower than 0.05 (0.03) deg needs T >~ 8e3 (2e4), far beyond the
T ~ 1/rho = 19 (10) at which the closure's theta_a and the absorbing-state band match (Sec. 6).

### 8.4 Density at the line (p_2 = 0.45)

Stationary rho_bond above the line (mean over the last third of each run, drift over that window
per 1e5 steps; the rows with L 2048 and T 3e5 or L 1024 and T 1e6 are new, the others from Sec. 4):

| p_1 | dp_1 (p_1* = 0.18535) | L | T | rho (last third) | drift per 1e5 | file |
|---|---|---|---|---|---|---|
| 0.18535 | +0.00000 | 2048 | 3e+05 | 0.1486 | -0.0208 | p2_0.45_p1_0.18535_L2048_s1.txt |
| 0.18538 | +0.00003 | 1024 | 1e+06 | 0.0765 | -0.0075 | p2_0.45_p1_0.185375_L1024_T1e6_s1.txt |
| 0.18538 | +0.00003 | 2048 | 3e+05 | 0.1660 | +0.0055 | p2_0.45_p1_0.185375_L2048_s1.txt |
| 0.18540 | +0.00005 | 1024 | 1e+06 | 0.1914 | +0.0032 | p2_0.45_p1_0.18540_L1024_T1e6_s1.txt |
| 0.18540 | +0.00005 | 1024 | 1e+06 | 0.1926 | -0.0062 | p2_0.45_p1_0.18540_L1024_T1e6_s2.txt |
| 0.18540 | +0.00005 | 2048 | 3e+05 | 0.1823 | -0.0010 | p2_0.45_p1_0.18540_L2048_s1.txt |
| 0.18545 | +0.00010 | 1024 | 1e+06 | 0.2020 | +0.0004 | p2_0.45_p1_0.18545_L1024_T1e6_s1.txt |
| 0.18545 | +0.00010 | 1024 | 1e+06 | 0.2047 | -0.0033 | p2_0.45_p1_0.18545_L1024_T1e6_s2.txt |
| 0.18545 | +0.00010 | 2048 | 3e+05 | 0.2044 | +0.0051 | p2_0.45_p1_0.18545_L2048_s1.txt |
| 0.18550 | +0.00015 | 1024 | 1e+06 | 0.2104 | -0.0017 | p2_0.45_p1_0.18550_L1024_T1e6_s1.txt |
| 0.18550 | +0.00015 | 1024 | 1e+06 | 0.2110 | -0.0020 | p2_0.45_p1_0.18550_L1024_T1e6_s2.txt |
| 0.18550 | +0.00015 | 1024 | 3e+05 | 0.2091 | +0.0005 | p2_0.45_p1_0.1855_L1024.txt |
| 0.18550 | +0.00015 | 2048 | 3e+05 | 0.2127 | +0.0077 | p2_0.45_p1_0.18550_L2048_s1.txt |
| 0.18550 | +0.00015 | 2048 | 6e+04 | 0.2101 | -0.0296 | p2_0.45_p1_0.1855_s1.txt |
| 0.18555 | +0.00020 | 2048 | 3e+05 | 0.2148 | -0.0032 | p2_0.45_p1_0.18555_L2048_s1.txt |
| 0.18565 | +0.00030 | 1024 | 1e+06 | 0.2246 | +0.0012 | p2_0.45_p1_0.18565_L1024_T1e6_s1.txt |
| 0.18565 | +0.00030 | 1024 | 1e+06 | 0.2261 | -0.0007 | p2_0.45_p1_0.18565_L1024_T1e6_s2.txt |
| 0.18565 | +0.00030 | 2048 | 3e+05 | 0.2252 | +0.0016 | p2_0.45_p1_0.18565_L2048_s1.txt |
| 0.18575 | +0.00040 | 2048 | 6e+04 | 0.2301 | +0.0066 | p2_0.45_p1_0.18575_s1.txt |
| 0.18600 | +0.00065 | 1024 | 3e+05 | 0.2416 | +0.0019 | p2_0.45_p1_0.186_L1024.txt |
| 0.18600 | +0.00065 | 2048 | 6e+04 | 0.2426 | +0.0098 | p2_0.45_p1_0.186_s1.txt |
| 0.18700 | +0.00165 | 1024 | 3e+05 | 0.2677 | -0.0000 | p2_0.45_p1_0.187_L1024.txt |
| 0.18700 | +0.00165 | 2048 | 6e+04 | 0.2674 | +0.0010 | p2_0.45_p1_0.187_s1.txt |
| 0.18800 | +0.00265 | 2048 | 6e+04 | 0.2814 | +0.0001 | p2_0.45_p1_0.188_s1.txt |
| 0.19000 | +0.00465 | 2048 | 6e+04 | 0.3013 | +0.0021 | p2_0.45_p1_0.190_s1.txt |

Time series of the L 1024 runs to 1e6:

| run | 1e+02 | 1e+03 | 1e+04 | 3e+04 | 1e+05 | 2e+05 | 3e+05 | 5e+05 | 1e+06 | drift over the last third (per 1e5 steps) | absorbed at |
|---|---|---|---|---|---|---|---|---|---|---|---|
| p2_0.45_p1_0.185375_L1024_T1e6_s1 | 0.2982 | 0.2465 | 0.2089 | 0.1765 | 0.1768 | 0.1742 | 0.1269 | 0.0787 | 0.0367 | -0.0075 | - |
| p2_0.45_p1_0.18540_L1024_T1e6_s1 | 0.2986 | 0.2463 | 0.2171 | 0.1994 | 0.1920 | 0.1735 | 0.1752 | 0.1603 | 0.2000 | +0.0032 | - |
| p2_0.45_p1_0.18540_L1024_T1e6_s2 | 0.2979 | 0.2452 | 0.2106 | 0.2067 | 0.1996 | 0.1796 | 0.1756 | 0.1868 | 0.1925 | -0.0062 | - |
| p2_0.45_p1_0.18545_L1024_T1e6_s1 | 0.2987 | 0.2438 | 0.2159 | 0.2047 | 0.2048 | 0.2122 | 0.2116 | 0.1938 | 0.1906 | +0.0004 | - |
| p2_0.45_p1_0.18545_L1024_T1e6_s2 | 0.2991 | 0.2515 | 0.2206 | 0.2033 | 0.2068 | 0.1969 | 0.2069 | 0.1994 | 0.1957 | -0.0033 | - |
| p2_0.45_p1_0.18550_L1024_T1e6_s1 | 0.2995 | 0.2483 | 0.2212 | 0.2170 | 0.2112 | 0.2131 | 0.2101 | 0.2048 | 0.2056 | -0.0017 | - |
| p2_0.45_p1_0.18550_L1024_T1e6_s2 | 0.2995 | 0.2489 | 0.2192 | 0.2160 | 0.2143 | 0.2173 | 0.2258 | 0.2072 | 0.2029 | -0.0020 | - |
| p2_0.45_p1_0.18565_L1024_T1e6_s1 | 0.2979 | 0.2550 | 0.2312 | 0.2257 | 0.2245 | 0.2280 | 0.2238 | 0.2262 | 0.2223 | +0.0012 | - |
| p2_0.45_p1_0.18565_L1024_T1e6_s2 | 0.2997 | 0.2485 | 0.2286 | 0.2226 | 0.2271 | 0.2293 | 0.2266 | 0.2259 | 0.2098 | -0.0007 | - |
| p2_0.50_p1_0.17062_L1024_T1e6_s1 | 0.3054 | 0.2543 | 0.2241 | 0.2101 | 0.2081 | 0.2036 | 0.2212 | 0.1977 | 0.1927 | +0.0002 | - |

Fits (`python3 mfine.py fit`; one point per p_1, the L 2048 run; 0.185375 excluded as not
stationary; rms = root-mean-square residual):

```
# excluded as not stationary (dp_1 < 0.5e-4): (0.00003, 0.1660 still decaying)
# points used (one per p_1: longest L 2048 run, seed-averaged; L 1024 if no L 2048 run): (0.00005, 0.1823, L 2048), (0.00010, 0.2044, L 2048), (0.00015, 0.2127, L 2048), (0.00020, 0.2148, L 2048), (0.00030, 0.2252, L 2048), (0.00040, 0.2301, L 2048), (0.00065, 0.2426, L 2048), (0.00165, 0.2674, L 2048), (0.00265, 0.2814, L 2048), (0.00465, 0.3013, L 2048)
# fits of rho(dp_1), dp_1 = p_1 - p_1*, p_1* = 0.18535 fixed
  dp_1 <= 1e-03: n = 7
    log   rho = A/ln(B/dp):        A = 2.0016  B = 2.236e+00   rms = 0.00329
    power rho = rho0 + C dp^beta:  rho0 = -0.2304  C = 0.683  beta = 0.05   rms = 0.00242
  dp_1 <= 5e-03: n = 10
    log   rho = A/ln(B/dp):        A = 2.2498  B = 7.313e+00   rms = 0.00336
    power rho = rho0 + C dp^beta:  rho0 = 0.0009  C = 0.525  beta = 0.11   rms = 0.00236

# p_1* free (grid 0.18515 .. 0.185397), points with p_1 - 0.18535 <= 1e-3
  log  : best p_1* = 0.18539  rms = 0.00138  A = 2.9320 B = 1.174e+02
  power: best p_1* = 0.18540  rms = 0.00108  rho0 = 0.1689 C = 0.701 beta = 0.30
```

Reading:
- The stationary density is 0.185(5), 0.204, 0.213, 0.215, 0.225 at delta p_1 = 0.5, 1, 1.5, 2, 3e-4
  (L 2048 over 1e5-3e5; the L 1024 runs to 1e6 agree within 0.005: 0.191/0.193, 0.202/0.205,
  0.210/0.211, 0.225/0.226), then 0.230, 0.243, 0.267, 0.281, 0.301 at 4e-4, 6.5e-4, 1.65e-3, 2.65e-3,
  4.65e-3 (Sec. 4). Drifts over the last third are within +-0.008 per 1e5 steps in every stationary
  run, but on L 2048 the density wanders by +-0.008 at delta p_1 = 0.5e-4 (0.175-0.192 between 2e5
  and 3e5), which sets the +-0.005 on the value there. On the second cut: 0.200 at +0.5e-4 (0.17062,
  L 2048 1e5-3e5; L 1024 to 1e6: 0.193), 0.21 at +1e-4, 0.228 at +2e-4.
- With p_1* fixed at the velocity reversal 0.18535 (7 points, delta p_1 <= 1e-3): A/ln(B/delta p_1)
  with A = 2.00, B = 2.24 fits to rms 0.0033 (it gives 0.187, 0.200, 0.208, 0.215, 0.224, 0.246 at
  0.5, 1, 1.5, 2, 3, 6.5e-4), and rho_0 + C delta p_1^beta to rms 0.0024 only with rho_0 = -0.23 and
  beta = 0.05, i.e. by imitating a logarithm; with rho_0 >= 0 imposed the best power fit is rho_0 ~ 0,
  beta = 0.11 (rms 0.0024 over delta p_1 <= 4.65e-3). So if the density line is the velocity
  reversal, the data prefer a density that vanishes there (logarithmically or as a small power), not
  a jump.
- With p_1* free, both forms push the line to the top of the allowed range, 0.18539 (log: A = 2.93,
  B = 117, rms 0.0014) and 0.18540 (power: rho_0 = 0.169, beta = 0.30, rms 0.0011), i.e. to where the
  active phase actually loses stability (Sec. 8.2), 0.4-0.5e-4 above the velocity reversal; there
  a jump of 0.17 with beta ~ 0.3 fits marginally better than a logarithmic vanishing.
- The rms differences (0.0033 vs 0.0024 at fixed p_1*, 0.0014 vs 0.0011 free) are below the
  scatter of the stationary densities near the line (+-0.003-0.005), so **the data cannot decide**
  between a logarithmic vanishing and a jump. What they establish: the active phase keeps a
  stationary density of 0.185(5) up to 0.5e-4 above the reversal, is still declining at 2.5e-5 above
  it (0.16-0.17 on L 2048 at 2e5-3e5, gone on L 1024 by 1e6) and coarsens away at the reversal; and
  the fall from 0.204 to 0.185 over the last 0.5e-4 is as large as the fall from 0.225 to 0.204
  over the preceding 2e-4, i.e. rho(delta p_1) steepens toward the line as both forms require,
  with the end of the stationary branch 0.5e-4 above the velocity reversal.

### 8.5 Numbers in `paper/main.tex` that would change (Appendix D and Table IV, `tab:molsim`)

Applied to `paper/main.tex` on 10-30 (see `paper/CHANGES_2026-10-30.md`, Addendum 3; the fine scan is the inset of Fig. 10(c), `paper/figs/fig_molofsky.py`). Appendix D (`app:molofsky`), paragraph "Simulation":
- Protocol sentence: "planar fronts from a strip of the active phase ($L_x = 1024$ and 2048,
  $L_y = 4096$, $6\times10^4$ steps, three to eight seeds, steps of $2.5\times10^{-4}$ in $p_1$
  near the reversal)" -> add "$L_x = 2048$ and 4096, $10^5$ steps, four seeds, steps of $10^{-4}$
  (and $0.5\times10^{-4}$) within $3\times10^{-4}$ of the reversal"; bulk runs "($L = 256$--2048,
  $6\times10^4$--$10^6$ steps)" still true (L 2048 now to $3\times10^5$, L 1024 to $10^6$).
- "$p_1^* = 0.18535(10)$ and $0.17057(10)$" -> "$0.18535(3)$ and $0.17057(3)$"; "slope $\approx 14$"
  -> "$\approx 15$ and $14$"; "$-3.41\pm0.02$" -> "$-3.41(1)$" ($-3.413(7)$ and $-3.398(3)$); the angles
  $106.33^\circ$ and $106.40^\circ$ and the extrapolations to $p_2 = 1$ and $0$ are unchanged.
- "The active phase is stationary down to $1.5\times10^{-4}$ above $p_1^*$, where its interface
  density is $0.21$, constant from $10^4$ to $3\times10^5$ steps at $L = 1024$" -> "stationary down to
  $0.5\times10^{-4}$ above $p_1^*$, where its interface density is $0.185(5)$, constant from $10^5$ to
  $3\times10^5$ steps at $L = 2048$ and to $10^6$ at $L = 1024$; $2.5\times10^{-5}$ above $p_1^*$ it
  is still declining at $3\times10^5$ steps, and at $p_1^*$ itself it coarsens away".
- "the density falls toward the line without a plateau, from $0.30$ at $5\times10^{-3}$ above it,
  roughly as $1/\ln(1/\Delta p_1)$, so whether it jumps at the line or vanishes there
  logarithmically ... our data do not decide" -> keep the verdict; "from $0.30$ at $5\times10^{-3}$
  to $0.19$ at $0.5\times10^{-4}$; $A/\ln(B/\Delta p_1)$ and $\rho_0 + C\,\Delta p_1^\beta$ fit
  equally well (rms $1$--$3\times10^{-3}$), and with $p_1^*$ free both put the end of the active branch
  $0.5\times10^{-4}$ above the reversal, where the stability test puts it".
- "Below $p_1^*$ the active phase has no metastable branch: $2.5\times10^{-4}$ below the reversal its
  interface density declines from the start ... and the time to absorption grows with $L$,
  $4.8\times10^4$ and $1.3\times10^5$ steps at $L = 256$ and 512 on the first cut and $1.6\times10^5$
  and $3.1\times10^5$ on the second" -> "at the reversal and $1$--$3\times10^{-4}$ below it ...
  $0.5$--$0.9\times10^5$ and $1.3$--$1.5\times10^5$ steps at $L = 256$ and 512 on the first cut,
  $0.5$--$1.4\times10^5$ and $1.1$--$2.0\times10^5$ on the second, and $2.4\times10^5$ / $4.6\times10^5$
  and $2.3\times10^5$ / $8\times10^5$ at $p_1^*$ itself".
- "No stationary active slab with a retreating front exists within $2.5\times10^{-4}$ of $p_1^*$, i.e.
  within $0.3^\circ$ and $0.15^\circ$ on the two cuts, where the closure at $T = 16$--$64$ puts the
  ray $f'(1/2) = 0$ some $0.6^\circ$--$1.5^\circ$ beyond the front" -> "within $0.5\times10^{-4}$ of
  $p_1^*$, i.e. within $0.05^\circ$ and $0.03^\circ$: the active phase loses stability between the
  reversal and $0.5\times10^{-4}$ above it, on the invading side, where the closure puts the ray
  $f'(1/2) = 0$ $1.5^\circ$ ($T = 16$), $0.65^\circ$ ($64$) and $0.15^\circ$ ($1024$) beyond the front".
- "while the front and the bound on $\theta_{1/2}$ want $T \gtrsim 300$" -> "$T \gtrsim 10^4$" (the
  bound on $\theta_{1/2}$ is below the closure's $T = 1024$ width).
- Main text, Discussion ("bound the active-state branch to less than a fifth of its predicted
  width at $T = 16$") -> "less than a thirtieth"; Introduction/abstract sentences ("within
  $10^{-3}$ of the transition line", "a band $2\times10^{-3}$ wide") unchanged.

Table IV (`tab:molsim`):
- "$p_1^*$: $0.18535(10)$ / $0.17057(10)$" -> "$0.18535(3)$ / $0.17057(3)$".
- "active phase stable for $\Delta p_1 \ge$: $1.5\times10^{-4}$ / $1.8\times10^{-4}$" ->
  "$0.5\times10^{-4}$ / $0.5\times10^{-4}$".
- "coarsens for $\Delta p_1 \le$: $-1.0\times10^{-4}$ / $-0.7\times10^{-4}$" -> "$2.5\times10^{-5}$ / $0$"
  (it coarsens at the reversal itself on both cuts, and $2.5\times10^{-5}$ above it on the first).
- "i.e. $\theta_{1/2} - \theta_{\rm front} <$: $0.3$ / $0.15$" -> "$|\theta_{1/2} - \theta_{\rm front}| \le$:
  $0.05$ / $0.03$" (with $\theta_{1/2} \le \theta_{\rm front}$); the closure row $1.49$, $0.65$, $0.15$ is
  now exceeded at every listed $T$.
- "$\rho$ at $\Delta p_1 = 1.5\times10^{-4}$: $0.21$ / $0.22$" -> add "$\rho$ at $\Delta p_1 =
  0.5\times10^{-4}$: $0.185(5)$ / $0.20$" (keep the $1.5\times10^{-4}$ row: $0.213$ / $0.22$).
- Fig. caption (c): the fine-scan points could replace or join the $2.5\times10^{-4}$ points
  (`fine.png` panel (a) has them; triangles there are $L_x = 4096$).

Sec. 7 above ("Resolution reached: 2.5e-4 ...", "nowhere within 2.5e-4 of the line", and the draft
manuscript sentences with $2.5\times10^{-4}$) is superseded by the numbers of this section:
resolution 0.5e-4 in p_1 (0.05 / 0.03 deg), velocities +-1-3e-4 at the crossing from 4 seeds at
Lx 2048 and 2 at Lx 4096 over 1e5 steps, active-phase stability from L 256 to 2048 over 3e5-1e6
steps at 1e-4 steps.

## 9. Seed-averaged density of the active phase at the line: jump or continuous vanishing? (`bulk_dens/`)

Question (continuing Sec. 8.4): on the cut p_2 = 0.45, where the front reverses at p_1* = 0.18535(3),
does the stationary interface density rho_bond of the active phase jump at the transition line or
vanish there continuously, logarithmically (A/ln(B/dp_1), as at a generalized-voter transition) or
as a small power? Sec. 8.4 could not decide from single runs: both forms fit to rms 0.002-0.003 with
p_1* fixed, and with p_1* free both put the end of the branch 0.4-0.5e-4 above the reversal, where
the active phase loses stability (Secs. 8.2-8.3). The target here was +-0.002 per point from seed
averages at the same dp_1 and below, within 50-80 core-hours. Throughout dp_1 = p_1 - 0.18535.

Short answer: the data do not decide, and they show why a density measurement cannot. The active
phase is stationary down to dp_1 = 0.5e-4, where its density is 0.188(1) at L = 2048 and 0.190(1) at
L = 4096; at 0.4e-4 it no longer settles (the L 2048 seeds wander between 0.15 and 0.19 for 1.2e6
steps, mean 0.173(3); L 4096 declines slowly to 0.176(3) by 9e5 steps), at 0.3e-4 it declines, at the
reversal it coarsens (Sec. 8). Along the stationary branch the density is linear in ln dp_1 (0.021(5) per
e-fold from 3e-4 down to 0.5e-4, no steepening), which a logarithm, a power dp_1^0.10 and a jump
with beta = 0.3 all reproduce. A continuous vanishing anchored at the central value of the
reversal is excluded by the two points below 0.5e-4, which lie 0.01-0.03 under every such form;
freeing p_1* moves the line to 0.18536-0.18539 for every form, within the error of the reversal and
onto the stability edge. So either the line is 1-4e-5 above 0.18535 and the density vanishes there
so steeply that it is still 0.17 one part in 1e5 above it, or the branch ends at finite density
0.17-0.19: a jump for every practical purpose, and the two cannot be told apart without the line
itself to better than 1e-5.

### 9.1 Protocol and cost

Runs (`bulk_dens/make_runs.py` writes the lists, chains `run_dens.sh`, `run_big.sh`, `run_ext.sh`,
`run_missing.sh`, `run_one.sh`, `run_ext2.sh`, log `run_dens.log`; invocation as in Sec. 8:
`--molofsky p1,0.45 --nbhd vn --init random:0.5 --stop-absorbed --every 2000 --ppd 10`, checkpoints
every 5e4 steps, 12-14 niced single-thread jobs on the 8-core machine, 1.0-1.2e9 cell updates/s per
job under load; dnsim v1.4, selftest passed):

| set | p_1 (dp_1 x 1e4) | L | seeds | T | purpose |
|---|---|---|---|---|---|
| (i) main | 0.18539, 0.18540, 0.18542, 0.18545, 0.18550, 0.18555, 0.18565 (0.4, 0.5, 0.7, 1, 1.5, 2, 3) | 2048 | 11-18 | 6e5 | seed-averaged rho(dp_1) |
| (iii) edge | 0.18538 (0.3) | 2048 | 11-14 | 6e5 | below the first stationary point of Sec. 8 |
| (ii) size | 0.18540, 0.18545 (0.5, 1) | 4096 | 11-12 | 3e5 -> 6e5 | does the density move with L |
| side | 0.18539 (0.4) | 4096 | 11-12 | 6e5 -> 9e5 | same, at the point where the L 2048 seeds split |
| extension | 0.18538, 0.18539, 0.18540 | 2048 | 11-18 (11-14 at 0.18538) | 6e5 -> 1.2e6 | are the near-edge points stationary or decaying |

Seeds 11-18 are new (Secs. 4 and 8 used 1-4). The continued runs were resumed in place from their
checkpoints with a larger --T (`--resume`; the output stays one continuous series on the 2000-step
grid, checked for duplicate or missing times; the wall time of the first segment is kept in
`bulk_dens/walltimes_main.txt`, the `# done` line of the output covers the last segment only).

Mishap: at 17:00 the main chain's list `runs.txt` was rewritten by a budget check that re-ran
`make_runs.py` in its then-default mode while xargs was still reading the list lazily (GNU xargs
4.10 reads 8 kB at a time), so 16 of the 64 main runs (seeds 17-18 at every p_1, seed 16 at
0.18555 and 0.18565) dropped out of the stream, and one command arrived cut off before its
`--out` (0.18555 seed 16): that run went to completion with its series on /dev/null and only its
checkpoint saved, which `--resume` then took for a finished run. The 15 skipped runs were run
afterwards from `runs_missing.txt` (`run_missing.sh`); the orphan checkpoint was deleted and the run
redone (`run_one.sh`). Every output's header was checked against its file name (p_1, p_2, L, seed,
`nbhd vn`, `init random:0.5`): no run was lost or duplicated, only delayed. `make_runs.py` without
arguments is now a dry run, and no chain's list is rewritten while it runs (`run_ext2.txt` was
written by hand once).

Cost: 65.6 job-hours (sum of the wall times of every run segment, `bulk_dens/cost.py`); wall clock
15:36-21:29 on 30 Oct 2026 on the 8 cores with hyperthreading (5.9 h x 8 = 47 core-hours; a
job-hour at 12 jobs on 8 cores is 2/3 of a core-hour). Within the 50-80 budget.

Analysis (`mdens.py`; every table of this section is in `bulk_dens/analysis.md`, regenerated by
`sh bulk_dens/analysis.sh`, which also redraws `dens.png` from `plot_dens.py`):
- Windows: the protocol value of a point is the mean over the last half of each run (3e5-6e5 for the
  6e5 runs, 6e5-1.2e6 for the continued ones, 4.5e5-9e5 for the L 4096 side runs), after checking
  that the two quarters of that half agree ("Q2 - Q1", with the drift = slope over the half per 1e5
  steps); the fixed windows 1.5e5-3e5, 3e5-6e5 and 6e5-1.2e6 are also tabulated so that L 2048 and
  L 4096, and early and late times, can be compared like with like. All on the 2000-step output grid.
- Errors: per run std x sqrt(2 tau / n) with tau the integrated autocorrelation time (Sokal
  windowing, c = 5); per point the larger of the seed-to-seed standard error and the combined
  per-run error, with a floor of 0.001. The Sec. 4 single runs (L 2048, 6e4 steps) carry +-0.003.
- Status of a point: 'absorbing' if a seed reached an absorbing state; 'split' if the seed means
  spread by more than 3 within-run standard deviations (the seeds are not sampling one state);
  'drifting' if the two quarters differ by more than 2 standard errors or the drift exceeds 2.5
  se; otherwise 'stationary'. Only stationary points enter the fits.
- Fits: weighted least squares; A/ln(B/dp_1); rho_0 + C dp_1^beta with rho_0 >= 0, and with
  rho_0 = 0; a jump, rho_0 + C dp_1^beta with beta fixed at 0.3 (and 0.5); p_1* fixed at 0.18535,
  at 0.18532 and 0.18538 (its error), and free on a 2e-6 grid up to 1e-6 below the lowest point;
  ranges: new points only (dp_1 <= 3e-4), new + Sec. 4 up to 1e-3, all points to 4.65e-3. The rho_0
  profile is the range of rho_0 with dchi^2 <= 1 (<= 4) when beta and C are refitted.

### 9.2 Run by run: within 0.5e-4 of the reversal the seeds stop sampling one state

Block means over 1e5 steps of every run with dp_1 <= 0.5e-4 (`mdens.py blocks 0.5e-4`;
`dens.png` (c) shows the same runs smoothed over 2e4 steps):

| run | status | 0-1e5 | 1-2e5 | 2-3e5 | 3-4e5 | 4-5e5 | 5-6e5 | 6-7e5 | 7-8e5 | 8-9e5 | 9-10e5 | 10-11e5 | 11-12e5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p2_0.45_p1_0.18538_L2048_s11 | done | 0.202 | 0.182 | 0.168 | 0.171 | 0.165 | 0.151 | 0.148 | 0.142 | 0.145 | 0.135 | 0.121 | 0.127 |
| p2_0.45_p1_0.18538_L2048_s12 | done | 0.201 | 0.186 | 0.174 | 0.163 | 0.178 | 0.175 | 0.182 | 0.189 | 0.182 | 0.179 | 0.182 | 0.180 |
| p2_0.45_p1_0.18538_L2048_s13 | done | 0.199 | 0.186 | 0.173 | 0.167 | 0.168 | 0.172 | 0.175 | 0.172 | 0.173 | 0.161 | 0.168 | 0.152 |
| p2_0.45_p1_0.18538_L2048_s14 | done | 0.203 | 0.182 | 0.153 | 0.144 | 0.142 | 0.137 | 0.120 | 0.139 | 0.138 | 0.135 | 0.128 | 0.134 |
| p2_0.45_p1_0.18539_L2048_s11 | done | 0.196 | 0.176 | 0.157 | 0.164 | 0.167 | 0.158 | 0.151 | 0.154 | 0.155 | 0.169 | 0.170 | 0.171 |
| p2_0.45_p1_0.18539_L2048_s12 | done | 0.202 | 0.183 | 0.187 | 0.176 | 0.168 | 0.171 | 0.175 | 0.174 | 0.173 | 0.163 | 0.170 | 0.183 |
| p2_0.45_p1_0.18539_L2048_s13 | done | 0.203 | 0.184 | 0.182 | 0.178 | 0.176 | 0.165 | 0.171 | 0.169 | 0.170 | 0.173 | 0.181 | 0.176 |
| p2_0.45_p1_0.18539_L2048_s14 | done | 0.207 | 0.183 | 0.172 | 0.166 | 0.161 | 0.154 | 0.150 | 0.152 | 0.158 | 0.169 | 0.169 | 0.170 |
| p2_0.45_p1_0.18539_L2048_s15 | done | 0.200 | 0.179 | 0.183 | 0.190 | 0.188 | 0.187 | 0.178 | 0.171 | 0.166 | 0.179 | 0.187 | 0.190 |
| p2_0.45_p1_0.18539_L2048_s16 | done | 0.201 | 0.187 | 0.184 | 0.190 | 0.188 | 0.186 | 0.178 | 0.170 | 0.177 | 0.174 | 0.176 | 0.174 |
| p2_0.45_p1_0.18539_L2048_s17 | done | 0.203 | 0.181 | 0.183 | 0.180 | 0.186 | 0.190 | 0.186 | 0.171 | 0.167 | 0.169 | 0.177 | 0.177 |
| p2_0.45_p1_0.18539_L2048_s18 | done | 0.202 | 0.180 | 0.176 | 0.175 | 0.179 | 0.180 | 0.182 | 0.184 | 0.188 | 0.180 | 0.187 | 0.186 |
| p2_0.45_p1_0.18539_L4096_s11 | done | 0.202 | 0.179 | 0.180 | 0.185 | 0.182 | 0.174 | 0.172 | 0.172 | 0.169 | - | - | - |
| p2_0.45_p1_0.18539_L4096_s12 | done | 0.205 | 0.187 | 0.184 | 0.184 | 0.178 | 0.178 | 0.180 | 0.181 | 0.178 | - | - | - |
| p2_0.45_p1_0.18540_L2048_s11 | done | 0.203 | 0.194 | 0.196 | 0.185 | 0.191 | 0.189 | 0.185 | 0.189 | 0.192 | 0.190 | 0.187 | 0.180 |
| p2_0.45_p1_0.18540_L2048_s12 | done | 0.207 | 0.194 | 0.184 | 0.181 | 0.173 | 0.176 | 0.179 | 0.185 | 0.183 | 0.186 | 0.189 | 0.193 |
| p2_0.45_p1_0.18540_L2048_s13 | done | 0.204 | 0.186 | 0.179 | 0.182 | 0.186 | 0.185 | 0.191 | 0.194 | 0.191 | 0.179 | 0.193 | 0.188 |
| p2_0.45_p1_0.18540_L2048_s14 | done | 0.207 | 0.183 | 0.173 | 0.172 | 0.184 | 0.191 | 0.185 | 0.195 | 0.188 | 0.193 | 0.187 | 0.170 |
| p2_0.45_p1_0.18540_L2048_s15 | done | 0.201 | 0.181 | 0.192 | 0.185 | 0.190 | 0.189 | 0.186 | 0.188 | 0.193 | 0.194 | 0.195 | 0.193 |
| p2_0.45_p1_0.18540_L2048_s16 | done | 0.203 | 0.177 | 0.173 | 0.175 | 0.181 | 0.184 | 0.188 | 0.198 | 0.186 | 0.189 | 0.187 | 0.188 |
| p2_0.45_p1_0.18540_L2048_s17 | done | 0.205 | 0.194 | 0.193 | 0.194 | 0.193 | 0.191 | 0.193 | 0.189 | 0.186 | 0.189 | 0.187 | 0.180 |
| p2_0.45_p1_0.18540_L2048_s18 | done | 0.207 | 0.190 | 0.189 | 0.191 | 0.181 | 0.178 | 0.179 | 0.183 | 0.195 | 0.197 | 0.196 | 0.190 |
| p2_0.45_p1_0.18540_L4096_s11 | done | 0.205 | 0.189 | 0.186 | 0.189 | 0.191 | 0.189 | - | - | - | - | - | - |
| p2_0.45_p1_0.18540_L4096_s12 | done | 0.205 | 0.191 | 0.192 | 0.186 | 0.191 | 0.191 | - | - | - | - | - | - |

- **0.18540 (0.5e-4): stationary at both sizes.** All eight L 2048 seeds stay between 0.17 and
  0.20 for 1.2e6 steps (block means 0.170-0.198), with excursions of +-0.01 lasting 1-2e5 steps;
  the seed average is 0.1845 over 3e5-6e5 and 0.1882 over 6e5-1.2e6 (the rise is the recovery of
  the three seeds that dipped to 0.172-0.181 around 3e5-5e5), and over the last half the eight seed
  means lie within 0.186-0.191. Both L 4096 runs sit at 0.186-0.192 from 1e5 steps on, 0.1896(10)
  over 3e5-6e5, drift +0.001(1) per 1e5, and agree with each other to 0.0001: the late L 2048
  window and L 4096 agree, 0.188(1) and 0.190(1).
- **0.18539 (0.4e-4): no longer one state.** The L 2048 seeds wander between 0.15 and 0.19 with
  excursions of +-0.02 lasting 2-4e5 steps: s11 falls from 0.196 to 0.151 (6-7e5) and climbs back to
  0.171 (11-12e5); s14 0.207 -> 0.150 -> 0.170; s15 holds 0.19 to 6e5, dips to 0.166 at 8-9e5 and
  returns to 0.190; s17 holds 0.18-0.19 to 7e5, then 0.167-0.177; s18 stays at 0.175-0.188
  throughout; none trends to zero within 1.2e6 steps. The seed average falls window by window,
  0.1785 (1.5e5-3e5) -> 0.1759 (3e5-6e5) -> 0.1727 (6e5-1.2e6), the eight means over the last half
  spread from 0.161 to 0.185 (sd 0.008, 3.3 within-run standard deviations: 'split'), and within
  that half the second quarter is above the first by 0.006(3). The two L 4096 runs agree with each
  other to 0.001 in the windows to 6e5 (self-averaging, as at 0.18540) and then drift apart as one
  of them keeps declining: blocks 0.185/0.184 (3-4e5) -> 0.182/0.178 -> 0.174/0.178 (5-6e5) ->
  0.172/0.180 -> 0.172/0.181 -> 0.169/0.178 (8-9e5); 0.1801(14) over 3e5-6e5 with drift -0.004(1)
  per 1e5, 0.1761(33) over the last half 4.5e5-9e5 (seed means 0.173 and 0.179, 0.006 apart where
  the L 4096 pairs at 0.18540 and 0.18545 agree to 0.0002). Over the matched window 3e5-6e5 the
  L 4096 density is 0.004 above the L 2048 one (0.180(1) vs 0.176(4)); over the last halves they
  agree (0.176(3) vs 0.173(3)), and both are below the 0.188 of the point above.
- **0.18538 (0.3e-4): declining.** Two of the four L 2048 seeds fall steadily, 0.20 -> 0.17
  (2-3e5) -> 0.15 (5-6e5) -> 0.127 and 0.134 (11-12e5); s13 reaches 0.15-0.17 (0.167 over the last
  half, still falling by 0.013 between its quarters); s12 dips to 0.163 at 3-4e5 and recovers to
  0.18. The seed average falls monotonically, 0.172 -> 0.161 -> 0.155 over the three windows and
  0.145 at 1.2e6 steps (`dens.png` (d)), with no sign of levelling off; no seed reaches an absorbing
  state by 1.2e6 steps (at L 1024 the Sec. 8.4 run at this p_1 was at 0.037 by 1e6: the smaller
  system empties faster, as the lifetimes of Sec. 8.2 do).
- From 0.18542 (0.7e-4) up every seed is stationary from 1.5e5 steps on, the seeds agree within
  0.001-0.004 and the quarters within their errors (table below).

So the active phase is stationary for dp_1 >= 0.5e-4 (to 1.2e6 steps at L 2048, 6e5 at L 4096),
marginal at 0.4e-4 (no decay at L 2048 within 1.2e6 steps, but +-0.02 excursions that do not
average out across seeds, and a slow decline at L 4096), and declining at 0.3e-4. If the decline
at 0.3e-4 is the loss of the phase and not a slow approach to a lower density, the stability edge of
Sec. 8.3 (0 < p_1(theta_half) - p_1* <= 5e-5) narrows to 3e-5 < p_1(theta_half) - p_1* <= 5e-5
(|theta_half - theta_front| between 0.03 and 0.05 deg), with the 0.4e-4 point on the edge; the edge
and the reversal, 0.18535(3), still coincide within the latter's error.

### 9.3 Seed-averaged density

The protocol value (mean over the last half of each run) with its error, the seed-to-seed sd, the
difference of the two quarters, the drift, the autocorrelation time, the means in the fixed windows
(number of runs covering the window in brackets), the status and the seed means
(`bulk_dens/sec9_tables.py`; the full tables per window are in `analysis.md`):

| p_1 | dp_1 (1e-4) | L | seeds | T | rho (last half) | err | seed sd | Q2 - Q1 | drift /1e5 | tau | 1.5e5-3e5 | 3e5-6e5 | 6e5-1.2e6 | status | seed means (last half, sorted) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.18538 | 0.3 | 2048 | 4 | 1.2e+06 | 0.1545 | 0.0121 | 0.0242 | -0.0087 +- 0.0039 | -0.0023 +- 0.0013 | 42837 | 0.1717 (4) | 0.1611 (4) | 0.1545 (4) | split | 0.132 0.136 0.167 0.182 |
| 0.18539 | 0.4 | 2048 | 8 | 1.2e+06 | 0.1727 | 0.0028 | 0.0078 | +0.0062 +- 0.0028 | +0.0018 +- 0.0007 | 41259 | 0.1785 (8) | 0.1759 (8) | 0.1727 (8) | split | 0.161 0.162 0.173 0.173 0.175 0.175 0.179 0.185 |
| 0.18539 | 0.4 | 4096 | 2 | 9e+05 | 0.1761 | 0.0033 | 0.0047 | -0.0015 +- 0.0026 | -0.0010 +- 0.0011 | 30326 | 0.1812 (2) | 0.1801 (2) | - | stationary | 0.173 0.179 |
| 0.18540 | 0.5 | 2048 | 8 | 1.2e+06 | 0.1882 | 0.0010 | 0.0020 | -0.0001 +- 0.0021 | -0.0001 +- 0.0007 | 37395 | 0.1853 (8) | 0.1845 (8) | 0.1882 (8) | stationary | 0.186 0.186 0.187 0.187 0.189 0.189 0.190 0.191 |
| 0.18540 | 0.5 | 4096 | 2 | 6e+05 | 0.1896 | 0.0010 | 0.0001 | +0.0019 +- 0.0022 | +0.0013 +- 0.0012 | 18633 | 0.1896 (2) | 0.1896 (2) | - | stationary | 0.190 0.190 |
| 0.18542 | 0.7 | 2048 | 8 | 6e+05 | 0.1941 | 0.0013 | 0.0038 | +0.0009 +- 0.0011 | +0.0003 +- 0.0007 | 13833 | 0.1954 (8) | 0.1941 (8) | - | stationary | 0.188 0.189 0.194 0.195 0.196 0.196 0.198 0.198 |
| 0.18545 | 1.0 | 2048 | 8 | 6e+05 | 0.2035 | 0.0010 | 0.0010 | +0.0007 +- 0.0010 | -0.0001 +- 0.0007 | 12790 | 0.2032 (8) | 0.2035 (8) | - | stationary | 0.202 0.202 0.203 0.203 0.204 0.204 0.204 0.205 |
| 0.18545 | 1.0 | 4096 | 2 | 6e+05 | 0.2017 | 0.0010 | 0.0002 | +0.0018 +- 0.0011 | +0.0011 +- 0.0003 | 19263 | 0.2040 (2) | 0.2017 (2) | - | drifting | 0.202 0.202 |
| 0.18550 | 1.5 | 2048 | 8 | 6e+05 | 0.2103 | 0.0010 | 0.0011 | +0.0001 +- 0.0004 | -0.0002 +- 0.0002 | 5321 | 0.2110 (8) | 0.2103 (8) | - | stationary | 0.208 0.209 0.209 0.210 0.211 0.211 0.211 0.212 |
| 0.18555 | 2.0 | 2048 | 8 | 6e+05 | 0.2174 | 0.0010 | 0.0004 | +0.0005 +- 0.0003 | +0.0002 +- 0.0002 | 7152 | 0.2173 (8) | 0.2174 (8) | - | stationary | 0.217 0.217 0.217 0.217 0.218 0.218 0.218 0.218 |
| 0.18565 | 3.0 | 2048 | 8 | 6e+05 | 0.2253 | 0.0010 | 0.0005 | -0.0000 +- 0.0003 | +0.0001 +- 0.0002 | 5431 | 0.2251 (8) | 0.2253 (8) | - | stationary | 0.225 0.225 0.225 0.225 0.225 0.226 0.226 0.226 |

- The +-0.002 target is met from 0.5e-4 up (errors 0.0010-0.0013, seed sd 0.0004-0.004); at 0.4e-4
  and 0.3e-4 the seed sd is 0.008 and 0.024 and the error means little, since the seeds are not
  sampling one state.
- Against the Sec. 8.4 single runs (3e5 steps): 0.1823 -> 0.188(1) at 0.18540 (the single run sat
  in a dip; the L 1024 runs of Sec. 8.4, 0.191-0.193 to 1e6 steps, were closer), 0.2044 -> 0.2035
  (0.18545), 0.2127 -> 0.2103 (0.18550), 0.2148 -> 0.2174 (0.18555), 0.2252 -> 0.2253 (0.18565): the
  single runs scattered by +-0.003 about the seed averages, as their quoted error said.
- L 2048 against L 4096 in the same window, 3e5-6e5: 0.1845(20) vs 0.1896(10) at 0.18540, which the
  late L 2048 window resolves (0.1882(10)); 0.2035(10) vs 0.2017(10) at 0.18545; 0.1759(44) vs
  0.1801(14) at 0.18539, and 0.1727(28) vs 0.1761(33) over the last halves. No size dependence
  beyond 0.002 where the phase is stationary.
- Sec. 4 points used in the fits (L 2048, 6e4 steps, +-0.003): 0.2301 (0.18575, dp_1 = 4e-4),
  0.2426 (0.186, 6.5e-4), 0.2674 (0.187, 1.65e-3), 0.2814 (0.188, 2.65e-3), 0.3013 (0.190, 4.65e-3).

### 9.4 Fits

Stationary points only: the six new points 0.5e-4 .. 3e-4 (L 2048, last half) and, for the 1e-3
range, the two Sec. 4 points at 4e-4 and 6.5e-4 (`mdens.py fit`; chi^2 per degree of freedom, rms
of the residuals in units of rho; full output with every parameter in `analysis.md`):

| p_1* | range (points) | A/ln(B/dp_1) | C dp_1^beta | rho_0 + C dp_1^beta, rho_0 >= 0 | rho_0 + C dp_1^0.3 | rho_0 + C dp_1^0.5 | rho_0 with dchi^2 <= 1 (<= 4) |
|---|---|---|---|---|---|---|---|
| 0.18535 | <= 3e-4 (6) | 1.31, rms 0.0010 (A 2.06, B 2.7) | 0.76, 0.0008 (beta 0.100) | 1.01, 0.0008 (rho_0 0.004(7), beta 0.10) | 2.19, 0.0012 (rho_0 0.137(2)) | 4.92, 0.0018 (rho_0 0.165(2)) | [0, 0.070] ([0, 0.124]) |
| 0.18535 | <= 1e-3 (8) | 1.18, 0.0015 (A 2.11, B 3.4) | 0.53, 0.0008 (beta 0.100) | 0.64, 0.0008 (rho_0 0, beta 0.10) | 2.31, 0.0024 (rho_0 0.140(2)) | 6.22, 0.0041 (rho_0 0.168(1)) | [0, 0.054] ([0, 0.106]) |
| 0.18535 | all (11) | 2.01, 0.0019 (A 2.31, B 9.1) | 0.55, 0.0013 (beta 0.102) | 0.54, 0.0011 (rho_0 0.027(4), beta 0.115) | 7.53, 0.0042 | 28, 0.0081 | [0, 0.060] ([0, 0.084]) |
| 0.18532 | <= 1e-3 (8) | 3.41, 0.0027 (A 1.74, B 0.74) | 1.92, 0.0019 (beta 0.122) | 2.30, 0.0019 (rho_0 0) | 4.24, 0.0032 (rho_0 0.125(3)) | 8.30, 0.0047 | - |
| 0.18538 | <= 1e-3 (8) | 0.93, 0.0012 (A 2.87, B 95) | 1.85, 0.0019 (beta 0.072) | 0.51, 0.0007 (rho_0 0.147(2), beta 0.25) | 0.53, 0.0010 (rho_0 0.157(2)) | 3.32, 0.0032 (rho_0 0.177(1)) | - |

p_1* free (2e-6 grid up to 0.185398; best value, dchi^2 <= 1 range, chi^2 at the best value and at
0.18535):

| form | <= 3e-4 (6 points) | <= 1e-3 (8 points) |
|---|---|---|
| A/ln(B/dp_1) | 0.185370 [0.185360, 0.185376], 2.2 vs 5.3 (A 2.53, B 21) | 0.185370 [0.185362, 0.185374], 2.2 vs 7.1 (A 2.53, B 21) |
| C dp_1^beta | 0.185362 [0.185350, 0.185370], 2.1 vs 3.0 (beta 0.090) | 0.185362 [0.185352, 0.185368], 2.2 vs 3.2 (beta 0.090) |
| rho_0 + C dp_1^beta | 0.185362 [0.185350, 0.185390], 2.1 vs 3.0 (rho_0 0.001(7)) | 0.185362 [0.185352, 0.185386], 2.1 vs 3.2 (rho_0 0.030(5), beta 0.105) |
| rho_0 + C dp_1^0.3 | 0.185382 [0.185374, 0.185388], 2.5 vs 8.8 (rho_0 0.158(2)) | 0.185384 [0.185378, 0.185388], 2.9 vs 13.9 (rho_0 0.160(2)) |

With the 0.4e-4 point added as if stationary, at 0.177(4) (an error wide enough to cover both
sizes; `points_edgevar.txt`): at p_1* = 0.18535 the logarithm gives chi^2/dof 1.66 (rms 0.0032),
the pure power 0.96 (0.0026), the beta = 0.3 jump 2.80 (0.0039) on the 1e-3 range; free, the
logarithm wants 0.185372 [0.185368, 0.185378] (3.8 vs 11.6 at 0.18535), the pure power 0.185366
[0.185356, 0.185370] (4.3 vs 6.7), the free-rho_0 power runs to the top of the grid as a jump
(rho_0 0.166(1), beta 0.32), the beta = 0.3 jump 0.185386 [0.185382, 0.185388] (3.4 vs 19.6).

- **With p_1* fixed at the reversal, 0.18535, the stationary branch does not discriminate.** The
  logarithm (chi^2/dof 1.2-1.3), the pure power with beta = 0.10 (0.5-0.8) and a jump of 0.14 with
  beta = 0.3 (2.2-2.3) all describe the six (eight) points at rms 0.0008-0.0024 against errors of
  0.001; only the beta = 0.5 jump fails (chi^2/dof 5-6). The rho_0 profile allows [0, 0.05-0.07] at
  dchi^2 <= 1 and up to 0.11-0.12 at dchi^2 <= 4. Over the whole range to 4.65e-3 the power
  (rho_0 = 0.027(4), beta 0.115) fits (chi^2/dof 0.5) and the logarithm strains (2.0), but the Sec. 4
  points at 1.65e-3 and above come from 6e4-step single runs. Sec. 8.4's best values, A = 2.0,
  B = 2.2 and beta = 0.05 .. 0.30, are what the seed averages give again (A = 2.06-2.11,
  B = 2.7-3.4; beta = 0.10).
- **Moving p_1* within its error changes the preferred form, not the quality.** At 0.18532 nothing
  fits well (chi^2/dof >= 1.9); at 0.18538 the jump wins (rho_0 = 0.15-0.16 with beta 0.25-0.3,
  chi^2/dof 0.5) and the logarithm still fits (0.9) with B = 95, i.e. a form that is 0.18 one part
  in 1e5 above the line.
- **p_1* free: every form puts the line above 0.18535.** On the stationary points alone the
  logarithm prefers 0.185370 [0.185360, 0.185376] (dchi^2 = 3.1-4.9 against 0.18535), the pure
  power 0.185362 [0.185350, 0.185370] (dchi^2 = 1), the beta = 0.3 jump 0.185382 [0.185374,
  0.185388] (dchi^2 = 6-11); the improvements are 1-3 sigma, and the power with free rho_0 is flat
  between 0.18535 and 0.18539. The grey line in `dens.png` (a) is the free logarithm.
- **The points below 0.5e-4 are what exclude a continuous vanishing at 0.18535.** Every form
  anchored there that fits the branch predicts 0.184-0.187 at 0.4e-4 and 0.179-0.183 at 0.3e-4
  (logarithm 0.186 / 0.181, power 0.185 / 0.179, beta = 0.3 jump 0.187 / 0.183). Measured: 0.173(3)
  at L 2048 and 0.176(3) at L 4096 over the last halves (0.4e-4, if the state is stationary at all;
  0.180(1) at L 4096 over 3e5-6e5), 0.155(12) and falling at 0.3e-4: 3-5 sigma below at 0.4e-4 at
  either size, 2 sigma at 0.3e-4 with its large error. Entered with the deliberately wide error above, the 0.4e-4 point alone pulls
  the free line to 0.18537(1) for the logarithm (dchi^2 = 6-8 against 0.18535) and 0.18537(1) for
  the pure power (dchi^2 = 2.2-2.4). If instead the state is not stationary there, the branch ends
  above it, which excludes the continuous forms anchored at 0.18535 just as well: they need a
  stationary phase all the way to the line.
- **Steepening.** The measured slope d rho / d ln dp_1 over the five segments from 0.5e-4 to 3e-4
  is 0.018, 0.026, 0.017, 0.025, 0.019 (+-0.005 each): constant, 0.021(5) per e-fold, with no
  steepening toward the line (`mdens.py slopes`). Neither continuous form requires one here: with
  A = 2.1 and B = 3.4 the logarithm's slope is A/ln^2(B/dp_1) = 0.016-0.019 across this range and
  steepens only below ~1e-5, and the power with beta = 0.10 has a constant slope in ln dp_1 by
  construction. The steepening that does appear is below 0.5e-4: 0.188 -> 0.173-0.180 (0.4e-4) ->
  0.155 (0.3e-4), i.e. 0.04-0.07 per e-fold, where the states are not stationary. So the
  "steepening toward the line" that Sec. 8.4 read into the single runs (0.204 -> 0.185 from 1e-4 to
  0.5e-4) is not in the seed averages on the stationary branch (0.2035 -> 0.188, the same slope as
  above), and where it occurs it is the branch giving way.

### 9.5 Verdict

The data do not decide between a jump and a logarithmic vanishing, and they show why a density
measurement cannot. The active phase keeps an interface density of 0.188(1) at the last point where
it is stationary, 0.5e-4 above the reversal (L 2048, 8 seeds, 6e5-1.2e6 steps; 0.190(1) at
L 4096), is at 0.17-0.18 and no longer settling at 0.4e-4 (L 2048 seeds wander between 0.15 and
0.19 for 1.2e6 steps, mean 0.173(3); L 4096 declines slowly to 0.176(3)), declines at 0.3e-4 (two
of four seeds to 0.13 by 1.2e6 steps, seed average 0.145 and falling) and coarsens at the reversal
(Sec. 8). Along the stationary branch, from 0.5e-4 to 3e-4, the density is linear in ln dp_1,
0.021(5) per e-fold, with no steepening, and a logarithm A/ln(B/dp_1), a power dp_1^0.10 and a jump
rho_0 + C dp_1^0.3 with rho_0 = 0.14 anchored at the reversal all fit it (chi^2/dof 1.2, 0.5, 2.3);
if the line is at 0.18535 and the branch continued to it, rho_0 <= 0.07 (dchi^2 <= 1), if it is at
0.18538, rho_0 = 0.15-0.16 is preferred. What the data exclude is a continuous vanishing anchored at
the central value of the reversal that continues through the points at 0.4 and 0.3e-4: every such
form predicts 0.18 there, 0.01-0.03 above what is measured, and with p_1* free every form moves the
line to 0.18536-0.18539, within the error of the reversal and onto the stability edge of Sec. 8. Two
readings remain, and they are not distinguishable here or by any affordable run: (a) the line is
1-4e-5 above 0.18535 and the density vanishes there continuously but so steeply (B = 20-100 in the
logarithm) that it is still 0.17 one part in 1e5 above it; (b) the active branch ends at finite
density, 0.17-0.19, between 0.3e-4 and 0.5e-4 above the reversal, a jump. The bound, in either
reading: wherever the active phase is stationary its interface density is >= 0.188(1) (0.5e-4 above
the reversal, both sizes), whatever happens to it happens within 5e-5 of the reversal, starting from
>= 0.17, and if the vanishing is continuous the line is at 0.18536-0.18539, not at 0.18535. Telling
(a) from (b) needs the line to better than 1e-5 (the front reversal gives 3e-5) and stationary
densities within 1e-5 of it, where these runs are not stationary at 1.2e6 steps. For the paper the
statement is unchanged in substance and sharper in numbers: the density is 0.188(1) at the last
stationary point, the branch is stable for dp_1 >= 0.5e-4, declines at 0.3e-4, and the two forms
remain indistinguishable, now with the bound above.

### 9.6 Sentences in `paper/main.tex` that the result would change

Applied to `paper/main.tex` on 10-30, late evening, by the paper session (see `paper/CHANGES_2026-10-30.md`, Addendum 4): the density paragraph of Appendix D, the protocol sentence, and the Table IV rows and caption. The list below is as handed over, before that edit.

Appendix D, paragraph beginning "On the side where the active phase invades, the transition has the
signatures of Secs.~\ref{sec:bulk} and \ref{sec:drops}.":
- "The active phase is stationary down to $0.5\times10^{-4}$ above $p_1^*$, where its interface
  density is $0.185(5)$, constant from $10^5$ to $3\times10^5$ steps at $L = 2048$ and to $10^6$
  steps at $L = 1024$; $2.5\times10^{-5}$ above $p_1^*$ it is still declining at $3\times10^5$
  steps, and at $p_1^*$ itself it coarsens away." -> "The active phase is stationary down to
  $0.5\times10^{-4}$ above $p_1^*$, where its interface density is $0.188(1)$ at $L = 2048$ (eight
  seeds, $6\times10^5$ to $1.2\times10^6$ steps) and $0.190(1)$ at $L = 4096$; $0.4\times10^{-4}$
  above $p_1^*$ it no longer settles (the $L = 2048$ seeds wander between $0.15$ and $0.19$ for
  $1.2\times10^6$ steps, $L = 4096$ declines slowly to $0.176(3)$ by $9\times10^5$ steps),
  $0.3\times10^{-4}$ above it
  declines (two of four seeds to $0.13$ by $1.2\times10^6$ steps), and at $p_1^*$ itself it coarsens
  away."
- "The density falls toward the line without a plateau, from $0.30$ at $5\times10^{-3}$ above it to
  $0.204$ at $10^{-4}$ and $0.185$ at $0.5\times10^{-4}$, steepening as the line is approached;
  $A/\ln(B/\Delta p_1)$ and $\rho_0 + C\,\Delta p_1^{\beta}$ fit it equally well, and with $p_1^*$
  left free both place the end of the active branch where the stability test places it,
  $0.5\times10^{-4}$ above the reversal." -> "The density falls toward the line without a plateau,
  from $0.30$ at $5\times10^{-3}$ above it to $0.204(1)$ at $10^{-4}$ and $0.188(1)$ at
  $0.5\times10^{-4}$, linearly in $\ln\Delta p_1$ ($0.021(5)$ per e-fold from $3\times10^{-4}$ down,
  with no steepening); $A/\ln(B/\Delta p_1)$, $C\,\Delta p_1^{0.10}$ and a jump
  $\rho_0 + C\,\Delta p_1^{0.3}$ with $\rho_0 = 0.14$ fit the stationary branch equally well, and
  with $p_1^*$ left free all of them place the line $1$-$4\times10^{-5}$ above the reversal, within
  its error, where the stability test places the end of the active branch; the seed-averaged
  densities $0.3$ and $0.4\times10^{-4}$ above the reversal lie $0.01$-$0.03$ below every continuous
  form anchored at it."
- "So whether the density jumps at the line or vanishes there logarithmically, as at a generalized
  voter transition, our data do not decide." -> "So whether the density jumps at the line or
  vanishes there logarithmically, as at a generalized voter transition, our data do not decide; they
  bound the room for a continuous vanishing to the last $5\times10^{-5}$ above the reversal, starting
  from a density of at least $0.17$, with the line then $1$-$4\times10^{-5}$ above the front
  reversal's central value."

Table IV (`tab:molsim`), p_2 = 0.45 column (the p_2 = 0.50 column is untouched):
- "active phase stable for $\Delta p_1\ge$ & $0.5\times10^{-4}$ & $0.5\times10^{-4}$" -> unchanged
  (now to $1.2\times10^6$ steps at $L = 2048$, eight seeds, and $6\times10^5$ at $L = 4096$).
- "\quad coarsens for $\Delta p_1\le$ & $0.25\times10^{-4}$ & $0$" -> "$0.3\times10^{-4}$" (declines
  without reaching an absorbing state by $1.2\times10^6$ steps at $L = 2048$; the $0.4\times10^{-4}$
  point is marginal).
- "\quad i.e. $|\theta_{1/2} - \theta_{\rm front}|\le$ & $0.05$ & $0.03$" -> unchanged as a bound
  ($0.03$-$0.05$ if the lower end is quoted).
- "$\rho$ at $\Delta p_1 = 0.5\times10^{-4}$ & $0.185(5)$ & $0.20$" -> "$0.188(1)$" ($0.190(1)$ at
  $L = 4096$).
- "$\rho$ at $\Delta p_1 = 1.5\times10^{-4}$ & $0.21$ & $0.22$" -> "$0.210(1)$".
- caption: "The active phase loses stability between the reversal and the first point above it, so
  the row for $|\theta_{1/2} - \theta_{\rm front}|$ is a bound." -> "The active phase loses stability
  between $0.3$ and $0.5\times10^{-4}$ above the reversal at $p_2 = 0.45$ (between the reversal and
  the first point above it at $p_2 = 0.50$), so the row for $|\theta_{1/2} - \theta_{\rm front}|$ is
  a bound."
