# Round 2 (2026-10-31): the experiments left open after the referee round

Everything here was designed and run on 2026-10-31 (afternoon) after Addendum 6; the aim is to close the
three items the changelog lists as "still open" (orientation test with direct front tracking; the 2x2-seed
initial condition of Sec. IV C; wetting with more seeds and Ly 2048) and the referee's two-droplet
coexistence test (report of 2026-10-29, M4), plus two cheap supporting measurements (roughness of the
wetting-layer edges, asked for in Sec. IV D itself; the correlation length of the active phase at p = 1/2,
report M5).  HPC: Torch, `cs` partition (128-core Xeon 8592+ nodes), account torch_pr_1342_general.
Torch's home directory is over its file quota (touch fails), so the binary, scripts, task lists, logs and
outputs all live under `/scratch/mo2273/dn/{bin2,campaign2,logs2,runs/R2_*}`; `/scratch` itself is at
99.5% of its space quota (k6-packing 89 GB, ssu ~4.9 TB, dn 5.6 GB), so the snapshot dumps are kept small.

## 1. Simulator: dnsim v1.5 (`dnsim.c` here; patch_v15.py applies it to the v1.4 source)

- `--pdir a,b`: the `--profile` is binned along the normal of a tilted band, s = (a y - b x) mod L
  (Lx = Ly = L; (1,0) is the row profile), L bins of L cells, 1/sqrt(a^2+b^2) apart in normal distance.
- `--snap file --snap-every n [--snap-from t]`: full-lattice records (u64 t, then the Lx/64 x Ly words,
  checkpoint layout) every n steps and at T; truncated to t <= t0 on resume like the profile.
- `--init hole:p0:R`: random:p0 outside a disk of radius R, the background inside (a dead droplet in the
  active phase).
- Checks: selftest passed (local gcc 15, Torch gcc 13 in the apptainer image, -march=x86-64-v3); the time
  series of an E16 wall run (v1.3, Torch) is reproduced bit for bit on both machines (md5 of the t <= 100
  lines); `--pdir 1,0` gives a profile identical to the row profile; on a (2,1) band of tilt.py the bin sums
  equal the time-series totals at every record and the band occupies bins [0, L/2); `hole:0.5:40` has the
  expected alive fraction; the snapshot at T equals `--save`'s checkpoint; a run interrupted at T = 55 and
  resumed to 100, and a run resumed from an older checkpoint over newer outputs, both reproduce the straight
  run's time series, profile and snapshots exactly (up to the extra record at the interruption time).

## 2. Experiments (make_tasks2.py; shard lists under runs/R2_*/tasks on Torch; camp2.slurm, 32 runs/shard)

| exp | runs | what | est core-h | jobs |
|---|---|---|---|---|
| R2_orient | 360 | bands (1,0), (2,1), (1,1) on L 4096 from tilt.py checkpoints (numpy seed 100000 + L + 100(3a+b) + s, 8 seeds), p = 0.1, 0.5, 0.9, q/q* = 0.99, 0.995, 1, 1.005, 1.01 (Table I q*), T 2e5, --pdir; snapshots (every 4e4 from 4e4) for the x = 1 bands | 177 | 21539742 (12 shards) |
| R2_orient2k | 120 | the same at L 2048, p = 0.5, T 1e5 | 7 | 21539743 (4) |
| R2_drop2 | 32 | p = 0.5, L 4096, T 2e5: disk:0.5:1000 in the dead state and hole:0.5:1000 in random:0.5, q/q* = 0.994, 0.997, 1.003, 1.006, 4 seeds, with row profiles | 16 | 21539744 |
| R2_dropbulk | 8 | bulk L 2048 at the four q, 2 seeds (iota(t) reference) | 1 | 21539746 |
| R2_wall | 52 | p = 0.5 walls at the E16 q values: Lx 1024 x Ly 2048 at x = 0.95 ... 0.995 (4 seeds, T 2e6; 3e6 at 0.995); Lx 1024 x Ly 4096 at x = 0.995, 0.9975 (4 seeds, 3e6); Lx 4096 x Ly 1024 at x = 0.98, 0.99, 0.995: seeds 3, 4 new to 2e6 and seeds 1, 2 continued from their E16 checkpoints (original seeds, so the continuation is what the run would have produced) to 2e6; Lx 4096 x Ly 2048 at 0.99 (2 seeds); Lx 256 x Ly 2048 at 0.99, 0.995 (4 seeds); bare wall q = 0 at Ly 2048 (2 seeds); snapshots every 1e5 from 1e5 | 25 | 21539747 (2 shards) + 21539785 (the 6 continued runs, which the first job skipped because their copied outputs still carried "# done") |

Local (this machine, 10 slots at nice 10, `runs_local.txt`, `run_local.log`):
- R2_seed22: 2x2 block, full (disk:1:0.75) and half-filled (disk:0.5:0.75, the E14 initial condition),
  L 1024, T 3e4, 4000 replicas, at the E14 q values x = 1.05, 1.1, 1.2, 1.4, 2.0 (of 0.0459).
- R2_corr: bulk L 2048 at x = 1.01, 1.02, 1.05, 1.1 (of 0.0459), 2 seeds, T 3e5, snapshots every 1e4 from
  1e5: two-point function and correlation length of the active phase.

## 3. Analysis scripts

orient2.py (direct front tracking on the tilted profiles: outer-envelope and integrated widths, normal
velocity by least squares for t >= t0, crossings with fit and seed-bootstrap errors), drop2.py (droplet and
hole radii from the alive fraction and from the row profile; growth-law prediction), wall2.py (stationarity
windows, Gamma(delta), fits with and without offset, psi = 1/3 alternative), rough.py + snap.py (per-column
edges of the wetting layer from the snapshots: width distribution, edge roughness sigma, height-difference
correlation C(r) and the local exponent), seed22.py, corr.py (FFT two-point function, chi, xi).
reduce.sh runs the Torch-side reductions into /scratch/mo2273/dn/figdata/R2/.

## 4. Results

(filled in as the runs finish)

### 4.1 R2_seed22 (local, done 15:45; seed22_table.txt)

| q/q* | full 2x2 block | half-filled block (E14 initial condition) |
|---|---|---|
| 1.05 | 0.0003(3) | 0.0000 |
| 1.10 | 0.0032(8) | 0.0003(3) |
| 1.20 | 0.011(2) | 0.0030(9) |
| 1.40 | 0.042(3) | 0.014(2) |
| 2.00 | 0.202(4) | 0.061(3) |

4000 replicas each, L 1024, T 3e4 (every survivor has escaped by t ~ 1e3).  The paper's 0.014 at 1.4 q*
was the half-filled block (E14, 1000 replicas, 14/1000), confirmed; its 0.046 at 2 q* (46/1000) is 2 sigma
below the new 0.061(3).  The full block, the initial condition the text named, survives three times as
often: 0.042(3) at 1.4 q* (the referee's own 0.044(3)) and 0.20 at 2 q*.  Applied to Sec. IV C and the
Fig. 4 caption; fig_drop.py now takes the R = 0.75 points from R2_seed22 (half) and draws the full block
as stars; R_c values unchanged (2.60 ... 41.8, D = 0.199, dR = 1.60).

### 4.2 R2_orient2k (Torch, done 15:52; figdata/R2/orient2k_*.txt|csv): L 2048, p = 0.5, 8 seeds, T 1e5

Outer-envelope tracking, fit t >= 1e4: x* = 0.9990(2), 0.9987(2), 0.9991(1) for (1,0), (2,1), (1,1)
(fit error; seed bootstrap the same), slopes q* dv/dq = 0.275, 0.285, 0.296, chi2/dof 9-19/3 (seed scatter
larger than the s.e.m. suggests at some q).  The three orientations agree to 0.04%; all three sit 0.1%
below the Table I strip value 0.04596 (v = +0.00035-0.0004 at x = 1 for every orientation).  Integrated
width: 0.9998(2), 0.9971(2), 0.9994(2) but chi2 = 84 for (2,1): the integrated estimator is the worse one
for tilted bands (the interior density of a tilted band is not that of the (1,0) band at the same time).

### 4.3 R2_corr (local; corr.py): preliminary, 18 snapshots per q

xi_exp = 15.2 (fit 4 <= r <= 40) / 18.0 (10 <= r <= 60) at x = 1.009; 12.7 / 14.9 at 1.019; 9.7 / 11.0 at
1.049; 7.0 / 11.0 at 1.099.  C(1) = 0.40, C(5) = 0.20, C(10) = 0.12, C(20) = 0.057, C(40) = 0.019 at
x = 1.009: a fast drop over the first lattice spacing (the active phase is half interface sites), then a
decay with xi ~ 15-18.  This is the direct measurement behind the paper's xi ~ (chi/2 pi)^1/2 ~ 15 (Sec. IV B)
and rules out xi ~ chi^1/2 ~ 40.

Final (16:10; 42 snapshots per q, 410 time-series samples; corr_table.txt): xi_exp = 15.0 / 17.8 (fits over
4-40 / 10-60) at x = 1.009, 12.7 / 14.8 at 1.019, 9.4 / 10.0 at 1.049, 7.0 / 9.0 at 1.099; second-moment
length xi_2 = (sum r^2 C / 4 sum C)^1/2 over r <= 100: 20.6, 16.9, 12.2, 8.0.  chi_L = L^2 <m^2> from the
time series: 498(106), 216(45), 118(26), 69(16); the truncated 2D sum of C(r) over r <= 200 gives 389, 237,
120, 61.  These agree with E13 (bulk05.py, L 2048, 4 seeds): chi = 464 at 1.008 q*, 234 at 1.018, 117 at
1.048, 63 at 1.098 (and 1080 at 0.998 q*, the value the paper rounds to "1300 at the transition").  With
C(1) = 0.40, chi = 2 pi C_0 xi^2 gives xi ~ 13-15 for chi ~ 460-500, consistent with the direct fits.

### 4.4 R2_orient (Torch, 360 runs done 16:38; figdata/R2/orient_*.txt|csv): L 4096, T 2e5, 8 seeds

Outer-envelope tracking on the normal profile, fit t >= 1e4.  The five-point lines have chi2/dof up to
70/3 because the 0.99 q* points retreat faster than linear: at delta = 0.01 the band interior
(8e6 sites) nucleates the absorbing phase within the 2e5-step window (bulk lifetime 7e4 at 0.988 q*,
Sec. IV B), most at (1,0) whose interior is the widest.  Crossings from the four points x >= 0.995
(fit error / seed bootstrap of 1000 resamples agree):

| p | (1,0) | (2,1) | (1,1) | spread | mean vs Table I q* | local slope q* dv/dq |
|---|---|---|---|---|---|---|
| 0.1 | 0.99711(35) | 0.99711(10) | 0.99711(22) | 0.000% | -0.29% | 0.17 |
| 0.5 | 0.99913(10) | 0.99886(11) | 0.99910(13) | 0.027% | -0.10% | 0.28 |
| 0.9 | 1.00529(26) | 1.00524(24) | 1.00474(15) | 0.055% | +0.51% | 0.18 |

Three-point fits (0.995-1.005) move these by <= 0.0005; the five-point fits by <= 0.0008 (p 0.1) and
<= 0.0003 otherwise; fit windows from 2e4 or 5e4 by <= 0.0005.  The integrated-width estimator agrees
for (1,0) and (1,1) but is biased for (2,1) (chi2 up to 220: the interior density of a tilted band is not
the (1,0) one).  L 2048 at p = 0.5 (Sec. 4.2): 0.9990, 0.9987, 0.9991, i.e. no size dependence beyond
0.02%.  So: orientation dependence < 1e-3 (was: < 4e-3 with the interface-count method, E15, whose
crossings sat 0.1-1.0% above these); the strip thresholds of Table I differ from the band values by
-0.1% (p 0.5), -0.3% (0.1), +0.5% (0.9), the systematic of the method.  The local slopes at q* are
steeper than Table I's secants over q* < q <= 1.065 q* (0.22 at p 0.5, 0.14 at 0.1, 0.11 at 0.9): v(q)
is concave above q*, as the strip campaign found at p >= 0.99.  Applied to Sec. IV A, Fig. 2 (new
fig_orient from campaign/round2/fig_orient_r2.py), Appendix E.

### 4.5 R2_drop2 / R2_dropbulk (Torch, done 16:38; figdata/R2/drop2.txt): two-droplet test, L 4096, T 2e5

dR/dt from the row profiles (time-series estimator agrees within its larger noise), fit t >= 1e4,
mean +- s.e.m. over 4 seeds (1e-3 sites/step); predicted = +-v - D/R with v = 0.278 (x - 0.999), D = 0.20:

| q/q* | active disk | predicted | dead hole | predicted |
|---|---|---|---|---|
| 0.994 | -1.80(53) | -1.59 | +2.20(17) | +1.19 |
| 0.997 | -1.08(13) | -0.76 | +0.72(10) | +0.36 |
| 1.003 | +0.87(5) | +0.91 | -1.22(12) | -1.31 |
| 1.006 | +1.73(6) | +1.75 | -1.86(7) | -2.15 |

At every q one droplet grows and the other shrinks: no generic-coexistence window at the 0.3% level
(both would shrink).  Above q* the rates match the growth law to 0.1e-3 (disk) and 0.1-0.3e-3 (hole);
below q* the hole opens 0.4-1.0e-3 faster than the law, with the largest seed scatter: the metastable
active background (delta = 0.003-0.006, bulk lifetime ~ 1e5-5e5 steps) starts to nucleate dead patches
within the run, which add to the dead area.  New Fig. (fig_drop2) and paragraph in Sec. IV A.

### 4.6 R2_wall (Torch, 52 runs done 17:30; figdata/R2/wall2.txt, figdata/R2_wall/*.csv, rough_*.csv, whist_*.csv)

Stationarity: at delta = 0.0065 the successive 1e6-step windows of the Ly 4096 runs (T 3e6) give seed means
of I = 93.0, 95.0, 93.4 (2% apart); at delta = 0.0040 117.7, 119.6, 118.7.  The E16 Lx 1024 x Ly 1024 runs
at 0.0065 relax from 96.4 ([3e5, 1e6]) to 94.5 ([5e5, 1e6]): the published Gamma = 88 was 2-3% high.
Bare wall I_0 = 8.71 at Ly 2048 (8.72 at 1024).

Gamma = I - I_0 (late half of each run; E16: [3e5, 1e6]; runs cut where the wall has drifted so far that
|alive - 1/2| >= 0.2, fig_wall05.py):

| delta | Lx 1024, Ly 1024 (E16) | Lx 1024, Ly 2048 | Lx 1024, Ly 4096 | others |
|---|---|---|---|---|
| 0.0514 | 27.64(15) | 27.40(9) | | Lx 256: 27.2; Lx 4096: 27.4 |
| 0.0315 | 36.78(11) | 37.05(20) | | 37.5 / 36.9 |
| 0.0215 | 45.87(22) | 45.50(31) | | 44.8 / (box too small) |
| 0.0115 | 64.87(46) | 64.38(1.03) | | Lx 256 Ly 2048: 60.6(4); Lx 4096 Ly 2048: 64.40(44) |
| 0.0065 | 87.68(1.03) | 88.56(2.24), n 3 | 85.33(48) | Lx 256 Ly 2048: 86.7(4.3) |
| 0.0040 | | | 109.13(1.79) | |

Lx and Ly independence to 3% (Lx 256 is 3-6% low at delta <= 0.0115: the wall wanders most there).  Local
exponent of the Lx 1024 series: 0.55 (0.021-0.011), 0.53 (0.011-0.0065, using 87.2), 0.46-0.51
(0.0065-0.0040, using 87.2 or 85.3).  Fits over the 7 points Lx 1024, Ly >= 2048 (delta <= 0.051): free
power law psi = 0.55; A delta^-1/3 + B needs B = -32 and misses the ends by 2-3% (26.7 vs 27.4 at 0.051;
105.6 vs 109.1 at 0.004).  The q* uncertainty (0.1%) is 25% of delta at 0.004.

Wandering walls: the wetted wall drifts along y with a persistent velocity of order 1e-4 rows/step in some
runs (alive fraction 0.52 -> 0.89 over 1.8e6 steps in Lx 1024 x Ly 2048 s2 at 0.0065, then absorbed into
the all-alive state; 0.52 -> 0.04 over 2e6 steps in Lx 4096 x Ly 1024 s1 at 0.011).  Absorbed or collapsed
before T: Lx 1024 x Ly 2048 s2 (0.0065), Lx 256 x Ly 2048 s4 (0.0065), and 4 of the 6 continued/new
Lx 4096 x Ly 1024 runs at delta <= 0.021 (the layer edges' excursions, rms 70-100 at Lx 4096, reach the
other wall 512 rows away).  Those points are dropped from Fig. 6; the Ly 4096 runs (delta 0.0065, 0.004)
are all stationary to 3e6.  The alive fraction is Z2-symmetric on average; which slab wins is random.

Roughness (rough.py; per-column lowest/highest interface rows within Ly/4 of each wall; t >= 1e5):
C(r) = <[h(x+r) - h(x)]^2> ~ 25 r for r <= 50 at every delta and Lx, i.e. zeta_loc(2-8) = 0.50-0.52,
zeta_loc(14-38) = 0.44-0.49; the free fronts of the (1,0) bands at q* (L 4096, t 4e4-2e5): zeta 0.48-0.50
(2-8), 0.46-0.50 (14-38), the same amplitude.  Edge rms sigma at delta 0.011: 28 (Lx 256), 50 (1024),
87 (4096, ~ the layer thickness 100); at Lx 1024, Ly 2048: 35, 41, 42, 50, 55 for delta 0.051 ... 0.0065
and 62 at 0.004 (Ly 4096).  Layer thickness per column: mean 46, 61, 73, 99, 131, 168 with sd 24, 31,
35, 46, 61, 76 (sd/mean ~ 0.5), exponential tails (Fig. 7(b)).  Bare wall (q = 0): sigma 17-18, zeta 0.5.
Collision length of the edges xi_par ~ w^2/25 ~ 400-1200 for w = 100-170: comparable to or larger than
Lx, so the entropic (walks-and-walls, psi = 1/3) regime Lx >> xi_par is not reached; Gamma does not change
when sigma reaches w (Lx 4096).  Applied to Sec. IV D (Gamma list, stationarity to 3e6, sizes, psi 0.5-0.55
to delta 0.004, the 1/3 + offset alternative, the wandering walls, the roughness paragraph), Fig. 6
(fig_wall05.py: Ly series, fits, drift cut) and the new Fig. 7 (fig_rough.py).
