# Campaign: is the symmetric absorbing transition in noisy Day & Night first order?

Model: D&N (B3678/S34678) with α-asynchronous updates plus a symmetric minority flip.
A cell with exactly j = 2 disagreeing neighbours flips with probability q, so f(0) = 0 and
both all-dead and all-alive stay absorbing. Coarsening (ordered) below q*, a disordered
active phase above.

Competing predictions: GV/voter class (ρ ~ 1/ln t at q_c, β = 0); split Ising + DP
transitions (the z=8 case Droz et al. 2003 left unresolved); or first order (the pilot
showed a start-dependent window and a front velocity changing sign near q* ≈ 0.0185).

| exp | what | measures | est core-h* |
|---|---|---|---|
| E0_map, E0_fine | random start, L=1024, T=2e4, α ∈ {0.125..1}, j ∈ {1,2}, wide q | q*(α); done locally 10-14 | 3 |
| E1_front | strip vs all-dead, α=0.5, 11 q in [0.016, 0.026], Lx ∈ {2048, 8192}, Ly=4096, T=1e5 | front velocity v(q): linear sign change ⇒ first order | 40 |
| E2_bulk | random start, α=0.5, 9 q in [0.017, 0.023], L = 512…8192, T=1e5 | steady ρ on the disordered side vs L (jump or → 0), lifetimes, \|m\| | 68 |
| E3_alpha | E1 protocol at α ∈ {0.125, 0.25, 0.75, 1.0}, q ≈ 0.85–1.4 q*(α), T ∝ 1/α | does the jump shrink toward small α (tricritical/GV point)? | 38 |

\*at 7e8 updates/s/core; update with the Torch benchmark.

## Run
    python3 campaign/make_tasks.py --rate <measured> --only E1_front E2_bulk E3_alpha   # root $SCRATCH/dn
    sbatch --array=0-<nshards-1> --export=ALL,EXP=E1_front campaign/campaign.slurm   # etc.
    python3 campaign/summarize.py $SCRATCH/dn/runs/E1_front   # inside the container (numpy)

Runs checkpoint about every 10 min and resume, so resubmitting an array is always safe:
finished runs are skipped and the rest continue from their checkpoints with byte-identical
results. Outputs go to `$SCRATCH/dn/runs/<EXP>/out/<name>.{txt,prof,ckpt,log}` (about 900 files):
/home is at 98% of its 30k file quota (mostly ~/heesch), so nothing per-run goes there.

## Local findings feeding the design (10-14, one seed, L=1024, T=2e4)
q*(α) ≈ 0.04·α: 0.0055 (α=0.125), 0.0085 (0.2), 0.0105 (0.25), 0.014 (0.375), 0.0175 (0.5),
0.024 (0.625), 0.037 (0.75), ~0.08 (1.0). The jump looks smaller at small α, but small α
also means fewer updates per step, hence E3's T ∝ 1/α.

## Phase 1 results (Torch, 10-14; 263 runs, ~2 h wall on 4 nodes; full tables in phase1_results.txt)
Use `fronts.py` for front velocities. It tracks the two strip edges in the profile and fits
only while the interior is still disordered and the strip has not yet filled the box.
summarize.py's `v` column is biased for q < q* (the strip interior orders) and for large q
(the strip fills the box).

At α = 0.5 the transition is first order:
- **Front velocity:** v(q) changes sign linearly at q* = 0.01837 with slope ≈ 4.5–5 rows/step
  per unit q. Lx = 2048 and 8192 agree. Values: −0.0019 (0.018), +0.0007 (0.0185),
  +0.0028 (0.019), +0.0046 (0.0195), +0.0062 (0.020), +0.0089 (0.021), +0.0135 (0.023).
- **Active-phase density:** the active phase stays finite right down to q*. Bulk ρ_bond is
  0.226 at q = 0.0185 and is flat in L from 512 to 8192 and in t from 1e4 to 1e5.
- **Metastability below q*:** the active phase decays by nucleation, with lifetime
  independent of L (512–8192). τ_half = 1.7e4 (q = 0.017), 3e4 (0.0175), 9e4 (0.018).

E3 shows the same signature at every α from 0.125 to 1. q* ≈ 0.0062, 0.0108, 0.0184, 0.038
and 0.084 for α = 0.125, 0.25, 0.5, 0.75 and 1. The reduced slope q*·(dv/dq)/α is about
0.13–0.18 at every α, and the interior interface fraction at q* is 0.66–0.9. So the jump does
not shrink at small α and there is no tricritical point in this range.

Context (closest prior work):
- Rodrigues, Chatelain, Tomé & de Oliveira, J. Stat. Mech. 2015 P01035 (arXiv 1501.01116)
  studied z=6 symmetric-absorbing spin models. They report a discontinuous
  paramagnetic–absorbing line ending at the voter point in mean field and on the 3D cubic
  lattice. Their 2D triangular-lattice runs did not settle the order of that transition.
- Castellano, Muñoz & Pastor-Satorras (PRE 80, 041129, 2009) found a bistable regime in mean
  field and argued that fluctuations remove it in spatially extended systems.
- So 2D first-order coexistence (a front velocity that changes sign, and nucleation lifetimes
  independent of L) in a Life-like CA appears new. The obvious next question is where along
  the path from D&N to the voter model the first-order line ends.

## Phase 2 (submitted 10-14 ~10:00; ~136 core-h)
Local pilots, 10-14 (scratchpad, one seed each):
- **D&N→voter path.** Mix voter steps with probability p into α=0.5 D&N (`--voter p`), with the
  j=2 flip q. q*(p) is non-monotonic: 0.018 (p=0), ~0.036 (0.1), ~0.047 (0.25–0.5), 0.026
  (0.75), 0.0122 (0.9), 0.0065 (0.95), 0.0027 (0.98). So q* ≈ 0.13 (1−p) near the voter
  point. q* rises at first because voter steps flip j=5 cells, which D&N never does, and
  that smooths interfaces.
- **Droplets at α=0.5** (L=1024, T=2e4, 4 seeds). A disordered disk in the dead state survives
  only above a critical radius: R_c ≈ 32–64 (q=0.019), 16–32 (0.020, 0.021), 8–16 (0.023).

| exp | job | runs | what |
|---|---|---|---|
| E4_voter | 19805158 (0-1) | 120 | strip fronts, p ∈ {0.5, 0.75, 0.9, 0.95, 0.98, 0.99}, 10 q per p over 0.7–1.5 q*(p), Lx 2048 × Ly 4096, T 1e5, 2 seeds |
| E4_voter_bulk | 19805159 (0-1) | 180 | random start, same (p, q), L 2048 ×2 seeds and 4096 ×1, T 1e5: does the active density at q* shrink to 0 as p → 1 (end of the first-order line), and where? |
| E5_drop | 19805160 (0) | 1120 | disk R0 on a √2 grid around R_c at q ∈ {0.019, 0.0195, 0.020, 0.021, 0.023}, 32 seeds, T 5e4: P_surv(R0) → R_c(q), check R_c ∝ 1/(q − q*) |

Analysis: `fronts.py runs/E4_voter runs/E4_voter_bulk` and `drops.py runs/E5_drop`.

### Phase 2 results (complete 10-14 16:10; full tables in phase2_results.txt)
All 1420 runs finished (E4_voter 120, E4_voter_bulk 180, E5_drop 1120). E4_voter_bulk timed out at first: voter-mixed runs in the active
phase run at 1.6e8 updates/s, not 6e8. The remaining 122 runs were finished with
`make_tasks.py --remaining --rate 1.6e8` as job 19822043 (0-2, 1.6–2 h per shard).
Phase 2 total ≈ 220 core-h.

**E5_drop** (32 seeds). P_surv(R0) = 1/2 at R_c ≈ 35, 25, 19, 13.5 and 9 for q = 0.019, 0.0195,
0.020, 0.021 and 0.023. R_c·v0(q), with v0 from the E1 flat fronts, is 0.098, 0.115, 0.118,
0.120 and 0.122. This is consistent with a constant D ≈ 0.115 in dR/dt = v0 − D/R: a
curvature-limited critical nucleus, which is the first-order signature. R_c is not simply
∝ 1/(q − q*), because v0(q) is not linear over this range.

**E4_voter fronts** (Lx 2048). v(q) changes sign linearly at every p:

| p | q* | q*·dv/dq |
|---|---|---|
| 0.5 | 0.0459 | 0.21 |
| 0.75 | 0.0276 | 0.22 |
| 0.9 | 0.0119 | ~0.12 |
| 0.95 | 0.0062 | ~0.10 |
| 0.98 | 0.0025 | ~0.05 |
| 0.99 | 0.00127 | ~0.03 |

**E4_voter_bulk** (L 2048 × 2 seeds, L 4096 × 1 seed, T 1e5). The bulk agrees with the
front q* at every p: the active phase decays below q* and holds above it. L = 2048 and
4096 agree everywhere, both the plateaus and τ_half, which is 3e3–7e4 depending on the
distance below q*. Active density just above q*. Drift is the change from
t ∈ [2e4, 4e4] to [8e4, 1e5]:

| p | q/q* | ρ_active | drift | below q* |
|---|---|---|---|---|
| 0 (E2) | 1.007 | 0.226 | flat | τ_half 9e4 at 0.98 q* |
| 0.5 | 1.002 | 0.290 | −0.001 | collapses (ρ 0.03 at 0.95 q*) |
| 0.75 | 1.036 | 0.316 | 0.000 | metastable at 0.989 q* (0.27 → 0.24) |
| 0.9 | 1.025 | 0.291 | 0.000 | τ_half 4.7e4 at 0.974 q* |
| 0.95 | 1.048 | 0.289 | −0.001 | slow decay at 0.995 q* (0.268 → 0.257) |
| 0.98 | 1.028 / 1.08 | 0.258 / 0.276 | −0.006 / −0.002 | 0.245 → 0.219 at 0.972 q* |
| 0.99 | 1.024 / 1.126 | 0.255 / 0.272 | −0.005 / 0.000 | 0.248 → 0.232 at 0.969 q* |

So there is a finite jump at every p, with no sign of ρ_active → 0 at q*+. It is about
0.29 for p = 0.5–0.95 and shrinks to about 0.25–0.26 at p = 0.98–0.99. Within about 3%
above q*, the p ≥ 0.98 plateaus still drift slowly (−0.005 over the run), at roughly
half the rate of pure-voter coarsening. Below q*, single seeds decay gradually rather
than plateau-then-drop, as expected when many nuclei grow at once in L ≥ 2048 (the
lifetimes are L-independent). |m| ≈ 0.01 in the active phase at every p, so it is
paramagnetic.

Pure-voter reference (local, L 1024): at p = 1, q = 0, ρ = 0.30, 0.27, 0.26 at
t = 1e3, 1e4, 3e4. At p = 1 with q = 0.0003, 0.001 and 0.003, ρ_ss = 0.28, 0.30 and 0.33.

**Reading.** The transition is first order along the whole D&N→voter path up to p = 0.99,
with no tricritical point. The jump tends toward the voter model's slowly (logarithmically)
vanishing density as p → 1 and q* → 0. That fits a first-order line that ends at the voter
point, as in the mean-field and 3D results of Rodrigues et al. 2015. The open caveat: a
logarithmically vanishing jump and a short GV segment very close to p = 1 cannot be told
apart at T = 1e5. That needs T ≳ 1e6 bulk runs at p = 0.98–0.995 within a few % of q*.

## Phase 3 (10-14 evening): faster kernel, then long runs near the voter point
**Kernel v1.2.** A single uniform per cell replaces v1.1's one random stream per flip
class. Voter-mixed tables have about 8 classes, so voter-mixed runs are about 4× faster
(1.6e8 → about 6e8 updates/s on Torch), and plain D&N is about 20% faster. The law of the
dynamics is unchanged. Validation is in ../README.md. Phase 1–2 data are v1.1; Phase 3
data are v1.2.

**E6_long** (job 19840430, 3 × 12 cores; about 70 core-h): bulk runs, L 2048, T = 1e6 (10×
Phase 2), 2 seeds. p = 0.98 (q* 0.00252) and 0.99 (q* 0.00127) at q/q* = 0.97, 1.0, 1.01,
1.02, 1.035, 1.05, 1.1, 1.2. Plus p = 1 references: the pure voter model (q = 0), and the
voter model with the j=2 flip at q = 0.0003.
**E6_long_L4096** (job 19840432, 4 cores, 12 h limit; about 31 core-h): one seed at L = 4096,
at 1.02 q* for p = 0.98 and 0.99, at 1.05 q* for p = 0.99, and the pure voter model.
**E6_long_p995** (submitted once q*(0.995) is known from the local strip scan): p = 0.995 on
the same grid, plus L = 4096 at 1.02 q*.

What separates the two scenarios: in a first-order transition, every q > q* settles at a
finite plateau by t = 1e6, and runs just below q* decay. If a GV segment opens near p = 1,
ρ at q* follows the voter model's slow ~1/ln t decay (the p = 1 reference), and the
plateau shrinks continuously toward 0 as q → q*+.

Naming note: E4 run names round q to 5 decimals, while --jflip got the full value (name
q0.00136 = q 0.001365). That shifts the q/q* in the Phase 2 tables by ≤ 0.4%. E6 names
carry 7 digits.

### Phase 3 results: complete (L 2048 10-14 20:05; L 4096 10-15 01:12)
All 52 L = 2048 runs finished: E6_long in 1.4–1.6 h per shard, E6_long_p995 in 1.4 h. The
measured v1.2 rate on cs nodes is 7.7e8 updates/s for voter-mixed L 2048 and 6.1e8 for
L 4096. Tables: `python3 longruns.py $SCRATCH/dn/runs` and `fronts.py` on E6_long*.

q*(0.995) from the local strip scan: v = −0.0016, −0.0001, +0.0007 at q = 0.00058, 0.00061,
0.00064 (2 seeds, Lx = Ly = 2048, T 1e5), so q* ≈ 0.00062. dv/dq ≈ 26, about the same as at
p = 0.98 (20) and 0.99 (24), so q* ≈ 0.125 (1 − p) all the way to the voter point.

Bulk thresholds at T = 1e6 (bracketed by the grid; q/q* relative to the grid's q*):

| p | decays / collapses | stationary plateau (lowest q) | ρ there | drift 1–2e5 → 7–10e5 |
|---|---|---|---|---|
| 0.98 | 0.97 (absorbed ~1e6), 1.00 mixed | 1.02 | 0.261 | −0.004, +0.007 |
| 0.99 | 0.97 (absorbed by 7–9e5), 1.00 | 1.01 | 0.244 | −0.010, +0.006 |
| 0.995 | 0.97, 1.00, 1.01, 1.02 | 1.035 (≈ 1.02 × the bulk q*) | 0.236 | +0.004, −0.006 |
| 1 (voter, q = 0) | — | — | 0.235 → 0.213 | −0.017 (log decay) |

The p = 0.995 bulk threshold is ≈ 0.00063–0.00064, 2–3% above the strip-scan value. The
pure voter at L 4096 follows 1/ρ = 0.184 ln t + 2.13 for t ≥ 1e3 (final, T = 1e6). At L 2048 it already
fluctuates strongly at 1e6 (the two seeds give −0.076 and +0.008). The voter model with the
j=2 flip at q = 0.0003 is stationary at ρ = 0.277.

Reading:
- The active state holds a stationary plateau of ρ ≈ 0.24–0.26 down to within 1–2% of
  the threshold, at every p up to 0.995.
- Just below the threshold the active state is metastable for ~1e5–5e5 steps and then
  collapses to the absorbing state.
- The plateau decreases slowly toward the threshold: at p = 0.99, 1/ρ_ss ≈ 3.60, 3.75,
  3.87, 3.96, 4.02, 4.11 at δ = q/q* − 1 = 0.2, 0.1, 0.05, 0.035, 0.02, 0.01, roughly
  linear in ln(1/δ).
  A finite jump fits this, and so does a jump that vanishes logarithmically (GV-like:
  ρ_ss ≈ the voter density at a crossover time t× ~ δ^−0.8). Getting δ ≥ 0.01 cannot tell
  these apart; this is the intrinsic limit of density measurements this close to a
  logarithmic point.
- The magnetisation correlation length grows toward the threshold. rms m over t ≥ 3e5 is
  0.011–0.03 at 1.2 q* and 0.025–0.08 at ~1.02 q*, i.e. ξ ~ rms m · L ≈ 20–60 → 50–160
  cells. The L 4096 runs test whether ξ stays finite (rms m should halve).
- ρ does not depend on L (L 4096 plateaus 0.259, 0.247, 0.258 vs 0.262, 0.249, 0.259).

Final L 4096 (all four to T = 1e6, one seed each). χ_L = L²⟨m²⟩ and U4 over t ≥ 2e5;
τ_int is the integrated autocorrelation time of m, n_eff = (8e5)/(2 τ_int); the χ error is
χ/√n_eff (a lower bound; √(2/n_eff) for Gaussian m). `python3 chi_fss.py $SCRATCH/dn/runs`.

| p, q/q* | L 2048 χ (2 seeds) | L 4096 χ | ratio | L 4096 U4 | τ_int 2048 → 4096 | n_eff 4096 |
|---|---|---|---|---|---|---|
| 0.98, 1.02 | 3.1k, 2.9k | 2.5k ± 0.6k | 0.8 | 0.05 | 4–5e4 → 2.5e4 | 16 |
| 0.99, 1.02 | 6.2k, 5.6k | 21.7k ± 12k | 3.7 ± 2 | 0.48 | 4.5–6.4e4 → 1.3e5 | 3 |
| 0.99, 1.05 | 4.6k, 3.3k | 6.4k ± 1.8k | 1.6 ± 0.5 | 0.22 | 3–5e4 → 3.1e4 | 13 |
| 0.995, 1.02 | 23k, 24k | 34k ± 20k | 1.4 | 0.39 | 6–8e4 → 1.3e5 | 3 (decaying) |

- p = 0.98: χ flat, τ flat. Finite ξ ≈ 50 cells; first order.
- p = 0.99, 1.02 q*: the ×3.7 jump survived to T = 1e6, but it is one slow swing. m goes
  from −0.046 (t 2–3e5) to +0.042 (7e5) with ρ steady at 0.245–0.252, so n_eff ≈ 3 and the
  ratio is 1.3σ from 1. τ_int also doubled. Consistent with Ising-like growth (×3.4) and
  with one large fluctuation; it needs more seeds (E7 tier 1).
- p = 0.99, 1.05 q*: ×1.6 ± 0.5, between paramagnet (×1) and Ising (×3.4).
- p = 0.995, 1.02 q*: below the bulk threshold at L 4096 (ρ 0.239 → 0.179, still falling
  at 1e6). Its m drifts +0.05 → −0.075 as it coarsens, so χ there is not an equilibrium value.
- Pure voter at L 4096: m wanders ±0.035 over 1e6 steps, the same size as the p = 0.99,
  1.02 q* swings. Near q* the magnetisation is nearly voter-conserved, so τ_m is long and one
  run cannot measure χ.

## Phase 4 (designed 10-15): Binder cumulant and spreading. E7 tier 1 done 10-15 15:21; E8b done 10-15 19:10; E8a done 10-15 21:04 (verdict: first order up to p = 0.995)
Why: a referee check (10-15) turned up Azizi & Pleimling, PRE 102, 022112 (2020), arXiv 2007.10933.
In their 2D model with two symmetric absorbing states, the Ising and DP lines merge into a GV line
that ends away from the voter point. The Phase 1–3 density data near p = 1 cannot separate three
scenarios: (i) first order up to the voter point, (ii) a short GV segment, (iii) a split into an
Ising transition plus a DP transition, with a magnetised active phase in between. Two experiments
separate them: magnetisation finite-size scaling (E7) for (iii), and time-dependent spreading
(E8) for (ii).

### Calibration (local, 10-15, scratchpad p4cal/, v1.2, T 4e5)
χ_L = L²⟨m²⟩ over the active part of each run (t ≥ 2e4, ρ > 0.15). L 512 is local; L 2048 is
E6 (2 seeds); L 4096 is E6, one seed, partial (t ≈ 7.5e5). The final L 4096 values (T 1e6,
t ≥ 2e5) are in "Phase 3 results": 2.5k, 21.7k, 6.4k for the first three rows.

| p, q/q* | L 512 | L 2048 | L 4096 | reading |
|---|---|---|---|---|
| 0.98, 1.02 | 1.9k | 2.8k, 2.7k | 2.3k | flat: finite ξ |
| 0.99, 1.02 | 4.7k, 7.6k | 6.2k, 5.2k | 19.5k | flat to 2048, then ×3.4 (one run, n_eff ≈ 3) |
| 0.99, 1.05 | — | 5.2k, 2.8k | 8.6k | ×2 |
| 0.995, 1.05 | 3.8k | 6.9k, 4.1k | — | |

For comparison, an Ising critical point gives ×3.4 per doubling of L (L^1.75), a paramagnet ×1.
- U4 ≈ 0 (−0.08 to +0.18) at L 512, so m is Gaussian there.
- τ_m is 2–8×10⁴ steps at L 512, the same as at L 2048, so it does not grow with L.
- **L 256 is too small:** all five L 256 runs reached an absorbing state within 3×10⁴–2.6×10⁵
  steps. L 512 runs survived 4e5.

### E7_binder: magnetisation FSS (no kernel change; ~250 core-h)
Bulk random start, `--every 500`, T = 1e6 (far points 5e5), discard t < 2e5. Observables per
(p, q, L), pooled over seeds and restricted to the active window (ρ above 0.8 × its plateau):
- χ_L = L²⟨m²⟩;
- U4 = 1 − ⟨m⁴⟩/(3⟨m²⟩²);
- the histograms P(m) and P(m, ρ);
- τ_m;
- block errors with blocks ≥ 4 τ_m.

Target: n_eff ≥ 25 independent samples per point at L 4096 (±28% on χ) and ≥ 50 at L ≤ 1024.
U4's statistical error is about 1.6/√n_eff, so U4 is a secondary check and χ_L is primary.

| points | L 512 | L 1024 | L 2048 | L 4096 | core-h |
|---|---|---|---|---|---|
| near: (0.98, 1.02) (0.99, 1.01) (0.99, 1.02) (0.99, 1.05) (0.995, 1.05) | 8 seeds | 8 | 4 new (+2 E6) | 4 new (+1 E6 where it exists; else 5) | ~41 each, ~215 |
| far: (0.98, 1.2) (0.99, 1.2) (0.995, 1.2); control (0.9, 1.03) | 4 | 4 | 2 new | 2 | ~10 each, ~40 |

Rates assumed: 8e8 updates/s at L ≤ 1024, 7.7e8 at L 2048, 6.1e8 at L 4096. Wall time is
7.6 h per L 4096 run single-threaded, so give L 4096 tasks 2 threads. The L 4096 seeds are 85%
of the cost. The cheaper tier-1 version runs only (0.99, 1.02) and (0.99, 1.05): about 80 core-h,
and it answers the question raised by the one partial L 4096 run.

**Tier 1 submitted 10-15 12:12** (Matias approved "E7"; read as the tier 1 that was proposed):
- `E7_binder`, job 19950228 (cs602, 16 cores, 3 h limit): 40 runs at L 512 / 1024 (s1–s8) and
  L 2048 (s3–s6) at p = 0.99, q = 0.0012954 (1.02 q*) and 0.0013335 (1.05 q*). Est. 19.7 core-h,
  ~1.5 h wall.
- `E7_binder_L4096`, job 19950229 (cs610, 16 cores, 6 h limit): L 4096 s2–s5 at the same two q,
  8 runs × 2 threads (`THREADS=2`, new in run_shard.sh). Est. ~3.6 h wall, ~58 core-h.
- Measured: 1.33e9 updates/s per 2-thread L 4096 run, 7.9e8 per L 2048 run. Output is
  bit-identical for 1 and 2 threads (checked locally, L 1024 T 3000).
- Resubmit after a timeout with the same command; runs resume from their checkpoints:
  `sbatch -J dn-e7 --array=0 -t 06:00:00 --export=ALL,EXP=E7_binder_L4096,THREADS=2 campaign.slurm`
- Analysis: `python3 binder.py $SCRATCH/dn/runs` pools E6 and E7 by (p, q, L): χ with
  seed-scatter errors, U4 with jackknife errors, ⟨|m|⟩, τ_m, n_eff, and χ_L/χ_{L/2}.
  With E6, tier 1 gives 8 / 8 / 6 / 5 runs at L 512 / 1024 / 2048 / 4096 per point.

### E7 tier 1 results (10-15; both jobs COMPLETED rc 0: 19950228 in 1.5 h, 19950229 in 3.15 h)
All 48 runs reached T = 1e6, except 4 of the 8 L 512 runs at 1.02 q*, which were absorbed at
t = 2.4e5–8.6e5. L 4096 ran at 1.48e9 updates/s per 2-thread run. `python3 binder.py`, pooled
with E6 (t ≥ 2e5, active window):

| q/q* | L | runs | χ_L ± scatter | χ_L/χ_{L/2} | U4 ± jackknife | ⟨\|m\|⟩ | per-run χ |
|---|---|---|---|---|---|---|---|
| 1.02 | 512 | 7 of 8 (4 absorbed) | 3.1k ± 0.5k | — | 0.17 ± 0.06 | 0.089 | 1.5k–4.9k |
| 1.02 | 1024 | 8 | 5.3k ± 0.6k | 1.69 | 0.17 ± 0.08 | 0.059 | 3.3k–7.9k |
| 1.02 | 2048 | 6 | 8.3k ± 1.5k | 1.58 | −0.15 ± 0.25 | 0.035 | 4.1k–14.6k |
| 1.02 | 4096 | 5 | 6.4k ± 3.6k | 0.76 | −0.16 ± 0.20 | 0.015 | 1.3k, 3.6k, 6.2k, 6.7k, 21.7k (E6) |
| 1.05 | 512 | 8 | 2.7k ± 0.2k | — | 0.18 ± 0.03 | 0.083 | 2.2k–4.0k |
| 1.05 | 1024 | 8 | 3.8k ± 0.4k | 1.42 | 0.00 ± 0.11 | 0.048 | 1.6k–5.7k |
| 1.05 | 2048 | 6 | 4.4k ± 0.4k | 1.16 | −0.29 ± 0.11 | 0.025 | 3.4k–5.7k |
| 1.05 | 4096 | 5 | 3.5k ± 0.8k | 0.81 | 0.02 ± 0.07 | 0.012 | 2.0k–6.4k |

Reading:
- **No Ising ordering at p = 0.99 for q ≥ 1.02 q\*.** χ_L saturates by L ≈ 2048 at both points:
  ξ = √χ ≈ 65–90 cells, and ⟨|m|⟩ ∝ 1/L for L ≥ 1024. Ising growth from L 2048 would give
  ≈ 28k at L 4096; 4 of the 5 runs there are ≤ 6.7k. U4 stays in [−0.3, 0.2], with no drift
  toward 0.61 or 2/3. The E6 ×3.7 was one run's slow m swing.
- Earlier growth (512 → 2048) is the approach to saturation, since L 512 ≈ 6–8 ξ. At 1.02 q*,
  L 512 is also absorbing-prone. Negative U4 at L 2048 comes from slow m wandering (~7
  independent samples per run), not a trend.
- **What is left of scenario (iii):** an Ising line would have to sit between the bulk
  threshold (≈ 1.00–1.01 q*) and 1.02 q*.
- **Still open, scenario (ii) GV:** ξ grows toward q*. At L 2048 (E6 + E7), χ ≈ 0.84k, 2.1k,
  4.4k, 8.3k at δ = 0.2, 0.1, 0.05, 0.02, roughly χ ∝ δ^−1 (ξ ∝ δ^−1/2). δ = 0.01 gives 2.8k,
  but from 2 seeds only. A finite ξ at q* (first order) versus ξ → ∞ (GV) is not settled by
  E7 tier 1. That is what E8 (spreading, domain wall) and the δ = 0.01 point of the full E7
  grid would test.

### E8_spread: time-dependent spreading (kernel v1.3 needed for E8a; ~100–150 core-h)
**E8a, seed spreading.**
- Start: an all-dead lattice (L 2048 torus) with a 2×2 live block at the centre (`--init disk:1:0.75`).
- Per replica, at log-spaced t up to T = 1e5:
  - survival s(t) = [any interface cell];
  - interface cells I(t);
  - live cells A(t);
  - R²(t), the mean squared distance of interface cells from the seed.
- Aggregate:
  - P_s(t) = ⟨s⟩ ~ t^−δ;
  - N(t) = ⟨I⟩ over all replicas ~ t^η;
  - R²(t) over survivors ~ t^z;
  - local slopes over factors of 2–4 in t.
- Grid: p ∈ {0.9, 0.98, 0.99, 0.995} × q/q* ∈ {0.95, 1.0, 1.02, 1.05, 1.1, 1.2}. Validation and
  reference points: p = 1, q = 0 (pure voter, the GV point) and p = 1, q = 3e-4.
- Replicas are adaptive: run until 300 replicas survive to T, or 3e6 replicas.
- T = 1e5 is the practical ceiling. Voter-like survival is P_s ~ ln t / t, so going to T = 1e6
  would cost ~100× more. Near p = 1 the crossover away from voter behaviour happens at
  10³–2×10⁴ steps (E6), so T = 1e5 gives about one decade after it.

Kernel v1.3 (about 200 lines):
- `--box`: step only the bounding box of interface cells plus one cell. This is exact: cells
  outside have j = 0 and f(0) = 0, so they never flip. The RNG counter is (t, y, x), so the
  result is bit-identical to a full step. The box grows when activity reaches its edge and is
  re-scanned to shrink every 32 steps. If it comes within 2 cells of the torus seam, the replica
  is flagged.
- `--reps N`: run N replicas in one process. Replica i uses Philox key (seed, i) and clears only
  its own box. The accumulators (Σs, ΣI, ΣI², ΣA, ΣI·R², Σ box area) are written at the output
  times, plus a list of per-replica death times. Checkpoint every 1000 replicas.
- Selftest: box stepping vs full stepping, bit-identical over 200 voter-mixed seeds.
- Validation run: pure voter, where ⟨A(t)⟩ = 4 exactly (a martingale) and t·P_s/ln t → const
  (Bramson & Griffeath).

Cost: a dying replica costs ~10⁶–10⁷ updates; a survivor to 1e5 costs ~2×10¹⁰. That is ≤ 8
core-h per point and 26 points, ~100–150 core-h, since most points hit the survivor target or die
out early. Full-lattice stepping at L 2048 would cost 20–60× more (several thousand core-h).

**E8b, domain wall** (works with v1.2; `--box` would make it ~4× cheaper).
- Start: `--init strip:1:0.5 --bg 0`, i.e. two flat walls between all-alive and all-dead.
  Lx 1024 × Ly 4096, T 1e6, 2 seeds.
- Grid: p ∈ {0.9, 0.98, 0.99, 0.995} × q/q* ∈ {0.98, 1.0, 1.02, 1.05}. That is 32 runs, ~48 core-h
  at v1.2.
- From `--profile`: the active-layer width W(t) and its interior density ρ_in(t).
- No active phase exists at t = 0, unlike the E4 strips. So this asks whether an active layer
  forms at all, and whether it broadens ballistically or diffusively.

### E8 submitted (10-15; Matias: "go for E8")
**E8b_wall**, job 19985325 (cs631, 32 cores, 4 h limit), v1.2 binary:
- 32 runs, as designed: p ∈ {0.9, 0.98, 0.99, 0.995} × q/q* ∈ {0.98, 1.0, 1.02, 1.05} × 2 seeds.
- Lx 1024 × Ly 4096, T 1e6, `--every 10000 --ppd 20 --profile`. q* is the front value
  (`QSTAR_FRONT`: 0.0119, 0.00252, 0.00127, 0.00062).
- Est. 31 core-h; the absorbing bulk runs at ~4e9 updates/s early on.
- Analysis: `python3 wall.py $SCRATCH/dn/runs` (summary) or `wall.py --run PREFIX`.

**E8a_spread**, job 19986484 (26 cores, 12 h limit), kernel v1.3 (`bin/dnsim_v13`, `--reps`):
- 26 points: the same 4 p × q/q* ∈ {0.95, 1.0, 1.02, 1.05, 1.1, 1.2}, plus p = 1 at q = 0
  and q = 3e-4.
- L 2048, T 1e5, `--init disk:1:0.75`, `--ppd 20`.
- Each run stops at 300 replicas surviving to T (100 for x ≥ 1.1 and for voter + noise) or at
  3e6 replicas.
- Est. ~40–80 core-h, packed as 3 core-h each.
- Analysis: `python3 spread.py --summary $SCRATCH/dn/runs` or `spread.py FILE`.
- A run killed at the time limit resumes from its checkpoint, with identical results. To extend
  a finished run, resubmit it with a larger `--reps` or `--target`.

Kernel v1.3 validation (local + Torch build): see dn/README.md.
- Selftest: box stepping is bit-identical to full stepping over 200 seeds.
- Normal mode is byte-identical to v1.2.
- The Torch and local builds give the same md5 in both modes.
- Pure voter, 200,000 replicas: ⟨A(t)⟩ = 4 within |z| < 0.94 at all 38 times; z_eff 0.98–0.99;
  δ_eff → 0.94 at t ~ 3e3.

E8a calibration (local, L 2048, T 1e5; scratchpad e8cal):
- Near threshold, a 2×2 seed follows the voter survival curve to t ~ 10³ and then dies:
  - p 0.98, 1.02 q*: 0 of 1000 alive at T;
  - p 0.995, 1.05 q*: 0 of 1000;
  - p 0.99, q*: 1 of 2000.
- p 0.9 at q*: δ_eff rises past 1 (1.2–1.4), i.e. faster than voter decay, as a first-order
  transition needs (a seed must reach R_c).
- P_s(1e5) ≈ 0.7–1.7% at 1.2 q*.
- Cost: 0.1–0.4 ms per dying replica, 20–50 s per survivor; one p 0.9, 1.2 q* survivor
  escaped at t ≈ 6e4.

### E8b results (job 19985325 COMPLETED rc 0 at 19:10, 1.5 h wall; all 32 runs to T = 1e6)
`python3 wall.py RUNS_DIR`. Two new columns, computed per wall (2 walls per run):
- **W**: the distance between the outermost rows whose interface density (9-row running mean)
  exceeds 0.3. This is the slab width. The window is 1024 rows either side of each wall, so W
  saturates at 2047.
- **rc**: the interface density over those W rows, gaps included.

Values are means over walls and seeds. v = (dW/dt)/2 is fitted over t 1e5–1e6 while
W < 1640.

| p | q/q* | W 1e4 | W 1e5 | W 3e5 | W 1e6 | W(1e6)/W(1e5) | rc 1e5 | rc 1e6 | v (rows/step) |
|---|---|---|---|---|---|---|---|---|---|
| 0.9 | 0.98 | 127 | 167 | 175 | 200 | 1.2 | 0.45 | 0.48 | 0.0 |
| 0.9 | 1.00 | 147 | 246 | 308 | 264 | 1.1 | 0.52 | 0.54 | 0.0 |
| 0.9 | 1.02 | 168 | 581 | 1131 | 2047 | (full) | 0.60 | 0.71 | 1.4e-3 |
| 0.9 | 1.05 | 204 | 1371 | 2047 | 2047 | (full) | 0.69 | 0.74 | 4.8e-3 |
| 0.98 | 0.98 | 189 | 375 | 366 | 450 | 1.2 | 0.53 | 0.48 | 1.0e-4 |
| 0.98 | 1.00 | 202 | 481 | 718 | 1206 | 2.5 | 0.54 | 0.49 | 3.0e-4 |
| 0.98 | 1.02 | 196 | 670 | 1143 | 2047 | (full) | 0.53 | 0.66 | 1.4e-3 |
| 0.98 | 1.05 | 192 | 854 | 2040 | 2047 | (full) | 0.58 | 0.68 | 3.0e-3 |
| 0.99 | 0.98 | 203 | 502 | 784 | 672 | 1.3 | 0.50 | 0.51 | 0.3e-4 ± 0.7e-4 |
| 0.99 | 1.00 | 181 | 486 | 618 | 1153 | 2.4 | 0.52 | 0.50 | 3.0e-4 |
| 0.99 | 1.02 | 200 | 670 | 1150 | 2046 | (full) | 0.55 | 0.59 | 1.0e-3 |
| 0.99 | 1.05 | 194 | 626 | 1485 | 2047 | (full) | 0.54 | 0.65 | 2.0e-3 |
| 0.995 | 0.98 | 189 | 537 | 823 | 820 | 1.5 | 0.48 | 0.46 | 1.9e-4 |
| 0.995 | 1.00 | 191 | 574 | 672 | 1157 | 2.0 | 0.52 | 0.45 | 2.4e-4 |
| 0.995 | 1.02 | 174 | 513 | 795 | 1464 | 2.9 | 0.48 | 0.48 | 5.6e-4 |
| 0.995 | 1.05 | 218 | 537 | 1128 | 1994 | 3.7 | 0.47 | 0.53 | 1.0e-3 |

References: local runs, v1.3, one seed each, same geometry (scratchpad e8bref/):
- pure voter wall (p = 1, q = 0): W = 69, 192, 699, 1324 at t = 1e3, 1e4, 1e5, 3e5, then
  1400–2000 to 1e6 (it reaches the 2047 window by 4e5). rc = 0.60, 0.53, 0.50, then 0.37–0.54
  (noisy, no clear trend);
- p = 0.99, q = 0 (deep in the ordered phase): W = 102–116 from t = 1e4 to 1e6, a thin rough wall.

Reading:
- **Voter-like start.** For t ≲ 10⁴, every point with p ≥ 0.98 gives the same layer as the pure
  voter wall: W ∝ t^½, W ≈ 190 at t = 1e4, independent of q.
- **Then q takes over, at every p.** Let t_x be the time from which W(1.05 q*)/W(0.98 q*)
  stays above 1.5: ≈ 8e3 (p 0.9), 4e4 (0.98), 3e5 (0.99), 4e5 (0.995). That is roughly
  ∝ (1 − p)^−1.3, but with 2 seeds the individual values are rough.
- **Below q\* the wall becomes a stationary active slab.** At 0.98 q*, W grows only ×1.2–1.5
  from 1e5 to 1e6. From 3e5 on it stops: p 0.99 goes 784 → 672 and p 0.995 goes 823 → 820,
  while the pure voter wall keeps widening (699 at 1e5, 1324 at 3e5, 1400–2000 at 1e6). The
  slab is still 4–8 times wider than the q = 0 wall (≈ 100). Its interior keeps the
  active-phase density: rc = 0.46–0.51, against 0.47–0.58 for the bulk active phase at q* (E6). So the
  boundary between the two absorbing phases is wetted by a layer of the active phase. The
  layer is thicker the closer p is to 1 (W(1e6) = 200, 450, 670, 820 rows for p = 0.9, 0.98,
  0.99, 0.995). An active strip shrinks below q* (E4, v < 0), but here the slab cannot vanish,
  because the two absorbing phases must stay separated.
- **Above q\* the slab edges are fronts, with the E4 speeds.** dv/d(q/q*) from 1.02 to
  1.05 q* is 0.11, 0.053, 0.033, 0.015 for p = 0.9, 0.98, 0.99, 0.995. The E4 front slopes
  q*·dv/dq are ~0.12, ~0.05, ~0.03 (0.9, 0.98, 0.99). So the slab grown from a wall
  between two absorbing phases is the same active phase as a prepared active strip. At p ≥ 0.98
  the fronts are still slow: v·t ≲ 10³ rows by 1e6, comparable to the √t spreading. So
  W(1e6)/W(1e5) is still ≈ 3 at 1.02–1.05 q*, and the purely ballistic regime lies beyond
  this lattice and time.
- **At q\* itself (p ≥ 0.98)** the slab still widens slowly: ×2.0–2.5 per decade, v ≈ 2.4–3e-4.
  Either the wall threshold sits ~1% below the front q*, or the width creeps by fluctuations
  at coexistence. Both fit first order. A GV point would need W ∝ t^½ to continue indefinitely
  at some q_c. Between 1.0 and 1.02 q*, W grows at about that rate over this window (as fast as
  the voter wall at 1.02 q*). So E8b alone cannot exclude a GV point in that 2% gap, even though
  the bulk there is active with finite ξ (E7). E8a addresses that gap directly.
- **What E8b does not settle.** rc does not separate the scenarios here. The voter wall's
  interior shows no clear thinning over 1e4–1e6 either (0.53, then 0.37–0.54 in one noisy run),
  because a 1/ln t decay is too slow. The evidence is W. The first-order behaviour appears
  only after t_x, and t_x grows toward p = 1. Extrapolating, p = 0.999 would need ≳ 5e6 steps, so E8b says nothing about p > 0.995.

### E8a results (job 19986484 COMPLETED rc 0 at 21:04, 3.1 h wall on 26 cores; all 26 points to T = 1e5)
`python3 spread.py --brief $SCRATCH/dn/runs` (one line per point) or `--summary` (full tables).
Each point ran until 300 replicas survived to T (100 for x ≥ 1.1 and for voter + noise) or until
the 3e6-replica cap. δ_eff = −dln P_s/dln t over a factor 4 in t, so 3e4 is the last time the
window fits; ± is the jackknife over groups. "÷ voter" is P_s(1e5) over the pure-voter value.

| p | q/q* | replicas | alive at 1e5 | P_s(1e5) | ÷ voter | δ_eff 3e3 | 1e4 | 3e4 | η_eff 3e4 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | q = 0 | 1.91M | 300 | 1.57e-04 | 1 | 0.93 | 0.93 | 0.93 ± 0.02 | −0.08 |
| 1 | q = 0.0003 | 0.05M | 100 | 2.12e-03 | 13.5 | 0.66 | 0.44 | 0.11 ± 0.02 | 1.73 |
| 0.9 | 0.95 | 3.00M | 0 | 0.00e+00 | 0 | 2.78 | — | — | — |
| 0.9 | 1.00 | 3.00M | 0 | 0.00e+00 | 0 | 1.92 | 3.17 | — | — |
| 0.9 | 1.02 | 3.00M | 3 | 1.00e-06 | 0.00635 | 1.70 | 1.95 | 1.64 ± 0.28 | −1.16 |
| 0.9 | 1.05 | 3.00M | 226 | 7.53e-05 | 0.479 | 1.35 | 1.10 | 0.41 ± 0.03 | 1.34 |
| 0.9 | 1.10 | 0.15M | 100 | 6.82e-04 | 4.33 | 0.86 | 0.30 | 0.04 ± 0.01 | 2.12 |
| 0.9 | 1.20 | 0.02M | 100 | 4.08e-03 | 25.9 | 0.34 | 0.03 | −0.00 ± 0.00 | — |
| 0.98 | 0.95 | 3.00M | 0 | 0.00e+00 | 0 | 1.31 | 1.74 | 4.04 | −3.78 |
| 0.98 | 1.00 | 3.00M | 10 | 3.33e-06 | 0.0212 | 1.22 | 1.43 | 2.19 ± 0.11 | −1.44 |
| 0.98 | 1.02 | 3.00M | 54 | 1.80e-05 | 0.114 | 1.16 | 1.30 | 1.54 ± 0.07 | −0.73 |
| 0.98 | 1.05 | 3.00M | 291 | 9.70e-05 | 0.616 | 1.07 | 1.08 | 0.96 ± 0.03 | 0.19 |
| 0.98 | 1.10 | 0.25M | 100 | 4.04e-04 | 2.57 | 0.99 | 0.80 | 0.41 ± 0.04 | 1.10 |
| 0.98 | 1.20 | 0.06M | 100 | 1.61e-03 | 10.2 | 0.79 | 0.42 | 0.06 ± 0.02 | 1.89 |
| 0.99 | 0.95 | 3.00M | 23 | 7.67e-06 | 0.0487 | 1.13 | 1.35 | 1.92 ± 0.08 | −1.23 |
| 0.99 | 1.00 | 3.00M | 103 | 3.43e-05 | 0.218 | 1.07 | 1.19 | 1.38 ± 0.03 | −0.59 |
| 0.99 | 1.02 | 3.00M | 173 | 5.77e-05 | 0.366 | 1.04 | 1.12 | 1.28 ± 0.03 | −0.41 |
| 0.99 | 1.05 | 2.53M | 300 | 1.19e-04 | 0.754 | 1.02 | 1.02 | 0.97 ± 0.02 | −0.04 |
| 0.99 | 1.10 | 0.36M | 100 | 2.76e-04 | 1.75 | 0.97 | 0.93 | 0.64 ± 0.07 | 0.59 |
| 0.99 | 1.20 | 0.11M | 100 | 9.39e-04 | 5.96 | 0.79 | 0.65 | 0.35 ± 0.05 | 1.41 |
| 0.995 | 0.95 | 3.00M | 69 | 2.30e-05 | 0.146 | 1.04 | 1.18 | 1.54 ± 0.04 | −0.93 |
| 0.995 | 1.00 | 3.00M | 144 | 4.80e-05 | 0.305 | 1.01 | 1.12 | 1.34 ± 0.03 | −0.59 |
| 0.995 | 1.02 | 3.00M | 168 | 5.60e-05 | 0.356 | 1.02 | 1.09 | 1.26 ± 0.04 | −0.49 |
| 0.995 | 1.05 | 3.00M | 262 | 8.73e-05 | 0.555 | 0.97 | 1.03 | 1.13 ± 0.03 | −0.28 |
| 0.995 | 1.10 | 0.58M | 100 | 1.74e-04 | 1.1 | 0.97 | 0.95 | 0.92 ± 0.06 | 0.07 |
| 0.995 | 1.20 | 0.26M | 100 | 3.86e-04 | 2.45 | 0.90 | 0.83 | 0.63 ± 0.05 | 0.56 |

Reading:
- **Validation.** Pure voter: ⟨A⟩ = 4.02 at 1e4 (a martingale, exactly 4), δ_eff = 0.93 flat
  over 3e3–3e4, and t·P_s/ln t = 1.37 at 1e5, as P_s ~ ln t / t requires. Voter plus the j=2
  flip (q = 3e-4) is active: P_s is flat at 2.1e-3 from 3e4 on.
- **Seeds die faster than voter seeds where the bulk active phase is stable.** Phase 3 found
  stationary bulk plateaus from 1.02, 1.01 and 1.035 q* at p = 0.98, 0.99, 0.995. At those q,
  and up to a few % above, δ_eff is above 1 and still rising at the end of the window:
  p 0.98, 1.02 q*: 1.16 → 1.30 → 1.54 ± 0.07; p 0.99, 1.02 q*: 1.04 → 1.12 → 1.28 ± 0.03;
  p 0.995, 1.05 q*: 0.97 → 1.03 → 1.13 ± 0.03 (t = 3e3, 1e4, 3e4). η_eff is negative there
  (−0.73, −0.41, −0.28), so the interface population is shrinking too.
  At a GV point δ_eff stays at its critical value (→ 1 from below, like the voter curve), at a
  DP point at 0.45, and above either it falls. A local exponent above 1 and rising means the seed
  is subcritical: it has to grow past a critical droplet before it can spread. That is a
  nucleation barrier, i.e. first order. This covers the 1.00–1.02 q* gap that E8b left open at
  p 0.98–0.99, and at p 0.995 it extends to 1.05 q*.
- **The separatrix moves toward q\* with time.** x where δ_eff crosses the voter value:

  | p | t = 3e3 | 1e4 | 3e4 | P_s(1e5) = voter |
  |---|---|---|---|---|
  | 0.9 | 1.093 | 1.061 | 1.037 | 1.057 |
  | 0.98 | 1.132 | 1.076 | 1.052 | 1.060 |
  | 0.99 | 1.121 | 1.101 | 1.055 | 1.062 |
  | 0.995 | 1.158 | 1.115 | 1.096 | 1.091 |

  Just above it the seed overshoots and then turns toward a plateau: at 1.05 q* δ_eff peaks at
  1.11 (p 0.98, t ≈ 1.3e4) and 1.04 (p 0.99, t ≈ 6e3), then falls to 0.96–0.97 by 3e4. First
  order predicts this drift: above threshold the critical droplet is finite, so every seed that
  reaches it survives, and the time to do so (∝ R_c²) grows as q → q*. The crossing therefore
  tends to the threshold as t → ∞.
- **Above the separatrix, P_s plateaus.** At 1.2 q*: 4.1e-3 (p 0.9) and 1.6e-3 (p 0.98), with
  η_eff → 1.9–2.2 (compact growth). At p 0.99 and 0.995, P_s at 1.2 q* is still falling at
  1e5 (δ_eff 0.35, 0.63).
- **Below and at q\*** seeds die outright at p 0.9 (0 of 3e6 at 1e5 for x ≤ 1.0) and at p 0.98,
  0.95 q*. Elsewhere δ_eff climbs to 1.3–4.
- **The barrier weakens toward p = 1.** P_s(1e5)/voter at 1.02 q* is 0.006, 0.11, 0.37, 0.36 for
  p = 0.9, 0.98, 0.99, 0.995. The same trend as E8b's t_x.

Caveats:
- The reference is the p = 1 voter curve. A GV critical curve at p < 1 would have the same
  asymptotic form but possibly a different amplitude and corrections. That is why the argument
  rests on δ_eff exceeding 1 and rising, not on P_s/voter < 1.
- The argument needs the Phase 3 bulk thresholds (T = 1e6, a 1–2% grid). At p 0.98–0.99 the
  decisive points (1.02 q*) sit at the lowest stationary plateau. The cleanest point is p 0.995,
  1.05 q*, 1.5% above its bulk threshold.
- P_∞ > 0 just above the bulk threshold is predicted (the seeds should eventually plateau) but not
  shown: at 1.02 q* and p ≥ 0.98 the window ends while δ_eff is still rising.
- Nothing here about p > 0.995.

### E8 verdict (10-15)
Walls (E8b) and seeds (E8a) both point to first order up to p = 0.995. Below q* a wall becomes a
stationary active slab, and above q* its edges move at the E4 front speeds. Just above the bulk
threshold, single seeds die faster than voter seeds with an accelerating decay, which is a
nucleation barrier. A GV point would need a q_c where the bulk is active and seeds decay no
faster than voter seeds, and none of the 16 points between q* and 1.2 q* at p 0.98–0.995 shows
that. Together with E7 tier 1 (no Ising ordering at p 0.99), scenarios (ii) and (iii) are
disfavoured up to p = 0.995. The voter-like transient grows toward p = 1 (t_x ∝ (1 − p)^−1.3 in
E8b), so p > 0.995 stays open.

### What each scenario predicts
| | E7: χ_L, U4 | E8a: seed | E8b: wall |
|---|---|---|---|
| (i) first order to p → 1 | χ_L flat for L ≥ 1024 at every δ ≥ 0.01; U4 → 0; ξ(δ) stays finite as δ → 0 | no power-law window at any q: P_s plunges (δ_eff grows) for q ≤ q*; P_∞ > 0 above q* with an essential-singularity onset (a seed must reach R_c ∝ 1/(q − q*)) | W ≈ W₀ + 2v(q)t with v linear through 0; W bounded at q*; ρ_in flat at the plateau |
| (ii) GV segment | no ordered-active window; χ_∞(δ) diverges as δ → 0 | δ_eff → 1 (P_s ~ ln t / t), z → 1, ⟨A⟩ ≈ const at q_c | W ~ t^1/2 and ρ_in ~ 1/ln t at q_c |
| (iii) split Ising + DP | χ_L ~ L^1.75 at q_I above the bulk threshold, U4 crossing ≈ 0.61, then U4 → 2/3 with ρ > 0 (magnetised active phase) | DP exponents at q_DP: δ = 0.451, η = 0.230, z = 1.13 | layer grows ballistically above q_DP (not decisive) |

Analysis scripts to write: binder.py (χ, U4, P(m, ρ), τ_m, block errors), spread.py (P_s, N, R²,
local slopes, P_∞(q)), wall.py (W(t), ρ_in(t) from the profiles).

Suggested order:
1. E7 tier 1 (~80 core-h, no code). It settles the p = 0.99 χ jump.
2. Kernel v1.3 with local validation.
3. E8a.
4. The rest of E7.
5. E8b.

## Theory: pair approximation (10-16; `theory/`, details in theory/README.md)
PA and site MF are bistable at every q ≥ 0 and p < 1. The absorbing states are always linearly
stable, and the disordered active state is also stable, even for pure D&N at q = 0. So the
transition is discontinuous at the mean-field level, but homogeneous PA has no q*: only fronts
select the phase. A slow-manifold front criterion gets q* ∝ (1−p) near the voter point, but the
synchronous map is low by a factor of 2.0 (continuous time 5.5) for p ≥ 0.9, and it fails
qualitatively at p ≲ 0.5. The PA ρ* ≈ 0.43–0.50 is about 2× the measured ρ_bond. Next options:
perturbation theory around the voter point, a 2×2 plaquette closure.

## Theory: voter-point perturbation theory (10-16; theory/README.md)
Near p = 1, f = α·f_R + κ·f₂ exactly (κ = q/(1−p)), computed from coalescing Moore walks
(`theory/cwalk.c`, `vpt.py`). Prediction: q* = c·α·(1−p) with c = 0.26 ± 0.015, no fit
parameters. It matches E4/E6 (c = 0.238–0.254 at α = 0.5). E9_alpha (Phase 5, 4 seeds,
p 0.98 and 0.99, α 0.25–1) confirms q* ∝ α (exponent 1.00 ± 0.02) with c = 0.248–0.258. The coexistence
window shrinks like powers of 1/ln(1/ε), so the first-order line weakens logarithmically and
ends at the voter point. Literature checked 10-16 (theory/LITERATURE.md): no prior 2D first-order result found.

## Phase 5 (submitted 10-17): pre-submission reinforcement (Matias: "start on items 2-4")
Three experiments, all in make_tasks.py, account torch_pr_1342_general, partition cs:

| exp | job | runs | what | est core-h |
|---|---|---|---|---|
| E9_alpha | 20180445 (16 c, 4 h) | 120 | q* ∝ α: E4 strip fronts (Lx 2048 × Ly 4096, T 1e5) at α ∈ {0.25, 0.5, 1}, p ∈ {0.98, 0.99}, q = x·2α·q*(p, α 0.5) with x ∈ {0.9 … 1.1}, 4 seeds. α 0.5 is rerun so all three share one protocol | 35 |
| E7_binder_rest | 20180442 (16 c, 4 h) | 100 | the rest of the E7 grid at L 512–2048 (seeds continue E6's): near points (0.98, 1.02), (0.99, 1.01), (0.995, 1.05), T 1e6; far (0.98/0.99/0.995, 1.2) and control (0.9, 1.03), T 5e5 | 38 |
| E7_binder_rest_L4096 | 20180443 (44 c, 8 h, THREADS=2) | 22 | the same points at L 4096 (5 seeds near, 4 where E6 has one; 2 far) | 113 |
| E8a_long | 20180446 (48 c; cancelled at 2:38, see restructure below) → 20187060 (17 c, 8 h), dnsim_v13 | 48 | E8a to T = 1e6 on L 8192 (no escapes before 1e6): (0.98, 1.02/1.05), (0.99, 1.02/1.05), (0.995, 1.02 = below-threshold control / 1.05). 8 sub-runs per point (k1–k8, own seeds), each ≤ 375k replicas or 5 survivors to 1e6 ((0.98, 1.05): T 3e5, 10 survivors) | 126 used + ~85 (≤ 136 at the 8 h limit) |

- **E8a_long question:** does P_s plateau after 1e5 just above the bulk threshold? Voter / GV
  gives P_s(1e6)/P_s(1e5) ≈ 0.12 (ln t / t); a first-order point above threshold gives a ratio
  near 1, and the below-threshold control should fall below 0.12. Local smoke test (v1.3,
  L 8192, T 1e6, 4000 replicas at p 0.99, 1.05 q*): 3.4 ms per replica, 2.2e9 updates/s.
- **Analysis:**
  - `python3 spread.py --exp E8a_long --brief $SCRATCH/dn/runs` pools the k sub-runs (stacked
    jackknife groups) and adds the 1e6/1e5 ratio.
  - `python3 binder.py $SCRATCH/dn/runs` picks up E7_binder_rest* automatically.
  - `python3 fronts.py $SCRATCH/dn/runs/E9_alpha` gives v(q) per (p, α), pooled over seeds.
- **Resubmit** a timed-out array with the same sbatch line; all runs resume from checkpoints.

### E9_alpha result (10-17): q* ∝ α
Zero crossing of the seed-mean v(q), from a linear fit through x = 0.95, 1.0, 1.05. The errors
are a seed bootstrap only. The α 0.5 crossings reproduce E4 to 0.2% (p 0.98) and 1.7% (p 0.99),
so allow about ±2% systematic. Crossing fit: `fronts.py` table, then a 3-point linear fit per (p, α).

| p | α | predicted 2α·q*(0.5) | measured q* | ratio | c = q*/(α(1−p)) |
|---|---|---|---|---|---|
| 0.98 | 0.25 | 0.001260 | 0.001274 ± 0.000009 | 1.011 | 0.255 |
| 0.98 | 0.5 | 0.002520 | 0.002515 ± 0.000019 | 0.998 | 0.252 |
| 0.98 | 1 | 0.005040 | 0.005074 ± 0.000030 | 1.007 | 0.254 |
| 0.99 | 0.25 | 0.000635 | 0.000621 ± 0.000018 | 0.977 | 0.248 |
| 0.99 | 0.5 | 0.001270 | 0.001292 ± 0.000007 | 1.017 | 0.258 |
| 0.99 | 1 | 0.002540 | 0.002553 ± 0.000009 | 1.005 | 0.255 |

- Fitting q* ∝ α^β over α 0.25 → 1 gives β = 0.997 (p 0.98) and 1.02 (p 0.99).
- c = 0.248–0.258, mean 0.254. That sits on κ_front = 0.2553 and at the lower edge of the
  perturbation-theory bracket [0.255, 0.265] (theory/README.md).
- A 5-point fit gives c = 0.252–0.260, within 2% of the 3-point values.
- All 4 seeds agree on the sign of v at x = 0.9 and 1.1 everywhere, and at 1.05 except one seed
  at (p 0.99, α 0.25) with v = −0.00002. That point is the noisiest (its x = 1.0 mean is above
  its x = 1.05 mean), hence its larger error. Below x = 1 the strip interior often orders before
  the fit window ends (t_order column), so v there is noisier.
- This replaces the 1-seed local α check (c = 0.244, 0.255): same answer, now with 4 seeds,
  two p values and a common protocol.

### E8a_long restructure (10-17 12:35)
- Job 20180446 held 48 cores, but 31 sub-runs were done after 2.5 h. The (0.98, 1.05) sub-runs
  could not have finished in the 12 h limit:
  - Survivors there spread ballistically (v = 0.0027), so one survivor to 1e6 costs about
    (2v)²T³/3 updates. All 8 were still on their first survivor after 2.6 h.
  - The docstring's "~25 min per survivor" assumed diffusive growth (box area ∝ t), which was
    wrong for x = 1.05.
  - Letting the job run would have charged ~576 core-h.
- Cancelled at 2:38. Two changes:
  - (0.98, 1.05) moves to T 3e5 with 10 survivors per sub-run (E8A_LONG_SHORT, names `_T3e5_k*`).
    Its partial T1e6 files are in `E8a_long/abandoned_T1e6/`.
  - The other 9 unfinished sub-runs resume from their between-replica checkpoints. Their
    in-progress replicas are replayed, losing about 30 core-h.
- Resubmitted as job 20187060 (`dn-e8l2`, 17 cores, 8 h):
  `sbatch -J dn-e8l2 --array=0 -c 17 -t 08:00:00 --export=ALL,EXP=E8a_long,DNSIM=$HOME/dn/bin/dnsim_v13 campaign.slurm`.
- `spread.py --brief` now treats T ≥ 3e5 as long. The ratio column is P_s(T)/P_s(1e5); voter
  gives 0.37 at T = 3e5.
- Snapshots of a running reps file cover only the replicas finished so far. The replica in
  progress is excluded, and it is more likely a survivor, so P_s from a running sub-run is biased
  low by up to one survivor per sub-run. Use only finished sub-runs.

### E8a_long interim (10-17, 31 of 48 sub-runs done)
| p | q/q* | vs bulk threshold | reps | alive at T | P_s 1e5 → 3e5 → 1e6 | d_eff at 3e4, 1e5, 2.5e5 | η_eff same t | P_s(1e6)/P_s(1e5) |
|---|---|---|---|---|---|---|---|---|
| 0.98 | 1.02 | at / ≤ 2% above | 3.0M | 9 | 2.07e-5 → 4.3e-6 → 3.0e-6 | 1.51, 1.39, 1.21 | −0.56, −0.43, +0.04 | 0.145 ± 0.045 |
| 0.99 | 1.02 | ~1% above | 3.0M | 3 | 4.4e-5 → 6.0e-6 → 1.0e-6 | 1.36, 1.68, 2.10 | −0.43, −0.82, −1.22 | 0.023 ± 0.013 |
| 0.995 | 1.02 | below (control) | 3.0M | 0 | 6.2e-5 → 8.0e-6 → 0 | 1.24, 1.51, 2.26 | −0.41, −0.77, −1.91 | 0 |
| 0.995 | 1.05 | ~1.5% above | 2.86M (7/8 done) | 14 | 9.5e-5 → 2.1e-5 → 4.9e-6 | 1.10, 1.28, 1.34 | −0.23, −0.38, −0.33 | 0.051 ± 0.014 |

Reading (interim):
- **The plain ratio test doesn't work this close to threshold.** The README predicted a ratio
  near 1 above threshold, but that needs the survivors at 1e5 to be past the critical droplet
  already. Within ~1–2% of the bulk threshold they are not. R_c ∝ 1/(q − q*) is large there, so
  P_s keeps falling after 1e5, faster than the voter curve.
- **The turnover has started at (0.98, 1.02).** δ_eff peaks and falls (1.51 → 1.21), and η_eff
  crosses zero, so the interface population has stopped shrinking. P_s(1e6)/P_s(3e5) = 0.69,
  against 0.33 for voter. This is the overshoot-then-plateau that nucleation predicts, now on the
  1e5–1e6 scale. Only 9–13 survivors are behind it.
- **At ~1% above threshold the decay is still accelerating by 1e6.** At (0.99, 1.02), δ_eff
  reaches 2.1. A seed below R_c dies faster than a voter seed, which fits a large R_c.
  - The model is monotone here (p/8 ≥ (1−p)α), so P_s(t; q) ≥ P_s(t; q_c) for q ≥ q_c.
  - A GV curve decays with δ_eff ≤ 1. So δ_eff > 1, if sustained, needs q_c above this q or no
    GV point at all.
  - E6 has the bulk active from 1.01 q*. A continuous transition would put q_c at or below that,
    which leaves the barrier as the explanation.
  - One decade with δ_eff rising from 1.4 to 2.1 is suggestive, not a proof.
  - P_s(1e6) = 1.0e-6, against 1.9e-5 on the p = 1 voter curve.
- **The control dies out completely:** 0 of 3e6 at 1e6, with δ_eff 2.3.
- **The order of the turnover times follows q − q\*:** 5.0e-5 at (0.98, 1.02), 3.1e-5 at
  (0.995, 1.05), 2.5e-5 at (0.99, 1.02). The x = 1.05 points at p 0.98 and 0.99 (1.3e-4 and
  6.4e-5) are still running.
- P_∞ > 0 is still not shown at any point within 2% of threshold. That needs t ≫ R_c²/D, beyond
  1e6. The x = 1.05 points are the test of an actual plateau.

### E7 rest final (10-17 13:00; jobs 20180442 and 20180443 COMPLETED)
`binder.py` (pooled with E6 and E7 tier 1; p 0.9 added to its QSTAR). χ_L for L = 512, 1024,
2048, 4096, with x2 = χ_4096/χ_2048:

| p | q/q* | χ_L | x2 at 4096 | U4 at 4096 | τ_m (L = 1024 → 4096) |
|---|---|---|---|---|---|
| 0.9 | 1.03 | 323, 353, 313, 412 | 1.32 | 0.05 ± 0.08 | 5–6e3 |
| 0.98 | 1.02 | 2612, 3667, 2537, 2888 | 1.14 | 0.05 ± 0.04 | 4.6e4 → 3.9e4 |
| 0.98 | 1.2 | 401, 499, 448, 441 | 0.98 | 0.03 ± 0.01 | 0.7–1.5e4 |
| 0.99 | 1.01 | 4176, 8469, 14797 ± 8933, 12213 ± 2675 | 0.83 | 0.19 ± 0.07 | 5.3e4 → 7.0e4 |
| 0.99 | 1.02 | 3128, 5291, 8340, 6357 | 0.76 | −0.16 ± 0.20 | 5.0e4 → 5.5e4 |
| 0.99 | 1.05 | 2660, 3773, 4370, 3546 | 0.81 | 0.02 ± 0.07 | 4.5e4 → 2.9e4 |
| 0.99 | 1.2 | 797, 1092, 924, 846 | 0.92 | −0.04 ± 0.06 | 1–1.5e4 |
| 0.995 | 1.05 | 4164, 8142, 10430 ± 2236, 18632 ± 4362 | 1.79 | 0.23 ± 0.08 | 6.2e4 → 8.0e4 |
| 0.995 | 1.2 | 2025, 2554, 2939, 2521 | 0.86 | 0.10 ± 0.03 | 2.0e4 → 2.5e4 |

Ising references: x2 = 3.36, U4* = 0.61 → 2/3 when ordered.

Reading:
- **At every p ≤ 0.99 point and at all far points, χ_L is flat from L 1024–2048 on**, with U4 → 0.
  So m has a finite correlation length and there is no ordered active window. That includes
  δ = 0.01 at p 0.99.
- **τ_m barely changes with L** at fixed (p, q): at most +35% from 1024 to 4096, and it falls at
  half the points. An Ising critical point would give τ ∝ L^2.17, i.e. ×20 over that range. m acts like a slow
  Ornstein–Uhlenbeck mode with an L-independent relaxation time. That fits a weakly broken voter
  conservation law: τ_m and χ grow as p → 1, but not with L.
- **(0.995, 1.05) is the one point not clearly flat.** x2 = 1.79 at 4096. The per-run χ scatter
  is large (τ_m 8e4 against a 8e5 window, n_eff 29), so roughly 1.8 ± 0.6. That is ~1.3σ from 1
  and ~2.6σ below Ising. U4 is 0.23 ± 0.08 at both 2048 and 4096 (Ising: 0.61 rising to 2/3).
  ⟨|m|⟩ falls with L (0.042 → 0.028). Early and late halves of the window agree (χ 18.1k vs
  19.2k at 4096), so it is not still equilibrating.
  - Verdict: no sign of Ising ordering, but at p 0.995 χ alone cannot exclude ξ ≳ 1000 near
    threshold. The p 0.995 first-order claim rests on E8 (walls and seeds), not E7.
- **The p 0.9 control** is flat at χ ≈ 350.

### E8a_long results (10-17 20:40; 46 of 48 sub-runs done)
- Job 20187060 hit its 8 h limit with (0.99, 1.05) k1 and k3 each one survivor short of 5. They
  are resubmitted as 20206748 (`dn-e8l3`, 2 cores, 5 h). k3 has one escaped replica (a cluster
  reached the L 8192 box edge before 1e6, counted as a survivor); this was the only escape.
- The (0.99, 1.05) row below uses the 6 finished sub-runs only, because running ones are biased
  low. All 8 give 0.442 ± 0.058.

| p | q/q* | vs bulk threshold | reps | alive at T | P_s 1e5 → 3e5 → 1e6 | d_eff 3e4, 1e5, 2.5e5 | η_eff same t | P_s(T)/P_s(1e5) | voter |
|---|---|---|---|---|---|---|---|---|---|
| 0.98 | 1.02 | ≤ 2% above | 3.0M | 9 | 2.07e-5 → 4.3e-6 → 3.0e-6 | 1.51, 1.39, 1.21 | −0.56, −0.43, +0.04 | 0.145 ± 0.045 | 0.12 |
| 0.98 | 1.05 | ~3–5% above | 1.32M (T 3e5) | 80 | 8.19e-5 → 6.07e-5 | 1.02, 0.60, – | 0.10, 1.00, – | 0.741 ± 0.051 | 0.37 |
| 0.99 | 1.02 | ~1% above | 3.0M | 3 | 4.4e-5 → 6.0e-6 → 1.0e-6 | 1.36, 1.68, 2.10 | −0.43, −0.82, −1.22 | 0.023 ± 0.013 | 0.12 |
| 0.99 | 1.05 | ~4% above | 0.58M (6/8) | 30 | 1.09e-4 → 6.37e-5 → 5.17e-5 | 0.98, 0.77, 0.40 | −0.08, 0.51, 1.33 | 0.476 ± 0.070 | 0.12 |
| 0.995 | 1.02 | below (control) | 3.0M | 0 | 6.2e-5 → 8.0e-6 → 0 | 1.24, 1.51, 2.26 | −0.41, −0.77, −1.91 | 0 | 0.12 |
| 0.995 | 1.05 | ~1.5% above | 2.86M | 15 | 9.6e-5 → 2.1e-5 → 5.2e-6 | 1.10, 1.28, 1.32 | −0.22, −0.36, −0.28 | 0.055 ± 0.014 | 0.12 |

Reading:
- **P_∞ > 0 is now shown at 1.05 q\* for p 0.98 and 0.99,** i.e. 3–5% above the bulk threshold.
  - P_s levels off at ≈ 6e-5 and ≈ 5e-5. From 3e5 to 1e6 at p 0.99 it falls only 0.81×, against
    0.33× for voter.
  - d_eff falls toward 0 and η_eff climbs toward 2 (compact ballistic growth): survivors past R_c
    simply spread.
  - E8a showed these points overshooting δ_eff > 1 at 1e4. So the full sequence holds: barrier
    (decay faster than voter), then plateau.
- **Within 1–2% of the bulk threshold the barrier phase lasts past 1e6:**
  - (0.99, 1.02) and (0.995, 1.05) still decay faster than voter.
  - (0.98, 1.02) has begun to turn over: δ_eff falls and η_eff ≈ 0.
- **P_∞ is very steep in δ:** between 1.02 and 1.05 q* it rises from below 3e-6 (an upper bound,
  still falling) to 5–6e-5, at least ×20. That fits an onset controlled by a critical droplet,
  not a power law. Three points per p are too few to fit the form.
- **The control** (0.995, 1.02, below the bulk threshold) dies out completely, 0 of 3e6.
- **For the paper:** seeds show a nucleation barrier followed by P_∞ > 0 at 1.05 q* (p ≤ 0.99).
  Closer to threshold only the barrier is visible within 1e6. P_∞ at p 0.995 is not shown.

Phase 5 cost, actually charged (sacct CPUTimeRAW): E9 23.4, E7_binder_rest 39.2,
E7_binder_rest_L4096 135.8, E8a_long 126.2 + 136.0 (+ ≤ 10 for 20206748) = **~470 core-h**.
20187060 also idled: the 8 T3e5 sub-runs finished in ~1.2 h, but the job held 17 cores for 8 h
(~55 core-h idle). Next time, put short and long runs in separate jobs.
- **Final (10-17 23:40, job 20206748 COMPLETED, 48/48 sub-runs):** (0.99, 1.05) with all 8
  sub-runs: 750,839 reps, 40 alive at 1e6. P_s 1.17e-4 → 6.39e-5 → 5.33e-5; d_eff 0.95, 0.81,
  0.43; η_eff −0.05, 0.53, 1.42. **P_s(1e6)/P_s(1e5) = 0.455 ± 0.057** (voter 0.12). This
  matches the 6-sub-run value. E8a_long total charged: 126.2 + 136.0 + 5.8 = 268 core-h (Phase 5 total ≈ 466).

## Phase 6 (submitted 10-17 ~23:55): asynchrony check + P_∞(δ) (Matias: "queue up something overnight")
**Why.** Referees will ask whether first order is an artefact of synchronous updating. The GV
and voter classifications are for continuous-time models.

**Method: lazy updating, no kernel change.**
- Each cell updates with probability r per step: the flip table is scaled by r and passed via
  `--ftab` (`make_tasks.lazy_ftab`).
- r → 0 is continuous time, and t·r is physical time. The voter part becomes a lazy voter model,
  whose dual is lazy coalescing walks.
- Checked locally: at r = 1, `--ftab` output is bit-identical to `--rule/--jflip`. The lazy
  table runs ~12% slower.

**Local pilot** (scratch `lazy/pilot`, Lx 256 × Ly 2048, 1 seed, p 0.98, T 3e4/r):
- v changes sign near 1.11 q*(r=1) at r 0.5, and about 1.15–1.25 at r 0.25 and 0.125.
- So asynchrony raises q* by roughly 10–20%. The grids below span 1.0–1.36 × q*(r=1).

| exp | job | runs | what | est core-h |
|---|---|---|---|---|
| E10_lazy_front | 20215052 (2 × 16 c, 4 h) | 112 | E4 strips (Lx 2048 × Ly 4096) to physical t 1e5. (p, r) = (0.98, 0.5/0.25/0.125) and (0.99, 0.25). 7 q from 1.0 to 1.36 q*(r=1), 4 seeds | 78 (+12%) |
| E10_lazy_bulk | 20215053 (10 c, 6 h) | 10 | E6-style bulk (random init, L 2048) to physical t 1e6 (T 4e6) at r 0.25, p 0.98. q = 1.06–1.30 q*(r=1), 2 seeds | 31 (+12%) |
| E8a_pinf | 20215054 (24 × 1 c, 6 h) | 24 | E8a_long protocol at p 0.98, x = 1.03, 1.04 (T 5e5) and 1.07 (T 3e5). 8 sub-runs × ≤ 375k reps / 10 survivors. One sub-run per array task, so no idle cores | 48 |

All jobs use `DNSIM=$HOME/dn/bin/dnsim_v13`. Resubmit timed-out arrays with the same line;
checkpoints resume.

**Analysis:**
- `python3 fronts.py $SCRATCH/dn/runs/E10_lazy_front`: q*(r) from the v sign change. fronts.py
  now parses `r` in run names, and v is per step.
- E10_lazy_bulk: ρ(t) tables (fronts.py bulk mode; its times are steps, so divide by r for
  physical time). Look for the lowest stationary plateau and its ρ.
- `python3 spread.py --exp E8a_pinf --brief $SCRATCH/dn/runs`.

**Theory.** `theory/cwalk.c` gained a lazy dual (env `CWALK_R`; r = 1 is bit-identical
single-threaded). Running `CWALK_R=r python3 vpt.py run --T <grid ×1/r>` gives the
parameter-free c(r) to compare with E10 (theory/vpt_run_r*.txt).

**Next if E10 confirms q*(r) and a finite plateau:** seeds and walls at r 0.25 (E8a/E8b
protocols at the measured q*(r)) for the full first-order signature under asynchrony.

### Phase 6 update (10-18 ~00:20): lazy theory vs pilot 2
- **Lazy perturbation theory** (`CWALK_R=r vpt.py run`, 1e6 trials, T to 2^16/r; theory/vpt_run_r0.25.txt,
  vpt_run_r0.125.txt). κ_a, κ_front and κ_th3 converge, and agree, much faster than at r = 1:
  - c(0.25) = 0.208 ± 0.004 and c(0.125) = 0.201 ± 0.004, against 0.255–0.265 at r = 1.
  - So the theory predicts **q\*(r) ≈ 0.8 q\*(1)**.
  - Lazy walks coalesce far more: P(≥ 3 classes) is 0.23 at physical time 65536 vs 0.53 at r = 1.
    Synchronous walks can swap places without meeting.
- **Pilot 2** (local, Lx 512 × Ly 2048, 2 seeds, physical T 3e4; scratch `lazy/pilot2`):
  - r = 1 crosses at x ≈ 1.03, so this geometry reads ~3% high.
  - r = 0.25: v per step = −0.0029, −0.0025, −0.0002, −0.0002, +0.0009, +0.0026 at x = 0.7, 0.8,
    0.9, 1.0, 1.2, 1.4. The crossing is ≈ 1.0–1.05, flat near zero over 0.9–1.0.
  - Pilot 1 (Lx 256, 1 seed) had put it at 1.15–1.25; the narrow strip misled.
  - **So the simulation says q\*(0.25) ≈ q\*(1)** and the theory says 0.8. That is a ~25% discrepancy
    if E10 confirms it.
  - Implication: the r = 1 agreement (c 0.26 vs 0.252) does not carry over to r < 1. Either it is
    partly fortuitous, or finite-ε corrections differ strongly between synchronous and lazy
    dynamics. In 2D they are powers of 1/ln(1/ε), and the walks' approach to T → ∞ differs a lot.
- **Extensions, so E10 brackets both answers:**
  - E10_lazy_front_lo, job 20215662 (16 c, 4 h): x = 0.76, 0.82, 0.88, 0.94, 3 seeds, same
    (p, r) points, 48 runs.
  - E10_lazy_bulk_lo, job 20215663 (4 c, 6 h): x = 0.94, 1.0.
  - Analyse both front dirs together: `fronts.py $SCRATCH/dn/runs/E10_lazy_front $SCRATCH/dn/runs/E10_lazy_front_lo`.
  - Overnight total ≈ 220 core-h.

### Phase 6 results (10-18 morning)
All five jobs completed (E10_lazy_bulk finished last, at ~5.5 h of its 6 h limit).

**E10 lazy fronts: q\*(r).** Physical v = (per-step v)/r. A linear fit of the seed v against x over
x ∈ [0.88, 1.12] gives the crossing x\* = q\*(r)/q\*(1). The x = 0.76 and 0.82 fronts collapse early
(t_fitmax 1–3e5 steps, npts 26–72), so their v is poorly determined. Including x = 0.82 moves x\* up
by 0–0.05, which is the systematic. Fit script: scratch `lazy/cross.py`.

| p | r | x\* (0.88–1.12 fit) | x\* (0.82–1.06) | c = q\*/(α(1−p)) | theory c(r) |
|---|---|---|---|---|---|
| 0.98 | 1 | 1 (E4/E9) | | 0.252 | 0.255–0.265 |
| 0.98 | 0.5 | 0.864 ± 0.010 | 0.903 ± 0.012 | 0.218 | 0.223 ± 0.003 |
| 0.98 | 0.25 | 0.858 ± 0.016 | 0.869 ± 0.010 | 0.216–0.219 | 0.208 ± 0.004 |
| 0.98 | 0.125 | 0.870 ± 0.012 | 0.923 ± 0.024 | 0.219 | 0.201 ± 0.004 |
| 0.99 | 0.25 | 0.945 ± 0.018 | 0.927 ± 0.014 | 0.236–0.240 | 0.208 ± 0.004 |

- **The fronts keep the first-order form under asynchrony.** v(q) crosses zero linearly with a
  finite slope. Behind the front, ρ_b = 0.64–0.67 at every (r, x), including the retreating fronts
  at x 0.76–0.82, where the active phase is metastable and loses.
- **q\* shifts down, as the theory says, but by less.** At p 0.98, q\*(r) ≈ 0.865 q\*(1) at all three r,
  so the shift saturates by r = 0.5. The theory gives 0.80 (r 0.25) and 0.78 (r 0.125).
- **So the pilot-2 "25% discrepancy" is gone.** It came from Lx 512 strips and T 3e4. At p 0.98 the
  theory is 4–8% low under asynchrony, against 3% high at r = 1.
- **Open: p dependence.** At p 0.99 the shift is only −5%, so c rises from 0.218 to 0.24 between
  ε = 0.02 and 0.01. The synchronous c is p-independent to 2% (E9), and the theory is the ε → 0
  limit, so this runs the wrong way. Either the p 0.99 fronts are not converged at physical T 1e5
  (εT = 1e3; at 0.94 the fronts moved < 300 cells), or the lazy finite-ε corrections are large.
  Test: p 0.99 (and 0.995) at r 0.25 to physical T 3e5 with 4+ seeds at x 0.88–1.0.

**E10 lazy bulk (r 0.25, p 0.98, L 2048, physical t to 1e6 = 4e6 steps, 2 seeds).** Bond density ρ
(column 3), seed means. The last column averages steps 2e6–4e6 (physical 5e5–1e6).

| x | ρ at physical 2.5e4 | 2.5e5 | 5e5–1e6 mean | seeds (5e5–1e6) |
|---|---|---|---|---|
| 0.94 | 0.173 | 0.181 | 0.176 | 0.178, 0.174 |
| 1.00 | 0.180 | 0.182 | 0.184 | 0.184, 0.184 |
| 1.06 | 0.187 | 0.192 | 0.187 | 0.187, 0.188 |
| 1.12 | 0.193 | 0.192 | 0.191 | 0.191, 0.191 |
| 1.18 | 0.193 | 0.194 | 0.194 | 0.194, 0.194 |
| 1.24 | 0.198 | 0.194 | 0.197 | 0.197, 0.197 |
| 1.30 | 0.198 | 0.199 | 0.200 | 0.200, 0.200 |

- None absorbed. At x 0.94, only ~8% above q\*(0.25), the active state holds a flat plateau
  ρ ≈ 0.176 over a 40× range of time. Over that range the synchronous voter fit
  (1/ρ = 0.184 ln t + 2.13) loses 15%. The lazy-voter constants are not measured.
- ρ falls smoothly and by only 12% from x 1.30 to 0.94, with no sign of vanishing toward q\*(r).
  That is the synchronous E6 picture (ρ 0.26 at 1.02 q\*).
- The absolute ρ is lower than in the synchronous runs. Lazy updating removes the simultaneous
  two-cell swaps, so it carries fewer interfaces.
- Not yet run under asynchrony: bulk below q\*(r), to show the collapse (the fronts at 0.76–0.82
  already show the active phase losing), and the lazy-voter ρ(t) baseline.

**E8a_pinf (p 0.98, seeds, P_s(T)/P_s(1e5) ± jackknife; voter: 0.37 at 3e5, 0.23 at 5e5, 0.12 at 1e6).**

| q/q\* | reps | alive at T | P_s(3e5) | ratio | reading |
|---|---|---|---|---|---|
| 1.02 (E8a_long) | 3.0e6 | 9 at 1e6 | 4.3e-6 | 0.145 ± 0.045 (1e6) | still decaying at 1e6; turnover begun |
| 1.03 | 3.0e6 | 22 at 5e5 | 8.3e-6 | 0.218 ± 0.041 (5e5) | voter-like decay to 5e5 |
| 1.04 | 2.17e6 | 74 at 5e5 | 3.6e-5 | 0.561 ± 0.045 (5e5) | above voter, still decaying |
| 1.05 (E8a_long) | 1.32e6 | 80 at 3e5 | 6.1e-5 | 0.741 ± 0.051 (3e5) | levelling |
| 1.07 | 4.7e5 | 80 at 3e5 | 1.7e-4 | 0.930 ± 0.028 (3e5) | plateau: P_∞ ≈ 1.7e-4 |

- P_∞ turns on steeply above q\*: ≈ 1.7e-4 at 1.07, ≈ 6e-5 at 1.05, ≲ 3.5e-5 at 1.04, and no plateau
  by 5e5 at 1.03. The time at which survival levels off grows rapidly as δ → 0.
- This fits nucleation: a critical radius R_c ∝ 1/δ, reached by voter-like spreading in t_c ∝ R_c²,
  so P_∞ ∝ δ² up to logs. From 1.05 to 1.07 the ratio is ≥ 2.8, against (7/5)² = 2.0 if δ = x − 1. It
  matches if the seed threshold sits ~1–2% above the front q\* (δ = x − 1.02 gives 2.8). That is within
  the q\* systematic, so do not quote an exponent from three points.

**Cost.** E10f 121.7, E10fl 40.0, E10bl 20.0, E10b ≈ 55, E8p 40.6: **≈ 277 core-h** against the 220
estimate. The lazy runs ran at 8.7–9.5e8 cell-updates/s per core, against the 1.34e9 assumed
(1.5e9/1.12). For later lazy jobs use `--rate 9e8`.

**Next.**
1. Seeds (E8a protocol) and walls (E8b) at r 0.25, p 0.98, at the measured q\*(0.25) = 0.865 q\*(1),
   to get the full first-order signature under asynchrony.
2. The p 0.99 lazy convergence test (above).
3. Done: lazy theory at r 0.5 gives c = 0.223 ± 0.003 (theory/vpt_run_r0.5.txt), within 2% of the
   simulation. The theory's r-dependence is right from r 1 to 0.5. Below 0.5 the theory keeps
   falling (0.208, 0.201) while the simulation stays flat (0.218 ± 0.004). The theory's drop from
   r 0.5 to 0.125 is 10%, against 0.7 ± 1.8% in the simulation, a ~4σ difference. Unexplained.

## Phase 7 (submitted 10-18 morning): E11, the first-order signature under asynchrony (Matias: "Go for it")
All three are 1-core array tasks (one run each, 2G), `DNSIM=$HOME/dn/bin/dnsim_v13`, estimated at
`--rate 9e8`. q\*(0.5) = 0.864 × 0.00252 = 0.0021773 (`make_tasks.QSTAR_LAZY`). r 0.5 is used because
E10 showed the q\* shift already saturated there, and it costs half as much as r 0.25.
Checked locally: `--reps` runs with `--ftab`, and at r = 1 the output is bit-identical to `--rule/--jflip`.

| exp | job | runs | what | est core-h |
|---|---|---|---|---|
| E11_lazy_seed | 20249854 (24 × 1 c, 8 h) | 24 | E8a seeds at r 0.5, p 0.98, L 2048, physical t 1e5 (T 2e5). x = 0.95, 1.0, 1.05, 1.1, 1.2 × q\*(0.5), plus the lazy voter (p 1, q 0). 4 sub-runs per point, each ≤ 750k replicas or 75 survivors (25 at x ≥ 1.1) | 72 (rough) |
| E11_lazy_wall | 20249855 (10 × 1 c, 5 h) | 10 | E8b walls at r 0.5, p 0.98, Lx 1024 × Ly 4096, physical t 1e6 (T 2e6). x = 0.95, 1.0, 1.05, 1.1 × q\*(0.5), plus the lazy voter wall, 2 seeds | 26 |
| E11_lazy_front_p99 | 20249856 (20 × 1 c, 6 h) | 20 | E10 strips at p 0.99, r 0.25 to physical t 3e5 (T 1.2e6). x = 0.84, 0.88, 0.92, 0.96, 1.0 × q\*(0.99, r 1), 4 seeds. Is the p 0.99 crossing (0.93–0.95 at 1e5) converged? | 62 |

**Analysis.** spread.py and wall.py now read lazy runs: physical time t·r, and q/q\* relative to
q\*(r). Both still reproduce the E8a_pinf and E8b tables exactly.
- `python3 spread.py --exp E11_lazy_seed --brief $SCRATCH/dn/runs`. Compare δ_eff against the lazy
  voter row, not the synchronous voter.
- `python3 wall.py $SCRATCH/dn/runs --exp E11_lazy_wall`. W and rc against the lazy voter wall. The
  RHO_HALF of 0.3 assumes an active interface density of ~0.6; lazy fronts give ρ_b 0.64–0.67, so it holds.
- `python3 fronts.py $SCRATCH/dn/runs/E11_lazy_front_p99` (v per step; physical v = v/r). Also fit
  v over an early and a late window to see whether the crossing drifts with time.

### Phase 7 results (10-18, 16:45): first order survives asynchronous updating
All 54 tasks completed. Charged 95.0 core-h, against the 160 estimate: survivors and walls were
cheaper than assumed. Seeds took 0.4–1.4 h per sub-run, walls 1.3–2.8 h, fronts 2–3 h.

**E11_lazy_seed (r 0.5, p 0.98; spread.py, physical t; q/q\* relative to q\*(0.5)).**

| x | reps | alive at 1e5 | P_s(1e5) | P_s/voter | δ_eff at 3e3, 1e4, 3e4 | η_eff at 3e4 |
|---|---|---|---|---|---|---|
| lazy voter | 2.14M | 300 | 1.40e-4 | 1 | 0.92, 0.94, 0.87 | −0.02 |
| 0.95 | 3.0M | 53 | 1.77e-5 | 0.13 | 1.03, 1.19, 1.57 | −0.93 |
| 1.00 | 3.0M | 206 | 6.87e-5 | 0.49 | 0.98, 1.03, 1.09 | −0.31 |
| 1.05 | 1.67M | 300 | 1.80e-4 | 1.29 | 0.92, 0.89, 0.78 | +0.18 |
| 1.10 | 0.27M | 100 | 3.73e-4 | 2.7 | 0.87, 0.83, 0.51 | +0.73 |
| 1.20 | 0.10M | 100 | 1.00e-3 | 7.1 | 0.70, 0.57, 0.25 | +1.30 |

- At and below q\*(r), δ_eff is above 1 and rising, and η_eff is negative: the nucleation barrier,
  as in synchronous E8a (p 0.98, 1.02 q\*: 1.16 → 1.30 → 1.54).
- Above q\*(r), P_s turns toward a plateau: at 1.2, P_s(3e4) = 1.10e-3 and P_s(1e5) = 1.00e-3.
- P_s(1e5) equals the voter value at x ≈ 1.04 (log interpolation between 1.0 and 1.05). The
  synchronous value is 1.060.

**E11_lazy_wall (r 0.5, p 0.98; wall.py --exp E11_lazy_wall, physical t).** W | rc:

| x | 1e4 | 1e5 | 3e5 | 1e6 | a_W (1e5–1e6) | v |
|---|---|---|---|---|---|---|
| 0.95 | 168, 0.46 | 350, 0.43 | 400, 0.42 | 352, 0.41 | 0.07 | 0 ± 8e-6 |
| 1.00 | 164, 0.46 | 594, 0.42 | 613, 0.41 | 1038, 0.44 | 0.26 | 2.7e-4 |
| 1.05 | 170, 0.49 | 660, 0.41 | 1726, 0.42 | 2045 (full), 0.48 | 0.97 | 2.7e-3 |
| 1.10 | 193, 0.48 | 957, 0.45 | 2044 (full), 0.52 | 2032, 0.54 | — | 3.0e-3 |
| lazy voter | 182, 0.44 | 550, 0.42 | 1097, 0.41 | 1186, 0.40 | 0.32 | 2.5e-4 |

- **Below q\*(r), a stationary wetting slab.** W stops at ≈ 350–400 from 1e5 to 1e6, while the lazy
  voter wall doubles. This is synchronous E8b at 0.98 q\*.
- **At q\*(r)** the slab creeps at the voter-wall rate. The synchronous runs showed the same: ×2.5
  per decade, v 3e-4.
- **Above q\*(r)** the fronts are ballistic. At 1.10 the wall v (3.0e-3) matches the E10 front v at
  the same q (2.7e-3). At 1.05 the wall v is ~2× the front's (1.3–1.5e-3), because over 1e5–3e5
  the slab edges also spread voter-like (the voter wall widens ~1.4e-3 over the same interval).
- rc is 0.41–0.54, below the synchronous 0.46–0.71. That is the same lower lazy interface density
  seen in the bulk.

**E11_lazy_front_p99 (p 0.99, r 0.25, to physical t 3e5).** fronts.py (one fit from t 2000 steps to
the end) gives per-step v −0.00038, +0.00008, +0.00016, +0.00049, +0.00060 at x 0.84, 0.88, 0.92,
0.96, 1.0. The crossing is ≈ 0.87, against 0.93–0.95 from E10 at 1e5. Physical v in windows
(`vwin.py`, seed mean ± sem of per-run fits):

| x | [5e2, 1e5] | [2.5e4, 1e5] | [1e5, 3e5] |
|---|---|---|---|
| 0.84 | −0.0036 ± 0.0008 | −0.0040 ± 0.0007 | −0.0006 ± 0.0010 |
| 0.88 | −0.0017 ± 0.0007 | −0.0019 ± 0.0012 | +0.0009 ± 0.0006 |
| 0.92 | +0.0001 ± 0.0002 | +0.0009 ± 0.0003 | +0.0021 ± 0.0013 |
| 0.96 | +0.0012 ± 0.0003 | +0.0012 ± 0.0001 | +0.0021 ± 0.0008 |
| 1.00 | +0.0023 ± 0.0006 | +0.0025 ± 0.0007 | +0.0020 ± 0.0004 (strip filling) |

- **The p 0.99 anomaly was a finite-time effect.** Fronts just above q\*(r) retreat for ~1e5
  (physical) before advancing, so the crossing in a ≤ 1e5 window sits at ≈ 0.92. In [1e5, 3e5] it
  is ≈ 0.85–0.88, the same as p 0.98 (0.86–0.87).
- So c_sim(r 0.25) ≈ 0.86 × 0.254 ≈ 0.22 at both p, against the theory's 0.208.
- At p 0.98 the transient should be roughly half as long (∝ 1/ε). There the [2.5e4, 1e5] window
  agrees with the full fit (x 0.88: +0.0001 at r 0.25, +0.0008 at r 0.5). The p 0.98 x\* could
  still drift down by ≲ 0.02 asymptotically.

**Verdict.** Under asynchronous updating (r 0.5 for seeds and walls, r 0.25 for bulk and p 0.99
fronts, r 0.125–0.5 for p 0.98 fronts), all four first-order signatures reappear at a shifted
q\*(r) ≈ 0.86 q\*(1):
1. linear front velocity crossing, with a metastable active phase;
2. finite bulk plateau;
3. a seed nucleation barrier, then P_∞ > 0;
4. a stationary wetting slab below q\*(r) and ballistic slabs above.

The first-order result is not an artefact of synchronous updating. The parameter-free theory
predicts the shift's sign and most of its size (c 0.26 → 0.208, against 0.252 → 0.22 measured).

## Phase 8 (submitted 10-27 ~20:00): the clean model, flagship p = 0.5 (Matias: "go for it", aim PRE / J. Stat. Mech.)
**Why.** At p = 0 the model has exponentially many frozen configurations (all j ∈ {0,1,4,5};
entropy ≈ 0.32 nats/site), dynamically irrelevant for q > 0 but enough to make "two absorbing
states" false. For p > 0 every j ≥ 1 flips with probability ≥ p j/8, and on the Moore lattice
every non-uniform configuration has a cell with j ≥ 3 (so any rule with w(j) > 0 for j ≥ 3
has exactly two absorbing states). The paper's flagship therefore moves to p = 0.5 (q* = 0.0459,
ρ jump 0.29), which until now had only fronts (E4, 2 seeds) and bulk (E4_voter_bulk). Phase 8
gives p = 0.5 the full first-order signature, fills the phase diagram between p = 0 and 0.5,
tests front isotropy, measures the wetting of the absorbing–absorbing wall, and checks a second
two-absorbing-state rule. Details of the in-house checks: research_notes/.../in_house_checks.md.

Local pilots (10-27, scratchpad p8/):
- q*(p) from strips (Lx 512, 1 seed): ≈ 0.021 (p 0.01), 0.021 (0.02), 0.029 (0.05), 0.049 (0.25);
  p 0.1 refined (Lx 1024, 2 seeds): v = −0.0018 at 0.0365, +0.0011 at 0.037 → q* ≈ 0.0368,
  q*·dv/dq ≈ 0.22.
- Droplets at p 0.5 (L 2048, 200 replicas, T 2e4): P_surv = 1/2 near R ≈ 20 (x 1.05), 11 (1.1),
  8 (1.2), 4 (1.4). 2×2 seeds: 0 of 200 survive at every x ≤ 1.4.
- Walls at p 0.5 (Lx 1024 × Ly 1024, T 1e5, 1 seed): interface cells per unit wall length
  I ≈ 8.5 (q = 0), 18 (x 0.8), 26 (0.9), 35 (0.95), 50 (0.98), 67–72 (0.99); at x = 1.0 I grows
  without bound (30 → 176 over 1e3–1e5). Excess Γ = I − I(q=0) ≈ 9.5, 17.5, 26.5, 41.5, 60:
  roughly Γ ∝ δ^−0.5 over δ = 0.2 → 0.01 (single seed, not converged near δ = 0.01).
- E4_voter_bulk at p 0.5 (L 2048, 4096): τ_half ≈ 2800 (x 0.8), 4000 (0.9), 8900 (0.95),
  identical at both L (multi-droplet regime); plateau 0.290 at 1.002 q*.
- Tilted bands (tilt.py; L 1024, 1 seed, T 1.5e4): v = +0.0087, +0.0088, +0.0088 at x 1.05 for
  (a, b) = (1, 0), (1, 1), (2, 1). `--load` + `--resume` is bit-identical after an interruption.

| exp | job | runs | what | est core-h |
|---|---|---|---|---|
| E12_phase | 21159521 | 133 | E4 strips at p 0.01, 0.02, 0.05, 0.1, 0.25 (x 0.9–1.1, 3 seeds) and p 0.5 (x 0.96–1.04, 4 seeds) | 40 |
| E12_phase_bulk | 21159522 | 50 | bulk L 2048, T 1e5 at the new p, x 0.97–1.1, 2 seeds: the jump vs p | 8 |
| E13_bulk05 | 21159523 (0-2) | 1020 | p 0.5 bulk FSS: L 64–4096 (32…2 seeds), x 0.9–1.1, T 1e6 (L ≤ 1024), 3e5, 2e5; plateau vs L, lifetimes vs L (single- vs multi-droplet), χ_L | 86 |
| E14_drop05 | 21159525 | 48 | p 0.5 critical droplets (--reps, disk:0.5:R), x 1.02–2.0, 500–1000 replicas per (x, R) | 9 (uncertain) |
| E15_orient | 21159526 | 252 | tilted bands (a, b) = (1,0), (2,1), (1,1), L 2048, T 1e5, p 0.1 / 0.5 / 0.9, x 0.97–1.03, 4 seeds (`tilt.py`) | 38 |
| E16_wall | 21159527 | 80 | p 0.5 walls, x 0.8–0.995, Lx 256 / 1024 / 4096 × Ly 1024, T 1e6, + bare wall; p 0.9 walls x 0.9–0.99 | 21 |
| E17_eta_front / _bulk | 21159528 / 21159530 | 21 / 28 | single-noise rule w = α[j ∈ {3,6,7,8}] + η[j ∈ {1,2,4,5}], α 0.5, η = 0.0125 × (0.9–1.1): strips (3 seeds) and bulk (L 1024, 2048) | 12 |

All with `DNSIM=$HOME/dn/bin/dnsim_v13`. E15 checkpoints: `tilt.py make` (seed 100 s + 7), in
`$SCRATCH/dn/runs/E15_orient/init/`. Analysis: `fronts.py` (E12, E17), `tilt.py table` (E15),
`wall.py` (E16; QSTAR needs p 0.5), new scripts for E13 (lifetimes) and E14 (P_surv(R)).

Added the same evening:
- `E16_wall_ref` (job 21160700, 2 runs): the synchronous pure-voter wall (p 1, q 0) in the E8b
  geometry (Lx 1024 × Ly 4096, T 1e6), the reference E8b only had from a local one-seed run.
- `E18_lazy05_front` / `_bulk` (jobs 21163452 / 21163461, 24 + 10 runs, ~18 core-h): asynchrony at
  the flagship. Local pilot (Lx 512, 2 seeds, physical t 1.5e4): at p 0.5, r 0.25 the per-step v is
  −0.0034, +0.0018, +0.0039 at x = 0.85, 0.90, 0.95 × q*(1), so q*(0.25) ≈ 0.88 q*(1) (p 0.98: 0.86).
  Strips to physical 5e4 (3 seeds, x 0.84–0.96), bulk L 2048 to physical 1e5 (x 0.86–0.96).
- Early results (10-27 ~20:30): E17 η* ≈ 0.0129–0.0130 (strip v −0.0022, +0.0020, +0.0054 at
  η 0.01275, 0.013125, 0.01375; bulk ρ 0.252 at 0.013125, half-lives L-independent below). E12 bulk:
  stationary ρ just above q* = 0.245 (p 0.01), 0.256 (0.02), 0.257 (0.05), 0.288 (0.1), 0.307 (0.25).
  E14: R_c (P_surv = 1/2) ≈ 43, 21, 12.5, 7.9, 4.8, 2.6 at x = 1.02, 1.05, 1.1, 1.2, 1.4, 2.0; R_c·v
  0.18–0.26. E15 p 0.9: crossings x* = 1.0106(7), 1.0104(3), 1.0114(5) for (1,0), (2,1), (1,1).
- Correction: the E8a "monotone coupling" argument (P_s(t; q) nondecreasing in q for p ≥ 0.8) does
  not hold. Attractiveness orders configurations, but raising q adds flips in both directions at
  j = 2 (births of dead sites with 2 live neighbours and deaths of live sites with 2 dead ones), so
  the q-dynamics are not stochastically ordered. δ_eff > 1 and rising stays as physical evidence of
  a nucleation barrier, not as a proof.

### Phase 8 results (10-27, 21:30; all jobs COMPLETED except one E14 run, x 1.02 R 90, still adding replicas)
Charged ≈ 165 core-h (project total ≈ 1750). Figures and tables: `paper/figs/*.py`, manuscript
`paper/main.tex` (all numbers below are in it).
- **E12 transition line** (strip crossing, 3–4 seeds): q* = 0.02070 (p 0.01), 0.02280 (0.02), 0.02874
  (0.05), 0.03680 (0.1), 0.04932 (0.25), 0.04596 ± 0.00005 (0.5; E4 + E12). q*·dv/dq through q* over
  (1, 1.065] q*: 0.10, 0.13, 0.13, 0.14, 0.19, 0.22. Bulk just above: ρ+ = 0.245, 0.256, 0.257, 0.288,
  0.307, 0.302. q* rises steeply off p = 0 (0.0184 → 0.0207 at p 0.01), peaks near p 0.25–0.5.
- **E13 p 0.5 bulk FSS** (q* = 0.04596): stationary ρ independent of L (256–4096) to 1e-4: 0.3022–0.3023
  at 1.008 q*, 0.3111–0.3113 at 1.018 q*; at 0.998 q* ρ = 0.284 for 1e6 steps at all L ≥ 512 (within
  the q* uncertainty, so metastable or at coexistence). Half-life below q*: 3.4e3, 8.0e3, 1.6e4, 2.6e4,
  7.0e4 at 0.899, 0.949, 0.969, 0.979, 0.988 q*, L-independent for L ≥ 256 (multi-droplet). L ≤ 128 die
  sooner (L ~ ξ: whole-lattice fluctuations), L 64 even above q*; no L ≥ 256 run at ≥ 1.008 q* absorbed.
  τ ~ δ^−1.4 over δ 0.012–0.1 (power law rms 5% vs CNT exp(B/3δ)δ^−2/3 13%, B ≈ 0.06): R_c ≈ 70 at
  δ 0.012 is only 2–3 ξ, so the asymptotic essential singularity is out of reach. χ_L flat in L ≥ 512 at
  every x: 63 (1.098), 117, 250, 430 (1.008), ≈ 1300 (0.998): ξ grows toward q* but stays ≲ 40.
- **E14 droplets**: R_c = 2.6, 4.8, 7.9, 12.5, 21, 42 at q/q* = 2.0, 1.40, 1.20, 1.10, 1.05, 1.02;
  R_c = D/v + ΔR with D = 0.20, ΔR = 1.6 (R_c·v 0.21–0.26). 2×2 seeds survive with 0.014 (1.4 q*), 0.046 (2 q*).
- **E15 orientation**: zero crossings for (1,0), (2,1), (1,1) agree to 0.37% (p 0.1), 0.16% (0.5),
  0.15% (0.9). Ordering (1,0) > (1,1) > (2,1) everywhere, matching the 1/√(a²+b²) bias from the
  slow relaxation of the band interior (tilt.py uses ι_in fixed at t = 1000). A bulk-drift correction
  (paper/figs/orient_corr.py) over-corrects (band interiors do not nucleate like the bulk), so the
  paper quotes 0.4% as an upper bound on any orientation dependence.
- **E16 wetting at p 0.5**: bare wall I0 = 8.7 interface sites per unit length. Stationary excess
  Γ = I − I0 = 10.0, 17.7, 27.6, 36.8, 45.1, 63.5, 88 at δ = 1 − q/q* = 0.20, 0.10, 0.051, 0.031, 0.021,
  0.011, 0.0065; same for Lx 256 / 1024 / 4096 within 3%; stationary from 1e5 to 1e6. Local exponent
  → 0.55 for δ ≤ 0.05 (0.47–0.50 if q* is 0.15% lower): ψ = 0.5 ± 0.06, vs 1/3 for equilibrium complete
  wetting (2D, short-range). At δ 0.0065, Lx 256 runs got absorbed (walls met). Centred layer profile
  peaks at 0.53 (δ 0.0065) vs 0.74 bulk interface-site fraction.
- **E16_wall_ref**: the synchronous voter wall keeps widening (W ≈ 2000 by 1e6), unlike the 0.98 q*
  slabs of E8b.
- **E17 η rule**: η* = 0.01295 (strip), ρ = 0.252 at 1.013 η*, half-lives L-independent below.
- **E18 lazy at p 0.5** (r 0.25): q*(r)/q*(1) = 0.867 ± 0.007 (p 0.98: 0.86); bulk ρ = 0.211 stationary
  at 1.014 q*(r) to physical t 1e5, slow decay at 0.99 q*(r).
