# The lineage closure applied to the five-cell rule of Molofsky et al. (10-30)

    cc -O3 -march=native -fopenmp -o cwalk5 cwalk5.c
    python3 vpt.py selftest                                                   # z = 8 and z = 4 identities
    python3 vpt.py run --z 4 --trials 1000000 --seed 2 --T 0,1,2,...,65536    # -> patterns5.txt (15 s on 16 cores)
    python3 vpt.py molofsky --load patterns5.txt --boot 300 --csv molofsky.csv # -> vpt_run_molofsky.txt
    python3 plot_molofsky.py                                                  # -> molofsky.png

Files: `cwalk5.c` (five lazy coalescing walks), `patterns5.txt` (10^6 trials, seed 2),
`patterns5_1e7.txt` (10^7 trials, seed 3, cross-check), `vpt_run_molofsky.txt`,
`vpt_run_molofsky_1e7.txt`, `molofsky.csv`, `molofsky_1e7.csv`, `molofsky.png`. `vpt.py` now takes
the neighbourhood size `z` (default 8); `python3 vpt.py run --load patterns.txt` reproduces
`vpt_run_th3.txt` byte for byte.

## 1. The rule

Molofsky, Durrett, Dushoff, Griffeath & Levin, Theor. Popul. Biol. 55, 270 (1999): square lattice,
synchronous update, each site's next state depends on the number n of species-1 cells among the
five cells of its von Neumann neighbourhood including itself: P(next = 1 | n) = p_n with
p_0 = 0, p_5 = 1, p_3 = 1 − p_2, p_4 = 1 − p_1. In the language of the manuscript (flip
probability as a function of the number j of disagreeing neighbours among the four): a site in
state 0 with j disagreeing neighbours has n = j and flips with probability p_j; a site in state 1
has n = 5 − j and flips with probability 1 − p_{5−j} = p_j. So

    w(j) = p_j,   j = 0..4:   w(0) = 0,  w(1) = p_1,  w(2) = p_2,  w(3) = 1 − p_2,  w(4) = 1 − p_1,

the same for both states: the rule is exactly Z2-symmetric with two absorbing states. The voter
point is (p_1, p_2) = (1/5, 2/5), where w(j) = j/5. Writing w(j) = j/5 + g(j),

    g = (0, δ_1, δ_2, −δ_2, −δ_1),     δ_1 = p_1 − 1/5,   δ_2 = p_2 − 2/5,

since g(3) = 1 − p_2 − 3/5 = −δ_2 and g(4) = 1 − p_1 − 4/5 = −δ_1. The symmetry p_j + p_{5−j} = 1
together with j/5 + (5 − j)/5 = 1 gives g(j) + g(5 − j) = 0 for j = 1..4 (not g(j) + g(4 − j)).
The perturbation has two parameters, (δ_1, δ_2), and the physical region is the whole square
0 ≤ p_1, p_2 ≤ 1, so every direction around the voter point is allowed.

Their mean-field map (their eq. 1) is linear in the p_n, so its increment is exactly

    x_{t+1} − x_t = Σ_{m=1}^{4} C(5,m) g(m) x^m (1 − x)^{5−m}
                  = x(1−x)(1−2x) [5 δ_1 (1 − x + x²) + 10 δ_2 x(1−x)]  ≡  f_MF(x),

with J(0) = 1 + f_MF'(0) = 5 p_1 (their eq. 2) and J(1/2) = 1 + f_MF'(1/2) = (15 − 15 p_1 − 10 p_2)/8
(their eq. 3). Their "multiple outcome" regime (0, 1/2 and 1 all stable) is p_1 < 1/5 and
15 p_1 + 10 p_2 > 7: the wedge δ_1 < 0, δ_2 > −1.5 δ_1, whose tip is the voter point, between
the directions 90° and 123.7° measured from the +δ_1 axis. Their "phase separation" wedge is the
opposite one (p_1 > 1/5, 15 p_1 + 10 p_2 < 7).

## 2. The dual

At the voter point the probability that a site becomes 1 is n/5, which is exactly "copy a
uniformly chosen cell among the five, itself included". Since each site's new state depends only
on its own neighbourhood and its own randomness, the law of the automaton is the same whichever
way the per-site randomness is realised, and the voter-point automaton is the synchronous
five-cell voter model: ξ_{t+1}(x) = ξ_t(x + e_x), e_x uniform in {0, ±x̂, ±ŷ}, independent over
sites and times. Tracing a site's state backward gives a lazy walk that stays with probability
1/5 and steps to each of its four neighbours with probability 1/5; walks of different sites use
independent choices while apart and the same choice once they occupy the same site, so they
coalesce on contact and never split. `cwalk5.c` runs five such walks from the origin and its four
neighbours (uniform move drawn from a 32-bit word, bias 5/2^32) and records, at each checkpoint
T, the histogram of lineage-class patterns in the format of `cwalk.c`. The T = 1 histogram agrees
with exact enumeration of the 5^5 move combinations (all seven patterns within 2.1σ at 10^6
trials; P(all five distinct after one step) = 1458/3125 = 0.46656).

## 3. Reaction function: exact structure

With z = 4 the pattern polynomials of `vpt.py` give, for g antisymmetric as above,

    f(u) = δ_1 G_1(u) + δ_2 G_2(u),     G_i(u) = h(u) [a_i + c_i u(1−u)],    h(u) = u(1−u)(1−2u).

Checked in `vpt.py selftest` (z = 8 and z = 4): g linear in j gives f ≡ 0; one- and two-class
patterns cancel; three-class patterns (n_1, n_2) give h(u)[g(n_1) + g(n_2) − g(n_1 + n_2)];
four-class patterns (n_1, n_2, n_3) give h(u)[g(n_1) + g(n_2) + g(n_3) − g(n_1 + n_2 + n_3)];
the only five-class pattern (1,1,1,1), all five lineages distinct, gives exactly f_MF, the
mean-field quintic. Hence, with P_5(T) the probability of five classes and P(pattern) the
others,

    c_1 = −5 P_5,   c_2 = 10 P_5,   a_1 = 5 P_5 + b_1,   a_2 = b_2,
    b_1 = 2P_{11} + P_{21} + 2P_{31} + P_{22} + 3P_{111} + 3P_{211},
    b_2 = −P_{11} + 2P_{21} − P_{31} + 2P_{22} + P_{111} + P_{211}

(brackets: (1,1) → 2δ_1 − δ_2, (2,1) → δ_1 + 2δ_2, (3,1) → 2δ_1 − δ_2, (2,2) → δ_1 + 2δ_2,
(1,1,1) and (2,1,1) → 3δ_1 + δ_2). At T = 0, P_5 = 1 and f = f_MF; `molofsky_coeffs` recovers
(a_1, c_1, a_2, c_2) = (5, −5, 0, 10) and, at every T, the quotient G_i/h is a + c u(1−u) to
rounding (column `resid`).

For a direction (δ_1, δ_2) = ρ(cos θ, sin θ) write a = δ_1 a_1 + δ_2 a_2 and c = δ_1 c_1 + δ_2 c_2.
Then

    f'(0) = a,       f'(1/2) = −(a + c/4)/2,       ∫_0^{1/2} f du = (6a + c)/192,

using h'(1/2) = −1/2, ∫_0^{1/2} h = 1/32 and ∫_0^{1/2} h·u(1−u) = 1/192. So the plane around the
voter point is cut by three rays through it (two opposite rays per line; the relevant ones
are those with c > 0):

- θ_a: a = 0, the absorbing state loses linear stability (ergodic for a > 0, i.e. θ < θ_a);
- θ_front: 6a + c = 0, the front between u = 0 and u = 1/2 is stationary; the active state
  invades for 6a + c > 0;
- θ_half: a + c/4 = 0, the active state loses stability (clustering for θ > θ_half).

Both states are stable for θ_a < θ < θ_half, the closure's bistable sector, and the front
reverses inside it (at a = 0, 6a + c = c > 0; at a + c/4 = 0, 6a + c = 2a < 0). Its width is set
by the five-class weight: Δθ ≈ c(θ_a)/(4|∇a|) ≈ 2.6 P_5/a_1 rad (1.92° at T = 64 against 1.94°
measured). Because c ∝ P_5 carries both signs as θ turns, one of the two rays of the line a = 0
always has c > 0: **a rule with two free parameters always has a bistable direction at finite T.**
The sign question of the manuscript decides the matter only for a one-parameter family (as for
D&N at fixed α, or for a single transect of Molofsky's plane). The sector logic was checked by
evaluating f on the polynomials at T = 64: θ = 100° ergodic; 105.7° and 106.5° bistable with the
active state invading; 106.9° and 107.0° bistable with the absorbing state invading; 107.6° and
112° clustering.

## 4. Numbers

10^6 trials (`patterns5.txt`); angles in degrees from the +δ_1 axis; b_1, b_2 as above; `th_3cl` is
the ray a = 0 from three-class patterns only (the analogue of `kappa_th3`, the finite-T form of
Cox & Perkins' Θ_3 = 0). Mean field is the T = 0 row.

| T | P(≥3 cl) | P_5/P(≥3) | b_1 | b_2 | θ_a | θ_front | θ_half | width | th_3cl |
|---|---|---|---|---|---|---|---|---|---|
| 0 (MF) | 1 | 1 | 0 | 0 | 90.00 | 111.80 | 123.69 | 33.69 | – |
| 1 | 0.991 | 0.470 | 1.432 | 0.474 | 97.18 | 110.34 | 117.27 | 20.09 | 107.26 |
| 2 | 0.970 | 0.311 | 1.718 | 0.552 | 99.71 | 109.51 | 114.61 | 14.91 | 104.73 |
| 4 | 0.925 | 0.190 | 1.797 | 0.560 | 101.82 | 108.63 | 112.13 | 10.30 | 104.06 |
| 8 | 0.852 | 0.112 | 1.698 | 0.516 | 103.35 | 107.87 | 110.17 | 6.82 | 104.08 |
| 16 | 0.763 | 0.066 | 1.506 | 0.451 | 104.39 | 107.33 | 108.82 | 4.43 | 104.40 |
| 32 | 0.671 | 0.040 | 1.294 | 0.384 | 105.02 | 106.94 | 107.92 | 2.90 | 104.66 |
| 64 | 0.585 | 0.025 | 1.097 | 0.323 | 105.40 | 106.69 | 107.33 | 1.94 | 104.89 |
| 128 | 0.508 | 0.017 | 0.928 | 0.272 | 105.67 | 106.55 | 106.99 | 1.33 | 105.14 |
| 256 | 0.440 | 0.011 | 0.786 | 0.230 | 105.82 | 106.43 | 106.73 | 0.91 | 105.28 |
| 512 | 0.382 | 0.0074 | 0.669 | 0.195 | 105.96 | 106.38 | 106.59 | 0.63 | 105.45 |
| 1024 | 0.333 | 0.0052 | 0.573 | 0.167 | 106.03 | 106.33 | 106.48 | 0.45 | 105.57 |
| 2048 | 0.292 | 0.0036 | 0.494 | 0.144 | 106.05 | 106.26 | 106.36 | 0.32 | 105.61 |
| 4096 | 0.256 | 0.0026 | 0.428 | 0.124 | 106.09 | 106.25 | 106.33 | 0.24 | 105.70 |
| 8192 | 0.226 | 0.0019 | 0.373 | 0.109 | 106.14 | 106.25 | 106.31 | 0.17 | 105.79 |
| 16384 | 0.200 | 0.0013 | 0.328 | 0.095 | 106.08 | 106.16 | 106.20 | 0.12 | 105.75 |
| 32768 | 0.178 | 0.0011 | 0.289 | 0.083 | 106.05 | 106.12 | 106.15 | 0.10 | 105.75 |
| 65536 | 0.158 | 0.0008 | 0.256 | 0.074 | 106.07 | 106.12 | 106.15 | 0.08 | 105.80 |

Bootstrap errors (300 multinomial resamples): the angles ±0.01° at T = 1 rising to ±0.13° at
T = 65536; the width ±0.016° at T = 1 falling to ±0.007° at T = 65536. The 10^7 run
(`vpt_run_molofsky_1e7.txt`) agrees everywhere within these errors, with errors three times
smaller (at T = 65536: θ_front = 106.14 ± 0.04°, width 0.07 ± 0.002°).

Slopes dp_2/dp_1 at T = 65536 (10^7 run): −3.47 (θ_a), −3.46 (θ_front), −3.45 (θ_half); from
three-class patterns alone −3.71 and still drifting toward the full rays. The line through the
voter point with slope −3.46 meets p_2 = 1 at p_1 = 0.026 and p_2 = 0 at p_1 = 0.316.

Width of the sector in the (p_1, p_2) plane at the closure's own matching distance ρ = 1/T
(perpendicular width ρ·Δθ): 0.045 (T = 4), 0.015 (8), 0.005 (16), 0.0016 (32), 0.0005 (64),
0.0002 (128).

## 5. Comparison with Molofsky et al.

- **Mean field.** The T = 0 row is their Fig. 3 exactly: the multiple-outcome wedge between 90°
  and 123.7° (p_1 < 1/5, 15 p_1 + 10 p_2 > 7). The closure's sector lies inside that wedge at
  every T (θ_a ≥ 90°, θ_half ≤ 123.7°), on the side δ_1 < 0, δ_2 > 0: small p_1, large p_2,
  where they looked. It shrinks from the full wedge to a line at ≈ 106°.
- **Their spatial boundary.** They report the clustering/ergodic boundary running "from about
  (0.024, 1) at the top edge, down through the linear voter point (0.2, 0.4)", and the dyad rule
  3 p_1 + p_2 > 1 for spread, "in excellent agreement" with it. Through the voter point: the chord
  to (0.024, 1) has direction 106.4°, the dyad line 108.4°. The closure's T → ∞ transition line is
  at 106.1° (slope −3.46 against −3.41 and −3): it reproduces the direction of their empirical
  boundary to within a degree with no fitted parameter. Extended linearly, it meets the top edge at
  p_1 = 0.026 (theirs: 0.024) and the bottom edge at p_1 = 0.316 (theirs: clustering at 0.27,
  phase separation at 0.31, ergodic at 0.35; dyad rule 1/3). The extension is a linear
  extrapolation to O(1) distances where the closure is not controlled; take it as suggestive. The
  tilt itself (a_2 = b_2 > 0: p_2 helps the rare species, which in mean field it cannot) is the
  dyad effect: a rare 1 has a neighbouring 1 with a probability of order one because the two
  lineages have not separated.
- **Their transects.** p_2 = 0.4 with p_1 ∈ [0, 0.2] is the direction 180°: clustering at every T
  (extinction time growing toward the voter point, their Fig. 8 top). p_1 = 0.2 with
  p_2 ∈ [0.5, 1] is 90°: a = b_2 δ_2 > 0, ergodic at every T ≥ 1 (their Fig. 8 bottom, with
  clustering near p_2 = 0.5). Neither transect crosses the sector, which starts at 97° already at
  T = 1.
- **Their search for multiple outcomes.** They started "barely ergodic systems for large sampling
  values of p_2" from random configurations of density 0.01 and saw "small isolated configurations
  of species 1, such as the dyads, nucleate linearly growing patches of the symmetric
  equilibrium". Those parameters are at distance ρ ≈ 0.1–0.6 from the voter point, where the
  matching time is T ≈ 2–10 and the closure is at best qualitative. Where it does apply, it
  predicts exactly what they saw: on the active-invades side of the sector (Case 4A) the absorbing
  state is linearly stable but only metastable against seeds below a critical droplet of radius
  R_c ∼ (D T)^{1/2} × O(1), D = 1/5; at ρ = O(1) that is a site or two, and their own monad/dyad
  rules (a single 1 dies on average for 5 p_1 < 1, a dyad grows for 3 p_1 + p_2 > 1) describe a
  nucleus of two sites. A density-0.01 initial condition on a 512² lattice holds ≈ 50 dyads and
  nucleates at once. Their result therefore rules out a *robust* multiple-outcome region at O(1)
  distance from the voter point, not a window of the size the closure predicts near it: a width
  Δp ≈ 5×10^{-3} at ρ = 1/16 and 5×10^{-4} at ρ = 1/64, invisible to any parameter grid and
  undetectable by seeding at density 0.01.

## 6. Comparison with D&N (Moore neighbourhood, `patterns.txt`)

| T | Moore P(≥5 cl)/P(≥3) | five-cell P_5/P(≥3) | ratio | D&N (κ_a − κ_front)/κ_front | five-cell Δs/\|s\| |
|---|---|---|---|---|---|
| 16 | 0.70 | 0.066 | 10 | 1.13 | 0.30 |
| 64 | 0.43 | 0.025 | 17 | 0.54 | 0.13 |
| 256 | 0.25 | 0.011 | 23 | 0.28 | 0.059 |
| 1024 | 0.15 | 0.0052 | 29 | 0.16 | 0.029 |
| 4096 | 0.090 | 0.0026 | 34 | 0.098 | 0.015 |
| 16384 | 0.056 | 0.0013 | 42 | 0.060 | 0.0080 |
| 65536 | 0.036 | 0.0008 | 42 | 0.039 | 0.0051 |

(s = dp_2/dp_1 of the rays; Δs/|s| is the relative width of the sector in the slope, the analogue
of the relative width of [κ_front, κ_a] in κ = q/(1 − p). Both are parametrisation-dependent;
the five-lineage fractions are not.) With five walkers starting on a plus, all five must stay
distinct to give a five-class pattern; with nine on a 3×3 block there are many ways to keep five
classes. The five-lineage weight is 10–40 times smaller for the five-cell rule, and its window is
4–8 times narrower in relative terms. Both close in the same way, as a power of 1/ln T (the width
falls by a factor ≈ 0.5–0.65 per factor 4 in T at large T).

## 7. Answer

**The closure predicts a bistable window for the Molofsky rule near its voter point.** It lies in
the direction of small p_1 and large p_2, inside the mean-field multiple-outcome wedge, around a
transition line through the voter point of slope dp_2/dp_1 ≈ −3.5 that coincides with their
empirical clustering/ergodic boundary (dyad line, slope −3; their top-edge point, slope −3.4).
On the ergodic side of the line, within the sector, the absorbing state is linearly stable and
the active phase wins by invasion (Case 4A); on the other side the active state is stable and
the absorbing phase invades. The window is narrow: 4.4° wide at T = 16, 1.9° at T = 64, 0.45° at
T = 1024, 0.08° at T = 65536, four to eight times narrower (relative) than D&N's at the same T,
because the five-lineage weight is 10–40 times smaller. The closure does not, and for any
two-parameter rule cannot, predict the absence of a window: the sign of the five-lineage term is
favourable in some direction at every finite T. What it predicts for Molofsky's rule is where the
window is and how thin, and it is consistent with their negative search, which probed O(1)
distances with one- and two-site seeds. It does not explain the difference between the two rules
by a sign; it explains it by a weight.

**Caveats.** (i) The closure is heuristic: first order in g, voter statistics at scale T, no
limit theorem in d = 2. (ii) The matching T ∼ 1/ε is undetermined up to an O(1) factor, and here
ε is a direction-dependent magnitude ρ = (δ_1² + δ_2²)^{1/2}; the angular widths are quoted
against T, the conversion to Δp uses ρ = 1/T. (iii) The constants are for the synchronous model
(lazy von Neumann dual); the continuous-time version would differ at the 5–10% level seen for
D&N. (iv) The agreement of the asymptotic line with their boundary at the edges of the square is
an extrapolation to ε = O(1). (v) Nothing in Secs. 1–8 is a simulation of the Molofsky rule (Sec. 9 reports one, run
afterwards); the test that would settle it is the manuscript's own: planar fronts and critical droplets for their rule
at ρ ≈ 0.02–0.05 along the line p_2 − 2/5 = −3.46 (p_1 − 1/5), resolving Δp ≈ 10^{-3} across
it, where the closure says both phases are stable and the front reverses.

## 8. For the manuscript (superseded: the simulation of Sec. 9 and this closure were written into `paper/main.tex`, Appendix D, Table IV and Fig. 10, on 10-30; see `paper/CHANGES_2026-10-30.md`, Addendum 2)

*Update (Sec. 9).* The simulation has since been run. The last sentence of the draft below
("Their search ... could not have resolved it; planar fronts and critical droplets ... would") is
outdated, and the first ("predicts a bistable window there as well") needs the qualification that
only the absorbing-state half of the window was found. A merged draft that folds in the
simulation is in `campaign/molofsky_sim/NOTES.md`, Sec. 7. The two other Discussion sentences
flagged at the end of this section still need the change described there (sign → weight).

Replacement for "The first test should be the rule of Molofsky \emph{et al.}~\cite{molofsky1999},
whose five-cell neighborhood admits five-lineage patterns: the closure predicts whether a bistable
window exists there, where the simulations found none." in the Discussion:

> We have applied the closure to the five-cell rule of Molofsky \emph{et al.}~\cite{molofsky1999},
> whose only five-lineage pattern is the mean-field one. It predicts a bistable window there as
> well, in the direction of small $p_1$ and large $p_2$ in which they searched, around a transition
> line through the voter point of slope $dp_2/dp_1\approx-3.5$ that coincides with their empirical
> boundary $3p_1+p_2=1$; but the window is four to eight times narrower than for D\&N at the same
> $T$ ($4^\circ$ wide at $T=16$ and $0.5^\circ$ at $T=1024$), because five lineages stay distinct ten
> to forty times more rarely than nine. Their search, at distances of order one from the voter
> point and with seeds of one or two sites, could not have resolved it; planar fronts and critical
> droplets within $10^{-3}$ of the predicted line would.

Two other sentences now contradict the result and should change with it. In the paragraph "Why is
the bistability not washed out?": "the difference lies in the rule, and the closure of
Sec.~\ref{sec:theory} says where: in the sign of the five-lineage terms, which the first-derivative
closure of Ref.~\cite{vazquez2008} discards" → the sign is favourable for their rule too; the
closure locates the difference in the *weight* of the five-lineage patterns (an order of magnitude
smaller for a five-cell neighbourhood) and in the size of the critical nucleus at the parameters
they studied (a pair of sites). And in the last paragraph: "What matters near the voter point is
the sign of the contribution of five or more coalescing lineages to the reaction term, which the
rule's nonlinearity in $j$ controls; for the D\&N rule it opens a bistable window" → for a
one-parameter family the sign decides; in a two-parameter plane there is always a bistable
direction, and the rule controls its width and orientation. The lineage criterion is a criterion
for the width of the window, not for its existence.

## 9. Direct simulation (added after the side session of 10-30; `campaign/molofsky_sim/`)

The test named in Sec. 7(v) was run in a side session: the rule was added to a copy of the
simulator (`campaign/molofsky_sim/dnsim.c`, `--molofsky p1,p2`, diff in `dnsim_v14.patch`; the
root `dnsim` is untouched). Checked here again: `./dnsim --selftest` passes (von Neumann counts
against a naive count, the rule's table and its Z2 symmetry, flip frequencies, box stepping), and
on the Moore path a voter-mixed strip run (Lx 256 × Ly 256, T 500, seed 5) gives the same
time-series lines and the same binary profile as the root binary (only the version header
differs). Runs: planar fronts at p_2 = 0.45 and 0.50 (Lx 1024 and 2048 × Ly 4096, T 6e4, steps of
2.5e-4 near the crossing, 3–8 seeds), 2 × 2-seed spreading (2e4 replicas to 1e4), active disks
(≤ 400 replicas per radius to 2e4), bulk stability from random 1/2 (L 256–2048, 6e4–1e6 steps).
Tables, figure (`fronts.png`) and the full verdict: `campaign/molofsky_sim/NOTES.md`.

Against Secs. 4 and 7 (angles from the +δ_1 axis at the voter point; p_1 = 1/5 + (p_2 − 2/5) cot θ):

| | closure (Sec. 4) | simulation |
|---|---|---|
| transition line through the voter point | θ_front = 107.3° (T = 16) → 106.1° (T → ∞) | v = 0 at p_1* = 0.18535(10) for p_2 = 0.45 and 0.17057(10) for p_2 = 0.50: θ = 106.33°, 106.40°, slope dp_2/dp_1 = −3.41 ± 0.02 |
| line extended to the edges of the square | (0.026, 1) and (0.316, 0) | (0.0235, 1) and (0.318, 0); Molofsky: (0.024, 1); clustering at (0.27, 0), ergodic at (0.35, 0) |
| absorbing state linearly stable up to θ_a | 103.4° (T = 8), 104.4° (T = 16), 105.0° (T = 32) | 104.0–104.8°: the mean population of a 2 × 2 seed decays up to 1.5–2e-3 above p_1* |
| critical droplet on the invading side (Case 4A) | R_c = D/v, D = 1/5 | R_c v = 0.16–0.19 at five points (R_c = 20, 10, 27, 12 at v = 0.0090, 0.0180, 0.0059, 0.0155) |
| stable active phase with a retreating front (Case 4B), θ_front − θ_half | 1.5° (T = 16), 0.6° (T = 64), 0.15° (T = 1024) | none within 2.5e-4 of p_1* (0.3° at p_2 = 0.45, 0.15° at 0.50): just below p_1* the active phase coarsens away, its lifetime growing with L (4.8e4, 1.3e5, > 6e4 steps at L = 256, 512, 2048) |
| width of the bistable band, perpendicular to the line, at ρ ≈ 0.05 | 3.9e-3 (T = 16), of which 1.3e-3 on the clustering side | 1.7e-3, all on the ergodic side |

Reading.

1. **The line is right.** The measured line is 0.2–0.3° from the closure's asymptotic one and on
   the chord to Molofsky's top-edge point; extended, it reproduces their boundary at both edges of
   the square. It coincides with the closure's θ_front at T ≈ 250–1000, not at T ∼ 1/ρ ≈ 10–20,
   where the closure's front ray is 1–1.5° too high.
2. **The absorbing-state half of the window (4A) is right**, with the predicted critical nucleus
   D/v and about the width the closure gives at T ≈ 10–40: 2e-3 in p_1 at ρ = 0.05.
3. **The active-state half (4B) is absent at the resolution reached.** The active phase loses
   stability where the front stops, within 2.5e-4 in p_1, where the closure puts θ_half 1.4e-3
   (T = 16) to 1.4e-4 (T = 1024) further. The quintic form f = h(u)[a + c u(1 − u)] fixes the front
   reversal at two thirds of the way from θ_a to θ_half (Sec. 3); the rule puts it at the window's
   edge. No single T fits all three rays: θ_a wants T ≈ 10–40, θ_front T ≈ 250–1000, and the bound
   on θ_half needs T ≳ 300–1000 if the closure's form held at all.

**Revised answer.** For Molofsky's rule the bistable window is the metastability of the absorbing
state alone, on the ergodic side of a sharp transition line through the voter point in the
direction the closure predicts. The closure locates the line and that metastability, and it
explains why their search at O(1) distance saw neither (the critical nucleus is D/v, a few cells
once v ≈ 0.05); it overstates the stability of the active phase beyond the front reversal. The
difference from D&N, whose active phase stays metastable far past the reversal, is still a weight
of the five-lineage patterns, not a sign; but what the weight controls is the extent of the
active-state branch, not the existence of the absorbing-state branch. Sec. 7's statement "on the
other side the active state is stable and the absorbing phase invades" is not borne out at the
2.5e-4 level. A finer scan (1e-4 steps, Lx 4096) would tighten the bound by a factor of 2–3.
