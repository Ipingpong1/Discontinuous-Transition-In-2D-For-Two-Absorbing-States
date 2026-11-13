# Pair approximation for noisy D&N (10-16)

    python3 pa.py selftest
    python3 pa.py point --p 0 --q 0.02 --alpha 0.5      # fixed points, stability, absorbing rate
    python3 pa.py compare --csv compare.csv             # vs simulated q* (table in compare.txt)
    python3 pa.py curve --alpha 0.5 --csv curve_a05.csv # q_front along p
    python3 plot_compare.py                             # -> pa_vs_sim.png

## Model in PA form
D&N is self-complementary: births at k ∈ {3,6,7,8} and deaths at k ∈ {0,1,2,5} are both
j ∈ {3,6,7,8} in terms of j = # disagreeing Moore neighbours. So the whole dnsim model is
one flip table

    w(j) = (1−p)·α·[j ∈ {3,6,7,8}] + p·j/8 + q·[j = 2]

and sits in the generalized-voter family (Vázquez & López, PRE 2008).

**Closure.** Bethe lattice with z = 8. Given a site's state, its neighbours are i.i.d.
and disagree with probability θ_s = ρ/(2x_s), where ρ is the fraction of disagreeing
bonds. Variables (m, ρ). Two versions:
- **Continuous time:** rate w(j), one time unit = one dnsim step.
- **Synchronous map** (the actual dnsim update): a bond changes state iff exactly one end
  flips. The two ends' flips are independent given the bond.

Site MF is the same with ρ slaved to 2x(1−x). Equations are in the `pa.py` docstring.

**Selftest:**
- The voter limit conserves m, gives ρ* = (z−2)/(2(z−1))·(1−m²), and leaves the absorbing
  state marginal (θ₀ → 6/7).
- The ρ equation matches the (2/z)Σ w(j)(z−2j) form.
- The map → the ODE for small w.
- D&N is self-complementary.

## Results
1. **The absorbing states are always linearly stable for p < 1.** Analytically, a dilute
   minority grows at rate z·w(1)·θ − F₀(θ). With w(1) = p/8 this is
   −[(1−p)α·P_rule(θ) + q·P₂(θ)] < 0. Numerically the rate is −0.5 = −(1−p)α at p = 0.
2. **The disordered active state is also stable at every q ≥ 0, every p < 1, every α.**
   Checked at α = 0.125–1 and p = 0–0.99. The symmetric fixed point has ρ* ≈ 0.43–0.50, is
   stable in both ρ and m, and exists even for pure D&N at q = 0. An unstable ordered saddle
   (m ≈ 0.5, ρ ≈ 0.27 at p = 0) separates the two basins. There is no ordered-active
   attractor, so the PA has no Ising transition. Site MF agrees: m = 0 is stable at q = 0.
3. **So PA/MF are bistable everywhere and have no q\*.** Homogeneous stability never
   selects a phase. In 2D the phase is selected by which phase invades the other, which is
   what the simulations measure (front velocity, droplets, walls). A Bethe lattice has no
   curvature-driven coarsening. That is why it keeps D&N active at q = 0, where 2D
   coarsens.
4. **The front criterion is only a rough guide.** For a bistable reaction–diffusion front
   with constant diffusion, the absorbing phase invades iff I = ∫₀¹ ṁ(m, ρ_s(m)) dm > 0,
   with ρ on its stable slow branch. ρ is the fast variable: λ_ρ ≈ −0.7…−0.9 against
   λ_m ≈ −0.2. Taking q_front as the zero of I(q) (`compare.txt`, `pa_vs_sim.png`):
   - **Near the voter point it gets the scaling right but not the slope.** For p ≥ 0.9,
     q* ∝ (1−p) with the measured ratio q*_sim/q_front = 2.0 for the synchronous map (5.5 in
     continuous time). The linearity alone is automatic: near p = 1 everything is first order
     in (1−p) and q. The constant ratio says the PA shape of the perturbation is right while
     the amplitude is off by 2×.
   - **Away from the voter point it fails.** It gives a hump in q*(p), but peaking at
     p ≈ 0.6 rather than the measured 0.25–0.5, and 4–10× too low. It gives q_front = 0
     (active wins even with no noise) for p ≈ 0.05–0.2 at α = 0.5. At p = 0 it misses
     q* ∝ α: the ratio is 3–90 and not monotone.
   - **ρ\* is about 2× the measured ρ_bond** (0.43–0.50 against 0.23–0.29). The Bethe
     closure ignores the triangles of the Moore neighbourhood.

## What this means for the paper
- **Usable as stated:**
  - The PA agrees that the transition is discontinuous at the mean-field/Bethe level: the
    absorbing and active states coexist, with no continuous (DP/GV/Ising) route. This
    matches the mean-field picture of Rodrigues et al. 2015.
  - The PA also makes clear why 2D needs front-based numerics: a q* exists only through
    interface motion.
- **Not usable:** the PA as a quantitative predictor of q*(p, α).
- **Next steps, if a theory section needs more:**
  1. **A perturbative treatment around the voter point.** Write w = j/8 + (1−p)·δw_rule +
     q·δw_2 and compute the drift of the Al Hammal et al. Langevin coefficients to first
     order. This is the controlled regime, and the PA already shows q* ∝ (1−p) with an O(1)
     amplitude error. The slope 0.125 is the number to explain.
  2. **A 2×2-plaquette cluster approximation**, which includes Moore triangles. It should fix
     ρ* and may move the hump.
  3. **An interface-level argument at p = 0:** the sign of the front velocity from wall
     roughening versus curvature, to get q* ∝ α.

# Voter-point perturbation theory (10-16)

    python3 vpt.py selftest
    python3 vpt.py run --trials 1000000 --seed 2 --T 0,1,2,...,65536 --boot 200   # builds cwalk; ~1 min on 12 cores
    python3 vpt.py run --load patterns.txt                                        # re-analyse (output in vpt_run.txt)
    python3 plot_vpt.py                                                           # -> vpt_kappa.png

## Idea
Near p = 1 write w(j) = j/8 + ε·g(j), with ε = 1 − p, κ = q/(1−p) and

    g(j) = α·R(j) − j/8 + κ·[j = 2],   R(j) = [j ∈ {3,6,7,8}]

This is a voter-model perturbation in the sense of Cox, Durrett & Perkins (Astérisque 349,
2013, arXiv:1103.1676; their eq. 1.16 is this f). Their theorems are for d ≥ 3 only; the 2D
use here is heuristic. The 2D analogues use asymptotic coalescence probabilities (Cox & Perkins
2008, 2023). See LITERATURE.md. To first order in ε, the local density obeys a reaction–diffusion equation
∂u = (diffusion) + ε·f(u), with

    f(u) = E_u[(1 − ξ₀)·g(# live nbrs) − ξ₀·g(# dead nbrs)].

The expectation is over the voter state. By duality, the colours of the origin and its 8
Moore neighbours are those of coalescing random walks (`cwalk.c`: the discrete-time Moore
walk, dual to the synchronous voter update); each surviving lineage is live with probability
u. Exact properties, all checked by `vpt.py selftest`:
- Any g linear in j gives f ≡ 0, so the voter part and the −j/8 term drop out.
- Patterns with 1 or 2 lineage classes cancel exactly by the up–down symmetry.
- A 3-class pattern, whose non-origin classes hold n₁ and n₂ neighbours, contributes
  u(1−u)(1−2u)·[g(n₁) + g(n₂) − g(n₁+n₂)].

So **f = α·f_R + κ·f₂ exactly, and every threshold is κ\* = α × (a pure number)**. That
gives the prediction **q\* = c·α·(1−p)** near the voter point.

Two thresholds from f:
- **κ_a**, where f′(0) = 0: the absorbing state loses stability.
- **κ_front**, where ∫₀^½ f du = 0: the front between u = 0 and the active state u = ½ stops.

For κ_front < κ < κ_a both u = 0 and u = ½ are stable (f′(½) < 0 throughout). That is a
first-order coexistence window. In 2D all lineages eventually coalesce, so the pattern
weights depend on the walk time T ~ 1/ε. Once 3-class patterns dominate (T → ∞), the two
thresholds merge; the 4+-class terms that open the window fade only like powers of 1/ln T.

## Results (1M walks, `vpt_run.txt`, `vpt_kappa.png`)

| T | κ_a/α | κ_front/α | P(3 classes) / P(≥3) |
|---|---|---|---|
| 64 | 0.355 | 0.231 | 0.21 |
| 1024 | 0.289 | 0.249 | 0.48 |
| 65536 | 0.265 | 0.255 | 0.72 |

`kappa_th3` (3-class patterns only; the finite-T version of Cox & Perkins' Θ₃ = 0, see
LITERATURE.md) is 0.308, 0.279 and 0.269 at the same T, falling toward the same limit.

Bootstrap errors are about ±0.002. The two thresholds close in from opposite sides, so
**κ\*/α → 0.26 ± 0.015 as T → ∞**. Extrapolations in 1/ln T scatter between 0.21 and 0.27
and aren't reliable; the bracket is the honest error.

**Comparison with the simulations:**
- **α = 0.5, p = 0.9–0.995** (E4/E6 fronts): q*/((1−p)α) = 0.238, 0.248, 0.250, 0.254,
  0.248. These lie inside the theoretical window [κ_front, κ_a] for every T ≳ 100 and within
  4% of the T → ∞ value. The near-constancy of q*/(1−p) is the predicted linear law. The
  slight rise with p (0.238 at p = 0.9) fits κ_front(T) increasing with T ~ 1/ε.
- **α test at p = 0.99** (new, local, `alpha_check/`, `alpha_check_fronts.txt`; strip
  fronts, Lx 1024 × Ly 4096, T 6e4, 1 seed per q). The prediction q*(α) = 2α·q*(0.5) =
  2α·0.00127 had no free parameter:
  - α = 0.25: v changes sign between q = 0.000572 (−0.0019) and 0.000635 (+0.0012; that
    strip ordered at 5.1e4), so q* ≈ 0.00061 and κ*/α ≈ 0.244. Predicted 0.000635.
  - α = 1: v = −0.0003 at 0.00254 and +0.0056 at 0.00279, so q* ≈ 0.00255 and
    κ*/α ≈ 0.255. Predicted 0.00254.
  - So **q\* ∝ α holds within ~4% over a factor 4 in α.** (The PA predicted nothing here.)
- **E9_alpha on Torch (10-17; campaign README, Phase 5).** This is the same test with E4
  geometry (Lx 2048 × Ly 4096, T 1e5), 4 seeds, 5 q per (p, α), and p ∈ {0.98, 0.99}. Results:
  - c = q*/(α(1−p)) = 0.255, 0.252, 0.254 (p 0.98) and 0.248, 0.258, 0.255 (p 0.99) for
    α = 0.25, 0.5, 1.
  - The exponent in q* ∝ α^β is β = 1.00 ± 0.02.
  - Systematics are about ±2%: the α 0.5 rerun vs E4, and 3- vs 5-point fits.
  - The mean c = 0.254 sits on κ_front(T = 65536) = 0.2553.

## What this gives the paper
- **A controlled theory of the phase boundary near the voter point.** q* = c·α·(1−p) with
  c = 0.26 ± 0.015 from coalescing-walk statistics, with no fit parameters. It matches
  E4/E6 and the new α runs to a few %.
- **A mechanism for the first-order character fading toward p = 1.** In the reaction–
  diffusion limit the coexistence window comes only from 4+-lineage terms, which are
  suppressed by powers of 1/ln(1/ε). So the theory predicts a first-order line that becomes
  weak only logarithmically and ends at the voter point (tricritical at p = 1). This fits
  the slowly shrinking jump (E4/E6) and the growing voter-like transient (E8b,
  t_x ∝ (1−p)^−1.3), and it explains why p > 0.995 can't be settled numerically.
- **Caveats:**
  - The rigorous version of this limit is proved for d ≥ 3. In d = 2 the T ~ 1/ε
    identification carries log corrections (cf. Cox & Perkins, AAP 2008; Cox & Perkins,
    arXiv:2310.00812, which covers only the 2D coexistence side).
  - Here f is the mean-field reaction term of the rescaled equation. Whether 2D noise
    destroys the window at small ε is exactly the GV-versus-first-order question, which
    E7/E8 test numerically.
  - The α runs are 1 seed each on a smaller lattice than E4.
- **Literature:** checked 10-16; see LITERATURE.md. Molofsky et al. 1999 found *no*
  bistable regime in their 2D (von Neumann + self) model, because tiny seeds nucleate.
  Cox–Durrett–Perkins leave the bistable cases (3, 4B) open even in d ≥ 3. No prior 2D
  first-order result was found.

## Lazy (asynchronous) dual, 10-18
`cwalk.c` reads `CWALK_R` (default 1; at r = 1 the output is bit-identical single-threaded). Each
lineage jumps with probability r per step. This is the dual of dnsim with the flip table scaled by
r (campaign README, Phase 6).

| r | κ_a/α at T → ∞ | κ_front/α | κ_th3/α | P(≥ 3 classes) at physical T 65536 |
|---|---|---|---|---|
| 1 | 0.2652 | 0.2553 | 0.2692 | 0.53 |
| 0.5 | 0.223 ± 0.003 (physical T 65536) | 0.221 | 0.226 | 0.29 |
| 0.25 | 0.208 ± 0.004 | 0.206 | 0.209 | 0.23 |
| 0.125 | 0.201 ± 0.004 | 0.200 | 0.205 | 0.21 (T 524288 steps) |

- The prediction is q*(r) ≈ 0.8 q*(1) for r ≤ 0.25. That is the continuous-time
  (Cox–Perkins-type) limit.
- A local pilot instead found q*(0.25) ≈ q*(1). E10 (campaign Phase 6 results, 10-18) superseded it:
  - p 0.98: q*(r)/q*(1) = 0.86–0.87 at r = 0.5, 0.25, 0.125, so c_sim = 0.218 against the theory's
    0.208 (r 0.25) and 0.201 (r 0.125).
  - p 0.99, r 0.25: 0.93–0.95, so c_sim = 0.236–0.240.
- **Status.** The theory gets the sign of the shift and most of its size at p 0.98 (−20% predicted,
  −14% measured). Its error is +3% at r = 1 and −4 to −8% at r ≤ 0.25. The p 0.99 lazy value is
  further off and runs the wrong way for an ε → 0 limit, so check it for convergence before using it.
- r 0.5 (added 10-18): theory 0.223 against the simulation's 0.218. The theory and the simulation
  agree on the shift from r 1 to 0.5. Below 0.5 the simulation is flat while the theory falls
  another 10% (~4σ).
- p 0.99 resolved (Phase 7, E11_lazy_front_p99 to physical 3e5): the high value was a ~1e5
  retreating transient. Late-window crossing ≈ 0.85–0.88, so c_sim(r ≤ 0.5) ≈ 0.22 at both p
  (theory 0.201–0.223).
- For the paper: present c as the leading-order prediction, with an O(5–10%) error budget that
  covers the updating scheme. Do not present the 3% match at r = 1 as exact.
