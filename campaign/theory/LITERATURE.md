# Literature check (10-16)

"Verified" means the abstract, publisher page or full text was read. Unverified items are
marked.

## Same model family (symmetric flips depending on local counts, both states absorbing)
- **Molofsky, Durrett, Dushoff, Griffeath & Levin**, Theor. Pop. Biol. 55, 270 (1999),
  doi:10.1006/tpbi.1998.1404. Full text read.
  - Model: synchronous PCA on the 5-cell neighbourhood (von Neumann + self) with parameters
    p₁, p₂. The voter point is (0.2, 0.4).
  - 2D phases: clustering, ergodic (linear-rate invasion fronts), a tiny phase-separation
    region, and a quasi-periodic stripe phase. The clustering/ergodic boundary runs from
    (0.024, 1) through the voter point.
  - **They found no bistable ("multiple outcome") regime in 2D:** "small isolated
    configurations … have a propensity to nucleate linearly growing patches … there is no
    multiple outcome regime of the spatial model." The order of the transition is not
    discussed.
  - → D&N with the Moore neighbourhood is a counterexample: E5/E8a show seeds below R_c die.
- **Cox & Durrett**, "Nonlinear voter models", in *Random Walks, Brownian Motion and
  Interacting Particle Systems* (Birkhäuser 1991), pp. 189–202. Citation verified, paper not
  read. It introduces the 5-site discrete family and calls the Levin line 3p₁+p₂ = 1 the
  empirical boundary.
- Liggett, Ann. Probab. 22, 764 (1994): coexistence in threshold voter models. Durrett &
  Steif, Ann. Probab. 21, 232 (1993): fixation for thresholds > n/2.
- **Schweitzer & Behera**, EPJ B 67, 301 (2009), arXiv:cond-mat/0307742, "Nonlinear voter
  models: the transition from invasion to coexistence". Pair approximation, m = 4. (Earlier
  mis-cited here as Cox & Durrett 2016.)

## Voter-model perturbation theory (the method in vpt.py)
- **Cox, Durrett & Perkins**, Astérisque 349 (2013), arXiv:1103.1676.
  - f(u) = ⟨(1−ξ₀)h₁ − ξ₀h₀⟩_u via coalescing walks (eqs. 1.16, 1.26). Proved for **d ≥ 3
    only**.
  - For the continuous-time Molofsky model:
    f₁ = u(1−u)(1−2u)[b₁(1−u+u²) + b₂u(1−u)], with four cases from the signs of f′(0) and
    f′(½).
  - **Cases 3 and 4B (bistable, which is our regime): "our method yields no result."** The
    words first order and tricritical do not appear.
  - In the bistable case the sign of the wave speed is the sign of ∫f (eq. 1.28).
- Durrett, EJP 19 (2014), arXiv:1406.5876: spatial games with weak selection, d ≥ 3, cubic
  reaction terms, bistable cases decided by wave speed.
- Cox & Durrett, SPA 126, 2388 (2016), arXiv:1511.04713: evolutionary games on the torus,
  d ≥ 3.
- **2D work.**
  - Cox & Perkins, AAP 18, 747 (2008), arXiv:0804.0158: 2D Lotka–Volterra renormalization.
  - Cox, Merle & Perkins, EJP 15, 1190 (2010): 2D coexistence region with a thorn shape,
    tip at the voter point.
  - **Cox & Perkins, arXiv:2310.00812 (2023; EJP 30, 2025): complete convergence for 2D voter
    perturbations including the q-voter.** Read §1–5 on 10-16; see the section below.
- **Not found:** any computation of a first-order boundary near the 2D voter point.

## Cox & Perkins 2023 in detail (read 10-16; lit/2310.00812_cox_perkins2023.pdf)
- **Setting.** Continuous-time spin flips c_ε = c_vm + ε·c*, with asymptotic rate
  r^s(A) = g₁(1_A) (eq. 1.18): the O(ε) rate for a 0 to flip when the neighbours in A are 1.
  For us r^s(A) = g(|A|) = αR(|A|) − |A|/8 + κ[|A| = 2] (the −|A|/8 part drops out, Remark 5.4).
- **Θ₃ (eqs. 1.10, 1.19)** = Σ_A r^s(A)[Θ⁺(A) − Θ⁻(A)], where Θ± are sums of K₃(A₁,A₂,A₃),
  the constants in P(three given groups of walkers from N̄ have not met by t, but each group
  has coalesced) ~ K₃/(log t)³ (Prop. 1.6). **This is exactly the 3-class part of our
  f′(0)**: by Lemmas 5.1–5.2 it reduces to E_π[g(n₁) + g(n₂) − g(n₁+n₂)] over the K₃-weighted
  random 3-partition of N̄, i.e. our 3-class bracket. So **Θ₃ = 0 is our threshold at T → ∞**,
  where κ_a and κ_front merge.
  - Finite-T estimate: `vpt.py` now prints `kappa_th3/alpha` (3-class patterns only). It
    falls slowly: 0.308 (T = 64), 0.279 (1024), 0.272 (8192), 0.269 (65536). It must reach the
    same limit as κ_a (0.265, falling) and κ_front (0.255, rising), so the bracket
    [0.255, 0.265] stands. At T = 65536 the 4+-class terms still carry 28% of the weight.
  - Strict subadditivity of r^s gives Θ₃ > 0 (Cor. 5.3). D&N's rule is far from that:
    R(1) + R(2) − R(3) = −1. The j = 2 noise is what restores it, hence a κ*.
- **Theorems.**
  - Thm 1.15: low-density rescaling (ε_N = (log N)³/N) converges to super-Brownian motion with
    drift Θ₃. No monotonicity or cancellativity needed; exact up–down symmetry (ours) suffices.
  - Cor. 1.17: monotone + Θ₃ > 0 ⇒ a finite seed survives with positive probability, for
    ε small. Our w(j) is monotone in j iff p/8 ≥ (1−p)α (p ≥ 0.8 at α = 0.5; p ≥ 0.89 at α = 1).
  - Thm 1.9 (complete convergence with coexistence) additionally needs cancellativity, which
    is unchecked for D&N (the |N| = 8 check for the q-voter needed Maple, Appendix §11).
  - **Θ₃ < 0 (Remark 1.12) is only covered by the SBM limit with negative drift**: sparse 1s
    die. Nothing on large droplets, fronts, bistability or the transition order; for q > 1
    (founder control) "the tools … do not seem to be currently available".
- **What this means for us.**
  1. Cite C&P as the rigorous 2D framework and Θ₃ as the leading-order (T → ∞) form of our
     threshold. Our new content: the numerical value κ*/α = 0.26 for D&N/Moore, the finite-T
     corrections (4+-class terms) that open the [κ_front, κ_a] coexistence window, and the
     simulation check.
  2. At leading order in 2D there is a single threshold and no window. So the first-order
     character we see near p = 1 is a sub-leading effect, suppressed by powers of 1/log(1/ε).
     That matches "first order, weakening toward the voter point", and nothing rigorous
     contradicts it.
  3. Caveat: C&P (and CDP) are continuous time; dnsim is synchronous. Our cwalk uses the
     matching discrete-time dual, so the constants are for the synchronous model. Cor. 1.17
     applies, strictly, only to the continuous-time version of our rates.
- Their refs worth knowing: Agarwal–Simper–Durrett EJP 2021 (q-voter, d = 3 torus);
  Handjani AAP 1999 (threshold voter CCT, d ≥ 2); Sturm–Swart AAP 2008 (rebellious voter =
  1D Molofsky); Nowak, Stoń & Sznajd-Weron, Sci. Rep. 11, 6098 (2021) (discontinuous
  transitions in a noisy multi-state q-voter, mean-field).

## Physics: two symmetric absorbing states (all verified)
- Al Hammal, Chaté, Dornic & Muñoz, PRL 94, 230601 (2005): Langevin theory; GV or a split
  into Ising + DP. **d = 2 is the upper critical dimension of GV**: in d = 3 the direct
  transition is "fully first-order, in the sense that the density of interfaces now also
  jumps discontinuously to zero", while in 2D it is marginal. So our claim, in their words,
  is that D&N's 2D transition is *fully* first order (interface density jumps), not GV.
- Droz, Ferreira & Lipowski, PRE 67, 056108 (2003), "Splitting the voter criticality".
- Vazquez & López, PRE 78, 061127 (2008): microscopic rates → Langevin; GV/Ising/DP.
- **Castellano, Muñoz & Pastor-Satorras, PRE 80, 041129 (2009): the mean-field bistable
  regime of the q-voter model is "washed away" by fluctuations in 2D.** Full text read 10-16.
  - Their mean-field three-stable-state regime has the active state *at the origin* φ = 0
    of a single magnetization variable. The argument: "fluctuations can take the system from
    the origin to any of the absorbing states but not the other way around", so on large
    scales the origin plays no role and a single GV transition remains.
  - 2D evidence: one model (q = 4 on a square lattice, L = 5000) at one point
    (ε = 1/4): 1/ρ ∝ log t, the voter C(r,t) collapse, a linear exit probability, and a
    measured drift ⟨v(φ)⟩ renormalized to zero. Not a general theorem.
  - The response for D&N: the argument is zero-dimensional (well-mixed). In space the phase is
    picked by which phase invades the other. D&N has a front with a definite velocity of
    either sign (E1/E3), critical droplets (E5) and an L-independent active-state density;
    the active state is not "the origin" but has its own interface density (0.23), set by the
    rule. The PA shows the same: ρ is a fast variable with its own stable branch.
- Rodrigues, Chatelain, Tomé & de Oliveira, J. Stat. Mech. P01035 (2015): Ising + DP lines
  meeting at the voter point. Phase width is linear in the distance from the voter point,
  "possibly with logarithmic corrections in 2D". Full text read 10-16.
  - **They draw a "discontinuous" P–A line ending at the voter point, in 3D and on a 2D
    triangular lattice.** Closest prior art; must be cited and distinguished.
  - Their criterion is the magnetization jumping 0 → 1 with finite susceptibility, which
    any P→A transition has (the GV transition too). It was analysed on the cubic lattice
    only; the 2D diagram is "the same topology", from L ≤ 60, with no interface-density,
    front or nucleation analysis.
  - So in 2D they do not establish a *fully* first-order (Al Hammal sense) transition. Our
    E1–E8 data do: ρ_bond jumps to 0.23–0.29, fronts reverse sign, droplets have R_c, and
    lifetimes are L-independent.
- Azizi & Pleimling, PRE 102, 022112 (2020): Ising and DP lines merge into GV away from the
  voter point.
- Russell & Blythe, PRL 106, 165702 (2011): noise-induced voter → Ising coarsening crossover
  in 2D.
- Castellano & Pastor-Satorras, PRE 86, 051123 (2012): the 2D GV class.
- Llabrés, San Miguel & Toral, Phys. Rev. Research 8, 013015 (2026): a discontinuous
  "Modified GV" transition with a tricritical point, but only when the noise removes the
  absorbing states, and apparently mean-field.
- 1D precedent: Bagnoli, Boccara & Rechtman, PRE 63, 046116 (2001). A first-order line
  between the two absorbing states (not absorbing vs active), meeting DP lines.

## Life-like
- No noisy, stochastic or asynchronous Day & Night work found.
- Rendell (J. Cellular Automata): deterministic B34678/S3678, where islands grow or shrink
  depending on size.
- α-asynchronous Life (Blok & Bergersen, PRE 59, 3876, 1999): continuous, DP-like
  (unverified).

## Assessment (revised 10-16 after reading Cox & Perkins, Castellano, Al Hammal, Rodrigues)
- **Claim (a), first order in 2D, still appears novel, but phrase it precisely:** a *fully*
  first-order transition (interface density jumps; fronts, nucleation and droplets) in a 2D
  model with two symmetric absorbing states.
  - Prior art to distinguish: Rodrigues et al. 2015 draw a 2D discontinuous P–A line
    (magnetization-jump criterion only, L ≤ 60). Al Hammal et al. 2005 predict fully first
    order only for d ≥ 3, with 2D marginal. Castellano et al. 2009 and Molofsky et al. 1999 say
    bistability is washed out in 2D.
  - Answer to the washout argument: it is well-mixed. In 2D the phase is chosen by front
    motion, and D&N's active phase has its own rule-set interface density.
- **Claim (b), q\* = c·α·(1−p):**
  - The rigorous 2D framework exists: Cox & Perkins 2023 (Θ₃ from K₃ constants), on the
    coexistence side only.
  - Θ₃ = 0 is our T → ∞ threshold. Our contribution is its value for D&N/Moore (0.26), the
    finite-T window that makes the leading-order single threshold into a narrow first-order
    window, and the numerical check.
  - Rigorous results say nothing about Θ₃ < 0 beyond low-density extinction.
- **Still unread:** Cox & Durrett 1991 (book chapter; not online). Cox–Merle–Perkins 2010 and
  Sturm–Swart 2008 are downloaded but unread; low priority.
