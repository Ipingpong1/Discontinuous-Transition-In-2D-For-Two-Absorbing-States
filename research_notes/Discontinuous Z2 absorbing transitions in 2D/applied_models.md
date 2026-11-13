# Applied-model literature (opinion dynamics, spatial ecology, evolutionary games): two symmetric absorbing states on 2D lattices, and whether any has a discontinuous active–absorbing transition

Scope and legend. This file covers sociophysics, spatial ecology and evolutionary-game models. It does not repeat /home/matias/Desktop/silly/dn/campaign/theory/LITERATURE.md (Molofsky 1999, Cox–Durrett–Perkins, Cox–Perkins 2023, Al Hammal 2005, Rodrigues 2015, Azizi–Pleimling 2020, Russell–Blythe 2011, Bagnoli 2001, Vazquez–López 2008). It only adds new quotes or details for Castellano et al. 2009.

How each entry is marked:
- **Geometry**: [2D] = 2D lattice result, [NET] = network or random-graph result, [MF] = mean-field result. Only [2D] items can threaten the novelty claim.
- **Read**: FT = full text, or the relevant sections of it, read by me. ABS = abstract only. SEC = second-hand (a search summary or another paper's description).
- **Relevance**: H/M/L.
- **THREAT**: flags an item that bears on the claim that a *fully* first-order (interface-density-jumping) active→absorbing transition in a 2D Z2-symmetric two-absorbing-state model is new.

Bottom line:
- No applied-model paper found establishes such a transition, using interface-density jump, metastability/nucleation and L-independent lifetimes.
- Two published 2D items do use the words "first order" or "discontinuous" for a transition into two symmetric absorbing (sets of) states. The user must cite and distinguish both:
  - **Szabó & Sznaider (2004)**, with the Szabó & Fáth (2007) review, on the cyclic four-species "defensive alliances" model. They report a sudden drop of activity and a front velocity that vanishes linearly.
  - **Vazquez, Castelló & San Miguel (2010)**, on the Abrams–Strogatz volatility model on a square lattice. They call it "first order" from a deterministic Ginzburg–Landau potential, at the voter point.
- Separately, there is a terminology trap: de Oliveira et al. 1993, Dornic et al. 2001 and Castellano–Fortunato–Loreto 2009 call the 2D voter point itself a "first-order" (β = 0) critical point.

## Q1. q-voter and nonlinear voter models on 2D lattices without noise (both consensus states absorbing)

### Takeaway
Every noise-free 2D-lattice nonlinear-voter study found has its order–disorder transition at, or controlled by, the voter point, and is classified as generalized-voter (GV). This includes the q-voter, the Dornic kinetic-Ising family, majority/minority, symmetric Lotka–Volterra and Molofsky-type rules. The labels "first order"/"discontinuous" appear, but only in the β = 0 magnetization-jump sense: Vazquez–Castelló–San Miguel 2010 use the deterministic potential, and Llabrés et al. 2025 call GV itself "discontinuous". No study reports an interface-density jump, metastability or nucleation. The Sznajd-Weron/Timpanaro/Mobilia lines of work contain no noise-free 2D-lattice q-voter study. Their 2D or lattice work always includes independence, anticonformity, zealots or mass media (no absorbing states), or is 1D or on networks.

### Cited Findings
- **Castellano, Muñoz & Pastor-Satorras**, "The non-linear q-voter model", PRE 80, 041129 (2009); arXiv:0907.1775; doi:10.1103/PhysRevE.80.041129. [2D]+[MF].
  - Read: FT (local copy).
  - Model: square lattice L = 5000, random sequential. A voter consults q neighbours (with repetition) and flips if they agree; otherwise it flips with probability ε. Both consensus states are absorbing.
  - Details not in LITERATURE.md:
    - (i) On the 2D lattice "for several values of q > 1. In all cases, we find a transition separating a paramagnetic phase for large ε from a ferromagnetic one at low ε". Only q = 4 at ε = 1/4 was analysed in detail ("Analog results are found for other values of q (data not shown)").
    - (ii) In the no-repetition variant in d = 2, "for q < 4 no qualitative change occurs: there is a GV transition… Things change only in the special case q = 4, where all neighbors of a site are considered. The ordered phase disappears and for any value of ε a disordered state is reached." So the von Neumann q = 4 no-repetition rule has no ordered phase at all.
    - (iii) "For q = 2 and ε = 0 the q-voter model practically coincides with Sznajd model… Sznajd model is just a ferromagnetic model in its ordered phase."
  - Key washout quotes: "fluctuations renormalize the deterministic potential so that the transition is in the GV class for any value of q", and "b < 0 renormalizes on large scales to b = 0".
  - Evidence: one q, one ε, a single large L. H. Counter-position (not a threat). — [Source](https://arxiv.org/abs/0907.1775)
- **Dornic, Chaté, Chave & Hinrichsen**, "Critical coarsening without surface tension: the voter universality class", PRL 87, 045701 (2001); arXiv:cond-mat/0101202; doi:10.1103/PhysRevLett.87.045701. [2D].
  - Read: FT.
  - Model: square lattice, 4 nearest neighbours. A Z2-symmetric single-spin-flip "kinetic Ising" family with flip probability r_{s,h}, local field h ∈ {−4,−2,0,2,4}, and no bulk noise (r₄ = 0 ⇒ two absorbing states). Sizes up to 16384². Pair-update and m-conserving families were also studied.
  - Claims:
    - "We claim that this situation is generic: the 4-parameter space of the general kinetic Ising model above without bulk noise (space defined by r4 = 0) is divided into a disordered and an ordered region separated by a codimension-1 critical manifold of voter-like rules."
    - "This suggests to conjecture that all critical Z2-symmetric rules without bulk noise form a codimension-1 'voter-like' manifold separating order from disorder, characterized by the logarithmic decay of both ρ and m."
    - d = 2 "appears as the upper critical dimension where voter-like phase-ordering can be induced by the sole presence of interfacial noise."
  - Terminology point: citing de Oliveira, Mendes & Santos, J. Phys. A 26, 2317 (1993), they write that along the zero-bulk-noise line up to the voter point "the (static) magnetization ms jumps from 0 to ±1. Thus, the order parameter exponent β = 0, and, for such a first-order d = 2 nonequilibrium critical point…". Their conservative estimate is "β = 0, γ = 1, and ν = 1/2, with appropriate logarithmic corrections".
  - Evidence: simulations in several 2D planes of the z = 4 family only, with random initial conditions, coarsening, bubble and line experiments.
  - H. This is the most explicit counter-conjecture to the D&N claim. D&N is a Z2-symmetric rule without bulk noise (w(0) = 0), but it lies outside the explored family in three ways: z = 8; a Day & Night component that is non-subadditive and makes w(j) non-monotone whenever p/8 < (1−p)α; and synchronous updates. Not a threat. — [Source](https://arxiv.org/abs/cond-mat/0101202)
- **Vazquez, Castelló & San Miguel**, "Agent based models of language competition: macroscopic descriptions and order–disorder transitions", J. Stat. Mech. P04007 (2010); arXiv:1002.1251; doi:10.1088/1742-5468/2010/04/P04007. [2D]+[NET]+[MF].
  - Read: FT of §2.1 and §4.
  - Model: Abrams–Strogatz model with volatility a, P(s→−s) = ½(1−sv)(σ₋ₛ)^a, where σ₋ₛ is the fraction of opposite neighbours. Square lattice (4 neighbours), random sequential. v = 0 is Z2-symmetric. Both consensus states are absorbing, and a = 1 is the linear voter model.
  - Claims for the square lattice:
    - "the behavior of these models in square lattices is very different to their behavior in fully connected or complex networks… the mean-field approach that works well in fully connected networks gives poor results in lattices."
    - From a coarse-grained deterministic time-dependent Ginzburg–Landau equation: "The order-disorder nonequilibrium transition at a = 1 is reminiscent of the well known Ising model transition… An important difference is that the transition is here first order, since the low volatility stable sates φ = ±1 appear discontinuously at a = 1."
    - "Complete ordering for a > 1 is achieved through domain coarsening driven by surface tension". For a = 0.5 "the system reaches an active disordered stationary state" (128×128 snapshots).
    - "The coexistence-dominance transition in the ASM for v = 0 is at ac ≃ 1.0, as in fully connected and complex networks".
  - Evidence: a deterministic potential (no noise), snapshots and an applet. There is no analysis of interface density, fronts, droplets or finite-size scaling.
  - H. **PARTIAL THREAT (terminology)**: a published "first order" claim for a 2D Z2 nonlinear voter model with two absorbing states. Distinguishable: the transition sits exactly at the voter point, uses the magnetization-jump criterion, and the same group later classifies these models as GV (next item). — [Source](https://arxiv.org/abs/1002.1251)
- **Llabrés, San Miguel & Toral**, "Universality of noise-induced transitions in nonlinear voter models", arXiv:2505.11358 (2025); doi:10.1103/6xdk-2cjy. Journal per project notes: Phys. Rev. Research 8, 013015 (2026), not re-verified here. [MF] for the absorbing case; [2D]/[NET] for noisy models.
  - Read: FT of intro, §II.A and parts of §V and §VI.
  - Absorbing case:
    - "By mapping several models with symmetric absorbing states onto a canonical model… we confirm that they exhibit a Generalized Voter (GV) transition."
    - The analysis "focuses exclusively on the mean-field limit… (complete graph)".
    - GV is described as "characterized by a discontinuous transition from m = 0 to m = ±1."
    - Models mapped: nonlinear (Abrams–Strogatz-type) voter, nonlinear partisan voter, q-voter, Sznajd, and voter with aging.
  - Noisy case: the discontinuous "Modified GV" and the tricritical point exist only once noise removes the absorbing states. Tricritical scaling was studied on complete graphs and random networks. On the d = 2 lattice only the continuous Ising line was checked ("supporting their inclusion in the 2D Ising universality class").
  - M–H. Not a threat. It supports the point that "discontinuous" in this community means the β = 0 jump. — [Source](https://arxiv.org/abs/2505.11358)
- **Ramirez, Vazquez, San Miguel & Galla**, "Ordering dynamics of nonlinear voter models", PRE 109, 034307 (2024); arXiv:2310.20548. [NET]+[MF].
  - Read: ABS + text search; the journal reference was verified from the reference list of arXiv:2507.11521.
  - Model: rate ∝ (fraction of neighbours)^q, multi-state, on complete graphs and ER graphs, with a pair approximation. No 2D lattices.
  - Claim: "For q>1… The ordering to consensus is driven by deterministic drift"; for q < 1 "noise-driven succession of metastable states". L. — [Source](https://arxiv.org/abs/2310.20548)
- **Schweitzer & Behera**, "Nonlinear voter models: the transition from invasion to coexistence", EPJ B 67, 301 (2009); arXiv:cond-mat/0307742. [2D]+[MF]+PA.
  - Read: ABS + partial FT (local).
  - Model: m = 4 square-lattice stochastic CA, with Molofsky-type parameters (α₁, α₂).
  - Claims: three regimes, "complete invasion, random coexistence, and… correlated coexistence". Also: "the separatrix is divided into two pieces by the linear voter point, (0.2,0.4). Above that point, the separatrix is very narrow, but below the voter point it has in fact a certain extension in the parameter space". On the invasion side the system "converges into a 'frozen' state".
  - The order of the transition is not characterized. M. Not a threat. — [Source](https://arxiv.org/abs/cond-mat/0307742)
- **Lambiotte & Redner**, "Dynamics of vacillating voters", J. Stat. Mech. L10001 (2007); arXiv:0710.0914. [2D]+[MF].
  - Read: FT.
  - Model: square lattice. The flip probability is "0, 1/2, 5/6, and 1… when the number of anti-aligned neighbors is 0, 1, 2, and ≥ 3", a concave nonlinear voter with both consensus states absorbing.
  - Claims: "C1 reaches the stationary value 0.31, so that domains of opposite opinions coexist… and only a rare macroscopic fluctuation allows consensus to be reached"; "In two and higher dimensions, domains slowly anti-coarsen". This shows a genuine 2D active phase with finite interface density in a noise-free nonlinear voter model.
  - There is no tuning parameter, so no transition is studied. M. — [Source](https://arxiv.org/abs/0710.0914)
- **Lambiotte & Redner**, "Dynamics of non-conservative voters", EPL 82, 18007 (2008); arXiv:0712.0364. Read: ABS. 1D only ("a family of opinion formation models in one dimension"). L. — [Source](https://arxiv.org/abs/0712.0364)
- **Mobilia & Redner**, "Majority versus minority dynamics: phase transition in an interacting two-state spin system", PRE 68, 046106 (2003); arXiv:cond-mat/0306061. [all d].
  - Read: ABS + text search.
  - Model: groups of size G = 3 adopt the majority with probability p and the minority with probability 1−p. Both consensus states are absorbing.
  - Claim: "For group size G=3, there is a phase transition at p_c=2/3 in all spatial dimensions… For p=p_c… the system obeys classical voter model dynamics." The transition is pinned at the voter point. M. Not a threat. — [Source](https://arxiv.org/abs/cond-mat/0306061)
- **Stauffer, Castelló, Eguíluz & San Miguel**, "Microscopic Abrams–Strogatz model of language competition", Physica A 374, 835 (2007); arXiv:physics/0603042. [1D/2D/3D].
  - Read: ABS + text search.
  - Biased case (status s ≠ 1/2) with a = 1.31 and a = 1: one language dominates. There is no symmetric transition analysis. L. — [Source](https://arxiv.org/abs/physics/0603042)
- **Castelló, Eguíluz & San Miguel**, "Ordering dynamics with two non-excluding options: bilingualism in language competition", NJP 8, 308 (2006); arXiv:physics/0609079. [2D]+[NET].
  - Read: ABS.
  - Model: A/B/AB states; the only absorbing states are consensus A or B.
  - Claim: "Agents in the AB state define the interfaces, changing the interfacial noise driven coarsening of the voter model to curvature driven coarsening. We argue that this change in the coarsening mechanism is generic for perturbations of the voter model dynamics." No active phase. L–M. — [Source](https://arxiv.org/abs/physics/0609079)
- **Sznajd-Weron group and related q-voter lattice work**: nothing noise-free in 2D was found.
  - Przybyła, Sznajd-Weron & Tabiszewski, PRE 84, 031117 (2011), arXiv:1107.0366: 1D exit probability (ABS). — [Source](https://arxiv.org/abs/1107.0366)
  - Jędrzejewski, Sznajd-Weron & Szwabiński, "Mapping the q-voter model: from a single chain to complex networks", Physica A (2016), doi:10.1016/j.physa.2015.11.005, arXiv:1501.05091. Square lattices, but for the "modified q-voter model" (ABS; not noise-free per its framing). L. — [Source](https://arxiv.org/abs/1501.05091)
  - Timpanaro & Prado, PRE 89, 052808 (2014) (1D + large networks; SEC). Timpanaro & Galam, "An analytical expression for the exit probability of the q-voter model in one dimension", PRE 92, 012807 (2015), arXiv:1408.2734 (1D; SEC). **No "Timpanaro & Galesi" paper was found; the co-authors are Prado and Galam.** — [Source](https://arxiv.org/abs/1408.2734)
  - Mobilia, PRE 92, 012803 (2015). Mellor, Mobilia & Zia, EPL 113, 48001 (2016) and PRE 95(1) (2017), arXiv:1610.06092. These have zealots (no absorbing states) and are mean-field (SEC). L. — [Source](https://arxiv.org/abs/1610.06092)
  - Abramiuk-Szurlej, Lipiecki, Pawłowski & Sznajd-Weron, Sci. Rep. 11 (2021), doi:10.1038/s41598-021-97155-0. Generalized anticonformity (no absorbing states), on ER graphs (SEC). L. — [Source](https://arxiv.org/abs/2107.03338)
  - Woloszyn, Stauffer & Kułakowski, "Phase transition in Nowak–Sznajd opinion dynamics", Physica A 378, 453 (2007), arXiv:physics/0609178. Square lattice, but "with probabilities beta and gamma the opinions flip due to mass-media advertising" (no absorbing states). It finds "sharp phase transitions" (ABS). L. — [Source](https://arxiv.org/abs/physics/0609178)
  - Jędrzejewski & Sznajd-Weron, "Statistical physics of opinion formation: is it a SPOOF?", C. R. Physique 20, 244 (2019), arXiv:1903.04786. Read: text search. It has no 2D-lattice discontinuity statement; the q-voter discontinuities discussed ("discontinuous for q > 5") are mean-field. L. — [Source](https://arxiv.org/abs/1903.04786)
- **Coevolving nonlinear voter models** (networks, not lattices).
  - Jędrzejewski et al., PRE 102, 042313 (2020), arXiv:2004.08927: "Discontinuous phase transitions and hysteresis loops between both active phases are also possible" (adaptive networks, PA + MC; ABS).
  - Min & San Miguel, Sci. Rep. 7, 12864 (2017), arXiv:1707.08111.
  - Kim, Bae, Hur & Jeong, arXiv:2607.13616 (2026): random geometric graph with link deletion (ABS).
  - L. — [Source](https://arxiv.org/abs/2004.08927)
- **Physics-side items found incidentally** (for the stat-phys researcher): 2D interacting monomer–dimer model with two symmetric absorbing states.
  - Nam, Park, Kim & Lee, J. Stat. Mech. L06001 (2011), arXiv:1103.1222: "two nearby continuous transitions" (Z2 non-Ising + DP).
  - Park, PRE 85, 041140 (2012), arXiv:1202.2201, and Park, J. Stat. Mech. P10009 (2015), arXiv:1507.05537: argue the order–disorder transition is Ising.
  - All continuous. ABS only. — [Source](https://arxiv.org/abs/1202.2201)

### Inferences
- The explored 2D parameter families (z = 4 kinetic Ising, q-voter with q ≤ 4, Abrams–Strogatz with a, majority/minority G = 3, Molofsky 5-cell) all have rates that are monotone in the number of disagreeing neighbours, or close to linear. Their transitions are located at or emanate from the voter point. D&N's Day & Night component is strongly non-subadditive (R(1)+R(2)−R(3) = −1 per LITERATURE.md) and makes w(j) non-monotone unless p/8 ≥ (1−p)α. With z = 8 and synchronous updates, D&N sits outside every family that has actually been scanned. Dornic et al.'s "generic GV" is a conjecture from z = 4 data and does not cover it.
- "First order" in this literature (de Oliveira 1993; Dornic 2001; RMP 2009; Vazquez 2010; Llabrés 2025) consistently means the β = 0 magnetization jump at a critical (diverging-ξ, ν = 1/2) GV point. The D&N write-up must define first order by the interface-density jump, finite correlation length, metastability and nucleation (Al Hammal's "fully first order"). It should explicitly contrast Vazquez–Castelló–San Miguel's "first order" at a = 1.
- The vacillating-voter result shows that a stable, L-independent active phase with finite nearest-neighbour disorder (C1 = 0.31) exists in noise-free 2D nonlinear voters. So the D&N active phase is not exotic by itself; the novelty must rest on the discontinuity of the transition into it.

### Gaps
- The arXiv API was rate-limited (HTTP 429) part-way through, so searches for 2016–2026 Sznajd-Weron-group 2D papers were done by web search only. A noise-free 2D q-voter paper from that group could have been missed.
- de Oliveira, Mendes & Santos, J. Phys. A 26, 2317 (1993) was not read (known only via Dornic et al.).
- The Google Scholar / Semantic Scholar citation graphs of Castellano et al. 2009 and Dall'Asta & Castellano 2007 were not traversed systematically. Coverage relies on keyword searches plus the reference lists of the 2024–2026 reviews (arXiv:2507.11521; Llabrés 2025).

## Q2. Voter variants that acquire effective surface tension or bistability (noise-reduced, inertia/memory/aging, majority rule, Sznajd, Ising–voter interpolations, threshold)

### Takeaway
In 2D, memory, inertia, persistence, confidence, noise reduction and an intermediate (bilingual) state all generically turn voter (interfacial-noise) coarsening into Ising-like, curvature-driven coarsening with stripes. They produce ordering, not a new active phase. Where mean-field predicts a long-lived mixed state (confident voter), the 2D lattice instead shows surface tension and stripe states. None of these papers reports a discontinuous active–absorbing transition on a lattice.

### Cited Findings
- **Dall'Asta & Castellano**, "Effective surface-tension in the noise-reduced voter model", EPL 77, 60005 (2007); arXiv:cond-mat/0612186. [2D]+[3D]+[MF].
  - Read: ABS + text search.
  - Model: a spin changes only when a counter reaches a threshold r; consensus is absorbing.
  - Claims: "The numerical results for low-dimensional lattices show a scaling behavior in good agreement with usual Cahn-Allen curvature-driven coarsening". The interface density follows the t^(−1/2) slope "expected for curvature-driven" coarsening on a 2D lattice of L = 5000.
  - No active phase. L–M. — [Source](https://arxiv.org/abs/cond-mat/0612186)
- **Stark, Tessone & Schweitzer**, "Decelerating microdynamics can accelerate macrodynamics in the voter model", PRL 101, 018701 (2008); arXiv:0711.1133. [2D]+[NET]+[MF].
  - Read: ABS + text search.
  - Inertia grows with persistence time and saturates below 1 ("Choosing νs < 1 avoids trivial frozen states"). On 2D regular lattices (among others), consensus times decrease with inertia for µ < µ*.
  - No transition to an active phase. L. — [Source](https://arxiv.org/abs/0711.1133)
- **Volovik & Redner**, "Dynamics of confident voting", J. Stat. Mech. P04003 (2012); arXiv:1111.3883. [2D]+[MF].
  - Read: ABS.
  - MF: "a population of size N is quickly driven to a mixed state and remains close to this state before consensus is eventually achieved in a time of the order of ln N."
  - 2D: "two distinct times — one that scales linearly with N and another that appears to scale as N^{3/2}… stripe states arise from an effective surface tension".
  - Lattice-vs-MF contrast: the MF mixed state does not become a 2D active phase. M. — [Source](https://arxiv.org/abs/1111.3883)
- **Latoski, Dantas & Arenzon**, "Opinion inertia and coarsening in the Persistent Voter model", PRE 109, 054115 (2024); arXiv:2403.02483. Read: ABS. Transient, dynamically induced zealots give "an effective surface tension and the PVM has a crossover from a fluctuation-driven dynamics, as in the VM, to a curvature-driven one, as in the Ising Model at low temperature". L–M. — [Source](https://arxiv.org/abs/2403.02483)
- **Peralta, Khalil & Toral**, "Ordering dynamics in the voter model with aging", Physica A, art. 122475, doi:10.1016/j.physa.2019.122475; arXiv:1904.07128. [MF] only ("studied at the mean-field level"). Read: ABS. Anti-aging gives coexistence; aging gives consensus or a frozen state. L. — [Source](https://arxiv.org/abs/1904.07128)
- **Lambiotte, Saramäki & Blondel**, "Dynamics of latent voters", PRE 79, 046107 (2009); arXiv:0811.1464. [MF]+[1D]. Read: ABS. Latency drives toward zero magnetization; the majority rule with latency gives "a rich phenomenology". No 2D work. L. — [Source](https://arxiv.org/abs/0811.1464)
- **Krapivsky & Redner**, "Dynamics of majority rule in two-state interacting spin systems", PRL 90, 238701 (2003); arXiv:cond-mat/0303182. [d-dim lattices]+[MF]. Read: ABS. "On finite-dimensional lattices… the consensus time fluctuates strongly between realizations and grows as a dimension-dependent power of N." Consensus is always reached; there is no active phase. L. — [Source](https://arxiv.org/abs/cond-mat/0303182)
- **Chen & Redner**, "Majority rule dynamics in finite dimensions", PRE 71, 036101 (2005); arXiv:cond-mat/0408219. Read: ABS. "Repeated application of this update step leads to a coarsening mosaic of spin domains and ultimate consensus… with the longer time scale arising from realizations in which spins organize into coherent single-opinion bands." L. — [Source](https://arxiv.org/abs/cond-mat/0408219)
- **Lanchier & Neufer**, "Stochastic dynamics on hypergraphs and the spatial majority rule model", J. Stat. Phys. (2013), doi:10.1007/s10955-012-0543-5 (volume not verified); arXiv:1301.0151. [2D, rigorous + simulations]. Read: ABS. With n×n block updates, "opinion 1 wins when n is even while the system clusters when n is odd". L. — [Source](https://arxiv.org/abs/1301.0151)
- **Drouffe & Godrèche**, "Coarsening and persistence in a class of stochastic processes interpolating between the Ising and voter models", J. Phys. A 32, 249 (1999); arXiv:cond-mat/9807356. [2D].
  - Read: ABS + text search.
  - A two-temperature family (interfacial and bulk noise); "Special lines… correspond to the Ising model, voter model and majority vote model". On the zero-bulk-noise line (absorbing states), the voter point separates coarsening from the disordered region.
  - Only coarsening and persistence are studied (θ ≈ 0.22); the transition order is not studied. L. — [Source](https://arxiv.org/abs/cond-mat/9807356)
- **Corberi, dello Russo & Smaldone**, "Ordering kinetics with long-range interactions: interpolating between voter and Ising models", J. Stat. Mech. 093206 (2024); arXiv:2404.06917. Read: title/metadata only. Ordering kinetics only. L. — [Source](https://arxiv.org/abs/2404.06917)
- **Lipowski & Lipowska**, "Ising-doped voter model on complete graph", PRE 105, 024119 (2022); arXiv:2110.08800. [MF]+[NET]. Read: ABS. A noisy Ising/voter mixture with no absorbing states. L. — [Source](https://arxiv.org/abs/2110.08800)
- **Sznajd on the square lattice**: Castellano et al. 2009 state that Sznajd (q = 2, ε = 0) "is just a ferromagnetic model in its ordered phase" (FT; see Q1). Llabrés et al. 2025 map Sznajd onto the GV canonical form in MF. No active phase or discontinuous active–absorbing transition was found for the noise-free 2D Sznajd model. — [Source](https://arxiv.org/abs/0907.1775)
- **Reviews' framing**:
  - Starnini et al., "Opinion dynamics: statistical physics and beyond", arXiv:2507.11521v2 (May 2026). Read: relevant sections. "When noise-reduction is introduced, the surface tension emerges and drives the ordering dynamics (Dall'Asta and Castellano, 2007). Similarly, the coarsening of voter models with inertia or multiple states can also be driven by surface tension (Dall'Asta and Galla, 2008; Stark et al., 2008)." For q-voter extensions with anticonformity or independence, "consensus is no longer an absorbing state". — [Source](https://arxiv.org/abs/2507.11521)
  - Redner, "Reality-inspired voter models: a mini-review", C. R. Physique 20, 275 (2019); arXiv:1811.11888. Read: text search. No 2D discontinuity statements found. L. — [Source](https://arxiv.org/abs/1811.11888)
- **Threshold and conformity models**: covered in LITERATURE.md (Liggett 1994; Durrett–Steif 1993; Handjani 1999). Nothing new on the order of the transition was found. Noisy symmetric-threshold models (e.g. Nowak & Sznajd-Weron, PRE 101, 052316 (2020); SEC) have no absorbing states.

### Inferences
- The generic 2D effect of non-voter perturbations on the ordering side is surface tension (Allen–Cahn coarsening and stripes). This matches D&N's coarsening below q*. It means the absorbing side of D&N has a well-defined interface tension, which is a prerequisite for droplet and nucleation physics. In the GV scenario that tension vanishes at the transition (Lavrentovich–Nelson measure σ/k_BT ∼ |α|^0.61 → 0 at the voter point; see Q3). The D&N evidence of a finite R_c and interface tension at q* is therefore a sharp discriminator from GV.
- The confident-voter case is the clearest opinion-dynamics example of a mean-field metastable mixed state that does *not* survive as a 2D active phase.

### Gaps
- The 2D lattice behaviour of the noise-free aging voter (Peralta et al.) and of the latent voter was not found. Both are MF/1D only in what I read.
- Stauffer, Sousa & de Oliveira (2000, IJMPC; Sznajd on the square lattice) was not retrieved.

## Q3. Spatial ecology: founder control / bistable competition / mutualism on lattices; two-species Allee effects

### Takeaway
Spatial ecology's standard heuristic (Durrett–Levin 1994; Durrett 2009 "Case 2"; Lanchier 2015) is that mean-field bistability does *not* survive in space. Instead the state with the positive travelling-wave speed wins, so outcomes switch where the front speed changes sign. That is the same mechanism D&N invokes. Symmetric two-type models in 2D, however, put the coexistence (active) transition at the voter point and find voter- or GV-type behaviour:
- Neuhauser–Pacala Lotka–Volterra on the diagonal; Cox–Perkins 2008.
- Molofsky 1999.
- Symmetric mutualism with a 2D frontier (Lavrentovich & Nelson 2014).
- Pigolotti et al. 2013 in 2D, not characterized.

The one ecology model with two symmetric absorbing (sets of) states and a reported discontinuous transition in 2D is the cyclic four-species "defensive alliance" model (Szabó & Sznaider 2004). It is listed under Q4 with the evolutionary-game material because it comes from the same group and review.

### Cited Findings
- **Durrett**, "Coexistence in stochastic spatial models" (Special invited paper), Ann. Appl. Probab. 19, 477–496 (2009); doi:10.1214/08-AAP590. [mostly d = 2 IPS + heuristics].
  - Read: FT of §1–2.
  - Quotes:
    - "one goal of this paper is to explain the idea of Durrett and Levin (1994) that one can determine whether coexistence happens in the stochastic spatial model by examining properties of the mean-field ODE."
    - "With ecological problems in mind, we will usually take d = 2."
    - "Case 2. Two locally attracting fixed points. In this case, the limiting behavior of the ODE depends on the initial density. However, this is not the expected behavior for the particle system, and in this case the outcome of competition is dictated by the behavior of the PDE."
    - For the sexual-reproduction (single-absorbing) example: "the sign of c = the sign of ∫₀^ρ₂ g(y) dy".
  - M–H (supports the D&N mechanism; not a threat). — [Source](https://sites.math.duke.edu/~rtd/wald/waldsurvey.pdf)
- **Durrett & Levin**, "The importance of being discrete (and spatial)", Theor. Pop. Biol. 46, 363–394 (1994). Read: SEC (search summary + Durrett 2009 + Lanchier 2015). It compares mean-field ODEs, patch models, reaction–diffusion and IPS on three examples; the discrete stochastic models diverge from the differential-equation models in the third. M. — [Source](https://collaborate.princeton.edu/en/publications/the-importance-of-being-discrete-and-spatial/)
- **Durrett & Levin**, "Spatial aspects of interspecific competition", Theor. Pop. Biol. 53, 30–43 (1998). Read: SEC/ABS. Proves that a species with a competitive advantage over all others takes over; competition–dispersal trade-offs; a cyclic three-type case. It is not about symmetric bistability. L. — [Source](https://sites.math.duke.edu/~rtd/reprints/paper96.pdf)
- **Neuhauser & Pacala**, "An explicitly spatial version of the Lotka–Volterra model with interspecific competition", Ann. Appl. Probab. 9, 1226–1259 (1999). [Z^d].
  - Read: SEC (abstract via search) + Cox–Perkins 2008 intro.
  - Per the abstract, local competition shrinks the founder-control region, shrinks the coexistence region and produces spatial segregation where the ODE predicts coexistence.
  - Cox & Perkins (AAP 18, 747 (2008), arXiv:0804.0158, FT of intro) write: "Neuhauser and Pacala have conjectured coexistence holds along the diagonal α0 = α1 for α0 < 1… but Theorem 1.3 and the above Conjecture suggest that for d = 2 their conjecture is quite delicate near α0 = 1", where d = 2 is "the critical dimension". The symmetric model is Z2-symmetric with two absorbing states; the active/coexistence boundary on the diagonal sits at the voter point α = 1. The order is not characterized. M. — [Source](https://arxiv.org/abs/0804.0158)
- **Lavrentovich & Nelson**, "Asymmetric mutualism in two- and three-dimensional range expansions", PRL 112, 138102 (2014); arXiv:1309.0273. [2D frontier = d = 2+1].
  - Read: FT.
  - Model: generation-by-generation (synchronous) stepping-stone/Domany–Kinzel-type updates on stacked triangular lattices. Two cell types; fixation of either is absorbing; symmetric when α = β.
  - Claims:
    - "When d = 2 + 1, two-dimensional domains at the voter model point α = β = 0 lack a surface tension and readily 'dissolve'… these dynamics allow for a mutualistic phase for all α = β > 0, with a remarkable pinning of the corner of the 'wedge' of mutualism… to the origin."
    - "When r = 0, we find a voter model transition at s = sc(r = 0) = 0, and any s > 0 pushes the system into a mutualistic phase with a non-zero steady-state domain interface density."
    - For α = β < 0: "dynamics similar to a kinetic Ising model… interface density decay ⟨ρ⟩ ∼ t^(−1/2)… an effective line tension σ". "As we approach the voter model point (α → 0−), σ/kBT vanishes with an apparent power law σ/kBT ∼ |α|^0.61."
  - H as counter-evidence: a 2D synchronous PCA with two symmetric absorbing states and an active (mutualistic) phase gives a voter-point transition, not first order. Not a threat. — [Source](https://arxiv.org/abs/1309.0273)
- **Pigolotti, Benzi, Perlekar, Jensen, Toschi & Nelson**, "Growth, competition and cooperation in spatial population genetics", Theor. Popul. Biol. (2013; volume/pages not verified); arXiv:1208.4973. [1D]+[2D].
  - Read: FT of §4.
  - 2D: "Mutualism in 2d is computationally challenging and, to the best of our knowledge, has not been studied systematically… along the cut εA = εB [the heterozygosity] displays a transition from a regime in which it seems to decay logarithmically… to a regime in which fixation does not occur." The order is not characterized.
  - 1D: "the dashed diagonal line below the cusp… is like a first order phase transition". This is a jump between the two *absorbing* outcomes, not active↔absorbing. M. — [Source](https://arxiv.org/abs/1208.4973)
- **Korolev & Nelson**, "Competition and cooperation in one-dimensional stepping stone models", PRL 107, 088103 (2011); arXiv:1104.2464. [1D]. Read: ABS. "mutualism is generically lost via a directed percolation process, with a phase diagram strongly influenced by an exceptional DP2 transition." L (1D). — [Source](https://arxiv.org/abs/1104.2464)
- **Dall'Asta, Caccioli & Beghè**, "Fixation–coexistence transition in spatial populations", EPL 101, 18003 (2013). Read: SEC/ABS. The transition between fixation and coexistence is shown "for a one-dimensional system" via a modified voter model. L (1D). — [Source](https://epljournal.edpsciences.org/10.1209/0295-5075/101/18003)
- **Molofsky, Bever & Antonovics**, "Coexistence under positive frequency dependence", Proc. R. Soc. B 268, 273–277 (2001); doi:10.1098/rspb.2000.1355. [2D].
  - Read: ABS + methods via text search.
  - Model: 100×100 grid, generation (synchronous) updates, five-square neighbourhood, positive frequency dependence (a nonlinear voter of majority type). The neutral case is "equivalent to the voter model".
  - Claim: "This coexistence arises from the generation of banding patterns that will be stable over long time-periods". This is frozen or long-lived stripes, not an active phase. L. — [Source](https://beverlab.ku.edu/sites/beverlab/files/documents/56.pdf)
- **Borile, Maritan & Muñoz**, "The effect of quenched disorder in neutral theories", J. Stat. Mech. P04032 (2013); arXiv:1304.0965. Read: ABS. In nonlinear voter models, quenched random-field disorder "does not allow the separation of the up-down (Ising-like) symmetry breaking from the active-to-absorbing phase transition in low-dimensional systems… as a consequence of the well-known Imry-Ma argument". M (relevant if disorder is added to D&N). — [Source](https://arxiv.org/abs/1304.0965)

### Inferences
- In the ecology language, D&N is a "Case 2" system in which the bistable pair is (active mixed state, absorbing consensus) rather than (species 1, species 2). Durrett–Levin and Lanchier predict that the front-speed sign selects the phase, giving a first-order switch where the speed is zero. The 2D GV literature (Dornic, Castellano) asserts that near the voter point fluctuations destroy this, and every symmetric 2D ecology model found (NP Lotka–Volterra, Lavrentovich–Nelson, Molofsky 1999) does pin the transition at the voter point. So the ecology heuristic supports D&N's mechanism. The physics of 2D symmetric absorbing states is what makes it non-trivial and potentially novel.
- Rigorous results (Cox–Durrett–Perkins; Cox & Perkins 2008/2023) treat d = 2 as the critical dimension. No ecology paper computes a first-order boundary near the 2D voter point.

### Gaps
- No two-species Allee-effect lattice model with two symmetric absorbing states was found. This was searched only briefly (an arXiv query was blocked by rate limiting). Single-absorbing-state Allee/facilitation lattice models (e.g. Kéfi et al. TPB 2007; Villa Martín et al. PNAS 2015) were not covered (adjacent scope).
- The full texts of Neuhauser & Pacala 1999 and Durrett & Levin 1994 were not read; the Project Euclid PDF fetch returned HTML.
- Priority-effect metacommunity models (patch-occupancy, not IPS) were not surveyed.

## Q4. Evolutionary games and cyclic-dominance models on lattices with imitation/invasion updates (absorbing homogeneous states): any discontinuous transition between consensus/absorbing and a mixed active state?

### Takeaway
Discontinuous transitions between absorbing or homogeneous phases and active mixed phases are routine on 2D square lattices in multi-strategy spatial games. The Szabó/Szolnoki/Perc school locates them by front-invasion ("stripe") experiments, the same diagnostic as D&N's front-velocity sign change. Most of these games lack Z2 symmetry. The exception is **the even-species cyclic predator–prey model**:
- It has two symmetric "defensive alliance" absorbing sets.
- The activity drops "suddenly" to zero and the invading front's velocity vanishes linearly, v_inv ∼ (µcr − µ).
- The Szabó–Fáth review says "The Monte Carlo data indicates a first-order phase transition".
- The original 2004 paper instead says the transition "is equivalent to those found in spreading models with two equivalent absorbing states" (the voter class).

This is the closest prior art found for "2D, two symmetric absorbing states, first-order, linear front velocity" and must be cited and distinguished. For 2×2 Z2-symmetric coordination games under imitation, no study of a discontinuous transition into a mixed active phase on a lattice was found. Lattice results show a front- or cluster-advance-selected winner (risk dominance) instead of mean-field bistability.

### Cited Findings
- **Szabó & Sznaider**, "Phase transition and selection in a four-species cyclic predator-prey model", PRE 69, 031911 (2004); arXiv:q-bio/0310017 (v1 titled "…four-species cyclic Lotka-Volterra model"). [2D].
  - Read: FT (arXiv v1).
  - Model: square lattice (L = 2000–4000), random sequential updates, four species with cyclic dominance (1→2→3→4→1) and vacancies at concentration ρ₀. Individuals migrate to empty neighbouring sites or invade prey. "the mixed states containing only two neutral species (1+3 or 2+4) are also absorbing states and will be denoted as D13 and D24. In these stationary states the ratio of the two species remains constant." The two alliances are exchanged by the cyclic symmetry.
  - Claims:
    - Abstract: "Above the treshold, a symmetry breaking ordering occurs via growing domains containing only two neutral species inside… According to our Monte Carlo simulations the observed phase transition is equivalent to those found in spreading models with two equivalent absorbing states although the present model has continuous sets of absorbing states".
    - Text: "the invasion activity (or Ppp) decreases and drops suddenly to zero at ρ0 = ρcr. The discontinuous transition is accompanied with enhanced fluctuations (in all the quantities we studied) and a critical slowing down."
    - "The dynamical systems with two equivalent absorbing states represent a curious universality class (named after the voter model)… whose general features are consistent to those described above."
    - Front experiments: "the C state invades the territories of defensive alliances for low concentration of vacant sites. The average invasion velocity decreases monotonously with ρ0 and becomes zero at ρ0 = ρcr."
  - Evidence: simulations only, with no finite-size scaling of the jump, no droplet/nucleation tests and no interface-density analysis. The authors themselves are ambiguous between "discontinuous" and "voter class".
  - H. **PARTIAL THREAT** (2D, two symmetric absorbing sets, a discontinuity claimed, a sign-changing front velocity). Distinguishing points:
    - The absorbing states are continuous families (fixed composition, still mixing), not frozen consensus configurations.
    - The symmetry is the cyclic Z4 → Z2 alliance exchange, not spin-flip.
    - The active phase is cyclic invasion (a rock–paper–scissors-type phase) with vacancies.
    - The order was never established quantitatively, and the abstract calls it voter-class equivalent. — [Source](https://arxiv.org/abs/q-bio/0310017)
- **Szabó & Fáth**, "Evolutionary games on graphs", Phys. Rep. 446, 97 (2007); arXiv:cond-mat/0607344. [2D, review].
  - Read: FT of §8.1, the §6 PD passage and text search.
  - §8.1 ("A four-species cyclic predator-prey model with local mixing"): square lattice, neutral pairs exchange sites with probability µ, and "µcr = 0.02662(2) (Szabó, 2005)".
    - "The Monte Carlo data indicates a first-order phase transition at µ = µcr."
    - "Above the threshold value of µ, simulations display a domain growth process, which is analogous to the one described by the kinetic Ising model below its critical temperature".
    - "this system has four homogeneous stationary states, two continuous sets with two-species, and one symmetric stationary solution with four-species".
    - "For µ < µcr the four-species association expands, but its average invasion velocity vanishes linearly as the critical mixing rate is approached (Szabó and Sznaider, 2004), vinv ∼ (µcr − µ)."
    - On approximations: "The traditional mean-field approximation cannot account for site exchange… The two- and four-site approximations also fail… A qualitatively correct prediction can be achieved when the generalized mean-field approximation is performed on a 3 × 3 cluster."
  - In the section on the PD on the square lattice with Fermi imitation: "the four-site approximation predicts a first-order phase transition, i.e., a sudden jump to ρ = 0 at b = bc2, while the more accurate nine-site approximation predicts linear extinction" (the true lattice transition is DP). This is a cautionary example of spurious first order from cluster approximations.
  - H. **PARTIAL THREAT**: an explicit "first-order" statement plus linear front velocity for a 2D model whose inactive phases are two symmetric absorbing sets. Distinguish as above. — [Source](https://arxiv.org/abs/cond-mat/0607344)
- **Szabó, Szolnoki & Sznaider**, "Segregation process and phase transition in cyclic predator-prey models with even number of species", arXiv:0710.3371 (2007; journal reference not verified). [2D].
  - Read: ABS + FT of the transition paragraphs.
  - Model: n = 4, 6, 8 species, square lattice, site exchange with probability X, no vacancies.
  - Claims:
    - For n = 4 the four-species state survives "if X < Xc(4) = 0.02662(2)… At the critical value of X a sudden change occurs in both quantities".
    - For n = 6 "The predator-prey probability (ppp) drops suddenly to zero at X = Xc(6) = 0.00654(2)".
    - For n = 8 there is "blocking of this segregation process within a range of X".
    - "more accurate determination of the critical values and systematic analysis of the corresponding phase transitions exceed our computing capacity."
  - H (same family; the order is still not established). — [Source](https://arxiv.org/abs/0710.3371)
- **Helbing, Szolnoki, Perc & Szabó**, "Punish, but not too hard: how costly punishment spreads in the spatial public goods game", NJP 12, 083005 (2010); arXiv:1007.0431. [2D].
  - Read: FT of §3.1–3.2.
  - Model: square lattice, Fermi imitation, four strategies (C, D, PC, PD); homogeneous phases are absorbing; no Z2 symmetry.
  - Claims:
    - "Dashed blue lines indicate discontinuous first-order phase transitions, while solid red lines represent continuous, second-order phase transitions."
    - "The system always leaves the punishment-free state via a discontinuous first-order phase transition, while the transition between the mixed D+PC phase and the PC-only phase is continuous", and the latter "agrees with the directed percolation universality class".
  - M. Shows that first-order active/absorbing transitions are common in 2D multi-strategy games. Not Z2, so no direct threat. — [Source](https://arxiv.org/abs/1007.0431)
- **Lanchier**, "Evolutionary games on the lattice: payoffs affecting birth and death rates", Ann. Appl. Probab. 25, 1108–1154 (2015); arXiv:1302.0069. [Z^d, d ≥ 1, rigorous + simulations].
  - Read: ABS + FT of §1.
  - Claims:
    - "the inclusion of local interactions induces a reduction of the coexistence region of the replicator equation and the presence of a dominant strategy that wins even when starting at arbitrarily low density in the region where the replicator equation displays bistability."
    - "simple heuristic arguments looking at the interface between two adjacent blocks of the two strategies suggest that the most selfish strategy… always wins… The fact that bistability in the mean-field model results in the presence of a strong type in the interacting particle system has already been observed for a number of models" (citing Durrett 2009 and Durrett–Levin 1994).
    - Coexistence occurs in a "thorn-shaped region starting at the bifurcation point".
  - M (lattice-vs-MF statement; not a threat). — [Source](https://arxiv.org/abs/1302.0069)
- **Roca, Cuesta & Sánchez**, "Effect of spatial structure on the evolution of cooperation", PRE 80, 046106 (2009); arXiv:0806.1649. [2D lattices k = 4, 8 + networks]. Read: ABS + text search. Stag-hunt cooperation is promoted on lattices. The explanation is cluster growth and "the advance conditions of a flat in[terface]". There is no study of the order of the transition. L–M. — [Source](https://arxiv.org/abs/0806.1649)
- **Roca, Cuesta & Sánchez**, "Promotion of cooperation on networks? The myopic best response case", EPJ B (2009), doi:10.1140/epjb/e2009-00189-0; arXiv:0901.0355. Read: ABS + text search. "nontrivial results differing from the well-mixed population in the case of coordination games on lattices, which we explain in terms of the formation of spatial clusters and the conditions for their advancement". L–M. — [Source](https://arxiv.org/abs/0901.0355)
- **Szabó, Szolnoki, Varga & Hanusovszky**, "Ordering in spatial evolutionary games for pairwise collective strategy updates", PRE 82, 026110 (2010); arXiv:1008.1063. Read: ABS. "In the case of the anti-coordination game this system behaves alike an anti-ferromagnetic kinetic Ising model". The update allows random strategy pairs, so there are no absorbing states. L. — [Source](https://arxiv.org/abs/1008.1063)
- **Hauert**, "Effects of space in 2×2 games", Int. J. Bifurcation Chaos 12, 1531–1548 (2002). Citation verified; **not read**. Roca et al. (2009) call it a "preliminary and inconclusive attempt" (SEC). — [Source](https://personal.math.ubc.ca/~hauert/publications/reprints/hauert_ijbc02.pdf)

### Inferences
- The cyclic even-species family is the most dangerous prior art for the *headline* phrasing, more so than any opinion-dynamics model. It was described in 2D with "discontinuous"/"first-order" wording and with a front velocity vanishing linearly at the threshold, which is exactly D&N's E1/E3 diagnostic. The defensible novelty for D&N becomes:
  - frozen Z2 consensus states under spin-flip symmetry (not neutral two-species mixtures);
  - a nonlinear-voter (imitation-only) rule;
  - first order established quantitatively: interface-density jump, L-independent metastable lifetimes, critical droplets and nucleation barriers, wetting;
  - persistence of the first-order character close to the voter point (p = 0.995).
- The Szabó school's "stripe/invasion of competing solutions" method is an established 2D first-order locator. D&N should cite it (Szabó & Sznaider 2004; Szabó & Fáth 2007 §8.1, Eq. (163)) as methodological precedent rather than present front-velocity reversal as a new diagnostic.

### Gaps
- No study of a Z2-symmetric 2×2 game (a = d, b = c) under pure imitation on a 2D lattice that maps the coordination→anti-coordination (ordered→active) transition and its order was found. Cox & Durrett (d ≥ 3) and Cox–Perkins 2023 (2D, coexistence side) are the closest; both are in LITERATURE.md.
- The source of µcr = 0.02662(2) ("Szabó, 2005") and whether any later paper re-examined the order of the four-species transition with finite-size scaling were not verified. A citation-graph check of Szabó & Sznaider 2004 is recommended.

## Q5. Explicit lattice-vs-mean-field statements: do mean-field discontinuities survive on 2D lattices?

### Takeaway
The statements split by model class:
- Inertial majority vote (Z2, no absorbing states): the MF and network discontinuity does **not** survive on 2D lattices with k = 4, 8 (Moore) or 12. It reappears only for k = 20 (up to fourth neighbours), without hysteresis but with bimodal distributions and 1/N shifts.
- q-voter with absorbing states (Castellano 2009; Dornic 2001; Llabrés 2025 for MF): bistable or "three-stable-state" regimes are washed out in 2D to a GV transition.
- Spatial ecology and games (Durrett 2009; Lanchier 2015): MF bistability turns into front-speed-selected single outcomes. This implies first-order-type switching rather than washout.
- Cluster and pair approximations can produce spurious first-order transitions (Szabó–Fáth: four-site vs nine-site; Abramiuk-Szurlej: PA vs simulations on ER graphs).

### Cited Findings
- **Encinas, Harunari, de Oliveira & Fiore**, Sci. Rep. 8, 9338 (2018) (volume and article number as given in the assignment); doi:10.1038/s41598-018-27240-4; arXiv:1802.05256, titled "Fundamental ingredients for the emergence of discontinuous phase transitions in the majority vote model". [2D]+[NET]+[MF].
  - Read: FT.
  - Model: majority vote with inertia θ (flip probability depends on the local field including the own spin). Random-regular networks and 2D square lattices with interactions to first (k = 4), first+second (k = 8, i.e. Moore), first–third (k = 12) and first–fourth (k = 20) neighbours. Noisy, so no absorbing states.
  - MF: "MFT predicts a continuous phase transition for k = 4 irrespective the value of θ"; for k = 8 "the coexistence line arises only when θ > 1/3 and is very tiny (ff − fb is about 2·10−4)".
  - 2D, k = 4, 8, 12: "Thereby, there is an important difference between random and regular structures: The phase transitions are continuous irrespective the inertia" (Binder-cumulant crossing at the Ising value U₀* ≈ 0.61).
  - 2D, k = 20: "the phase transition becomes discontinuous… Contrary to the random complex case, hysteresis is absent… and the order-parameter distribution exhibits a bimodal shape". Pseudo-transition points scale as N^(−1); "The phase coexistence occurs for θ > 1/3, larger than θ > 3/13 (RR structure)".
  - Conclusions: "discontinuous transitions in the MV also manifest in low dimensional regular topologies. Also, its finite size behavior… is identical to that exhibited by discontinuous phase transitions into absorbing states", and "For random regular networks, we found that a minimum neighborhood is k = 7, whereas about k = 20 are required for changing the order of transition in bidimensional lattices."
  - H for the washout argument. On a Moore (z = 8) 2D lattice, the inertial-MV mean-field discontinuity is washed out. D&N claims first order *at z = 8*, so the user should explain why D&N differs (rule non-monotonicity, absorbing states, synchronous updates). Not a threat (no absorbing states). — [Source](https://arxiv.org/abs/1802.05256)
- **Chen, Shen, Zhang, Li, Hou & Kurths**, inertial majority vote, arXiv:1609.00469 (titled "Explosive phase transition in a majority-vote model with inertia"). [NET]+[MF].
  - Read: ABS + text search.
  - Topology: ER, random-regular and scale-free networks only; no lattices.
  - Claim: "the order-disorder phase transition is changed from a usual continuous type to a discontinuous or an explosive one when the inertia is above an appropriate level… strong hysteresis… a disordered phase and two symmetric ordered phases are coexisting".
  - **Citation conflict**: Encinas et al., Harunari et al. and Chen et al. (Chaos 2017) cite it as PRE 95, 042304 (2017). The arXiv journal-ref field says PRE 97, 062304 (2018) (doi:10.1103/PhysRevE.97.062304). Unresolved; the published title "First-order phase transition in a majority-vote model with inertia" is unverified. L–M. — [Source](https://arxiv.org/abs/1609.00469)
- **Harunari, de Oliveira & Fiore**, "Partial inertia induces additional phase transition in the majority vote model", PRE 96, 042305 (2017); arXiv:1706.07465. [NET]+[MF]. Read: ABS. Explosive transitions in homogeneous and heterogeneous networks; "the classification of phase transition depends on the inertia and the lattice topology". L. — [Source](https://arxiv.org/abs/1706.07465)
- **Chen, Shen, Zhang & Kurths**, "Large deviation induced phase switch in an inertial majority-vote model", Chaos 27, 081102 (2017); arXiv:1708.06052. [MF/WKB]. Read: ABS. Switching times are exponential in N within the hysteresis region; the coexistence line is where forward and backward actions intersect. L. — [Source](https://arxiv.org/abs/1708.06052)
- **Castellano, Muñoz & Pastor-Satorras (2009)**, [2D] vs [MF] (FT): "This third (so far unexplored) scenario… is rationalized and found to be specific of mean-field, i.e. fluctuations are explicitly shown to wash it out in spatially extended systems", and "as fluctuations can take the system from the origin to any of the absorbing states but not the other way around, effectively, on sufficiently large scales, the stable state at the origin plays no role, and the system exhibits a single ordering transition in the generalized-voter class." — [Source](https://arxiv.org/abs/0907.1775)
- **Castellano, Fortunato & Loreto**, "Statistical physics of social dynamics", Rev. Mod. Phys. 81, 591 (2009); arXiv:0710.3256. Read: §III.A.
  - "It has been argued (Dornic et al., 2001) that voter behavior is generically observed at order-disorder non-equilibrium transitions, driven by interfacial noise, between dynamically symmetric absorbing states. This symmetry may be enforced either by an up-down symmetry of the local rules or by global conservation of the magnetization. The universal exponents associated to the transition are β = 0, and ν = 1/2 in all dimensions, while γ = 1/2 for d = 1, and γ = 1 for d > 2 with logarithmic corrections at the upper critical dimension d = 2."
  - Also: "As the strength of bulk noise is decreased, the length grows and the voter first-order transition occurs for zero noise." — the "first-order" label applied to the voter point.
  - H (terminology and counter-position). — [Source](https://arxiv.org/abs/0710.3256)
- **Vazquez, Castelló & San Miguel (2010)**, [2D] vs [MF]/[NET]: "the mean-field approach that works well in fully connected networks gives poor results in lattices… pair approximation does not provide a good enough description of the dynamics in lattices either", yet "The coexistence-dominance transition in the ASM for v = 0 is at ac ≃ 1.0, as in fully connected and complex networks". — [Source](https://arxiv.org/abs/1002.1251)
- **Volovik & Redner (2012)**, [MF] vs [2D] (ABS): a quick MF approach to a long-lived mixed state, versus 2D "stripe states [that] arise from an effective surface tension". — [Source](https://arxiv.org/abs/1111.3883)
- **Durrett (2009)**, [MF] vs [IPS, d = 2]: "Two locally attracting fixed points… this is not the expected behavior for the particle system, and in this case the outcome of competition is dictated by the behavior of the PDE." — [Source](https://sites.math.duke.edu/~rtd/wald/waldsurvey.pdf)
- **Lanchier (2015)**, [MF] vs [Z^d]: "the presence of a dominant strategy that wins even when starting at arbitrarily low density in the region where the replicator equation displays bistability." — [Source](https://arxiv.org/abs/1302.0069)
- **Szabó & Fáth (2007)**, cluster approximations vs [2D] MC: the four-site approximation predicts a spurious "first-order phase transition, i.e., a sudden jump to ρ = 0" in the lattice PD, whereas nine-site and MC give linear (DP) extinction. In the cyclic model, MF and small-cluster approximations "fail" and only 3×3 is qualitatively right. — [Source](https://arxiv.org/abs/cond-mat/0607344)
- **Abramiuk-Szurlej et al. (2021)**, [NET] (SEC from abstract): discontinuities in the generalized-anticonformity q-voter can persist on ER graphs with small mean degree. Pair approximation predicts discontinuous transitions that simulations do not show when q_a ≥ q_c − 1. — [Source](https://arxiv.org/abs/2107.03338)
- **Llabrés, San Miguel & Toral (2025)**, [MF] vs [2D]: the discontinuous MGV line and the tricritical point were analysed on complete graphs and ER networks. On d = 2 lattices only the continuous Ising line was confirmed for the noisy nonlinear voters. — [Source](https://arxiv.org/abs/2505.11358)
- **Lavrentovich & Nelson (2014)**, [2D frontier]: the symmetric mutualistic transition in d = 2+1 is a voter-model transition pinned at zero benefit, versus a finite-benefit threshold in d = 1+1 (DP2 point at α = β = αc = 0.1242(5)). — [Source](https://arxiv.org/abs/1309.0273)

### Inferences
- The most useful contrast for the D&N paper is Encinas et al.:
  - In the inertial majority vote, a Z2 model whose MF already predicts a discontinuity for k = 8, the discontinuity does *not* survive on the 2D Moore lattice and needs k ≈ 20.
  - D&N claims a robust first-order transition at k = 8 *with* absorbing states, where theory (Al Hammal; Dornic; Castellano) predicts even stronger washout (2D is GV's upper critical dimension).
  - Either D&N's ingredients (a non-subadditive Day & Night component, the j = 2 term, synchronous updates, absorbing states) produce a genuinely different mechanism, or the first-order signals should be re-checked for very large L and long times. Encinas et al.'s criteria (absence or presence of hysteresis, bimodality, N^(−1) shift of pseudo-transition points) are a natural extra test.
- The ecology and game "front-speed selection" statements (Durrett; Lanchier; Szabó) and the opinion-dynamics "fluctuation washout" statements (Castellano; Dornic) are not contradictory. The washout argument concerns the *magnetization* of a mixed state at the origin. The front argument concerns phases with their own finite correlation lengths. D&N's response in LITERATURE.md (an active phase with its own interface density, set by the rule) is consistent with this reading.

### Gaps
- No study was found that tests whether the q-voter's MF discontinuity (with independence, q > 5) survives on a 2D square lattice. Web searches returned only MF and random-graph work.
- Vieira & Crokidakis, Physica A 450, 30 (2016) (majority vote with independence on square lattices; cited in SPOOF) was not read.
- de Oliveira, da Luz & Fiore, PRE 92, 062126 (2015) (finite-size scaling for discontinuous absorbing transitions, which Encinas et al. invoke) was not read here (adjacent scope).
