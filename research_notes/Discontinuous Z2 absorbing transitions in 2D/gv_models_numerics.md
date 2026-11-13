# Microscopic lattice models with two Z2-symmetric absorbing states: numerics, transition types, and prior first-order claims (focus d = 2)

Scope and method:
- Researcher: simulation/numerics literature. Done 2026-10-27.
- Read levels used below:
  - FULL means the arXiv full text was downloaded and the relevant sections read.
  - PARTS means targeted sections were read.
  - ABSTRACT means only the abstract or metadata was read.
- Forward citations came from the Semantic Scholar Graph API and the OpenAlex API. Google Scholar was not accessed. The arXiv API returned HTTP 429 (rate-limited), so arXiv full-text search was not possible. Keyword sweeps used OpenAlex title/abstract search (20 queries, 707 records).
- Project context assumed:
  - The model is Day & Night (B3678/S34678) on the Moore neighbourhood (z = 8), mixed with voter steps and j = 2 noise.
  - In terms of j disagreeing neighbours, D&N flips iff j ∈ {3,6,7,8}. This is a non-monotone rule.

---

## Key question 1 — Definitive answer: has a discontinuous (first-order) active–absorbing transition with two Z2-symmetric absorbing states been reported in d = 2?

### Takeaway
**Not in the strong ("fully first-order") sense. No 2D short-range lattice model with exactly two Z2-related homogeneous absorbing configurations has been reported with a jump in the interface/activity density backed by front, nucleation or metastability evidence.**

There are, however, four kinds of prior "discontinuous/first-order" statements that the project must cite and distinguish:
- **(a) Lipowski & Droz 2002** called d = q = 2 "most likely … first order". The evidence was a ρ jump. Its own dynamic Monte Carlo (MC) disagreed, and the authors withdrew the claim in 2003: the model is in the voter (GV) class.
- **(b) Rodrigues et al. 2015** label a 2D (triangular-lattice) P–A line "discontinuous". The jump is the voter-type magnetization jump. Azizi & Pleimling 2020 nevertheless describe it as "a line of discontinuous phase transitions … replacing the generalized voter critical line".
- **(c) Guisoni, Loscar & Girardi 2013** (abstract only) report "two lines of irreversible first-order phase transition between the two absorbing phases and the active phase". The model is a four-species cyclic predator–prey lattice model. This is the closest genuine precedent, but its absorbing "phases" are neutral-pair alliances, not two single configurations.
- **(d) The GV point itself** is routinely called "first-order critical point", "first-order d = 2 nonequilibrium critical point", "abrupt" or "discontinuous", because the magnetization jumps 0 → ±1 with β = 0.

Direct counter-evidence also exists:
- Al Hammal et al. 2005: 2D Langevin.
- Azizi & Pleimling 2020: a microscopic model at b = −8 is still GV.
- Castellano et al. 2009 (q-voter) and Molofsky et al. 1999.
- Droz et al. 2003: the LD Potts model on the Moore z = 8 neighbourhood is GV.

### Cited Findings

**Partial precedents that claim discontinuity in 2D:**
- **Lipowski & Droz 2002 (2D, q = 2, square lattice, z = 4)** — FULL.
  - On the steady state: "now the jump in ρ is almost twice smaller than previously [q = 3]. We should be also aware of the fact that an observed jump might be a finite-size effect."
  - On their reading: "In our opinion, from the steady-state and time-dependent measurements of ρ the most likely possibility is that the transition is of first order and the plateau will develop … but only at a larger time scale."
  - Their spreading runs disagree: P(t) ∼ t^−δ′ with δ′ = 0.90(2), and η = 0. The abstract concludes "Simulations are inconclusive, however, in the two-dimensional case for q = 2." — [arXiv:cond-mat/0112130](https://arxiv.org/abs/cond-mat/0112130); [PRE 65, 056114](https://doi.org/10.1103/PhysRevE.65.056114)
- **The same authors withdrew this in 2003.**
  - Review: "Such an unusual behaviour explains the numerical difficulties in studying models of this universality class [19]. The nonequilibrium Potts model for q = 2 and on square lattice with nearest neighbour interactions also belongs to the voter universality class [22]." — [Droz & Lipowski, Braz. J. Phys. 33, 526 (2003), arXiv:cond-mat/0308419](https://arxiv.org/abs/cond-mat/0308419)
  - DFL 2003 footnote: "We thank H. Chaté for drawing our attention to the voter model." — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **Rodrigues, Chatelain, Tomé & de Oliveira 2015** (2D triangular lattice and 3D cubic lattice, both z = 6) — FULL (local copy).
  - "As to the PA line, that separates the P and A phases, it is a discontinuous transition line that ends at the voter point."
  - "Along the PA line the jump in the order parameter is a constant and equals the jump observed in the voter model." — [arXiv:1501.01116](https://arxiv.org/abs/1501.01116); [J. Stat. Mech. P01035](https://doi.org/10.1088/1742-5468/2015/01/P01035)
- **Azizi & Pleimling 2020 read Rodrigues et al. as a genuinely discontinuous 2D line:** "the model of Rodrigues et al. only has two critical lines as a line of discontinuous phase transitions is replacing the generalized voter critical line. In addition, the Ising and directed percolation critical lines end exactly at the point where the classical voter model is realized." — [arXiv:2007.10933](https://arxiv.org/abs/2007.10933)
- **Guisoni, Loscar & Girardi 2013** — ABSTRACT only. The full text is paywalled; the CONICET repository PDF returned HTTP 401.
  - "By applying a gradient in the invasion rate parameter one can see, in the same simulation, the presence of two symmetric absorbing phases, composed by neutral pairs, and an active phase that includes all four species."
  - "We show that the model presents two lines of irreversible first-order phase transition between the two absorbing phases and the active phase. Depending on the value of the system parameters, these lines can converge into a triple point, which is the beginning of a first-order irreversible line between the two absorbing phases, or end in two critical points belonging to the directed percolation universality class. Standard simulations for some characteristic values of the parameters confirm the order of the transitions as determined by the gradient method." — [PRE 88, 022133](https://doi.org/10.1103/PhysRevE.88.022133); [CONICET record hdl:11336/7499](https://ri.conicet.gov.ar/handle/11336/7499)

**"Discontinuous" used for the GV point itself (magnetization jump only):**
- de Oliveira, Mendes & Santos 1993 (2D square lattice): "The phase diagram in the parameter space displays a critical line that terminates at a first-order critical point." — [J. Phys. A 26, 2317](https://doi.org/10.1088/0305-4470/26/10/006)
- Dornic et al. 2001: "the (static) magnetization ms jumps from 0 to ±1. Thus, the order parameter exponent β = 0, and, for such a first-order d = 2 nonequilibrium critical point, the susceptibility and the correlation length exponents should satisfy γ = 2ν." — [arXiv:cond-mat/0101202](https://arxiv.org/abs/cond-mat/0101202)
- Krause, Böttcher & Bornholdt 2012 (2D kinetic Ising with suppressed bulk noise): "we characterize the abrupt phase transition". — [PRE 85, 031126](https://doi.org/10.1103/PhysRevE.85.031126)
- Krause & Bornholdt 2012: "the discontinuous generalized voter model phase transition". — [PRE 86, 056106](https://doi.org/10.1103/PhysRevE.86.056106)

**Counter-evidence (2D: GV, not fully first order):**
- **Al Hammal, Chaté, Dornic & Muñoz 2005** (2D Langevin numerics) — FULL.
  - "For b < b∗, a unique transition is observed at a = aGV, across which the steady-state magnetization jumps from zero to ±1. On the other hand, the density of interfaces goes continuously to zero as a → aGV from below (not shown)."
  - Only for d = 3: "indeed it is 'fully first-order', in the sense that the density of interfaces now also jumps discontinuously to zero. We note that this is in agreement with d = 2 being the upper critical dimension for the GV class."
  - The 2D b < 0 example shown is b = −0.2, with D = 0.5 and σ² = 0.8. — [arXiv:cond-mat/0411693](https://arxiv.org/abs/cond-mat/0411693)
- **Azizi & Pleimling 2020** (2D microscopic cell model with tunable Langevin a, b) — PARTS.
  - At the GV point: ρ ∼ 1/ln t, N_f ∼ 1/ln t (η_GV = 0), δ_GV ≈ 1.
  - "The same behavior is encountered for b < 0. This is illustrated … for the case b = −8 for which the phase transition is found to take place at a = −4.132(2)." — [arXiv:2007.10933](https://arxiv.org/abs/2007.10933); [PRE 102, 022112](https://doi.org/10.1103/PhysRevE.102.022112)
- **Droz, Ferreira & Lipowski 2003**, LD Potts model with z = 8 (the Moore neighbourhood): "either the model belongs to the voter universality class or there is only extremely small splitting beyond the resolution of our simulations." — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **Dornic et al. 2001 conjecture:** "This suggests to conjecture that all critical Z2-symmetric rules without bulk noise form a codimension-1 'voter-like' manifold separating order from disorder, characterized by the logarithmic decay of both ρ and m." — [arXiv:cond-mat/0101202](https://arxiv.org/abs/cond-mat/0101202)
- **Castellano, Muñoz & Pastor-Satorras 2009:** the mean-field bistable regime of the 2D q-voter model is washed out, leaving a GV transition (project notes; full text read 10-16). — [arXiv:0907.1775](https://arxiv.org/abs/0907.1775); [PRE 80, 041129](https://doi.org/10.1103/PhysRevE.80.041129)
- **Molofsky, Durrett, Dushoff, Griffeath & Levin 1999:** no "multiple outcome regime" in the 2D synchronous PCA (project notes; full text read). — [Theor. Pop. Biol. 55, 270](https://doi.org/10.1006/tpbi.1998.1404)

**Context: first-order absorbing transitions are rare in low d.** Ódor's review: "First order phase transitions have rarely been seen in low dimensions." — [Ódor, Rev. Mod. Phys. 76, 663 (2004), arXiv:cond-mat/0205644](https://arxiv.org/abs/cond-mat/0205644)

### Inferences

**Claim (a) survives only in a narrowed form.** A defensible version is: *a fully first-order (Al Hammal sense) active–absorbing transition in a d = 2 short-range lattice model whose only absorbing states are the two Z2-related homogeneous configurations, demonstrated by interface-density jump, front-velocity sign change, critical droplets and L-independent metastable lifetimes.*
- I found no paper that meets all of these.
- Guisoni et al. 2013 is the main item that could be read as pre-empting the broad phrase "first first-order transition with two symmetric absorbing states in 2D". It must be cited and distinguished on four points:
  - Its absorbing phases are alliance families of configurations. In the parent model, Szabó & Sznaider describe them as "continuous sets of absorbing states".
  - Neutral-pair exchange continues inside them, so they are not frozen.
  - The evidence is a gradient-method interface analysis plus "standard simulations".
  - The lattice dimension is not stated in the abstract. 2D is probable (the parent model is on a square lattice) but unverified.

**Avoid ambiguous vocabulary.** "Discontinuous" alone is ambiguous in this literature: the magnetization at the GV point is "discontinuous" by convention (β = 0).
- The project should state explicitly that the *activity/interface density* ρ jumps, and contrast this with Al Hammal 2005's statement that in 2D ρ "goes continuously to zero".

**The historical pitfall is exactly LD 2002.** GV's slower-than-power-law decay mimicked a ρ jump.
- The project's front velocity, droplets and lifetimes are qualitatively stronger than anything in LD 2002, Rodrigues 2015, or (from its abstract) Guisoni 2013.
- Showing the ρ-jump and droplet data at several L, with log-time plots ruling out ρ ∝ 1/ln t at q*, would pre-empt the obvious "GV in disguise" objection.

**The strongest objection is Azizi–Pleimling's b = −8 GV result.** It sits at a = −4.132 ≈ b/2. That reading is my inference: it would be the mean-field coexistence point of the single-field potential. If so, even a deep mean-field bistability gives GV in a 2D microscopic model.
- The project's rebuttal must rest on D&N's active phase not being the φ = 0 state of a single-field Langevin. It has its own rule-generated interface density (0.23–0.29), effectively a second slow field.

**Neighbourhood size is not the explanation.** In the Lipowski–Droz family, enlarging the range splits GV into Ising + DP (z = 12) and gives GV at z = 8. The D&N Moore model is first order at the same z = 8.
- So the qualitative difference must come from the rule structure. The LD Metropolis rule is monotone in j (derivation below); D&N's is non-monotone. Synchronous versus random-sequential updating is not the cause, since the project checked asynchronous updates.

### Gaps
- **Guisoni et al. 2013 full text is not accessible** (APS paywall; repository bitstream restricted). Still unknown:
  - the exact rates;
  - the lattice dimension;
  - whether the first-order lines are supported by hysteresis, bimodality, metastability or fronts;
  - whether the alliance phases are frozen.
- Highest-priority item to obtain (library or authors).
- Henkel, Hinrichsen & Lübeck (2008), *Non-Equilibrium Phase Transitions Vol. 1*: not accessed, so any book statement on first order with two absorbing states in 2D is unverified.
- Google Scholar "cited by" was not accessed. Semantic Scholar/OpenAlex coverage is incomplete (Semantic Scholar returned 0 citers for LD 2002; OpenAlex returned 17).

---

## Key question 2 — Droz, Ferreira & Lipowski 2003: full content (z = 4, z = 8, z = 12, 3D) and any first-order statement

### Takeaway
Full text read (arXiv v1).
- The model is the q = 2 "absorbing-state Potts" (modified Metropolis) model.
- z = 4 is in the voter class.
- **z = 8 (up to 2nd neighbours, i.e. Moore) is "voter class or … extremely small splitting beyond the resolution"**, with no data shown.
- z = 12 splits into Ising + DP.
- The 3D nearest-neighbour model also splits.
- There is **no first-order claim for q = 2**. The only "discontinuous" remark refers to q ≥ 3.

### Cited Findings
- **Citation.**
  - M. Droz, A. L. Ferreira, A. Lipowski, published as "Splitting the voter Potts model critical point" (Crossref title; arXiv title "Splitting the voter criticality"). Phys. Rev. E 67, 056108 (2003), DOI 10.1103/PhysRevE.67.056108, arXiv:cond-mat/0301291 (v1, 16 Jan 2003). Read: FULL. — [arXiv](https://arxiv.org/abs/cond-mat/0301291); [DOI](https://doi.org/10.1103/PhysRevE.67.056108)
  - Relevance: **High**, because it uses the same lattice and the same z = 8 as D&N. Threat: **none** (no first order for q = 2).
- **Model.** Spins σ = 0, 1.
  - Dynamics: "(i) Select randomly the site i and its possible final state. (ii) If ∆E < z update the site i with the probability min{1, e^{−∆E/T}}". So a site whose z neighbours all agree with it cannot change.
  - "the unit of time is defined as a single (on average) update of every site", i.e. random-sequential.
  - Measured: density of active sites ρ (a site is active if a neighbour differs), magnetization m, variance s, Binder cumulant U, and spreading P(t) and N(t). — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **(A) z = 4 (square lattice, nearest neighbours).** The variance s diverges where ρ vanishes, at Tc = 1.7585 (L = 500), with s⁻¹ linear, i.e. γ = 1. This is voter class. — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **z = 8 (second neighbours = Moore).** "First we studied the model with interactions up to second-nearest neighbours (z = 8). It turned out, that in this case either the model belongs to the voter universality class or there is only extremely small splitting beyond the resolution of our simulations." No figure or numbers are given. — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **(B) z = 12 (up to third neighbours).**
  - s diverges at T_I, clearly above Tc. The Binder crossing is near the 2D Ising value U = 0.6107 (t ∼ 5·10⁶ MCS, L = 100–300).
  - The absorbing transition is at Tc = 4.7380(5) with β = 0.61(4), η = 0.25(3), δ = 0.44(3). The DP values are 0.584, 0.230 and 0.451. Statics used L = 400–500; spreading used L = 5000 and 2·10⁴ runs.
  - Interpretive line: "the z = 12 model also has two absorbing states and is driven only by interfacial noise. Nevertheless, the critical behaviour in this case is much different from the voter model … simple criterions most likely cannot be used for classifying nonequilibrium critical points." — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **(C) 3D, simple cubic, nearest neighbours (z = 6).** "here again the absorbing and magnetic phase transitions are separated" (L = 60). 3D Ising + (3+1) DP is expected, but "Confirmation … is left for the future." — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **Motivation.** A pair approximation showed that "for sufficiently large coordination number z a structure of the solution qualitatively changes in a way that clearly indicates two separate transitions" (ref. [19], "in preparation"). — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **The only first-order/discontinuous statement** concerns q > 2: "there are some indications that for d > 1 models with more than two symmetric absorbing states the phase transition should be discontinuous [9]", citing LD 2002. — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **On GV's β:** "Although the phase transition in the voter universality class is continuous, the decay of the order parameter ρ upon approaching the critical point is slower than any power law … Formally such a decay might be described with β = 0." — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **Summary sentence:** "models with two absorbing states in d > 1 dimensions might exhibit two transitions where the first one breaks the symmetry and the second one brings the model into an absorbing state. In the voter model and some related models these two transitions coincide." — [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- **Later confirmations.**
  - Azizi, Stidham & Pleimling 2018 confirm the z = 12 split with dynamic and aging quantities. — [arXiv:1809.08366](https://arxiv.org/abs/1809.08366)
  - Azizi & Pleimling 2020: "In numerical studies of kinetic Ising models [19, 22, 23] the different scenarios predicted by the Langevin equation are not found for a fixed range of the spin-spin interactions. Only when modifying the interaction range can one observe a change in scenario." — [arXiv:2007.10933](https://arxiv.org/abs/2007.10933)

### Inferences
- **Derivation of the DFL rule** (mine, from the stated dynamics). For q = 2 the only possible final state is the other spin. With j disagreeing neighbours, ΔE = z − 2j, so the flip probability is w(j) = min{1, e^{−(z−2j)/T}} for j ≥ 1, and w(0) = 0.
  - At z = 8: w(1) = e^{−6/T}, w(2) = e^{−4/T}, w(3) = e^{−2/T}, and w(j ≥ 4) = 1.
  - This is a **monotone, majority-like rule** with thermal minority flips.
  - D&N instead has w ∝ [j ∈ {3,6,7,8}]. It is non-monotone: it does not flip at j = 4, 5 and does flip at j = 3.
- **The DFL z = 8 result is the most direct "control" for the project.** The same lattice and neighbourhood, with a monotone rule, gives GV (or an unresolvable split), never first order.
- This supports the project's framing that first order comes from the D&N rule structure (e.g. R(1) + R(2) − R(3) < 0, i.e. the subadditivity violation in the project's notes), not from z = 8 per se.

### Gaps
- DFL show no z = 8 data, so how close z = 8 is to a split is unknown.
- The promised pair-approximation paper ("[19] in preparation") was not found among the 27 Semantic Scholar / 32 OpenAlex citers. It may never have appeared.

---

## Key question 3 — The Lipowski–Droz program (q = 2 vs q ≥ 3 in d = 2), Lipowski 2000, and follow-ups

### Takeaway
In d = 2 (square lattice, z = 4, random-sequential modified Metropolis):
- **q = 3 and q = 4 are discontinuous.** This is clear and was later reconfirmed with a jump in both M and ρ.
- **q = 2 was first called inconclusive / "most likely first order", then reclassified as GV.**
- A larger range (z = 12) turns even q = 3 into two continuous transitions: 2D Potts(3) + DP.

Lipowski 2000 is about infinitely many absorbing states with a gauge symmetry, not Z2. It is irrelevant to the novelty question except as a "first-order with power-law singularity" curiosity.

### Cited Findings
- **Lipowski & Droz, "Phase transitions in nonequilibrium d-dimensional models with q absorbing states".** Phys. Rev. E 65, 056114 (2002), DOI 10.1103/PhysRevE.65.056114, arXiv:cond-mat/0112130 (v2, 26 Feb 2002). FULL. Relevance **High**. Threat: **weak/historical** (claim withdrawn).
  - **Model.** q-state Potts Metropolis on the square lattice (nearest neighbours) and the simple cubic lattice. "when all neighbours of a given site are in the same state as this site, then this site cannot change its state".
  - **d = 2, q = 3.** "For T < 1.237 it is virtually impossible to reach an active steady state value of ρ, which suggests that the model undergoes a discontinuous phase transitions." Supporting evidence:
    - an infinitely long plateau of ρ(t) at T = 1.237 (L = 500);
    - P(t) and N(t) bending with no power law (up to 10⁸ runs).
  - "Thus, as a summary, for the d = 2, q = 3 case our model undergoes a discontinuous phase transition."
  - **d = 2, q = 4:** "the jump is larger and the first-order character of the transition is even more transparent. We expect that such a behaviour persists for all q > 3."
  - **d = 2, q = 2:** Tc = 1.7585(3). The steady-state ρ shows a jump about half the q = 3 one, and ρ(t) has no clear plateau. A local slope gives δ′ = 0.07(2), possibly 0. Spreading gives δ′ = 0.90(2) and η = 0. Their reading is "most likely … first order". They suggest d = q = 2 lies "close to or at the crossing point" of non-mean-field, mean-field and first-order regions.
  - **d = 3, q = 2:** "continuous, transition with β = 1 (mean-field)" (L = 40, 60).
  - They propose a (q, d) diagram "very similar to that of the equilibrium Potts model". — [arXiv:cond-mat/0112130](https://arxiv.org/abs/cond-mat/0112130); [DOI](https://doi.org/10.1103/PhysRevE.65.056114)
- **Droz & Lipowski, "Some recent developments in models with absorbing states".** Braz. J. Phys. 33, 526–532 (2003), DOI 10.1590/S0103-97332003000300014, arXiv:cond-mat/0308419. PARTS. Relevance **High** (the authors' consolidated view).
  - "(ii) d = 2, q = 3 In this case (nearest-neighbour interactions) there is a clear evidence of the discontinuous phase transition. In particular, the steady-state density of active sites ρ has a discontinuous behaviour and the time dependent ρ(t) develops a plateau at a critical point."
  - d = 2, q = 2 is voter class (quoted in KQ1).
  - "(iii) d = 3, q = 2 … increasing the dimensionality up to d = 3 also results in two separate phase transitions [22]". This revises the 2002 β = 1 claim. — [arXiv:cond-mat/0308419](https://arxiv.org/abs/cond-mat/0308419); [DOI](https://doi.org/10.1590/S0103-97332003000300014)
- **Azizi, Stidham & Pleimling, "Dynamic critical properties of non-equilibrium Potts models with absorbing states".** J. Stat. Mech. (2018) 103208, DOI 10.1088/1742-5468/aae2dd, arXiv:1809.08366. PARTS (intro, models, q = 3 sections, summary). Relevance **High** for q ≥ 3; Medium for q = 2. Threat: none for q = 2.
  - **q = 3, z = 4:** "Both order parameters display a jump at the same temperature" (Fig. 5, L = 400). The text gives T ≈ 1.273 in one place and T = 1.237 elsewhere; the latter matches LD 2002, so 1.273 is likely a typo.
  - Two-time autocorrelations show no dynamical scaling, like the equilibrium Potts model with q = 7 or 8.
  - **q = 3, z = 12:** "the single discontinuous transition encountered for only nearest-neighbor interactions is replaced by two continuous phase transitions when extending the interaction range to third nearest neighbors: a symmetry-breaking phase transition belonging to the universality class of the two-dimensional q = 3 Potts model and an absorbing phase transition belonging to the Directed Percolation universality class."
  - **q = 2:** z = 4 is GV and z = 12 splits. They also simulate the KBB GV model. — [arXiv:1809.08366](https://arxiv.org/abs/1809.08366); [DOI](https://doi.org/10.1088/1742-5468/aae2dd)
- **Lipowski & Droz, "Dynamics-dependent criticality in models with q absorbing states".** Phys. Rev. E 66, 016106 (2002), arXiv:cond-mat/0110404. ABSTRACT. 1D only. Relevance Low.
  - "For q = 2, as expected, the model belongs to the parity-conserving universality class. For q = 3 the critical behavior depends on the dynamics of the model. Under a certain dynamics it remains generically in the active phase, which is also the feature of some other models with three absorbing states. However, a modified dynamics induces a parity-conserving phase transition." — [DOI](https://doi.org/10.1103/PhysRevE.66.016106)
- **Lipowski, "First-order transition with power-law singularity in models with absorbing states".** Phys. Rev. E 62, 4401 (2000), arXiv:cond-mat/0002378. FULL. Relevance **Low**. Threat: none.
  - Model: bond variables w_ij ∈ (−0.5, 0.5) on the 2D square lattice and a 1D triangular ladder. A site is active if the product of its bonds is < r. Updates are random-sequential.
  - There are **infinitely many absorbing states** with a "gauge-like symmetry", not Z2.
  - "jump of the order parameter coincides with its power-law singularity". The transition is at r_c = 0 exactly, with ρ = ρ0 + A r^β: in 2D, ρ0 = 0.359 and β = 0.58(1); in 1D, ρ0 = 0.314827 and β = 0.66(3). — [arXiv:cond-mat/0002378](https://arxiv.org/abs/cond-mat/0002378); [DOI](https://doi.org/10.1103/PhysRevE.62.4401)
- **Later Lipowski-group work, on random graphs (mean-field-like).** Relevance Low.
  - "Analogous phase transitions for models with Metropolis dynamics are discontinuous" (heat-bath gives mean-field Ising). — [Lipowski, Ferreira & Lipowska, Entropy 25, 1615 (2023), arXiv:2312.02398](https://doi.org/10.3390/e25121615)
  - Ising–voter model: "For the Metropolis dynamics, the phase transition is most likely discontinuous, at least for not too many antivoters" (abstract snippet via OpenAlex). — [PRE 111, 024317 (2025)](https://doi.org/10.1103/PhysRevE.111.024317)
- **Lipowski, Lipowska & Ferreira, "Phase transition and power-law coarsening in Ising-doped voter model".** Phys. Rev. E 96, 032145 (2017), arXiv:1706.03575. ABSTRACT. Relevance Low (Ising agents presumably remove absorbing states; unverified).
  - "even a very small fraction (∼ 1%) of the Ising agents drastically changes the behaviour of the Voter model … a logarithmically slow coarsening of the Voter model (d = 2), or its active steady state (d = 3), change into an Ising-type power-law coarsening." — [arXiv:1706.03575](https://arxiv.org/abs/1706.03575)

### Inferences
- In this model family, discontinuity in 2D is tied to **q ≥ 3 with short range**. Increasing the range (z = 12) removes it even for q = 3.
- D&N shows the opposite trend: first order at q = 2 with z = 8. The project can present this as a qualitatively new route to first order, driven by the shape of the rule rather than by the number of absorbing states.
- The LD Metropolis-type models on random graphs become discontinuous (Lipowski 2023/2025). That is a mean-field-level hint that Metropolis-like or threshold-like rule structure, not detailed balance per se, controls the order. It does not bear on 2D.

### Gaps
- Pair-approximation results for large z (promised in DFL 2003 ref. [19]) were not found.
- No LD-model study at z = 8 with q = 3 was found.

---

## Key question 4 — Map of microscopic models with two symmetric absorbing states (2D focus; 1D and 3D for context)

### Takeaway
Every 2D microscopic model with exactly two homogeneous Z2-symmetric absorbing states that I could find reports either GV (direct transition, log decays) or a split into Ising + DP. There is a still-unresolved Ising vs non-Ising dispute in the interacting monomer–dimer models.

"Discontinuous" in 2D appears only in three places:
- the magnetization-jump sense (de Oliveira 1993; Dornic 2001; KBB 2012; Rodrigues 2015);
- one alliance-type cyclic predator–prey model (Guisoni 2013, abstract only);
- one withdrawn claim (LD 2002).

In 1D, first-order lines occur only *between* the two absorbing states, plus front-velocity jumps. In 3D, Langevin predicts "fully first-order", which is preliminary. Microscopic 3D models show a split or a "discontinuous" P–A line in the magnetization sense.

### Cited Findings

#### Summary map (R = relevance H/M/L; T = threat to novelty)

| Model / paper | d, lattice, z, update | Absorbing states | Transition reported | Evidence | Read | R / T |
|---|---|---|---|---|---|---|
| General z = 4 kinetic Ising with no bulk noise ([Dornic et al. 2001](https://arxiv.org/abs/cond-mat/0101202)) | 2D square, z = 4 | 2, homogeneous | Codimension-1 GV manifold; GV point is a "first-order d = 2 nonequilibrium critical point" (magnetization jump) | MC, several planes, L up to 16384² | FULL | H / none (the conjecture concerns critical rules) |
| [de Oliveira, Mendes & Santos 1993](https://doi.org/10.1088/0305-4470/26/10/006) | 2D square, z = 4 | 2 on the zero-bulk-noise line | Ising line ending at a "first-order critical point" (voter) | Pair approximation + MC | ABSTRACT | M / none |
| [Drouffe & Godrèche 1999](https://doi.org/10.1088/0305-4470/32/2/003) | 2D, 2 "temperatures" (interfacial and bulk noise) | 2 when bulk noise = 0 | Coarsening; persistence exponent | MC | ABSTRACT | L–M / none |
| LD Potts q = 2, z = 4 ([LD 2002](https://arxiv.org/abs/cond-mat/0112130); [DFL 2003](https://arxiv.org/abs/cond-mat/0301291)) | 2D square, random-sequential Metropolis | 2 | 2002: "inconclusive / most likely first order"; 2003: GV | MC, L ≤ 600 | FULL | H / weak (withdrawn) |
| LD Potts q = 2, **z = 8 (Moore)** ([DFL 2003](https://arxiv.org/abs/cond-mat/0301291)) | 2D square, z = 8, random-sequential | 2 | GV, or a split too small to resolve | MC (no data shown) | FULL | **H** / none |
| LD Potts q = 2, z = 12 ([DFL 2003](https://arxiv.org/abs/cond-mat/0301291); [Azizi et al. 2018](https://arxiv.org/abs/1809.08366)) | 2D, z = 12 | 2 | Split: Ising (Binder) + DP | MC | FULL / PARTS | H / none |
| Hinrichsen 1997 2D model ([cond-mat/9608065](https://arxiv.org/abs/cond-mat/9608065)) | 2D, 80×80 | 2 | δ = 0.9(1), η = 0.00(5), z = 1.0(1), "roughly … mean-field"; later understood as GV | MC (small) | PARTS | M / none |
| 2D generalized contact process ([Lee & Vojta 2011](https://arxiv.org/abs/1010.3298)) | 2D square | 2 | Generic transition GV; σ = 0 line "not critical" (trivially) | MC up to 20000², t = 3·10⁶ | PARTS | M / none |
| Interacting monomers ([Park 2012](https://arxiv.org/abs/1202.2201)) | 2D | 2 (sublattice) | Split: Ising order–disorder + absorbing transition | MC | ABSTRACT | M / none |
| Interacting monomer–dimer with nearest-neighbour repulsion ([Nam et al. 2011](https://doi.org/10.1088/1742-5468/2011/06/L06001); [Nam et al. 2014](https://doi.org/10.1088/1742-5468/2014/08/P08011); [Park 2015](https://arxiv.org/abs/1507.05537)) | 2D square | 2 (symmetric) | Split; non-Ising (Nam) vs Ising (Park) | MC | ABSTRACT | M / none |
| Nonlinear voter models ([Vazquez & López 2008](https://arxiv.org/abs/0810.1436)) | 2D square | 2 | GV / Ising / DP | Langevin from rates + MC | ABSTRACT (+skim) | H / none |
| q-voter ([Castellano et al. 2009](https://arxiv.org/abs/0907.1775)) | 2D square | 2 | GV; mean-field bistability "washed away" | MC, L = 5000 | project FULL | H / counter-evidence |
| Three GV-class models ([Castellano & Pastor-Satorras 2012](https://arxiv.org/abs/1207.4333)) | 2D | 2 | GV; non-trivial exit probability | MC | PARTS | M / none |
| KBB kinetic Ising ([Krause et al. 2012](https://doi.org/10.1103/PhysRevE.85.031126)) | 2D square | 2 | GV, "abrupt" | Fokker–Planck, mean field, MC | ABSTRACT | M / terminology |
| Russell–Blythe cell model ([PRL 106, 165702](https://doi.org/10.1103/PhysRevLett.106.165702)) | 2D grid of L² cells, N spins each | 2 | Noise-induced Ising-like → voter-like coarsening | MC + Langevin | PARTS (via Azizi–Pleimling) | M / none |
| Azizi–Pleimling cell model ([2020](https://arxiv.org/abs/2007.10933)) | 2D cells | 2 | Full Al Hammal diagram; GV line incl. **b = −8** | MC | PARTS | **H / counter-evidence** |
| [Rodrigues et al. 2015](https://arxiv.org/abs/1501.01116) | 2D triangular and 3D cubic, z = 6, continuous time | 2 | Ising + DP lines ending at the voter point; "discontinuous" P–A line (jump equals voter jump) | Mean field + MC | FULL | **H / partial** |
| [Molofsky et al. 1999](https://doi.org/10.1006/tpbi.1998.1404) PCA | 2D, 5-cell neighbourhood, synchronous | 2 | Clustering / ergodic; no bistable regime | MC | project FULL | H / counter-evidence |
| Four-species cyclic predator–prey with vacancies ([Szabó & Sznaider 2004](https://arxiv.org/abs/q-bio/0310017)) | 2D square | 2 symmetric alliance families (continuous sets) | Symmetry breaking "equivalent to … spreading models with two equivalent absorbing states" | MC | ABSTRACT | M / none |
| **Four-species cyclic predator–prey with neutral-pair exchange ([Guisoni, Loscar & Girardi 2013](https://doi.org/10.1103/PhysRevE.88.022133))** | Lattice (2D probable, unverified) | **"two symmetric absorbing phases, composed by neutral pairs"** | **"two lines of irreversible first-order phase transition between the two absorbing phases and the active phase"**; triple point; DP endpoints | Gradient method + "standard simulations" | ABSTRACT | **H / partial–moderate** |
| AB monomer–monomer reaction ([Evans & Ray 1993](https://doi.org/10.1103/PhysRevE.47.1018)) | 2D, no diffusion | 2 poisoned states | Voter-like slow poisoning | MC + analysis | ABSTRACT | L–M / none |
| Monomer–monomer surface reaction ([Fichthorn, Gulari & Ziff 1989](https://doi.org/10.1103/PhysRevLett.63.1527)) | 2D lattice MC | Poisoned states (details unverified) | "noise-induced transition to bistability" of a deterministically monostable model | MC | ABSTRACT | L–M / none |
| Bacterial-competition lattice model ([Shimaya & Takeuchi 2019](https://arxiv.org/abs/1804.09895)) | 3D channel (2D cross-section) | 2 (no mutation) | "consistent with the generalized voter class … Ising … not ruled out" | MC | ABSTRACT + conclusion | L–M / none |
| Linear Glauber = noisy voter ([de Oliveira 2003](https://doi.org/10.1103/PhysRevE.67.066101)) | d-dimensional hypercubic | 2 at zero noise | Exact: β = 0, γ = 1, ν = 1/2 for d ≥ 2, log corrections at d_c = 2 | Exact | ABSTRACT | L–M / none |
| Noisy nonlinear voter models ([Llabrés et al. 2026](https://arxiv.org/abs/2505.11358)) | Mean field, 2D lattice, networks | 2 without noise; none with noise | Noise-free: GV; noisy: Ising + discontinuous "Modified GV" + tricritical point | Mean field + FSS | PARTS | M / none (discontinuity only without absorbing states) |

#### 2D details beyond the table
- **Dornic, Chaté, Chave & Hinrichsen 2001.** "Critical coarsening without surface tension: the universality class of the voter model", PRL 87, 045701, DOI 10.1103/PhysRevLett.87.045701, arXiv:cond-mat/0101202. FULL.
  - Abstract: "We argue that voter-like growth is generically observed at order-disorder nonequilibrium transitions solely driven by interfacial noise between dynamically symmetric absorbing states."
  - Rule space: "the 4-parameter space of the general kinetic Ising model above without bulk noise (space defined by r4 = 0)". This is a z = 4 square-lattice family with flip probabilities r_{s,h}.
  - Exponents: "a conservative interpretation of the data is that β = 0, γ = 1, and ν = 1/2, with appropriate logarithmic corrections".
  - On dimension: "d = 2, which appears as the upper critical dimension where voter-like phase-ordering can be induced by the sole presence of interfacial noise." — [arXiv](https://arxiv.org/abs/cond-mat/0101202)
  - The conjecture is stated only for **critical** rules and only for z = 4 rule space, so it does not logically exclude first-order rules.
- **Al Hammal, Chaté, Dornic & Muñoz 2005.** "Langevin description of critical phenomena with two symmetric absorbing states", PRL 94, 230601, arXiv:cond-mat/0411693. FULL.
  - Equation: ∂φ/∂t = (aφ − bφ³)(1 − φ²) + D∇²φ + σ√(1 − φ²) η.
  - Mean field: "b ≤ 0: unique GV transition … For b < 0, the φ⁴ term is present in the potential, but it is not stabilizing … we expect to have a situation similar to the one for b = 0, i.e. a unique GV transition."
  - 2D numerics: b* ≈ 0.50. For b > b*, Ising (Binder U* ≈ 0.61) then DP. For b < b*, GV, with ρ continuous (quoted in KQ1).
  - 3D: "fully first-order" (preliminary). Relevance **High** (defines "fully first-order"); threat none. — [arXiv](https://arxiv.org/abs/cond-mat/0411693)
- **Hinrichsen 1997.** "Stochastic lattice models with several absorbing states", PRE 55, 219, arXiv:cond-mat/9608065. PARTS.
  - The 2D two-absorbing-state model on 80×80: "λc = 0.99(1), δ = 0.9(1), η = 0.00(5), and z = 1.0(1). These results agree roughly with the mean-field exponents … Therefore, we conjecture that 1 < dc ≤ 2 is the upper critical dimension of systems with two absorbing states."
  - Geometry: "active sites in two-dimensional models with two absorbing states are arranged in fractal 'lines' along the boundaries of inactive islands." — [arXiv](https://arxiv.org/abs/cond-mat/9608065)
- **Lee & Vojta 2011.** "Generalized contact process with two symmetric absorbing states in two dimensions", PRE 83, 011114, arXiv:1010.3298. PARTS.
  - Abstract: "The critical behavior of the generic transition is compatible with the generalized voter universality class … In contrast, the transition at zero domain-boundary activation rate is not critical."
  - Body: the σ = 0 line "is thus not a true critical line. It only appears critical because the stationary density ρ_st (trivially) vanishes with σ". So it is **not** a first-order jump. — [arXiv](https://arxiv.org/abs/1010.3298)
  - 1D companion paper: Lee & Vojta, PRE 81, 061128 (2010), arXiv:1003.1308 (ABSTRACT). — [DOI](https://doi.org/10.1103/PhysRevE.81.061128)
- **Barghathi & Vojta.** PRL 109, 170603 (2012) and PRE 93, 022120 (2016), arXiv:1511.04725. ABSTRACTS.
  - "random-field disorder fails to destroy the nonequilibrium phase transition of the one- and two-dimensional generalized contact process … In the two-dimensional case, we map the dynamics onto that of the well studied low-temperature random-field Ising model." Relevance Low. — [PRE 93, 022120](https://doi.org/10.1103/PhysRevE.93.022120); [PRL](https://doi.org/10.1103/PhysRevLett.109.170603)
- **Park 2012** (PRE 85, 041140, arXiv:1202.2201), ABSTRACT.
  - "a model of two-dimensional interacting monomers which has two symmetric absorbing states and exhibits two kinds of phase transition".
  - Conclusion: "the order-disorder transition in two-dimensional models with two symmetric absorbing states is of the Ising universality class, contrary to the recent claim [Nam et al.]", and "the Binder cumulant could be a misleading guide to the critical point in these systems". — [arXiv](https://arxiv.org/abs/1202.2201)
- **Nam, Park, Kim & Lee 2011** (J. Stat. Mech. L06001, arXiv:1103.1222), title "Double transitions, non-Ising criticality and the critical absorbing phase in an interacting monomer–dimer model on a square lattice". Title and Crossref metadata only. — [DOI](https://doi.org/10.1088/1742-5468/2011/06/L06001)
- **Nam, Kim & Lee 2014** (J. Stat. Mech. P08011), ABSTRACT.
  - "two nearby continuous transitions: the Z2 symmetry-breaking order-disorder transition and the absorbing transition with directed percolation criticality".
  - Non-Ising values β ≃ 0.149(2), η ≃ 0.30(1), ν ≃ 1.37(2), and 1/z ≃ 0.466(5). — [DOI](https://doi.org/10.1088/1742-5468/2014/08/P08011)
- **Park 2015** (J. Stat. Mech. P10009, arXiv:1507.05537), ABSTRACT: the critical relaxation agrees with Ising, "unlike the recent claim by Nam et al". — [arXiv](https://arxiv.org/abs/1507.05537)
- **Vazquez & López 2008.** "Systems with two symmetric absorbing states: relating the microscopic dynamics with the macroscopic behavior", PRE 78, 061127, arXiv:0810.1436. ABSTRACT plus skim; the full text is local in campaign/theory/lit.
  - "Starting from the microscopic dynamics on a square lattice, we derive a Langevin equation … the macroscopic behavior only depends on the first derivatives of the spin-flip probabilities … the three types of transitions commonly observed in these systems—generalized voter, Ising and directed percolation."
  - I found no first-order or discontinuous discussion; grep found only "first order in Δφ" expansions. — [arXiv](https://arxiv.org/abs/0810.1436)
- **Castellano & Pastor-Satorras 2012.** PRE 86, 051123, arXiv:1207.4333. PARTS.
  - "three different spin models belonging to the generalized voter class for ordering dynamics in two dimensions … for nonlinear voter models the exit probability … assumes a nontrivial shape".
  - GV point of their nonlinear voter (NLV) model: ε_c = 0.3996(4). — [arXiv](https://arxiv.org/abs/1207.4333)
- **Szabó & Sznaider 2004.** "Phase transition and selection in a four-species cyclic predator-prey model", PRE 69, 031911, arXiv:q-bio/0310017. ABSTRACT.
  - Square lattice with vacancies. "Above the threshold, a symmetry breaking ordering occurs via growing domains containing only two neutral species inside … the observed phase transition seems to be equivalent to those found in spreading models with two equivalent absorbing states although the present model has continuous sets of absorbing states with different portions of the two neutral species." — [arXiv](https://arxiv.org/abs/q-bio/0310017)
- **Guisoni, Loscar & Girardi 2013.** PRE 88, 022133 (Aug 2013). Authors from the repository record: Nara Cristina Guisoni, Ernesto Selim Loscar, Mauricio Girardi.
  - Its 10 OpenAlex citers (Pleimling-group cyclic-game papers and others) do not follow up the first-order claim. — [DOI](https://doi.org/10.1103/PhysRevE.88.022133); [record](https://ri.conicet.gov.ar/handle/11336/7499)
- **Krause, Böttcher & Bornholdt 2012** (PRE 85, 031126, arXiv:1212.4782), ABSTRACT: "a kinetic Ising model with suppressed bulk noise, which is a prominent representative of the generalized voter model phase transition … we characterize the abrupt phase transition. The system shows nonequilibrium dynamics in the presence of absorbing states." — [DOI](https://doi.org/10.1103/PhysRevE.85.031126)
- **Evans & Ray 1993** (PRE 47, 1018), ABSTRACT: "The two-dimensional monomer-monomer (AB) surface reaction model … For equal reactant adsorption rates … simulations reveal the same form of slow poisoning … the same as that found in simulations of the two-dimensional voter model." — [DOI](https://doi.org/10.1103/PhysRevE.47.1018)
- **Fichthorn, Gulari & Ziff 1989** (PRL 63, 1527), ABSTRACT: "a noise-induced transition to bistability in a Monte Carlo simulation of a model heterogeneous catalytic chemical reaction which is deterministically monostable … a consequence of the interaction of noise with the spatial degrees of freedom." — [DOI](https://doi.org/10.1103/PhysRevLett.63.1527)
- **Oliveira & Dickman 2011** (PRE 84, 011125), ABSTRACT. This is **not** Z2-absorbing: there is a single vacuum plus sublattice ordering.
  - "Monte Carlo simulations in two dimensions verify these predictions qualitatively, except for a first-neighbor creation rate of zero. (In the latter case one of the phase transitions is discontinuous.)" Relevance Low. — [DOI](https://doi.org/10.1103/PhysRevE.84.011125)
- **Other 2D items.**
  - Bordogna & Albano 2011 (PRE 83, 046111): interfaces in a 2D generalized voter model, including "between an ordered phase (consensus) and a disordered one (δ = 1/2)". ABSTRACT. Low. — [DOI](https://doi.org/10.1103/PhysRevE.83.046111)
  - Martínez-García, Vázquez, López & Muñoz 2012 (PRE 85, 051125): temporal Griffiths phases in the Ising and GV classes. ABSTRACT. Low. — [DOI](https://doi.org/10.1103/PhysRevE.85.051125)
  - Kumar & Dasgupta 2020 (PRE 102, 052111): a 2D Ising model with modified Metropolis/Glauber rules that break detailed balance. The transition is continuous with 2D Ising exponents, and the abstract does not mention absorbing states. ABSTRACT. Low. — [DOI](https://doi.org/10.1103/PhysRevE.102.052111)
  - Ai, Levine & Chaikin, arXiv:2610.03936 (2 Oct 2026), "An Absorbing State Model with an Emergent Ising Transition". ABSTRACT summary. Repulsive needles on a square lattice give an absorbing-state transition with infinitely many absorbing configurations; a checkerboard state melts with 2D Ising exponents. Low. — [arXiv](https://arxiv.org/abs/2610.03936)

#### 1D context (all "two symmetric absorbing states" or related)
- Grassberger, Krause & von der Twer, "A new type of kinetic critical phenomenon", J. Phys. A 17, L105 (1984), DOI 10.1088/0305-4470/17/3/003. ABSTRACT.
  - 1D cellular automaton with kinks between ordered states; n ∼ t^−α with α = 0.27 ± 0.08.
  - "Possible extensions to two or more dimensions … are discussed." — [DOI](https://doi.org/10.1088/0305-4470/17/3/003)
- Menyhárd & Ódor, J. Phys. A 29, 7739 (1996): 1D non-equilibrium kinetic Ising models (NEKIM), PC class. Citation verified only through the reference lists of DFL 2003 and Al Hammal 2005; **not read**. — [ref. list in arXiv:cond-mat/0411693](https://arxiv.org/abs/cond-mat/0411693)
  - No Ódor- or Menyhárd-authored 2D Z2-absorbing simulation turned up in any sweep. Ódor's RMP reviews 2D Z2 spin rules (citing Drouffe–Godrèche and Dornic) as Ising/voter-class. — [arXiv:cond-mat/0205644](https://arxiv.org/abs/cond-mat/0205644)
- Hwang, Kwon, Park & Park, "Critical phenomena of nonequilibrium dynamical systems with two absorbing states", PRE 57, 6438 (1998), arXiv:cond-mat/9712259. ABSTRACT. 1D interacting monomer–dimer models, PCAs and NEKIM.
  - "a continuous phase transition … which belongs to the universality class of the models with the parity conservation. However, when we break the symmetry … the system goes back to the conventional directed percolation universality class." — [arXiv](https://arxiv.org/abs/cond-mat/9712259)
- Kwon, Hwang & Park, PRE 59, 4949 (1999). ABSTRACT. In 1D, with a symmetry-breaking field, "the spreading velocity of this driven interface shows a discontinuous jump at criticality", explained by "a finite relaxation time in the absorbing phase". Low–Medium, as a front-velocity analogue. — [DOI](https://doi.org/10.1103/PhysRevE.59.4949)
- Bagnoli, Boccara & Rechtman, PRE 63, 046116 (2001). ABSTRACT. 1D PCA: "two second-order and one first-order transition lines that meet at a bicritical point". The first-order line lies between the two absorbing states. — [DOI](https://doi.org/10.1103/PhysRevE.63.046116)
- Bagnoli, Franci & Rechtman, PRE 71, 046108 (2005), arXiv:cond-mat/0405604. ABSTRACT. "a long-range one-dimensional cellular automaton with two symmetric absorbing states … extremely rich phase diagram". — [arXiv](https://arxiv.org/abs/cond-mat/0405604)
- Monetti, PRE 58, 144 (1998). ABSTRACT. A 1D monomer–dimer model with lateral interactions: "first-order irreversible phase transition between two absorbing states". — [DOI](https://doi.org/10.1103/PhysRevE.58.144)
- Dantas, Ticona & Stilck, Braz. J. Phys. 35, 536 (2005), arXiv:cond-mat/0311152. ABSTRACT. In 1D, "the critical line ends at a point where the transition is discontinuous", i.e. the voter point. — [DOI](https://doi.org/10.1590/s0103-97332005000300028)
- Lee & Vojta 2010: 1D generalized contact process, with a special point that "shares some characteristics with a multicritical point". — [DOI](https://doi.org/10.1103/PhysRevE.81.061128)
- Park 2020, "Order-parameter critical exponent of absorbing phase transitions in one-dimensional systems with two symmetric absorbing states", PRE 101, 052114, arXiv:2002.03143. Title only. — [arXiv](https://arxiv.org/abs/2002.03143)
- Kockelkoren & Chaté, "Absorbing phase transitions of branching-annihilating random walks", PRL 90, 125701 (2003), arXiv:cond-mat/0208497. ABSTRACT. A **1D** classification of mA → (m+k)A, nA → (n−l)A into "Four universality classes"; nothing on 2D Z2 models. Low. — [arXiv](https://arxiv.org/abs/cond-mat/0208497)
- Ha, Huse & Samajdar, "Absorbing state transitions with discrete symmetries", arXiv:2502.08702 (v2, Jan 2026). PARTS.
  - 1D only: "models with multiple absorbing states related by discrete symmetries, such as Z2 … and Z3 or S3".
  - A nonlocal variant's surface transition "exhibits a first-order-like behavior". Low for 2D. — [arXiv](https://arxiv.org/abs/2502.08702)
- Hinrichsen's review states that DP2 (directed Ising) and PC "coincide in one spatial dimension". — [Hinrichsen, Adv. Phys. 49, 815 (2000), arXiv:cond-mat/0001070](https://arxiv.org/abs/cond-mat/0001070)

#### 3D context
- LD 2002: 3D q = 2 was first called "continuous … β = 1". DFL 2003 then found separate transitions (L = 60). — [arXiv:cond-mat/0112130](https://arxiv.org/abs/cond-mat/0112130); [arXiv:cond-mat/0301291](https://arxiv.org/abs/cond-mat/0301291)
- Rodrigues et al. 2015, cubic lattice: Ising + DP lines and a "discontinuous" P–A line ending at the voter point, with the same topology as mean field. — [arXiv:1501.01116](https://arxiv.org/abs/1501.01116)
- Al Hammal et al. 2005: 3D Langevin with b ≲ 0 is "fully first-order" (preliminary). — [arXiv:cond-mat/0411693](https://arxiv.org/abs/cond-mat/0411693)
- Lipowski et al. 2017: the 3D voter model has an "active steady state" that Ising dopants change into Ising-type coarsening. — [arXiv:1706.03575](https://arxiv.org/abs/1706.03575)

### Inferences
- **The 2D landscape is binary:** GV, or Ising + DP. Every route to a split found so far is continuous: longer range, cell models with b > b*, or the interacting monomer / monomer–dimer models.
- **Interface density as a discriminator.** No model with two homogeneous absorbing states reports a non-zero ρ at the transition from the active side. Al Hammal 2005 states the opposite for 2D (ρ → 0 continuously). D&N's ρ = 0.23–0.29 down to q* is therefore the cleanest discriminating observable.
- **Guisoni 2013 differs in structure.** The four-species models have (i) two alliance phases, each a continuous family of absorbing configurations, and (ii) dynamics inside the absorbing phase (vacancy hopping or neutral-pair exchange).
  - The DP endpoints and the triple point in Guisoni et al. look more like generic "competing phases" (ZGB-like) first-order physics than like a Z2 Langevin system.
  - This is a plausible argument for why it does not contradict Al Hammal in 2D, and why D&N is still a different, cleaner case.
- **Rodrigues et al. are the only microscopic 2D Z2 model besides D&N whose P–A line is labelled discontinuous** (their criterion and lattice sizes are in the Gaps).
- **The two published readings of it should not be conflated.** Azizi & Pleimling's characterization ("a line of discontinuous phase transitions is replacing the generalized voter critical line") is stronger than Rodrigues et al.'s own wording ("jump … equals the jump observed in the voter model").

### Gaps
- Rodrigues et al.'s maximum lattice size (L ≤ 60, per the project notes) and their exact discontinuity criterion were not re-verified in this pass.
- Sanservino, López, Albano & Monetti, "Study of a monomer-monomer reaction system with short-range repulsive interactions", Eur. Phys. J. B 40, 305–312 (2004), DOI 10.1140/epjb/e2004-00271-1.
  - It cites Dornic 2001, DFL 2003 and Hinrichsen 1997, so it is probably a 2D two-absorbing-state model.
  - **The abstract could not be retrieved**: the Springer page redirected to a login, and Semantic Scholar/OpenAlex have no abstract. Its transition type is unknown, and it should be checked. — [DOI](https://doi.org/10.1140/epjb/e2004-00271-1)
- Fichthorn et al. 1989: whether the model's poisoned states are strictly absorbing (desorption) was not verified.
- Russell & Blythe 2011 was read only through grep and Azizi–Pleimling's description.

---

## Key question 5 — q ≥ 3 symmetric absorbing states in 2D: are discontinuous transitions reported?

### Takeaway
**Yes. A 2D discontinuous transition with q ≥ 3 symmetric absorbing states is established prior art.**
- The non-equilibrium Potts model at q = 3 and 4 with z = 4 shows it (LD 2002; Azizi et al. 2018, with a jump in both M and ρ_a).
- Naming-game models on the 2D square lattice show it, with infinitely many label-symmetric absorbing consensus states.

So "first order with q ≥ 3 symmetric absorbing states in 2D" is not a novel direction per se. What would be new is a q ≥ 3 D&N-like rule with front, nucleation and droplet characterization, or a demonstration that extended range does **not** split it. In the LD model, z = 12 splits q = 3 into Potts(3) + DP.

Two items suggested as 2D q = 3 examples are not:
- Bassler & Browne's three-species monomer–monomer model is **1D**.
- Most "multi-state q-voter" discontinuous transitions are mean-field or network results **with noise that removes the absorbing states**.

### Cited Findings
- **LD 2002.** 2D q = 3: "discontinuous phase transition"; q = 4: "the jump is larger and the first-order character of the transition is even more transparent. We expect that such a behaviour persists for all q > 3." — [arXiv:cond-mat/0112130](https://arxiv.org/abs/cond-mat/0112130)
- **Droz & Lipowski 2003 review:** "d = 2, q = 3 … clear evidence of the discontinuous phase transition". — [arXiv:cond-mat/0308419](https://arxiv.org/abs/cond-mat/0308419)
- **Azizi, Stidham & Pleimling 2018.**
  - q = 3, z = 4: "Both order parameters display a jump at the same temperature".
  - q = 3, z = 12: "replaced by two continuous phase transitions": 2D q = 3 Potts + DP. — [arXiv:1809.08366](https://arxiv.org/abs/1809.08366)
- **Crokidakis & Brigatti** (J. Stat. Mech. P01019, 2015, arXiv:1412.2994), ABSTRACT.
  - "we study on a 2-dimensional square lattice a recent version of the Naming Game … The model displays a nonequilibrium phase transition at a critical point σc ≈ 25.6, which separates an absorbing consensus state from an active fragmented state … The finite-size scaling analysis of our simulations strongly suggests that the phase transition is discontinuous." — [arXiv](https://arxiv.org/abs/1412.2994)
- **Brigatti & Hernández** (PRE 94, 052308, 2016, arXiv:1609.02869), ABSTRACT + model section.
  - "a non-equilibrium phase transition between an absorbing ordered consensus state, which occurs for small noise, and a disordered phase … We observe a scaling behavior typical of a discontinuous transition."
  - Setting: "simulated on a regular 2D square lattice … agents communicate only with their four nearest neighbors." — [arXiv](https://arxiv.org/abs/1609.02869)
- **Neto & Brigatti** (PRE 101, 022112, 2020, arXiv:1911.07294), ABSTRACT: in the 2D Naming Game "a nonequilibrium low-dimensional system with different absorbing states … the spatial quenched disorder does not destroy the discontinuous transition." — [arXiv](https://arxiv.org/abs/1911.07294)
- **Bassler & Browne 1996** (PRL 77, 4094, arXiv:cond-mat/9609191), ABSTRACT + text.
  - "a reactive steady state bordered by three equivalent unreactive phases … The transition from the reactive to a saturated phase shows directed percolation critical behavior. Each pair of these reactive-saturated phase boundaries join at a bicritical point where the universal behavior is in the even branching annihilating random walk class."
  - The text says one species "saturates the chain", so it is 1D.
  - The 1998 follow-up (J. Phys. A 31, 6309) states "We study the one dimensional three species monomer-monomer reaction model" and "the transitions between adsorbing phases are first-order". — [PRL arXiv](https://arxiv.org/abs/cond-mat/9609191); [1998 DOI](https://doi.org/10.1088/0305-4470/31/30/001)
- **1D q = 3 is dynamics-dependent** (LD 2002b: generically active, or PC under modified dynamics). Ha–Huse–Samajdar (1D Z3/S3) is likewise 1D only. — [PRE 66, 016106](https://doi.org/10.1103/PhysRevE.66.016106); [arXiv:2502.08702](https://arxiv.org/abs/2502.08702)
- **Noisy multi-state models (absorbing states removed):**
  - 2D S-state block voter with noise: "continuous transition for S ≤ 4 to a discontinuous transition for S ≥ 5" (OpenAlex abstract snippet). — [Araújo, Filho & Moreira, arXiv:1609.08573](https://arxiv.org/abs/1609.08573)
  - Mean-field multi-state noisy q-voter: discontinuous (project notes). — [Nowak, Stoń & Sznajd-Weron, Sci. Rep. 11, 6098 (2021)](https://doi.org/10.1038/s41598-021-85361-9)
- **Szabó & Szolnoki 2002** (three-state cyclic voter with Potts energy, 2D, PRE 65, 036115): studies pattern formation (loops vs spirals). No absorbing-transition order claim. ABSTRACT. — [DOI](https://doi.org/10.1103/PhysRevE.65.036115)

### Inferences
- **For the "richer direction" (q ≥ 3):** the baseline expectation in 2D is already first order for short range and monotone Metropolis rules. A q ≥ 3 D&N generalization would be incremental unless it shows something qualitatively new. Candidates:
  - first order that survives larger range (contrast Azizi 2018's z = 12 split);
  - an interplay between Z2 and Z3 sectors;
  - wetting or three-phase-contact phenomena with an active slab.
- **Naming-game models are a different symmetry class:** infinitely many labels and memory or inventory variables. They are useful as "2D first-order into symmetric absorbing states" precedents but are not Potts-like.

### Gaps
- No 2D q = 3 study with D&N-type non-monotone rules was found.
- The ordering of the Lipowski–Droz (q, d) diagram for q ≥ 3 at z = 8 is untested.

---

## Key question 6 — Forward-citation sweep (Dornic 2001, DFL 2003, LD 2002, Rodrigues 2015, Azizi–Pleimling 2020, and other seeds)

### Takeaway
I swept about 410 unique Semantic Scholar citing records over 12 seeds and about 400 OpenAlex citing records over 12 seeds. **No citing paper reports a fully first-order transition between a disordered active phase and two Z2 homogeneous absorbing states in 2D.**

The sweep did surface the relevant 2D models that were missing from the project notes:
- Lee & Vojta 2011;
- Park 2012 / Nam et al. 2011, 2014 / Park 2015;
- Krause et al. 2012;
- Guisoni et al. 2013 (found via keyword sweep, not citations);
- Szabó & Sznaider 2004;
- Neto & Brigatti 2020;
- Azizi et al. 2018;
- Barghathi & Vojta 2012/2016;
- Martínez-García et al. 2012;
- Shimaya & Takeuchi 2019.

### Cited Findings
- **Semantic Scholar citing-paper counts** (Graph API), 409 unique citers overall: Dornic 2001: 286; Al Hammal 2005: 84; Kockelkoren & Chaté 2003: 48; Hinrichsen 1997: 41; Bassler & Browne 1996: 29; DFL 2003: 27; Hwang 1998: 18; Russell & Blythe 2011: 12; Azizi & Pleimling 2020: 6 (including off-topic records); Rodrigues 2015: 4; LD 2002: 0 (coverage failure); Lipowski 2000: 0. — [Semantic Scholar API](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1103/PhysRevE.67.056108)
- **OpenAlex cited-by counts:** LD 2002: 17; Lipowski 2000: 9; DFL 2003: 32; Rodrigues 2015: 5; Azizi–Pleimling 2020: 5; Castellano–PS 2012: 13; Vazquez–López 2008: 72; Russell–Blythe 2011: 18; Hinrichsen 1997: 107; Azizi 2018: 1; Al Hammal 2005: 113; Guisoni 2013: 10. — [OpenAlex API](https://api.openalex.org/works/doi:10.1103/PhysRevE.65.056114)
- **Citers of Rodrigues 2015:**
  - Azizi et al. 2018, Azizi & Pleimling 2020, and Lipowski et al. 2023 (random graphs).
  - Fernandes, Silva & Bernardi 2018, a ZGB model with CO desorption (single absorbing state).
  - A 2019 São Paulo thesis (Wada).
  - None extends the 2D "discontinuous" line. — [OpenAlex](https://api.openalex.org/works/doi:10.1088/1742-5468/2015/01/P01035)
- **Citers of Azizi–Pleimling 2020:** Lipowski et al. 2023 plus unrelated records (nanofluids, percolation polynomials). No follow-up on Z2 first order. — [OpenAlex](https://api.openalex.org/works/doi:10.1103/PhysRevE.102.022112)
- **Citers of LD 2002** (OpenAlex) include:
  - Neto & Brigatti 2020 (2D naming game, discontinuous; see KQ5);
  - Kumar & Dasgupta 2020 (2D, continuous Ising);
  - Martínez-García et al. 2012;
  - Lipowska & Lipowski 2014 (naming game);
  - Syiem & Kharwanlang 2023 (1D three-state clock model with an absorbing state, hysteresis);
  - Noh & Park 2004 (1D).
  - None is a 2D Z2 first-order report. — [OpenAlex](https://api.openalex.org/works/doi:10.1103/PhysRevE.65.056114)
- **Citers of DFL 2003** include Lee & Vojta 2011, Park 2012, Oliveira & Dickman 2011, Barghathi & Vojta 2016, Castellano 2009, Castellano–PS 2012, Azizi 2018/2020, Rodrigues 2015, Llabrés 2026, Khaleque & Sen 2015 (Ising on networks), Timonin & Chitov 2014 (1D kinetic contact process with parallel update; DP), and Szabó & Sznaider 2004. None is first order in 2D with Z2. — [OpenAlex](https://api.openalex.org/works/doi:10.1103/PhysRevE.67.056108)
- **Citers of Dornic 2001** (286 in Semantic Scholar) are dominated by opinion-dynamics and network voter papers. Physics-relevant non-network items are already covered above.
  - One borderline item: "Critical discontinuous phase transition in the threshold contact process" (Silva & Oliveira 2011). This is a single-absorbing-state contact process, not Z2. Title only. — [Semantic Scholar API](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1103/PhysRevLett.87.045701/citations)
- **Recent citers (2023–2026)** are 1D or not Z2 absorbing:
  - Llabrés, San Miguel & Toral, PRR 8, 013015 (2026; Crossref issued 9 Jan 2026), arXiv:2505.11358.
  - Villegas & Meloni, arXiv:2610.18506 (2026): "a directed-percolation line, a tricritical point, and a first-order region whose spinodal is the explosive transition". Single absorbing state, mean-field/network. Methodologically notable: "In finite dimensions the location of a first-order transition is set by the stalling of fronts", which is the project's criterion.
  - Saif 2023 (arXiv:2308.07196, 1D).
  - Burridge 2025.
  - Ha–Huse–Samajdar 2025/26 (1D). — [Llabrés DOI](https://doi.org/10.1103/6xdk-2cjy); [Villegas–Meloni](https://arxiv.org/abs/2610.18506); [Saif](https://arxiv.org/abs/2308.07196)

### Inferences
- The citation graph around the GV and Langevin papers is well populated with GV and split results and is silent on 2D fully first-order Z2 transitions. This supports (but cannot prove) the novelty of the narrowed claim.
- Guisoni 2013 sits outside this citation cluster: it was not found by citation, only by keyword. Searches limited to the GV citation graph would miss it, so it deserves explicit discussion.

### Gaps
- Google Scholar "cited by" (usually larger counts, e.g. for Dornic 2001) was not accessed.
- arXiv full-text search was blocked (HTTP 429), so arXiv-only preprints that never got a DOI and use non-standard wording could have been missed.
- Semantic Scholar returned 0 citers for LD 2002 and Lipowski 2000. Only OpenAlex was used for those.

---

## Key question 7 — Keyword and research-group sweep (first-order / tricritical / Z2-absorbing; Brazilian, Polish, Hungarian, French, Spanish groups)

### Takeaway
The keyword sweep found 74 physics records whose abstracts combine "absorbing/voter" with "first-order/discontinuous/tricritical/bistable". Almost all are:
- (i) mean-field or network q-voter models with noise (no absorbing states);
- (ii) 1D models;
- (iii) single-absorbing-state 2D models;
- (iv) GV-point terminology.

The only 2D lattice items claiming first-order active–absorbing transitions with symmetric absorbing phases are **Guisoni et al. 2013** (four-species cyclic, alliance phases) and the **naming-game** models (infinitely many absorbing states).

By group, no Brazilian, Polish, Hungarian, French or Spanish work was found that reports a 2D fully first-order Z2 transition. Per-group details are under Cited Findings.

### Cited Findings
- **Sweep method:** OpenAlex title/abstract search, 20 queries, 707 records. Queries included "symmetric absorbing" (53 hits), "two absorbing states" (96), "double absorbing" (39), "generalized voter" (36), voter + "first-order" (186), voter + tricritical (5), "absorbing states" + Potts (8), "multiple absorbing states" (30), "three absorbing states" (9), and "interfacial noise" (18). — [OpenAlex search API](https://api.openalex.org/works?filter=title_and_abstract.search:%22symmetric%20absorbing%22)
- **"Tricritical" + voter.** Neither item below has absorbing states in its first-order regime, and neither is a 2D Z2-absorbing first-order result.
  - Llabrés et al. 2026: "a continuous Ising transition and a discontinuous transition we call Modified Generalized Voter (MGV). These two transition lines converge at a tricritical point". The noise "eliminates the absorbing states". 2D-lattice FSS is used for the Ising line; the noise-free models with symmetric absorbing states "exhibit a Generalized Voter (GV) transition". — [arXiv:2505.11358](https://arxiv.org/abs/2505.11358)
  - Kim, Park & Kahng, arXiv:2610.16053 ("A Shy Voter Model"): the abstract snippet mentions a tricritical point; the lattice and absorbing structure were **not verified**. — [arXiv](https://arxiv.org/abs/2610.16053)
- **Brazilian groups.**
  - de Oliveira and Tomé: OMS 1993, linear Glauber 2003, Rodrigues et al. 2015.
  - Fiore, two strands, neither a Z2 two-absorbing-state transition in 2D:
    - Discontinuous absorbing transitions in "two simple lattice models originally presenting discontinuous absorbing transitions" (Pianegonda & Fiore, J. Stat. Mech. (2015) P08018; ABSTRACT; the number and symmetry of absorbing states were not checked). — [DOI](https://doi.org/10.1088/1742-5468/2015/08/P08018)
    - Discontinuous Z2 order–disorder transitions without absorbing states, in the "inertial majority vote model" under temporal disorder (Encinas & Fiore, PRE 103, 032124 (2021); mean field and random-regular graphs; ABSTRACT). — [DOI](https://doi.org/10.1103/PhysRevE.103.032124)
  - Stilck: Dantas et al. 2005, 1D.
  - Brigatti and Crokidakis: 2D naming game, discontinuous.
  - Dickman: Oliveira & Dickman 2011, single absorbing state.
  - None reports a 2D fully first-order Z2 transition. — [Rodrigues 2015](https://arxiv.org/abs/1501.01116); [de Oliveira 2003](https://doi.org/10.1103/PhysRevE.67.066101); [Crokidakis–Brigatti 2015](https://arxiv.org/abs/1412.2994)
- **Polish groups.**
  - Sznajd-Weron: many discontinuous q-voter results (Nyczka et al. 2012; Abramiuk et al. 2019; Nowak & Sznajd-Weron 2020/2022; Lipiecki & Sznajd-Weron 2025). They are **mean-field or networks, with independence or anticonformity noise that removes the absorbing states**.
  - Lipowski (Poznań): the LD/DFL program covered above.
  - Sources: [Nyczka, Sznajd-Weron & Cisło, PRE 86, 011105 (2012)](https://doi.org/10.1103/PhysRevE.86.011105); [Abramiuk et al., Entropy 21, 521 (2019)](https://doi.org/10.3390/e21050521)
- **Hungarian groups.**
  - Ódor and Menyhárd: 1D NEKIM/PC and the RMP review, with no 2D Z2 first-order work found.
  - Szabó (Budapest): four-species cyclic model, GV-like.
  - Ódor's RMP: "First order phase transitions have rarely been seen in low dimensions." — [arXiv:cond-mat/0205644](https://arxiv.org/abs/cond-mat/0205644); [Szabó–Sznaider](https://arxiv.org/abs/q-bio/0310017)
- **French group (Chaté, Dornic).**
  - Dornic 2001 (GV conjecture); Al Hammal 2005 (Langevin, 3D fully first-order); Kockelkoren & Chaté 2003 (1D only).
  - Benitez, Duclut, Chaté, Delamotte, Dornic & Muñoz 2016, PRL 117, 100601: exact Langevin equations for reaction–diffusion. Azizi–Pleimling cite it as a reaction–diffusion system "described by (1)" whose physical parameter range cannot reach the different lines (ABSTRACT).
  - Sources: [Al Hammal 2005](https://arxiv.org/abs/cond-mat/0411693); [Benitez et al. 2016](https://doi.org/10.1103/PhysRevLett.117.100601); [Azizi–Pleimling](https://arxiv.org/abs/2007.10933)
- **Spanish groups.**
  - Muñoz: Al Hammal 2005; Castellano 2009; Martínez-García 2012; Al Hammal, Bonachela & Muñoz 2006 ("Absorbing state phase transitions with a non-accessible vacuum", J. Stat. Mech. P12007; title only).
  - San Miguel and Toral (IFISC): Llabrés 2023 partisan voter (mean field; continuous and discontinuous noise-induced transitions); Llabrés 2026.
  - Vázquez and López: Vazquez–López 2008.
  - None reports a 2D fully first-order Z2 transition. — [Castellano 2009](https://arxiv.org/abs/0907.1775); [Llabrés 2023 PRE 108, 054106](https://doi.org/10.1103/PhysRevE.108.054106)
- **Other groups found.**
  - Argentina (La Plata/CONICET: Loscar, Guisoni, Girardi; Albano; Monetti): **Guisoni 2013 (first order)**, Bordogna & Albano 2011, Monetti 1998 (1D), and Sanservino et al. 2004 (abstract unavailable).
  - Korea (Park, Kwon, Hwang, Nam, Lee): 1D directed-Ising/PC; 2D splits.
  - USA: Vojta (2D generalized contact process); Pleimling (Azizi).
  - Japan: Shimaya & Takeuchi. — [Guisoni 2013](https://doi.org/10.1103/PhysRevE.88.022133); [Lee–Vojta 2011](https://arxiv.org/abs/1010.3298)

### Inferences
- The only realistic novelty threats are Guisoni et al. 2013 and Rodrigues et al. 2015, plus the withdrawn LD 2002 and GV-terminology confusion. The project's paper should:
  - cite all four (Guisoni 2013, Rodrigues 2015, LD 2002, and the GV-terminology papers) explicitly;
  - define "fully first-order" operationally: ρ jump, front reversal, critical droplets, L-independent lifetimes;
  - state why each prior case does not meet that definition. Guisoni: alliance families with internal dynamics, gradient-method evidence. Rodrigues: magnetization-only criterion, small L. LD 2002: GV log decay.
- Sznajd-Weron and IFISC "discontinuous" voter results are a frequent source of reviewer confusion. Note explicitly that they occur only when noise destroys the absorbing states, and mostly in mean field or on networks.

### Gaps
- Sanservino et al. 2004 (2D monomer–monomer with repulsion) is unresolved.
- Kim, Park & Kahng 2026 (shy voter) is unverified.
- Llabrés et al. 2026's 2D-lattice behaviour of the MGV line was only partly read: it is unclear whether the MGV first-order line is seen on the 2D lattice or only in mean field.
- The Henkel–Hinrichsen–Lübeck book chapter on GV/DP2 was not accessed.
