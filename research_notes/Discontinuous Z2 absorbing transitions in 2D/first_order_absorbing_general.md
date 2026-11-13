# Discontinuous (first-order) absorbing-state transitions in d = 1, 2: no-go arguments, established cases, coexistence, disorder, nucleation/wetting, diagnostics, and what they imply for a 2D Z2 claim

Compiled 2026-10-27. These notes add to `campaign/theory/LITERATURE.md` (they do not repeat Molofsky 1999, Cox–Durrett–Perkins, Cox–Perkins 2023, Al Hammal 2005, Droz 2003, Vazquez–López 2008, Castellano 2009, Rodrigues 2015, Azizi–Pleimling 2020, Llabrés 2026 or Bagnoli–Boccara–Rechtman 2001, which are covered there).
**Read-depth tags:** [FT] full text read; [PT] targeted sections of the full text read (arXiv/OSTI PDF); [ABS] abstract only; [SUM] full text read only through a fetch-tool summary, so the quotes are approximate; [CIT] citation checked against a later paper's reference list, paper not seen. **Evidence:** Strong / Moderate / Weak / Disputed. **Relevance** to the D&N project: High / Medium / Low. All DOIs were checked against Crossref or arXiv metadata unless marked unverified.

## Q1. Hinrichsen's 1+1-D no-go argument: content, scope, counterexamples and extensions, and status in 2D

### Takeaway
Hinrichsen's argument, which was never published (arXiv only), says that a *fluctuating* ordered phase in 1+1 D cannot be stabilised by short-range domain-wall forces. It is accepted as a heuristic. Most later 1D first-order claims in simple short-range single-species models either turned out to be long crossovers (the triplet-creation model) or rely on one of the listed loopholes: long-range interactions, effective long-range forces as in wetting, a large local state space, or conservation. One exception I could not resolve: Ódor–Szolnoki (1996) report 1D first order for their range-4 synchronous CA with y = 4, 5, and I found no re-examination of it. The argument does not extend to 2D. In 2D, curvature-driven shrinkage of islands (line tension) or Toom-like anisotropic erosion supplies the missing mechanism, and first-order absorbing transitions with short-range rules are well established. They are not guaranteed, though: strong noise, limited diffusion or quenched disorder can still turn a mean-field first-order transition continuous in 2D.

### Cited Findings
- **Hinrichsen 2000** — H. Hinrichsen, "First-order transitions in fluctuating 1+1-dimensional nonequilibrium systems", arXiv:cond-mat/0006212 (2000), unpublished (Park 2009 calls it "his unpublished work"). [FT]. Models: triplet creation 3X→4X, X→Ø; annihilation/fission 2X→3X, 2X→Ø; X+Y→2Y, Y→X; all d = 1.
  - Claim: "first-order phase transitions in 1+1-dimensional nonequilibrium systems with fluctuating ordered phases are impossible, provided that there are no additional conservation laws, long-range interactions, macroscopic currents, or special boundary conditions … The apparent first-order behavior turns out to be a transient phenomenon, crossing over to a continuous transition after very long time."
  - Mechanism: a minority island is a pair of domain walls with only an exponentially decaying attraction. The walls' random walk gives a small but finite probability of unbinding, so "the ordered phase is only marginally stable, irrespective of the interaction range r₀".
  - On d ≥ 2: "both mechanisms for the elimination of minority islands – surface tension and anisotropic propagation velocities – can only be implemented in at least two spatial dimensions" and "these arguments are not valid in higher dimensions where the dynamics of domain walls may depend on the local curvature".
  - Proposed test: a stationary minority-island size distribution that decays exponentially and then saturates at a finite constant is "a hallmark for the presence of a short-range force between domain walls".
  - Triplet-creation model at D = 0.9: Dickman–Tomé's spinodals λ₋ = 10.12(1), λ₊ = 10.30(2) bracket Hinrichsen's λc = 10.145(10). Seeds show a CDP-like window (δ = 1/2, η = 0) from about 10³ to 10⁵ steps, then cross over to DP between 10⁵ and 10⁷. Dickman–Tomé "could only go up to 22000 time steps".
  - Evidence: Strong (simulations plus a heuristic argument; not a theorem). Relevance: Medium (it sets the d = 1 baseline and names the 2D mechanisms). [Source](https://arxiv.org/abs/cond-mat/0006212)
- **Dickman & Tomé 1991** — R. Dickman, T. Tomé, "First-order phase transition in a one-dimensional nonequilibrium model", Phys. Rev. A 44, 4833 (1991), doi:10.1103/PhysRevA.44.4833. [ABS]. Triplet creation with hopping, d = 1. Claim: "continuous, and in the directed percolation class, for low diffusion rates, and first order for sufficiently high diffusion rates", while 2X→3X "remains continuous". Evidence: Disputed and since overturned (see the next three entries). Relevance: Medium, as the classic false positive. [Source](https://doi.org/10.1103/PhysRevA.44.4833)
- **Cardozo & Fontanari 2006** — G. O. Cardozo, J. F. Fontanari, "Revisiting the nonequilibrium phase transition of the triplet-creation model", Eur. Phys. J. B 51, 555 (2006), doi:10.1140/epjb/e2006-00246-2, arXiv:cond-mat/0607691. [ABS]. Supported a discontinuous transition at high diffusion but admitted "Strong crossover effects … make it difficult, if not impossible, to establish the exact location of the tricritical point." Evidence: Disputed. Relevance: Low. [Source](https://doi.org/10.1140/epjb/e2006-00246-2)
- **Ódor & Dickman 2009** — G. Ódor, R. Dickman, "On the absorbing-state phase transition in the one-dimensional triplet creation model", J. Stat. Mech. (2009) P08024, doi:10.1088/1742-5468/2009/08/P08024, arXiv:0907.0112. [ABS]. The order parameter "appears to be a discontinuous function … no evidence of a stable interface between active and absorbing phases is found", yet via a mapping to modified CDP "the transition is continuous, despite the seemingly discontinuous phase transition suggested by studies of finite systems". Evidence: Strong. Relevance: Medium (the absence of a stable interface is the tell-tale; contrast D&N, which has stationary fronts). [Source](https://doi.org/10.1088/1742-5468/2009/08/P08024)
- **Park 2009** — S.-C. Park, "Absence of the discontinuous transition in the one-dimensional triplet creation model", Phys. Rev. E 80, 061103 (2009), doi:10.1103/PhysRevE.80.061103, arXiv:0903.3436. [ABS]. Finds the 1D TCM in DP "even for larger diffusion constant than the suggested tricritical point". The phase boundary is "well described by the crossover from the mean field to the directed percolation". Evidence: Strong. Relevance: Medium. [Source](https://doi.org/10.1103/PhysRevE.80.061103)
- **1D counterexamples via the loopholes.**
  - **Ginelli, Hinrichsen, Livi, Mukamel, Politi 2005**, "Directed percolation with long-range interactions: Modeling nonequilibrium wetting", Phys. Rev. E 71, 026121, doi:10.1103/PhysRevE.71.026121, arXiv:cond-mat/0407635 [ABS]. The activation rate at the edge of an inactive island of length l is 1 + a l^(−σ): "for σ>1 the transition is continuous and belongs to … directed percolation, while for 0<σ<1, the transition becomes first order." [Source](https://doi.org/10.1103/PhysRevE.71.026121)
  - **Fiore & de Oliveira 2013**, "Robustness of first-order phase transitions in one-dimensional long-range contact processes", Phys. Rev. E 87, 042101, doi:10.1103/PhysRevE.87.042101, arXiv:1306.3163 [ABS]: "for all examples presenting a single absorbing state, a discontinuous transition is possible for small σ". [Source](https://doi.org/10.1103/PhysRevE.87.042101)
  - **Fiore & Landi 2014**, Phys. Rev. E 90, 032123, doi:10.1103/PhysRevE.90.032123 [ABS]. For intermediate σ, diffusion "causes a suppression of the phase coexistence leading to a continuous transition belonging to the DP universality class". [Source](https://doi.org/10.1103/PhysRevE.90.032123)
  - **Ginelli, Hinrichsen, Livi, Mukamel, Torcini 2006**, "Contact processes with long range interactions", J. Stat. Mech. P08008, doi:10.1088/1742-5468/2006/08/P08008 [ABS]. A different long-range CP (creation rate decaying with distance to the nearest particle) is continuous with varying exponents, so the loophole is model-specific. [Source](https://doi.org/10.1088/1742-5468/2006/08/P08008)
  - **Hinrichsen, Livi, Mukamel, Politi 2000**, "First-order phase transition in a (1+1)-dimensional nonequilibrium wetting process", Phys. Rev. E 61, R1032, doi:10.1103/PhysRevE.61.R1032, arXiv:cond-mat/9906039 [ABS]: "the wet (unbound) and the non-wet (pinned) states coexist and are both thermodynamically stable in a domain of the dynamical parameters". Ginelli et al. 2005 attribute this to effective long-range interactions. [Source](https://doi.org/10.1103/PhysRevE.61.R1032)
  - **Basu, Basu & Mohanty 2010**, "Discontinuous absorbing state transition in (1+1) dimension", arXiv:1010.2007 (journal publication not verified) [ABS]. N valve channels per site; discontinuous "in the large N limit … in all (d+1) dimensions". [Source](https://arxiv.org/abs/1010.2007)
  - Evidence: Moderate to Strong. Relevance: Low to Medium. These are the loopholes, and the D&N model uses none of them.
- **Rigorous side, positive-rates problem.** P. Gács, "Reliable cellular automata with self-organization", J. Stat. Phys. 103, 45 (2001), doi:10.1023/A:1004823720305, arXiv:math/0003117 [ABS]. Keeping a bit forever in a positive-rates PCA "is nontrivial … there is a solution in 2 dimensions". The paper constructs a (hierarchical, extremely complex) 1D nonergodic PCA, so a 1D no-go cannot be a theorem for arbitrary local rules. Companion: L. F. Gray, "A reader's guide to Gács's 'positive rates' paper", J. Stat. Phys. 103, 1 (2001), doi:10.1023/A:1004824203467 [CIT, title only]. Relevance: Low. [Source](https://doi.org/10.1023/A:1004823720305)
- **Status in the reviews.** Ódor, "Universality classes in nonequilibrium lattice systems", Rev. Mod. Phys. 76, 663 (2004), doi:10.1103/RevModPhys.76.663 [PT, Sec. II]:
  - "First order phase transitions have rarely been seen in low dimensions … fluctuation induced second ordered phase transitions are likely to appear."
  - Restates Hinrichsen's hypothesis.
  - "In two dimensions there are certain stochastic cellular automata for which systematic cluster mean-field techniques combined with simulations made it possible to prove first order transitions (Ódor and Szolnoki, 1996) firmly."
  - [Source](https://doi.org/10.1103/RevModPhys.76.663)
- **Lübeck 2006** [PT] on d = 1: "the first-order line and the tricritical point collapse for D = 1 … (this agrees with the more general argument that first-order transitions can not occur in fluctuating one-dimensional systems…)", while "In D ≥ 2 the tricritical point separates a line of DP-like transitions from a line of first-order transitions". [Source](https://doi.org/10.1007/s10955-005-9000-z)
- **Windus & Jensen 2007** [ABS]: 2A→3A, A→Ø with diffusion is "continuous in 1+1 dimensions and of first-order in higher dimensions". [Source](https://doi.org/10.1088/1751-8113/40/10/005)
- **Weissmann & Shnerb 2014** — H. Weissmann, N. M. Shnerb, "Stochastic desertification", EPL 106, 28004 (2014), doi:10.1209/0295-5075/106/28004, arXiv:1401.0160. [ABS+PT]. Their simulations are on a **1d lattice of L = 10000 sites**. A bistable model with an absorbing state has a continuous transition under adiabatic driving, and "the front velocity renormalizes to zero at the extinction transition", leaving a finite region "where domain walls are stable". Evidence: Moderate (1D only). Relevance: Medium, as a warning that front velocities can be pinned at zero over a finite interval. [Source](https://doi.org/10.1209/0295-5075/106/28004)
- **2D does not guarantee first order: Villa Martín, Bonachela, Levin & Muñoz**, "Eluding catastrophic shifts", PNAS 112, E1828 (2015), doi:10.1073/pnas.1414708112. [SUM of the PMC full text].
  - Spatial Langevin model with an Allee (cubic) term, mainly 2D.
  - With baseline noise and D = 1 the 2D transition is truly discontinuous with bistability and hysteresis.
  - Stronger demographic noise (σ² = 4) or limited diffusion (D = 0.1) makes it continuous in the DP class; an RG flow places a tricritical point near D ≈ 0.9.
  - Any quenched spatial heterogeneity smooths it in 2D but not in 3D.
  - In 1D it is always continuous.
  - Diagnostics: decay and spreading experiments and interface ("Maxwell point") runs.
  - Evidence: Moderate to Strong. Relevance: High (fluctuations can kill first order in 2D). [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC4403222/)
- **Sampaio Filho et al. 2018**, "The symbiotic contact process: phase transitions, hysteresis cycles, and bistability", Phys. Rev. E 98, 062108, doi:10.1103/PhysRevE.98.062108 [ABS]. Mean field and random graphs show hysteresis and bistability, but "on a regular square lattice, we find no signs of bistability or hysteretic behavior … the transition in two dimensions is rather continuous". This is a 2D case where mean-field first order fails. Relevance: Medium. [Source](https://doi.org/10.1103/PhysRevE.98.062108)
- **Grassberger, Chen, Ghanbarnejad & Cai 2016**, "Phase transitions in cooperative coinfections: Simulation results for networks and lattices", Phys. Rev. E 93, 042316, doi:10.1103/PhysRevE.93.042316 [ABS]. These are SIR-type, not steady-state absorbing models: "No first order transitions exist on trees and on 2-d lattices with local contacts", but they do exist "on 2-d lattices with sufficiently long-ranged contacts" and for d ≥ 4. Relevance: Low to Medium. [Source](https://doi.org/10.1103/PhysRevE.93.042316)

### Inferences
- The status of the no-go maps onto D&N as follows: 1D no; 2D yes in principle, given a working island-elimination mechanism. D&N's measured R_c·v ≈ const is direct evidence for Hinrichsen's "surface tension" mechanism (1/R shrinkage), which is exactly the ingredient the argument says is missing in 1D. This is the strongest conceptual link between the project and the no-go literature.
- "2D, short-range, so first order is allowed" is not sufficient. The PNAS 2015 and symbiotic-CP results show that strong fluctuations can still renormalise a first-order transition into a continuous one in 2D. Near the voter point (p → 1) fluctuations are largest, so a tricritical point on D&N's q*(p) line is a plausible scenario to rule in or out (compare Lübeck's 2D first-order line ending at finite p).

### Gaps
- Hinrichsen's 2000 paper was never refereed. I found no formal theorem in either direction for 2D absorbing systems with Moore-type rules.
- Only the abstract of Dickman & Tomé 1991 was read.
- Grassberger 1982 (Z. Phys. B 47, 365, doi:10.1007/BF01313803) was not read and no abstract was available. Liu–Wang–Evans 2018 characterise it as finding a *continuous* transition for one realisation of Schlögl's second model.
- Prakash & Nicolis 1997 (J. Stat. Phys. 86, 1289, doi:10.1007/BF02183624) was not read either.
- Ódor–Szolnoki's 1D y = 4, 5 first-order claim has not, to my knowledge, been re-examined against Hinrichsen's argument.

## Q2. Established 2D first-order absorbing transitions (one or several asymmetric absorbing states, and q-symmetric cases) and the evidence used

### Takeaway
2D first-order absorbing transitions are firmly established in surface-reaction models (the ZGB CO-poisoning transition), in "type-2 Schlögl", quadratic or threshold contact processes, in synchronous PCAs (Bidaux–Boccara–Chaté; Ódor–Szolnoki) and in tricritical-DP models. The evidence that convinced the field combines:
- hysteresis and constant-coverage metastability loops;
- interface (stripe) stationarity;
- nucleation and metastable-decay kinetics (KJMA);
- bimodal quasi-stationary distributions with finite-size shifts ∝ 1/V.

For symmetric multi-absorbing-state models in 2D the record is thinner:
- q = 3 (Lipowski–Droz 2002) was judged discontinuous;
- q = 2 (two Z2 absorbing states) was explicitly *inconclusive*;
- other 2D Z2 models (generalised contact process, interacting monomer–dimer) show continuous GV, Ising or DP transitions.

### Cited Findings
#### ZGB family (two asymmetric absorbing states: one DP and one first-order transition)
- **Ziff, Gulari & Barshad 1986**, "Kinetic phase transitions in an irreversible surface-reaction model", Phys. Rev. Lett. 56, 2553, doi:10.1103/PhysRevLett.56.2553. [ABS]
  - "undergo both first- and second-order kinetic phase transitions".
  - Values from Hinrichsen's review [PT]: O-poisoning for y < y₁ ≃ 0.389 (DP) and CO-poisoning for y > y₂ ≃ 0.525 (first order); "In real catalytic reactions, however, only the discontinuous transition at y = y₂ can be observed."
  - Evidence: Strong. Relevance: Medium.
  - [Source](https://doi.org/10.1103/PhysRevLett.56.2553); [review](https://doi.org/10.1080/00018730050198152)
- **Meakin & Scalapino 1987**, "Simple models for heterogeneous catalysis: Phase transition-like behavior in nonequilibrium systems", J. Chem. Phys. 87, 731, doi:10.1063/1.453570. [ABS]
  - ZGB on hexagonal, 1D and strip lattices: for 1D and width-2 strips no reactive steady state exists; for widths ≥ 3 the reactive window grows with width.
  - Also presents new results for the monomer–monomer C+D reaction.
  - Relevance: Low to Medium. [Source](https://doi.org/10.1063/1.453570)
- **Evans & Miesch 1991**, "Characterizing kinetics near a first-order catalytic-poisoning transition", Phys. Rev. Lett. 66, 833, doi:10.1103/PhysRevLett.66.833. [ABS]
  - Epidemic analysis of an empty patch in the CO-poisoned phase. Exponents "vary strongly with reaction rate".
  - "The increase in epidemic survival probability with initial patch size is smooth … rather than indicating a critical size for survival."
  - Evidence: Moderate. Relevance: High for signatures 4 and 5. [Source](https://doi.org/10.1103/PhysRevLett.66.833)
- **Brosilow & Ziff 1992**, "Effects of A desorption on the first-order transition in the A-B₂ reaction model", Phys. Rev. A 46, 4534, doi:10.1103/PhysRevA.46.4534. [ABS]
  - Above a desorption probability Pc ≈ 0.039 there is no first-order transition.
  - The coexistence interface is non-fractal below Pc and has fractal dimension D = 1.46(5) at Pc.
  - Relevance: Medium (absorbing-breaking noise ends the first-order line). [Source](https://doi.org/10.1103/PhysRevA.46.4534)
- **Ziff & Brosilow 1992**, "Investigation of the first-order phase transition in the A-B₂ reaction model using a constant-coverage kinetic ensemble", Phys. Rev. A 46, 4630, doi:10.1103/PhysRevA.46.4630. [ABS]
  - The constant-coverage ensemble "allows the metastability loop and spinodal point … to be carefully examined"; y₂ = 0.52560(1).
  - Evidence: Strong. Relevance: High (method). [Source](https://doi.org/10.1103/PhysRevA.46.4630)
- **Tomé & Dickman 1993**, "Ziff-Gulari-Barshad model with CO desorption: An Ising-like nonequilibrium critical point", Phys. Rev. E 47, 948, doi:10.1103/PhysRevE.47.948. [ABS]. With CO desorption k the first-order line "terminates at a critical point k_c", in the Ising class. Relevance: Medium. [Source](https://doi.org/10.1103/PhysRevE.47.948)
- **Evans & Ray 1994**, "Interface propagation and nucleation phenomena for discontinuous poisoning transitions in surface-reaction models", Phys. Rev. E 50, 4302, doi:10.1103/PhysRevE.50.4302. [ABS]
  - Reactive/poisoned interface: "Kardar-Parisi-Zhang equation, but with the nonlinearity and kinetic surface tension vanishing on approaching the transition".
  - Nuclei: "divergence and 'sharpening' of the critical size … The deviation from mean-field divergence approaching the transition is related to the vanishing of the kinetic surface tension. The sharpening is related to the reduced influence of fluctuations."
  - Evidence: Strong. Relevance: High (closest methodological precedent for D&N's fronts and droplets). [Source](https://doi.org/10.1103/PhysRevE.50.4302)
- **Machado, Buendía & Rikvold 2005**, "Decay of metastable phases in a model for the catalytic oxidation of CO", Phys. Rev. E 71, 031603, doi:10.1103/PhysRevE.71.031603, arXiv:cond-mat/0411335. [ABS+PT]
  - Finite-size scaling of fluctuations and the fourth-order cumulant shows a first-order transition below a critical CO desorption rate.
  - Lifetimes follow KJMA nucleation and growth; there are "two distinct regimes … multidroplet or single-droplet decay", and "the crossover … approaches the coexistence curve logarithmically with system size".
  - Evidence: Strong. Relevance: High for signature 3. [Source](https://doi.org/10.1103/PhysRevE.71.031603)
  - Note: the four-author **Machado, Buendía, Rikvold & Ziff 2005**, Phys. Rev. E 71, 016120, doi:10.1103/PhysRevE.71.016120 [ABS], is a *different* paper. It covers periodic CO-pressure forcing and a dynamic phase transition in the Ising class. Hysteresis and constant coverage belong to Ziff & Brosilow 1992. [Source](https://doi.org/10.1103/PhysRevE.71.016120)
- **Loscar & Albano 2003**, "Critical behaviour of irreversible reaction systems", Rep. Prog. Phys. 66, 1343, doi:10.1088/0034-4885/66/8/203 [ABS]. Review with "special attention … to … first-order IPTs". Relevance: Medium. [Source](https://doi.org/10.1088/0034-4885/66/8/203)

#### Monomer–monomer and voter-type catalytic models (two absorbing states, no separate active phase)
- **Frachebourg & Krapivsky 1996**, "Exact results for kinetics of catalytic reactions", Phys. Rev. E 53, R3009, doi:10.1103/PhysRevE.53.R3009, arXiv:cond-mat/9508123. [ABS]. In the reaction-controlled limit the dimer–dimer model *is* the voter model: "power law decay for D<2 and a slow logarithmic decay in two dimensions", with the relevance to monomer–monomer discussed. Relevance: Medium. [Source](https://doi.org/10.1103/PhysRevE.53.R3009)

#### Quadratic contact process, Schlögl's second model, sexual reproduction, threshold contact processes (single absorbing state)
- **Liu, Guo & Evans 2007**, "Quadratic contact process: Phase separation with interface-orientation-dependent equistability", Phys. Rev. Lett. 98, 050601, doi:10.1103/PhysRevLett.98.050601. [ABS; numbers from Guo's 2009 thesis, PT]
  - "discontinuous phase transition … but equistability between these states depends on the orientation of the separating interface … phase coexistence over a finite region of their two-dimensional parameter space."
  - Evidence: Strong. Relevance: High. [Source](https://doi.org/10.1103/PhysRevLett.98.050601)
- **Guo, Liu & Evans 2007**, "Generic two-phase coexistence, relaxation kinetics, and interface propagation in the quadratic contact process: Simulation studies", Phys. Rev. E 75, 061129, doi:10.1103/PhysRevE.75.061129. [ABS]
  - Covers rapid poisoning above an effective spinodal, "nucleation-mediated poisoning in the metastable region, the dynamics of poisoned droplets", and orientation-dependent planar interfaces.
  - Companion: **Guo, Evans & Liu 2008**, "… Analytic studies", Physica A 387, 177, doi:10.1016/j.physa.2007.09.002 [title and metadata only].
  - [Source](https://doi.org/10.1103/PhysRevE.75.061129); [Source](https://doi.org/10.1016/j.physa.2007.09.002)
- **Guo, Liu & Evans 2009**, "Schloegl's second model for autocatalysis with particle diffusion: Lattice-gas realization exhibiting generic two-phase coexistence", J. Chem. Phys. 130, 074106, doi:10.1063/1.3074308 [ABS]. Also X. Guo, PhD thesis, Iowa State 2009, OSTI 1342563 [PT].
  - With hopping h > 0: "p(f)<p(e). This strict inequality contrasts a postulate of Durrett".
  - Thesis numbers at h = 0: p_eq(S=1) = p_e = 0.09443(3); p_eq(S→0 or ∞) = 0.0869(4); effective spinodal ≈ 0.0995(5).
  - The orientation dependence "quickly weaken[s] with increasing h".
  - Evidence: Strong. Relevance: High. [Source](https://doi.org/10.1063/1.3074308); [thesis](https://www.osti.gov/biblio/1342563)
- **Wang, Liu & Evans 2012**, "Schloegl's second model for autocatalysis on hypercubic lattices: Dimension dependence of generic two-phase coexistence", Phys. Rev. E 85, 041109, doi:10.1103/PhysRevE.85.041109. [ABS]
  - Coexistence for p_f(d) < p < p_e(d) "spanned by the orientation-dependent stationary points for planar interfaces", with p_e(f) ≈ 0.2113765 + c_e(f)/d and Δc ≈ 0.014.
  - "The dRDEs display artificial propagation failure, absent in the stochastic model due to fluctuations."
  - Relevance: High. [Source](https://doi.org/10.1103/PhysRevE.85.041109)
- **Wang, Liu & Evans 2015**, "Discontinuous non-equilibrium phase transition in a threshold Schloegl model for autocatalysis: Generic two-phase coexistence and metastability", J. Chem. Phys. 142, 164105, doi:10.1063/1.4918908 (OSTI 1227298). [PT, Secs. 1 and 3]. N = 2 threshold contact process on the square lattice. Details under Q3 and Q5. Evidence: Strong. Relevance: High. [Source](https://www.osti.gov/biblio/1227298)
- **da Silva & de Oliveira 2011**, "Critical discontinuous phase transition in the threshold contact process", J. Phys. A 44, 135002, doi:10.1088/1751-8113/44/13/135002. [ABS]
  - Same N = 2 model, run with borders "permanently occupied". Claims a discontinuous transition "at a definite value … in accordance with the Gibbs phase rule" that "exhibits critical behavior".
  - Disputed by Wang–Liu–Evans 2015 on both counts (orientation-dependent coexistence; non-zero line tension). Evidence: Disputed. Relevance: High as a methodological cautionary tale. [Source](https://doi.org/10.1088/1751-8113/44/13/135002)
- **Fiore 2014**, "Minimal mechanism leading to discontinuous phase transitions for short-range systems with absorbing states", Phys. Rev. E 89, 022104, doi:10.1103/PhysRevE.89.022104, arXiv:1312.2408. [ABS+PT]
  - Contact-process variants on the square lattice (z = 4) "in which at least two adjacent particles … are required to create" offspring.
  - Evidence used: bimodal order-parameter distributions and "finite scale analysis".
  - "Unlike the original CP, the transition becomes … discontinuous for dimensions larger than 1."
  - Evidence: Moderate. Relevance: Medium. [Source](https://doi.org/10.1103/PhysRevE.89.022104)
- **Pianegonda & Fiore**, two papers. Relevance: Low to Medium.
  - J. Stat. Mech. (2015) P08018, doi:10.1088/1742-5468/2015/08/P08018 [ABS]: in two lattice models diffusion "does not change the order of the transition, regardless its strength", in partial contrast with Langevin results. [Source](https://doi.org/10.1088/1742-5468/2015/08/P08018)
  - Physica A 451, 349 (2016), doi:10.1016/j.physa.2016.01.047, arXiv:1603.08835 [ABS]: "limited competitive rates suppress the phase coexistence". [Source](https://doi.org/10.1016/j.physa.2016.01.047)
- **Windus & Jensen 2007**, "Phase transitions in a lattice population model", J. Phys. A 40, 2287, doi:10.1088/1751-8113/40/10/005, arXiv:cond-mat/0612447 [ABS]. First order in 2+1 D, located through "the presence of phase coexistence" in density histograms. This model is the pure first-order reference used by Villa Martín et al. 2014. Related: **Windus & Jensen 2007**, "Allee effects and extinction in a lattice model", Theor. Popul. Biol. 72, 459, doi:10.1016/j.tpb.2007.07.006 [CIT]. [Source](https://doi.org/10.1088/1751-8113/40/10/005)
- **Varghese & Durrett 2013**, "Phase transitions in the quadratic contact process on complex networks", Phys. Rev. E 87, 062819, doi:10.1103/PhysRevE.87.062819 [ABS]. Discontinuous with bistability on random-regular and Erdős–Rényi graphs; continuous on power-law graphs. Not 2D lattices. Relevance: Low. [Source](https://doi.org/10.1103/PhysRevE.87.062819)
- **Rigorous sexual-reproduction process.** Relevance: Low to Medium.
  - **H.-N. Chen 1992**, "On the stability of a population growth model with sexual reproduction on Z²", Ann. Probab. 20, 232, doi:10.1214/aop/1176989926 [ABS]. In the symmetric Z² model the vacant state is stable against small initial densities but "unstable under the second kind of perturbation" (spontaneous birth β > 0). [Source](https://doi.org/10.1214/aop/1176989926)
  - **Chen 1994**, Ann. Probab. 22, 1195, doi:10.1214/aop/1176988600 [ABS]: results depend on the rule geometry (type III vs IV) and on dimension. [Source](https://doi.org/10.1214/aop/1176988600)
  - **Noble 1992**, Ann. Probab. 20, 724, doi:10.1214/aop/1176989802 [ABS]: rapid diffusion. [Source](https://doi.org/10.1214/aop/1176989802)
  - **Durrett & Neuhauser 1994**, Ann. Probab. 22, 289, doi:10.1214/aop/1176988861 [ABS]: fast-stirring limit leading to reaction–diffusion. [Source](https://doi.org/10.1214/aop/1176988861)
  - **Durrett 1999**, "Stochastic spatial models", SIAM Rev. 41, 677, doi:10.1137/S0036144599354707 [ABS]: the Durrett–Levin idea that behaviour is predicted from the mean-field ODE. Per Guo et al. 2009 this is the source of the "Durrett postulate" p_f = p_e. [Source](https://doi.org/10.1137/S0036144599354707)
- **Schlögl 1972**, "Chemical reaction models for non-equilibrium phase transitions", Z. Phys. 253, 147, doi:10.1007/BF01379769 [CIT]: the origin of the models. [Source](https://doi.org/10.1007/BF01379769)

#### Synchronous PCAs (closest to D&N's synchronous outer-totalistic rule)
- **Bidaux, Boccara & Chaté 1989**, "Order of the transition versus space dimension in a family of cellular automata", Phys. Rev. A 39, 3094, doi:10.1103/PhysRevA.39.3094. [ABS]
  - "At odds with usual beliefs on two-state automata with one absorbing phase, first-order transitions are indeed exhibited as soon as d>1 … For d=1, the transition is continuous, but with critical exponents different from those of directed percolation."
  - Evidence: Moderate. Relevance: High (2D synchronous PCA, first order).
  - Correction: **Jensen 1991**, "Universality class of a one-dimensional cellular automaton", Phys. Rev. A 43, 3187, doi:10.1103/PhysRevA.43.3187 [ABS]. Time-dependent simulations give 1D exponents "in excellent agreement with directed percolation", showing that steady-state fits had produced spurious non-DP exponents.
  - [Source](https://doi.org/10.1103/PhysRevA.39.3094); [Source](https://doi.org/10.1103/PhysRevA.43.3187)
- **Ódor & Szolnoki 1996**, "Directed-percolation conjecture for cellular automata", Phys. Rev. E 53, 2231, doi:10.1103/PhysRevE.53.2231, arXiv:cond-mat/9506139. [ABS+PT]
  - Rule: s(t+1) = X (Prob(X=1) = p) if y ≤ Σ ≤ 6, else 0 (the BBC rule family).
  - 2D: y = 1 is DP; "in the cases of y>1 the transitions become first order".
  - Method: "Since it is not easy to prove first order transition by simulation, we applied the generalized mean-field method" (cluster-size convergence of the jump).
  - 1D: y = 4, 5 are also reported discontinuous (range-4 rule).
  - Evidence: Moderate. Relevance: High. [Source](https://doi.org/10.1103/PhysRevE.53.2231)
- **Ódor, Boccara & Szabó 1993**, Phys. Rev. E 48, 3168, doi:10.1103/PhysRevE.48.3168 [ABS]. In a 1D PCA, mixing changes the order "from second order to first order", with a tricritical point. Relevance: Low. [Source](https://doi.org/10.1103/PhysRevE.48.3168)

#### Tricritical DP (first-order lines in 2D)
- **Lübeck 2006**, "Tricritical directed percolation", J. Stat. Phys. 123, 193, doi:10.1007/s10955-005-9000-z, arXiv:cond-mat/0608339. [PT, Sec. 6]. 2D tricritical contact process (pair-reaction terms added).
  - Superheating limit: "in contrast to a continuous phase transition the transition point p_o exhibits no systematic system size dependence".
  - Supercooling limit p_u obtained with a small external field h → 0, giving "a small but finite hysteresis, i.e., the two phases coexist between p_o < p < p_u".
  - Stripe test (one orientation): at q = 0.97 the outcome is ~50/50 at p = 0.2217, and "this interval decreases with increasing system size L"; Δp/p_c ≈ 0.009.
  - "the first-order line ends for q → 1 at a finite p value in contrast to the mean field phase diagram."
  - Warning: "it is notoriously difficult to distinguish a continuous phase transition with a small but finite value of the exponent β_t from a discontinuous jump."
  - Evidence: Strong. Relevance: High. [Source](https://doi.org/10.1007/s10955-005-9000-z)
- **Grassberger 2006**, "Tricritical directed percolation in 2+1 dimensions", J. Stat. Mech. P01004, doi:10.1088/1742-5468/2006/01/P01004, arXiv:cond-mat/0510428. [ABS+PT]
  - Generalised Domany–Kinzel in 2+1 D: "If … P_k increases quickly with k so that the formation of inactive holes surrounded by active sites is suppressed, the transition is first order."
  - p_c(q) is fixed by "the condition that neither phase will win" for a straight interface.
  - Cluster growth is "similar to the growth of a droplet in a usual first order phase transition, except that now one of the two phases is not fluctuating."
  - Tricritical exponents agree "rather poor[ly]" with the ε = 3 − d expansion.
  - Evidence: Strong. Relevance: High. [Source](https://doi.org/10.1088/1742-5468/2006/01/P01004)
- **Jo & Kahng 2020**, "Tricritical directed percolation with long-range interaction in one and two dimensions", Phys. Rev. E 101, 022121, doi:10.1103/PhysRevE.101.022121 [ABS]. Relevance: Low. [Source](https://doi.org/10.1103/PhysRevE.101.022121)

#### Experiments and other 2D settings
- **Néel, Rondini, Turzillo, Mujica & Soto 2014**, "Dynamics of a first-order transition to an absorbing state", Phys. Rev. E 89, 042206, doi:10.1103/PhysRevE.89.042206, arXiv:1401.1839. [ABS]
  - Quasi-2D vibrated granular layer: first order "with long metastable periods, followed by rapid transitions driven by homogeneous nucleation".
  - A Ginzburg–Landau model with multiplicative noise "explains the macroscopic size of the critical nuclei".
  - Relevance: Medium. [Source](https://doi.org/10.1103/PhysRevE.89.042206)
- **Maire, Plati, Smallenburg & Foffi 2025**, J. Stat. Mech. 123206, doi:10.1088/1742-5468/ae23bc [ABS]: the granular absorbing transition can be "both continuous and discontinuous" depending on parameters. Relevance: Low. [Source](https://doi.org/10.1088/1742-5468/ae23bc)
- **Bertin, Andrix & Le Godais 2024**, J. Stat. Phys. 191, 113, doi:10.1007/s10955-024-03327-3 [ABS]: the laminar–turbulent transition treated as a discontinuous absorbing transition, with a "faster-than-exponential increase of the turbulence lifetime". Relevance: Low. [Source](https://doi.org/10.1007/s10955-024-03327-3)

#### Symmetric multi-absorbing-state models in 2D (most relevant prior art beyond LITERATURE.md)
- **Lipowski & Droz 2002**, "Phase transitions in nonequilibrium d-dimensional models with q absorbing states", Phys. Rev. E 65, 056114, doi:10.1103/PhysRevE.65.056114, arXiv:cond-mat/0112130. [PT]. Nonequilibrium Potts model (Metropolis rule with sites frozen when all neighbours agree); q symmetric absorbing states; square lattice; L = 500.
  - d = 2, q = 3: "undergoes a discontinuous phase transition". Spreading P(t) and N(t) bend away from power laws.
  - d = 2, q = 2: an apparent jump in ρ "almost twice smaller" than for q = 3, but "an observed jump might be a finite-size effect".
    - ρ(t) shows no clear plateau, and δ′ = 0.07(2) "might be equal to zero".
    - Spreading gives P(t) ~ t^−0.90(2) over three decades with N(t) ≈ const.
    - Verdict: "Simulations are inconclusive … in the two-dimensional case for q = 2 … close to or at the crossing point of lines separating three different types of phase transitions".
    - Below Tc, ρ(t) ~ t^−1/2 because "active sites are located mainly at the domain walls".
  - Ódor's review summarises this as suggesting first order for q > 2 in d > 1.
  - Evidence: Moderate (q = 3), Inconclusive (q = 2). Relevance: High, the direct precedent for the D&N claim.
  - [Source](https://doi.org/10.1103/PhysRevE.65.056114). 1D companion: Lipowski & Droz, Phys. Rev. E 66, 016106 (2002), doi:10.1103/PhysRevE.66.016106 [ABS]. [Source](https://doi.org/10.1103/PhysRevE.66.016106)
- **Lee & Vojta 2011**, "Generalized contact process with two symmetric absorbing states in two dimensions", Phys. Rev. E 83, 011114, doi:10.1103/PhysRevE.83.011114, arXiv:1010.3298. [ABS+PT]. The generic transition is GV ("symmetry-breaking and absorbing transitions coincide"); "the transition at zero domain-boundary activation rate is not critical". Domain walls between the two inactive states "can hop, branch and annihilate". Relevance: High. [Source](https://doi.org/10.1103/PhysRevE.83.011114)
- **2D interacting monomer–dimer with two symmetric absorbing states.** Relevance: Medium (shows how fragile cumulant-based locating is in Z2-absorbing models).
  - **Nam, Park, Kim & Lee 2011**, J. Stat. Mech. L06001, doi:10.1088/1742-5468/2011/06/L06001 [ABS]: "two nearby continuous transitions", a non-Ising order–disorder transition, and a "critical" absorbing phase. [Source](https://doi.org/10.1088/1742-5468/2011/06/L06001)
  - Rebutted by **Park 2012**, Phys. Rev. E 85, 041140, doi:10.1103/PhysRevE.85.041140, and **Park 2015**, J. Stat. Mech. P10009, doi:10.1088/1742-5468/2015/10/P10009 [ABS]: Ising class, and "the Binder cumulant is not an efficient quantity to locate the order-disorder transition point of the model with two symmetric absorbing states". [Source](https://doi.org/10.1088/1742-5468/2015/10/P10009)
- **Dornic, Chaté, Chave & Hinrichsen 2001**, "Critical coarsening without surface tension: The universality class of the voter model", Phys. Rev. Lett. 87, 045701, doi:10.1103/PhysRevLett.87.045701 [ABS]: "voter-like growth is generically observed at order-disorder nonequilibrium transitions solely driven by interfacial noise between dynamically symmetric absorbing states." Relevance: High (the null hypothesis that D&N must beat). [Source](https://doi.org/10.1103/PhysRevLett.87.045701)
- **Hinrichsen 1997**, "Stochastic lattice models with several absorbing states", Phys. Rev. E 55, 219, doi:10.1103/PhysRevE.55.219 [ABS]: 1D, n = 2 gives the BAW-even (parity-conserving) class without explicit parity conservation. Relevance: Low. [Source](https://doi.org/10.1103/PhysRevE.55.219)
- **1D PCAs with two absorbing states.**
  - **Bagnoli, Boccara & Palmerini 1997**, arXiv:cond-mat/9705171 (publication not verified) [ABS]: "a first-order phase transition between the fully occupied state and the empty state", with DP lines elsewhere. [Source](https://arxiv.org/abs/cond-mat/9705171)
  - **Bagnoli, Franci & Rechtman 2005**, Phys. Rev. E 71, 046108, doi:10.1103/PhysRevE.71.046108 [ABS]: an extended-range 1D version. [Source](https://doi.org/10.1103/PhysRevE.71.046108)
  - Relevance: Low to Medium (these are transitions between the two absorbing states, not active vs absorbing).

### Inferences
- The accepted 2D first-order cases all have **a stationary planar interface with a definite velocity of either sign, plus metastability**. The 1D false positive (TCM) lacked a stable interface. On this axis D&N resembles the accepted cases, not the false positives.
- The closest prior art for "two symmetric absorbing states, 2D" is Lipowski–Droz 2002 (q = 2: inconclusive, with voter-like spreading δ ≈ 0.9, η ≈ 0). A defensible novelty statement is therefore: *the first 2D model with two Z2-symmetric absorbing states shown to have a fully first-order active/absorbing transition, using interface, nucleation and wetting diagnostics rather than a finite-size jump.* It settles a case the (q, d) phase-diagram literature left open. This also needs the Rodrigues 2015 distinction already in LITERATURE.md.
- D&N belongs to the BBC and Ódor–Szolnoki class of synchronous outer-totalistic PCAs with upper thresholds, where 2D first order was found early. Jensen's 1991 correction of BBC's 1D exponents is a warning about steady-state fits.

### Gaps
- Full texts not read: Liu–Guo–Evans PRL 2007 (the numbers come from Guo's thesis), BBC 1989, Evans & Ray 1994, Evans & Miesch 1991, Meakin & Scalapino 1987, Varghese & Durrett 2013, Chen 1992 and 1994, Durrett 1999.
- The caller's "da Silva & Fiore" and "Fiore, Landi & de Oliveira" were not found as such. What exists: da Silva & de Oliveira 2011; Fiore & Landi 2014; Fiore, Harunari, Fernández Noa & Landi, Phys. Rev. E 104, 064123 (2021), "Current fluctuations in nonequilibrium discontinuous phase transitions", doi:10.1103/PhysRevE.104.064123 (title only) ([Source](https://doi.org/10.1103/PhysRevE.104.064123)).
- Ódor's 2003 triplet and quadruplet first-order results were seen only as cited in his review.

## Q3. Generic two-phase coexistence (GTPC) and orientation-dependent equistability: what follows if q* depends on front orientation

### Takeaway
In nonequilibrium systems the coexistence point is defined by a stationary planar interface. That point can depend on interface orientation; this is generic for nearest-neighbour square-lattice threshold and quadratic contact processes. When it does, the literature (Toom; Bennett–Grinstein; Muñoz et al.; Liu, Guo, Wang & Evans) says the single transition point becomes a **finite coexistence region [min_S q*(S), max_S q*(S)]**. Inside it, sufficiently large droplets of *either* phase become polygons bounded by the "losing" orientations and shrink at constant speed, so both phases are stable. The threshold for survival of finite seeds (p_f) then differs from the stability limit of the bulk phase (p_e), contradicting Durrett's postulate p_f = p_e.

How large the window is depends strongly on the rule and neighbourhood:
- ~6–11% of p_eq for nearest-neighbour rules;
- for the **Moore neighbourhood** in the threshold-2 Schlögl model, KMC finds **no detectable orientation dependence** (p_eq = 0.14681(1) "for all S"), even though the pair approximation predicts ~6%.
- The window shrinks with diffusion, with absorbing-breaking noise and with neighbourhood size, and vanishes as N → ∞, where a Gibbs-like unique p_eq re-emerges.

I found no GTPC study for a model with Z2-symmetric absorbing states. All Z2-symmetric GTPC examples (Toom NEC, the totally asymmetric kinetic Ising model, the sweep rule) have no absorbing states.

### Cited Findings
- **Toom 1980** — A. L. Toom, "Stable and attractive trajectories in multicomponent systems", in *Multicomponent Random Systems*, eds. R. L. Dobrushin & Ya. G. Sinai, Adv. Probab. Relat. Topics 6 (Dekker, 1980), pp. 549–575. [CIT; citation taken from Hinrichsen 2000 ref. 2 and Muñoz et al. 2005 ref. 1]. Hinrichsen describes it: "the interfaces in Toom's model propagate in a preferred direction … this propagation velocity depends on the orientation of the domain walls … minority islands quickly assume a triangular shape and begin to shrink with constant velocity … they coexist in a whole region of the parameter space." [Source](https://arxiv.org/abs/cond-mat/0006212)
- **Bennett & Grinstein 1985**, "Role of irreversibility in stabilizing complex and nonergodic behavior in locally interacting discrete systems", Phys. Rev. Lett. 55, 657, doi:10.1103/PhysRevLett.55.657. [ABS]. "Irreversibility stabilizes certain locally interacting discrete systems against the nucleation and growth of a most-stable phase … over a set of positive measure in the parameter space". Relevance: High (conceptual). [Source](https://doi.org/10.1103/PhysRevLett.55.657)
- **Grinstein, Jayaprakash & He 1985**, "Statistical mechanics of probabilistic cellular automata", Phys. Rev. Lett. 55, 2527, doi:10.1103/PhysRevLett.55.2527 [ABS]. Continuous transitions in up–down-symmetric PCA are expected in the kinetic Ising class. Relevance: Low to Medium. [Source](https://doi.org/10.1103/PhysRevLett.55.2527)
- **Muñoz, de los Santos & Telo da Gama 2005**, "Generic two-phase coexistence in nonequilibrium systems", Eur. Phys. J. B 43, 73, doi:10.1140/epjb/e2005-00029-3, arXiv:cond-mat/0411578. [ABS+PT]. Interfacial models (a KPZ interface bound to a wall; a depinning model with stress overshoots).
  - Two ingredients: "dynamical asymmetry" and "a robust mechanism for the elimination of islands". In Toom's model, islands "shrink … at a velocity roughly independent of their radius … compared with the Ising model … where a droplet of radius r shrinks with a velocity proportional to 1/r".
  - "at least one of the boundaries [of the coexistence region] has to be discontinuous".
  - Relevance: High (conceptual). [Source](https://doi.org/10.1140/epjb/e2005-00029-3)
- **de los Santos, Telo da Gama & Muñoz 2003**, "Nonequilibrium wetting transitions with short range forces", Phys. Rev. E 67, 021607, doi:10.1103/PhysRevE.67.021607 [ABS]: "the finite coexistence region survives in the thermodynamic limit". [Source](https://doi.org/10.1103/PhysRevE.67.021607)
- **Wang, Liu & Evans 2015** (J. Chem. Phys. 142, 164105) [PT] is the most explicit recipe.
  - Definition: "it is convenient to define an equistability point, p = p_eq, corresponding to a stationary planar interface."
  - Measured: "p_e = max_S p_eq(S) = p_eq(S=1) = 0.35178(2) is significantly above p_eq(S=2) = 0.3497(1) and p_f = min_S p_eq(S) = p_eq(S=0 or ∞) = 0.3161(1)."
  - Droplet shapes: "for p₋ < p < p₊, an embedded vacuum droplet will evolve to a diamond shape and then shrink, and an embedded populated droplet will evolve to a square shape and then shrink."
  - Shrink law: "droplets of either the populated or vacuum state embedded in the other state shrink at constant rate".
  - Pitfall: "(This feature is obscured starting with a random configuration [23].)"
  - Method: constant-concentration (Ziff–Brosilow) runs started from strips; slopes S other than 0, 1, ∞ use rectangular periodic cells of aspect ratio S.
  - Practical rule: "one can simply stipulate that embedded droplets must shrink at constant rate in order for a state to be assigned as stable."
  - Evidence: Strong. Relevance: **High**. [Source](https://www.osti.gov/biblio/1227298)
- **Guo et al. (QCP, Durrett-type rates)**, Guo thesis 2009 [PT].
  - GTPC for p_eq(S=∞) ≈ 0.0869 < p < p_eq(S=1) ≈ 0.0944.
  - "we show that p_f = p_eq(S=0) and that p_f < p_e. This strict inequality contrasts a postulate of Durrett".
  - The orientation dependence "quickly weaken[s] with increasing h" (hop rate).
  - [Source](https://www.osti.gov/biblio/1342563); [Source](https://doi.org/10.1063/1.3074308)
- **Liu, Wang & Evans 2018**, "Discontinuous phase transitions in nonlocal Schloegl models for autocatalysis: Loss and reemergence of a nonequilibrium Gibbs phase rule", Phys. Rev. Lett. 121, 120603, doi:10.1103/PhysRevLett.121.120603 (OSTI 1477196). [FT]. Constant-concentration strip runs on a 1024² lattice; p_eq(S) per neighbourhood:

  | Neighbourhood | p_eq range over S | Relative spread |
  |---|---|---|
  | NN, N = 4 | 0.080480(8) at S = 1/6 to 0.084941(40) at S = 1 | ~5% |
  | 4NN (r = √5 a), N = 8 | 0.151121(2) at S = 0 to 0.151037(2) at S = 1 | ~0.06% |
  | NN + 3NN (r = a, 2a), N = 8 | 0.147948(13) to 0.148134(10) | ~0.13% |
  | **Moore (NN + 2NN, r = a, √2 a), N = 8** | **0.146802(2) to 0.146808(1), no trend** | **< 1e-4** |
  | N = 12 (r = a, 2a, 3a) | 0.168721(12) to 0.169000(11) | ~0.17% |

  - On the Moore case: "For N = 8 with NN & 2NN sites (r = a, √2a) and for N>8, it is difficult to quantify any weak S-dependence of p_eq(S), although plausibly this persists."
  - "In the simulations … with random initial conditions, p_eq reflects an average over S of p_eq(S)."
  - Large-N limit: "while interface profiles shape depends on orientation as N→∞, we find a unique equistability point … We thus establish a Gibbs phase rule for these non-equilibrium models as N→∞."
  - I confirmed the Moore-column identification visually from the PDF table.
  - Evidence: Strong. Relevance: **High**, since D&N also uses the Moore neighbourhood. [Source](https://www.osti.gov/biblio/1477196); [Source](https://doi.org/10.1103/PhysRevLett.121.120603)
- **Liu, Shen & Evans 2025**, "Trends in spatial correlations, two-phase coexistence, and criticality in a class of type-2 Schloegl models for autocatalysis", Phys. Rev. E 111, 064134, doi:10.1103/y54y-5dvx (OSTI 2571583). [PT incl. Table I]
  - KMC p_eq for N = 8, d = {1, √2} (Moore): "0.14681(1) for all S". The pair approximation gives 0.122 (S = 1) to 0.130 (S = ∞).
  - N = 4: 0.0798 (S = ∞) to 0.08494 (S = 1).
  - "the width, Δp, of the region of generic 2PC decreases quickly with increasing ε > 0 [spontaneous creation] … vanish as a critical point ε_c ≈ 0.005" for N = 4.
  - ε_c ≈ 0.0043 for the Moore neighbourhood (0.0041 in the Fig. 5 caption), against the mean-field value 1/27, attributed to "the compact short-range choice".
  - Fluctuation FSS crossings sit close to the 2D-Ising value for N = 4. For the Moore neighbourhood they cross lower (more mean-field-like), and the authors only "anticipate" the Ising class.
  - Evidence: Strong. Relevance: **High**. [Source](https://doi.org/10.1103/y54y-5dvx)
- **Shen, Liu & Evans 2023**, Phys. Rev. E 107, 034104, doi:10.1103/PhysRevE.107.034104 [ABS]. Pair-approximation LDEs reproduce the orientation dependence within 15%. Kinks on near-vertical interfaces are stationary at p_eq(min). [Source](https://doi.org/10.1103/PhysRevE.107.034104)
- **Wang, Liu & Evans 2020**, "Extended families of critical and stationary droplets for nonequilibrium phase transitions in spatially discrete bistable systems", Phys. Rev. E 101, 022803, doi:10.1103/PhysRevE.101.022803 (OSTI 1602867). [ABS+intro]. In discrete deterministic LDEs, propagation failure plus orientation-dependent propagation give "entire families of stationary droplets", and "there can exist a regime of generic two-phase coexistence where arbitrarily large droplets of either phase always shrink". Relevance: High. [Source](https://doi.org/10.1103/PhysRevE.101.022803)
- **GTPC in Z2-symmetric systems, all without absorbing states.** Relevance: Medium.
  - **Godrèche & Luck 2017**, "Generic phase coexistence in the totally asymmetric kinetic Ising model", J. Stat. Mech. 073208, doi:10.1088/1742-5468/aa79b0 [ABS]: "The gist of their [Bennett–Grinstein] argument relies on the dynamics of interfaces and droplets". [Source](https://doi.org/10.1088/1742-5468/aa79b0)
  - **Ray, Laflamme & Kubica 2024**, Phys. Rev. E 109, 044141, doi:10.1103/PhysRevE.109.044141 [ABS]: the "sweep rule" maps "a region of two coexisting stable phases in the (T,h) plane". [Source](https://doi.org/10.1103/PhysRevE.109.044141)
- **Durrett & Levin 1994**, "The importance of being discrete (and spatial)", Theor. Popul. Biol. 46, 363, doi:10.1006/tpbi.1994.1032 [ABS]. Spatial stochastic models can disagree with mean field. Via Durrett 1999 and Weissmann–Shnerb this is the origin of the "the sign of the wave speed decides" (Maxwell-point) rule. Relevance: Medium. [Source](https://doi.org/10.1006/tpbi.1994.1032)

### Inferences
- **Direct answer to the key question.** If q* depends on orientation, the literature says the stationary first-order point becomes a **GTPC window**:
  - lower edge q_lo = min_S q*(S): below it absorbing droplets grow in every orientation, so the active phase is only metastable;
  - upper edge q_hi = max_S q*(S): above it active droplets grow in every orientation;
  - in between, both phases are stable against large droplets, which facet onto the losing orientations and shrink at constant speed (not ∝ 1/R).
- The "true transition point" is then not unique:
  - a fully active initial state stays active down to q_lo, up to nucleation from below (spinodal-limited);
  - a finite seed survives only above q_hi (the Evans-group p_f);
  - homogeneous random starts give an orientation-averaged pseudo-point.
- Practical consequence for D&N: measuring q*(S) for S = 0, 1/2, 1 and checking that droplets stay rounded (not polygonal) near q* is the single most important missing test. LWE18 and LSE25 suggest the Moore neighbourhood *may* be nearly isotropic, but that is for a monotone threshold rule. D&N's non-monotone rule (flip classes j ∈ {3,6,7,8}), its synchronous updating and its voter admixture are untested, so this cannot be assumed.
- Terminology: the [κ_front, κ_a] "coexistence window" in LITERATURE.md (front threshold vs the f′(0) low-density threshold) is a mean-field bistability or spinodal window. It is *not* GTPC, which requires both phases to be stable against arbitrarily large droplets. Reports should keep the two apart.
- If GTPC were found, the stronger result would be a 2D Z2 absorbing-state model with generic coexistence. I found no precedent, so it would be a natural "richer adjacent problem".

### Gaps
- No study found of orientation dependence for Z2-symmetric *absorbing* models, for synchronous PCAs with non-monotone rules, or for voter-perturbed rules.
- No published Wulff-type or kinetic-shape construction for stochastic absorbing-state droplets was found. The Evans group describes shapes qualitatively (diamond or square).
- The Liu–Guo–Evans PRL 2007 full text was not read.

## Q4. Quenched (and temporal) disorder

### Takeaway
By analogy with Imry–Ma, Aizenman–Wehr and Hui–Berker, Villa Martín, Bonachela & Muñoz (2014) argue and simulate that quenched disorder turns 2D first-order *absorbing* transitions continuous, giving disordered-DP criticality with a Griffiths phase. The stated reason is that the active phase fluctuates even though the absorbing phase does not. PNAS 2015 agrees for 2D (any heterogeneity smooths the transition) but not for 3D. Three caveats:
- Voronoi–Delaunay topological disorder does *not* destroy the ZGB first-order transition;
- temporal disorder does not forbid first order;
- random fields that bias between *two symmetric absorbing states* do not destroy the 1D GCP transition.

For a Z2 model, the relevant disorder (random-bond-like or random-field-like) and the expected outcome are open questions.

### Cited Findings
- **Aizenman & Wehr 1989**, "Rounding of first-order phase transitions in systems with quenched disorder", Phys. Rev. Lett. 62, 2503, doi:10.1103/PhysRevLett.62.2503. [ABS]. "in 2D quenched randomness results in the elimination of discontinuities in the density of the thermodynamic variable conjugate to the fluctuating parameter … 2D random-bond Potts models where a phase transition persists, but ceases to be first order." Evidence: Strong (rigorous, equilibrium). Relevance: Medium. [Source](https://doi.org/10.1103/PhysRevLett.62.2503)
- **Hui & Berker 1989**, "Random-field mechanism in random-bond multicritical systems", Phys. Rev. Lett. 62, 2507, doi:10.1103/PhysRevLett.62.2507 [ABS]: "under bond randomness, the phase transitions of q-state Potts models are second order for all q at dimensionality d≤2." [Source](https://doi.org/10.1103/PhysRevLett.62.2507)
- **Villa Martín, Bonachela & Muñoz 2014**, "Quenched disorder forbids discontinuous transitions in nonequilibrium low-dimensional systems", Phys. Rev. E 89, 012145, doi:10.1103/PhysRevE.89.012145, arXiv:1402.1406. [PT: intro, model, conclusions]. Pure 2D model: the Windus–Jensen pair-creation contact process with diffusion.
  - Pure-model evidence: an abrupt change in decay from full occupation at "0.070 < p_d_thr < 0.075"; bistability against initial density; mean survival time with a "multi-droplet" regime where it "ceases to depend on system size"; seeded spreading with N_s, R² ~ t² and P_s → const in the active phase.
  - Disordered model (site-dependent rates): a Griffiths phase and activated scaling at p_d ≈ 0.07652, compatible with DP plus quenched disorder.
  - Conjecture: "first-order phase transitions cannot appear in low-dimensional disordered systems with an absorbing state … the other phase is active and subject to fluctuation effects."
  - Exception: "the only possibility … would be … a discontinuous phase transition between two different fluctuation-less states".
  - Their text also says "systems exhibiting a first-order transition are bistable only at exactly the transition point". This is the conventional, Gibbs-like view, in tension with the Evans group's GTPC for closely related QCPs.
  - Evidence: Moderate (simulation plus heuristic). Relevance: High (robustness claims). [Source](https://doi.org/10.1103/PhysRevE.89.012145)
- **Caveat on VBM14's supporting example.** They cite Buendía & Rikvold 2012 as ZGB "in the presence of catalytic impurities (a fraction of inert sites)" whose first-order transition becomes continuous under quenched disorder. The cited paper's abstract concerns *gas-phase contaminants that adsorb*: "if the impurities can not desorb … the first order transition … disappear[s] … This is quite different from the behavior of a system in which the surface presents a fixed percentage of impurities."
  - G. M. Buendía & P. A. Rikvold, Phys. Rev. E 85, 031143 (2012), doi:10.1103/PhysRevE.85.031143, arXiv:1109.1315. [ABS]
  - So it is probably not a quenched-disorder test (judged from the abstract only). [Source](https://doi.org/10.1103/PhysRevE.85.031143)
- **de Oliveira, Alves & Ferreira 2016**, "Continuous and discontinuous absorbing-state phase transitions on Voronoi-Delaunay random lattices", Phys. Rev. E 93, 012110, doi:10.1103/PhysRevE.93.012110, arXiv:1509.05231. [ABS]. "the VD disorder does not change the nature of its [ZGB] discontinuous transition", consistent with Barghathi & Vojta 2014 on irrelevant topological disorder; "disorder in nonequilibrium APT may, under certain conditions, be irrelevant for the phase coexistence." Relevance: Medium. [Source](https://doi.org/10.1103/PhysRevE.93.012110)
  - Also: Barghathi & Vojta, "Phase transitions on random lattices: How random is topological disorder?", Phys. Rev. Lett. 113, 120602 (2014), doi:10.1103/PhysRevLett.113.120602 [title only]. [Source](https://doi.org/10.1103/PhysRevLett.113.120602)
- **Villa Martín, Bonachela, Levin & Muñoz 2015** (PNAS) [SUM]. In 2D "any amount of spatial heterogeneity induces a smooth transition", with a Griffiths phase. In 3D, disorder plus demographic noise "does not suffice". [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC4403222/)
- **Temporal disorder.** Relevance: Low to Medium.
  - **de Oliveira & Fiore 2016**, "Temporal disorder does not forbid discontinuous absorbing phase transitions in low-dimensional systems", Phys. Rev. E 94, 052138, doi:10.1103/PhysRevE.94.052138, arXiv:1603.08742 [ABS]. Three 2D models keep "bistability around the coexistence point and common finite-size scaling behavior with the inverse of the system volume", with no temporal Griffiths phase in d = 2. [Source](https://doi.org/10.1103/PhysRevE.94.052138)
  - **Encinas & Fiore 2021**, "Influence of distinct kinds of temporal disorder in discontinuous phase transitions", Phys. Rev. E 103, 032124, doi:10.1103/PhysRevE.103.032124 [title only]. [Source](https://doi.org/10.1103/PhysRevE.103.032124)
- **Barghathi & Vojta 2012**, "Random fields at a nonequilibrium phase transition", Phys. Rev. Lett. 109, 170603, doi:10.1103/PhysRevLett.109.170603, arXiv:1206.1878. [ABS+PT]
  - 1D generalized contact process with two symmetric absorbing states and random-field disorder (locally favouring one absorbing state): "the phase transition … persists", and "The ultraslow dynamics in the symmetry-broken phase is described by a Sinai walk of the domain walls between two different absorbing states".
  - "Although our explicit results are for one dimension, we expect our main conclusion to hold in higher dimensions, too."
  - Evidence: Moderate. Relevance: **High** for any Z2 disorder extension. [Source](https://doi.org/10.1103/PhysRevLett.109.170603)
- **Borile, Maritan & Muñoz 2013**, "The effect of quenched disorder in neutral theories", J. Stat. Mech. P04032, doi:10.1088/1742-5468/2013/04/P04032, arXiv:1304.0965 [ABS]. Voter and nonlinear-voter models (two symmetric absorbing states) with random-field-like disorder: coexistence "strongly favored by disorder". Relevance: Medium. [Source](https://doi.org/10.1088/1742-5468/2013/04/P04032)
- **Fytas, Theodorakis & Malakis 2017**, "Interfacial adsorption in two-dimensional pure and random-bond Potts models", Phys. Rev. E 95, 032126, doi:10.1103/PhysRevE.95.032126 [ABS]. This is the equilibrium analogue of disorder acting on the wetting signature: adsorption at "disordered induced" continuous transitions. Relevance: Medium. [Source](https://doi.org/10.1103/PhysRevE.95.032126)

### Inferences
- For the clean D&N lattice, disorder is not a flaw in the evidence. It does bound any claim of *physical robustness*: VBM14 predicts that site disorder in q (or α, p) would round the 2D transition into disordered-DP or GV-type criticality.
- In a Z2 model two disorder types act differently:
  - random-bond-type disorder (local preference active vs absorbing) is covered by VBM14's argument, because the active phase fluctuates;
  - random-field-type disorder (local preference 0 vs 1) acts on fluctuationless absorbing domains whose walls carry an active wetting slab. That combines the VBM14 and Barghathi–Vojta mechanisms, and I found no prior treatment.
- This is a well-posed "richer adjacent problem".

### Gaps
- No study found of quenched disorder at a first-order transition with *two symmetric* absorbing states in 2D.
- The VBM14 heuristic is not a theorem. VD-lattice results show exceptions.
- Barghathi & Vojta 2014 (PRL 113, 120602) was not read.

## Q5. Nucleation, spreading and wetting phenomenology

### Takeaway
Nucleation in absorbing-state first-order transitions is described with classical tools:
- planar-front velocity v ∝ Δp;
- critical droplet R_c ∝ σ_eff/Δp (Allen–Cahn: dR/dt = v − D/R);
- KJMA (Avrami) decay, C(t) ≈ C_m exp[−(t/τ)³] in 2D;
- an essential singularity, τ_nuc ≈ b δp^(−2/3) exp[c/(3δp)];
- single-droplet vs multidroplet regimes, with the crossover moving logarithmically with L.

Spreading from small seeds is a weak discriminator at first-order points. In the ZGB without diffusion, survival rises *smoothly* with patch size (no sharp critical size), and the critical size sharpens only when fluctuations are reduced. Wetting of the interface between two ordered phases by the disordered phase is textbook equilibrium physics at 2D q > 4 Potts transitions (adsorption W ∝ t^(−1/3)). I found **no nonequilibrium study of an active phase wetting the interface between two absorbing phases**, although activity sitting on domain walls between absorbing domains is documented in Z2 models.

### Cited Findings
- **Critical droplets and line tension.** Wang–Liu–Evans 2015 [PT].
  - "Rc ~ σ_eff/δp, where again σ_eff is an effective line tension". Critical vacuum droplets were measured with constant-concentration runs that stabilise them.
  - "C(t) ≈ C_m exp[-(t/τ_nuc)³], where τ_nuc ∝ (V_grow)^-2/3 (J_nuc)^-2/3" and "τ_nuc ≈ b (δp)^-2/3 exp[c_nuc/(3δp)]".
  - c_nuc ≈ 0.125 (threshold model) and ≈ 0.024 (Durrett model), hence "non-zero effective line tension".
  - "Note that droplet fate can be somewhat ambiguous given the fluctuating background in contrast to droplets embedded in an absorbing state".
  - Relevance: High. [Source](https://www.osti.gov/biblio/1227298)
- **Wang–Liu–Evans 2020** [ABS]. In the continuum "there is a unique critical droplet … with size diverging as p → p_eq". Discrete systems can instead have "entire families of stationary droplets". [Source](https://doi.org/10.1103/PhysRevE.101.022803)
- **Evans & Ray 1994** [ABS]: the critical size diverges and "sharpens", with deviations from mean-field divergence "related to the vanishing of the kinetic surface tension". [Source](https://doi.org/10.1103/PhysRevE.50.4302)
- **Evans & Miesch 1991** [ABS]: "The increase in epidemic survival probability with initial patch size is smooth … rather than indicating a critical size for survival." [Source](https://doi.org/10.1103/PhysRevLett.66.833)
- **Metastable decay.**
  - Machado–Buendía–Rikvold 2005 [ABS+PT]: KJMA; "multidroplet or single-droplet decay"; crossover ∝ log L. [Source](https://doi.org/10.1103/PhysRevE.71.031603)
  - **Rikvold, Tomita, Miyashita & Sides 1994**, "Metastable lifetimes in a kinetic Ising model: Dependence on field and system size", Phys. Rev. E 49, 5080, doi:10.1103/PhysRevE.49.5080 [ABS]: the lifetime is "an exponential in the inverse field … The field marking the crossover between the weak-field regime … single droplet, and the intermediate-field regime … finite density of droplets, vanishes logarithmically with system size." [Source](https://doi.org/10.1103/PhysRevE.49.5080)
  - **Meerson & Sasorov 2011**, "Extinction rates of established spatial populations", Phys. Rev. E 83, 011129, doi:10.1103/PhysRevE.83.011129 [ABS]: with a strong Allee effect, extinction "can be mapped into the overdamped limit of the theory of homogeneous nucleation due to Langer". [Source](https://doi.org/10.1103/PhysRevE.83.011129)
- **Spreading at first-order and tricritical points.**
  - Grassberger 2006 [PT]: "compact clusters with only few small holes and a slowly moving rough surface"; growth "similar to the growth of a droplet". [Source](https://doi.org/10.1088/1742-5468/2006/01/P01004)
  - **Windus & Jensen 2008**, "Cluster geometry and survival probability in systems driven by reaction-diffusion dynamics", New J. Phys. 10, 113023, doi:10.1088/1367-2630/10/11/113023 [ABS]: "a linear relationship for the probability of survival against initial cluster size at the tricritical point only". [Source](https://doi.org/10.1088/1367-2630/10/11/113023)
  - Hinrichsen 2000 [FT]: a CDP-like transient (δ = 1/2, η = 0) masquerading as first order in 1D. [Source](https://arxiv.org/abs/cond-mat/0006212)
  - Lipowski–Droz 2002 [PT]: 2D q = 2 gives P(t) ~ t^−0.90(2) and N ≈ const, i.e. voter-like. [Source](https://doi.org/10.1103/PhysRevE.65.056114)
  - VBM14 [PT]: in the active phase N_s, R² ∝ t² (ballistic) and P_s → const. [Source](https://doi.org/10.1103/PhysRevE.89.012145)
- **Seed survival vs bulk stability.** The Durrett postulate p_f = p_e fails when GTPC is present (Guo et al. 2009: p_f = p_eq(S=0) < p_e = p_eq(S=1)). [Source](https://doi.org/10.1063/1.3074308)
- **Wetting, equilibrium analogue.** Relevance: High (the template for signature 6).
  - **Selke & Pesch 1982**, Z. Phys. B 47, 335 [CIT: confirmed as a reference in arXiv:1610.08224; title and abstract not seen]. [Source](https://arxiv.org/abs/1610.08224)
  - **Selke & Huse 1983**, "Interfacial adsorption in planar Potts models", Z. Phys. B 50, 113, doi:10.1007/BF01304093 [metadata only]. [Source](https://doi.org/10.1007/BF01304093)
  - **Carlon, Iglói, Selke & Szalma 1999**, "Interfacial adsorption in two-dimensional Potts models", J. Stat. Phys. 96, 531, doi:10.1023/A:1004542105635, arXiv:cond-mat/9902182 [ABS+PT].
    - At the first-order transition (q = 16, MC + DMRG): "W ∝ t^(−1/3) on approach to the bulk critical temperature"; on "phenomenological grounds" one expects ω = 1/3.
    - Earlier q = 20 simulations had suggested ω "significantly larger than 1/3"; "the asymptotic regime is very narrow, and corrections to scaling are quite pronounced".
    - [Source](https://doi.org/10.1023/A:1004542105635)
  - **Fytas, Malakis, Selke & Shchur 2015**, "Interfacial adsorption in Potts models on the square lattice", Eur. Phys. J. B 88, 204, doi:10.1140/epjb/e2015-60326-8 [ABS]: at *continuous* transitions "the spread of the interfacial adsorption profiles is observed to increase linearly with the lattice size at the bulk transition point". [Source](https://doi.org/10.1140/epjb/e2015-60326-8)
- **Nonequilibrium wetting.** These concern a growing interface and a substrate, not an active slab between absorbing phases. Relevance: Medium.
  - **Hinrichsen, Livi, Mukamel & Politi 1997**, "Model for nonequilibrium wetting transitions in two dimensions", Phys. Rev. Lett. 79, 2710, doi:10.1103/PhysRevLett.79.2710 [ABS]: at the equilibrium point "complete wetting with critical exponents γ = 1/3 for the diverging interface width", with "a crossover to different exponents … related to a KPZ-type nonlinearity" away from it. [Source](https://doi.org/10.1103/PhysRevLett.79.2710)
  - **HLMP 2003**, "Wetting under nonequilibrium conditions", Phys. Rev. E 68, 041606, doi:10.1103/PhysRevE.68.041606 [ABS]. [Source](https://doi.org/10.1103/PhysRevE.68.041606)
  - HLMP 2000 (first order with coexistence) and de los Santos et al. 2003, see Q1 and Q3.
- **Activity on walls between absorbing domains** (the D&N analogue). Relevance: High.
  - Lee & Vojta 2011 [PT]: domain walls "can hop, branch and annihilate". [Source](https://doi.org/10.1103/PhysRevE.83.011114)
  - Lipowski–Droz 2002 [PT]: below Tc, "active sites are located mainly at the domain walls" and ρ ~ t^(−1/2) (curvature-driven coarsening). [Source](https://doi.org/10.1103/PhysRevE.65.056114)
  - Barghathi & Vojta 2012: domain walls between absorbing states perform Sinai walks under random fields. [Source](https://doi.org/10.1103/PhysRevLett.109.170603)
- **New diagnostics (preprints).**
  - **Lei & Ni 2026**, "Metastable hyperuniformity at discontinuous absorbing transitions", arXiv:2605.13561 [ABS]: S(k→0) ~ k^1.2 "appears only near the metastable regime", in both 2D and 3D. [Source](https://arxiv.org/abs/2605.13561)
  - **Lazarides & Pizzi 2026**, "Apparent bistability from weak long-range interactions", Phys. Rev. B 113, 224116, doi:10.1103/6c4l-kv48 [ABS]: R_c ~ h^(−1/(α−d)) makes "apparent bistability" possible with weak long-range couplings. Not applicable to short-range D&N. [Source](https://doi.org/10.1103/6c4l-kv48)

### Inferences
- **Signature 6 is generic as stated.** A finite-width active layer on a wall between two absorbing phases below the transition is expected in *any* Z2 absorbing model, continuous ones included (active domain walls, curvature-driven coarsening). What would make it a first-order signature, by analogy with q > 4 Potts complete wetting:
  - the slab interior carries the *bulk metastable active density* (a genuine third phase);
  - the width diverges as q → q*⁻ with a complete-wetting law, e.g. ∝ (q* − q)^(−1/3) if the 2D equilibrium analogy holds, logarithmic in mean field, possibly KPZ-modified out of equilibrium, rather than as the critical correlation length ξ⊥;
  - beyond q*, the edges move with exactly the front velocity v(q) of signature 1.
- I found no nonequilibrium theory fixing the exponent. A measurement would be new. This is another candidate "adjacent problem".
- **Signatures 4 and 5.** The ZGB precedent (no sharp critical size without diffusion) shows that a well-defined R_c is *not* automatic even at an accepted first-order point. D&N's R_c·v ≈ const is therefore informative, but the sharpness of the survival-vs-size curve (its width as a function of R) should be reported. R_c should also be checked for isotropy, or measured per orientation.
- Under the Allen–Cahn law, R_c·v ≈ const amounts to a measurement of the effective "kinetic line tension" D. Its non-vanishing at q* is precisely what Evans & Ray link to sharp nucleation. A vanishing D would signal a crossover toward criticality or a tricritical point.

### Gaps
- No paper found on interfacial adsorption or wetting by an active phase between two absorbing phases. No nonequilibrium prediction for the divergence exponent of such a slab.
- Selke & Pesch 1982 content not seen; Selke & Huse 1983 abstract not available.

## Q6. Methodology: accepted diagnostics for weak first order vs continuous, and known false positives and negatives

### Takeaway
Diagnostics accepted in this literature:

| Diagnostic | Typical source |
|---|---|
| Hysteresis / metastability loops | Constant-coverage (constant-concentration) ensembles |
| Superheating and supercooling limits (the latter via a vanishing field) | Lübeck 2006 |
| Bimodal QS order-parameter distributions with equal-area or susceptibility-maximum pseudo-transitions shifting ∝ 1/V | de Oliveira–da Luz–Fiore |
| Stripe and interface stationarity at several orientations | Evans group |
| Nucleation kinetics (KJMA, essential singularity, single- vs multidroplet) | Machado et al.; Rikvold et al. |
| Cluster mean-field convergence of the jump | Ódor–Szolnoki |

False positives come from long crossovers (CDP-like transients, the 1D TCM), from finite-size jumps (Lipowski–Droz q = 2), from homogeneous starts hiding orientation effects, and from steady-state exponent fits (BBC 1D). False negatives come from weak first order looking continuous ("notoriously difficult"). In Z2 two-absorbing-state models, Binder cumulants are unreliable for locating transitions.

### Cited Findings
- **de Oliveira & Dickman 2005**, "How to simulate the quasistationary state", Phys. Rev. E 71, 016129, doi:10.1103/PhysRevE.71.016129, arXiv:cond-mat/0407797 [ABS]: the QS method is "about an order of magnitude more efficient" at criticality. Relevance: High (tool). [Source](https://doi.org/10.1103/PhysRevE.71.016129)
- **de Oliveira, da Luz & Fiore 2015**, "Generic finite size scaling for discontinuous nonequilibrium phase transitions into absorbing states", Phys. Rev. E 92, 062126, doi:10.1103/PhysRevE.92.062126, arXiv:1510.08707 [ABS+PT].
  - "quantities such as response functions, cumulants, and equal area probability distributions all scale with the volume"; "|λ_V − λ₀| goes with 1/V".
  - Five lattice models, including the ZGB with a bimodal QS distribution.
  - Evidence: Moderate to Strong. Relevance: **High** (the standard FSS test D&N could add). [Source](https://doi.org/10.1103/PhysRevE.92.062126)
- **Ziff & Brosilow 1992** (constant coverage: metastability loop, spinodal, precise y₂), see Q2. [Source](https://doi.org/10.1103/PhysRevA.46.4630)
- **Fiore & de Oliveira 2007**, "Contact process with long-range interactions: A study in the ensemble of constant particle number", Phys. Rev. E 76, 041103, doi:10.1103/PhysRevE.76.041103 [ABS]: in the conserved ensemble "the presence of a 'loop' characterizes the first-order transition", and the cluster fractal dimension decides the order. [Source](https://doi.org/10.1103/PhysRevE.76.041103)
- **Lübeck 2006**: superheating limit with no systematic L dependence; supercooling via field h → 0; stripe probability that "decreases with increasing system size"; "notoriously difficult to distinguish … small but finite … β_t from a discontinuous jump". [Source](https://doi.org/10.1007/s10955-005-9000-z)
- **Ódor RMP 2004** [PT]: "By simulations it is quite difficult to decide whether a transition is really discontinuous. The order parameter of weak first order transitions – where the jump is small – may look very similar to continuous transitions. The hysteresis … is a demanding task to measure." He lists debated cases (Dickman & Tomé 1991; Hinrichsen 2000b, 2001a; Lipowski 1999; Lipowski & Lopata 1999; Szolnoki 2000; Tomé & de Oliveira 1989). [Source](https://doi.org/10.1103/RevModPhys.76.663)
- **Machado et al. 2005** (FSS of fluctuations and the fourth-order cumulant; KJMA lifetimes). **Rikvold et al. 1994** (single vs multidroplet). **Wang–Liu–Evans 2015** (Avrami fits, τ_nuc scaling, constant-concentration strips per orientation, constant-concentration-stabilised critical droplets). See Q5.
- **Ódor & Szolnoki 1996**: generalised mean-field cluster methods because "it is not easy to prove first order transition by simulation". [Source](https://doi.org/10.1103/PhysRevE.53.2231)
- **Windus & Jensen 2007**: density histograms showing coexistence. **Fiore 2014**: bimodal distributions. [Source](https://doi.org/10.1088/1751-8113/40/10/005)
- **Known false positives and controversies.**
  - 1D TCM: Dickman–Tomé 1991 → Hinrichsen 2000, Ódor–Dickman 2009, Park 2009. [Source](https://doi.org/10.1103/PhysRevE.80.061103)
  - BBC 1D non-DP exponents → Jensen 1991, DP. [Source](https://doi.org/10.1103/PhysRevA.43.3187)
  - Threshold CP: "critical discontinuous" transition at a unique point (da Silva & de Oliveira 2011) vs GTPC with finite line tension (Wang–Liu–Evans 2015). [Source](https://www.osti.gov/biblio/1227298)
  - 2D q = 2 nonequilibrium Potts: jump possibly finite-size (Lipowski–Droz 2002). [Source](https://doi.org/10.1103/PhysRevE.65.056114)
  - 2D interacting monomer–dimer: non-Ising claim (Nam et al. 2011) vs Ising (Park 2012, 2015), with "the Binder cumulant is not an efficient quantity to locate the order-disorder transition point of the model with two symmetric absorbing states". [Source](https://doi.org/10.1088/1742-5468/2015/10/P10009)
  - Mean-field first order that is continuous in 2D: symbiotic CP (Sampaio Filho et al. 2018); Allee model at strong noise or low diffusion (PNAS 2015). [Source](https://doi.org/10.1103/PhysRevE.98.062108)
- **Robustness test via absorbing-breaking noise.** In ZGB (CO desorption) and in the Schlögl models (spontaneous creation ε), the first-order line or GTPC region ends at an Ising-class critical point (Tomé–Dickman 1993; Brosilow–Ziff 1992; Liu–Shen–Evans 2025). [Source](https://doi.org/10.1103/PhysRevE.47.948); [Source](https://doi.org/10.1103/y54y-5dvx)

### Inferences
- In a Z2 two-absorbing-state model the extra hazard is the 2D GV class, which is marginal (d_c = 2) with logarithmic decays: ρ ~ 1/ln t, and spreading δ ≈ 0.9 with η ≈ 0, as in Lipowski–Droz. Slow log decays can masquerade as plateaus. D&N's density data should therefore be shown as 1/ρ vs ln t at q*, alongside the QS-ensemble bimodality and 1/L² shift tests.
- The cleanest accepted closure is that three independent estimates coincide within errors, extrapolated as L^(−2): the front-stationarity point (for each orientation), the QS equal-area or susceptibility-maximum point, and the seed-survival threshold. A gap between them signals GTPC; continuous vanishing of the jump signals approach to criticality.

### Gaps
- I found no accepted single "smoking gun"; the field relies on convergence of several diagnostics.
- Not read: Henkel–Hinrichsen–Lübeck (2008) and Marro–Dickman (1999) books; Lübeck's IJMPB 18, 3977 (2004) review (doi:10.1142/S0217979204027748, citation only).

## Q7. Assessment: which of the project's six signatures this literature treats as conclusive, and which known pitfalls may be missing

### Takeaway
Signature 1 (front-velocity sign change) and signature 4 (critical droplets with R_c·v ≈ const) are the kind of evidence that settled the 2D ZGB, QCP and TDP cases. They are strong, **except that signature 1 was measured for only one orientation**. In this literature that is the best-known way such a claim can be incomplete: the transition can turn out to be a GTPC window. Signature 2 is standard but needs the QS/FSS and log-decay checks. Signatures 3, 5 and 6, as stated, are *consistent with* first order but do not discriminate it from a continuous or GV transition. Each becomes discriminating with a specific additional measurement.

### Cited Findings
1. **Front velocity changes sign linearly at q*, one orientation only.**
   - Planar-interface stationarity is the accepted nonequilibrium "Maxwell" criterion (Lübeck 2006; Grassberger 2006; Evans & Ray 1994; Liu–Wang–Evans 2018: "phase coexistence or equistability is usually determined by stationarity of a planar interface"). [Source](https://doi.org/10.1103/PhysRevLett.121.120603)
   - Missing: the orientation scan. Square-lattice NN rules show 5–11% windows (Wang–Liu–Evans 2015; Guo 2009); the Moore-neighbourhood threshold model shows none detectable (Liu–Shen–Evans 2025). [Source](https://www.osti.gov/biblio/1227298); [Source](https://doi.org/10.1103/y54y-5dvx)
   - Also check: no zero-velocity plateau. Pinning or "renormalisation to zero" is seen in 1D (Weissmann–Shnerb) and propagation failure in discrete deterministic equations (Wang–Liu–Evans 2012, 2020); stochastic fluctuations usually remove it. [Source](https://doi.org/10.1209/0295-5075/106/28004)
   - Also check: front roughness and KPZ coefficients, since in ZGB the kinetic surface tension vanishes at the transition (Evans & Ray 1994). [Source](https://doi.org/10.1103/PhysRevE.50.4302)
   - Rating: strong, **conditional on isotropy**.
2. **Finite, L-independent interface density down to q*.**
   - Accepted when combined with metastability or hysteresis and FSS. Superheating limits show no systematic L shift (Lübeck); QS pseudo-points should shift ∝ 1/V (de Oliveira–da Luz–Fiore 2015). [Source](https://doi.org/10.1103/PhysRevE.92.062126)
   - Pitfalls: weak first order vs small β (Lübeck; Ódor); GV logarithmic decays; the 2D q = 2 precedent judged inconclusive on exactly this kind of data (Lipowski–Droz 2002). [Source](https://doi.org/10.1103/PhysRevE.65.056114)
   - Rating: supportive.
3. **L-independent metastable lifetimes below q*.**
   - L-independence is the *multidroplet* regime (VBM14; Machado et al. 2005; Rikvold et al. 1994). A continuous transition also gives L-independent decay times for L ≫ ξ⊥, so L-independence alone does not discriminate.
   - Discriminating forms: ln τ ∝ 1/Δq, via τ ≈ b δp^(−2/3) exp[c/(3δp)] (Wang–Liu–Evans 2015); an Avrami (t³) decay shape; a crossover to single-droplet (L-dependent) statistics as Δq → 0 at fixed L, moving logarithmically with L. [Source](https://doi.org/10.1103/PhysRevE.71.031603); [Source](https://doi.org/10.1103/PhysRevE.49.5080)
   - Rating: weak as stated, strong with these fits.
4. **Critical droplets with R_c·v ≈ const.**
   - Classical-nucleation-consistent (R_c ~ σ_eff/δp). It is the "finite line tension" evidence Wang–Liu–Evans used against "critical discontinuous" claims, and it directly demonstrates Hinrichsen's d ≥ 2 island-elimination mechanism. [Source](https://www.osti.gov/biblio/1227298)
   - Pitfalls: R_c may not be sharp in strongly fluctuating systems (Evans & Miesch 1991; sharpening with reduced fluctuations, Evans & Ray 1994). [Source](https://doi.org/10.1103/PhysRevLett.66.833)
   - Under GTPC, droplets facet and shrink at constant speed (Wang–Liu–Evans 2015; Muñoz et al. 2005). [Source](https://doi.org/10.1140/epjb/e2005-00029-3)
   - Rating: strong if the survival curve is sharp and droplets stay rounded.
5. **2×2 seeds decay faster than voter seeds, then plateau (P∞ > 0).**
   - P∞ > 0 above threshold holds for continuous transitions too, and near q* small seeds survive rarely in both scenarios.
   - Informative comparisons: early-time decay against the voter law (δ ≈ 0.9–1, η ≈ 0, Lipowski–Droz); and the **seed-survival threshold q_f vs the front threshold q*(S)**. The Durrett postulate p_f = p_e holds only without GTPC (Guo et al. 2009). [Source](https://doi.org/10.1063/1.3074308)
   - Rating: weak as stated, strong as a q_f-vs-q* consistency test.
6. **Wall between the absorbing phases wetted by a finite active slab below q*; edges become moving fronts above.**
   - Active domain walls between absorbing domains are generic in Z2 models (Lee–Vojta 2011; Lipowski–Droz 2002), so the existence of the slab is not diagnostic. [Source](https://doi.org/10.1103/PhysRevE.83.011114)
   - Discriminating form: the q > 4 Potts analogue, complete wetting by the third phase with W ∝ t^(−1/3) at first-order Tc, against adsorption growing ∝ L at continuous Tc (Carlon et al. 1999; Fytas et al. 2015). D&N would need the slab-width divergence law and a bulk-like interior density. [Source](https://doi.org/10.1023/A:1004542105635); [Source](https://doi.org/10.1140/epjb/e2015-60326-8)
   - Rating: suggestive; potentially novel if quantified.
- **Pitfalls not listed among the six signatures:**
  - **(a) Orientation dependence and GTPC** (Q3), including that homogeneous starts average over orientations ("obscured starting with a random configuration"). [Source](https://www.osti.gov/biblio/1227298)
  - **(b) Fluctuation-induced softening near the voter point.** A tricritical point on the q*(p) line is plausible: 2D first order destroyed by strong noise or low mobility (PNAS 2015); 2D TDP first-order lines ending at finite parameters (Lübeck 2006). Track the jump magnitude along the line. [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC4403222/)
  - **(c) Robustness to small symmetric absorbing-breaking noise.** In ZGB and Schlögl models the first-order line ends at an Ising-class critical point (Tomé–Dickman 1993; Liu–Shen–Evans 2025). For Z2 models compare the "Modified GV" line in Llabrés et al. 2026 (LITERATURE.md). [Source](https://doi.org/10.1103/PhysRevE.47.948)
  - **(d) Binder cumulants** are unreliable for locating transitions in two-symmetric-absorbing-state models (Park 2015). [Source](https://doi.org/10.1088/1742-5468/2015/10/P10009)
  - **(e) Quenched-disorder robustness** if any physical-relevance claims are made (VBM14). [Source](https://doi.org/10.1103/PhysRevE.89.012145)
  - **(f) Novelty phrasing** must cite Lipowski–Droz 2002 (2D q = 3 first order; q = 2 inconclusive) as well as the LITERATURE.md items. [Source](https://doi.org/10.1103/PhysRevE.65.056114)

### Inferences
- **Narrow option.** "First microscopic demonstration of a *fully* first-order active/absorbing transition with two Z2-symmetric absorbing states in 2D" is defensible and fills the q = 2, d = 2 gap left by Lipowski–Droz. It needs (i) fronts at ≥ 2–3 orientations with matching q*, (ii) QS/FSS bimodality with 1/L² shifts, and (iii) lifetime fits showing an essential singularity. Without (i), a referee familiar with the Evans-group work could argue the transition is a GTPC window.
- **Richer adjacent problems with no precedent found:**
  - GTPC (orientation-dependent equistability) in Z2-symmetric absorbing systems, where droplets of both absorbing phases and the active phase interplay;
  - nonequilibrium interfacial adsorption and complete wetting of absorbing–absorbing walls by the active phase, with the divergence law as a new observable;
  - quenched random-field vs random-bond disorder at a Z2 first-order absorbing transition, between the VBM14 and Barghathi–Vojta mechanisms;
  - the tricritical end of the line approaching the voter (GV, d_c = 2) point.
- The first two reuse D&N's existing front and slab machinery directly.

### Gaps
- The assessment relies on the project summary provided. I did not inspect the D&N data, so whether the six measurements already include the extra checks (orientation, FSS shifts, Avrami fits, slab-width law) is unknown.
- No literature found on synchronous-versus-asynchronous effects on orientation dependence.
