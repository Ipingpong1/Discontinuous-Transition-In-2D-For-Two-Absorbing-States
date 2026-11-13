# Rigorous probability literature on clustering, coexistence and bistability in two-type spin systems with two traps (emphasis d = 2)

**Scope.** These notes cover 1975 to October 2026, including arXiv preprints.

**Reading-level codes:**
- **FT** — full text (sections named), read in this research, by me or by a delegated sub-agent.
- **ABS** — abstract only.
- **REV** — zbMATH or MathSciNet review.
- **SEC** — known only through a citing paper.

**Relevance:** H / M / L.

**Threat** means the source already proves or conjectures a first-order clustering/coexistence transition in a 2D two-trap Z2 model, or computes the bistable window near the 2D voter point.

**Access notes:**
- Local full texts are in `/home/matias/Desktop/silly/dn/campaign/theory/lit/txt/`.
- Project Euclid blocked automated access. Several Ann. Probab. / AAP primary texts were therefore read only at ABS/REV level; this is marked on each entry.
- Quotes are verbatim unless marked as paraphrase.

## 1. Nonlinear and threshold voter models: what is proved about clustering vs coexistence, and about the type of transition?

### Takeaway
In d = 2 the rigorous theory of threshold and nonlinear voter models covers only the two ends of the family:
- **Threshold 1, any range:** coexistence, with a complete convergence theorem (Liggett 1994; Handjani 1999).
- **Thresholds above about half the neighbourhood:** fixation (Durrett–Steif 1993).
- **In between:** nothing rigorous.
- **The 2D Molofsky family:** known only from simulation.
  - The clustering/coexistence boundary runs through the voter point, close to the "Levin/dyad line".
  - The spatial model shows no "multiple outcome" regime.

No source proves or conjectures a discontinuous clustering/coexistence transition. The one careful numerical study of the type of such a transition in this class (the 1D rebellious voter model, Swart–Vrbenský 2010) finds it continuous.

### Cited Findings

#### Cox & Durrett (1991), "Nonlinear voter models"
- **Citation:** in R. Durrett and H. Kesten (eds.), *Random Walks, Brownian Motion, and Interacting Particle Systems* (Festschrift for F. Spitzer), Progress in Probability 28, Birkhäuser; MR1146446; zbMATH 0825.60053. — [zbMATH](https://zbmath.org/?q=an:0825.60053)
- **Page range — CONFLICT:**
  - 189–201 per zbMATH and per Sturm–Swart 2008, ref. [7]. — [Sturm–Swart](https://arxiv.org/abs/math/0701555)
  - 203–213 per Cox–Durrett–Perkins 2013, ref. [4], and Cox–Perkins 2023, ref. [5]. — [CDP](https://arxiv.org/abs/1103.1676); [CP](https://arxiv.org/abs/2310.00812)
  - The project's LITERATURE.md has 189–202.
- **Read:** not read; there is no online copy, and the zbMATH review text is withheld. Content below is SEC.
- **Model:** nonlinear voter models, where the flip rate a_ℓ depends only on the number ℓ of disagreeing neighbours in N. The threshold voter model is a special case.
  - Cox–Perkins 2023 define the class "As defined in [5]" and take from it (Prop. 3.2) the matrix criterion for a nonlinear voter model to be cancellative. — [CP 2023](https://arxiv.org/abs/2310.00812)
  - Cox–Durrett–Perkins: "General models of this type were introduced and studied in [4]." — [CDP](https://arxiv.org/abs/1103.1676)
- **Statements, as reported by citing papers:**
  - Liggett 1994: "Cox and Durrett conjectured that it [threshold-1 voter] coexists in all cases except N = d = 1." — [Liggett 1994](https://doi.org/10.1214/aop/1176988729)
  - Per the zbMATH review of Handjani 1999, Cox–Durrett also conjectured complete convergence. — [Handjani](https://doi.org/10.1214/aop/1022677260)
- **UNVERIFIED:** the project note that this paper "introduces the 5-site discrete family and calls the Levin line 3p₁+p₂ = 1 the empirical boundary". The verified sources for that line are Molofsky et al. 1999 (dyad rule), Durrett's web survey and CDP 2013, Fig. 3 (all below).
- **Relevance:** M–H (origin of the class and of the cancellativity criterion). **Threat:** none known.

#### Molofsky, Durrett, Dushoff, Griffeath & Levin (1999), "Local frequency dependence and global coexistence"
- **Citation:** Theor. Popul. Biol. 55, 270–282; doi:10.1006/tpbi.1998.1404. — [DOI](https://doi.org/10.1006/tpbi.1998.1404)
- **Read:** FT (local PDF; §§2–4 re-checked for this note).
- **Model:** synchronous PCA on Z², with a 5-cell neighbourhood (site plus 4 nearest neighbours).
  - P(ξₙ₊₁(x) = 1) = p_k when k cells of the neighbourhood are 1.
  - p₀ = 0 and p₅ = 1 (two traps); p₁ = 1 − p₄ and p₂ = 1 − p₃ (Z2 symmetry).
  - The voter point is (0.2, 0.4).
- **Simulation results:**
  - "In every instance, however, small isolated configurations of species 1—such as the dyads described in the next section—have a propensity to nucleate linearly growing patches of the symmetric equilibrium. We conclude that, in contrast to the mean-field model, there is no multiple outcome regime of the spatial model." — [DOI](https://doi.org/10.1006/tpbi.1998.1404)
  - "the boundary separating the clustering and ergodic regions runs from about (0.024, 1) at the top edge, down through the linear voter point (0.2, 0.4), before a bifurcation encloses the tiny phase separation." — [DOI](https://doi.org/10.1006/tpbi.1998.1404)
  - "the multiple outcome region (found when p1 is small and p2 is relatively large) disappears in the explicit spatial model, making coexistence possible under a wider range of conditions than is found in the mean-field approximation." — [DOI](https://doi.org/10.1006/tpbi.1998.1404)
  - "if the monad rule is violated but the dyad rule is satisfied (e.g., if p1 = 0.1 and p2 = 0.5), most introductions will fail. A few isolated monads, however, stochastically will duplicate, getting beyond the initial threshold and allowing the dyad rule to become operational." — [DOI](https://doi.org/10.1006/tpbi.1998.1404)
  - At (0.31, 0): "two distinct stable mixed equilibria dominated by 0's and 1's, respectively" (phase separation). — [DOI](https://doi.org/10.1006/tpbi.1998.1404)
- **Not discussed:** the order of the clustering/ergodic transition. A keyword search for "first order", "discontinu", "abrupt", "jump" and "hysteres" finds nothing relevant.
- **Relevance:** H. This is the direct 2D precedent of the same model class.
- **Threat:** L–M. It reports the opposite of bistability: nucleation from dyads, and no multiple-outcome regime.

#### Durrett, web survey page "Nonlinear Voter Model" (undated)
- **Read:** FT (sub-agent).
- "Empirically the boundary between coexistence and clustering is given by the Levin line: 3 p[1] + p[2] = 1." — [Durrett survey](https://sites.math.duke.edu/~rtd/survey/survc5.html)
- Open question on the "Colorado Springs" model (p₁ = 0.3, p₂ = 0): "It would be very interesting to show that a small perturbation of this model… has two translation invariant stationary distribution." — [Durrett survey](https://sites.math.duke.edu/~rtd/survey/survc5.html)
- **Relevance:** M. **Threat:** L.

#### Liggett (1994), "Coexistence in threshold voter models"
- **Citation:** Ann. Probab. 22(2), 764–802; doi:10.1214/aop/1176988729. Correction: Ann. Probab. 27 (1999), 2120–2121 (per sub-agent). — [DOI](https://doi.org/10.1214/aop/1176988729)
- **Read:** ABS.
- **Model:** threshold-1 voter model on Z^d, range N. Whether the range uses the sup-norm, which would make N = 1 in d = 2 the Moore neighbourhood, is UNVERIFIED.
- **Theorem:** "This system is known to cluster if N = d = 1, and to coexist if N ≥ 4 in one dimension and if N is reasonably large in other dimensions. Cox and Durrett conjectured that it coexists in all cases except N = d = 1. In this paper, we prove this conjecture… part of the proof is computer aided." — [DOI](https://doi.org/10.1214/aop/1176988729)
- **Relevance:** H (the coexisting end of the family in d = 2). **Threat:** none.

#### Durrett & Steif (1993), "Fixation results for threshold voter systems"
- **Citation:** Ann. Probab. 21(1), 232–247; doi:10.1214/aop/1176989403. — [DOI](https://doi.org/10.1214/aop/1176989403)
- **Read:** ABS and REV.
- **Theorem (per review):** "if 0∈N and τ>(|N|−1)/2, then starting from any initial configuration, lim ξ_t(x) exists (i.e., the system 'fixates' at x) almost surely for all x." — [DOI](https://doi.org/10.1214/aop/1176989403)
- **Picture from the abstract:** "There is a critical value, θ_c, so that if τ = θn with θ > θ_c and n is large then most sites never flip, while for θ ∈ (1/2, θ_c) and n large, the limiting state consists mostly of large regions of points of the same type. In d = 1, θ_c ≈ 0.6469076 while in d > 1, θ_c = 3/4." How much of the d > 1 part is proved is UNVERIFIED. — [DOI](https://doi.org/10.1214/aop/1176989403)
- **Relevance:** H (the fixating end). **Threat:** none.

#### Andjel, Liggett & Mountford (1992), "Clustering in one-dimensional threshold voter models"
- **Citation:** Stoch. Proc. Appl. 42, 73–90; doi:10.1016/0304-4149(92)90027-N. — [DOI](https://doi.org/10.1016/0304-4149(92)90027-N)
- **Read:** ABS.
- **Result:** in d = 1 with threshold equal to half the neighbourhood, "We prove that clustering occurs for all N ≥ 1". — [DOI](https://doi.org/10.1016/0304-4149(92)90027-N)
- **Relevance:** M. **Threat:** none.

#### Handjani (1999), "The complete convergence theorem for coexistent threshold voter models"
- **Citation:** Ann. Probab. 27(1), 226–245; doi:10.1214/aop/1022677260. — [DOI](https://doi.org/10.1214/aop/1022677260)
- **Read:** ABS and REV.
- **Theorem:** "in these coexistent cases, η_t ⇒ αδ_0 + βδ_1 + (1 − α − β)η_∞^{1/2}", with α = P(τ₀ < ∞) and β = P(τ₁ < ∞). This holds in every d ≥ 2 and every range. It is the template "two traps plus one non-trivial law" later used by Cox–Perkins. — [DOI](https://doi.org/10.1214/aop/1022677260)
- **Relevance:** M–H. **Threat:** none.

#### Not read
- Liggett, "Clustering and coexistence in threshold voter models", NATO ASI Ser. C 396 (Kluwer 1993), 403–410; zbMATH 0869.60086. — [zbMATH](https://zbmath.org/?q=an:0869.60086)
- Liggett, *Stochastic Interacting Systems* (Springer 1999), threshold-voter chapter. Its notes on thresholds between 1 and n/2 were not accessed.
- Both may hold statements about intermediate thresholds.

#### Dayue Chen
- No paper on threshold voter models was found in zbMATH author searches. The closest is "The consensus times of the majority vote process on a torus", J. Stat. Phys. 86 (1997) 779–802. — [DOI](https://doi.org/10.1007/BF02199120)
- Related: X. Xue, threshold voter models on trees (JSP 146, 2012) and on tori (JSP 159, 2015); titles only (sub-agent).
- **Relevance:** L.

#### Steif (1994), "The threshold voter automaton at a critical point"
- **Citation:** Ann. Probab. 22(3), 1121–1139.
- **Read:** ABS.
- **Content:** a deterministic automaton in d = 1. "in some sense, there is a discontinuity both from the left and from the right at this critical value". This is a transition between fixation regimes, not clustering/coexistence. — [DOI](https://doi.org/10.1214/aop/1176988597)
- **Relevance:** L.

#### Sowers (1999)
- **Citation:** J. Funct. Anal. 169, 421–455.
- **Read:** REV.
- **Content:** on Z², a threshold voter model in which a site flips if at least 2 of its 4 nearest neighbours disagree. Under diffusive rescaling the interfaces converge to mean-curvature flow. This is a rigorous 2D coarsening statement for a symmetric majority-type two-trap rule. — [DOI](https://doi.org/10.1006/jfan.1999.3477)
- **Relevance:** L–M.

#### Lanchier and co-authors (rigorous lattice models with mean-field bistability)
- **Lanchier & Neufer (2013), "Stochastic dynamics on hypergraphs and the spatial majority rule model"** — J. Stat. Phys. 151, 21–45; read FT §§1–2.
  - Theorem 4: "Assume that d = 2 and n = 2. Then, opinion 1 wins" (ties go to type 1).
  - Conjecture 1: "Clustering occurs when n is odd".
  - Heuristic: majority-type rules behave "similar to that of threshold voter models with threshold parameter slightly smaller than half of the number of neighbors". Simulations "suggest that… interface dynamics follow motion by mean curvature." — [arXiv:1301.0151](https://arxiv.org/abs/1301.0151)
- **Lanchier (2015), "Evolutionary games on the lattice: payoffs affecting birth and death rates"** — Ann. Appl. Probab. 25(3), 1108–1154; read ABS and introduction.
  - Proves "the presence of a dominant strategy that wins even when starting at arbitrarily low density in the region where the replicator equation displays bistability". — [arXiv:1302.0069](https://arxiv.org/abs/1302.0069)
- **Evilsizor & Lanchier, best-response dynamics** — EJP 19 (2014); doi 10.1214/EJP.v19-3126.
  - Under best-response dynamics the favoured strategy wins in every d. — [arXiv:1310.5192](https://arxiv.org/abs/1310.5192)
- **Evilsizor & Lanchier, death-birth updating** — EJP 21 (2016); doi 10.1214/16-EJP4380. — [arXiv:1506.07233](https://arxiv.org/abs/1506.07233)
- **Foxall & Lanchier, "…death and birth of the fittest"** — ALEA 14 (2017).
  - In the doubly-stable region "either one strategy wins or the system fixates to a configuration where both strategies are present". — [arXiv:1605.04037](https://arxiv.org/abs/1605.04037)
- **Overall:** asymmetric games, not Z2-symmetric two-trap models near the voter point. **Relevance:** M. **Threat:** none.
- **Correction to the assignment's list:** the death-birth paper is Evilsizor–Lanchier (EJP 21, 2016), not Foxall–Lanchier, and the best-response paper is EJP 19 (2014).

#### Sturm & Swart (2008), "Voter models with heterozygosity selection"
- **Citation:** Ann. Appl. Probab. 18(1), 59–99; doi:10.1214/07-AAP444; arXiv:math/0701555.
- **Read:** FT §§1–2.5 (local PDF).
- **Models:** neutral Neuhauser–Pacala model; affine voter model (voter mixed with threshold voter); 1D rebellious voter model.
- "It is believed, and has been proved in special cases, that coexistence occurs for sufficiently small competition parameters α or in high dimensions, while noncoexistence occurs in low dimensions and for large competition parameters." — [arXiv](https://arxiv.org/abs/math/0701555)
- **Lemma 1 (duality):** coexistence ⇔ the dual parity-preserving system survives; survival ⇔ the dual persists. — [arXiv](https://arxiv.org/abs/math/0701555)
- **Open question Q1:** "Is coexistence of X equivalent to survival of X?" — [arXiv](https://arxiv.org/abs/math/0701555)
- **Discussion:** "Very little is known about noncoexistence, except for the pure voter model in dimensions d = 1, 2 and the disagreement voter model … noncoexistence can be expected in dimension one only. In fact, it seems that physicists believe that for these and similar models, there is a critical dimension dc ≈ 4/3 such that only below dc there is a nontrivial phase transition between coexistence and noncoexistence." — [arXiv](https://arxiv.org/abs/math/0701555)
- **Numerics:** "For the rebellious voter model, numerical simulations suggest the existence of a critical value αc ≈ 0.5". — [arXiv](https://arxiv.org/abs/math/0701555)
- **Relevance:** M–H. **Threat:** none.

#### Swart & Vrbenský (2010), "Numerical analysis of the rebellious voter model"
- **Citation:** J. Stat. Phys. 140(5), 873–899; doi:10.1007/s10955-010-0021-x; arXiv:0911.1266.
- **Read:** FT of v3 (abstract, §1, §2.1, §3.3, §5), by sub-agent.
- **Model:** d = 1 only.
- **Simulation results:**
  - α_c = 0.510 ± 0.002 for the two-sided model; the one-sided model appears to have α_c = 1/2. — [arXiv](https://arxiv.org/abs/0911.1266)
  - "Our data strongly suggest that ρ(αc) = 0 = χ(αc), which implies that at criticality, the process exhibits neither coexistence nor interface stability." The proposed formula "implies that the order parameter critical exponent β for this model is one". So the transition is continuous. — [arXiv](https://arxiv.org/abs/0911.1266)
  - On the Neuhauser–Pacala model: "For d = 2, their model is supposed to exhibit noncoexistence only for α = 1". Also: "there exists no rigorous proof of noncoexistence… except in the trivial cases α = 1 (pure voter) or d = 1, N = 1." — [arXiv](https://arxiv.org/abs/0911.1266)
- **Relevance:** M. **Threat:** none (1D, continuous).

#### Sturm & Swart (2015), "A particle system with cooperative branching and coalescence"
- **Citation:** Ann. Appl. Probab. 25(3), 1616–1649; arXiv:1311.0417.
- **Read:** ABS and §§1.3–1.5.
- **Model:** d = 1 only.
- **Order of the transition:** not proved. First-order behaviour holds only for "the corresponding mean-field model".
- **On d ≥ 2:** minority-favouring models "in dimensions d ≥ 2 tend to have an invariant law in which both types are present for all values of the selection parameter, but in dimension one undergo a phase transition." — [arXiv](https://arxiv.org/abs/1311.0417)
- **Relevance:** M–L.

#### Swart, *A Course in Interacting Particle Systems* (lecture notes)
- **Citation:** arXiv:1703.10007 (v5, Nov 2025).
- **Read:** FT of the relevant sections (sub-agent).
- "It can be difficult to prove that a given phase transition is first or second order." — [arXiv](https://arxiv.org/abs/1703.10007)
- Exercise 3.3 asks, for the mean-field cooperative-branching model: "Is it first- or second order?" — [arXiv](https://arxiv.org/abs/1703.10007)
- **Relevance:** M.

#### Cox & Durrett (1995), "Hybrid zones and voter model interfaces"
- **Citation:** Bernoulli 1(4), 343–370.
- **Read:** ABS.
- **Content:** "We study the dynamics of hybrid zones in the absence of selection." This is the pure voter model, not a majority perturbation. — [DOI](https://doi.org/10.2307/3318488)
- **Relevance:** L.

### Inferences
- **Two endpoints, no middle.** Between threshold 1 (coexistence) and "more than half" (fixation), there is no 2D theorem. Lanchier–Neufer and Sowers suggest that majority-type rules near one half coarsen by mean curvature, which is the clustering phenomenology of a symmetric bistable rule.
- **Theory does not reach the D&N rule far from the voter point.** The project's rule interpolates between the voter model (p = 1) and Day & Night (p = 0), and D&N is non-monotone (flips at j = 3 but not at j = 4, 5).
- **Molofsky et al. is the main empirical counterpoint to 2D bistability in this class.** Their 2D boundary follows a small-seed criterion: the "dyad rule", the expected offspring of two adjacent 1's.
  - A boundary fixed by a small-seed criterion is what one expects from a continuous, seed-controlled transition, not from a front-speed (Maxwell) line.
  - The project's critical droplets (seeds below R_c die) are therefore a genuinely different 2D observation.
  - The project should state explicitly why D&N/Moore differs from the 5-cell rule.
  - It should also quantify the survival probability against seed size, since Molofsky saw "most introductions will fail" while a few monads succeed. Positive survival probability from small seeds is compatible with both pictures.
- **The only quantitative study of the transition type in this class is 1D and continuous** (Swart–Vrbenský). A 2D analogue study does not exist in the probability literature.

### Gaps
- Cox & Durrett 1991 is unread, and its page range conflicts between sources.
- Not read: Liggett's 1999 book chapter and the 1993 NATO paper (intermediate thresholds); Durrett–Steif beyond abstract and review.
- How much of the Durrett–Steif d > 1 "θ_c = 3/4" picture is proved is unknown.
- No Dayue Chen threshold-voter paper was found. The assignment's reference may be mistaken.

## 2. Voter-model perturbations: what is proved, and how does the project's calculation relate to it?

### Takeaway
**In d ≥ 3** there is a reaction–diffusion limit, and its reaction term classifies the phases (Cox–Durrett–Perkins 2013, "CDP").
- Coexistence is proved even in the tristable **Case 4A**: f′(0) < 0, but the 1/2-phase invades the traps.
- The symmetric bistable **Cases 3 and 4B** are explicitly left open: "our method yields no result".

**In d = 2** there is only a low-density super-Brownian limit with drift Θ₃, and coexistence/survival when Θ₃ > 0 (Cox–Perkins 2023/2025). In their words the reaction function "is not well-defined for d = 2", and the founder-control side is beyond "currently available" tools even in d ≥ 3.

No paper from 2010–2026 proves, disproves or explicitly conjectures a first-order clustering/coexistence transition or a bistable window near the 2D voter point.

There is conjecture-level prior art for the **structure**: a 4A/4B coexistence/clustering line, set by the front speed and running from the 2D voter point.
- CDP 2013, Fig. 3.
- Huang–Durrett 2021.
- Durrett's 2026 book draft, which states "we have clustering" in 4B.
- None of them calls that line discontinuous.

### Cited Findings

#### Cox, Durrett & Perkins (2013), "Voter model perturbations and reaction diffusion equations"
- **Citation:** Astérisque 349, vi+113 pp.; ISBN 978-2-85629-355-3; Zbl 1277.60004; MR3075759; arXiv:1103.1676. — [arXiv](https://arxiv.org/abs/1103.1676); [Numdam](https://www.numdam.org/item/AST_2013__349__R1_0/)
- **Read:** FT §§1.1–1.5 and §1.7.3 (local).
- **Dimension:** "throughout this work we assume d ≥ 3". The abstract: "converges to a solution of a reaction diffusion equation in dimensions d ≥ 3 … An important feature of our general results is that they do not require the process to be attractive." — [arXiv](https://arxiv.org/abs/1103.1676)
- **Nonlinear voter model (§1.5).** The model is a continuous-time version of Molofsky's, with 4 neighbours drawn from range L.
  - Reaction term: f₁ = b₁u(1−u)⁴ + b₂u²(1−u)³ − b₂u³(1−u)² − b₁u⁴(1−u), with b₁ = 4a(1) − a(4) and b₂ = 6a(2) − 4a(3).
  - f₁′(0) = b₁; f₁′(1/2) = −(6b₁ + 2b₂)/16; ∫₀^{1/2} f₁ = (5b₁ + b₂)/192. — [arXiv](https://arxiv.org/abs/1103.1676)
- **Case 3 (0 and 1 attracting, 1/2 unstable):** "In this bistable case the winner is dictated by the sign of the speed of the traveling wave, but by symmetry (recall (1.53)) the speed is 0. One would guess that clustering occurs in this case and there are only trivial stationary distributions, but our method yields no result." — [arXiv](https://arxiv.org/abs/1103.1676)
- **Case 4 (0, 1/2 and 1 all attracting):** "Symmetry implies c2 = −c1, but we can have c1 < 0 < c2 (Case 4A) … Hence we expect coexistence in Case 4A and all invariant distributions to have density near 1/2. If c1 > 0 > c2 (Case 4B) and L is large enough, there is a standing wave solution w0(x) of the PDE in d = 1 with w0(−∞) = 0, w0(∞) = 1 … and our method yields no result." — [arXiv](https://arxiv.org/abs/1103.1676)
- **Theorem 1.13 (d ≥ 3, L large, ε small):**
  - (a) Coexistence in Cases 1, 2 and 4A.
  - (b) In Cases 1 and 4A every stationary law ν with ν(ξ ≡ 0 or ξ ≡ 1) = 0 satisfies sup_x |ν(ξ(x)=1) − 1/2| ≤ η. — [arXiv](https://arxiv.org/abs/1103.1676)
- **Key remark:** "Case 4A is of particular interest as there is coexistence even though f′(0) < 0. Here the low density limit theorem in [7] shows convergence to super-Brownian motion with drift f′(0) < 0 … From this one might incorrectly guess (after an exchange of limits) that there is a.s. extinction of 0's for the nonlinear voter model, while our proof of coexistence will show that there is positive probability of survival of 0's even starting with a single 0." — [arXiv](https://arxiv.org/abs/1103.1676)
- **Case 2 (two interior stable states):** "Simulations in [36] … suggest that in this case there may be two nontrivial stationary distributions". — [arXiv](https://arxiv.org/abs/1103.1676)
- **Figure 3**, captioned "Conjectured phase diagram for the discrete time two-dimensional nonlinear voter model of [36]":
  - It has regions 1, 2 ("two equil.?"), 3 ("clustering"), 4A ("coexistence"), 4B and "periodic". The voter model sits at (0.2, 0.4).
  - Text: "There we will see that the division between 4A and 4B is dictated by the speed of traveling waves for the PDE. Here we have drawn the 'Levin line' 6p1 + 2p2 = 2 … Simulations suggest that the true boundary curve exits the square at (0.024, 1)". — [arXiv](https://arxiv.org/abs/1103.1676)
- **Figure 5**, "Phase diagram for the continuous time nonlinear voter model with large range in d ≥ 3": the lines b₁ = 0, 3b₁ + b₂ = 0 and 5b₁ + b₂ = 0 all emanate from the voter point. — [arXiv](https://arxiv.org/abs/1103.1676)
- The words "first order", "discontinu" and "tricritical" do not occur.
- **Relevance:** H.
- **Threat:** MEDIUM, at conjecture/heuristic level. It already identifies the tristable regime and the front-speed criterion for its coexistence/clustering boundary, proves the coexistence side in d ≥ 3, and draws a conjectured 2D diagram with a 4A/4B split at the voter point. It says nothing about the order of the transition and nothing rigorous in d = 2.

#### Cox & Perkins (2005), "Rescaled Lotka–Volterra models converge to super-Brownian motion"
- **Citation:** Ann. Probab. 33(3), 904–947. — [DOI](https://doi.org/10.1214/009117904000000973)
- **Read:** ABS.
- **Content:** the low-density limit is super-Brownian motion with drift, covering long-range models and nearest-neighbour models "for dimensions three and above". — [DOI](https://doi.org/10.1214/009117904000000973)
- **Relevance:** M.

#### Cox & Perkins (2007), "Survival and coexistence in stochastic spatial Lotka–Volterra models"
- **Citation:** PTRF 139, 89–142. — [DOI](https://doi.org/10.1007/s00440-006-0040-3)
- **Read:** REV and SEC.
- **Content:** in d ≥ 3, survival near (1,1), and coexistence on a cone to the left of (1,1) with slope set by m₀ = q₂/(q₂ + q₃) < 1. — [CMP 2010](https://doi.org/10.1214/EJP.v15-795)
- **Relevance:** M.

#### Cox & Perkins (2008), "Renormalization of the two-dimensional Lotka–Volterra model"
- **Citation:** Ann. Appl. Probab. 18(2), 747–812; doi:10.1214/07-AAP453; arXiv:0804.0158.
- **Read:** FT §1 (verified by me from arXiv).
- **Theorem 1.3 (d = 2):** survival from a single 1 holds when 1 − r(η) < α₀, α₁ < 1 + r(η), and
  - α₁ − 1 < (1 − η)(α₀ − 1) if α₀ ≥ 1;
  - α₁ − 1 < (1 + η)(α₀ − 1) if α₀ < 1.
  - So survival includes part of the mean-field founder-control quadrant (1 < α₁ < α₀). — [arXiv](https://arxiv.org/abs/0804.0158)
- **Conjecture:** "There is a continuous curve α1 = g(α0), tangent to α1 = α0 at (1, 1), so that survival fails above the curve for α0 close to 1. The corresponding conjecture has been proved for d ≥ 3 to the left of (1, 1) with the slope m0 in place of 1, but in fact is false to the right of (1, 1) (with slope m0) in this higher dimensional setting." — [arXiv](https://arxiv.org/abs/0804.0158)
- **On d = 2:** "For d = 2 in Theorem 1.3, we are in the critical case m0 = 1, and the above proof of coexistence just fails … Theorem 1.3 and the above Conjecture suggest that for d = 2 their conjecture is quite delicate near α0 = 1." — [arXiv](https://arxiv.org/abs/0804.0158)
- **Relevance:** H. **Threat:** none.

#### Cox, Merle & Perkins (2010), "Coexistence in a two-dimensional Lotka–Volterra model"
- **Citation:** EJP 15, paper 38, 1190–1266; doi:10.1214/EJP.v15-795.
- **Read:** FT §1 (local).
- **Theorem 1.2 (d = 2):** coexistence on the "thorn" |α₀ − α₁| ≤ ((K − η)/γ)(1 − α₀)/(log(1/(1−α₀)))², for α₀ near 1. — [DOI](https://doi.org/10.1214/EJP.v15-795)
- **Abstract:** "The result is delicate since coexistence fails for the two-dimensional voter model which corresponds to the tip of the thorn." — [DOI](https://doi.org/10.1214/EJP.v15-795)
- **Proposition 1.3 (2D coalescence asymptotics):** the probability that n walkers from distinct points have not met by t is ∼ K_n/(log t)^{n(n−1)/2}. So q₂(t) ∼ γ/log t and q₃(t) ∼ K/(log t)³. — [DOI](https://doi.org/10.1214/EJP.v15-795)
- **On α > 1:** "When α0 and α1 are both greater than 1 … this setting should favour large clusters of particles of the same type, and at least for d ≤ 2, invariant measures should remain degenerate." — [DOI](https://doi.org/10.1214/EJP.v15-795)
- **Conjecture 1.6:** extinction of ones just above the thorn (sharpness). "The above strengthens Conjecture 2 of Neuhauser and Pacala". — [DOI](https://doi.org/10.1214/EJP.v15-795)
- **Relevance:** H. **Threat:** none.

#### Cox & Perkins (2014), "A complete convergence theorem for voter model perturbations"
- **Citation:** Ann. Appl. Probab. 24(1), 150–197; doi:10.1214/13-AAP919; arXiv:1210.0830.
- **Read:** FT §1 (local).
- **Theorem 1.1 (d ≥ 2):** complete convergence for the symmetric LV(α), α ∈ (α_c, 1).
- **Theorem 1.2 (d ≥ 3):** complete convergence with coexistence for cancellative voter-model perturbations with f′(0) > 0. — [arXiv](https://arxiv.org/abs/1210.0830)
- "the condition f′(0) > 0 means there is a positive drift for the local density of 1's when the density of 1's is very small … this condition promotes coexistence. It also excludes voter models themselves for which the complete convergence theorem fails." — [arXiv](https://arxiv.org/abs/1210.0830)
- On the voter model itself: "If d ≤ 2, then there are exactly two extremal stationary distributions, δ0 and δ1, and hence no coexistence" (Holley–Liggett 1975, Ann. Probab. 3, 643–663). — [arXiv](https://arxiv.org/abs/1210.0830)
- **Relevance:** M.

#### Cox & Perkins (2023/2025), "A complete convergence theorem for the q-voter model and other voter model perturbations in two dimensions"
- **Citation:** Electron. J. Probab. 30 (2025), Paper No. 9, 86 pp.; doi:10.1214/25-EJP1267; arXiv:2310.00812 (v2, 22 Nov 2023); Zbl 1557.60255. Paper number and DOI were confirmed by two sub-agents via Crossref/zbMATH. — [arXiv](https://arxiv.org/abs/2310.00812); [DOI](https://doi.org/10.1214/25-EJP1267)
- **Read:** FT §§1–5 (local; §1 re-read here).
- **Theorem 1.1:** for |N| ≤ 8 and d = 2, 3, 4, complete convergence with coexistence holds for the q-voter model with q_c < q < 1.
  - Remark 1.3: the hypotheses "also hold for d = 2 when N is the unit sphere in Z2 in the L∞ norm (|N| = 8)", i.e. the Moore neighbourhood. — [arXiv](https://arxiv.org/abs/2310.00812)
- **Theorem 1.9:** "Assume for 0 < ε ≤ ε0, ξ[ε] is a cancellative and monotone finite range voter model perturbation in Z2, and Θ3 > 0. There is an ε1 > 0 such that for ε ∈ (0, ε1), the complete convergence theorem with coexistence (CCT) holds for ξ[ε]." — [arXiv](https://arxiv.org/abs/2310.00812)
- **Corollary 1.17:** "Assume … ξ[ε] is a monotone asymptotically symmetric finite range voter model perturbation in Z2, and Θ2 + Θ3 > 0." Then for small ε a single 1 survives with positive probability. No cancellativity is needed. — [arXiv](https://arxiv.org/abs/2310.00812)
- **No reaction function in d = 2:** "The reaction function f is defined in terms of the invariant measures for the voter model and so is not well-defined for d = 2. Therefore … we replace f′(0) with the drift Θ3 in a super Brownian motion low density limit theorem." — [arXiv](https://arxiv.org/abs/2310.00812)
- **Remark 1.12 (q > 1):** "Θ3 < 0 will have the opposite sign … the conclusion of Theorem 1.11 with drift Θ3 < 0". This is only the low-density super-Brownian limit with negative drift. — [arXiv](https://arxiv.org/abs/2310.00812)
- **Founder control:** "When q > 1 … we now expect a type of founder control, where one type or another will take over, with probabilities depending on the initial configuration. The tools for proving such a result for q > 1 and close enough to 1 in two dimensions, or even for d ≥ 3, do not seem to be currently available." — [arXiv](https://arxiv.org/abs/2310.00812)
- **Conjectures 1.4 and 1.5:** complete convergence with coexistence for every neighbourhood, every d ≥ 2, and every 0 < q < 1. — [arXiv](https://arxiv.org/abs/2310.00812)
- **Theorems 5.6–5.8:** 2D Lotka–Volterra, affine voter and geometric voter models all have Θ₃ > 0 and satisfy complete convergence near the voter point. Strict subadditivity of r^s gives Θ₃ > 0 (Cor. 5.3). — [arXiv](https://arxiv.org/abs/2310.00812)
- **Relevance:** H. **Threat:** none. It explicitly leaves Θ₃ ≤ 0 open.

#### Agarwal, Simper & Durrett (2021), "The q-voter model on the torus"
- **Citation:** EJP 26, Paper No. 118, 33 pp.; doi:10.1214/21-EJP682; arXiv:2006.04162.
- **Read:** FT §1 and §5 (local).
- **Model:** 3D torus.
- "If q > 1, then the process rapidly reaches fixation on one opinion." — [arXiv](https://arxiv.org/abs/2006.04162)
- **Theorem 2:** du/dt = ±c_k u(1−u)(1−2u)f_k(u), where f_k > 0 on [0,1]. So the q-voter model has only the Case 1/Case 3 dichotomy and no tristable Case 4. — [arXiv](https://arxiv.org/abs/2006.04162)
- **Conjecture:** persistence for time exp(n^β) when q < 1. — [arXiv](https://arxiv.org/abs/2006.04162)
- **Relevance:** M. **Threat:** none.

#### Durrett & Remenik — the 2D PDE limit for voter perturbations
- **Source:** Durrett & Remenik, "Evolution of dispersal distance", J. Math. Biol. 64 (2012) 657–666. I read the 18 May 2011 preprint.
- **Statement of a 2D PDE limit:**
  - "In d = 2 the probability two random walks do not coalesce by time t is O(1/log t), so we need to further speed up time to τε = ε−2 log(1/ε)".
  - The limit is stated as a KPP-type equation with m* computed from the potential kernel.
  - "The mathematical details are still being worked out by Durrett and Remenik (2011)."
  - Their reference reads: "Voter model perturbations in two dimensions. Manuscript in preparation." — [Durrett–Remenik preprint](https://sites.math.duke.edu/~rtd/evoldisp/evol_dispersal2.pdf)
- **Publication status:** no published version of "Voter model perturbations in two dimensions" was found. Cox–Merle–Perkins 2010 cite it as "Preprint. MR2538084"; that MR number is UNVERIFIED. — [CMP 2010](https://doi.org/10.1214/EJP.v15-795)
- **Relevance:** M (only the asymmetric, two-lineage case; unpublished).

#### Durrett (2014), "Spatial evolutionary games with small selection coefficients"
- **Citation:** EJP 19, Paper No. 121; doi:10.1214/EJP.v19-3621; arXiv:1406.5876.
- **Read:** FT §§2, 6 (sub-agent).
- **Dimension:** d ≥ 3.
- **Theorem 6.1:** "In case S2 [bistable], 1's take over if ū < 1/2, and 2's take over if ū > 1/2." The symmetric ū = 1/2 case is excluded. — [arXiv](https://arxiv.org/abs/1406.5876)
- "bistablity in the replicator equation disappears in the spatial model". — [arXiv](https://arxiv.org/abs/1406.5876)
- **Open problem:** "we cannot handle the case in which the replicator equation is bistable … we are not aware of results that prove the existence of traveling waves solutions, or methods to prove convergence of solutions to them". — [arXiv](https://arxiv.org/abs/1406.5876)
- **Relevance:** M–H.

#### Cox & Durrett (2016), "Evolutionary games on the torus with weak selection"
- **Citation:** Stoch. Proc. Appl. 126(8), 2388–2409; arXiv:1511.04713.
- **Read:** FT §§5, 7 (sub-agent).
- **Content:** d ≥ 3 only; bistable games are decided by ū versus 1/2. There are no d = 2 statements and nothing on the order of transitions. — [arXiv](https://arxiv.org/abs/1511.04713)
- **Relevance:** M.

#### Huang & Durrett (2021), "Motion by mean curvature in interacting particle systems"
- **Citation:** PTRF 181, 489–532; doi:10.1007/s00440-021-01082-0; arXiv:2006.03090.
- **Read:** abstract verified by me; FT §§1.2–1.3 by sub-agent.
- **Abstract:** "in a wide variety of examples when the RDE limit has a bistable reaction term and traveling waves have speed 0, one can run time faster and further rescale space to obtain convergence to motion by mean curvature. This opens up the possibility of proving that the sexual reproduction model with fast stirring has a discontinuous phase transition, and that in Region 2 of the phase diagram for the nonlinear voter model studied by Molofsky et al there were two nontrivial stationary distributions." — [arXiv](https://arxiv.org/abs/2006.03090)
- **In the text (sub-agent):**
  - Region 3: "we expect clustering".
  - Case 4B: "we expect clustering".
  - Conjecture 1.5: clustering for the symmetric Lotka–Volterra model at θ₀ = θ₁.
  - Voter-perturbation results are for d ≥ 3. — [arXiv](https://arxiv.org/abs/2006.03090)
- **Relevance:** H. **Threat:** conjecture-level (clustering in 4B). "Discontinuous" is used only for the fast-stirring sexual-reproduction model.

#### Hughes & Lin, "Generalized front propagation for spatial stochastic population models"
- **Citation:** arXiv:2408.02817 (v2 Nov 2025).
- **Read:** ABS and §7 (sub-agent).
- **Content:** the "rescaled nonlinear voter model with d ≥ 3" converges to generalized mean-curvature flow. No phase-transition statements. — [arXiv](https://arxiv.org/abs/2408.02817)
- **Relevance:** M.

#### Durrett, *Interacting Particle Systems: Ideas, Techniques, Applications* (book draft)
- **Citation:** draft PDF PASTA0902.pdf, linked from the PASTA page; copyright notice 2025; not peer-reviewed.
- **Read:** §7.5 "Molofsky, et al (1999)", verified by me.
- **Statements:**
  - "if v > 0 which occurs in Region 4B then the two waves move toward each other … there is a traveling wave that decreases from 1 to 0, and has velocity 0 by symmetry and we have clustering." — [draft PDF](https://sites.math.duke.edu/~rtd/PASTA/PASTA0902.pdf)
  - "The dyad heuristic is a simple computation that allows one to determine the location of the line separating 4A and 4B … This is > 2 if 3p1 + p2 > 1 which is a line containing the points (0, 1) and (0.2, 0.4), the voter model." — [draft PDF](https://sites.math.duke.edu/~rtd/PASTA/PASTA0902.pdf)
  - "In order to use this machinery we have to reformulated the problem in Zd with d ≥ 3." — [draft PDF](https://sites.math.duke.edu/~rtd/PASTA/PASTA0902.pdf)
  - "We are not able to prove results in Regions 2 and 3 for the model symmetric under interchange of 0 and 1, since we have a bistable PDE with wave speed 0." — [draft PDF](https://sites.math.duke.edu/~rtd/PASTA/PASTA0902.pdf)
  - The front matter carries the line "Nonergodic two dimensional nonlinear voter model with positive rates?" — [draft PDF](https://sites.math.duke.edu/~rtd/PASTA/PASTA0902.pdf)
- **Relevance:** H.
- **Threat:** conjecture-level. A 2D coexistence/clustering line through the voter point is conjectured and placed by a dyad heuristic. Its order is not discussed.

#### Ma & Durrett (2018)
- **Citation:** Ann. Appl. Probab. 28(5), 2896–2921.
- **Read:** ABS.
- **Content:** "In d=2 … there are similar phase transitions, with coexistence occurring when (1+2θ)/(1+θ)<λ<(1+2θ)/θ". The 2D range is c√(log N); there is no bistability. — [DOI](https://doi.org/10.1214/17-AAP1378)
- **Relevance:** L–M.

#### Low-relevance rigorous works (no 2D lattice bistability or order-of-transition content)
- Yu-Ting Chen: Ann. Appl. Probab. 23 (2013); 28 (2018); SPA 157 (2023). These are on regular and finite graphs. — [arXiv:1107.1513](https://arxiv.org/abs/1107.1513)
- Huo & Durrett: latent voter model (SPA 2017) and zealot voter model (AAP 2019).
- Newman, Ravishankar & Schertzer, EJP 22 (2017): 1D voter perturbations.
- Bessonov & Durrett (2017): 2D quadratic contact process edge speeds.
- Köppl, Lanchier & Mercer, EJP 30 (2025): lattice games.

#### Citation sweep (sub-agent, Semantic Scholar and OpenAlex)
- Papers citing each source: CDP 2013, 40; Molofsky 1999, 93; Neuhauser–Pacala 1999, 173; Cox–Perkins 2023/2025, 1 (a Rev. Mod. Phys. 2026 review, for terminology only).
- Of these, only the works listed above treat 2D voter perturbations, the q-voter model or bistable cases.
- The sweep found no follow-up on the Θ₃ < 0 side in 2D. — [Semantic Scholar API](https://api.semanticscholar.org/graph/v1/paper/arXiv:1103.1676/citations)

### Inferences

**1. How the project's model fits Cox–Perkins.**
- The rule w(j) = j/8 + ε g(j) is, in continuous time, a nonlinear voter model in the Cox–Durrett/Cox–Perkins sense, with asymmetric part r^a = 0 and symmetric part r^s(A) = g(|A|). So Θ₂ = 0, and Θ₃ is the leading drift.
- **Contrast with asymmetric models.** In asymmetric 2D models such as Lotka–Volterra with α₀ ≠ α₁, the pair drift γ*(θ₀ − θ₁), of order 1/log, dominates. It decides survival even inside the mean-field founder-control quadrant (CP 2008, Thm 1.3).
- Z2 symmetry is therefore exactly the situation in which the sub-leading 3-lineage term, and the 4-lineage terms beyond it, decide the phase.

**2. Which theorems apply to the project's rule (continuous-time version only):**
- **Theorem 1.15** (super-Brownian limit with drift Θ₃): applies. It needs only asymptotic symmetry and finite range.
- **Corollary 1.17** (single-seed survival): applies where Θ₃ > 0. It needs monotonicity, which holds near p = 1 (p/8 ≥ (1−p)α).
- **Theorem 1.9** (complete convergence): needs cancellativity.
  - **My computation (not from the literature).** I applied Cox–Perkins' Prop. 3.2, from Cox–Durrett 1991, to |N| = 8:
    - Criterion: solve a_ℓ = Σ_k α_k M(k, ℓ), with M(k, j) = Σ_{odd i ≤ j∧k} C(j, i) C(8 − j, k − i), and require α_k ≥ 0.
    - M is invertible (det M = 2²⁸), so α is unique.
    - The test reproduces the known cases. The voter model gives α = (1/8, 0, …, 0). The threshold voter model gives all α_k = 1/128. The q-voter model gives α ≥ 0 for q = 0.9 and 0.99, and negative entries for q = 1.01 and 1.1.
    - For the project's rule, α_k for k ≥ 2 equals ε·(M^{-T}g)_k. This is negative at k = 2, 5 and 6 for α = 0.5 and 1 and every κ/α in [0, 0.5] tested. Example: α = 1, κ = 0.26 gives (M^{-T}g)₂, ₅, ₆ ≈ −0.047, −0.066, −0.047.
    - The full w(j) at p = 0.95 also fails.
    - Script: `/tmp/claude-1000/-home-matias-Desktop-silly-dn/7f0d9a3b-c595-4ffd-a22c-cdbe8fd1d247/scratchpad/canc/canc.py` (scratchpad; may not persist).
  - So Theorem 1.9 cannot be invoked through the known criterion. Whether some other cancellative representation exists is not settled; I think it unlikely, by symmetrizing over permutations of N, but have not proved it.
- **Discrete time.** All the theorems are continuous-time. The synchronous model has no rigorous counterpart, and the discrete-time 2D Molofsky model appears only in conjectured diagrams.

**3. Leading order versus sub-leading order in 2D.**
- In d ≥ 3, the tristable region 4 is a wedge with non-zero opening at the voter point. The front-speed line 5b₁ + b₂ = 0 is distinct from the line b₁ = 0 where f′(0) changes sign (CDP Fig. 5).
- In d = 2 at leading order there is a single threshold, Θ₃ = 0. This parallels CP 2008's "critical case m0 = 1", where the survival line and the diagonal coincide to leading order, and the coexistence thorn opened only through the next (log)⁻³ term (CMP 2010).
- By CMP's K_n/(log t)^{n(n−1)/2}, 4-class terms (∼ (log t)⁻⁶) are suppressed relative to 3-class terms (∼ (log t)⁻³) by about (log t)⁻³. With t ∼ ε⁻¹ polylog(1/ε), this suggests a window of relative width shrinking roughly like (log 1/ε)⁻³. That is consistent with "a first-order line weakening logarithmically toward the voter point" (heuristic; constants unknown).
- The project's own figure, that 4+-class terms still carry 28% of the weight at T = 65536, shows the asymptotic regime is very far away.
- A referee could argue that the 2D transition is "continuous at leading order". The defensible statement is "first-order at fixed small ε, with a window closing logarithmically".

**4. Prior art in d ≥ 3.** CDP Theorem 1.13(b) proves density ≈ 1/2 for every non-trivial stationary law on the 4A side. Huang–Durrett and Durrett's 2026 draft expect clustering in 4B. Together they imply a discontinuity across the 4A/4B (front-speed) line in d ≥ 3. That reading is my inference; none of these authors draws it.

**5. What would be new.** The project would be the first to:
- identify that line explicitly as a first-order transition;
- show that it survives in d = 2 as a sub-leading window around Θ₃ = 0;
- give its location for D&N/Moore (κ*/α ≈ 0.26).
No published Θ₃ value exists for any Moore-neighbourhood rule other than the sign for strictly subadditive rates (CP Cor. 5.3). D&N's R is far from subadditive.

### Gaps
- There is no rigorous 2D reaction–diffusion (hydrodynamic) limit for symmetric voter perturbations. Durrett–Remenik's 2D work appears never to have been published.
- Final EJP versions were not compared with the arXiv texts.
- The CMP 2010 citation "MR2538084" for Durrett–Remenik is unverified.
- Whether the D&N-based rule is cancellative by some route other than Prop. 3.2 is open.

## 3. Neuhauser–Pacala (1999) and founder control: d = 2 versus d ≥ 3; is there a rigorous front-speed criterion or a proved discontinuity in d = 2?

### Takeaway
**Neuhauser–Pacala's conjectures:**
- Conj. 1: coexistence on the diagonal below 1.
- Conj. 2: in the mean-field founder-control quadrant (α₀, α₁ > 1), the spatial model picks a winner off the diagonal, rather than the outcome depending on initial conditions.

**Status near the voter point:**
- **d ≥ 3:** Conj. 2 is proved via the sign of ∫f, i.e. the wave speed (CDP Thm 1.9).
- **d = 2:** partially. CP 2008 proves survival of the favoured type's single seed inside the founder-control quadrant.

**The symmetric founder-control diagonal (α₀ = α₁ > 1):**
- This is the analogue of CDP Case 3 and of the project's clustering side.
- It has no theorem in any d. CMP 2010 expect degenerate invariant measures "at least for d ≤ 2"; CP 2025 say the tools are missing.

There is no rigorous 2D front-speed criterion for a symmetric bistable case, and no proved discontinuity.

Lotka–Volterra has no stable mixed state in its founder-control region, so its bistability is "0 versus 1" (Case 3), not "absorbing versus active" (Case 4) as in the project.

### Cited Findings

#### Neuhauser & Pacala (1999), "An explicitly spatial version of the Lotka–Volterra model with interspecific competition"
- **Citation:** Ann. Appl. Probab. 9(4), 1226–1259; doi:10.1214/aoap/1029962871; MR1728561.
- **Read:** ABS only (full text blocked). The conjectures are known through quotation by CMP 2010, CP 2008 and CDP 2013.
- **Abstract:** "Local competitive interactions between species result in (1) a reduction of the parameter region where coexistence occurs in the classical model, (2) a reduction of the parameter region where founder control occurs in the classical model, and (3) spatial segregation of the two species in parts of the parameter region where the classical model predicts coexistence." — [UMN Experts abstract](https://experts.umn.edu/en/publications/an-explicitly-spatial-version-of-the-lotka-volterra-model-with-in/)
- **Conjecture 1** (quoted by CMP): "For d ≥ 2 coexistence holds for LV(α0, α1) whenever α0 = α1 < 1." Neuhauser–Pacala proved it for small α. — [CMP 2010](https://doi.org/10.1214/EJP.v15-795)
- **Conjecture 2** (paraphrased by CDP): "Conjecture 2 of [37] states that 1's take over for α ∈ Λ01, α0 > 1 and 0's take over for α ∈ Λ00, α0 > 1." The exact wording and any dimension restriction are UNVERIFIED. — [CDP](https://arxiv.org/abs/1103.1676)
- **Relevance:** H. **Threat:** none.

#### CDP (2013), Lotka–Volterra section
- **Heuristic:** "On R4 ∪ R5, u* ∈ (0, 1) is an unstable fixed point, while 0 and 1 are attracting fixed points for the ODE. This is case 2 of Durrett and Levin [19], so we expect the winner of the competition to be predicted by the direction of movement of the speed of the decreasing traveling wave … Symmetry dictates that the speed is 0 when θ0 = θ1 … r has the same sign as ∫₀¹ f(u) du". — [arXiv](https://arxiv.org/abs/1103.1676)
- **Theorem 1.9** (d ≥ 3, finite-range p, α near (1,1)): 0's take over on Λ₀^η, 1's take over on Λ₁^η. "Theorem 1.9 establishes this result asymptotically as α gets close to (1, 1) at least for d ≥ 3." — [arXiv](https://arxiv.org/abs/1103.1676)
- This is the only rigorous front-speed criterion in a bistable voter perturbation (d ≥ 3, asymmetric).

#### Cox & Perkins (2008), 2D
- Theorem 1.3 and the d = 2 conjecture are quoted in §2 above.
- The phrase "false to the right of (1, 1) (with slope m0) in this higher dimensional setting" refers to the d ≥ 3 result that the survival boundary in the founder-control quadrant has front-speed slope 1, not the low-density slope m₀. — [arXiv](https://arxiv.org/abs/0804.0158)

#### Cox, Merle & Perkins (2010)
- For α₀, α₁ > 1: "at least for d ≤ 2, invariant measures should remain degenerate" (intuition, not a theorem). — [DOI](https://doi.org/10.1214/EJP.v15-795)

#### Cox & Perkins (2014, 2025)
- Complete convergence for the symmetric LV(α) with α < 1 near 1 in d ≥ 2 (2014, Thm 1.1; 2025, Thm 5.6). — [arXiv 2014](https://arxiv.org/abs/1210.0830); [arXiv 2025](https://arxiv.org/abs/2310.00812)
- Nothing for α > 1.

#### Huang & Durrett (2021)
- Conjecture 1.5: clustering for the symmetric Lotka–Volterra model on the line θ₀ = θ₁ (d ≥ 3 setting; sub-agent FT). — [arXiv](https://arxiv.org/abs/2006.03090)

#### Durrett's "Case 2" principle (heuristic)
- **Durrett 2009:** "Case 2. Two locally attracting fixed points. In this case, the limiting behavior of the ODE depends on the initial density. However, this is not the expected behavior for the particle system, and in this case the outcome of competition is dictated by the behavior of the PDE." — [Durrett 2009, arXiv:0906.2293](https://arxiv.org/abs/0906.2293)
- **Durrett 1999:** "in the stochastic spatial model, there is one stronger equilibrium that is the winner starting from any positive initial density. To determine the stronger equilibrium, we can start with one half-plane in each equilibrium and watch the direction of movement of the front." — [SIAM Rev. 41 (1999)](https://doi.org/10.1137/S0036144599354707)
- The principle originates in Durrett & Levin, Theor. Popul. Biol. 46 (1994) 363–394 (ABS only). — [DOI](https://doi.org/10.1006/tpbi.1994.1032)

### Inferences
- **The two criteria coincide in 2D at leading order.** In d ≥ 3 the Lotka–Volterra survival boundary to the right of (1,1) is the front-speed line (slope 1), distinct from the low-density slope m₀ < 1 (CP 2008 remark; CDP Thm 1.9). In d = 2, m₀ = 1, so the two coincide at leading order. This is the Lotka–Volterra analogue of the project's statement that κ_a and κ_front merge at Θ₃ = 0 as T → ∞. It is a good structural parallel to cite.
- **Founder control is ruled out off the symmetric line.** The Durrett–Levin principle, rigorous in the cases above, says that in infinite volume bistability turns into "the front decides". The project's claim fits this: the active phase wins where its front advances, and the absorbing phase wins otherwise. The new element is that one of the bistable states is an active, non-trivial phase, and that the switch happens at fixed small ε in d = 2.
- **No 2D front-speed theorem exists for any bistable voter perturbation.** The 2D results are low-density drift results only.

### Gaps
- The Neuhauser–Pacala full text was not read (Project Euclid blocked). The exact Conjecture 2 wording and its dimension restriction are unverified.
- No 2D result on the symmetric founder-control line exists in any source found.

## 4. Fast stirring and critical attractive systems: bistable reaction terms, the wave-speed sign, and whether "critical seeds die while the full configuration survives" is compatible with first order

### Takeaway
**Fast stirring.** Bistable reaction terms are handled rigorously in the limit:
- The critical value converges to the zero of the wave speed (∫f = 0), not to the ODE saddle-node. For the sexual-reproduction process β_c → 4.5, not 4.
- Combined with the limiting equilibrium density (ρ₂ → 2/3), this gives a density jump in the ε → 0 limit only. Nothing is proved at fixed ε in d ≥ 2.
- The first proofs of a discontinuity at fixed (large) stirring appeared in 2026, and only in d = 1, for one-trap pair and triplet creation processes (Durrett, arXiv:2608.29878).

**Durrett's 2002 Memoir.** Its coexistence theorem needs a Lyapunov function that blows up at the boundary, so it excludes bistable cases.

**Bezuidenhout–Gray (1994).**
- For finite-range attractive one-trap systems, the set of parameters where a finite seed survives is open, so finite seeds die at criticality.
- The upper invariant density is right-continuous, so if it jumps it is positive at λ_c. "Critical seed dies, full configuration survives" is therefore exactly the first-order signature, and is compatible with this theory.
- Their hypothesis, per the review, excludes two-trap voter-type models.

### Cited Findings

#### Durrett & Neuhauser (1994), "Particle systems and reaction-diffusion equations"
- **Citation:** Ann. Probab. 22(1), 289–333; doi:10.1214/aop/1176988861.
- **Read:** REV and SEC. Durrett's reprint is an image-only scan.
- **Review:** "The authors show for five examples that there exist non-trivial limits, which coincide with the critical value at which the speed of the traveling wave changes its sign." — [DOI](https://doi.org/10.1214/aop/1176988861)
- **Durrett 2009, Theorem 7:** "If we introduce fast stirring at rate ε−2 then in the sexual reproduction model, βc → 4.5 as ε → 0". The wave-speed sign is that of ∫₀^{ρ₂} g; at β = 4.5 the roots are 0, 1/3, 2/3. — [Durrett 2009](https://arxiv.org/abs/0906.2293)
- **Durrett 1999 (SIAM):** δ_e(ε) → 2/9, and "the critical value for survival from a finite set δf(ε) → 2/9". — [SIAM Rev.](https://doi.org/10.1137/S0036144599354707)
- **Dimension:** the bistable propagation result is "first proved in d = 1 by Fife and McLeod … See page 296 and the appendix of [20] for the extension to d > 1" (CDP). — [CDP](https://arxiv.org/abs/1103.1676)
- **Secondhand, UNVERIFIED in Durrett–Neuhauser itself:** Chatterjee–Durrett 2013 report that "Based on simulations they conjectured that the phase transition was continuous." — [Chatterjee–Durrett](https://arxiv.org/abs/1011.1567)
- **Relevance:** H. **Threat:** none (one trap).

#### Noble (1992), "Equilibrium behavior of the sexual reproduction process with rapid diffusion"
- **Citation:** Ann. Probab. 20(2), 724–745; doi:10.1214/aop/1176989802.
- **Read:** REV.
- **Review:** "an upper bound for the critical value of λ is given by limsup_{ε↓0} λ_c(ε) ≤ 4.5 and it is also conjectured that liminf_{ε↓0} λ_c(ε) ≥ 4.5". — [DOI](https://doi.org/10.1214/aop/1176989802)
- Nobody calls it "first order"; the jump holds only in the ε → 0 limit.
- **Relevance:** H. **Threat:** none.

#### Durrett (2002), "Mutual invadability implies coexistence in spatial models"
- **Citation:** Mem. AMS 156, no. 740; doi:10.1090/memo/0740.
- **Read:** REV and SEC.
- **Main theorem (as stated in Durrett 2009, Thm 5):** "Suppose there is a function φ that (i) decreases along solutions of the mean-field ODE, (ii) is convex, and (iii) φ(u) → ∞ when mini ui → 0. Then there is coexistence in the model with fast stirring". Condition (iii) excludes bistable cases. — [Durrett 2009](https://arxiv.org/abs/0906.2293)
- **Relevance:** M–L.

#### Durrett (2009), "Coexistence in stochastic spatial models"
- **Citation:** Ann. Appl. Probab. 19(2), 477–496; doi:10.1214/08-AAP590; arXiv:0906.2293.
- **Read:** FT (by me).
- **Problem 4:** uniqueness of the non-trivial stationary distribution for sexual reproduction with β > 4.5 and fast stirring.
- **Problem 5, the catalyst / Ziff–Gulari–Barshad model on Z².** Both all-CO and all-O are absorbing; the model is not Z2-symmetric. "Simulations suggest that the density of O atoms in equilibrium drops to 0 discontinuously at the upper critical value, but proving this is a much harder problem." — [arXiv](https://arxiv.org/abs/0906.2293)
- **Relevance:** H (a 2D two-trap model with a conjectured discontinuity).

#### Durrett (1999), "Stochastic spatial models"
- **Citation:** SIAM Rev. 41(4), 677–718.
- **Read:** FT (sub-agent).
- Conjecture 3.1, for the quadratic contact process with stirring: "If the stirring rate ν > 0 then δf = δe". — [DOI](https://doi.org/10.1137/S0036144599354707)
- Physics simulations reported by Chatterjee–Durrett (Evans–Guo–Liu) find δ_f < δ_e, i.e. a discontinuous transition, at low hopping. This is in tension with Conjecture 3.1. — [Chatterjee–Durrett](https://arxiv.org/abs/1011.1567)
- A sub-agent flagged a possible sign slip in Durrett's "It is known in a fair amount of generality that δf ≥ δe" (UNVERIFIED).

#### Bezuidenhout & Gray (1994), "Critical attractive spin systems"
- **Citation:** Ann. Probab. 22(3), 1160–1194; doi:10.1214/aop/1176988599.
- **Read:** REV only.
- **Main theorem (review):** "For every r < ∞ in discrete and continuous time the set {(β, δ) : δ(N'_r) > 0, ξ^{0} is viable} is open in the space of parameters (β, δ) corresponding to r-range attractive systems". — [DOI](https://doi.org/10.1214/aop/1176988599)
- **Consequences:**
  - A finite seed dies at the critical point.
  - Nothing is said about the upper invariant measure, which is compatible with a jump.
  - The condition δ(N'_r) > 0 (a site whose neighbours are all occupied can still die) excludes two-trap voter-type models. This is the sub-agent's reading of the review.
- **Relevance:** M.

#### Bramson & Gray (1991), "A useful renormalization argument"
- **Citation:** in the Spitzer Festschrift (Progress in Probability 28), 113–152; zbMATH 0745.60103.
- **Read:** REV and SEC.
- **Content:** a simpler proof of Toom's theorem, and a necessary-and-sufficient criterion for survival of generic attractive models. — [zbMATH](https://zbmath.org/?q=an:0745.60103)
- **Durrett, *Ten lectures*:** "Corollary 4.7. Two dimensional threshold birth rates for the Moore neighborhood in two dimension are robust if θ ≤ 3". — [Ten lectures (scan)](https://sites.math.duke.edu/~rtd/reprints/paper85.pdf)
- **Relevance:** M.

#### Durrett (1986), "Some peculiar properties of a particle system with sexual reproduction"
- **Citation:** in *Stochastic Spatial Processes*, LNM 1212, 106–111; doi:10.1007/BFb0076241.
- **Read:** FT, via an OCR'd scan by a sub-agent; symbols reconstructed.
- **Model:** Z², one trap; Toom's NE birth rule.
- **Results:** λ_f = ∞ (finite seeds never survive), while λ_c ≤ 110 (Durrett–Gray contour argument).
- **Conjecture:** "ρ(λ) ≡ P(0 ∈ ξ^1_∞) is discontinuous at λc … 'general nonsense' (see Griffeath (1981)) implies ρ(λ) is right continuous and hence ρ(λc) > 0." — [scan](https://sites.math.duke.edu/~rtd/reprints/paper30.pdf); [DOI](https://doi.org/10.1007/BFb0076241)
- **Relevance:** H.
- **Threat:** none. This is a one-trap growth model, but it is an explicit 2D discontinuity conjecture.

#### Chatterjee & Durrett (2013), "A first order phase transition in the threshold θ ≥ 2 contact process on random r-regular graphs and r-trees"
- **Citation:** Stoch. Proc. Appl. 123(2), 561–578; arXiv:1011.1567.
- **Read:** ABS and introduction (sub-agent).
- **Result:** a proved first-order transition, but only on trees and random regular graphs.
- "However outside the physics literature, see da Silva and de Oliveira (2011), there are no results about the nature of the phase transition on Z^d." — [arXiv](https://arxiv.org/abs/1011.1567)
- Open question Q2: "Is ρ(λ) discontinuous at λc? If so, then soft results imply that P(ξ^1_∞(x) = 1) > 0 when λ = λc." — [arXiv](https://arxiv.org/abs/1011.1567)
- **Related:** Fontes & Schonmann, PTRF 141 (2008) 513–541, treat threshold contact processes on trees. — [arXiv:math/0603109](https://arxiv.org/abs/math/0603109)

#### Durrett (2026), "Discontinuous phase transitions in the one-dimensional pair and triplet creation processes"
- **Citation:** arXiv:2608.29878 (v1 30 Aug 2026, v2 9 Sep 2026); no journal reference.
- **Read:** ABS verified by me; FT §§2.2–2.4 by sub-agent.
- **Content:** proves the transitions are discontinuous in both models on Z when the stirring rate is large. One trap; d = 1. — [arXiv](https://arxiv.org/abs/2608.29878)
- **Relevance:** M (a methodological template, with a forbidden-region argument after Chatterjee–Durrett).

#### Durrett, book draft (2026), on sexual reproduction and autocatalytic models
- "It seems difficult to prove that these models have a discontinuous phase transition, especially since the conclusion may be false." — [draft PDF](https://sites.math.duke.edu/~rtd/PASTA/PASTA0902.pdf)

#### Huang & Durrett (2021), on the 2D fast-stirring sexual-reproduction model
- Per sub-agent FT: "Theorem 1.2 suggests that there is a discontinuous phase transition but does not rule out the possibility that the phase transition could be continuous for any ǫ > 0." — [arXiv](https://arxiv.org/abs/2006.03090)

#### Front-motion limits (no traps)
- Katsoulakis & Souganidis, Arch. Rational Mech. Anal. 127 (1994) 133–157: generalized front evolution (metadata only). — [DOI](https://doi.org/10.1007/BF00377658)
- Funaki & Tsunoda, J. Stat. Phys. 177 (2019) 183–208: motion by mean curvature from Glauber–Kawasaki dynamics (ABS). — [arXiv](https://arxiv.org/abs/1812.10182)
- **Relevance:** L–M.

### Inferences
- **Answer to the Bezuidenhout–Gray question: yes, the two are compatible.**
  - For attractive one-trap systems, survival from a finite seed is an open condition, so seeds die at λ_c.
  - The upper-invariant density is right-continuous, so a jump means ρ(λ_c) > 0.
  - Hence "critical seed dies, full configuration survives" is precisely what a first-order transition looks like.
  - The 2D Toom-NE sexual-reproduction process is a rigorous extreme case: λ_f = ∞ while λ_c ≤ 110.
- **For two-trap symmetric systems, Bezuidenhout–Gray does not apply.**
  - There is no maximal non-trap initial state.
  - The natural analogue is the limit from Bernoulli(1/2), the "odd upper invariant law" of a cancellative system (Sturm–Swart, Lemma 1).
  - The relevant open question is Sturm–Swart's Q1 (coexistence ⇔ survival).
  - CDP's Case 4A remark shows that survival from a single seed and coexistence can both hold while the low-density drift is negative. So a first-order transition need not break Q1 except possibly at the critical point.
- **"The phase is set by the wave-speed sign"** is rigorous only:
  - in the fast-stirring limit (one-trap examples);
  - in d ≥ 3 voter perturbations with non-zero wave speed (CDP).
  Nothing in d = 2 makes it rigorous at fixed parameters.

### Gaps
- The primary texts of Durrett–Neuhauser, Noble, Bezuidenhout–Gray and the 2002 Memoir were not read (paywall or blocked). Their statements are from reviews and Durrett's later accounts.
- The claim that Durrett–Neuhauser conjectured continuity is secondhand.
- Bezuidenhout–Gray's own remarks on the upper invariant measure are unknown.

## 5. Must a 2D clustering/coexistence transition with two traps be continuous? Proved or conjectured discontinuities, and who conjectured what about the bistable case in d = 2

### Takeaway
**Nothing forces continuity.** No theorem anywhere says that a clustering/coexistence (or survival) transition in a 2D two-type or one-trap spin system must be continuous; all three literature sweeps came up empty.

**Nothing proves a discontinuity on Z² either.** Discontinuity is proved only on trees and random regular graphs, in fast-stirring limits, and (2026) on Z with large stirring.

**Conjectured 2D discontinuities** exist only for:
- one-trap growth models (Durrett 1986);
- the non-symmetric two-trap ZGB catalyst (Durrett 2009, Problem 5).

**For the symmetric bistable case in d = 2:**
- The probabilists' consensus conjecture is clustering on the founder-control or bistable side, made by CDP, Cox–Merle–Perkins, Cox–Perkins 2025, Huang–Durrett and Durrett's 2026 draft.
- A coexistence/clustering line through the 2D voter point is drawn as conjectured (CDP Fig. 3; Durrett 2026 draft Fig. 7.2).
- Nobody conjectures that this line is first order.

The closest explicit first-order claim is a physics one: Rodrigues et al. 2015 report a discontinuous P–A line ending at the voter point on 2D/3D lattices. It is contested by Azizi & Pleimling 2020.

### Cited Findings

#### No continuity theorem; difficulty of proving the order of a transition
- Chatterjee–Durrett 2013: "there are no results about the nature of the phase transition on Z^d" outside physics. — [arXiv](https://arxiv.org/abs/1011.1567)
- Swart: "It can be difficult to prove that a given phase transition is first or second order." — [arXiv](https://arxiv.org/abs/1703.10007)

#### Explicit 2D discontinuity conjectures (none Z2-symmetric with two traps)
- Durrett 1986: Toom-NE sexual reproduction, one trap. — [DOI](https://doi.org/10.1007/BFb0076241)
- Durrett 2009, Problem 5: ZGB catalyst, two non-symmetric traps. — [arXiv](https://arxiv.org/abs/0906.2293)
- Huang–Durrett 2021 aim at "proving that the sexual reproduction model with fast stirring has a discontinuous phase transition". — [arXiv](https://arxiv.org/abs/2006.03090)
- Durrett's 2026 draft cautions that such a discontinuity "may be false". — [draft PDF](https://sites.math.duke.edu/~rtd/PASTA/PASTA0902.pdf)

#### Conjectures on the symmetric bistable case
- **CDP 2013 (d ≥ 3 framework):** "One would guess that clustering occurs" in Case 3. 4B gets "no result". Fig. 3 is a conjectured 2D diagram with 4A coexistence / 4B / 3 clustering. — [arXiv](https://arxiv.org/abs/1103.1676)
- **Cox, Merle & Perkins 2010:** α₀, α₁ > 1 "should favour large clusters … at least for d ≤ 2, invariant measures should remain degenerate." — [DOI](https://doi.org/10.1214/EJP.v15-795)
- **Cox & Perkins 2025:** q > 1 should give "a type of founder control"; the tools are "not … currently available" in 2D or even d ≥ 3. — [arXiv](https://arxiv.org/abs/2310.00812)
- **Agarwal, Simper & Durrett 2021:** proved rapid fixation for q > 1 on the 3D torus. "Mean-field calculations suggest that there should be coexistence between opinions if q < 1 and clustering if q > 1." — [arXiv](https://arxiv.org/abs/2006.04162)
- **Huang & Durrett 2021:** clustering expected in Region 3 and Case 4B. — [arXiv](https://arxiv.org/abs/2006.03090)
- **Durrett, 2026 draft:** "we have clustering" in 4B. The dyad heuristic puts the 2D 4A/4B line at 3p₁ + p₂ = 1 through the voter point. — [draft PDF](https://sites.math.duke.edu/~rtd/PASTA/PASTA0902.pdf)
- **Sturm & Swart 2008, 2015; Swart & Vrbenský 2010:** for minority-favouring models, coexistence in d ≥ 2 for all parameters, so no non-trivial transition; "dc ≈ 4/3". — [arXiv:math/0701555](https://arxiv.org/abs/math/0701555); [arXiv:1311.0417](https://arxiv.org/abs/1311.0417); [arXiv:0911.1266](https://arxiv.org/abs/0911.1266)
- **Molofsky et al. 1999:** in 2D "there is no multiple outcome regime of the spatial model". — [DOI](https://doi.org/10.1006/tpbi.1998.1404)

#### Rigorous 2D two-trap examples on the clustering/fixation side (context)
- **Zero-temperature Glauber (majority) dynamics on Z²** (SEC via Morris and later papers):
  - Nanda–Newman–Stein: every site flips infinitely often from p = 1/2 in 2D.
  - Fontes–Schonmann–Sidoravicius: fixation to + for p near 1.
  - Morris, PTRF 149 (2011): p_c(Z^d) → 1/2 as d → ∞. — [arXiv:0809.0353](https://arxiv.org/abs/0809.0353)
- **U-voter dynamics** (both uniform states absorbing): Blanquicett (2020) proves p_c^vot(Z^d, U) < 1 for a wide class of U, i.e. fixation at + for initial densities near 1. — [arXiv:2003.02420](https://arxiv.org/abs/2003.02420)
- **Sowers 1999:** a 2D threshold voter interface converges to mean-curvature flow. — [DOI](https://doi.org/10.1006/jfan.1999.3477)

#### Physics (brief; read at the levels stated)
- **Rodrigues, Chatelain, Tomé & de Oliveira**, J. Stat. Mech. (2015) P01035; FT by a sub-agent.
  - "As to the PA line, that separates the P and A phases, it is a discontinuous transition line that ends at the voter point." — [arXiv:1501.01116](https://arxiv.org/abs/1501.01116)
  - The model is in the same class as the project's: 6-neighbour symmetric flips with two absorbing states, on 2D triangular and 3D cubic lattices.
- **Azizi & Pleimling**, PRE 102, 022112 (2020), report that in Rodrigues et al.'s model "a line of discontinuous phase transitions is replacing the generalized voter critical line" ("in stark contrast to the Langevin description"). Their own 2D model instead shows Ising and DP lines merging into a generalized-voter line away from the voter point. — [arXiv:2007.10933](https://arxiv.org/abs/2007.10933)
- **Castellano, Muñoz & Pastor-Satorras 2009:** the mean-field intermediate regime is "washed out" by fluctuations in spatially extended systems. — [arXiv:0907.1775](https://arxiv.org/abs/0907.1775)
- **Llabrés, San Miguel & Toral 2026:** mean-field abrupt switch with a metastability window. — [arXiv:2505.11358](https://arxiv.org/abs/2505.11358)

### Inferences
**Novelty verdict, rigorous side.** A rigorous result on a discontinuous 2D clustering/coexistence transition, or even on the Θ₃ < 0 side beyond the low-density limit, would be new. The field explicitly lacks the tools (CP 2025; Durrett 2014 open problem; Durrett 2026 draft; Chatterjee–Durrett 2013). The project's contribution is heuristic or numerical, and should be framed that way.

**Novelty verdict, conjecture and heuristic side.**
- The structure is anticipated: a tristable reaction term in which the active (1/2) phase and the traps are all locally stable, with the coexistence/clustering boundary at a zero front speed and the line starting at the voter point.
- Sources: CDP 2013 (d ≥ 3, with a conjectured 2D diagram), Huang–Durrett 2021 and Durrett's 2026 draft.
- Nobody states that this boundary is first order, analyses it in d = 2, or computes its location for any 2D rule.
- The project should present its claim as: "the CDP 4A/4B front-speed line is a first-order transition (interface density jumps; fronts reverse; critical droplets), and in d = 2 it survives near the voter point as a log-suppressed window around Cox–Perkins' Θ₃ = 0".
- It must also cite Rodrigues et al. 2015 (closest numerical prior art) and Azizi & Pleimling 2020 (counter-evidence).

**Main vulnerability.**
- At leading order in 2D (CP 2025) the low-density and front criteria coincide, so the transition is "continuous at leading order" in the sense of a single threshold.
- Molofsky's 2D observation that small seeds nucleate growing patches shows that a mean-field multiple-outcome regime can vanish in 2D.
- The project should therefore show:
  - the window survives at accessible ε;
  - the interface density of the active phase stays bounded away from 0 up to the front-speed line;
  - seeds below R_c have survival probability decaying with system-size-independent rates.

**Why the q-voter cannot be a counterexample.** Its reaction term is a single-signed cubic times a positive polynomial (Agarwal–Simper–Durrett, Thm 2), so it can never have a Case-4 window. The D&N-based rule can, because g(j) has mixed signs, with R strongly non-subadditive.

### Gaps
- The rigorous literature has no 2D symmetric nonlinear-voter study of the transition type, so the project has no rigorous comparator.
- Not found or not read: any statement in Liggett's 1999 book about the transition type for intermediate thresholds; Durrett's 1988 lecture notes; Durrett–Gray 1985 (an unpublished manuscript per Durrett's publication list).
- Physics sources are covered only briefly here; the project's LITERATURE.md and other researchers cover them.
