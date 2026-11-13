# Cellular-automaton literature for a noisy, asynchronous Day & Night model: self-complementary Life-like rules, α-asynchronous and probabilistic 2D CA, two-absorbing-state PCA, and frozen non-uniform configurations

Scope note for the report writer. Every entry gives the citation, what I actually read (full text, abstract only, or a secondary source), the model and update scheme, the claim (quoted where I saw the text), a relevance rating (High/Medium/Low), and a novelty-threat flag for the project's claim. The claim is a *fully* first-order transition in a 2D model with two Z2-symmetric absorbing states, from noisy α-asynchronous Day & Night (D&N). Items already in `campaign/theory/LITERATURE.md` (Molofsky 1999, Cox–Durrett–Perkins, Cox–Perkins 2023, Al Hammal 2005, Castellano 2009, Rodrigues 2015, Durrett–Steif 1993, Liggett 1994) are not repeated, only cross-referenced.

Notation used below: j = number of the 8 Moore neighbours that disagree with the cell. D&N flips iff j ∈ {3,6,7,8}. Anneal flips iff j ∈ {4,6,7,8}. Majority/"Vote" flips iff j ≥ 5. Minority (Moore, self included) flips iff j ≤ 4. I checked these translations against the rule strings and the CAM-book lookup tables (see Q1).

Citation corrections found while verifying the brief:
- Olejarz–Krapivsky–Redner is **PRL 109, 195702 (2012)**, not PRL 107 (2011).
- Spirin–Krapivsky–Redner "Freezing in Ising ferromagnets" is listed by arXiv as **PRE 65, 016119 (2002)**.
- The twisted-majority (ANNEAL) rule is credited by Toffoli & Margolus to **Vichniac (1986)**, not to Vichniac, Physica D 10 (1984). The 1984 paper covers voting rules and synchronous Ising CAs (details in Q1).
- Bagnoli–Boccara–Rechtman PRE 63, 046116 (2001) is **one-dimensional**.

---

## Q1. Is there prior work on Day & Night and its relatives (Anneal/twisted majority, other self-complementary outer-totalistic rules), on noisy versions of them, or on Z2-symmetric probabilistic outer-totalistic rules on the Moore neighbourhood?

### Takeaway
I found **no academic work on a stochastic, noisy or α-asynchronous Day & Night**, nor any phase-transition study of D&N; the closest prior art is in three places.
- **Descriptive deterministic studies.** Eppstein (2010) describes D&N and Anneal on random initial conditions: coarsening at density ½, shrinking minority clusters otherwise.
- **The noisy twisted majority, RAND-ANNEAL.** Toffoli & Margolus (1987) built it: a Z2-symmetric probabilistic voting rule that adds noise only at the marginal counts next to the majority threshold, "using noise rather than frustration" to thaw frozen majority boundaries. It was not analysed as a phase transition.
- **Vichniac (1984).** He catalogued deterministic voting rules on von Neumann and Moore neighbourhoods. To cure checkerboard artefacts in synchronous Ising CAs, he replaced the rule by a coin flip inside checkerboard domains.

These are the direct ancestors of "D&N + symmetric noise". None of them threatens the novelty claim.

### Cited Findings

**Day & Night: background, hobbyist and tertiary sources**
- D&N is B3678/S34678. It "was invented and named in 1997 by Nathan Thompson" and studied in detail by David I. Bell. The page calls it Wolfram "class 4", and says random initial conditions show "the formation and gradual coarsening of large day and night regions" with persistent activity at boundaries. Read: full page. Its only academic reference is Eppstein (2010). Relevance: Low (background). Threat: none. — [Wikipedia, Day and Night (cellular automaton)](https://en.wikipedia.org/wiki/Day_and_Night_(cellular_automaton))
- **Conflict on classification:** Eppstein says one can "view the behavior of Day & Night as belonging to Wolfram's Class II", against Wikipedia's Class 4. Eppstein also dates Bell's unpublished article to 1994, against 1997 elsewhere. — [Eppstein, arXiv:0911.2890](https://arxiv.org/abs/0911.2890); [Wikipedia](https://en.wikipedia.org/wiki/Day_and_Night_(cellular_automaton))
- LifeWiki returned HTTP 403 to my fetches; LifeWiki pages (OCA:Day & Night, Self-complementary) were not read. A search snippet attributes a Rule 110 unit cell in D&N (Turing completeness) to Peter Naszvadi (2020); this is hobbyist and unverified. — [LifeWiki Self-complementary (search result only)](https://conwaylife.com/wiki/Self-complementary)
- A NetLogo community model "Day and Night" exists (hobbyist). Not read. — [NetLogo community model](https://ccl.northwestern.edu/netlogo/models/community/Day%20and%20Night)

**Eppstein, D. (2010)**, "Growth and Decay in Life-Like Cellular Automata", in A. Adamatzky (ed.), *Game of Life Cellular Automata*, Springer, pp. 71–98, doi:10.1007/978-1-84996-217-9_6, arXiv:0911.2890.
- Read: full text of the D&N and Anneal paragraphs (§7). Model: deterministic synchronous outer-totalistic rules.
- On D&N: it "has the curious property that, if one exchanges the roles of live and dead cells, the rule's behavior is unchanged".
- From random fields with p < ½, "the live cells form clusters that gradually shrink, eventually forming a collection of scattered oscillators". At p = ½, "the clusters and live and dead cells gradually grow and merge", with cluster size and boundary curvature monotone in time.
- On Anneal (B4678/S35678): "Similarly to Day & Night, any pattern in this rule has the same behavior if all live cells are replaced by dead cells and vice versa". It forms "growing clusters of live and dead cells, with scattered oscillators". "Wolfram's classification is incapable of distinguishing this rule from Day & Night".
- B01245/S0125 is a "strobing version of Day & Night".
- Relevance: Medium. It is the only academic description of D&N's coarsening and its oscillator debris. Threat: none. — [arXiv:0911.2890](https://arxiv.org/abs/0911.2890)

**Toffoli, T. & Margolus, N. (1987)**, *Cellular Automata Machines: A New Environment for Modeling*, MIT Press (Series in Scientific Computation).
- Read: full text of §5.4 "Voting rules" and §8.3 "Voting rules, revisited" (searchable PDF on N. Margolus's MIT page).
- **MAJORITY.** 9SUM {0 0 0 0 0 1 1 1 1 1}, i.e. flips iff j ≥ 5. It gives "stable, quite fragmented domains" and "a frozen pattern of local alliances".
- **ANNEAL.** 9SUM {0 0 0 0 1 0 1 1 1 1}, i.e. flips iff j ∈ {4,6,7,8} (= B4678/S35678). It is credited to "Gerard Vichniac [67]", which is Vichniac 1986 (below). "By swapping the two table entries that are adjacent to the threshold ... one encourages reshuffling at the boundary between black and white domains, where the majority is marginal". The result is "a good model of surface tension, where boundaries behave as stretched membranes that exert a pull proportional to their curvature". Caveat in the book: "the ANNEAL rule had missed smoothing out a few rare spots: it's hard for a deterministic rule to guarantee 'smoothness' on the scale of only a few cells".
- **RAND-ANNEAL (§8.3).** A 5-cell (von Neumann plus self) majority, 6SUM {0 0 RAND −RAND 1 1}, with RAND = 1 with probability p (p = 1/32 in the book). "Here we proceed in an analogous way, this time using noise rather than frustration as a means of maintaining fluidity on the boundaries". "With a nondeterministic rule such as RAND-ANNEAL (with p not 0 or 1) such stuck patterns are eventually melted away". Fig. 8.2: "The frozen boundaries (a) produced by 6MAJ are thawed in (b) by the thermal noise provided by a probabilistic variant of this rule, namely RAND-ANNEAL". At p → 1 "the rule becomes a 5-neighbor version of ANNEAL".
- Update scheme: synchronous (CAM hardware).
- Relevance: **High**. It is the earliest Z2-symmetric probabilistic "noisy twisted majority" PCA that keeps both uniform states absorbing, and an explicit precedent for curing frozen majority states by noise placed only at marginal counts. Threat: none, since no phase-transition analysis is given. — [CAM book PDF](https://people.csail.mit.edu/nhm/cam-book.pdf)

**Vichniac, G. Y. (1984)**, "Simulating physics with cellular automata", Physica D 10, 96–116.
- Read: the OCR'd article inside the Computer History Museum scan of the Physica D 10 proceedings, §3.2 and §4.4–4.6.
- **§3.2, threshold voting rules.** Rules that output 1 if the count is above a threshold "always lead to the growth of clusters of one of the 'species' (0 or 1) until a stable configuration is reached". "The growth consists in filling concavities, and halts once the convex shapes are reached".
- Simple-majority rules are "self-conjugate" with critical initial concentration ½ ("percolation"). Biased thresholds give "nucleation".
- Moore entries: 4/8 threshold p_c = 0.333 and 3/8 threshold p_c = 0+.
- He also notes synchronous voting on a checkerboard background, whose "states oscillate between '0' and '1' both in space and time".
- **§4.6, fix for synchronous Ising CAs.** "Define now a new rule that follows the cellular-automaton Ising rule for aligned domains, but ignores it in the checkerboard domains of high energy. The new rule will flip a coin instead. Physically, it puts these regions at an infinite temperature. We observe then a growth of the T = 0 domains at the expense of the dislocations and the oscillatory monsters."
- **Twisted majority not found.** I did not find the twisted-majority/ANNEAL rule in this paper; searches of the OCR for "twist", "anneal" and "majority" found only the voting-rule and simulated-annealing passages. OCR gaps are possible.
- Relevance: **High**, as precedent: frozen convex clusters under deterministic majority, and a coin flip at checkerboard (j = 4-type) configurations to dissolve them. Threat: none. — [Physica D 10 scan (Computer History Museum)](https://tcm.computerhistory.org/CHMfiles/Cellular%20Automation%201984.pdf)

**Vichniac, G. (1986)**, "Cellular automata models of disorder and organization", in *Disordered Systems and Biological Organization* (Bienenstock et al., eds.), Springer, pp. 1–20.
- Read: citation only, via the CAM book's reference [67]. This is the source Toffoli & Margolus give for ANNEAL; its content is unverified. — [CAM book PDF, refs.](https://people.csail.mit.edu/nhm/cam-book.pdf)

**Bagnoli, F. & Mencarelli, L. (2026)**, "Coarsening and Bifurcations in Wide-Range Two-Dimensional Totalistic Cellular Automata", arXiv:2604.10082 (11 Apr 2026). No journal reference listed.
- Read: full text. Model: deterministic totalistic majority on a circular neighbourhood of radius R (K odd), serial or parallel update. Their "frustrated majority" maps uniform states to their opposite (F(0) = 1, F(K) = 0). It is not Vichniac's twist.
- Key quote: "In any majority model, there are at least two absorbing states ... there are also other stable configurations, formed by clusters of sites whose boundaries present a curvature related to the interacting radius ... These configurations are not forecast by a mean-field analysis, and are unstable under stochastic perturbation, even when these obeys the same symmetry law of the original rule, like in the voter model."
- Parallel update: "one can have oscillating checkerboard patterns for certain values of R".
- Initial-density bifurcation: "One can safely assume that the transition is first-order for large values of R".
- Relevance: Medium. It is the most recent (2026) statement that frozen non-uniform majority states are destroyed by symmetric stochastic (voter-like) perturbations. Threat: none (deterministic; bifurcation in initial density). — [arXiv:2604.10082](https://arxiv.org/abs/2604.10082)

**Kornyak, V. V. (2006)**, "Cellular Automata with Symmetric Local Rules", CASC 2006, LNCS 4194, pp. 240–250, arXiv:math-ph/0605040.
- Read: abstract. It covers rules with local permutation invariance (the "generalized Game of Life family") and the compatibility of several Life-family rules. The abstract does not name D&N. Relevance: Low. — [arXiv:math-ph/0605040](https://arxiv.org/abs/math-ph/0605040)

**Rendell (Journal of Cellular Automata)**, "A Life-Like Rule Where Random Islands Either Grow or Shrink". Already in LITERATURE.md. Citation page found; not re-read. — [Old City Publishing listing](https://www.oldcitypublishing.com/?p=22025)

### Inferences
- **What the predecessors share.** The D&N-family literature splits into deterministic descriptions (Eppstein; hobbyist pattern work) and the "voting-rule" line (Vichniac 1984 → ANNEAL 1986 → RAND-ANNEAL 1987). The second line already contains the project's design move: keep both uniform states absorbing, keep Z2 symmetry, and add randomness at the marginal counts to stop boundaries freezing. It never asked about phase transitions or the order of the active/absorbing transition. The project can cite it as the origin of "noisy self-complementary voting rules" without novelty risk.
- **RAND-ANNEAL in the project's notation (my own translation of the lookup table, z = 4).** Flip probabilities are w(0) = w(1) = 0, w(2) = p, w(3) = 1 − p, w(4) = 1. This is a von Neumann analogue of the project's family. The book's "stuck patterns are eventually melted away" is about small spots: straight interfaces still have j ≤ 1 on both sides and so stay frozen for 0 < p < 1. Toffoli–Margolus's noise therefore removes local frozen debris but not spanning flat stripes, the same situation as zero-temperature Ising (Q4).
- **Eppstein's deterministic D&N on random fields is a useful α = 1, noise-free baseline.** At density ≠ ½ the minority shrinks to "scattered oscillators", which are not absorbing under synchronous updating. At density ½ it coarsens. The project's model departs from this through α < 1 (no deterministic oscillators) and through noise.

### Gaps
- No stochastic, noisy or asynchronous D&N paper was found in arXiv, APS, Springer or the asynchronous-CA surveys (two independent searches, plus checking the reference lists of Fatès 2014, Regnault et al. 2010 and Vieira et al. 2021). Absence is not proof; Google Scholar full-text search was not available.
- Vichniac (1986) not read, so its exact ANNEAL definition and any noise discussion are unverified.
- Peña & Sayama, "Life Worth Mentioning: Complexity in Life-Like Cellular Automata" (reportedly Artificial Life, 2021) was **not found or verified**. It may contain D&N in a rule-space survey.
- I found no systematic study of Z2-symmetric probabilistic outer-totalistic rules on the Moore neighbourhood with general w(j), apart from the rigorous voter-perturbation framework already in LITERATURE.md (Cox–Perkins).

---

## Q2. What phase transitions are known in probabilistic and α-asynchronous 2D Life-like CA, and have any discontinuous transitions been reported?

### Takeaway
- **Asynchronous Life.** α-asynchronous Game of Life has a continuous transition at σ_c ≈ 0.906 (Blok & Bergersen 1999, verified). Its **static exponents are DP-like but its dynamic exponents are nonuniversal**, which the authors attribute to the "many possible absorbing states" (frozen still lifes). The absorbing ("frozen") phase is on the *synchronous* side.
- **Probabilistic Life.** Monetti & Albano (1995) reported **discontinuous (first-order) irreversible transitions** between an extinct (frozen) phase and a living phase.
- **Rule space.** Reia & Kinouchi (2014/2015) report a discontinuous absorbing transition across the space of outer-totalistic rules.
- **Asynchronous minority.** 2D minority shows an empirical α_c ≈ 0.57 (Moore) / 0.83 (von Neumann), below which frozen stripes or checkerboards take over (Regnault–Schabanel–Thierry).

All of these have a single absorbing phase or many frozen ones; none is Z2-symmetric with exactly two absorbing states. So they bound but do not threaten the novelty claim. Monetti–Albano is the most important precedent for "first-order absorbing transition in a 2D stochastic Life-like CA".

### Cited Findings

**Blok, H. J. & Bergersen, B. (1999)**, "Synchronous versus asynchronous updating in the 'game of Life'", Phys. Rev. E 59, 3876–3879, doi:10.1103/PhysRevE.59.3876.
- Read: **full text**, verified. The brief's "(verify: asynchronous Life, DP-like?)" is confirmed with a caveat on the dynamic exponents.
- Model: Conway's Life (B3/S23), square lattice, periodic, L = 32–512. Each site updates independently with probability σ per step (α-asynchronous); time is rescaled by σ. Random 50% initial state.
- Abstract: "a critical phase transition is observed ... the static exponents appear to belong to the directed percolation universality class in 2+1 dimensions. The dynamic exponents, however, are nonuniversal, as seen in other systems with multiple absorbing states."
- Phases: "for large σ the system eventually decays to a low-density frozen state. The one exception is σ = 1 which also allows simple periodic structures to survive". Small σ gives a fluctuating state, and "As σ → 0+ the steady-state is characterized by domains of alternating dead and live stripes."
- Order parameters: density and "activity (fraction of updated sites which flip states)".
- Numbers:
  - σ_c(∞) = 0.9060 ± 0.0004.
  - β = 0.595 ± 0.004 (density) and β′ = 0.595 ± 0.006 (activity), against DP β = 0.586 ± 0.014.
  - Spreading: δ = 0.25 ± 0.04 and ν_∥ = 0.93 ± 0.02, against DP δ = 0.451 and ν_∥ = 1.295. They comment: "The key ingredient for this to occur appears to be the possibility of the system settling into one of many possible absorbing states."
- Relevance: **High**. It is the closest analogue of the project's p = 0 line: an α-asynchronous Life-like rule with infinitely many frozen absorbing configurations. Threat: none (continuous, single active phase). — [Blok & Bergersen PDF](https://www.cs.ubc.ca/~rikblok/research/published/blok99.pdf)

**Monetti, R. A. & Albano, E. V. (1995)**, "Critical edge between frozen extinction and chaotic life", Phys. Rev. E 52, 5825–5831, doi:10.1103/PhysRevE.52.5825.
- Read: **abstract only** (APS page); the stochastic rules themselves were not visible.
- Claims: "by introducing stochastic rules in the GL it is found that ''life'' exhibits a very rich critical behavior". "Discontinuous (first-order) irreversible phase transitions (IPT's)" separate an extinct phase from a living steady state. "The edge of the frozen-chaotic transition coincides with that of the IPT's life extinction."
- Vieira et al. (2021) summarise the line: "Both continuous [13] and discontinuous [12] phase transitions can be observed, depending on the choice of probabilities", where [12] = Monetti & Albano and [13] = Nordfalk & Alstrøm.
- Relevance: **High** as prior art. Threat: **Low–Medium**. It is a first-order absorbing transition in a 2D stochastic Life-like CA, but with a single extinct absorbing phase, not Z2. It must be cited and distinguished. — [APS abstract](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.52.5825); [Vieira et al. arXiv:2101.07225](https://arxiv.org/abs/2101.07225)

**Nordfalk, J. & Alstrøm, P. (1996)**, "Phase transitions near the 'game of Life'", Phys. Rev. E 54, R1025–R1028, doi:10.1103/PhysRevE.54.R1025.
- Read: abstract only. Probabilities are introduced into Life's rules ("a hypercubic space of stochastic CA"); they "identify a phase-transition line separating 'Life' and 'Death'", which "ends very close to but not in the 'game of Life'". The abstract does not state the universality class; Vieira et al. classify it as continuous.
- Relevance: Medium. Threat: none. — [APS abstract](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.54.R1025)

**Schulman, L. S. & Seiden, P. E. (1978)**, "Statistical mechanics of a dynamical system based on Conway's game of Life", J. Stat. Phys. 19, 293–314, doi:10.1007/BF01011727.
- Read: Springer metadata and abstract summary via search only. Model: Life-like rules plus a stochastic element set by a "temperature". It reports a phase transition and a "thermodynamic" constant of motion. The order of the transition is not given in what I saw.
- Several IBM pages list the wrong journal; Springer and arXiv citations agree on J. Stat. Phys.
- Relevance: Low–Medium (historical first stochastic Life). Threat: none. — [Springer](https://link.springer.com/doi/10.1007/BF01011727)

**Vieira, A. P., Goles, E. & Herrmann, H. J. (2021)**, "Phase transitions in a conservative Game of Life", Phys. Rev. E 103, 012132, arXiv:2101.07225.
- Read: abstract and introduction. Model: Moore-neighbourhood Life with conservative, asynchronous, nonlocal swaps of "unsatisfied" live and dead cells.
- Claim: "a discontinuous phase transition between an inactive phase, in which the dynamics freezes after a finite time, and an active phase"; the second transition is continuous on the square lattice.
- Relevance: Medium (a discontinuous freezing transition in a 2D Life-like rule). Threat: none (conservative, not Z2). — [arXiv:2101.07225](https://arxiv.org/abs/2101.07225)

**Reia, S. M. & Kinouchi, O. (2014)**, "Conway's game of life is a near-critical metastable state in the multiverse of cellular automata", Phys. Rev. E 89, 052123, arXiv:1407.1006.
- Read: abstract. Mean-field maps for 6144 outer-totalistic rules. "The transition region in rule space seems to correspond to a nonequilibrium discontinuous absorbing phase transition". "LIFE is not at the 'border of chaos,' but thrives on the 'border of extinction.'"
- Relevance: Medium. Threat: Low (transition across rule space, single absorbing vacuum). — [arXiv:1407.1006](https://arxiv.org/abs/1407.1006)

**Reia, S. M. & Kinouchi, O. (2015)**, "Nonsynchronous updating in the multiverse of cellular automata", Phys. Rev. E 91, 042110 (DOI listed on arXiv), arXiv:1503.01350.
- Read: abstract. Compares synchronous, asynchronous and sequential updates over 6144 rules. "The first-order phase-transition in the multiverse of synchronous cellular automata ... seems to be robust not only to changes in update schema but also to different initial densities."
- Relevance: Medium. Threat: Low. — [arXiv:1503.01350](https://arxiv.org/abs/1503.01350)

**Roy, S., Paul, S. & Adak, S. (2024)**, "A Note on α-Asynchronous Life-Like Cellular Automata", ASCAT 2024 proceedings (Springer), pp. 192–203, doi:10.1007/978-3-031-56943-2_15.
- Read: abstract only. Studies Life-like rules under α-asynchrony, including "phase transition dynamics" and "their connection with second-order phase transition and first-order irreversible phase transition" for random rule sets. The abstract names no specific rules.
- Relevance: Medium. Threat: **unknown, possibly Low–Medium**. It is the only α-asynchronous Life-like rule survey found after 2010 that mentions first-order transitions. The full text should be checked for D&N or self-complementary rules. — [DTU Orbit record](https://orbit.dtu.dk/en/publications/a-note-on%CE%B1-asynchronous-life-like-cellular-automata/)

**Fatès, N. (2014)**, "A guided tour of asynchronous cellular automata", J. Cell. Autom. 9(5–6), 387–416, arXiv:1406.0792. Extends the Automata'13 invited paper, LNCS 8155, pp. 15–30.
- Read: **full text** (arXiv v2). Volume and pages confirmed from the author's publication list.
- On Blok & Bergersen: a transition "from a 'static' behaviour, where the system would settle on fixed points, to a 'living' behaviour ... labyrinth-like patterns", called "a second-order phase transition".
- On ECA under α-asynchrony (Fatès 2009): "no less than nine rules also displayed phase transitions ... directed percolation ... A unique rule, namely ECA 178, was shown to belong to another universality class, a fact that is explained by the symmetric role that is played by 0s and 1s in the transition rule."
- On 2D minority: it "will settle out on a checkerboard or on a stripe-like pattern depending on whether the rule is defined with the von Neumann or the Moore neighbourhood. To our knowledge, there is no mathematical explanation of this empirical observation."
- On full versus small asynchrony: "consider the minority rule in 2D with a von Neumann neighbourhood: with a fully asynchronous updating, the two uniform fixed points are not reachable from a non-uniform configuration [40] but this is not longer true if we allow a small degree of synchrony" ([40] = Fatès & Gerin, J. Cell. Autom. 4(4), 323–337, 2009).
  - Interpretation caveat (mine): all-0 and all-1 are *not* fixed points of minority, so "uniform fixed points" here presumably means the two checkerboard fixed points. Fatès & Gerin was not read.
- On Bersini & Detours (Artificial Life IV, 1994): their "stabilising effect" of asynchrony on Life, "freeze" onto labyrinths, was "intimately linked to finite-size effects".
- Relevance: **High** as a map of the field. Threat: none. It mentions no D&N, Vichniac, anneal or twisted majority. — [arXiv:1406.0792](https://arxiv.org/abs/1406.0792); [Fatès publications page](https://members.loria.fr/nazim.fates/publications.html)

**Fatès, N. (2009)**, "Asynchronism induces second order phase transitions in elementary cellular automata", J. Cell. Autom. 4(1), 21–38.
- Read: citation and claims via the 2014 survey, plus the arXiv abstract snippet (nlin/0703044): transitions "belonging to the directed percolation or to the parity conservation universality classes".
- My own decoding of the lookup table: ECA 178 flips the centre iff at least one neighbour disagrees. It is therefore the **1D Z2-symmetric, two-absorbing-state α-asynchronous** analogue of the project's model, and it gives a continuous non-DP (PC) transition.
- Relevance: Medium. Threat: none (1D). — [Fatès 2014 survey](https://arxiv.org/abs/1406.0792); [arXiv:nlin/0703044](https://arxiv.org/abs/nlin/0703044)

**Fatès, N. & Morvan, M. (2004)**, "Perturbing the Topology of the Game of Life Increases its Robustness to Asynchrony", arXiv:nlin/0405061. Proceedings venue not verified.
- Read: abstract and introduction. Studies "a phase transition from an 'inactive-sparse phase' to a 'labyrinth phase'" against synchrony rate and missing-link rate.
- Fatès (2010), "Does Life resist asynchrony?", in Adamatzky (ed.), *Game of Life Cellular Automata*, Springer, pp. 257–274, doi:10.1007/978-1-84996-217-9_14, is cited via the survey; not read.
- Relevance: Medium. Threat: none. — [arXiv:nlin/0405061](https://arxiv.org/abs/nlin/0405061); [Springer chapter page](https://www.springerprofessional.de/en/does-life-resist-asynchrony/1964024)

**Regnault, D., Schabanel, N. & Thierry, É. (2009)**, "Progresses in the analysis of stochastic 2D cellular automata: a study of asynchronous 2D Minority", Theor. Comput. Sci. 410(47–49), 4844–4855, doi:10.1016/j.tcs.2009.06.024, arXiv:0706.2479 (MFCS 2007 abstract version).
- Read: **full text of §1–2 and Fig. 1**.
- Model: 2D minority, von Neumann or Moore, α-asynchronous (each cell updates with probability α); α = 0 means one random cell per step.
- Empirical result (Fig. 1): for α below a threshold, "A stable pattern emerges (checkerboards and stripes for N- and M-neighborhoods respectively) and rapidly covers the whole configuration when both dimensions are even". α_c ≈ 0.83 (von Neumann) and **α_c ≈ 0.57 (Moore)**.
- Synchronous regime: it "converges to cycles of length 1 or 2"; from random initial conditions it gives "large flashing" regions. Asynchrony "creates a noise that progressively erodes the flashing homogenous large regions".
- Relevance: **High**. It is the 2D α-asynchronous totalistic CA whose ordered phase consists exactly of frozen stripes and checkerboards, with an α-driven transition. Threat: none. — [arXiv:0706.2479](https://arxiv.org/abs/0706.2479)

**Regnault, D., Schabanel, N. & Thierry, É. (2010)**, "On the analysis of 'simple' 2D stochastic cellular automata", DMTCS 12(2), 263–294, doi:10.46298/dmtcs.518.
- Read: **full text of §1–3**. Model: Moore-neighbourhood minority, fully asynchronous ("one single random cell updates at each time step").
- Definitions: "A cell is said active if its state would change if fired"; "A configuration c is stable if it remains unchanged under the dynamics, i.e., if all its cells are inactive."
- The Fig. 2 caption, comparing Life, "some random outer-totalistic" rules and Minority on the Moore neighbourhood: "Similar stripes emerge in asynchronous regime even if their synchronous behavior differ drastically." Also: "the stripes arise as well in many other asynchronous automata especially in totalistic cellular automata".
- They prove almost-sure convergence to stable stripe configurations (polynomial time under conditions on the initial state). Earlier work showed "the synchronous dynamics eventually leads to stable configurations or cycles between two conjugated configurations" (Goles & Martínez 1990).
- Relevance: **High** for frozen stripes and checkerboards in asynchronous Moore-neighbourhood CA. Threat: none. — [DMTCS 518](https://dmtcs.episciences.org/518)

**Rouquier, J.-B., Regnault, D. & Thierry, É. (2011)**, "Stochastic minority on graphs", Theor. Comput. Sci. 412(30), 3947–3963 (citation via the Fatès survey), arXiv:1011.5119.
- Read: abstract. Fully asynchronous minority on arbitrary graphs; convergence time "strongly dependent on the topology". Relevance: Low. — [arXiv:1011.5119](https://arxiv.org/abs/1011.5119)

**Bidaux, R., Boccara, N. & Chaté, H. (1989)**, "Order of the transition versus space dimension in a family of cellular automata", Phys. Rev. A 39, 3094–3105.
- Read: APS abstract. Model: probabilistic "two-state automata with one absorbing phase", rule chosen by mean-field theory, d = 1–4.
- Claims: "first-order transitions are indeed exhibited as soon as d>1"; "For d=1, the transition is continuous, but with critical exponents different from those of directed percolation."
- Relevance: **High** as context: first-order absorbing transitions in 2D PCA are established since 1989. Threat: Low (single absorbing state). This is likely the "Chaté" work most relevant to the brief. — [APS abstract](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.39.3094)

### Inferences
- **Two structural analogues for the p = 0 line.** In Blok–Bergersen the frozen (absorbing) phase sits at *high* synchrony and asynchrony creates the active labyrinth/stripe phase. In Regnault et al. it is the reverse: asynchrony (α < α_c) creates frozen stripes or checkerboards, and synchrony keeps the system flashing. The project's α-asynchronous D&N at p = 0 combines both mechanisms: frozen j ∈ {0,1,4,5} patterns, and α controlling whether deterministic synchronous oscillations survive. So the α dependence at p = 0 need not be monotone. That is worth checking numerically.
- **First-order transitions are not new in themselves.** Discontinuous transitions in stochastic or asynchronous Life-like CA exist: Monetti–Albano 1995, Vieira et al. 2021, and in rule space Reia–Kinouchi. The project's distinguishing elements must be Z2 symmetry with exactly two absorbing states, plus the *fully* first-order character (interface density jumps, fronts, droplets). "First-order in a 2D Life-like CA" alone would not be novel.
- **The 1D anchor.** Fatès's ECA 178 result is the natural anchor for "α-asynchronous self-complementary rule ⇒ symmetric absorbing-state class (PC in 1D)". The 2D analogue would generically be GV unless something makes the transition discontinuous.

### Gaps
- Full texts of Monetti & Albano 1995 and Nordfalk & Alstrøm 1996 not read: the exact stochastic rules, the noise placement, and how "frozen" was defined.
- Roy–Paul–Adak 2024 full text not read: whether D&N or Anneal appear, and what "first-order irreversible" means there.
- Fatès (2010), "Does Life resist asynchrony?", not read; it may survey other Life-like rules.
- No α-asynchronous study of any self-complementary 2D Life-like rule (D&N, Anneal, Vote) was found.

---

## Q3. What is known about 2D Z2-symmetric probabilistic CA with two absorbing states: Domany–Kinzel endpoints, compact DP, Bagnoli–Boccara–Rechtman-type models, Lübeck, Grassberger, Chaté & Manneville?

### Takeaway
- **Bagnoli–Boccara–Rechtman is 1D** (radius 1, three-input totalistic PCA). Its first-order line runs *between the two absorbing phases*, and the later Bagnoli work is also 1D (long-range) or deterministic 2D. No 2D two-absorbing-state PCA by the Bagnoli group was found.
- **Compact DP**, the DK endpoint where both uniform states are absorbing, is exactly solvable in 1D. Its density jumps (β = 0 "in all dimensions", Lübeck), and the two major reviews disagree on whether to call that "discontinuous": Hinrichsen says it is discontinuous; Lübeck says "misleading" and that it is really second order.
- **Hinrichsen (1997) 2D model.** It has two symmetric absorbing states and gives roughly mean-field exponents. Hinrichsen conjectures **1 < d_c ≤ 2** for two absorbing states.
- **Lipowski & Droz (2002) nonequilibrium Potts model, d = 2, q = 2.** It shows **a jump in the active density** but inconclusive dynamics. They place d = q = 2 near the point where continuous non-MF, MF and **first-order** regions meet. This is the closest prior hint of a discontinuous transition in a 2D Z2 two-absorbing-state model and is a **novelty concern to address explicitly**.
- **Grassberger (2005).** He shows that a 2+1D generalized DK model with a *single* absorbing state turns first order when hole formation is suppressed.

### Cited Findings

**Bagnoli, F., Boccara, N. & Rechtman, R. (2001)**, "Nature of phase transitions in a probabilistic cellular automaton with two absorbing states", Phys. Rev. E 63, 046116, doi:10.1103/PhysRevE.63.046116, arXiv:cond-mat/0002361.
- Read: **full text of §1–2** plus abstract.
- Model: "a one-dimensional, one-component, totalistic PCA with two absorbing states". The neighbourhood is the site and its 2 nearest neighbours, with p_s = P(x′ = 1 | σ = s). "The transition probabilities from an empty neighborhood is zero, and that from a completely occupied neighborhood is one." Synchronous.
- Abstract (v2 PDF): "two second-order and one first-order transition lines that meet at a bicritical point". Note: the arXiv listing page renders this as "tricritical point". This is a source inconsistency; the PDF and journal say bicritical.
- "Our model exhibits a first-order phase transition along the line of symmetry in the upper part of the phase diagram", shown via hysteresis with a small symmetry-breaking field. DP lines. At the bicritical point (p_t = 0.460(2)) the kinks are in the PC class (η = 0.00(2), δ = 0.285(5), z = 1.18(2)).
- Definition: "An absorbing state is represented by a set of configurations from which the system cannot escape". Construction: "if we assume that particles can be absorbed or desorbed only near the boundaries of a patch ... (when the neighborhood is not homogeneous), then the completely empty and occupied configurations are absorbing states."
- Relevance: High (1D precedent, already in LITERATURE.md). Threat: Low (1D; the first-order line is absorbing–absorbing, not active–absorbing). — [arXiv:cond-mat/0002361](https://arxiv.org/abs/cond-mat/0002361)

**Bagnoli, F., Franci, F. & Rechtman, R. (2005)**, "Phase transitions of extended-range probabilistic cellular automata with two absorbing states", Phys. Rev. E 71, 046108, arXiv:cond-mat/0405604.
- Read: abstract. A "long-range one-dimensional cellular automaton with two symmetric absorbing states", with ferromagnetic linear coupling competing with an antiferromagnetic nonlinear coupling. Relevance: Low–Medium. Threat: none. — [arXiv:cond-mat/0405604](https://arxiv.org/abs/cond-mat/0405604)

**Bagnoli, F., Franci, F. & Rechtman, R. (2002)**, "Opinion formation and phase transitions in a probabilistic cellular automaton with two absorbing states", LNCS 2493, pp. 249–258, arXiv:nlin/0511001.
- Read: skimmed the full text. One-dimensional, r = 1 nearest neighbours, with a first-order transition between the absorbing states (as in the 2001 paper). Relevance: Low. — [arXiv:nlin/0511001](https://arxiv.org/abs/nlin/0511001)

**Bagnoli, F. & Rechtman, R. (2014)**, "Phase Transitions of Cellular Automata", arXiv:1409.4284 (book chapter). Read: abstract and reference grep only. A review of asynchronism-induced, synchronization and chaotic transitions, citing the 2001 two-absorbing-state model. Relevance: Low. — [arXiv:1409.4284](https://arxiv.org/abs/1409.4284)

**Hinrichsen, H. (2000) review**, "Nonequilibrium critical phenomena and phase transitions into absorbing states", Adv. Phys. 49, 815, arXiv:cond-mat/0001070.
- Read: **full text of the CDP section**. At DK p₂ = 1: "there are two symmetric absorbing states ... CDP has a global Z2 symmetry".
- "In 1+1 dimensions CDP is fully equivalent to a zero temperature Ising model with Glauber dynamics or the voter model". Exponents: β = 0, β′ = 1, ν_∥ = 2, ν_⊥ = 1, δ = ½.
- "Because of the vanishing exponent β, the CDP transition is discontinuous."
- Relevance: High (definitions). Threat: none. — [arXiv:cond-mat/0001070](https://arxiv.org/abs/cond-mat/0001070)

**Lübeck, S. (2004/5) review**, "Universal scaling behavior of non-equilibrium phase transitions", arXiv:cond-mat/0501259. The journal reference (Int. J. Mod. Phys. B) was **not verified**.
- Read: **full text of §3.1 (DK) and Table 2**.
- CDP: "The associated critical exponent β = 0 differs in all dimensions from the directed percolation values". It "is characterized by a continuously vanishing percolation order parameter P_perc and algebraically diverging correlations lengths, indicating a second order phase transition. But due to the misleading discontinuous behavior of the steady state order parameter ρ_a the phase transition was often considered as first order."
- 2+1D: Lübeck simulates a generalized DK on a bcc lattice (site DP, L ≤ 2048) only for DP scaling functions. I found **no 2+1D study of the Z2-symmetric endpoint** in this review.
- Relevance: High (definitional conflict on "discontinuous" at symmetric endpoints). Threat: none. — [arXiv:cond-mat/0501259](https://arxiv.org/abs/cond-mat/0501259)
- **Conflict:** Hinrichsen calls the CDP transition "discontinuous"; Lübeck calls that "misleading" and says it is second order. — [Hinrichsen](https://arxiv.org/abs/cond-mat/0001070); [Lübeck](https://arxiv.org/abs/cond-mat/0501259)

**Hinrichsen, H. (1997)**, "Stochastic lattice models with several absorbing states", Phys. Rev. E 55, 219, arXiv:cond-mat/9608065.
- Read: **full text**.
- Models: generalized DK (model I) and generalized contact process (model II), with n "colours" of inactive state plus an active state, so each site has 3 states for n = 2.
- Abstract: "for n = 2 both models belong to the same universality class as branching annihilating walks with an even number of offspring" (PC in 1D).
- Design rule (C): "Boundaries between inactive domains of different colors are not allowed to stick to each other irreversibly." Otherwise "the colors would then be irrelevant ... and thus belong to directed percolation."
- **2D (model II, 80 × 80):** "active sites in two-dimensional models with two absorbing states are arranged in fractal 'lines' along the boundaries of inactive islands". λ_c = 0.99(1), δ = 0.9(1), η = 0.00(5), z = 1.0(1), "roughly ... mean-field". "We conjecture that 1 < d_c ≤ 2 is the upper critical dimension of systems with two absorbing states."
- On CDP: "On this line active clusters are not fractal but compact".
- Relevance: **High**. Threat: **Medium**. It is a prior 2D two-symmetric-absorbing-state model, continuous and roughly mean-field, consistent with "2D is marginal" (Al Hammal 2005). Distinguish it: three-state sites; no fronts or droplets; small lattice. — [arXiv:cond-mat/9608065](https://arxiv.org/abs/cond-mat/9608065)

**Lipowski, A. & Droz, M. (2002)**, "Phase transitions in nonequilibrium d-dimensional models with q absorbing states", Phys. Rev. E 65, 056114, arXiv:cond-mat/0112130.
- Read: **full text**.
- Model: Metropolis dynamics of the q-state Potts model with one change: "when all neighbours of a given site are in the same state as this site, then this site cannot change its state". This gives exactly q absorbing states for T > 0. Square lattice, random sequential.
- Abstract: "In two dimensions and q = 3 the model exhibits a discontinuous transition. For the three-dimensional case and q = 2 the model exhibits a continuous, transition with β = 1 (mean-field). Simulations are inconclusive, however, in the two-dimensional case for q = 2. We suggest that in this case the model is close to or at the crossing point of lines separating three different types of phase transitions."
- d = q = 2 details:
  - "the jump in ρ is almost twice smaller than previously [q = 3] ... an observed jump might be a finite-size effect", with T_c = 1.7585(3).
  - ρ(t) at T_c decays with δ′ = 0.07(2), possibly 0. P(t) shows δ′ = 0.90(2) and η = 0, "the same values" as Hinrichsen's 2D q = 2 model.
  - "the possibility that we have a discontinuous transition accompanied by some dynamical power-law characteristics, although exotic at first sight, cannot be ruled out."
- Activity definition: "A given site i is active when at least one of its neighbours is in a state different than i."
- Footnote [14]: models with infinitely many absorbing states are excluded; "some models of this kind might exhibit novel critical behaviour."
- Relevance: **High**. Threat: **Medium–High**. It is the closest prior claim of a *possibly discontinuous* transition in a 2D Z2-symmetric two-absorbing-state model. The project should cite it and say what it adds: a clear first-order signal (fronts, droplets, L-independent lifetimes), as opposed to Lipowski–Droz's "inconclusive" status. — [arXiv:cond-mat/0112130](https://arxiv.org/abs/cond-mat/0112130)

**Grassberger, P. (2005)**, "Tricritical directed percolation in 2+1 dimensions", arXiv:cond-mat/0510428. The journal reference (believed J. Stat. Mech. 2006) was **not verified**.
- Read: **full text of the abstract, §I–II and the first-order section**.
- Model: generalized DK in 2+1D, non-tilted square lattice, von Neumann upstream sites, wetting probabilities P_k with P₁ = p, P_{k+1} = P_k + (1 − P_k)pq. Single absorbing state.
- "If, however, P_k increases fast with k so that the formation of inactive holes surrounded by active sites is suppressed, the transition is first order. These two transition lines meet at a tricritical point."
- In 1+1D "the cross over to the first order transition, called 'compact DP' in this case, happens at the boundary of the DK phase diagram, just as the phase transition in the 1-d Ising model occurs at T = 0".
- At the first-order line "clusters starting with single sites survive with a stretched exponential probability", and the boundary "behaves ... like a generic fluctuating non-equilibrium (KPZ) surface".
- Relevance: **High** (mechanism: compact domains, suppressed holes ⇒ first order in 2D). Threat: Low (single absorbing state). This is the "Grassberger 2D" item; I found no Grassberger 2D CA with *two* absorbing states. — [arXiv:cond-mat/0510428](https://arxiv.org/abs/cond-mat/0510428)

**Domany, E. & Kinzel, W. (1984)**, PRL 53, 311, and **Essam, J. W. (1989)**, J. Phys. A 22, 4927 (CDP). Read: citations only, via Grassberger's reference list. — [Grassberger refs](https://arxiv.org/abs/cond-mat/0510428)

**Bidaux, Boccara & Chaté (1989)**: see Q2 (first order for d > 1, one absorbing phase). The Chaté–Manneville collective-behaviour CA papers (EPL 14, 409, 1991; Prog. Theor. Phys. 87, 1, 1992) concern synchronous totalistic CA with quasi-periodic collective behaviour. They were found only as citations and are not verified or read. — [search context: Gallas et al. follow-ups on IDEAS/RePEc](https://ideas.repec.org/a/eee/phsmap/v180y1992i1p19-41.html)

### Inferences
- **What "2D DK endpoint" means for this project.** The literature has no paper on the 2+1D DK Z2-symmetric endpoint (P₀ = 0, P_max = 1, P_k = 1 − P_{z−k}). The project's model *is* such an endpoint family. The 1D CDP lore says β = 0 always, so "the density jumps 0 → 1" cannot by itself establish first order. Lübeck says so explicitly, which supports the LITERATURE.md critique of Rodrigues et al.'s magnetization-jump criterion.
- **The closest 2D precedents are both marginal or inconclusive.** Lipowski–Droz report a jump with ambiguous power laws; Hinrichsen reports roughly MF exponents and d_c ≤ 2. Together with Al Hammal et al. (d_c = 2 for GV, already in LITERATURE.md), they set expectations: in 2D, the two-absorbing-state transition sits near a boundary between continuous (GV/MF-like) and first-order behaviour. The project's claim is therefore plausible but must be phrased as "fully first order, with fronts, droplets and an L-independent interface density". A density jump alone would not do.
- **Mechanism.** Grassberger's mechanism (first order when holes inside active regions are suppressed) parallels Bidaux–Boccara–Chaté. For D&N, the active phase is a compact, rule-stabilised mixed phase with its own interface density; its "hole-suppression" analogue is the stability of that mixed phase against invasion by uniform domains. This is a natural way to explain first-order behaviour in the D&N model.

### Gaps
- No 2D (2+1D) study of the Z2-symmetric DK/CDP endpoint line was found.
- No "Grassberger 2D CA with two absorbing states" paper was identified; the brief's reference may conflate Grassberger 2005 (single absorbing state) or Hinrichsen 1997.
- Lübeck's own 2+1D tricritical contact-process paper (cited by Grassberger as a 2005 preprint) was not located.
- Chaté & Manneville 1991/1992 not read.

---

## Q4. How do majority- and minority-like dynamics freeze into non-uniform absorbing or blocked configurations (stripes, checkerboards), how do noise or voter terms remove them, and does that change the universality or order of the transition?

### Takeaway
Freezing into non-uniform configurations is generic for deterministic or zero-temperature majority-type rules.
- **Zero-temperature Ising (Glauber).** It freezes into straight stripes with probability ≈ 1/3 on square tori; the probability equals a critical-percolation crossing probability.
- **Odd coordination numbers** freeze "into a disordered state".
- **Synchronous majority or minority** ends in fixed points or 2-cycles (Goles–Olivos/Martínez).
- **Asynchronous minority** freezes into checkerboards (von Neumann) or stripes (Moore).
- **Remedies in the literature:** randomise the marginal or tie counts (Toffoli–Margolus RAND-ANNEAL; Vichniac's coin flip in checkerboard regions), or add any symmetric stochastic or voter-like perturbation (Bagnoli–Mencarelli; Gács: islands then shrink in O(n²)).

Frozen states mostly change *final-state statistics* and relaxation times, not the coarsening universality (curvature-driven). But Hinrichsen's rule (C) implies that frozen interfaces between the two absorbing phases would make the Z2 labels irrelevant and push the transition toward DP.

### Cited Findings

**Spirin, V., Krapivsky, P. L. & Redner, S. (2001)**, "Fate of zero-temperature Ising ferromagnets", Phys. Rev. E 63, 036118, arXiv:cond-mat/0008019.
- Read: abstract and first page. "On the square lattice, a frozen two-stripe state is apparently reached approximately 1/4 of the time ... In greater than two dimensions ... the system typically ends up wandering forever within an iso-energy set of stochastically 'blinking' metastable states."
- Relevance: High. Threat: none. — [arXiv:cond-mat/0008019](https://arxiv.org/abs/cond-mat/0008019)

**Spirin, V., Krapivsky, P. L. & Redner, S. (2002)**, "Freezing in Ising ferromagnets", Phys. Rev. E 65, 016119 (arXiv journal-ref gives 2002), arXiv:cond-mat/0105037.
- Read: abstract and §I–II.
- "In two dimensions, the system either reaches a frozen stripe state with probability ≈ 1/3 or the ground state with probability ≈ 2/3." "For small but non-zero temperature, relaxation to the final state first proceeds by the formation of very long-lived metastable states."
- "Any stripe of width ≥ 2 which traverses the entire system is obviously stable at zero temperature."
- They avoid "the pathologies associated with odd-coordinated lattices, where stable convex islands of minority spins readily form and the system always freezes into a disordered state".
- Blinkers are "localized sets of spins which can flip ad infinitum without any energy cost". Newman–Stein type M: "a fraction of the spins flip a finite number of times, while the complementary fraction flip infinitely often."
- Relevance: **High** (definitions of frozen, blinking and metastable). Threat: none. — [arXiv:cond-mat/0105037](https://arxiv.org/abs/cond-mat/0105037)

**Barros, K., Krapivsky, P. L. & Redner, S. (2009)**, "Freezing into stripe states in two-dimensional ferromagnets and crossing probabilities in critical percolation", Phys. Rev. E 80, 040101, arXiv:0905.3521.
- Read: abstract and first page. The system "eventually converges to either a ground state or an infinitely long-lived metastable stripe state". The stripe probability is predicted analytically from percolation crossing probabilities. "Our approach generally applies to coarsening dynamics of non-conserved scalar fields in two dimensions."
- Relevance: Medium–High. Threat: none. — [arXiv:0905.3521](https://arxiv.org/abs/0905.3521)

**Olejarz, J., Krapivsky, P. L. & Redner, S. (2012)**, "Fate of 2D Kinetic Ising Ferromagnets and Critical Percolation Crossing Probabilities", Phys. Rev. Lett. **109**, 195702, doi:10.1103/PhysRevLett.109.195702, arXiv:1208.2944.
- Read: abstract. The probability of a stripe state winding a times horizontally and b times vertically equals the critical-percolation crossing probability P_{a,b}. This **corrects the brief's "PRL 107, 195702 (2011)"**.
- Olejarz et al. also have 3D zero-T papers in PRE 83 (2011), cited in the arXiv text; not verified.
- Relevance: Medium. Threat: none. — [arXiv:1208.2944](https://arxiv.org/abs/1208.2944)

**Godrèche, C. & Pleimling, M. (2018)**, "Freezing in stripe states for kinetic Ising models: a comparative study of three dynamics", J. Stat. Mech. (2018) 043209, arXiv:1801.07749.
- Read: abstract. It compares zero-T Glauber with reversible and irreversible limits of a directed Ising model, testing whether percolation predictions for frozen stripes carry over.
- Relevance: Medium. It is evidence that freezing statistics depend on the dynamics, not only on the stationary measure. — [arXiv:1801.07749](https://arxiv.org/abs/1801.07749)

**Drouffe, J.-M. & Godrèche, C. (1999)**, "Coarsening and persistence in a class of stochastic processes interpolating between the Ising and voter models", J. Phys. A 32, 249, arXiv:cond-mat/9807356.
- Read: abstract. A 2D square-lattice family with interfacial and bulk "temperatures". "Special lines in the plane of parameters correspond to the Ising model, voter model and majority vote model." Persistence θ ≈ 0.22 in the coarsening phase.
- Relevance: Medium–High. It is the standard 2D interpolation between majority (frozen at T = 0) and voter, i.e. the analogue of the project's p. Threat: none. — [arXiv:cond-mat/9807356](https://arxiv.org/abs/cond-mat/9807356)

**Goles, E. & Olivos, J. (1980)**, "Periodic behaviour of generalized threshold functions", Discrete Math. 30(2), 187–189, doi:10.1016/0012-365X(80)90121-1.
- Read: citation only (AMS reference list). The content, that synchronous symmetric threshold or majority networks reach only fixed points or 2-cycles, is attested by Regnault et al., citing Goles & Martínez (1990, Kluwer).
- A secondary source notes a conflict over which Goles–Olivos paper (1980 vs 1981) holds the proof. — [AMS ref list](https://www.ams.org/tran/1995-347-05/S0002-9947-1995-1297535-1/); [Regnault et al. DMTCS](https://dmtcs.episciences.org/518)

**Asynchronous minority (Regnault et al. 2009/2010; Fatès 2014)**: frozen checkerboards (von Neumann) and stripes (Moore); empirical α_c ≈ 0.83 / 0.57; uniform states unreachable under full asynchrony for von Neumann minority. See Q2. — [arXiv:0706.2479](https://arxiv.org/abs/0706.2479); [DMTCS 518](https://dmtcs.episciences.org/518); [Fatès survey](https://arxiv.org/abs/1406.0792)

**Toffoli & Margolus (1987)** and **Vichniac (1984)**: noise or coin flips at marginal counts or checkerboard regions thaw frozen boundaries. See Q1. — [CAM book](https://people.csail.mit.edu/nhm/cam-book.pdf); [Physica D 10 scan](https://tcm.computerhistory.org/CHMfiles/Cellular%20Automation%201984.pdf)

**Bagnoli & Mencarelli (2026)**: frozen curvature-limited clusters are "unstable under stochastic perturbation, even when these obeys the same symmetry law of the original rule, like in the voter model"; synchronous checkerboard oscillations. See Q1. — [arXiv:2604.10082](https://arxiv.org/abs/2604.10082)

**Gács, P. (2023/2024)**, "Probabilistic cellular automata with Andrei Toom", arXiv:2306.13226 (v3, Nov 2024).
- Read: **full text of §3 and §5**.
- "In the absence of noise the majority rule over the von Neumann neighborhood ... would not shrink an island, say a square of size n, at all. Simulations (and informal arguments) show that in unbiased noise it will be eliminated in about O(n²) steps."
- "If the majority rule over the von Neumann neighborhood is perturbed in a way that is, say, biased towards the 1's, then a large island will grow!"
- "It is not known whether the majority rule over the von Neumann neighborhood, perturbed in a small but unbiased noise, can be ergodic."
- Relevance: High (Q4 and Q6). — [arXiv:2306.13226](https://arxiv.org/abs/2306.13226)

**Hinrichsen (1997)** rule (C): sticking interfaces between different absorbing phases make the transition DP. See Q3. — [arXiv:cond-mat/9608065](https://arxiv.org/abs/cond-mat/9608065)

### Inferences
- **D&N's frozen set is a local, bulk phase.** In zero-T Ising, frozen states are *spanning* (system-wrapping straight stripes). They have finite probability on tori but do not stop bulk coarsening on Z², so they change final-state statistics, not universality. In D&N at p = 0, the frozen set j ∈ {0,1,4,5} contains *local* bulk patterns: checkerboards and width-2 diagonal stripes, with j = 4 everywhere. This is closer to the odd-coordination pathology SKR describe ("stable convex islands ... freezes into a disordered state") and to asynchronous minority's stripe/checkerboard phase than to Ising stripes. The p = 0 model has extra absorbing *phases* (the two checkerboard registries, and diagonal-stripe states), not just extra absorbing *configurations*. That can change the transition, not just the final-state statistics.
- **Universality consequence.** By Hinrichsen's rule (C), if the boundary between an all-0 and an all-1 domain can freeze, the Z2 labels lose their dynamical role and the active→absorbing transition is expected to be DP-like (statics). Blok–Bergersen's Life (many frozen states, DP statics, nonuniversal spreading) is the empirical template. So p = 0 (frozen states allowed) and p > 0 (exactly two absorbing states) may lie in different classes; p > 0 is then a relevant perturbation. This is a testable prediction, not a literature result.
- **Measurement pitfall.** Lipowski–Droz's "active site" (any disagreeing neighbour) and any interface-density order parameter count a frozen checkerboard as *maximally active* although it is absorbing. At p = 0 the order parameter must be the flip activity, as in Blok–Bergersen and Regnault's "active cell". Otherwise frozen states masquerade as an active phase with high interface density.
- **What is known about removing frozen states.** Noise restricted to marginal counts (Toffoli–Margolus) or to checkerboard regions (Vichniac) is the historical fix. Symmetric noise at all inhomogeneous neighbourhoods (BBR's construction; Lipowski–Droz's "only fully agreeing sites are blocked") is the absorbing-state physics fix. Neither literature reports a change of *order* of the ordering transition caused by removing frozen states. That question is open, and the project could address it.

### Gaps
- No paper found on frozen states specifically for Moore-neighbourhood zero-temperature majority (j ≥ 5 rule) with its checkerboard 2-cycles under synchronous update; only von Neumann and Ising-lattice results are verified.
- Nanda–Newman–Stein (2000) and Gandolfi–Newman–Stein (2000) on fixation on Z² were not read. A secondary statement (via a Morris paper) is that at p = ½ on Z² no vertex fixates almost surely; unverified here.
- No paper quantifies how the *order* of an absorbing-state transition changes when frozen non-uniform states are removed by a small symmetric perturbation.

---

## Q5. What is known about models with infinitely many absorbing states, and about Z2-symmetric models with exponentially many frozen configurations? Is there a framing in which the p = 0 frozen states are a feature?

### Takeaway
The standard theory, Muñoz–Grinstein–Dickman–Livi (1996), covers models with infinitely many absorbing states and **no additional symmetry**: static exponents in the DP/Reggeon field theory class, with spreading exponents that depend on the frozen background (Jensen's PCP; Dickman 1999). Asynchronous Life (Blok & Bergersen) is a 2D CA instance with exactly that signature. **I found no theory for Z2-symmetric models with exponentially many frozen absorbing configurations.** Lipowski & Droz explicitly excluded them, noting such models "might exhibit novel critical behaviour". That gives a natural framing: the p = 0 line is a member of the "many-absorbing-states" family (expected DP-like if continuous), p > 0 restores exactly two absorbing states (GV or first-order), and the crossover in p is itself new territory.

### Cited Findings

**Muñoz, M. A., Grinstein, G., Dickman, R. & Livi, R. (1996)**, "Critical behavior of systems with many absorbing states", Phys. Rev. Lett. 76, 451–454.
- Read: repository record and abstract via search; full text not read.
- The paper proposes Langevin equations for models "with an infinite number of absorbing states and no additional symmetries"; the RG analysis places all their static exponents in the Reggeon field theory (DP) class. A companion paper appeared in Physica D 103, 485 (1997), cited by Lipowski–Droz's coauthors (unverified).
- Relevance: High (framing). Threat: none. — [FLORE record](https://flore.unifi.it/handle/2158/213341)

**Jensen, I. (1993)**, "Critical behavior of the pair contact process", Phys. Rev. Lett. 70, 1465 (DOI listing found). **Dickman, R. (1999)**, arXiv:cond-mat/9909347.
- Read: secondary summaries only. PCP statics are DP. Spreading exponents δ and η vary with the background particle density; Dickman reports that "z and δ+η ... exhibit a systematic dependence on φ" and confirms a generalized hyperscaling relation.
- Relevance: Medium. — [APS DOI](https://link.aps.org/doi/10.1103/PhysRevLett.70.1465); [arXiv:cond-mat/9909347](https://arxiv.org/abs/cond-mat/9909347)

**Blok & Bergersen (1999)**: DP statics and nonuniversal dynamics in α-asynchronous Life, attributed to "many possible absorbing states" (see Q2). — [PDF](https://www.cs.ubc.ca/~rikblok/research/published/blok99.pdf)

**Lipowski & Droz (2002)**, footnote [14]: "In the present paper we do not consider models with infinitely many absorbing states. In some cases typical universality classes are reported like DP or PC but it is also possible that some models of this kind might exhibit novel critical behaviour." Their (q, d) diagram: first order for d = 2, q ≥ 3 ("the jump is larger and the first-order character of the transition is even more transparent" for q = 4). — [arXiv:cond-mat/0112130](https://arxiv.org/abs/cond-mat/0112130)

**Hinrichsen (1997)** rule (C): if absorbing domains of different colour stick irreversibly, the model is in DP. See Q3. — [arXiv:cond-mat/9608065](https://arxiv.org/abs/cond-mat/9608065)

**Regnault et al. (2009)**: in α-asynchronous Moore minority, frozen stripe configurations form the ordered phase below α_c ≈ 0.57. See Q2. — [arXiv:0706.2479](https://arxiv.org/abs/0706.2479)

**Vieira et al. (2021)**: a freezing (absorbing) transition in conservative Life that is discontinuous. See Q2. — [arXiv:2101.07225](https://arxiv.org/abs/2101.07225)

### Inferences
Three framings in which the p = 0 frozen states are a feature, each grounded in the cited work:
1. **"Many absorbing states" line versus "two absorbing states" bulk.** At p = 0 the model belongs to the Muñoz/Jensen/Blok–Bergersen family. A continuous transition there would be expected to have DP statics and background-dependent spreading exponents. For any p > 0 only all-0 and all-1 are absorbing, so the transition is GV-type or first order. Then p measures the distance from a "frozen-manifold" line, and the project could report how the transition changes along p.
2. **Competing ferro and antiferro absorbing phases.** Because D&N leaves j = 4 frozen, the checkerboard (two registries, exchanged by the Z2 complement plus a lattice shift) and diagonal-stripe states are *bulk* absorbing phases alongside all-0 and all-1. At p = 0 the model then has q > 2 symmetry-related ordered absorbing phases. In Lipowski–Droz's (q, d) diagram, d = 2 with q ≥ 3 is the first-order region. This is speculative, because their q counts permutation-equivalent ferro states. But it suggests frozen states could *favour* discontinuity at p = 0. Measuring the weight of checkerboard and stripe domains near the transition would test it.
3. **Asynchrony-induced ordering, as in Regnault's minority.** One can view the p = 0 frozen patterns as the "stripe/checkerboard phase" of an α-asynchronous self-complementary rule. α and p then act as two routes to melting it, and the user's "first-order" claim concerns the competition between the active mixed phase and the uniform phases once the frozen phase is melted.

Whichever framing is used, the claims must be stated separately for p = 0 (infinitely many absorbing states) and p > 0 (exactly two). An observable such as interface density cannot be used at p = 0 without the flip-activity check (Q4).

### Gaps
- No paper found on Z2-symmetric models with exponentially many frozen absorbing configurations (a "Z2 PCP"). Ódor's review (Rev. Mod. Phys. 2004) may cover this but was not checked.
- Muñoz et al. (1996) full text not read; the scope of "no additional symmetries" was not examined in detail.

---

## Q6. Are Toom's NEC rule, eroders, Gács's work and the positive-rates conjecture relevant to two phases coexisting over a parameter region?

### Takeaway
- **Toom.** His NEC rule (2D, majority of north, east and centre) is the classic proof that a noisy PCA with *all transition probabilities positive* can be non-ergodic. Its two phases coexist over a **finite area of parameter space** (generic bistability, against the Gibbs phase rule) because the deterministic rule and its dual are *eroders*: minority islands shrink at constant speed, so even biased noise cannot nucleate the other phase.
- **Isotropic majority** (von Neumann) is not an eroder: squares are frozen without noise, shrink only in O(n²) under unbiased noise, and *grow* under bias. So isotropic Z2 rules generically have coexistence only on the symmetry line.
- **Gács** extended non-ergodicity to 1D (counterexample to the positive-rates conjecture).
- **Relevance.** A coexistence *window* in the project's parameters is allowed in nonequilibrium systems (Toom; Bennett–Grinstein; Muñoz et al. 2005; Godrèche–Luck 2017). But in those works the window comes from anisotropy or irreversibility of interface motion. D&N is isotropic and non-monotone (Toom's eroder theorem needs monotonicity), so a finite window must come from a different mechanism, such as pinned or zero-velocity fronts.

### Cited Findings

**Gács (2023/2024)**, arXiv:2306.13226.
- Read: full text of §3 and §5.
- On positive rates: "the conjecture was considered that a probabilistic cellular automaton with no prohibited transitions ... is always ergodic ... Toom's best-known result is a refutation of this conjecture."
- The NEC rule is "a majority vote over three neighbors: north, east, self ... The novelty is just that the neighborhood over which the majority is taken is not (central-) symmetric—but this makes a big difference." An enclosing triangle shrinks each step, so "the island will shrink with constant speed. Small noise (even if it is biased!) can slow down this shrinking, but still leaves its speed constant".
- Eroders: "Toom called a deterministic cellular automaton an eroder if it eliminates any island". Toom's convex-hull criterion decides this for **monotonic** binary CA (Theorem 1). "Every binary eroder is stable" (Theorem 2, from Toom 1980).
- In 1D, ergodicity of monotonic two-state nearest-neighbour rules is due to Gray (1982, 1987). Gács (2001), J. Stat. Phys. 103, 45–267, is the 1D non-ergodic construction.
- Toom 1980 is cited there as "Stable and Attractive Trajectories in Multicomponent Systems", in *Multicomponent Systems* (Dobrushin, ed.), Advances in Probability 6, 549–575, Dekker.
- Relevance: **High** (conceptual). Threat: none. — [arXiv:2306.13226](https://arxiv.org/abs/2306.13226)

**Bennett, C. H. & Grinstein, G. (1985)**, "Role of irreversibility in stabilizing complex and nonergodic behavior in locally interacting discrete systems", Phys. Rev. Lett. 55, 657–660.
- Read: IBM record and abstract summary via search, plus the citation in an Astérisque reference list; full text not read.
- Irreversibility (NEC) lets two phases coexist nonergodically over a positive-measure region of parameter space; reversible analogues do not.
- Relevance: High. — [IBM record](https://research.ibm.com/publications/role-of-irreversibility-in-stabilizing-complex-and-nonergodic-behavior-in-locally-interacting-discrete-systems); [Astérisque 287 (Numdam)](https://www.numdam.org/item/AST_2003__287__71_0/)

**Muñoz, M. A., de los Santos, F. & Telo da Gama, M. M. (2005)**, "Generic two-phase coexistence in nonequilibrium systems", Eur. Phys. J. B 43, 73–79, doi:10.1140/epjb/e2005-00029-3, arXiv:cond-mat/0411578.
- Read: abstract and introduction.
- "Toom's cellular automaton, revealed that non-equilibrium systems may exhibit generic bistability, i.e. two-phase coexistence over a finite area of the (two-dimensional) phase diagram, in violation of the equilibrium Gibbs phase rule." They analyse "two interfacial models" (pinning–depinning, KPZ-type) with broad coexistence regions.
- Relevance: **High** for any claimed coexistence window. Threat: Low. — [arXiv:cond-mat/0411578](https://arxiv.org/abs/cond-mat/0411578)

**Godrèche, C. & Luck, J.-M. (2017)**, "Generic phase coexistence in the totally asymmetric kinetic Ising model", J. Stat. Mech. (2017) 073208, arXiv:1704.06120.
- Read: abstract. A 2D totally asymmetric kinetic Ising model whose zero-field stationary measure is the Ising one. Droplet shrinking velocity v(T, h) vanishes at boundary fields ±h_b(T), "which bound the coexistence region", with h_b vanishing near T_c with exponent Δ_h ≈ 3.2 ± 0.3.
- Relevance: Medium–High (a quantitative method: coexistence window = interval where droplet or front velocity has a fixed sign). — [arXiv:1704.06120](https://arxiv.org/abs/1704.06120)

**Grinstein, G., Jayaprakash, C. & He, Y. (1985)**, "Statistical mechanics of probabilistic cellular automata", Phys. Rev. Lett. 55, 2527, doi:10.1103/PhysRevLett.55.2527.
- Read: APS abstract (paraphrased there). Continuous ferromagnetic transitions in up-down-symmetric PCA (positive rates) belong to the kinetic Ising class even for irreversible rules. Periodic states seen in mean field may persist in 2D and 3D.
- Relevance: Medium. It is the baseline: Z2 PCA without absorbing states give Ising; with absorbing states, GV, DP2 or first order. — [APS abstract](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.55.2527)

**Marsan, Sablik & Törmä**, "A Perturbed Cellular Automaton with Two Phase Transitions for the Ergodicity" (J. Stat. Phys., per search). Read: search snippet only; journal details unverified. It describes a CA that is ergodic at small noise, non-ergodic at intermediate noise and ergodic again near 1, using Gács's construction. Relevance: Low. — [UTU record](https://research.utu.fi/converis/portal/detail/Publication/508988842)

### Inferences
- **Toom's theory does not apply to D&N directly.** The project's model is isotropic, Z2-symmetric by construction, and D&N is non-monotone (it flips at j = 3 but not at j = 4, 5). Toom's eroder theorem, which needs monotonicity, therefore does not apply. Toom/Bennett–Grinstein coexistence arises from *asymmetric* (constant-speed) erosion of either phase's islands.
- **Where a coexistence window in the project could come from.** It concerns the active mixed phase versus the uniform absorbing phases (the window [κ_front, κ_a] in LITERATURE.md). It would require that, over a parameter interval, neither phase's large droplets shrink. Muñoz–de los Santos–Telo da Gama show the interfacial routes: pinning, or KPZ-type nonlinearity with a front-orientation dependence. Godrèche–Luck give the measurement recipe: map v(params) for large droplets and define the window by v = 0 on both sides.
- **A caution against over-claiming.** In isotropic models, a finite-size "window" between front-reversal and nucleation thresholds is usually a metastability (hysteresis) window that shrinks with system size and time. True generic coexistence needs an erosion mechanism. The project should state which of the two it claims.

### Gaps
- Bennett & Grinstein (1985) and Toom (1980) full texts not read. Gray (1999) on continuous-time Toom stability was not read.
- No paper found on eroder-like properties of non-monotone isotropic rules such as D&N.

---

## Q7. Synthesis: how have others obtained exactly two absorbing states, and how could the D&N model be reformulated?

### Takeaway
The absorbing-state literature obtains exactly two absorbing states with one recipe: **every inhomogeneous neighbourhood has a positive flip probability; only fully agreeing cells are blocked.**
- BBR: "absorbed or desorbed only near the boundaries ... (when the neighborhood is not homogeneous)".
- Lipowski–Droz: Metropolis except "when all neighbours ... are in the same state".
- Toffoli–Margolus and Vichniac randomise only the marginal or checkerboard situations.

For the Moore lattice there is a short argument (my own, below) that **w(j) > 0 for all j ≥ 3 already suffices**. D&N at p = 0 fails only because w(4) = w(5) = 0. So adding a small symmetric flip probability ε at j ∈ {4, 5} (a RAND-ANNEAL-style marginal noise) restores exactly two absorbing states for all p, with minimal change to the rule.

### Cited Findings
- BBR construction of exactly two absorbing states (positive rates at all inhomogeneous neighbourhoods): see Q3. — [arXiv:cond-mat/0002361](https://arxiv.org/abs/cond-mat/0002361)
- Lipowski–Droz: "when all neighbours of a given site are in the same state as this site, then this site cannot change its state", which yields q absorbing states for T > 0. — [arXiv:cond-mat/0112130](https://arxiv.org/abs/cond-mat/0112130)
- Hinrichsen rule (C): interfaces between different absorbing phases must not stick irreversibly, or the transition reverts to DP. — [arXiv:cond-mat/9608065](https://arxiv.org/abs/cond-mat/9608065)
- RAND-ANNEAL: noise only at the two marginal counts; "stuck patterns are eventually melted away". — [CAM book](https://people.csail.mit.edu/nhm/cam-book.pdf)
- Vichniac's coin flip in checkerboard domains: "puts these regions at an infinite temperature", and T = 0 domains grow "at the expense of the dislocations and the oscillatory monsters". — [Physica D 10 scan](https://tcm.computerhistory.org/CHMfiles/Cellular%20Automation%201984.pdf)
- Symmetric stochastic perturbations destabilise frozen majority clusters "like in the voter model". — [arXiv:2604.10082](https://arxiv.org/abs/2604.10082)
- Under full asynchrony, von Neumann minority cannot reach its two "uniform" fixed points (presumably the checkerboards; see Q2) from non-uniform states, but small synchrony restores reachability. The update scheme therefore matters for which absorbing states are reachable. — [Fatès survey](https://arxiv.org/abs/1406.0792)

### Inferences
**Lemma (my derivation; not from the literature; worth a line in the paper).** On the Moore lattice, every non-uniform configuration contains a cell with j ≥ 3.
- Proof sketch. A non-uniform configuration has an orthogonally adjacent disagreeing pair x, y, because the nearest-neighbour graph is connected. x and y share four Moore neighbours, the two cells adjacent to both on each side. Each shared cell disagrees with exactly one of x, y. Counting also the mutual disagreement, j(x) + j(y) ≥ 2 + 4 = 6, so max(j(x), j(y)) ≥ 3.
- Consequence. In a synchronous PCA, a configuration is absorbing iff every cell has w(j) = 0. So **w(j) > 0 for all j ∈ {3, …, 8} implies the uniform states are the only absorbing states**, whatever w(1) and w(2) are.
- The condition is close to necessary:
  - w(3) = 0 freezes flat stripes (j ∈ {0, 3} everywhere).
  - w(4) = 0 freezes checkerboards (j = 4 everywhere).
  - w(6) = 0 freezes width-1 stripes (j = 6 everywhere).
  - Whether w(5) = 0 alone admits frozen configurations was not determined. It needs an exhaustive check on small tori, e.g. a SAT search for configurations with all j ∈ {0, 1, 2, 5}.
- **Caveat on closed classes.** The lemma rules out non-uniform absorbing *configurations*. It does not by itself rule out closed recurrent *sets* of non-uniform configurations, such as the stochastic "blinker" sets SKR find in d ≥ 3 zero-T Ising.
  - For p > 0 the voter term gives every j ≥ 1 a positive rate, which should make consensus reachable from any configuration on a finite torus.
  - At p = 0 with only j ≥ 2 or j ≥ 3 flips allowed, reachability of consensus from every configuration should be checked (proved or verified numerically).

**Concrete reformulations, in order of how little they change D&N:**
- (a) Add ε·[j ∈ {4, 5}] symmetric flips: RAND-ANNEAL-style marginal noise. This gives exactly two absorbing states for every p ≥ 0 when α(1 − p) > 0 or p > 0. ε → 0 recovers D&N, so results can be quoted as ε → 0⁺ limits.
- (b) Add ε·[j ≥ 1]: BBR/Lipowski–Droz-style positive rates at every inhomogeneous neighbourhood. This is the most standard choice in absorbing-state physics.
- (c) Keep the model but restrict all claims to p > 0, and treat p = 0 as the "many absorbing states" boundary line (Q5 framing 1).
- **Caution.** Flips at j = 5 push the rule toward majority/Ising-like curvature-driven dynamics (Toffoli–Margolus MAJORITY versus ANNEAL). ε should be reported and its effect on front velocities checked: does the sign of v change as ε → 0?

**Diagnostics to report, following the literature:**
- Flip activity, as opposed to interface density, as the order parameter at small p or ε (Blok–Bergersen; Regnault et al.).
- Absence of deterministic cycles. A deterministic cycle needs every cell to have w(j) ∈ {0, 1}. At p = 0, w(j) ≤ max(α, q), so for α < 1 and q < 1 no deterministic cycles exist. At α = 1, p = 0 the rule is deterministic on part of configuration space. D&N is not a threshold rule, so Goles–Olivos does not cap the period at 2.
- A demonstration that the first-order signal is not an artefact of a frozen-configuration "jump", given the CDP/Lipowski–Droz caveats on density jumps (Q3).

### Gaps
- Whether w(5) = 0 alone (with w(4) > 0) admits non-uniform frozen configurations on the Moore lattice is not settled here; it needs a small exhaustive or SAT search.
- No published model was found that adds marginal noise to D&N or Anneal and studies the resulting absorbing-state transition. The closest is the unanalysed RAND-ANNEAL.
