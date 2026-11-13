# In-house checks on the model (2026-10-27, local runs with dnsim v1.3)

These are not literature findings. They are quick computational checks run while the literature
search was in progress, so the report can address the user's concern (A): "the p = 0 flagship
has more than two absorbing states". All runs are single-seed, local, and meant as diagnostics,
not publication data. Scripts and checkpoints are in the session scratchpad (`absorb/`).

## 1. How many absorbing configurations does the p = 0 model have?
At p = 0 the flip table is w(j) = α·[j ∈ {3,6,7,8}] + q·[j = 2]. A configuration is absorbing
(frozen) iff every cell has j ∈ {0,1,4,5} (q > 0) or j ∈ {0,1,2,4,5} (q = 0).
- Transfer-matrix count on W × L tori: the number of frozen configurations grows exponentially
  with the area. Strip entropy ln λ_max / W ≈ 0.32 nats per site (W = 7, 8), against ln 2 = 0.69
  for all configurations. So the p = 0 model has exponentially many absorbing configurations,
  not just all-0, all-1, the 2 checkerboards and the 8 diagonal width-2 stripe states (all of
  which have j = 4 at every cell).
- Check: on the 3 × 3 torus (where every cell neighbours every other) the count is 254 =
  1 + C(9,4) + C(9,5) + 1, as it must be.

## 2. Are they dynamically relevant?
- Random start, p = 0, α = 0.5, L = 256, T = 1e6, q = 0.005, 0.012, 0.016, 0.018 (all below
  q* = 0.0184; 2 seeds each): every run ends exactly in all-0 or all-1. At q = 0 (pure
  α-asynchronous D&N) the run freezes with ≈1% interface cells of still-life debris in a
  uniform background.
- Active phase (q = 0.03–1, α = 0.5 and 1): checkerboard-like and diagonal-stripe-like 3 × 3
  blocks stay at chance level (0.2–0.6% and 1–1.5%).
- A prepared half-lattice of checkerboard or diagonal stripes, next to an active half or an
  all-dead half, is completely eroded within 1e4 steps at every q tested (0, 0.004 at α 0.125,
  0.01, 0.016, 0.02, 0.05, 0.2).
- A single flipped cell in a perfect 256 × 256 checkerboard or stripe pattern destroys the
  whole pattern within ~1e3 steps (11 of 12 seeds; in one the defect healed).
- Conclusion: the extra absorbing states are real but locally unstable, with zero basin of
  attraction. For q > 0 only all-0 and all-1 are attractors, so the p = 0 data are not wrong.
  But the statement "two absorbing states" is formally false at p = 0, and a referee would
  catch it.

## 3. Variants with exactly two absorbing states (any w(j) > 0 for all j ≥ 1)
- Any p > 0 (voter admixture, w(j) ≥ p·j/8) already qualifies. All Phase 2–7 runs
  (p = 0.5 … 0.995) are in this clean family.
- ε-regularised p = 0 (add flip probability ε on j ∈ {1,4,5}): with ε = 0.003 the horizontal
  front still reverses linearly. Per-step v = −0.0041, −0.0002, +0.0024, +0.0053 at
  q/q*(0) = 0.88, 0.92, 0.96, 1.00, so q*(ε) ≈ 0.92 q*(0), with a finite interior interface
  density (0.44–0.6) at the crossing. At ε = 0.01, q* ≈ 0.8 q*(0). q*(ε) is continuous as
  ε → 0, so the p = 0 results are the ε → 0+ limit of a clean model.
- Single-noise variant: w(j) = α·[j ∈ {3,6,7,8}] + η·[j ∈ {1,2,4,5}]. This has two parameters,
  exactly two absorbing states, and the noise acts only on interface cells. At α = 0.5 the
  front reverses between η = 0.012 (v ≈ −0.0027; the band interior decays) and η = 0.013
  (v ≈ +0.0013, interior interface fraction 0.60–0.67), and v = +0.0062 at 0.014. So η* ≈ 0.0125,
  with a large interface-density jump: a strongly first-order candidate flagship.

## 4. Orientation dependence of the front threshold (generic-coexistence check)
At p = 0.5, α = 0.5 (exactly two absorbing states; horizontal-front q* = 0.0459 from E4), horizontal
and diagonal (45°) bands of the same normal width were run (L = 1024, one seed). Normal front
speeds were v = −0.0056 / −0.0049 (x = q/q* = 0.98), +0.0007 / +0.0014 (1.00), +0.0041 / +0.0051
(1.02), +0.0098 / +0.0095 (1.04), horizontal / diagonal. Both reverse between x = 0.98 and 1.00.
So no orientation dependence of q* shows above ~1% at this point. A finite generic-coexistence
window narrower than ~1% is not excluded, and nothing was checked at p = 0 or near p = 1.
