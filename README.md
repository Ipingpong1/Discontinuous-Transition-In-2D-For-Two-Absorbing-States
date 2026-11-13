# dnsim — bit-packed stochastic Day & Night

Synchronous stochastic outer-totalistic CA on an Lx×Ly torus, 64 cells per word.
A cell in state s with k live neighbours flips with probability f[s][k]; built from

    f = (1-voter)·alpha·[rule flips] + voter·j/8 + Σ q_j'·[j = j']     (j = # disagreeing nbrs)

so `--rule B3678/S34678 --alpha 0.5 --jflip 2:0.02` is D&N with α-asynchrony plus a
symmetric minority flip. The resolved 18-entry table is printed in every output header.

    make && ./dnsim --selftest
    ./dnsim --L 1024 --T 20000 --alpha 0.5 --jflip 2:0.02 --seed 1 --out run.txt
    ./dnsim --Lx 1024 --Ly 2048 --T 20000 --alpha 0.5 --jflip 2:0.019 \
            --init strip:0.5:0.5 --bg 0 --out s.txt --profile s.prof
    python3 profile.py s.txt s.prof [rho_if]      # active width A(t) → front velocity

Output columns: `t alive_frac rho_bond flip_frac iface_frac` at log-spaced times
(`--ppd`, plus `--every n`). `--save/--load` checkpoints; the time continues, and the
RNG is counter-based (Philox4x32-10 keyed by seed, counter = t, row, word, draw), so a run
is bit-identical regardless of thread count or restarts. Each cell draws one uniform per
step and flips iff U < f of its class; all classes of a word share one bit-serial
comparison, so the cost does not grow with the number of distinct flip probabilities.

v1.2 (10-14 evening) replaced v1.1's one-stream-per-class sampler (counter word|value<<24)
with the shared comparison above. The law of the dynamics is unchanged and only the random
streams differ. Initial conditions are bit-identical. Phase 1–2 data are v1.1, kept as
`dnsim_v11` locally and `~/dn/bin/dnsim_v11` on Torch. Check: ρ_bond at p=0, q=0.023, L 256,
120 seeds each: v1.1 0.30742(12), v1.2 0.30727(11). Five voter-mixed sets agree within
|z| < 1.5.

v1.3 (10-15) adds `--reps`, the spreading mode for E8a. Normal-mode output is byte-identical
to v1.2: time series, profile and checkpoint, checked on bulk, D&N and strip runs. The
per-word update is now a shared always-inline function, at the same speed (−0.6%).

    ./dnsim --alpha 0.5 --voter 0.99 --jflip 2:0.0013 --L 2048 --T 100000 \
            --init disk:1:0.75 --reps 3000000 --target 300 --seed 7 \
            --out sp.txt --ckpt sp.ckpt --resume
    python3 campaign/spread.py sp.txt

- Each replica starts from `--init disk:p0:R` (the 2×2 block for disk:1:0.75) on an all-dead
  torus.
- Only the bounding box of the cells with a disagreeing neighbour is stepped. This is exact:
  every other cell has j = 0 and f(0) = 0, so it cannot flip, and the RNG counter depends only
  on (t, row, word). The box grows next to flips and is refitted every 32 steps.
- A replica whose box comes within 2 rows or 1 word of the lattice edge is stopped as
  *escaped*. It counts as alive from then on but adds nothing more to the sums.
- Replica i uses key (seed lo, seed hi ^ 0x9E3779B9·(i+1)) and adds into group i % 16.
- For every output time, `--out` gets per-group sums: replicas, surviving, escaped, I, I², A,
  A², Σr² over interface cells, and box area.
- The run stops after `--reps` replicas, or once `--target` replicas are alive at T.
- Checkpoints go between replicas (`--ckpt-every` counts replicas, default 1000, plus every
  10 min). The results do not depend on where a run was interrupted: two interruptions gave
  byte-identical data.

## Validation (2026-10-14)
- `--selftest`: Philox known-answer vectors; Bernoulli frequencies; bit-exact agreement
  with a naive per-cell implementation for 5 rules × 4 lattice shapes (incl. W=1, Ly=3);
  per-class flip frequencies for an 18-distinct-value table (z-variance over 560 repeats
  1.04 ± 0.06); restart and thread-count determinism. v1.2 adds a 134M-sample test of the
  shared sampler on a voter-like 8-class table: per-class frequencies and adjacent-lane
  independence, both |z| < 2.
- Reproduces the independent torch/GPU pilot: steady-state ρ in the disordered phase
  agrees within ~1% (q = 0.02: 0.266 vs 0.267 at L=512; 0.271 vs 0.273 at L=1024).
- v1.3 `--selftest` adds box stepping vs full stepping on 200 seeds (D&N + voter + j=2 noise;
  65 died, 102 escaped, 33 alive at t = 600). It compares grids, incremental live counts and
  measure_box against measure() every 25 steps: 0 mismatches. seed_disk equals init_grid's
  disk.
- Pure voter from a 2×2 block (20,000 replicas, T 1e4, L 1024): ⟨A(t)⟩ = 4.00 ± 0.01 → 3.3 ± 0.6
  (martingale, within errors); R² ∝ t (z_eff 0.93–1.01); δ_eff rises to ~1 and t·P_s/ln t falls
  slowly from 2.4 to 1.6 (Bramson–Griffeath log corrections). 0.54 ms per replica.

## Speed (Ryzen 7600X3D, AVX-512)
Near the transition (RNG-bound): 1.1e9 cell updates/s/core, 7.1e9 on 12 threads.
Deterministic D&N: 5.9e9/s/core.
v1.2, one core, L 1024, active phase: D&N q=0.0185 1.06e9 → 1.28e9; voter-mixed (p 0.5–0.99)
2.2–2.9e8 → 0.92–1.1e9. On the Torch login node the voter-mixed rate goes 1.09e8 → 4.33e8.

## Repository contents

This repository tracks the simulator, campaign and analysis source, HPC job
scripts and manifests, research notes, and manuscript sources and figures.
Large simulation datasets (`data/` and campaign run outputs), generated reports,
release archives, compiled executables, caches, and local tool settings remain
on the original workstation and are excluded by `.gitignore`.

Build the simulator with `make` and run its checks with `./dnsim --selftest`.
Reproducing existing campaign analyses and manuscript figures also requires the
corresponding simulation data; those datasets are not included in this source
checkout. Campaign output files are ignored by default; explicitly allow any
small curated dataset in `.gitignore` before adding it to version control.
