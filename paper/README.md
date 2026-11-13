# Manuscript: discontinuous transition into two symmetric absorbing states in 2D

Draft for Phys. Rev. E (revtex4-2). Red `\pending{...}` marks values that still have to be filled
in. Remove the macro before submission; `grep pending main.tex` should then return nothing.

## Build
No TeX is installed on this machine. Any TeX Live with revtex4-2 works, or Overleaf. Locally:

    conda create -p /path/to/texenv -c conda-forge --override-channels tectonic
    /path/to/texenv/bin/tectonic -X compile main.tex

## Figures
`figs/fig_*.py` write `figs/fig_*.pdf` and `.png`. Run them from `figs/`. They read:
- `../data/runs/<EXP>/out/*.txt`: time series mirrored from Torch (`$SCRATCH/dn/runs`, the
  `.txt` files only):

      rsync -a --include='*/' --include='*.txt' --exclude='*' --prune-empty-dirs \
          torch:/scratch/mo2273/dn/runs/<EXP> ../data/runs/

- `../data/figdata/`: reductions of the large `.prof` files, made on Torch inside the
  container by `scripts/extract_prof.py`:

      python3 extract_prof.py $SCRATCH/dn/runs $SCRATCH/dn/figdata fronts E1_front E4_voter ...
      python3 extract_prof.py $SCRATCH/dn/runs $SCRATCH/dn/figdata walls E8b_wall E16_wall ...

- `../campaign/theory/vpt_run_th3.txt`: the perturbation theory (`campaign/theory/vpt.py`); `../campaign/theory/molofsky_1e7.csv`: the same closure for the five-cell rule (`vpt.py molofsky`, see `campaign/theory/MOLOFSKY.md`). `../campaign/molofsky_sim/`: the direct simulation of that rule (`dnsim.c` v1.4 with `--molofsky`, front tables from `mfront.py table` and, for the fine scan across the reversal, `sec8_tables.py`; the seed-averaged densities at the line in `bulk_dens/`, `mdens.py`; all results in its `NOTES.md`), used for Table IV and Fig. 10(c).
- `../data/snap/*.ckpt`: dnsim checkpoints for the snapshot figure. They come from local runs at
  p = 0.5, α = 0.5, with q = 1.05, 0.95, 0.98 and 1.1 × 0.0459; the commands are in the
  caption of Fig. 1.

| figure | script | data |
|---|---|---|
| 1 snapshots | fig_snap.py | data/snap |
| 2 orientation / fronts | fig_orient.py | R2_orient (direct front tracking, L 4096; `data/figdata/R2/orient_runs.csv` from `campaign/round2/orient2.py`); the E15 interface-count version is kept as fig_orient_e15.py |
| 3 two-droplet test | fig_drop2.py (uses campaign/round2/drop2.py) | R2_drop2 (time series and row profiles in data/runs) |
| 4 bulk at p = 1/2 | fig_bulk05.py (bulk05.py) | E13_bulk05 |
| 5 critical droplets | fig_drop.py | E14_drop05, E4_voter / E12_phase fronts; R = 0.75 and the full 2×2 block from R2_seed22 (4000 replicas) |
| 6 wetting | fig_wall05.py | E16_wall and R2_wall (figdata; Ly 2048 and 4096 series, runs to 2–3e6 steps; the Lx 4096 seeds 1, 2 are the E16 runs continued) |
| 7 roughness of the layer edges | fig_rough.py | `data/figdata/R2/rough_*.csv`, `whist_*.csv` (campaign/round2/rough.py on the R2_wall and R2_orient snapshots) |
| 8 phase diagram | fig_phase.py (fronts_util.py) | E1, E4, E9, E12 fronts (p ≥ 0.99: the long strips of `../campaign/strip_p099/`, values in the script); E2, E4b, E6, E12b, E13 bulk and `../campaign/strip_p099/bulk` |
| 9 near the voter point | fig_voter.py | E6/E7 (χ table in the script, from `campaign/binder.py` with m = 2n − 1; output kept in `fig_voter_chi.txt`), E8a, E8a_long, E8b, E16_wall_ref; legends in (b) and (c) give q/q* relative to the Table I thresholds, the runs having been placed at multiples of 0.00127 (p = 0.99) and 0.00062 (p = 0.995) |
| 10 theory | fig_theory.py | campaign/theory |
| 11 single-noise rule | fig_eta.py | E17_eta_front, E17_eta_bulk |
| 12 Molofsky rule, closure and simulation (Appendix D) | fig_molofsky.py | campaign/theory/molofsky_1e7.csv; campaign/molofsky_sim/scan/table_p2_*.txt, scan_wide/table_p2_*.txt; inset: scan_fine/table.txt, scan_fine/wide/table.txt |

The bibliography is `refs.bib`. Every entry was checked against Crossref, arXiv or zbMATH; the
notes are in `refs_check.md`. Entries that nothing in `main.tex` cites were moved to
`refs_removed_2026-10-31.bib` (2026-10-31).
