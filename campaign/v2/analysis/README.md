# V2 numerical analysis

These modules leave all original scripts and inputs unchanged. Run from the project root:

```sh
python3 -B -m unittest campaign.v2.analysis.test_analysis -v
python3 -B -m campaign.v2.analysis.bulk --boot 1000
python3 -B -m campaign.v2.analysis.wetting --boot 1000
python3 -B -m campaign.v2.analysis.fronts data/runs/R2_orient/out data/runs/R2_orient2k/out --boot 500
python3 paper/v2/figs/fig_bulk05.py
python3 paper/v2/figs/fig_wall05.py
python3 paper/v2/figs/fig_orient.py
```

For new nested GPU outputs, the same readers accept `task/out.txt` and
`task/profile.bin`, using actual `# profile_dir` and model header values. Example:

```sh
python3 -B -m campaign.v2.analysis.fronts data/v2/gpu/runs/orientation --out campaign/v2/results/fronts_new --tmins 100000,200000,500000,1000000 --windows 100000:500000,500000:1000000,1000000:2000000
python3 -B -m campaign.v2.analysis.bulk data/v2/gpu/runs/flagship_bulk --out campaign/v2/results/bulk_new
```

Front reductions skip unfinished jobs by default. `--include-incomplete` is for
diagnostics and retains an explicit completion flag. Fixed-window front fits and
fits extending from different start times to the run end are reported separately.

Only NumPy is required for reductions/tests; figure scripts also require Matplotlib.
Outputs go to `campaign/v2/results/{bulk,wetting,fronts}`. Each contains per-run and
grouped tables plus JSON settings and SHA-256 input provenance. Figures read these
outputs, rather than recomputing selected subsets with separate definitions.

## Definitions and limits

* Bulk late means integrate the time series over `[max(1000,T/2),T]`. Every completed
  run contributes. Consensus is extended by its known absorbing state to T. Separate
  survivor-conditioned means are labelled as such; no record is filtered by density.
* A first recorded half-density crossing and complete consensus are different events.
  Both the preceding record and the first-below time are retained. Because activity
  is nonmonotone, an earlier unrecorded dip is possible: this pair is not a rigorous
  bracket on the first physical crossing. Product-limit curves use the first recorded
  crossing. Consensus is absorbing, so its observation bracket has no such ambiguity.
  Unfinished runs are listed but not
  mixed into completed-run summaries. RMST is an observed-horizon quantity, not an
  estimate of an uncensored mean lifetime. Event-free groups receive no bootstrap
  RMST interval because resampling would give misleading zero uncertainty.
* Fronts use the minimum and maximum qualifying bins of the smoothed profile. A hole
  inside the strip cannot stop this envelope. The old contiguous-component tracker is
  retained only as a diagnostic. Smoothing has the same physical normal length for
  each orientation. Run slopes, hole fractions, interior density, time-window choices,
  local interpolation, polynomial lack of fit and seed-bootstrap intervals are explicit.
  Crossing intervals condition on model/window; their spread is a separate systematic.
* Wetting fits use residuals in Gamma units for both candidate two-parameter models.
  Points sharing a q=0 baseline have correlated errors. Baselines and independent runs
  are resampled separately; q* is varied as a systematic range. Time windows begin at
  max(3e5,T/2), and terminate before a drifting wall consumes an absorbing slab
  (`abs(alive-.5)>=.2`). All truncations and retained fractions are reported. Such
  conditional wall windows must be checked against larger boxes.

The final operational p=.5 reference and its normalization-sensitivity range are
recorded in `results/coexistence_reference.json`. They summarize finite-horizon
front reversal measurements with explicit direction, window and tracker
dependence; the range is neither a pooled confidence interval nor a rigorous
bound on a limiting coexistence point. Small numbers of independent seeds limit
both standard errors and bootstrap coverage.

## Fixed GPU ensembles

Use the manifest-aware reduction for the new campaign:

```sh
python3 -m campaign.v2.analysis.gpu_campaign --boot 1000
python3 paper/v2/figs/fig_gpu_v2.py
```

The reduction merges `production_manifest_profiles50k.json`,
`production_manifest_largeL.json`, `production_manifest_near_large.json`,
`production_manifest_orientation5m.json`, `production_manifest_near_midpoint.json`,
and `production_manifest_axial10m.json`.
The exact-seed repair replaces one previously ineligible original task only
after the replay checks described below; it does not add another realization.
Every declared task appears in `results/gpu/status.json`; complete groups are
listed in `ensemble_eligibility.csv`. Publication summaries and curves require
all planned independent seeds. A completed text footer alone is insufficient:
text end, all profile times, terminal profile record, and `done.json` must agree.
Incomplete live copies and quota-interrupted files remain excluded. Front
crossing fits additionally require a valid geometrical slope for every declared
seed at the specified parameter and window; invalid geometries are listed in
`front_geometry.json`. No interior-density filter is applied; the near-voter
figure shows interior activity separately because a decaying interior can
invalidate interpretation of an envelope velocity as a stationary front.

`bulk_runs.csv`, `bulk_windows.csv`, `front_runs.csv`, and `front_records.csv`
retain completed individual-run diagnostics with an `ensemble_complete` flag.
The `bulk_summary`, `bulk_window_summary`, `bulk_curves`, `front_points`, and
`front_crossings` outputs contain only eligible complete ensembles. Common
absolute-time bulk windows distinguish size comparisons from changes in the
observation horizon. All final groups use four seeds except the L=4096
orientation groups, which use eight.

The figure normalization reference is centralized in
`results/coexistence_reference.json`. Updating it changes labels, never the
actual q values or the measured front velocities. Re-run the wetting analysis
with the same qstar before regenerating its figure. The two-droplet overlay now
interpolates the corrected seed-mean planar velocities at each actual q;
`results/droplets/two_drop_rates.csv` records the curvature comparison.

## Near-voter magnetic moments

`python3 -m campaign.v2.analysis.susceptibility` reads all completed E6/E7
trajectories, uses actual header probabilities, integrates each run from t=2e5
to its planned end, and then weights runs equally. An absorbing run has m²=1
and is retained. A common [2e5,5e5] window tests unequal legacy durations.
The explicitly labelled active-conditional diagnostic retains the old density
cutoff and reports retained time/run fractions; it is not the primary estimate.
Run-bootstrap intervals preserve the m²/m⁴ covariance in the Binder ratio.
Outputs: `campaign/v2/results/susceptibility/{runs,summary,ratios}.csv` and
`analysis.json` with hashes. `paper/v2/figs/fig_voter.py` reads the CSV directly.
Tests: `python3 -m unittest campaign.v2.analysis.test_susceptibility`.

## Completed bundle import after the scratch quota failure

The final scratch snapshot was frozen at 07:05:35 UTC on 2026-11-03. Do not
copy that source over resumed results. Runner v3 writes on node-local storage
and publishes atomic task archives under `/home/mo2273/dn_v2_stage/gpu/bundles/`.
Read-only synchronization and local import are:

```sh
rsync -a --exclude='*.tmp.*' torch:/home/mo2273/dn_v2_stage/gpu/bundles/ data/v2/gpu/bundles/
python3 -m campaign.v2.analysis.import_gpu_bundles
```

Only the completed `bundles` tree is eligible. The separate `checkpoints` tree
is never analysis input. The importer checks flat regular-file archive members,
manifest parameters, text/profile completion, and final-state size, time and
SHA256. Each task is staged before replacing its local directory; prior local
files remain in `data/v2/gpu/previous_runs/`. A conflicting completed final-state
digest is rejected. `bundle_imports.json` records every imported archive hash.
New wetting tasks are imported by the same helper but reduced separately by
`gpu_wetting.py`.

The legacy critical-droplet figure uses completed fixed-replica boundary-growth
counts, excludes the unfinished R90 batch, and reads corrected outer-envelope
E4/E12 strip velocities from `results/drop_fronts`. `droplets.py` interpolates
the 50% radius linearly in logit probability versus log radius; unbracketed or
ambiguous crossings are never extrapolated. `critical_radius.json` records
conditional binomial-bootstrap radius intervals, the joint replica/front-run
bootstrap for the curvature fit, and tracker/window sensitivity. These intervals
do not include the finite horizon, radius interpolation, or curvature-model
systematics. Recompute the profile reduction with:

```sh
python3 -m campaign.v2.analysis.fronts data/runs/E4_voter/out/p0.5_*.txt data/runs/E12_phase/out/p0.5_*.txt --out campaign/v2/results/drop_fronts --tmins 1000,2000,3000,10000 --thresholds 0.2,0.3,0.5
python3 paper/v2/figs/fig_drop.py
```

## Direct spatial correlations

`python3 -m campaign.v2.analysis.spatial_correlations` measures axis-averaged
spin correlations in the final configurations of complete four-seed flagship
bulk cohorts with L>=8192. The primary quantity is the raw spin product,
`s=2*sigma-1`: exact ensemble Z2 symmetry fixes the mean spin to zero. The
per-configuration mean-subtracted quantity is retained only as a diagnostic.
Absorbed configurations remain in the ensemble and contribute C(r)=1.
Uncertainty is a pointwise bootstrap of independent runs, not sites.

The packed-bit calculation is tested against independent site-array shifts,
including word-boundary shifts and ordered states. Every real grid must match
its completion SHA256/time and reproduce its recorded alive fraction and
nearest-neighbor bond density. These are equal-time, finite-horizon data;
one final snapshot per run alone does not establish time stationarity, an
asymptotic exponent, or a unique correlation length.

`fetch_frozen_grids.py` fetches only missing final.bin files for already-complete
flagship tasks from the frozen scratch source, never overwriting an existing
file, then checks the pre-existing local completion digest. It does not copy
old time series or profiles over resumed data. Its audit is in
`reports/v2/hpc/frozen_grid_fetch.json`.

`paper/v2/figs/fig_spatial_v2.py` waits for the declared large-size comparison
cohorts before plotting. Signed correlations and approximate pointwise 95%
run-bootstrap intervals are shown, including the long-distance fluctuations.

Bulk summaries retain `half_recording_lower_mean` and
`half_recording_upper_mean` (and the corresponding absorption fields). These
are mean observation brackets associated with the first recorded crossing; a run
right-censored at the common planned T contributes T to both bounds. A
half-density event already present at the t=1000 origin contributes [0,1000].
The primary means and bootstrap intervals use upper observation times.
Recording-resolution brackets and run-sampling intervals are distinct. Activity
is not monotone: an unobserved earlier dip can recover before the next record.
The preceding/current-record bracket therefore does not rigorously bound its
first physical threshold passage. The primary half-density event is explicitly
the first recorded crossing on the declared output grid. Consensus is an
absorbing event and does not have this ambiguity.

The local maintenance command `python3 -m campaign.v2.analysis.gpu_maintenance`
performs only read-only home-bundle synchronization and sequential local
import/reduction/figure updates every five minutes, until its explicit UTC
deadline. `data/v2/gpu/maintenance.jsonl` records each cycle. It never accesses
the frozen scratch source or submits jobs. `gpu_report.py` writes a compact
complete-ensemble numerical report; `results/gpu/figure_coverage.json` records
which declared parameter groups are still missing from each figure. Figures
are withheld until every panel has eligible data.

The separate `production_manifest_orientation5m.json` predeclares 36 independent
long-time front runs (nine four-seed groups). Their new initial and dynamical
seeds are disjoint from all earlier manifests. Both `fronts.summarize` and the
GPU reduction retain T in their grouping keys; the original orientation figure
explicitly selects T=2 million and `fig_gpu_fronts5m` selects T=5 million. The new figure compares
2--5m velocities, successive-window crossings, and paired 3.5--5 million minus 2--3.5 million
velocities. `front_stationarity.py` requires all four runs valid in both windows;
its paired crossing contrast also requires all three predeclared q groups and
resamples the same run indices across the two windows at each q. No early-window
observations are pooled across the original and new cohorts.

Before any 5-million-step worker began, its common q grid was amended to
.04593/.04594/.04595 to keep the normal(2,1) bands separated through 5 million steps. The
previous .04592/.04594/.04596 manifest is retained as
`production_manifest_orientation5m_prestart_v1.json`; the active manifest embeds
the old SHA256, fixed-grid amendment, geometric rationale and verification that all
four scheduled workers were still pending. The original 2-million-step grid was unchanged.

The schedule audit identified one frozen claimed completion,
`orientation_L16384x16384_ab1-1_q.04596_s2`, missing its t=1 million text record and t=950,000/1 million profile records.
It was already excluded for timeline mismatch.
Exact-seed repair job 21730749 uses the unchanged task ID, initial seed and
dynamical seed through an empty legacy source; frozen originals remain intact.
The importer checks the preserved evidence, all overlapping numerical records,
all preserved profile counts, and the original final-state digest before the
replacement can enter a cohort. `repair_verification.json` records the replay
result. The repaired trajectory exactly reproduces all 51 preserved text records,
all 50 preserved profile records, and the original final-state SHA256; it restores
all 52 scheduled records. Missing observations are regenerated by simulation,
never interpolated.

`python3 -m campaign.v2.analysis.legacy_reference_ratios` recomputes a historical
cache comparison after the shared p=.5 reference changes; it is not a paper
input and its output is explicitly marked superseded. It preserves
the absolute E18 crossing from the cached three-seed velocity table and records
that estimator's adaptive-window limitation. The old quoted uncertainty is
rescaled only; it is not reinterpreted as a newly validated confidence interval.

The subsequent `python3 -m campaign.v2.analysis.lazy05_audit` reanalyzes all
24 E18 raw profiles with the corrected envelope tracker, three declared seeds
per q, and fixed windows requiring valid geometry in every seed. Input time is
in update steps; physical time is 0.25 times that value, and physical velocity
is the per-step velocity divided by 0.25. Its results under `results/lazy05_fronts`
supersede a simple renormalization of the cached legacy tracker when reporting
the corrected absolute crossing. Late-window, threshold, and interior-density
diagnostics remain separate from the conditional run-bootstrap interval.

`production_manifest_near_midpoint.json` independently declares four new
L=8192, p=.995, q=.000642, T=20-million bulk trajectories. This parameter was
chosen between the earlier .000640 late-collapse and .000644 persistent
trajectories before any midpoint result existed. All four new seeds must
complete; absorptions are retained and no survivors from other q values are
pooled into the cohort. Worker submissions 21731854--21731857 preserve the
hard 11 UTC deadline and do not alter existing job priorities.

The wetting CLI defaults read the shared coexistence reference and systematic
half-range directly, preventing a later rerun from silently restoring the old
normalization. Explicit CLI values remain available for declared sensitivity
checks, and the reference plus analysis sources are included in provenance.

For delivery, maintenance stops when every declared primary/wall run is
imported, or at its explicit final local deadline. After the coordinated final
reference-dependent refresh, run `python3 -m campaign.v2.analysis.freeze_numerics`
to record numerical status, any withheld tasks, the chosen reference, and
SHA256 hashes of analysis outputs, figures, sources and import/repair provenance.
It refuses to run while the maintenance process is active and never regenerates
figures. Leave the recorded files unchanged during final compilation/packaging.

`mean_summary` enumerates all 256 ordered empirical-bootstrap resamples for
four-run means and paired means. This removes Monte Carlo variation between
identical four-run observables; it does not provide exact confidence coverage.
Other cohort sizes retain Monte Carlo sampling, and joint front-root and KM
bootstrap algorithms are unchanged. The existing random-stream advancement is
preserved so this refinement does not perturb unrelated sampled intervals.

`python3 -m campaign.v2.analysis.front5m_diagnostics` exports per-run late
interior ratios and all three threshold geometry/hole diagnostics for every
complete fixed 5M cohort, alongside the paired velocity contrasts. It writes
`reports/v2/front5m_diagnostics.json` without changing scientific estimates.

The final bounded `production_manifest_axial10m.json` independently declares
12 axial L8192 runs, four new seeds at each q=.04593/.04594/.04595, through
10 million steps. All 24 initial/dynamical seeds and their low32-bit keys are
disjoint from previous declarations. This check was declared after the5M study
resolved axial temporal drift; it does not reuse those trajectories. Its new
jobs alone have an11:05UTC ceiling; previous jobs keep their original11UTC
ceiling and priority. The local maintenance backstop is11:10UTC solely for
final completed-bundle import/reduction.

`axial10m_diagnostics.py` compares2–5M with5–10M within each new trajectory,
using all256 empirical bootstrap resamples for four-run paired mean changes.
Paired crossing differences require all three declared q groups, the same four
runs valid throughout both windows at each q, and unique in-grid roots. The
bootstrap reuses indices across the two windows within q and remains conditional
on both resampled roots existing; invalid or unbracketed cases are explicit.
Original T2M/T5M estimates and figures remain separate.

`python3 paper/v2/figs/fig_gpu_near20m.py` renders the largest near-voter
time-series and paired late-density-drift comparison only after all four declared
L8192, T20M cohorts are complete. Its pointwise bands enumerate all 256 ordered
four-run empirical bootstrap resamples. The existing four-panel near-voter figure
retains its original scope. A standalone caption is in
`paper/v2/gpu_near20m_caption.tex`.
