# Revised manuscript and reproducibility guide

`main.tex` is the authoritative v2 manuscript. The original paper and simulator
are preserved. `main.pdf` is the compiled reading copy; `refs.bib` contains the
references actually used. The `*_revision.tex` files and
`campaign/v2/assemble_manuscript.py` document the initial assembly and are not
the current manuscript. The assembler refuses to overwrite v2.

The central claim is a quantitative numerical coexistence result at
`p=alpha=1/2`: finite bulk activity, front reversal in several orientations,
metastable decay, and complementary droplet-growth tests in a binary model with
exactly two absorbing configurations. The introduction places this result in
the specific two-dimensional history and acknowledges earlier reports of
discontinuous transitions. The voter endpoint and asymptotic wetting law remain
separate open questions.

`evidence_index.json` maps the main scientific claims to their analysis scripts,
result tables and interpretation limits. Paths in that index refer to the
project root, including material outside the compact manuscript-source ZIP.

## What changed

- Bulk estimates now integrate each trajectory in time and weight independent
  runs equally, retaining absorption and censoring. Half-density crossing and
  complete consensus are distinct observables.
- Fronts use an outer-envelope tracker with equal physical smoothing across
  orientations. Local crossings, fitting windows, geometry loss, bootstrap
  uncertainty and systematic sensitivity are reported separately.
- Droplet success means reaching the finite-box boundary before extinction and
  the declared horizon. The critical radius is a 50% growth probability, with
  conditional uncertainty, rather than a strict minimum surviving seed.
- Wetting fits use the same adsorption units and shared-baseline covariance.
  The new larger-wall measurements resolve size and time effects and test
  predictions from the earlier fits. Neither the legacy effective one-half
  slope nor a one-third law is established as the asymptotic exponent.
- Magnetic fluctuations retain absorbed runs. Uncertainty at the closest
  near-voter points prevents a blanket exclusion of an intervening critical
  regime. Adaptively stopped legacy seed curves are identified as descriptive.
- The coalescing-walk calculation distinguishes exact finite-cutoff identities
  from the unproved matching step in two dimensions. Its cutoff dependence and
  correlated threshold uncertainty are explicit.
- Dense-lattice simulations use an independently validated CUDA backend; the
  CPU and GPU integer updates and counter-based random numbers agree bitwise.
  New scientific summaries require every seed in the declared ensemble.

The referee report is `reports/v2/referee_review.pdf`; the additional independent
audit is `reports/v2/independent_final_audit.txt`. These review the incoming
claims and record the corrections, rather than certify an asymptotic theorem.

## Completed evidence in this revision

- 440 declared GPU trajectories completed: 412 bulk/front trajectories in 88
  fixed ensembles, and 28 wall trajectories in eight ensembles. The largest
  bulk systems contain 16384² sites; near-voter runs extend to 20 million
  steps, and the independent axial-front check extends to 10 million steps.
- 176 CPU lifetime runs completed, including a second ensemble closer to the
  transition. All runs, including consensus events, enter their declared summaries.
- The independent fixed eight-million-seed confirmation completed all sixteen
  shards, with 19 survivors and no boundary escapes at one million steps.
  The separately launched sixteen-million-seed campaign remains excluded from
  publication pooling: two of eight shards completed, with 13,295,943 replicas
  checkpointed overall at the deadline. Its files are preserved in the separate
  recovery archive and `original16m_final_snapshot.json` receipt.
- The operational reference is `q*=0.04594`, with normalization sensitivity
  `[0.04592, 0.04596]`. This is not a pooled confidence interval or an
  asymptotic bound: the independent long-front data resolve temporal shifts.

`reports/v2/numerical_final_snapshot.json` records the final numerical and figure
hashes. CPU publication hashes are in
`reports/v2/hpc/recovery/cpu_publication_freeze.json`. The Slurm accounting and
utilization summaries are in `reports/v2/hpc/`: 58.4 GPU allocation hours and
99.0 CPU-only allocated core hours. No campaign jobs remain active.

## Build

From this directory:

```sh
make TECTONIC=/home/mo2273/miniforge3/envs/tex/bin/tectonic
```

On another machine, `make TECTONIC=tectonic` uses an installed Tectonic. Standard
LaTeX with REVTeX 4.2 and BibTeX also works. The supplied figure PDFs are the
publication inputs. A PDF-only rebuild does not rerun simulations or analyses.

From the project root, the analysis workflow and figure commands are documented
in `campaign/v2/analysis/README.md` and `campaign/v2/theory/README.md`. Run-level
tables, settings and input SHA256s are retained in `campaign/v2/results/` and
`campaign/v2/theory/results/`. Original legacy inputs are under `data/runs/`;
new fixed-manifest inputs are under `data/v2/`.

## HPC provenance and recovery

All production simulations run through Slurm on NYU Torch. The GPU worker
records utilization and stops sustained averages below 80%; sparse replicas and
lifetimes use small CPU jobs. Manifests, seeds, source hashes, job IDs, validation
outputs and recovery records are under `campaign/v2/hpc/` and `reports/v2/hpc/`.
`campaign/v2/hpc/collect_accounting.py` records Slurm allocation hours, including
failed and restarted attempts. `summarize_utilization.py --sync` collects the
staged workers' raw monitor logs and reports their sampled utilization. Those
samples include compute and in-run I/O; the post-startup means exclude build,
validation, inter-task staging, and each run's initial 20 seconds, so they are
not allocation-wide efficiency estimates.
All GPU campaigns completed before their 3 November 2026 deadlines
(11:00 UTC for the original campaigns and 11:05 UTC for the final fixed
twelve-run axial-front extension). The extension has initial and dynamical
seeds independent of all preceding cohorts. The CPU deadline is 11:25 UTC.
Do not blindly resubmit these dated manifests.

Two account-wide scratch quota failures required recovery. Active GPU files now
live on node-local storage, with atomic completed/checkpoint bundles in
`/home/mo2273/dn_v2_stage/gpu`. CPU scientific outputs live in
`/home/mo2273/dn_v2_stage/cpu`, reached through the preserved scratch symlink.
Incomplete checkpoints are resumable records, not completed scientific data.
The final old scratch GPU mirror was frozen at 07:05:35 UTC; later imports use
only completed home bundles and preserve previous local files.

Nine closed legacy campaign directories were moved only after every file was
byte-SHA verified against `/archive/mo2273/dn-v2-storage-recovery-20261103/`.
Their original scratch paths are symlinks. **Torch archive storage is mounted
on login nodes, not compute nodes**: restore any of these directories to a
writable compute-visible location before compute-node reanalysis. The transfer
records list every file and hash. Additional precautionary archive copies did
not replace their original scratch directories. A similarly verified archive
preserves the pip download cache removed to free home inodes.

The review copy omits author names and affiliations; these remain to be supplied
before submission. No public repository DOI or author metadata has been invented.
