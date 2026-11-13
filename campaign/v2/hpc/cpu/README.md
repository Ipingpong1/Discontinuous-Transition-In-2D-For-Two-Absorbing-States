# V2 CPU simulations on Torch

Submitted 2026-11-03 UTC, account `torch_pr_1342_general`, partition `cs`.

| Job | Work | Allocation |
|---|---|---|
| 21717999 | Lifetimes and near-line extension | 4 CPUs, 8 GiB, 6 hours |
| 21718000 | Fixed-replica seed tails | 4 CPUs, 8 GiB, 6 hours |

Both started at approximately 05:11:50 UTC. Each job runs at most four independent
single-threaded simulators. The combined allocation is eight CPUs. The runner
stops current work at **11:25 UTC** and does not start tasks with less than one
minute remaining; the outer deadline is 11:30 UTC. Checkpoints survive timeout.

Remote root: `/scratch/mo2273/dn/v2/cpu`.
The existing binary `/scratch/mo2273/dn/bin2/dnsim_v15` is read without modification.
Its SHA256 and every command are recorded in `manifest.json` and job events.

## Planned measurements

- Main lifetimes: p=alpha=0.5, q=0.046359, L=64,128,192,256, 32 seeds per size,
  T=10^7, random density-1/2 initial states, regular output every 1000 steps plus
  20 logarithmic outputs per decade, stop on observed absorption.
- Near-line extension after the main set: q=0.04605, L=192,256,384, 16 seeds per
  size, T=3x10^6, otherwise identical. These tasks may remain unstarted if the
  deadline is reached; the runner records this explicitly.
- Seed tails: p=0.99, alpha=0.5, q=0.0012954, L=8192, T=10^6, a fully live 2x2
  seed, eight independent shards of exactly 2,000,000 replicas each. Adaptive
  stopping is disabled (`--target 0`). Sixteen group sums are retained per shard.

Seeds are deterministic 64-bit values. Their high word comes from SHA256 of the
campaign/task identity; the low words are distinct, including between shards.
This avoids sharing the low RNG key word in the replica-key construction.

## Collection and interpretation

Fetch compact outputs without checkpoints:

    rsync -av --include='*/' --include='*.txt' --include='*.log' --include='*.jsonl' --exclude='*' torch:/scratch/mo2273/dn/v2/cpu/out/ data/v2/cpu/

The runner's event stream lives in remote `logs/*_events.jsonl`; Slurm logs are
`logs/*_<jobid>.slurm`. A `# done` record marks a completed simulation. A finite
run that reaches T with positive activity is right-censored, not an absorption
event. An absorption is first detected at an observation time; the previous
nonzero observation brackets the event, so time resolution must be retained in
the analysis. Do not replace censored lifetimes with T in an uncensored mean.

For seed survival, sum all group/shard numerators and denominators. Report the
number of actual survivors and escaped replicas separately. Unfinished shards retain a checkpointed replica count for progress reporting.
Primary inference waits for all eight fixed two-million-replica shards; no pooled
wall-clock-selected subset is treated as the final ensemble.

Checkpoint interval is 100,000 simulation steps for lifetimes and 1,000 replicas
for tails. Every task includes `--resume`; completed outputs are skipped. The
fixed deadline deliberately prevents blindly resubmitting these scripts on a
later date. A new explicitly budgeted campaign must update the deadline first.

No simulation was executed on a login node. Submission and provenance file
creation used the login node; all simulator processes run inside Slurm jobs.

## Quota recovery and audited statistics (2026-11-03)

Scratch-quota exhaustion stopped the original lifetime job 21717999 and seven
tail processes around 06:10 UTC. All checkpoints passed format/length validation.
Recovery jobs 21722717 (four lifetime workers) and 21723013 (three tail workers)
resumed the unchanged manifest. Atomic claims in `claims/tails_recovery` prevent
duplicate work; the surviving original shard 02 remained untouched.

At 06:48 UTC an additional four-worker tail job, 21724141, started shards 05–08.
The coordinator explicitly authorized up to twelve active single-threaded
simulations (eight tails and four lifetimes). The initial eight-process cap was
our conservative planning choice, not a user limit. Allocated CPU counts briefly
exceeded active simulations because old jobs retained idle allocated cores.

A second scratch-quota failure around 06:51 UTC stopped 21722717 and 21724141.
Original tail job 21718000 failed at 07:06 UTC. Job 21723013 retained only the
expensive live shard 03. At 07:12 UTC, all 71 stored CPU checkpoints were valid,
59 lifetime outputs had consistent completion footers, and all eight tail
checkpoints remained resumable. A complete local recovery snapshot was saved.

`relocate_cpu_home.py` can briefly stop only this job's simulator and runner,
verify that no scientific file is open, byte-verify a copy to home staging,
replace the old path by a symlink, and resume the same processes. The original
scratch tree is retained for already-open logs. Execution requires confirmed
home inode headroom from the coordinator. `compact_completed.py` packs only
validated completed-task checkpoints into a SHA-verified archive and removes
only their zero-byte logs; all raw scientific outputs remain separately readable.

`python3 -m campaign.v2.analysis.cpu_lifetimes` produces diagnostic progress
and final summaries only for wholly completed declared cohorts.
`paper/v2/figs/fig_lifetimes_v2.py` waits for all four 32-seed main cohorts.
The figure separates the first recorded half-density threshold crossing from
complete consensus and retains right censoring. Density is nonmonotone, so
observation brackets do not rigorously bound the first physical crossing.
Incomplete groups are not publication averages.

`python3 -m campaign.v2.analysis.seed_tails` keeps the old 3M and new fixed-N
samples separate. The original plan designated the 16M campaign as primary;
the predeclared 8M amendment below supersedes that designation. Each ensemble
requires all its fixed shards to finish. Selecting only shards that finish by a wall-clock deadline can select
against expensive surviving trajectories and is not accepted as final inference.
The original scientific seeds, binary, and 11:25 UTC deadline are unchanged.

At 07:17 UTC the compacted CPU tree was relocated to
`/home/mo2273/dn_v2_stage/cpu`; the old scratch path is now a symlink. All 117
files were byte-verified before the switch, and live shard 03 resumed in the
same process. New jobs 21725907 (four tail workers), 21725908 (three tail workers),
and 21725909 (four lifetime workers) use home Slurm logs. Tail recovery excludes
the still-running shard 03. The original scratch tree remains for old open
logs. Checkpoint compaction is repeated while the campaign runs to conserve
home inodes; the archive retains all completed checkpoint bytes and SHA256s.

## Independent confirmatory 8M amendment

At 07:38:23 UTC, before submission, `confirm8m/manifest.json` predeclared a
separate ensemble of sixteen 500000-replica shards, with adaptive stopping
disabled and new deterministic RNG low words disjoint from the original CPU
campaign and all mirrored E8 seed files. Physical parameters and the simulator
binary are unchanged. Jobs 21727007–21727010 each allocate four CPUs and
run four independent single-threaded shards; the total authorized active CPU
count is now 28, including the original eight tails and four lifetime workers.
All jobs retain the 11:25 UTC deadline.

The amendment responds to the long runtimes of rare trajectories in the
original two-million-replica shards. All original simulations continue unchanged.
The new 8M ensemble is the declared confirmatory analysis regardless of which
campaign finishes first. Its inference waits for all sixteen shards. The
original 16M ensemble is supplementary independent replication only if all
eight of its shards finish. No wall-clock-selected subsets are pooled.

New output root: `/home/mo2273/dn_v2_stage/cpu/confirm8m`; local scientific mirror:
`data/v2/cpu/confirm8m/tails`. The common analysis reports both ensembles
separately, and Fig. voter uses the new 8M only after its complete fixed ensemble
is available. Until then it retains the explicitly nominal legacy 3M curve.

After mirroring the four JSONL worker logs to
`reports/v2/hpc/confirm8m_logs`, run `python3 campaign/v2/hpc/cpu/audit_confirm8m.py`.
It verifies the frozen manifest hash, binary and runner hashes, exact task
commands and seeds, successful completion events, and the scientific reducer's
full header/count/schedule checks. Its report is
`reports/v2/hpc/confirm8m_completion_audit.json`. An incomplete ensemble is
explicitly ineligible for inference, and the audit exits with status 2.

The confirmatory ensemble completed at 10:32:46 UTC: all 8,000,000 planned
replicas, 19 survivors at T=1,000,000, and zero escapes. All four Slurm jobs
completed with exit status 0. The fixed-sample survival estimate is
2.375e-6, with exact 95% Clopper–Pearson interval [1.4299057e-6,3.7088542e-6].
The complete-cohort audit passed, and every one of 28,416 numerical output
rows matches its final checkpoint accumulator exactly. All sixteen final
checkpoints are preserved locally in
`data/v2/cpu/provenance/confirm8m_completed_checkpoints.tar.gz`, with a
separate verification receipt in `reports/v2/hpc/confirm8m_checkpoint_audit.json`.
The original 16M campaign is still separate; its incomplete state does not
affect the declared primary analysis.

## Repair of quota-induced observation gaps

A full regular-output-grid audit at 08:16 UTC found seven completed lifetime
files with 9–450 missing 1000-step records immediately before their resume
points. Valid final footers and checkpoints alone had not detected the missing
history. The affected tasks are L192 replicas12,13,15,16 and L256 replicas12,14,15
at q=0.046359. Final fixed-cohort inference now requires the complete regular
observation grid; the affected files remain excluded until repaired.

Job21729362, four additional CPUs, reruns exactly these seven original seeds
with the unchanged binary and all original physical parameters. Damaged files
are retained in a byte-verified archive. `repair_lifetimes/apply_repairs.py`
requires equality of every overlapping numeric record, final endpoint, and
final checkpoint byte, plus full regular coverage before canonical replacement.
The original main simulation job continues uninterrupted. Repair provenance
and both versions of the files are retained.

The seven repairs finished and were applied at approximately 08:43 UTC.
All 68,122 overlapping numerical records and all final checkpoint bytes matched
exactly; 1,280 missing records were restored. Both regular and logarithmic
schedules are now required. The canonical outputs have no rejected inputs.
The detailed mapping, hashes, and validation counts are in
`reports/v2/hpc/cpu_repairs_applied.json`; the original gapped files are also
preserved locally in `data/v2/cpu/provenance/gapped_originals.tar.gz`.

All 176 lifetime runs completed by 09:47:48 UTC: four main cohorts of 32
independent seeds and three near-line cohorts of 16. All seven complete cohorts
pass the full observation-grid audit. The 176 final checkpoints were copied to
`data/v2/cpu/provenance/lifetime_completed_checkpoints.tar.gz`; every archived
member was verified against its recorded SHA256. The main lifetime figure uses
only the predeclared 128-run campaign. Optional near-line results remain separate.
