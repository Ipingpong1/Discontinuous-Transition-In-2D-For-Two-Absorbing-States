# Overnight v2 campaign

The original manuscript, simulator and datasets are preserved. New simulations and
analysis are isolated under `campaign/v2`, `paper/v2`, and `data/v2`.

Resource policy: use small CPU jobs for sparse replicas and reductions. GPU
production must be preceded by CPU/GPU bitwise checks and a utilization benchmark;
target sustained GPU utilization above 80 percent. Do not use login nodes for
simulation. Record commands, seeds, source hashes, scheduler IDs and resource logs.

Scientific priorities:

1. Correct existing front tracking and quantify geometric/fit systematics.
2. Resolve finite-time drift near p=0.995 with longer bulk and strip runs.
3. Test the flagship coexistence point at higher precision and more orientations.
4. Improve intermediate-size lifetime distributions and seed-tail statistics.
5. Report closure cutoff sensitivity and its limitations separately from data.

Start: 2026-11-02 22:30 America/New_York. Target end: approximately 2026-11-03
06:30 America/New_York. Jobs should fit within that window where possible.
