# Dense CUDA backend

This is an isolated backend for `campaign/round2/dnsim.c`; the original source is never edited. `generate.py` copies its Model declaration and exact integer Philox, Bernoulli-table, neighbor-count, and word-update functions, records the source SHA256, and generates a C host retaining the original initialization, probability quantization, measurements, profiles, snapshot/checkpoint formats, and CLI.

Each CUDA thread advances one 64-site word. Moore and von Neumann neighborhoods use the same extracted update code as the CPU. Two resident lattice buffers alternate synchronously. A 256-step CUDA graph, instantiated for both buffer parities, advances a device time counter; shorter remainders use ordinary kernel launches. Grid data and the last step's flip count return to the host only at output, snapshot, or checkpoint events. RNG counters remain exactly (time, row, word).

Build with CUDA 12.8 or newer:

```sh
NVCC=/usr/local/cuda/bin/nvcc bash campaign/v2/hpc/gpu/build.sh
```

The default target is sm80 plus compute80 PTX; set CUDA_ARCH for another architecture. Set BUILD_DIR or CPU_SOURCE to override paths. The binary is `build/dnsim_gpu`. Set OMP_NUM_THREADS explicitly for the retained CPU measurements. CUDA_VISIBLE_DEVICES selects the GPU.

Validation:

```sh
python3 campaign/v2/hpc/gpu/test_host.py
python3 campaign/v2/hpc/gpu/validate.py --cpu /path/to/dnsim_cpu --gpu /path/to/dnsim_gpu --work /path/to/results
```

The first command checks the generated event-driven host and extracted update helpers using CPU emulation. It does **not** establish CUDA correctness. The second must run on a GPU before scientific production. It compares all numerical time-series rows and byte-for-byte profiles, snapshots, final saves, and checkpoints across Moore, Molofsky/VN, arbitrary probability classes, deterministic rules, one-word rows, tilted profiles, graph parity boundaries, resume, and the largest valid RNG times. Headers and timing differ intentionally.

Limitations:
- `--stop-absorbed` checks only scheduled time-series outputs. Its stopping time can therefore be later than the CPU's. This is stated in every output file. Use the original CPU binary for exact lifetime studies.
- `--reps` is rejected; sparse spreading remains on CPU.
- `--selftest` is rejected to avoid presenting the original CPU tests as GPU validation.
- One model/context per process; the model lives in CUDA constant memory.
- Throughput excludes initialization and graph creation, like the original timer excludes initialization. CPU measurement and transfer time remain included.
- CUDA errors are fatal and printed to stderr. Existing completed checkpoints remain usable.

## Production and storage recovery

Real H200 validation passed before production and is repeated at the beginning
of every worker. Initial scientific throughput was approximately 5--7e11 site
updates/s for large lattices. The utilization monitor samples NVIDIA utilization
with jitter, includes checkpoint/profile phases, and stops a worker if a mature
120-second rolling average remains below 80% for 30 seconds. A fixed ten-second
monitor in an early worker aliased with periodic I/O; its replacement and job
history are recorded in `submitted.json`.

The account's shared scratch quota filled twice, around 06:10 and 06:51 UTC on
2026-11-03. Old scratch-writing workers are superseded by
`production_stage.sbatch` and `production_runner_v3.py`. Do not resubmit the old
batch scripts. The audited staging design uses:

- Code: `/home/mo2273/dn_v2_stage/campaign/v2/hpc/gpu`.
- Active simulation files: the job's node-local `SLURM_TMPDIR`.
- Atomic completed archives: `/home/mo2273/dn_v2_stage/gpu/bundles/FAMILY/TASK.tar.gz`.
- Resumable archives: the separate `checkpoints/FAMILY/TASK.tar.gz` tree.
- Worker provenance, validation archive and utilization: `workers/JOB/`.
- Frozen pre-recovery inputs: `/scratch/mo2273/dn/v2/gpu/runs/`.

`storage_stage.py` verifies manifest headers, scheduled text/profile records,
global counts and checkpoint dimensions/times before resume. A damaged prefix
is restarted with its original fixed seed; the original input is preserved.
Full completion additionally requires the final grid and terminal record to
agree. Archives are fsynced and atomically published before local files and
claims are removed. Checkpoint archives never enter scientific analysis.

The staged pilot, 21725768, passed CPU/GPU validation and independent inspection
of its first million-step checkpoint and complete two-million-step result.
`submit_stage.py` then submitted 18 further workers, recorded in
`reports/v2/hpc/staged_submission_20261103.json`. The fixed manifests and seeds
were unchanged. All workers have a hard 11:00 UTC scientific deadline.

Run the meaningful local staging checks with:

```sh
python3 -m unittest discover -s campaign/v2/hpc/gpu -p 'test_storage_stage.py' -v
```

For read-only result synchronization and manifest-aware import, see
`campaign/v2/analysis/README.md`. A new campaign on another date requires an
explicitly updated time budget and manifest deadline; the current scripts are
intentionally not perpetual launchers.
