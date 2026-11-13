#!/usr/bin/env bash
#SBATCH --job-name=dn-v2-legacy-replay
#SBATCH --account=torch_pr_1342_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/home/mo2273/dn_v2_stage/cpu/legacy_replay/slurm_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec apptainer exec /share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif python3 /home/mo2273/dn_v2_stage/cpu/legacy_replay/runner.py
