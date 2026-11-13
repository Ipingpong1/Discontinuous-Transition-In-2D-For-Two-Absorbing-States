#!/usr/bin/env bash
#SBATCH --job-name=dn-v2-tails-recovery
#SBATCH --account=torch_pr_1342_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=05:00:00
#SBATCH --output=/home/mo2273/dn_v2_stage/cpu/logs/tails_recovery_%j.slurm
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNSIM_CPU_WORKERS=3 DNSIM_EXCLUDE_SHARDS=03
exec apptainer exec /share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif python3 /scratch/mo2273/dn/v2/cpu/tails_recovery_runner.py
