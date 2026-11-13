#!/usr/bin/env bash
#SBATCH --job-name=dn-v2-confirm8m
#SBATCH --account=torch_pr_1342_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/home/mo2273/dn_v2_stage/cpu/confirm8m/logs/confirm8m_%j.slurm
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNSIM_CPU_WORKERS=4 DNSIM_EXCLUDE_SHARDS=
exec apptainer exec --bind /home:/home --bind /scratch:/scratch /share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif python3 /home/mo2273/dn_v2_stage/cpu/confirm8m/runner.py
