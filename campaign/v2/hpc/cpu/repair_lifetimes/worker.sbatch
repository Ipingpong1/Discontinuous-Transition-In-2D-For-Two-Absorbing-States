#!/usr/bin/env bash
#SBATCH --job-name=dn-v2-life-repair
#SBATCH --account=torch_pr_1342_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/home/mo2273/dn_v2_stage/cpu/repair_lifetimes/logs/repair_%j.slurm
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec apptainer exec --bind /home:/home --bind /scratch:/scratch /share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif python3 /home/mo2273/dn_v2_stage/cpu/repair_lifetimes/runner.py lifetimes
