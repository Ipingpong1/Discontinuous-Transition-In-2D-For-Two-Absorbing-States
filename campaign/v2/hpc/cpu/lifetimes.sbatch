#!/bin/bash
#SBATCH --job-name=dn-v2-lifetimes
#SBATCH --account=torch_pr_1342_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/home/mo2273/dn_v2_stage/cpu/logs/lifetimes_%j.slurm
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
apptainer exec --bind /scratch:/scratch --bind /home:/home /share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif python3 /scratch/mo2273/dn/v2/cpu/runner.py lifetimes
