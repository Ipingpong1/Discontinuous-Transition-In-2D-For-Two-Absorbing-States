#!/usr/bin/env bash
#SBATCH --job-name=dn_v2_gpu_check
#SBATCH --account=torch_pr_1342_courant
#SBATCH --partition=h200_courant
#SBATCH --gres=gpu:h200:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/mo2273/dn/v2/logs/gpu_check_%j.log
set -euo pipefail
TASKROOT=/scratch/mo2273/dn/v2
IMG=/share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif
GPU_SOURCE=$TASKROOT/campaign/v2/hpc/gpu
export OMP_NUM_THREADS=2
mkdir -p "$TASKROOT/validation"
nvidia-smi --query-gpu=name,uuid,driver_version --format=csv
apptainer exec --nv --bind /scratch,/home "$IMG" bash "$GPU_SOURCE/build.sh"
apptainer exec --nv --bind /scratch,/home "$IMG" python3 "$GPU_SOURCE/validate.py" \
  --cpu /scratch/mo2273/dn/bin2/dnsim_v15 --gpu "$GPU_SOURCE/build/dnsim_gpu" \
  --work "$TASKROOT/validation/gpu_${SLURM_JOB_ID}"
printf 'VALIDATION_PASSED\n'
nvidia-smi --query-gpu=timestamp,utilization.gpu,utilization.memory,memory.used,power.draw \
  --format=csv --loop=2 > "$TASKROOT/validation/util_${SLURM_JOB_ID}.csv" &
MONITOR_PID=$!
trap 'kill "$MONITOR_PID" 2>/dev/null || true' EXIT
for LATTICE in 4096 8192; do
  apptainer exec --nv --bind /scratch,/home "$IMG" "$GPU_SOURCE/build/dnsim_gpu" \
    --L "$LATTICE" --T 50000 --alpha 0.5 --voter 0.995 --jflip 2:0.00064 \
    --seed 2026110201 --ppd 2 --every 10000 \
    --out "$TASKROOT/validation/benchmark_${SLURM_JOB_ID}_L${LATTICE}.txt"
done
printf 'BENCHMARK_DONE\n'
