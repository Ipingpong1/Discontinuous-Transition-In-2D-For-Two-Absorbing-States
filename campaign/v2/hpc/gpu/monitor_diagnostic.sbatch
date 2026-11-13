#!/usr/bin/env bash
#SBATCH --job-name=dn_v2_monitor_check
#SBATCH --account=torch_pr_1342_courant
#SBATCH --partition=h200_courant
#SBATCH --gres=gpu:h200:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:08:00
#SBATCH --output=/scratch/mo2273/dn/v2/logs/monitor_check_%j.log
set -euo pipefail
TASKROOT=/scratch/mo2273/dn/v2
IMG=/share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif
DIAG=$TASKROOT/validation/monitor_${SLURM_JOB_ID}
BIN=$TASKROOT/gpu/builds/21718937/dnsim_gpu
mkdir -p "$DIAG"
export OMP_NUM_THREADS=2
nvidia-smi --query-gpu=timestamp,uuid,utilization.gpu,utilization.memory,power.draw --format=csv --loop=2 > "$DIAG/host.csv" &
HOST_MONITOR=$!
apptainer exec --nv "$IMG" nvidia-smi --query-gpu=timestamp,uuid,utilization.gpu,utilization.memory,power.draw --format=csv --loop=2 > "$DIAG/container.csv" &
CONTAINER_MONITOR=$!
trap 'kill "$HOST_MONITOR" "$CONTAINER_MONITOR" 2>/dev/null || true' EXIT
apptainer exec --nv --bind /scratch,/home "$IMG" "$BIN" --L 8192 --T 500000 --seed 202611030701 --alpha .5 --voter .5 --jflip 2:.04592 --ppd 0 --every 50000 --out "$DIAG/bulk.txt"
apptainer exec --nv --bind /scratch,/home "$IMG" "$BIN" --L 8192 --T 200000 --seed 202611030702 --alpha .5 --voter .5 --jflip 2:.04592 --ppd 0 --every 5000 --pdir 1,1 --profile "$DIAG/profile.bin" --init strip:.5:.5 --out "$DIAG/profile.txt"
printf 'MONITOR_DIAGNOSTIC_DONE\n'
