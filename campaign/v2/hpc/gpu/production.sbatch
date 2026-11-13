#!/usr/bin/env bash
#SBATCH --job-name=dn_v2_gpu_prod
#SBATCH --account=torch_pr_1342_courant
#SBATCH --partition=h200_courant
#SBATCH --gres=gpu:h200:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=05:30:00
#SBATCH --signal=B:USR1@90
#SBATCH --output=/scratch/mo2273/dn/v2/logs/gpu_prod_%j.log
set -euo pipefail
FAMILY="${1:?Pass orientation, flagship_bulk, or near_voter}"
TASKROOT="${TASKROOT:-/scratch/mo2273/dn/v2}"
IMG="${IMG:-/share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif}"
GPU_SOURCE="$TASKROOT/campaign/v2/hpc/gpu"
export BUILD_DIR="$TASKROOT/gpu/builds/$SLURM_JOB_ID"
export CUDA_ARCH=90 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
mkdir -p "$BUILD_DIR" "$TASKROOT/gpu/workers/$SLURM_JOB_ID"
nvidia-smi --query-gpu=name,uuid,driver_version --format=csv
apptainer exec --nv --bind /scratch,/home "$IMG" bash "$GPU_SOURCE/build.sh"
apptainer exec --nv --bind /scratch,/home "$IMG" python3 "$GPU_SOURCE/validate.py" \
  --cpu /scratch/mo2273/dn/bin2/dnsim_v15 --gpu "$BUILD_DIR/dnsim_gpu" \
  --work "$TASKROOT/gpu/workers/$SLURM_JOB_ID/validation"
printf 'PRODUCTION_VALIDATION_PASSED\n'
exec apptainer exec --nv --bind /scratch,/home "$IMG" python3 "${DNSIM_RUNNER:-$GPU_SOURCE/production_runner_v2.py}" \
  --family "$FAMILY" --binary "$BUILD_DIR/dnsim_gpu" \
  --manifest "${DNSIM_MANIFEST:-$GPU_SOURCE/production_manifest.json}" --root "$TASKROOT/gpu"
