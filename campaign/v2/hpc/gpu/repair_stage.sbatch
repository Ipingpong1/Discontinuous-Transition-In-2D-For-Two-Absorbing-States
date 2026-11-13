#!/usr/bin/env bash
#SBATCH --job-name=dn_v2_gpu_repair
#SBATCH --account=torch_pr_1342_courant
#SBATCH --partition=h200_courant
#SBATCH --gres=gpu:h200:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:35:00
#SBATCH --signal=B:USR1@90
#SBATCH --output=/home/mo2273/dn_v2_stage/logs/gpu_prod_%j.log
set -euo pipefail
FAMILY="orientation"
TASKROOT="/scratch/mo2273/dn/v2"
STAGEROOT="/home/mo2273/dn_v2_stage/gpu"
IMG="/share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif"
GPU_SOURCE="/home/mo2273/dn_v2_stage/campaign/v2/hpc/gpu"
DNSIM_LOCAL_ROOT="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/dn-v2-${SLURM_JOB_ID}-XXXXXX")"
trap 'rm -rf -- "$DNSIM_LOCAL_ROOT"' EXIT
# The frozen source is never changed or reused by this exact-seed repair.
mkdir -p "$DNSIM_LOCAL_ROOT/empty_legacy"
export DNSIM_LOCAL_ROOT BUILD_DIR="$DNSIM_LOCAL_ROOT/build" CUDA_ARCH=90
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
mkdir -p "$BUILD_DIR" "$STAGEROOT/workers/$SLURM_JOB_ID"
df -h "$DNSIM_LOCAL_ROOT" "$STAGEROOT"
nvidia-smi --query-gpu=name,uuid,driver_version --format=csv
apptainer exec --nv --bind /scratch,/home,"$DNSIM_LOCAL_ROOT" "$IMG" bash "$GPU_SOURCE/build.sh"
apptainer exec --nv --bind /scratch,/home,"$DNSIM_LOCAL_ROOT" "$IMG" python3 "$GPU_SOURCE/validate.py" \
  --cpu /scratch/mo2273/dn/bin2/dnsim_v15 --gpu "$BUILD_DIR/dnsim_gpu" \
  --work "$DNSIM_LOCAL_ROOT/validation"
tar -czf "$STAGEROOT/workers/$SLURM_JOB_ID/validation.tar.gz" -C "$DNSIM_LOCAL_ROOT" validation
printf 'PRODUCTION_VALIDATION_PASSED\n'
apptainer exec --nv --bind /scratch,/home,"$DNSIM_LOCAL_ROOT" "$IMG" python3 "$GPU_SOURCE/production_runner_v3.py" \
  --family "$FAMILY" --binary "$BUILD_DIR/dnsim_gpu" \
  --manifest "${DNSIM_MANIFEST:-$GPU_SOURCE/production_manifest_gpu_repair.json}" \
  --root "$STAGEROOT" --local-root "$DNSIM_LOCAL_ROOT" --resume-root "$DNSIM_LOCAL_ROOT/empty_legacy"
