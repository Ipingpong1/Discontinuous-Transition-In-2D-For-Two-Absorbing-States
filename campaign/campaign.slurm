#!/bin/bash
#SBATCH -J dn-camp
#SBATCH -A torch_pr_1342_general
#SBATCH -p cs
#SBATCH -c 16
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH -o /home/mo2273/dn/logs/%x_%A_%a.log
# One shard of an experiment per array task (16 single-threaded runs at a time).
#   sbatch --array=0-$((NSHARDS-1)) --export=ALL,EXP=E1_front campaign.slurm
# Resubmitting the same array resumes unfinished runs from their checkpoints.
# Run data lives under $RUNROOT (default $SCRATCH/dn): /home is near its file-count quota.
IMG=/share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif
RUNROOT=${RUNROOT:-$SCRATCH/dn}
SHARD=$(printf "%s/runs/%s/tasks/shard_%03d.txt" "$RUNROOT" "$EXP" "$SLURM_ARRAY_TASK_ID")
[ -f "$SHARD" ] || { echo "no shard $SHARD"; exit 1; }
apptainer exec --bind /home,/scratch "$IMG" "$HOME/dn/campaign/run_shard.sh" "$SHARD"
