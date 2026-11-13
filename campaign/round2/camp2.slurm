#!/bin/bash
#SBATCH -A torch_pr_1342_general
#SBATCH -p cs
#SBATCH -c 32
#SBATCH --mem=16G
#SBATCH -t 10:00:00
#SBATCH -o /scratch/mo2273/dn/logs2/%x_%A_%a.log
# Round 2 (2026-10-31): one shard of an experiment per array task, 32 single-threaded dnsim v1.5 runs
# at a time (campaign/run_shard.sh; finished runs are skipped, unfinished ones resume).
#   sbatch -J r2_$EXP --array=0-$((NSHARDS-1)) --export=ALL,EXP=R2_orient camp2.slurm
# Everything lives under /scratch (home is over its file quota).
IMG=/share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif
RUNROOT=/scratch/mo2273/dn
SHARD=$(printf "%s/runs/%s/tasks/shard_%03d.txt" "$RUNROOT" "$EXP" "$SLURM_ARRAY_TASK_ID")
[ -f "$SHARD" ] || { echo "no shard $SHARD"; exit 1; }
export DNSIM=/scratch/mo2273/dn/bin2/dnsim_v15
apptainer exec --bind /home,/scratch "$IMG" /scratch/mo2273/dn/campaign2/run_shard.sh "$SHARD"
