#!/bin/bash
#SBATCH -J r2_bands
#SBATCH -A torch_pr_1342_general
#SBATCH -p cs
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o /scratch/mo2273/dn/logs2/%x_%j.log
# Initial bands for R2_orient / R2_orient2k: tilt.py make, numpy seed 100000 + L + 100 (3a + b) + s.
IMG=/share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif
D=/scratch/mo2273/dn/runs/R2_orient/init
mkdir -p $D
apptainer exec --bind /home,/scratch "$IMG" bash -c '
for L in 4096 2048; do for ab in "1 0" "2 1" "1 1"; do set -- $ab; for s in 1 2 3 4 5 6 7 8; do
  f='"$D"'/band_L${L}_ab$1-$2_s$s.ckpt
  [ -s $f ] && continue
  python3 /scratch/mo2273/dn/campaign2/tilt.py make $f --L $L --a $1 --b $2 --seed $((100000 + L + 100 * (3 * $1 + $2) + s)) && echo "made $(basename $f)"
done; done; done
ls -la '"$D"' | wc -l; echo BANDS_DONE'
