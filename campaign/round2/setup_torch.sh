#!/bin/bash
# Run on Torch after the staging tree has been unpacked into /scratch/mo2273/dn: out dirs, the E16
# Lx 4096 runs to be extended (copied with their checkpoints; their original seeds restored in the
# shard lines), then the jobs.
set -u
R=/scratch/mo2273/dn
for e in R2_orient R2_orient2k R2_drop2 R2_dropbulk R2_wall; do mkdir -p $R/runs/$e/out; done
for x in 0.044982 0.045441 0.0456705; do for s in 1 2; do
  n=p0.5_a0.5_j2_q${x}_Lx4096_Ly1024_s$s
  for ext in txt prof ckpt; do [ -f $R/runs/R2_wall/out/$n.$ext ] || cp $R/runs/E16_wall/out/$n.$ext $R/runs/R2_wall/out/$n.$ext; done
  seed=$(grep -o -- "--seed [0-9]*" $R/runs/E16_wall/out/$n.txt | head -1 | awk '{print $2}')
  [ -z "$seed" ] && seed=$(grep -m1 "^# Lx" $R/runs/E16_wall/out/$n.txt | awk '{print $8}')
  sed -i "s|\(R2_wall/out/$n\t.*--seed \)[0-9]*|\1$seed|" $R/runs/R2_wall/tasks/shard_*.txt $R/runs/R2_wall/manifest.tsv
  echo "extend $n from t=$(grep '^# done' $R/runs/R2_wall/out/$n.txt | awk '{print $4}') seed $seed"
done; done
grep -c "" $R/runs/R2_wall/tasks/shard_*.txt
cp /home/mo2273/dn/campaign/run_shard.sh /home/mo2273/dn/campaign/tilt.py $R/campaign2/
cp $R/campaign2/camp2.slurm $R/campaign2/bands.slurm $R/campaign2/ 2>/dev/null
cd $R/campaign2 || exit 1
BJ=$(sbatch --parsable bands.slurm); echo "bands job $BJ"
sbatch -J r2_R2_orient   --dependency=afterok:$BJ --array=0-11 --export=ALL,EXP=R2_orient   camp2.slurm
sbatch -J r2_R2_orient2k --dependency=afterok:$BJ --array=0-3  --export=ALL,EXP=R2_orient2k camp2.slurm
sbatch -J r2_R2_drop2    --array=0-0 --export=ALL,EXP=R2_drop2    camp2.slurm
sbatch -J r2_R2_dropbulk --array=0-0 --export=ALL,EXP=R2_dropbulk camp2.slurm
sbatch -J r2_R2_wall     --array=0-1 --export=ALL,EXP=R2_wall     camp2.slurm
squeue -u mo2273 -o "%.10i %.18j %.4t %.10M %.6D %R" | head -30
