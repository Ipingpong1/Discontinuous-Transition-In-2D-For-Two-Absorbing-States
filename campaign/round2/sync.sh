#!/bin/bash
# Mirror the round-2 results from Torch: time series of every R2_* run into data/runs/<EXP>/out (the
# project convention: *.txt only), the R2_drop2 row profiles (needed by fig_drop2.py), and the Torch-side
# reductions (figdata/R2: orient CSVs and tables, drop2/wall2 tables, rough CSVs; figdata/R2_wall: the
# extract_prof.py wall series and profiles for fig_wall05.py).
set -u
D=/home/mo2273/Desktop/dn/data
for e in R2_orient R2_orient2k R2_drop2 R2_dropbulk R2_wall; do
  mkdir -p $D/runs/$e/out
  rsync -a --include='*.txt' --include='manifest.tsv' --exclude='*' torch:/scratch/mo2273/dn/runs/$e/out/ $D/runs/$e/out/
  rsync -a torch:/scratch/mo2273/dn/runs/$e/manifest.tsv $D/runs/$e/ 2>/dev/null
done
rsync -a --include='*.prof' --exclude='*' torch:/scratch/mo2273/dn/runs/R2_drop2/out/ $D/runs/R2_drop2/out/
mkdir -p $D/figdata/R2 $D/figdata/R2_wall
rsync -a torch:/scratch/mo2273/dn/figdata/R2/ $D/figdata/R2/
rsync -a torch:/scratch/mo2273/dn/figdata/R2_wall/ $D/figdata/R2_wall/ 2>/dev/null
echo "synced $(date '+%F %T'): $(ls $D/runs/R2_*/out/*.txt | wc -l) time series, $(ls $D/figdata/R2 | wc -l) reduction files"
