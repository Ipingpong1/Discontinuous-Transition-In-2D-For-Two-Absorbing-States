#!/bin/bash
# Reductions on Torch (inside the apptainer image, numpy): run after the R2_* jobs; writes
# /scratch/mo2273/dn/figdata/R2/*.txt|csv.  Usage: bash reduce.sh [orient|drop|wall|rough|all]
R=/scratch/mo2273/dn; C=$R/campaign2; F=$R/figdata/R2; mkdir -p $F
what=${1:-all}
run() { echo "== $1 $(date +%T)"; }
if [ $what = orient ] || [ $what = all ]; then
  run orient; python3 $C/orient2.py $R/runs/R2_orient --csv $F/orient_runs.csv > $F/orient_outer.txt 2>&1
  python3 $C/orient2.py $R/runs/R2_orient --est integ > $F/orient_integ.txt 2>&1
  python3 $C/orient2.py $R/runs/R2_orient --t0 20000 > $F/orient_outer_t20k.txt 2>&1
  python3 $C/orient2.py $R/runs/R2_orient --t0 50000 > $F/orient_outer_t50k.txt 2>&1
  python3 $C/orient2.py $R/runs/R2_orient2k --csv $F/orient2k_runs.csv > $F/orient2k_outer.txt 2>&1
fi
if [ $what = drop ] || [ $what = all ]; then
  run drop; python3 $C/drop2.py $R/runs/R2_drop2 > $F/drop2.txt 2>&1
fi
if [ $what = wall ] || [ $what = all ]; then
  run wall; python3 $C/wall2.py $R/runs/R2_wall --e16 $R/runs/E16_wall --csv $F/wall_runs.csv > $F/wall2.txt 2>&1
  python3 /home/mo2273/dn/paper_scripts/extract_prof.py $R/runs $R/figdata walls R2_wall > $F/extract_walls.txt 2>&1
fi
if [ $what = rough ] || [ $what = all ]; then
  run rough; ls $R/runs/R2_wall/out/*.snap | sed 's/\.snap$//' | xargs -P 8 -n 1 sh -c 'python3 '"$C"'/rough.py "$0" --tmin 100000 --csv '"$F"'/rough_$(basename "$0").csv --hist '"$F"'/whist_$(basename "$0").csv > '"$F"'/rough_$(basename "$0").txt 2>&1'
  ls $R/runs/R2_orient/out/*ab1-0*.snap | sed 's/\.snap$//' | xargs -P 8 -n 1 sh -c 'python3 '"$C"'/rough.py "$0" --tmin 40000 --csv '"$F"'/rough_$(basename "$0").csv --hist '"$F"'/whist_$(basename "$0").csv > '"$F"'/rough_$(basename "$0").txt 2>&1'
fi
echo "REDUCE_DONE $what $(date +%T)"
