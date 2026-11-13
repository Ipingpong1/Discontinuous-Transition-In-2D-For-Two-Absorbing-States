#!/bin/sh
# Planar-front scan: strips of the active phase (random density 1/2, rows [0, Ly/2)) in the
# all-0 background, Molofsky rule at fixed p2, p1 over a list.  Usage:
#   sh scan.sh OUTDIR P2 "P1 P1 ..." "SEED SEED ..." T [LX LY EVERY]
# Runs 8 single-thread jobs at a time (8 physical cores); finished runs (with a "# done" line)
# are skipped, so a scan can be re-run to fill in.
D=$(cd "$(dirname "$0")" && pwd)
OUT=$1; P2=$2; P1S=$3; SEEDS=$4; T=$5; LX=${6:-1024}; LY=${7:-4096}; EVERY=${8:-2000}
mkdir -p "$OUT"
LIST=$OUT/runs_$(date +%H%M%S).txt
: > "$LIST"
for p1 in $P1S; do for s in $SEEDS; do
  f=$OUT/p2_${P2}_p1_${p1}_s${s}
  grep -q "^# done" "$f.txt" 2>/dev/null && continue
  echo "--molofsky $p1,$P2 --Lx $LX --Ly $LY --T $T --init strip:0.5:0.5 --bg 0 --stop-absorbed --seed $s --ppd 10 --every $EVERY --out $f.txt --profile $f.prof" >> "$LIST"
done; done
echo "$(wc -l < "$LIST") runs -> $OUT"
export OMP_NUM_THREADS=1
xargs -P 8 -I{} sh -c "$D/dnsim {} > /dev/null 2>&1" < "$LIST"
echo SCAN_DONE "$OUT"
