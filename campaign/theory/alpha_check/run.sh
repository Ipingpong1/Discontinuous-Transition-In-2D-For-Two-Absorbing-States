#!/bin/bash
# kappa* ∝ alpha check at p = 0.99: strip fronts at alpha 0.25 and 1, q around 2*alpha*q*(0.5)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=6
for a in 0.25 1; do
  for f in 0.8 0.9 1.0 1.1 1.25; do
    q=$(python3 -c "print(round(0.00127*2*$a*$f, 7))")
    n=p0.99_a${a}_j2_q${q}_Lx1024_Ly4096_s1
    [ -s $n.txt ] && grep -q '^60000 ' $n.txt && continue
    ../../../dnsim --alpha $a --voter 0.99 --jflip 2:$q --Lx 1024 --Ly 4096 --T 60000 \
      --init strip:0.5:0.5 --bg 0 --every 1000 --ppd 10 --stop-absorbed --seed 1 \
      --out $n.txt --profile $n.prof > $n.log 2>&1
    echo "$n done $(date +%T)"
  done
done
