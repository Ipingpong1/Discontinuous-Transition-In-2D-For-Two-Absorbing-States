#!/bin/sh
# Extra runs at the half-step points: strips at p1 = 0.18540 (p2 0.45) and 0.17062 (p2 0.50),
# bulk/lifetime runs there, and bulk runs 2.5e-5 above the line (0.185375).  Started with P=4 on
# the idle cores while run_fine.sh was in its last phases (raised with SIGUSR1 later).
D=$(cd "$(dirname "$0")" && pwd); cd "$D"
export OMP_NUM_THREADS=1
P=${P:-4}
echo "EXTRA_START $(date)" >> run_fine.log
xargs -P "$P" -I{} sh -c "nice -n 10 $D/dnsim {} > /dev/null 2>&1" < runs_fine_extra.txt
echo "EXTRA_DONE $(date)" >> run_fine.log
