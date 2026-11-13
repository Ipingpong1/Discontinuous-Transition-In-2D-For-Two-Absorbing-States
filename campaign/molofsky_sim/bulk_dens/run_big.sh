#!/bin/sh
# Two L 4096 runs at dp1 = 0.4e-4 (p1 = 0.18539) to 6e5 steps, started at 17:05 beside the main
# chain (14 jobs on 16 logical CPUs), nice 10, single thread.
D=$(cd "$(dirname "$0")/.." && pwd); cd "$D"
export OMP_NUM_THREADS=1
python3 bulk_dens/make_runs.py --big >> bulk_dens/run_dens.log
echo "BIG_START $(date)" >> bulk_dens/run_dens.log
xargs -P 2 -I{} sh -c "nice -n 10 $D/dnsim {} > /dev/null 2>&1" < bulk_dens/runs_big.txt
echo "BIG_DONE $(date)" >> bulk_dens/run_dens.log
