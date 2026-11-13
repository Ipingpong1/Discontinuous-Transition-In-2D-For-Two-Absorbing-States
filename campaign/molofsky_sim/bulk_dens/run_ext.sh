#!/bin/sh
# Extension chain: waits for the main chain's DONE line, saves the '# done' lines of the main runs
# (their wall times are lost when a run is continued in place), then continues the runs at
# dp1 = 0.3, 0.4, 0.5e-4 (L 2048) to 1.2e6 steps and the L 4096 runs at 0.5e-4 and 1e-4 to 6e5
# through --resume; 12 niced single-thread jobs.
D=$(cd "$(dirname "$0")/.." && pwd); cd "$D"
export OMP_NUM_THREADS=1
P=${P:-12}
until grep -q '^DONE' bulk_dens/run_dens.log; do sleep 60; done
grep -H '^# done' bulk_dens/p2_*.txt > bulk_dens/walltimes_main.txt
python3 bulk_dens/make_runs.py --ext >> bulk_dens/run_dens.log
echo "EXT_START $(date) P=$P $(wc -l < bulk_dens/runs_ext.txt) runs" >> bulk_dens/run_dens.log
xargs -P "$P" -I{} sh -c "nice -n 10 $D/dnsim {} > /dev/null 2>&1" < bulk_dens/runs_ext.txt
echo "EXT_DONE $(date)" >> bulk_dens/run_dens.log
