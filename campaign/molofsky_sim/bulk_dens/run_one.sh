#!/bin/sh
# The one run whose command the main chain received truncated (cut before --out, so its series went
# to /dev/null while its checkpoint was saved): 0.18555 seed 16, redone from scratch after deleting
# the orphan checkpoint.  Started 19:05 beside the other chains.
D=$(cd "$(dirname "$0")/.." && pwd); cd "$D"
export OMP_NUM_THREADS=1
L=$(grep 'p1_0.18555_L2048_s16' bulk_dens/runs_missing.txt)
echo "ONE_START $(date) $L" >> bulk_dens/run_dens.log
nice -n 10 $D/dnsim $L > /dev/null 2>&1
echo "ONE_DONE $(date)" >> bulk_dens/run_dens.log
