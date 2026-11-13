#!/bin/sh
D=$(cd "$(dirname "$0")" && pwd)
cd $D
while [ $(grep -c SCAN_DONE $D/run_all.log) -lt 2 ]; do sleep 10; done
sh scan.sh $D/scan 0.45 "0.18525 0.1855 0.18575" "1 2 3" 60000
export OMP_NUM_THREADS=1
xargs -P 8 -I{} sh -c "$D/dnsim {} > /dev/null 2>&1" < reps/runs.txt; echo REPS_DONE
xargs -P 8 -I{} sh -c "$D/dnsim {} > /dev/null 2>&1" < bulk/runs.txt; echo BULK_DONE
echo CHAIN2_DONE
