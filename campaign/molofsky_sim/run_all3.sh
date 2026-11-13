#!/bin/sh
D=$(cd "$(dirname "$0")" && pwd)
cd $D
while ! grep -q CHAIN2_DONE $D/run_all2.log; do sleep 15; done
export OMP_NUM_THREADS=1
xargs -P 8 -I{} sh -c "$D/dnsim {} > /dev/null 2>&1" < drop/runs.txt; echo DROP_DONE
xargs -P 8 -I{} sh -c "$D/dnsim {} > /dev/null 2>&1" < bulk_L/runs.txt; echo BULKL_DONE
sh scan.sh $D/scan 0.45 "0.185 0.18525 0.1855" "4 5 6 7 8" 60000
sh scan.sh $D/scan 0.50 "0.17025 0.1705 0.17075" "4 5 6 7 8" 60000
sh scan.sh $D/scan_wide 0.45 "0.18525 0.1855" "1 2 3" 60000 2048 4096
sh scan.sh $D/scan_wide 0.50 "0.1705 0.17075" "1 2 3" 60000 2048 4096
xargs -P 8 -I{} sh -c "$D/dnsim {} > /dev/null 2>&1" < bulk_long/runs.txt; echo BULKLONG_DONE
sh scan.sh $D/scan 0.45 "0.176 0.177 0.178 0.179 0.180 0.191 0.192 0.193 0.194 0.195 0.196" "1 2" 60000
sh scan.sh $D/scan 0.50 "0.160 0.161 0.162 0.163 0.164 0.165 0.176 0.177 0.178 0.179 0.180 0.181 0.182" "1 2" 60000
echo CHAIN3_DONE
