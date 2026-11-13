#!/bin/sh
D=$(cd "$(dirname "$0")" && pwd)
cd $D
sh scan.sh $D/scan 0.45 "0.181 0.182 0.183 0.184 0.185 0.186 0.187 0.188 0.189 0.190" "1 2 3" 60000
sh scan.sh $D/scan 0.50 "0.166 0.167 0.168 0.169 0.170 0.171 0.172 0.173 0.174 0.175" "1 2 3" 60000
export OMP_NUM_THREADS=1
xargs -P 8 -I{} sh -c "$D/dnsim {} > /dev/null 2>&1" < reps/runs.txt; echo REPS_DONE
xargs -P 8 -I{} sh -c "$D/dnsim {} > /dev/null 2>&1" < bulk/runs.txt; echo BULK_DONE
sh scan.sh $D/scan 0.45 "0.176 0.177 0.178 0.179 0.180 0.191 0.192 0.193 0.194 0.195 0.196" "1 2" 60000
sh scan.sh $D/scan 0.50 "0.160 0.161 0.162 0.163 0.164 0.165 0.176 0.177 0.178 0.179 0.180 0.181 0.182" "1 2" 60000
echo ALL_DONE
