#!/bin/sh
# run_one.sh "<dnsim args>": one dnsim v1.5 job at nice 10, single thread; log to <out>.log, end line to run_local.log.
D=/home/mo2273/Desktop/dn/campaign/round2
export OMP_NUM_THREADS=1
out=$(printf '%s\n' "$1" | sed -n 's/.*--out \([^ ]*\.txt\).*/\1/p')
if [ -f "$out" ] && grep -q '^# done' "$out"; then echo "SKIP $(basename "$out")" >> "$D/run_local.log"; exit 0; fi
# shellcheck disable=SC2086
nice -n 10 $D/dnsim $1 > "${out%.txt}.log" 2>&1
rc=$?
echo "end $(date '+%F %T') rc=$rc $(basename "$out") $(grep '^# done' "$out" 2>/dev/null | cut -c3-)" >> "$D/run_local.log"
