#!/bin/sh
# run_one.sh "<dnsim args>": one dnsim job at nice 10, single thread; stdout/stderr to <out>.log,
# end line (time, rc, the run's "# done" line) appended to run_all.log.  A lock <out>.running
# (holding the PID) is kept while the job runs, so that a second chain reaching the same line
# skips it instead of launching a second dnsim on the same checkpoint (added 10-31 04:50).
D=/home/mo2273/Desktop/dn/campaign/strip_p099
export OMP_NUM_THREADS=1
out=$(printf '%s\n' "$1" | sed -n 's/.*--out \([^ ]*\.txt\).*/\1/p')
lock="${out%.txt}.running"
if [ -e "$lock" ]; then echo "SKIP running $(basename "$out") (lock pid $(cat "$lock"))" >> "$D/run_all.log"; exit 0; fi
echo $$ > "$lock"
# shellcheck disable=SC2086  # $1 is the whitespace-separated dnsim option line
nice -n 10 /home/mo2273/Desktop/dn/dnsim $1 > "${out%.txt}.log" 2>&1
rc=$?
rm -f "$lock"
echo "end $(date '+%F %T') rc=$rc $(basename "$out") $(grep '^# done' "$out" 2>/dev/null | cut -c3-)" >> "$D/run_all.log"
