#!/bin/sh
# Third chain (10-31 02:10): runs_fill.txt, launched one run at a time into the slots the other chains
# free, keeping at most CAP dnsim processes on the machine (14 = the 12 + 2 of the first two chains);
# finished runs are skipped.
D=/home/mo2273/Desktop/dn/campaign/strip_p099
CAP=${CAP:-14}
echo "FILL_START CAP=$CAP $(date '+%F %T')" >> "$D/run_all.log"
while IFS= read -r line; do
    out=$(printf '%s\n' "$line" | sed -n 's/.*--out \([^ ]*\.txt\).*/\1/p')
    if [ -f "$out" ] && grep -q '^# done' "$out"; then echo "SKIP $(basename "$out")" >> "$D/run_all.log"; continue; fi
    while [ "$(pgrep -c -x dnsim)" -ge "$CAP" ]; do sleep 60; done
    "$D/run_one.sh" "$line" &
    sleep 30
done < "$D/runs_fill.txt"
wait
echo "FILL_DONE $(date '+%F %T')" >> "$D/run_all.log"
