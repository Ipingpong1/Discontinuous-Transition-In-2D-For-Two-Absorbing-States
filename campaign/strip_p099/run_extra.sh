#!/bin/sh
# Second chain (started 10-31 00:40 while run_all.sh was running): runs_extra.txt with P extra slots
# (2, so 14 jobs on the 16 logical CPUs), same run_one.sh; safe to re-run (finished runs skipped).
D=/home/mo2273/Desktop/dn/campaign/strip_p099
P=${P:-2}
echo "EXTRA_START P=$P $(date '+%F %T')" >> "$D/run_all.log"
cat "$D/runs_extra.txt" | while IFS= read -r line; do
    out=$(printf '%s\n' "$line" | sed -n 's/.*--out \([^ ]*\.txt\).*/\1/p')
    if [ -f "$out" ] && grep -q '^# done' "$out"; then echo "SKIP $(basename "$out")" >> "$D/run_all.log"; continue; fi
    printf '%s\n' "$line"
done | xargs -P "$P" -I{} "$D/run_one.sh" {}
echo "EXTRA_DONE $(date '+%F %T')" >> "$D/run_all.log"
