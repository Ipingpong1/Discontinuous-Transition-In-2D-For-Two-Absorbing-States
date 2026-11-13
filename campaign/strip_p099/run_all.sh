#!/bin/sh
# Chain for NOTES.md: runs_wide.txt (Lx 4096 strips), runs_inner.txt (strips at the three q nearest
# each crossing), runs_outer.txt (the outer q), runs_bulk.txt (bulk stationarity), as P single-thread
# jobs at nice 10 via run_one.sh (12 on the 8 cores + HT, 4 logical CPUs left free).  The lists are
# concatenated into one stream so that no slot idles between lists.  Finished runs are skipped
# (output has "# done"); unfinished ones resume from their checkpoint, so re-running is safe.
D=/home/mo2273/Desktop/dn/campaign/strip_p099
P=${P:-12}
echo "START P=$P $(date '+%F %T')" >> "$D/run_all.log"
cat "$D/runs_wide.txt" "$D/runs_inner.txt" "$D/runs_outer.txt" "$D/runs_bulk.txt" | while IFS= read -r line; do
    out=$(printf '%s\n' "$line" | sed -n 's/.*--out \([^ ]*\.txt\).*/\1/p')
    if [ -f "$out" ] && grep -q '^# done' "$out"; then echo "SKIP $(basename "$out")" >> "$D/run_all.log"; continue; fi
    printf '%s\n' "$line"
done | xargs -P "$P" -I{} "$D/run_one.sh" {}
echo "ALL_DONE $(date '+%F %T')" >> "$D/run_all.log"
