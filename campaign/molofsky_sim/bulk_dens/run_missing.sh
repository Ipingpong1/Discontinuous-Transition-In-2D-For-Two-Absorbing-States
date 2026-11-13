#!/bin/sh
# The 16 main runs the main chain never started (its list, runs.txt, was rewritten under the running
# xargs at 17:00, see make_runs.py).  Started at 18:50 beside the extension chain (12 jobs) and the two
# L 4096 side runs with P0 slots; once the extension chain has written EXT_DONE the slots are raised
# to P1 by sending SIGUSR1 to xargs (GNU xargs adds one parallel slot per SIGUSR1).
D=$(cd "$(dirname "$0")/.." && pwd); cd "$D"
export OMP_NUM_THREADS=1
P0=${P0:-6}; P1=${P1:-12}
python3 bulk_dens/make_runs.py --missing >> bulk_dens/run_dens.log
echo "MISS_START $(date) P=$P0 (-> $P1 after EXT_DONE) $(wc -l < bulk_dens/runs_missing.txt) runs" >> bulk_dens/run_dens.log
xargs -P "$P0" -I{} sh -c "nice -n 10 $D/dnsim {} > /dev/null 2>&1" < bulk_dens/runs_missing.txt &
XP=$!
( until grep -q '^EXT_DONE' bulk_dens/run_dens.log; do sleep 60; done
  i=$P0; while [ "$i" -lt "$P1" ]; do kill -USR1 "$XP" 2>/dev/null; i=$((i + 1)); sleep 1; done
  echo "MISS_SLOTS_RAISED $(date) to $P1" >> bulk_dens/run_dens.log ) &
wait "$XP"
echo "MISS_DONE $(date)" >> bulk_dens/run_dens.log
