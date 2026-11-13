#!/bin/bash
# Run one shard file (lines "<out prefix>\t<dnsim args>") with up to NCPU / THREADS runs at a
# time, THREADS OpenMP threads each (default 1; the output does not depend on it).  Finished
# runs (output has "# done") are skipped and unfinished ones resume from their checkpoint, so
# re-running a shard is always safe.
#   run_shard.sh SHARD_FILE        (env: SLURM_CPUS_PER_TASK, THREADS, DNSIM)
set -uo pipefail
SHARD=$1
NCPU=${SLURM_CPUS_PER_TASK:-1}
THREADS=${THREADS:-1}
NSLOT=$(( NCPU / THREADS > 0 ? NCPU / THREADS : 1 ))
BIN=${DNSIM:-$HOME/dn/bin/dnsim}
export OMP_NUM_THREADS=$THREADS
echo "host $(hostname) cpus $NCPU threads/run $THREADS shard $SHARD bin $BIN $(date '+%F %T')"
n=0 skip=0
while IFS=$'\t' read -r prefix args; do
    [[ -z $prefix || $prefix == \#* ]] && continue
    if [[ -f $prefix.txt ]] && grep -q '^# done' "$prefix.txt"; then skip=$((skip + 1)); continue; fi
    while (( $(jobs -rp | wc -l) >= NSLOT )); do wait -n; done
    n=$((n + 1))
    (
        # shellcheck disable=SC2086  # args are whitespace-separated dnsim options
        "$BIN" $args > "$prefix.log" 2>&1
        rc=$?
        echo "end $(basename "$prefix") rc=$rc $(date +%T) $(grep '^# done' "$prefix.txt" 2>/dev/null | cut -c3-)"
    ) &
done < "$SHARD"
wait
echo "SHARD_DONE launched $n skipped $skip $(date '+%F %T')"
