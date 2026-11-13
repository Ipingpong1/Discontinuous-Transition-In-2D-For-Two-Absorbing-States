#!/bin/sh
# Second extension chain (30 Oct 2026, 20:35): the two L 4096 side runs at p1 = 0.18539 from 6e5 to 9e5
# steps (is the decline of the 0.4e-4 point at L 4096 a decay or a fluctuation?), and seeds 17-18 at
# 0.18539 and 0.18540 (L 2048, run late by the missing-run chain) from 6e5 to 1.2e6, so that the
# 6e5-1.2e6 window has all 8 seeds there.  The list bulk_dens/runs_ext2.txt is written once, by hand
# (same command pattern as make_runs.py bulk()), before this chain starts, and is never rewritten.
D=$(cd "$(dirname "$0")/.." && pwd); cd "$D"; export OMP_NUM_THREADS=1
echo "EXT2_START $(date) P=6 $(wc -l < bulk_dens/runs_ext2.txt) runs" >> bulk_dens/run_dens.log
xargs -P 6 -I{} sh -c "nice -n 10 $D/dnsim {} > /dev/null 2>&1" < bulk_dens/runs_ext2.txt
echo "EXT2_DONE $(date)" >> bulk_dens/run_dens.log
