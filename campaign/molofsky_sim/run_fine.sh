#!/bin/sh
# Chain for NOTES.md Sec. 8: runs_fine.txt (strips s1-2, lifetimes, small bulks, strips s3-4,
# L 2048 bulks, L 1024 1e6 density runs) then runs_fine_wide.txt (Lx 4096 strips), P single-thread
# jobs at nice 10 (12 jobs on 8 cores + HT gave ~25 % more throughput than 8 in a 10 s test).
D=$(cd "$(dirname "$0")" && pwd); cd "$D"
export OMP_NUM_THREADS=1
P=${P:-12}
echo "START $(date)" >> run_fine.log
xargs -P "$P" -I{} sh -c "nice -n 10 $D/dnsim {} > /dev/null 2>&1" < runs_fine.txt
echo "MAIN_DONE $(date)" >> run_fine.log
xargs -P "$P" -I{} sh -c "nice -n 10 $D/dnsim {} > /dev/null 2>&1" < runs_fine_wide.txt
echo "WIDE_DONE $(date)" >> run_fine.log
