#!/bin/sh
# Chain for NOTES.md Sec. 9 (seed-averaged densities above the line): bulk_dens/runs.txt, P
# single-thread jobs at nice 10 (12 on the 8 cores + HT, as in run_fine.sh).  Re-runnable: the
# list generator leaves finished runs out and every run resumes from its checkpoint.
D=$(cd "$(dirname "$0")/.." && pwd); cd "$D"
export OMP_NUM_THREADS=1
P=${P:-12}
python3 bulk_dens/make_runs.py --main >> bulk_dens/run_dens.log
echo "START $(date) P=$P $(wc -l < bulk_dens/runs.txt) runs" >> bulk_dens/run_dens.log
xargs -P "$P" -I{} sh -c "nice -n 10 $D/dnsim {} > /dev/null 2>&1" < bulk_dens/runs.txt
echo "DONE $(date)" >> bulk_dens/run_dens.log
