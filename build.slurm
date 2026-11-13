#!/bin/bash
#SBATCH -J dn-build
#SBATCH -A torch_pr_1342_general
#SBATCH -p cpu_short
#SBATCH -c 32
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o /home/mo2273/dn/logs/build_%j.log
# Build dnsim (portable x86-64-v3 and AVX-512 x86-64-v4), run the selftest, a
# reproducibility check against a reference run, and a throughput benchmark.
IMG=/share/apps/images/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif
apptainer exec --bind /home "$IMG" bash -s <<'EOF'
set -euo pipefail
cd $HOME/dn
echo "host $(hostname)  $(lscpu | grep 'Model name' | sed 's/ \+/ /g')"
gcc -O3 -march=x86-64-v3 -fopenmp -Wall -Wextra -std=c11 -o bin/dnsim dnsim.c -lm
gcc -O3 -march=x86-64-v4 -fopenmp -Wall -Wextra -std=c11 -o bin/dnsim_v4 dnsim.c -lm
for b in dnsim dnsim_v4; do echo "== selftest $b"; OMP_NUM_THREADS=8 bin/$b --selftest | tail -3; done
echo "== reproducibility (data lines md5, compare with local)"
OMP_NUM_THREADS=32 bin/dnsim --L 512 --T 2000 --alpha 0.5 --jflip 2:0.02 --seed 5 | grep -v '^#' | md5sum
OMP_NUM_THREADS=3 bin/dnsim_v4 --L 512 --T 2000 --alpha 0.5 --jflip 2:0.02 --seed 5 | grep -v '^#' | md5sum
echo "== benchmark L=4096 T=300, alpha=0.5 jflip 2:0.02 (near transition)"
for b in dnsim dnsim_v4; do
  for th in 1 8 16 32; do
    printf "%-9s threads=%-3s " $b $th
    OMP_NUM_THREADS=$th OMP_PROC_BIND=close bin/$b --L 4096 --T 300 --alpha 0.5 --jflip 2:0.02 --ppd 1 | tail -1
  done
done
echo "== 32 independent 1-thread runs at once (ensemble mode), L=2048 T=300"
for i in $(seq 32); do OMP_NUM_THREADS=1 bin/dnsim --L 2048 --T 300 --alpha 0.5 --jflip 2:0.02 --seed $i --ppd 1 | tail -1 & done | awk '{s+=$NF; n++} END{printf "runs %d  mean per-run %.3g  aggregate %.3g updates/s\n", n, s/n, s}'
wait
echo BUILD_DONE
EOF
