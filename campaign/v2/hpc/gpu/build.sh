#!/usr/bin/env bash
set -euo pipefail
GPU_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-$GPU_DIR/build}"
CPU_SOURCE="${CPU_SOURCE:-$GPU_DIR/../../../round2/dnsim.c}"
CC="${CC:-gcc}"
NVCC="${NVCC:-nvcc}"
CUDA_ARCH="${CUDA_ARCH:-80}"
python3 "$GPU_DIR/generate.py" --source "$CPU_SOURCE" --out "$BUILD_DIR"
"$CC" -O3 -std=gnu11 -fopenmp -I"$GPU_DIR" -I"$BUILD_DIR" \
  -c "$BUILD_DIR/generated_host.c" -o "$BUILD_DIR/host.o"
"$NVCC" -O3 -std=c++14 -I"$GPU_DIR" -I"$BUILD_DIR" \
  -gencode "arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH}" \
  -gencode "arch=compute_${CUDA_ARCH},code=compute_${CUDA_ARCH}" \
  -c "$GPU_DIR/backend.cu" -o "$BUILD_DIR/backend.o"
"$NVCC" -Xcompiler=-fopenmp "$BUILD_DIR/host.o" "$BUILD_DIR/backend.o" \
  -lm -o "$BUILD_DIR/dnsim_gpu"
printf '%s\n' "$BUILD_DIR/dnsim_gpu"
