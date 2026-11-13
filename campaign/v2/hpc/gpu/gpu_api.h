#ifndef DNSIM_GPU_API_H
#define DNSIM_GPU_API_H
#include <stddef.h>
#include "generated_model.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct DnsimGpu DnsimGpu;
DnsimGpu *dnsim_gpu_init(const Model *model, size_t model_bytes, int W, int Ly,
                         const u64 *initial, u64 t0);
void dnsim_gpu_advance(DnsimGpu *gpu, u64 final_time);
void dnsim_gpu_download(DnsimGpu *gpu, u64 *state, u64 *last_flips);
void dnsim_gpu_finish(DnsimGpu *gpu);
#ifdef __cplusplus
}
#endif
#endif
