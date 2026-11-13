
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gpu_api.h"
#define __device__
#define __forceinline__ inline
#define __ffs(x) __builtin_ffs(x)
#include "generated_device.cuh"
struct DnsimGpu { Model model; int W, Ly; size_t bytes; u64 *a, *b, time, flips; };
DnsimGpu *dnsim_gpu_init(const Model *m, size_t n, int W, int Ly, const u64 *a, u64 t0) {
    if (n != sizeof(Model)) abort();
    DnsimGpu *g = calloc(1, sizeof *g);
    g->model = *m; g->W = W; g->Ly = Ly; g->bytes = (size_t)W * Ly * 8; g->time = t0;
    g->a = malloc(g->bytes); g->b = malloc(g->bytes); memcpy(g->a, a, g->bytes); return g;
}
void dnsim_gpu_advance(DnsimGpu *g, u64 final_time) {
    while (g->time < final_time) {
        ++g->time; g->flips = 0;
        for (int y = 0; y < g->Ly; ++y) {
            const u64 *U = g->a + (size_t)(y ? y - 1 : g->Ly - 1) * g->W;
            const u64 *C = g->a + (size_t)y * g->W;
            const u64 *D = g->a + (size_t)(y + 1 < g->Ly ? y + 1 : 0) * g->W;
            for (int x = 0; x < g->W; ++x) {
                u64 f = step_word(U,C,D,x,x ? x-1 : g->W-1,x+1 < g->W ? x+1 : 0,
                                  (u32)g->time,y,&g->model);
                g->b[(size_t)y*g->W+x] = C[x]^f; g->flips += __builtin_popcountll(f);
            }
        }
        u64 *swap = g->a; g->a = g->b; g->b = swap;
    }
}
void dnsim_gpu_download(DnsimGpu *g, u64 *state, u64 *flips) {
    memcpy(state, g->a, g->bytes); if (flips) *flips = g->flips;
}
void dnsim_gpu_finish(DnsimGpu *g) { free(g->a); free(g->b); free(g); }
