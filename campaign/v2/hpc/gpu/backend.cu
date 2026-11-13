/* Dense CUDA backend. Integer update helpers are generated from dnsim.c. */
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include "gpu_api.h"
#include "generated_device.cuh"

#ifndef DNSIM_GRAPH_BATCH
#define DNSIM_GRAPH_BATCH 256
#endif
static_assert(DNSIM_GRAPH_BATCH >= 2 && DNSIM_GRAPH_BATCH % 2 == 0,
              "graph batch must be a positive even number");

#define CUDA_CHECK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
    std::fprintf(stderr, "CUDA error %s:%d: %s: %s\n", __FILE__, __LINE__, #call, \
                 cudaGetErrorString(e_)); std::exit(1); } } while (0)

struct DnsimGpu {
    u64 *state[2], *base;
    u32 *flips, *host_flips;
    size_t words, bytes;
    int W, Ly, parity, blocks;
    u64 time;
    cudaStream_t stream;
    cudaGraph_t graph[2];
    cudaGraphExec_t executable[2];
};
/* One model per process; no floating-point arithmetic occurs in the kernels. */
__constant__ Model gpu_model;

__global__ static void update_words(const u64 *cur, u64 *next, int W, int Ly,
                                    const u64 *base, u32 offset, u32 *counts)
{
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t words = (size_t)W * Ly;
    if (i >= words) return;
    int y = (int)(i / W), x = (int)(i % W);
    const u64 *U = cur + (size_t)(y ? y - 1 : Ly - 1) * W;
    const u64 *C = cur + (size_t)y * W;
    const u64 *D = cur + (size_t)(y + 1 < Ly ? y + 1 : 0) * W;
    int xl = x ? x - 1 : W - 1, xr = x + 1 < W ? x + 1 : 0;
    u64 flip = step_word(U, C, D, x, xl, xr, (u32)(*base + offset), y, &gpu_model);
    next[i] = C[x] ^ flip;
    if (counts) counts[i] = (u32)__popcll((unsigned long long)flip);
}
__global__ static void increment_time(u64 *base, u32 count) { *base += count; }

extern "C" DnsimGpu *dnsim_gpu_init(const Model *model, size_t model_bytes,
                                    int W, int Ly, const u64 *initial, u64 t0)
{
    if (model_bytes != sizeof(Model) || W <= 0 || Ly < 3) {
        std::fprintf(stderr, "Invalid CUDA model layout or lattice\n"); std::exit(2);
    }
    DnsimGpu *g = new DnsimGpu{};
    g->W = W; g->Ly = Ly; g->words = (size_t)W * Ly;
    g->bytes = g->words * sizeof(u64);
    g->parity = 0; g->time = t0;
    const int threads = 128;
    if ((g->words + threads - 1) / threads > INT_MAX) {
        std::fprintf(stderr, "CUDA lattice exceeds one-dimensional launch range\n"); std::exit(2);
    }
    g->blocks = (int)((g->words + threads - 1) / threads);
    CUDA_CHECK(cudaStreamCreateWithFlags(&g->stream, cudaStreamNonBlocking));
    CUDA_CHECK(cudaMalloc(&g->state[0], g->bytes));
    CUDA_CHECK(cudaMalloc(&g->state[1], g->bytes));
    CUDA_CHECK(cudaMalloc(&g->base, sizeof(u64)));
    CUDA_CHECK(cudaMalloc(&g->flips, g->words * sizeof(u32)));
    g->host_flips = (u32 *)std::malloc(g->words * sizeof(u32));
    if (!g->host_flips) { std::fprintf(stderr, "Out of host memory\n"); std::exit(1); }
    CUDA_CHECK(cudaMemcpyToSymbol(gpu_model, model, sizeof(Model)));
    CUDA_CHECK(cudaMemcpy(g->state[0], initial, g->bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g->base, &t0, sizeof(u64), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(g->flips, 0, g->words * sizeof(u32)));
    for (int start = 0; start < 2; start++) {
        CUDA_CHECK(cudaStreamBeginCapture(g->stream, cudaStreamCaptureModeThreadLocal));
        for (u32 k = 1; k <= DNSIM_GRAPH_BATCH; k++) {
            int from = start ^ ((k - 1) & 1);
            update_words<<<g->blocks, threads, 0, g->stream>>>(
                g->state[from], g->state[from ^ 1], W, Ly, g->base, k,
                k == DNSIM_GRAPH_BATCH ? g->flips : nullptr);
        }
        increment_time<<<1, 1, 0, g->stream>>>(g->base, DNSIM_GRAPH_BATCH);
        CUDA_CHECK(cudaStreamEndCapture(g->stream, &g->graph[start]));
        CUDA_CHECK(cudaGraphInstantiateWithFlags(&g->executable[start], g->graph[start], 0));
    }
    int device;
    cudaDeviceProp properties;
    CUDA_CHECK(cudaGetDevice(&device));
    CUDA_CHECK(cudaGetDeviceProperties(&properties, device));
    std::fprintf(stderr, "CUDA dense backend: %s, %zu words, graph batch %d\n",
                 properties.name, g->words, DNSIM_GRAPH_BATCH);
    return g;
}
extern "C" void dnsim_gpu_advance(DnsimGpu *g, u64 final_time)
{
    if (final_time < g->time || final_time > UINT32_MAX) {
        std::fprintf(stderr, "Invalid CUDA target time\n"); std::exit(2);
    }
    u64 count = final_time - g->time;
    while (count >= DNSIM_GRAPH_BATCH) {
        CUDA_CHECK(cudaGraphLaunch(g->executable[g->parity], g->stream));
        count -= DNSIM_GRAPH_BATCH;
    }
    u32 tail = (u32)count;
    for (u32 k = 1; k <= tail; k++) {
        update_words<<<g->blocks, 128, 0, g->stream>>>(
            g->state[g->parity], g->state[g->parity ^ 1], g->W, g->Ly, g->base, k,
            k == tail ? g->flips : nullptr);
        g->parity ^= 1;
    }
    if (tail) increment_time<<<1, 1, 0, g->stream>>>(g->base, tail);
    CUDA_CHECK(cudaGetLastError());
    g->time = final_time;
}
extern "C" void dnsim_gpu_download(DnsimGpu *g, u64 *state, u64 *last_flips)
{
    CUDA_CHECK(cudaMemcpyAsync(state, g->state[g->parity], g->bytes,
                              cudaMemcpyDeviceToHost, g->stream));
    if (last_flips)
        CUDA_CHECK(cudaMemcpyAsync(g->host_flips, g->flips, g->words * sizeof(u32),
                                  cudaMemcpyDeviceToHost, g->stream));
    CUDA_CHECK(cudaStreamSynchronize(g->stream));
    if (last_flips) {
        u64 total = 0;
        for (size_t i = 0; i < g->words; i++) total += g->host_flips[i];
        *last_flips = total;
    }
}
extern "C" void dnsim_gpu_finish(DnsimGpu *g)
{
    CUDA_CHECK(cudaStreamSynchronize(g->stream));
    for (int i = 0; i < 2; i++) {
        CUDA_CHECK(cudaGraphExecDestroy(g->executable[i]));
        CUDA_CHECK(cudaGraphDestroy(g->graph[i]));
        CUDA_CHECK(cudaFree(g->state[i]));
    }
    CUDA_CHECK(cudaFree(g->base));
    CUDA_CHECK(cudaFree(g->flips));
    CUDA_CHECK(cudaStreamDestroy(g->stream));
    std::free(g->host_flips);
    delete g;
}
