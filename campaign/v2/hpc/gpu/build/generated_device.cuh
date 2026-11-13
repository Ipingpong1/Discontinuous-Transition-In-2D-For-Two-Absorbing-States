/* Generated from dnsim.c; SHA256 dd56be11b1e7865b52ec4607432eb5733f2306b4a714312690df801793f22784. Do not edit. */
#include "generated_model.h"
typedef struct { u32 ctr[4], key[2], buf[4]; int left; } rstream;
static __device__ __forceinline__ int gpu_ctz(u32 v) { return __ffs((int)v) - 1; }
static __device__ __forceinline__ void philox4x32_10(const u32 in[4], const u32 key[2], u32 out[4])
{
    u32 c0 = in[0], c1 = in[1], c2 = in[2], c3 = in[3];
    u32 k0 = key[0], k1 = key[1];
    for (int r = 0; r < 10; r++) {
        if (r) { k0 += 0x9E3779B9u; k1 += 0xBB67AE85u; }
        u64 p0 = (u64)0xD2511F53u * c0;
        u64 p1 = (u64)0xCD9E8D57u * c2;
        u32 n0 = (u32)(p1 >> 32) ^ c1 ^ k0;
        u32 n2 = (u32)(p0 >> 32) ^ c3 ^ k1;
        c1 = (u32)p1; c3 = (u32)p0; c0 = n0; c2 = n2;
    }
    out[0] = c0; out[1] = c1; out[2] = c2; out[3] = c3;
}

static __device__ __forceinline__ void rs_init(rstream *r, const u32 key[2], u32 a, u32 b, u32 c)
{
    r->ctr[0] = a; r->ctr[1] = b; r->ctr[2] = c; r->ctr[3] = 0;
    r->key[0] = key[0]; r->key[1] = key[1];
    r->left = 0;
}

static __device__ __forceinline__ u64 rs_next(rstream *r)
{
    if (!r->left) { philox4x32_10(r->ctr, r->key, r->buf); r->ctr[3]++; r->left = 2; }
    r->left--;
    int i = r->left ? 0 : 2;
    return ((u64)r->buf[i + 1] << 32) | r->buf[i];
}

static __device__ __forceinline__ u64 bern_table(const Model *M, const u64 *m, u64 lanes, rstream *r)
{
    u64 und = lanes, res = 0;
    for (int b = 31; b >= M->minlow; b--) {
        u64 R = rs_next(r), P = 0;
        for (u32 cs = M->bitcls[b]; cs; cs &= cs - 1) P |= m[gpu_ctz(cs)];
        res |= und & P & ~R;
        und &= ~(P ^ R);
        for (u32 cs = M->endcls[b]; cs; cs &= cs - 1) und &= ~m[gpu_ctz(cs)];
        if (!und) break;
    }
    return res;
}

static __device__ __forceinline__ void count_eq(const u64 *U, const u64 *C, const u64 *D, int x, int xl, int xr, u64 eq[9])
{
    u64 u = U[x], c = C[x], d = D[x];
    u64 uw = (u << 1) | (U[xl] >> 63), ue = (u >> 1) | (U[xr] << 63);
    u64 cw = (c << 1) | (C[xl] >> 63), ce = (c >> 1) | (C[xr] << 63);
    u64 dw = (d << 1) | (D[xl] >> 63), de = (d >> 1) | (D[xr] << 63);
    u64 su = uw ^ u ^ ue, cu = (uw & u) | (ue & (uw ^ u));
    u64 sd = dw ^ d ^ de, cd = (dw & d) | (de & (dw ^ d));
    u64 sm = cw ^ ce, cm = cw & ce;
    u64 b0 = su ^ sd ^ sm, y1 = (su & sd) | (sm & (su ^ sd));
    u64 z1 = cu ^ cd ^ cm, z2 = (cu & cd) | (cm & (cu ^ cd));
    u64 b1 = z1 ^ y1, y2 = z1 & y1;
    u64 b2 = z2 ^ y2, b3 = z2 & y2;
    u64 nb3 = ~b3;
    u64 t00 = ~b1 & ~b0, t01 = ~b1 & b0, t10 = b1 & ~b0, t11 = b1 & b0;
    u64 lo = nb3 & ~b2, hi = nb3 & b2;
    eq[0] = lo & t00; eq[1] = lo & t01; eq[2] = lo & t10; eq[3] = lo & t11;
    eq[4] = hi & t00; eq[5] = hi & t01; eq[6] = hi & t10; eq[7] = hi & t11;
    eq[8] = b3;
}

static __device__ __forceinline__ void count_eq_vn(const u64 *U, const u64 *C, const u64 *D, int x, int xl, int xr, u64 eq[9])
{
    u64 u = U[x], c = C[x], d = D[x];
    u64 cw = (c << 1) | (C[xl] >> 63), ce = (c >> 1) | (C[xr] << 63);
    u64 s1 = u ^ d, c1 = u & d, s2 = cw ^ ce, c2 = cw & ce;
    u64 b0 = s1 ^ s2, cc = s1 & s2;
    u64 b1 = c1 ^ c2 ^ cc, b2 = (c1 & c2) | (cc & (c1 ^ c2));
    u64 nb = ~b2;
    eq[0] = nb & ~b1 & ~b0; eq[1] = nb & ~b1 & b0; eq[2] = nb & b1 & ~b0; eq[3] = nb & b1 & b0;
    eq[4] = b2; eq[5] = eq[6] = eq[7] = eq[8] = 0;
}

static __device__ __forceinline__ u64
step_word(const u64 *U, const u64 *C, const u64 *D, int x, int xl, int xr, u32 t, int y, const Model *M)
{
    u64 eq[9], m[2 + MAXV] = {0};
    if (M->vn) count_eq_vn(U, C, D, x, xl, xr, eq);
    else count_eq(U, C, D, x, xl, xr, eq);
    u64 c = C[x];
    for (int k = 0; k < 9; k++) {
        m[M->slot[0][k]] |= ~c & eq[k];
        m[M->slot[1][k]] |= c & eq[k];
    }
    u64 flip = m[SLOT_ALWAYS], rnd = 0;
    for (int v = 0; v < M->nv; v++) rnd |= m[2 + v];
    if (rnd) {
        rstream r;
        rs_init(&r, M->key, t, (u32)y, (u32)x);
        flip |= bern_table(M, m + 2, rnd, &r);
    }
    return flip;
}
