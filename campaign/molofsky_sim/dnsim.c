/*
 * dnsim — bit-packed stochastic Life-like cellular automaton.
 *
 * Built for Day & Night (B3678/S34678) under Z2-symmetric, absorbing-state-preserving
 * noise, but the model is any synchronous outer-totalistic stochastic CA on an
 * Lx x Ly torus (Moore neighbourhood): a cell in state s with k live neighbours flips
 * with probability f[s][k].  The flip table is built from
 *
 *     f = (1 - voter) * alpha * [rule says flip]  +  voter * j/8  +  sum_j' q_j' [j == j']
 *
 * where j = number of DISAGREEING neighbours (k if s=0, 8-k if s=1), clipped to [0,1];
 * or given verbatim with --ftab.  For a self-complementary rule (D&N) every term
 * depends on j only, so the dynamics commute with complementation, and f(j=0) = 0
 * keeps all-0 and all-1 absorbing.
 *
 * Layout: 64 cells per uint64 word, bit b of word x in row y = cell (64x+b, y).
 * Randomness: counter-based Philox4x32-10 keyed by --seed, counter (t, y, x, i), so a
 * run is bit-reproducible for any thread count and restartable from a checkpoint.
 * Flip probabilities are quantised to multiples of 2^-32; each cell draws one uniform
 * 32-bit U per step and flips iff U < p of its class, exactly for the quantised value
 * (bit-serial comparison, all classes of a word at once).  v1.1 used one stream per
 * class (counter x|v<<24), which is ~4x slower for tables with many classes (voter
 * mixing); the dynamics are the same in law, the random streams differ.
 *
 * v1.4 (10-30, campaign/molofsky_sim): --nbhd vn selects the von Neumann neighbourhood
 * (k = live N, S, W, E neighbours, 0..4; table entries k >= 5 unused) and --molofsky p1,p2
 * builds the five-cell rule of Molofsky, Durrett, Dushoff, Griffeath & Levin, Theor. Popul.
 * Biol. 55, 270 (1999) on it.  The Moore path is unchanged.
 */
#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef uint64_t u64;
typedef uint32_t u32;

/* ------------------------------------------------------------------ RNG */

/* Philox4x32-10 (Salmon, Moraes, Dror, Shaw, SC'11), as in Random123. */
static inline void philox4x32_10(const u32 in[4], const u32 key[2], u32 out[4])
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

typedef struct { u32 ctr[4], key[2], buf[4]; int left; } rstream;

static inline void rs_init(rstream *r, const u32 key[2], u32 a, u32 b, u32 c)
{
    r->ctr[0] = a; r->ctr[1] = b; r->ctr[2] = c; r->ctr[3] = 0;
    r->key[0] = key[0]; r->key[1] = key[1];
    r->left = 0;
}

static inline u64 rs_next(rstream *r)
{
    if (!r->left) { philox4x32_10(r->ctr, r->key, r->buf); r->ctr[3]++; r->left = 2; }
    r->left--;
    int i = r->left ? 0 : 2;
    return ((u64)r->buf[i + 1] << 32) | r->buf[i];
}

/* Independent Bernoulli(p / 2^32) bits on the lanes set in `lanes` (other lanes are
 * garbage), 0 < p < 2^32, low = ctz(p).  Bit-serial comparison U < p, MSB first; stops
 * once every requested lane is decided or the lowest set bit of p has been passed
 * (after it, U == p-prefix can only mean U >= p).  Restricting to the requested lanes
 * cuts the expected number of draws from ~log2(64) to ~log2(popcount(lanes)). */
static inline u64 bern(u32 p, int low, u64 lanes, rstream *r)
{
    u64 und = lanes, res = 0;
    for (int b = 31; b >= low; b--) {
        u64 R = rs_next(r);
        if ((p >> b) & 1) { res |= und & ~R; und &= R; }
        else und &= ~R;
        if (!und) break;
    }
    return res;
}

/* ---------------------------------------------------------------- model */

#define MAXV 18
enum { SLOT_NEVER = 0, SLOT_ALWAYS = 1 };

typedef struct {
    double f[2][9];   /* requested flip probabilities */
    int slot[2][9];   /* SLOT_NEVER, SLOT_ALWAYS, or 2+v */
    int nv;
    u32 p[MAXV];      /* quantised probability of value v, units of 2^-32 */
    int low[MAXV];
    u32 bitcls[32];   /* classes v whose p has bit b set */
    u32 endcls[32];   /* classes v whose lowest set bit is b */
    int minlow;
    u32 key[2];
    int vn;           /* 1: von Neumann neighbourhood (k in 0..4), 0: Moore */
} Model;

static void model_quantise(Model *M)
{
    M->nv = 0;
    for (int s = 0; s < 2; s++)
        for (int k = 0; k < 9; k++) {
            double q = M->f[s][k];
            if (q < 0) q = 0;
            if (q > 1) q = 1;
            long long pq = llround(q * 4294967296.0);
            if (pq <= 0) {
                if (q > 0) fprintf(stderr, "warning: f[%d][%d]=%g rounds to 0\n", s, k, q);
                M->slot[s][k] = SLOT_NEVER;
            } else if (pq >= 4294967296LL) {
                if (q < 1) fprintf(stderr, "warning: f[%d][%d]=%.12g rounds to 1\n", s, k, q);
                M->slot[s][k] = SLOT_ALWAYS;
            } else {
                int v;
                for (v = 0; v < M->nv; v++) if (M->p[v] == (u32)pq) break;
                if (v == M->nv) { M->p[v] = (u32)pq; M->low[v] = __builtin_ctz((u32)pq); M->nv++; }
                M->slot[s][k] = 2 + v;
            }
        }
    memset(M->bitcls, 0, sizeof M->bitcls);
    memset(M->endcls, 0, sizeof M->endcls);
    M->minlow = 32;
    for (int v = 0; v < M->nv; v++) {
        for (int b = 0; b < 32; b++) if ((M->p[v] >> b) & 1) M->bitcls[b] |= 1u << v;
        M->endcls[M->low[v]] |= 1u << v;
        if (M->low[v] < M->minlow) M->minlow = M->low[v];
    }
}

/* Lane i of class v (lanes of m[v], disjoint) flips iff U_i < p[v], one uniform U_i per
 * lane: bit-serial comparison, MSB first, where the threshold bit-plane at bit b is the
 * union of the classes whose p has bit b set.  A lane leaves once U_i and p differ, or
 * after the lowest set bit of its p (then U_i >= p).  Expected draws ~log2(popcount(lanes))
 * + 2, independent of the number of classes. */
static inline u64 bern_table(const Model *M, const u64 *m, u64 lanes, rstream *r)
{
    u64 und = lanes, res = 0;
    for (int b = 31; b >= M->minlow; b--) {
        u64 R = rs_next(r), P = 0;
        for (u32 cs = M->bitcls[b]; cs; cs &= cs - 1) P |= m[__builtin_ctz(cs)];
        res |= und & P & ~R;
        und &= ~(P ^ R);
        for (u32 cs = M->endcls[b]; cs; cs &= cs - 1) und &= ~m[__builtin_ctz(cs)];
        if (!und) break;
    }
    return res;
}

/* f = (1 - voter) alpha [rule flips] + voter j/8 + extra[j], clipped to [0, 1]; Bm/Sm are
 * bit masks of the birth/survival counts. */
static void model_from_rule(Model *M, int Bm, int Sm, double alpha, double voter, const double extra[9])
{
    for (int s = 0; s < 2; s++)
        for (int k = 0; k < 9; k++) {
            int j = s ? 8 - k : k;
            double det = s ? !((Sm >> k) & 1) : (Bm >> k) & 1;
            double f = (1 - voter) * alpha * det + voter * j / 8.0 + extra[j];
            M->f[s][k] = f < 0 ? 0 : f > 1 ? 1 : f;
        }
}

static double model_eff(const Model *M, int s, int k)
{
    int sl = M->slot[s][k];
    return sl == SLOT_NEVER ? 0.0 : sl == SLOT_ALWAYS ? 1.0 : M->p[sl - 2] / 4294967296.0;
}

/* Five-cell rule of Molofsky et al. (1999): square lattice, synchronous; the next state of a
 * site is 1 with probability p_n, n = number of 1s among the site and its four nearest
 * neighbours, with p_0 = 0, p_5 = 1, p_3 = 1 - p_2, p_4 = 1 - p_1.  Over k = live von Neumann
 * neighbours: a 0 has n = k and flips with p_k; a 1 has n = k + 1 and flips with 1 - p_{k+1}
 * = p_{4-k}, so f[1][k] = f[0][4-k] (exact Z2 symmetry) and f(j = 0) = 0 for both states (both
 * uniform states absorbing).  Voter point (p_1, p_2) = (1/5, 2/5): f = j/5. */
static void model_molofsky(Model *M, double p1, double p2)
{
    double p[6] = { 0, p1, p2, 1 - p2, 1 - p1, 1 };
    M->vn = 1;
    for (int k = 0; k < 9; k++) {
        M->f[0][k] = k <= 4 ? p[k] : 0;
        M->f[1][k] = k <= 4 ? 1 - p[k + 1] : 0;
    }
}

/* Highest neighbour count of the model's neighbourhood. */
static inline int model_kmax(const Model *M) { return M->vn ? 4 : 8; }

/* ------------------------------------------------------------- lattice */

typedef struct { int Lx, Ly, W; u64 *a; } Grid;

static Grid grid_new(int Lx, int Ly)
{
    Grid g = { Lx, Ly, Lx / 64, NULL };
    g.a = aligned_alloc(64, ((size_t)g.W * Ly * 8 + 63) / 64 * 64);
    if (!g.a) { fprintf(stderr, "out of memory\n"); exit(1); }
    return g;
}

/* Bit-sliced neighbour count of the word at column-word x: k = b0 + 2 b1 + 4 b2 + 8 b3,
 * decoded into eq[k] (bit set where the count is exactly k). */
static inline void count_eq(const u64 *U, const u64 *C, const u64 *D, int x, int xl, int xr, u64 eq[9])
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

/* Von Neumann count k = N + S + W + E (0..4), bit-sliced: eq[k] set where the count is
 * exactly k; eq[5..8] = 0. */
static inline void count_eq_vn(const u64 *U, const u64 *C, const u64 *D, int x, int xl, int xr, u64 eq[9])
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

/* Flip mask of column-word x in row y at time t (rows U, C, D = y-1, y, y+1).  Depends only
 * on the 3x3 word neighbourhood and the counter (t, y, x), so stepping any subset of words
 * gives the same result as stepping the whole lattice. */
static inline __attribute__((always_inline)) u64
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

/* Cells of word x with at least one disagreeing Moore neighbour. */
static inline u64 iface_word(const u64 *U, const u64 *C, const u64 *D, int x, int xl, int xr)
{
    u64 u = U[x], c = C[x], d = D[x];
    u64 uw = (u << 1) | (U[xl] >> 63), ue = (u >> 1) | (U[xr] << 63);
    u64 cw = (c << 1) | (C[xl] >> 63), ce = (c >> 1) | (C[xr] << 63);
    u64 dw = (d << 1) | (D[xl] >> 63), de = (d >> 1) | (D[xr] << 63);
    return (c ^ uw) | (c ^ u) | (c ^ ue) | (c ^ cw) | (c ^ ce) | (c ^ dw) | (c ^ d) | (c ^ de);
}

/* One synchronous update at time t (t >= 1). Returns the number of flipped cells;
 * row_flips (optional) receives per-row counts. */
static u64 step(const Grid *g, u64 *restrict nxt, u32 t, const Model *M, u32 *row_flips)
{
    const int W = g->W, Ly = g->Ly;
    const u64 *cur = g->a;
    u64 total = 0;
#pragma omp parallel for schedule(static) reduction(+ : total)
    for (int y = 0; y < Ly; y++) {
        const u64 *U = cur + (size_t)(y ? y - 1 : Ly - 1) * W;
        const u64 *C = cur + (size_t)y * W;
        const u64 *D = cur + (size_t)(y + 1 < Ly ? y + 1 : 0) * W;
        u64 *N = nxt + (size_t)y * W;
        u64 rf = 0;
        for (int x = 0; x < W; x++) {
            int xl = x ? x - 1 : W - 1, xr = x + 1 < W ? x + 1 : 0;
            u64 flip = step_word(U, C, D, x, xl, xr, t, y, M);
            N[x] = C[x] ^ flip;
            rf += __builtin_popcountll(flip);
        }
        if (row_flips) row_flips[y] = (u32)rf;
        total += rf;
    }
    return total;
}

/* ------------------------------------------------- bounding-box stepping (--reps) */

/* Rows ya..yb and column-words xa..xb, inclusive and not wrapping.  Invariant kept by
 * Rep: every cell with a disagreeing neighbour lies in the box, so every cell outside has
 * j = 0 and cannot flip (f(0) = 0 for an absorbing model), and both buffers agree outside
 * the box.  Stepping the box alone is then bit-identical to stepping the whole lattice. */
typedef struct { int ya, yb, xa, xb; } Box;

typedef struct { u64 flips, born; int fy0, fy1, fx0, fx1; } BoxStep;

static BoxStep step_box(const Grid *g, u64 *restrict nxt, u32 t, const Model *M, Box b)
{
    const int W = g->W, Ly = g->Ly;
    const u64 *cur = g->a;
    u64 total = 0, born = 0;
    int fy0 = Ly, fy1 = -1, fx0 = W, fx1 = -1;
    long work = (long)(b.yb - b.ya + 1) * (b.xb - b.xa + 1);
#pragma omp parallel for schedule(static) reduction(+ : total, born) reduction(min : fy0, fx0) \
    reduction(max : fy1, fx1) if (work >= 16384)
    for (int y = b.ya; y <= b.yb; y++) {
        const u64 *U = cur + (size_t)(y ? y - 1 : Ly - 1) * W;
        const u64 *C = cur + (size_t)y * W;
        const u64 *D = cur + (size_t)(y + 1 < Ly ? y + 1 : 0) * W;
        u64 *N = nxt + (size_t)y * W;
        for (int x = b.xa; x <= b.xb; x++) {
            int xl = x ? x - 1 : W - 1, xr = x + 1 < W ? x + 1 : 0;
            u64 flip = step_word(U, C, D, x, xl, xr, t, y, M);
            N[x] = C[x] ^ flip;
            if (flip) {
                total += __builtin_popcountll(flip);
                born += __builtin_popcountll(flip & ~C[x]);
                if (y < fy0) fy0 = y;
                if (y > fy1) fy1 = y;
                if (x < fx0) fx0 = x;
                if (x > fx1) fx1 = x;
            }
        }
    }
    return (BoxStep){ total, born, fy0, fy1, fx0, fx1 };
}

/* Exact word bounding box of the cells with a disagreeing neighbour inside b; 0 if none. */
static int iface_bbox(const Grid *g, Box b, Box *out)
{
    const int W = g->W, Ly = g->Ly;
    const u64 *a = g->a;
    Box r = { Ly, -1, W, -1 };
    for (int y = b.ya; y <= b.yb; y++) {
        const u64 *U = a + (size_t)(y ? y - 1 : Ly - 1) * W;
        const u64 *C = a + (size_t)y * W;
        const u64 *D = a + (size_t)(y + 1 < Ly ? y + 1 : 0) * W;
        for (int x = b.xa; x <= b.xb; x++) {
            int xl = x ? x - 1 : W - 1, xr = x + 1 < W ? x + 1 : 0;
            if (iface_word(U, C, D, x, xl, xr)) {
                if (y < r.ya) r.ya = y;
                if (y > r.yb) r.yb = y;
                if (x < r.xa) r.xa = x;
                if (x > r.xb) r.xb = x;
            }
        }
    }
    *out = r;
    return r.yb >= 0;
}

typedef struct { u64 alive, iface; double r2; } BoxObs;

/* Live cells, interface cells and the sum over interface cells of their squared distance
 * from the seed centre (Lx/2, Ly/2), inside b (all live cells are inside b). */
static BoxObs measure_box(const Grid *g, Box b)
{
    const int W = g->W, Ly = g->Ly;
    const u64 *a = g->a;
    const double cx = g->Lx / 2.0, cy = Ly / 2.0;
    u64 alive = 0, iface = 0;
    double r2 = 0;
    for (int y = b.ya; y <= b.yb; y++) {
        const u64 *U = a + (size_t)(y ? y - 1 : Ly - 1) * W;
        const u64 *C = a + (size_t)y * W;
        const u64 *D = a + (size_t)(y + 1 < Ly ? y + 1 : 0) * W;
        const double dy = y + 0.5 - cy;
        for (int x = b.xa; x <= b.xb; x++) {
            int xl = x ? x - 1 : W - 1, xr = x + 1 < W ? x + 1 : 0;
            alive += __builtin_popcountll(C[x]);
            u64 m = iface_word(U, C, D, x, xl, xr);
            int n = __builtin_popcountll(m);
            iface += n;
            r2 += n * dy * dy;
            for (; m; m &= m - 1) {
                double dx = 64.0 * x + __builtin_ctzll(m) + 0.5 - cx;
                r2 += dx * dx;
            }
        }
    }
    return (BoxObs){ alive, iface, r2 };
}

typedef struct { u64 alive, bonds, iface; } Obs;

/* alive cells, disagreeing nearest-neighbour bonds (each bond once), cells with >= 1
 * disagreeing Moore neighbour; per-row alive/iface counts if requested. */
static Obs measure(const Grid *g, u32 *row_alive, u32 *row_iface)
{
    const int W = g->W, Ly = g->Ly;
    const u64 *a = g->a;
    u64 alive = 0, bonds = 0, iface = 0;
#pragma omp parallel for schedule(static) reduction(+ : alive, bonds, iface)
    for (int y = 0; y < Ly; y++) {
        const u64 *U = a + (size_t)(y ? y - 1 : Ly - 1) * W;
        const u64 *C = a + (size_t)y * W;
        const u64 *D = a + (size_t)(y + 1 < Ly ? y + 1 : 0) * W;
        u64 ra = 0, ri = 0, rb = 0;
        for (int x = 0; x < W; x++) {
            int xl = x ? x - 1 : W - 1, xr = x + 1 < W ? x + 1 : 0;
            u64 u = U[x], c = C[x], d = D[x];
            u64 uw = (u << 1) | (U[xl] >> 63), ue = (u >> 1) | (U[xr] << 63);
            u64 cw = (c << 1) | (C[xl] >> 63), ce = (c >> 1) | (C[xr] << 63);
            u64 dw = (d << 1) | (D[xl] >> 63), de = (d >> 1) | (D[xr] << 63);
            ra += __builtin_popcountll(c);
            rb += __builtin_popcountll(c ^ cw) + __builtin_popcountll(c ^ u);
            ri += __builtin_popcountll((c ^ uw) | (c ^ u) | (c ^ ue) | (c ^ cw) | (c ^ ce) |
                                       (c ^ dw) | (c ^ d) | (c ^ de));
        }
        if (row_alive) row_alive[y] = (u32)ra;
        if (row_iface) row_iface[y] = (u32)ri;
        alive += ra; bonds += rb; iface += ri;
    }
    return (Obs){ alive, bonds, iface };
}

/* --------------------------------------------------------- init / io */

static u64 random_word(const u32 key[2], int y, int x, u32 pq)
{
    if (pq == 0) return 0;
    rstream r;
    rs_init(&r, key, 0, (u32)y, (u32)x | (255u << 24));
    return bern(pq, __builtin_ctz(pq), ~0ULL, &r);
}

/* mode: random:p0 | strip:p0:frac | disk:p0:R   (outside region = bg) */
static void init_grid(Grid *g, const char *spec, int bg, const u32 key[2])
{
    char mode[16] = {0};
    double p0 = 0.5, par = 0;
    int n = sscanf(spec, "%15[^:]:%lf:%lf", mode, &p0, &par);
    if (n < 2) { fprintf(stderr, "bad --init %s\n", spec); exit(1); }
    int full = p0 >= 1.0;
    u32 pq = full ? 0 : (u32)llround(p0 * 4294967296.0);
    u64 bgw = bg ? ~0ULL : 0;
    for (int y = 0; y < g->Ly; y++)
        for (int x = 0; x < g->W; x++) {
            u64 w = full ? ~0ULL : random_word(key, y, x, pq), in;
            if (!strcmp(mode, "random")) in = ~0ULL;
            else if (!strcmp(mode, "strip")) in = y < par * g->Ly ? ~0ULL : 0;
            else if (!strcmp(mode, "disk")) {
                in = 0;
                double cy = g->Ly / 2.0, cx = g->Lx / 2.0;
                for (int b = 0; b < 64; b++) {
                    double dx = 64.0 * x + b + 0.5 - cx, dy = y + 0.5 - cy;
                    if (dx * dx + dy * dy <= par * par) in |= 1ULL << b;
                }
            } else { fprintf(stderr, "bad --init mode %s\n", mode); exit(1); }
            g->a[(size_t)y * g->W + x] = (w & in) | (bgw & ~in);
        }
}

/* init_grid's disk:p0:R on an all-dead lattice, touching only the words that meet the disk
 * (same cells, same random words).  Returns those words' box grown by one row and one word
 * on each side, which contains every cell with a disagreeing neighbour. */
static Box seed_disk(Grid *g, double p0, double R, const u32 key[2])
{
    const double cy = g->Ly / 2.0, cx = g->Lx / 2.0;
    int full = p0 >= 1.0;
    u32 pq = full ? 0 : (u32)llround(p0 * 4294967296.0);
    int y0 = (int)floor(cy - R - 1), y1 = (int)ceil(cy + R + 1);
    int x0 = (int)floor((cx - R - 1) / 64), x1 = (int)floor((cx + R + 1) / 64);
    if (y0 < 0) y0 = 0;
    if (y1 > g->Ly - 1) y1 = g->Ly - 1;
    if (x0 < 0) x0 = 0;
    if (x1 > g->W - 1) x1 = g->W - 1;
    for (int y = y0; y <= y1; y++)
        for (int x = x0; x <= x1; x++) {
            u64 in = 0;
            for (int b = 0; b < 64; b++) {
                double dx = 64.0 * x + b + 0.5 - cx, dy = y + 0.5 - cy;
                if (dx * dx + dy * dy <= R * R) in |= 1ULL << b;
            }
            u64 w = !in ? 0 : full ? ~0ULL : random_word(key, y, x, pq);
            g->a[(size_t)y * g->W + x] = w & in;
        }
    Box b = { y0 - 1, y1 + 1, x0 - 1, x1 + 1 };
    if (b.ya < 0) b.ya = 0;
    if (b.yb > g->Ly - 1) b.yb = g->Ly - 1;
    if (b.xa < 0) b.xa = 0;
    if (b.xb > g->W - 1) b.xb = g->W - 1;
    return b;
}

/* One replica of a spreading run: two buffers that agree outside the box (see Box). */
typedef struct { Grid g, nx; Box b; long long alive; } Rep;

enum { REP_LIVE = 0, REP_DEAD = 1, REP_ESCAPED = 2 };
#define REP_RESCAN 32

/* Seed an all-dead replica (both buffers) and fit the box.  Returns REP_DEAD if the seed
 * is empty. */
static int rep_seed(Rep *r, double p0, double R, const u32 key[2])
{
    seed_disk(&r->nx, p0, R, key);
    Box b = seed_disk(&r->g, p0, R, key);
    r->alive = (long long)measure_box(&r->g, b).alive;
    return !r->alive || !iface_bbox(&r->g, b, &r->b) ? REP_DEAD : REP_LIVE;
}

/* Advance the replica to time t.  REP_DEAD: every cell is dead.  REP_ESCAPED: the box came
 * within 2 rows or 1 word of the lattice edge, where it would have to wrap; the replica is
 * stopped there (the state at t is still exact). */
static int rep_step(Rep *r, u32 t, const Model *M)
{
    BoxStep s = step_box(&r->g, r->nx.a, t, M, r->b);
    u64 *tmp = r->g.a; r->g.a = r->nx.a; r->nx.a = tmp;
    r->alive += 2 * (long long)s.born - (long long)s.flips;
    if (!r->alive) return REP_DEAD;
    if (s.flips) {
        /* new disagreements appear only next to flipped cells */
        if (s.fy0 <= r->b.ya) r->b.ya = s.fy0 - 1;
        if (s.fy1 >= r->b.yb) r->b.yb = s.fy1 + 1;
        if (s.fx0 <= r->b.xa) r->b.xa = s.fx0 - 1;
        if (s.fx1 >= r->b.xb) r->b.xb = s.fx1 + 1;
        if (r->b.ya < 2 || r->b.yb > r->g.Ly - 3 || r->b.xa < 1 || r->b.xb > r->g.W - 2) return REP_ESCAPED;
    }
    if (t % REP_RESCAN == 0) {
        /* shrink: first make the buffers agree on the old box, then refit */
        const int W = r->g.W;
        for (int y = r->b.ya; y <= r->b.yb; y++)
            memcpy(r->nx.a + (size_t)y * W + r->b.xa, r->g.a + (size_t)y * W + r->b.xa,
                   8 * (size_t)(r->b.xb - r->b.xa + 1));
        Box nb;
        if (!iface_bbox(&r->g, r->b, &nb)) return REP_DEAD;   /* uniform; all dead, since alive < N */
        r->b = nb;
    }
    return REP_LIVE;
}

/* Zero both buffers on the box: the lattice is all dead again (outside the box it already is). */
static void rep_clear(Rep *r)
{
    const int W = r->g.W;
    for (int y = r->b.ya; y <= r->b.yb; y++) {
        memset(r->g.a + (size_t)y * W + r->b.xa, 0, 8 * (size_t)(r->b.xb - r->b.xa + 1));
        memset(r->nx.a + (size_t)y * W + r->b.xa, 0, 8 * (size_t)(r->b.xb - r->b.xa + 1));
    }
}

static const char MAGIC[8] = "DNSIM1\0";

static void save_grid(const Grid *g, u64 t, const char *path)
{
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); exit(1); }
    u64 hdr[3] = { (u64)g->Lx, (u64)g->Ly, t };
    fwrite(MAGIC, 1, 8, f); fwrite(hdr, 8, 3, f);
    fwrite(g->a, 8, (size_t)g->W * g->Ly, f);
    if (fclose(f)) { perror(path); exit(1); }
}

/* write path.tmp, then rename: a kill mid-write never leaves a torn checkpoint */
static void save_atomic(const Grid *g, u64 t, const char *path)
{
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    save_grid(g, t, tmp);
    if (rename(tmp, path)) { perror(path); exit(1); }
}

/* On resume from a checkpoint at t0: drop time-series lines with t > t0 (and the
 * "# done" line) that a killed run may have written after its last checkpoint. */
static void truncate_text(const char *path, u64 t0)
{
    FILE *f = fopen(path, "r");
    if (!f) return;
    char tmp[4096], line[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *o = fopen(tmp, "w");
    if (!o) { perror(tmp); exit(1); }
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#') {
            if (strncmp(line, "# done", 6)) fputs(line, o);
            continue;
        }
        if (strtoull(line, NULL, 10) > t0) break;
        fputs(line, o);
    }
    fclose(f);
    if (fclose(o) || rename(tmp, path)) { perror(path); exit(1); }
}

/* Same for the binary profile: keep whole records with t <= t0. */
static void truncate_profile(const char *path, int Ly, u64 t0)
{
    FILE *f = fopen(path, "rb");
    if (!f) return;
    long rec = 8 + 8L * Ly, keep = 0;
    fseek(f, 0, SEEK_END);
    long nrec = ftell(f) / rec;
    for (long i = 0; i < nrec; i++) {
        long long t;
        if (fseek(f, i * rec, SEEK_SET) || fread(&t, 8, 1, f) != 1 || (u64)t > t0) break;
        keep = i + 1;
    }
    fclose(f);
    if (truncate(path, keep * rec)) { perror(path); exit(1); }
}

static u64 load_grid(Grid *g, const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    char mg[8]; u64 hdr[3];
    if (fread(mg, 1, 8, f) != 8 || memcmp(mg, MAGIC, 8) || fread(hdr, 8, 3, f) != 3) {
        fprintf(stderr, "%s: not a dnsim checkpoint\n", path); exit(1);
    }
    if ((int)hdr[0] != g->Lx || (int)hdr[1] != g->Ly) {
        fprintf(stderr, "%s: size %llux%llu != %dx%d\n", path, (unsigned long long)hdr[0],
                (unsigned long long)hdr[1], g->Lx, g->Ly); exit(1);
    }
    if (fread(g->a, 8, (size_t)g->W * g->Ly, f) != (size_t)g->W * g->Ly) {
        fprintf(stderr, "%s: truncated\n", path); exit(1);
    }
    fclose(f);
    return hdr[2];
}

/* ------------------------------------------------------------ selftest */

static int naive_cell(const Grid *g, int X, int Y)
{
    X = (X % g->Lx + g->Lx) % g->Lx; Y = (Y % g->Ly + g->Ly) % g->Ly;
    return (g->a[(size_t)Y * g->W + X / 64] >> (X % 64)) & 1;
}

static int naive_k(const Grid *g, int X, int Y)
{
    int k = 0;
    for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++)
            if (dx || dy) k += naive_cell(g, X + dx, Y + dy);
    return k;
}

static int naive_k_vn(const Grid *g, int X, int Y)
{
    return naive_cell(g, X - 1, Y) + naive_cell(g, X + 1, Y) + naive_cell(g, X, Y - 1) + naive_cell(g, X, Y + 1);
}

static int selftest(void)
{
    int bad = 0;
    /* 1. Philox known-answer vectors (Random123 kat_vectors). */
    struct { u32 c[4], k[2], o[4]; } kat[3] = {
        {{0, 0, 0, 0}, {0, 0}, {0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8}},
        {{0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff}, {0xffffffff, 0xffffffff},
         {0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd}},
        {{0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344}, {0xa4093822, 0x299f31d0},
         {0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1}},
    };
    for (int i = 0; i < 3; i++) {
        u32 o[4];
        philox4x32_10(kat[i].c, kat[i].k, o);
        int ok = !memcmp(o, kat[i].o, sizeof o);
        printf("philox KAT %d: %s\n", i, ok ? "ok" : "FAIL");
        bad |= !ok;
    }
    /* 2. Bernoulli frequencies. */
    u32 key[2] = { 12345, 678 };
    double ps[] = { 0.5, 1.0 / 3, 0.02, 0.019, 0.999, 1e-4, 0.75 };
    for (size_t i = 0; i < sizeof ps / sizeof *ps; i++) {
        u32 pq = (u32)llround(ps[i] * 4294967296.0);
        u64 ones = 0, n = 0, bitcnt[64] = {0};
        const int NW = 1 << 20;
        for (int w = 0; w < NW; w++) {
            rstream r; rs_init(&r, key, 1, (u32)w, 7);
            u64 b = bern(pq, __builtin_ctz(pq), ~0ULL, &r);
            ones += __builtin_popcountll(b); n += 64;
            for (int j = 0; j < 64; j++) bitcnt[j] += (b >> j) & 1;
        }
        double p = pq / 4294967296.0, z = (ones - n * p) / sqrt(n * p * (1 - p)), zmax = 0;
        for (int j = 0; j < 64; j++) {
            double zj = fabs((bitcnt[j] - NW * p) / sqrt(NW * p * (1 - p)));
            if (zj > zmax) zmax = zj;
        }
        int ok = fabs(z) < 5 && zmax < 5.5;
        printf("bern p=%-8g freq=%.6f z=%+.2f max|z| over bit lanes=%.2f %s\n", p, (double)ones / n, z, zmax,
               ok ? "ok" : "FAIL");
        bad |= !ok;
    }
    /* 2b. bern_table on a voter-like 8-class table, random class per lane: per-class
     * frequency, and flip correlation between adjacent lanes (independent uniforms). */
    {
        Model M = {0};
        for (int k = 0; k < 9; k++)
            M.f[0][k] = M.f[1][8 - k] = 0.9 * k / 8 + 0.05 * (k == 3) + 0.0013 * (k == 2);
        model_quantise(&M);
        double cnt[MAXV] = {0}, fl[MAXV] = {0}, both = 0, eboth = 0, vboth = 0;
        const int NW = 1 << 21;
        for (int w = 0; w < NW; w++) {
            rstream rc, r;
            rs_init(&rc, key, 2, (u32)w, 0);
            u64 c0 = rs_next(&rc), c1 = rs_next(&rc), c2 = rs_next(&rc), m[MAXV];
            for (int v = 0; v < 8; v++)
                m[v] = (v & 1 ? c0 : ~c0) & (v & 2 ? c1 : ~c1) & (v & 4 ? c2 : ~c2);
            rs_init(&r, key, 3, (u32)w, 0);
            u64 res = bern_table(&M, m, ~0ULL, &r);
            double pl[64];
            for (int i = 0; i < 64; i++) {
                int v = (c0 >> i & 1) | (c1 >> i & 1) << 1 | (c2 >> i & 1) << 2;
                pl[i] = M.p[v] / 4294967296.0;
                cnt[v] += 1; fl[v] += res >> i & 1;
            }
            for (int i = 0; i < 63; i++) {
                double q = pl[i] * pl[i + 1];
                both += (res >> i & 1) & (res >> (i + 1) & 1); eboth += q; vboth += q * (1 - q);
            }
        }
        double zmax = 0;
        for (int v = 0; v < M.nv; v++) {
            double p = M.p[v] / 4294967296.0, z = (fl[v] - cnt[v] * p) / sqrt(cnt[v] * p * (1 - p));
            if (fabs(z) > zmax) zmax = fabs(z);
        }
        double zc = (both - eboth) / sqrt(vboth);
        int ok = M.nv == 8 && zmax < 5 && fabs(zc) < 5;
        printf("bern_table %d classes: max|z| per class=%.2f, adjacent-lane pair z=%+.2f %s\n", M.nv, zmax, zc,
               ok ? "ok" : "FAIL");
        bad |= !ok;
    }
    /* 3. Deterministic rules: bit-sliced step == naive per-cell step, exactly. */
    const char *rules[][2] = { {"3", "23"}, {"3678", "34678"}, {"36", "23"}, {"1357", "1357"}, {"0123", "8"} };
    int sizes[][2] = { {64, 64}, {128, 96}, {192, 7}, {64, 3} };
    for (size_t ri = 0; ri < sizeof rules / sizeof *rules; ri++)
        for (size_t si = 0; si < sizeof sizes / sizeof *sizes; si++) {
            Model M = {0}; M.key[0] = 1; M.key[1] = 2;
            for (int k = 0; k < 9; k++) {
                M.f[0][k] = strchr(rules[ri][0], '0' + k) != NULL;
                M.f[1][k] = strchr(rules[ri][1], '0' + k) == NULL;
            }
            model_quantise(&M);
            Grid g = grid_new(sizes[si][0], sizes[si][1]), h = grid_new(sizes[si][0], sizes[si][1]);
            int mism = 0;
            for (int trial = 0; trial < 20; trial++) {
                char spec[32]; snprintf(spec, sizeof spec, "random:%g", 0.1 + 0.04 * trial);
                u32 k2[2] = { (u32)trial, 99 };
                init_grid(&g, spec, 0, k2);
                step(&g, h.a, 1, &M, NULL);
                for (int Y = 0; Y < g.Ly; Y++)
                    for (int X = 0; X < g.Lx; X++) {
                        int s = naive_cell(&g, X, Y), k = naive_k(&g, X, Y);
                        int want = s ? strchr(rules[ri][1], '0' + k) != NULL : strchr(rules[ri][0], '0' + k) != NULL;
                        if (naive_cell(&h, X, Y) != want) mism++;
                    }
            }
            printf("deterministic B%s/S%s %dx%d: %s (%d mismatches)\n", rules[ri][0], rules[ri][1], g.Lx, g.Ly,
                   mism ? "FAIL" : "ok", mism);
            bad |= mism != 0;
            free(g.a); free(h.a);
        }
    /* 3b. Von Neumann neighbourhood: bit-sliced count == naive N+S+W+E count, exactly, for
     * deterministic 5-entry rules on the same lattice shapes. */
    {
        int vrules[][2] = { {0x02, 0x1c}, {0x18, 0x01}, {0x15, 0x0a}, {0x1f, 0x00}, {0x04, 0x1b} };
        int vmism = 0, vn_total = 0;
        for (size_t ri = 0; ri < sizeof vrules / sizeof *vrules; ri++)
            for (size_t si = 0; si < sizeof sizes / sizeof *sizes; si++) {
                Model M = {0}; M.key[0] = 1; M.key[1] = 2; M.vn = 1;
                for (int k = 0; k < 9; k++) {
                    M.f[0][k] = k <= 4 && ((vrules[ri][0] >> k) & 1);
                    M.f[1][k] = k <= 4 && !((vrules[ri][1] >> k) & 1);
                }
                model_quantise(&M);
                Grid g = grid_new(sizes[si][0], sizes[si][1]), h = grid_new(sizes[si][0], sizes[si][1]);
                for (int trial = 0; trial < 20; trial++) {
                    char spec[32]; snprintf(spec, sizeof spec, "random:%g", 0.1 + 0.04 * trial);
                    u32 k2[2] = { (u32)trial, 99 };
                    init_grid(&g, spec, 0, k2);
                    step(&g, h.a, 1, &M, NULL);
                    for (int Y = 0; Y < g.Ly; Y++)
                        for (int X = 0; X < g.Lx; X++) {
                            int s = naive_cell(&g, X, Y), k = naive_k_vn(&g, X, Y);
                            int want = s ? (vrules[ri][1] >> k) & 1 : (vrules[ri][0] >> k) & 1;
                            vmism += naive_cell(&h, X, Y) != want;
                            vn_total++;
                        }
                }
                free(g.a); free(h.a);
            }
        printf("deterministic von Neumann rules, 5 rules x 4 shapes (%d cells): %s (%d mismatches)\n", vn_total,
               vmism ? "FAIL" : "ok", vmism);
        bad |= vmism != 0;
    }
    /* 4. Stochastic table with 18 distinct classes: per-class flip frequency, exactness of 0/1. */
    {
        Model M = {0}; M.key[0] = 4; M.key[1] = 5;
        for (int s = 0; s < 2; s++)
            for (int k = 0; k < 9; k++) M.f[s][k] = 0.05 + 0.05 * (s * 9 + k);
        M.f[0][0] = 0; M.f[1][8] = 0; M.f[0][4] = 1; M.f[1][3] = 1;
        model_quantise(&M);
        Grid g = grid_new(128, 128), h = grid_new(128, 128);
        double cnt[2][9] = {{0}}, fl[2][9] = {{0}};
        u32 k2[2] = { 7, 7 };
        init_grid(&g, "random:0.5", 0, k2);
        for (u32 t = 1; t <= 400; t++) {
            step(&g, h.a, t, &M, NULL);
            for (int Y = 0; Y < g.Ly; Y++)
                for (int X = 0; X < g.Lx; X++) {
                    int s = naive_cell(&g, X, Y), k = naive_k(&g, X, Y);
                    cnt[s][k] += 1; fl[s][k] += naive_cell(&h, X, Y) != s;
                }
            /* re-randomise every 10 steps so all 18 classes stay populated */
            if (t % 10 == 0) { k2[0] = t; init_grid(&g, "random:0.5", 0, k2); }
            else { u64 *tmp = g.a; g.a = h.a; h.a = tmp; }
        }
        double zmax = 0;
        int exact_bad = 0;
        for (int s = 0; s < 2; s++)
            for (int k = 0; k < 9; k++) {
                double p = model_eff(&M, s, k);
                if (p == 0 || p == 1) { exact_bad += fl[s][k] != p * cnt[s][k]; continue; }
                double z = (fl[s][k] - cnt[s][k] * p) / sqrt(cnt[s][k] * p * (1 - p));
                if (fabs(z) > zmax) zmax = fabs(z);
            }
        int ok = zmax < 5 && !exact_bad;
        printf("stochastic 18-class table: max|z|=%.2f, exact-class violations=%d %s\n", zmax, exact_bad,
               ok ? "ok" : "FAIL");
        bad |= !ok;
        free(g.a); free(h.a);
    }
    /* 4b. Molofsky rule: table structure, and per-class flip frequencies against the naive
     * von Neumann count at a generic point (10 distinct classes), exactness of f(j=0) = 0. */
    {
        Model M = {0}; M.key[0] = 6; M.key[1] = 5;
        model_molofsky(&M, 0.2, 0.4);
        int tab_bad = 0;
        for (int k = 0; k <= 4; k++) tab_bad += fabs(M.f[0][k] - k / 5.0) > 1e-12 || fabs(M.f[1][k] - (4 - k) / 5.0) > 1e-12;
        model_molofsky(&M, 0.17, 0.53);
        for (int k = 0; k <= 4; k++) tab_bad += fabs(M.f[1][k] - M.f[0][4 - k]) > 1e-12;
        tab_bad += M.f[0][0] != 0 || M.f[1][4] != 0 || !M.vn;
        model_quantise(&M);
        Grid g = grid_new(128, 128), h = grid_new(128, 128);
        double cnt[2][5] = {{0}}, fl[2][5] = {{0}};
        u32 k2[2] = { 8, 8 };
        init_grid(&g, "random:0.5", 0, k2);
        for (u32 t = 1; t <= 400; t++) {
            step(&g, h.a, t, &M, NULL);
            for (int Y = 0; Y < g.Ly; Y++)
                for (int X = 0; X < g.Lx; X++) {
                    int s = naive_cell(&g, X, Y), k = naive_k_vn(&g, X, Y);
                    cnt[s][k] += 1; fl[s][k] += naive_cell(&h, X, Y) != s;
                }
            if (t % 10 == 0) { k2[0] = t; init_grid(&g, "random:0.5", 0, k2); }
            else { u64 *tmp = g.a; g.a = h.a; h.a = tmp; }
        }
        double zmax = 0;
        int exact_bad = 0;
        for (int s = 0; s < 2; s++)
            for (int k = 0; k <= 4; k++) {
                double p = model_eff(&M, s, k);
                if (p == 0 || p == 1) { exact_bad += fl[s][k] != p * cnt[s][k]; continue; }
                double z = (fl[s][k] - cnt[s][k] * p) / sqrt(cnt[s][k] * p * (1 - p));
                if (fabs(z) > zmax) zmax = fabs(z);
            }
        int ok = !tab_bad && M.nv == 4 && zmax < 5 && !exact_bad;
        printf("molofsky (0.17, 0.53) table + 4-value (10-class) von Neumann frequencies: %d classes, max|z|=%.2f, "
               "table errors=%d, exact-class violations=%d %s\n", M.nv, zmax, tab_bad, exact_bad, ok ? "ok" : "FAIL");
        bad |= !ok;
        free(g.a); free(h.a);
    }
    /* 5. Restart determinism: 0->40 in one go == 0->17, save/load, 17->40; thread count independence. */
    {
        Model M = {0}; M.key[0] = 9; M.key[1] = 9;
        for (int k = 0; k < 9; k++) { M.f[0][k] = (k == 3 || k >= 6) * 0.5; M.f[1][8 - k] = M.f[0][k]; }
        M.f[0][2] += 0.02; M.f[1][6] += 0.02;
        model_quantise(&M);
        Grid a = grid_new(256, 64), b = grid_new(256, 64), tmp = grid_new(256, 64);
        u32 k2[2] = { 3, 3 };
        init_grid(&a, "random:0.5", 0, k2);
        memcpy(b.a, a.a, (size_t)a.W * a.Ly * 8);
        for (u32 t = 1; t <= 40; t++) { step(&a, tmp.a, t, &M, NULL); u64 *x = a.a; a.a = tmp.a; tmp.a = x; }
#ifdef _OPENMP
        int nt = omp_get_max_threads();
        omp_set_num_threads(nt > 1 ? 1 : 2);
#endif
        for (u32 t = 1; t <= 17; t++) { step(&b, tmp.a, t, &M, NULL); u64 *x = b.a; b.a = tmp.a; tmp.a = x; }
        save_grid(&b, 17, "/tmp/dnsim_selftest.ckpt");
        memset(b.a, 0, (size_t)b.W * b.Ly * 8);
        u64 t0 = load_grid(&b, "/tmp/dnsim_selftest.ckpt");
        remove("/tmp/dnsim_selftest.ckpt");
        for (u32 t = (u32)t0 + 1; t <= 40; t++) { step(&b, tmp.a, t, &M, NULL); u64 *x = b.a; b.a = tmp.a; tmp.a = x; }
#ifdef _OPENMP
        omp_set_num_threads(nt);
#endif
        int ok = !memcmp(a.a, b.a, (size_t)a.W * a.Ly * 8);
        printf("restart + thread-count determinism: %s\n", ok ? "ok" : "FAIL");
        bad |= !ok;
        free(a.a); free(b.a); free(tmp.a);
    }
    /* 6. Box stepping (--reps) == full-lattice stepping, bit for bit, over 200 seeds of D&N +
     * voter + j=2 noise, on lattices small enough that some replicas escape; seed_disk ==
     * init_grid's disk; incremental alive count and measure_box == measure(). */
    {
        const int Bm = 1 << 3 | 1 << 6 | 1 << 7 | 1 << 8, Sm = 1 << 3 | 1 << 4 | 1 << 6 | 1 << 7 | 1 << 8;
        int mism = 0, n[3] = {0};
        for (int sd = 0; sd < 200; sd++) {
            Model M = {0}; M.key[0] = (u32)sd; M.key[1] = 77;
            double extra[9] = {0};
            extra[2] = (sd % 3) ? 0.02 : 0;
            model_from_rule(&M, Bm, Sm, 0.5, 0.9 + 0.05 * (sd % 3), extra);
            model_quantise(&M);
            int Lx = sd % 2 ? 256 : 512, Ly = sd % 2 ? 64 : 192;
            double p0 = sd % 4 ? 1.0 : 0.6, rad = sd % 4 ? 0.75 : 4 + sd % 7;
            Grid f = grid_new(Lx, Ly), tmp = grid_new(Lx, Ly);
            Rep r = { grid_new(Lx, Ly), grid_new(Lx, Ly), {0, 0, 0, 0}, 0 };
            memset(r.g.a, 0, (size_t)r.g.W * Ly * 8); memset(r.nx.a, 0, (size_t)r.g.W * Ly * 8);
            char spec[64]; snprintf(spec, sizeof spec, "disk:%.17g:%.17g", p0, rad);
            init_grid(&f, spec, 0, M.key);
            int st = rep_seed(&r, p0, rad, M.key);
            size_t bytes = (size_t)f.W * Ly * 8;
            mism += memcmp(f.a, r.g.a, bytes) != 0;
            for (u32 t = 1; t <= 600 && st == REP_LIVE; t++) {
                step(&f, tmp.a, t, &M, NULL);
                u64 *x = f.a; f.a = tmp.a; tmp.a = x;
                st = rep_step(&r, t, &M);
                if (st != REP_LIVE || t % 25 == 0) {
                    mism += memcmp(f.a, r.g.a, bytes) != 0;
                    Obs o = measure(&f, NULL, NULL);
                    if (st == REP_DEAD) { mism += o.alive != 0; continue; }
                    BoxObs bo = measure_box(&r.g, r.b);
                    mism += bo.alive != o.alive || bo.iface != o.iface || (long long)o.alive != r.alive;
                }
            }
            n[st]++;
            rep_clear(&r);
            for (size_t i = 0; i < bytes / 8; i++) mism += r.g.a[i] != 0 || r.nx.a[i] != 0;
            free(f.a); free(tmp.a); free(r.g.a); free(r.nx.a);
        }
        int ok = !mism && n[REP_DEAD] && n[REP_ESCAPED] && n[REP_LIVE];
        printf("box stepping == full stepping, 200 seeds (%d died, %d escaped, %d alive at t=600): %s (%d mismatches)\n",
               n[REP_DEAD], n[REP_ESCAPED], n[REP_LIVE], ok ? "ok" : "FAIL", mism);
        bad |= !ok;
    }
    /* 6b. Box stepping == full stepping for the Molofsky rule (von Neumann), 120 seeds. */
    {
        int mism = 0, n[3] = {0};
        for (int sd = 0; sd < 120; sd++) {
            Model M = {0}; M.key[0] = (u32)sd; M.key[1] = 78;
            model_molofsky(&M, 0.15 + 0.03 * (sd % 5), 0.3 + 0.1 * (sd % 3));
            model_quantise(&M);
            int Lx = sd % 2 ? 256 : 512, Ly = sd % 2 ? 64 : 192;
            double p0 = sd % 4 ? 1.0 : 0.5, rad = sd % 4 ? 0.75 : 3 + sd % 7;
            Grid f = grid_new(Lx, Ly), tmp = grid_new(Lx, Ly);
            Rep r = { grid_new(Lx, Ly), grid_new(Lx, Ly), {0, 0, 0, 0}, 0 };
            memset(r.g.a, 0, (size_t)r.g.W * Ly * 8); memset(r.nx.a, 0, (size_t)r.g.W * Ly * 8);
            char spec[64]; snprintf(spec, sizeof spec, "disk:%.17g:%.17g", p0, rad);
            init_grid(&f, spec, 0, M.key);
            int st = rep_seed(&r, p0, rad, M.key);
            size_t bytes = (size_t)f.W * Ly * 8;
            mism += memcmp(f.a, r.g.a, bytes) != 0;
            for (u32 t = 1; t <= 600 && st == REP_LIVE; t++) {
                step(&f, tmp.a, t, &M, NULL);
                u64 *x = f.a; f.a = tmp.a; tmp.a = x;
                st = rep_step(&r, t, &M);
                if (st != REP_LIVE || t % 25 == 0) {
                    mism += memcmp(f.a, r.g.a, bytes) != 0;
                    Obs o = measure(&f, NULL, NULL);
                    if (st == REP_DEAD) { mism += o.alive != 0; continue; }
                    BoxObs bo = measure_box(&r.g, r.b);
                    mism += bo.alive != o.alive || bo.iface != o.iface || (long long)o.alive != r.alive;
                }
            }
            n[st]++;
            rep_clear(&r);
            free(f.a); free(tmp.a); free(r.g.a); free(r.nx.a);
        }
        int ok = !mism && n[REP_DEAD] && n[REP_LIVE];
        printf("box stepping == full stepping, molofsky, 120 seeds (%d died, %d escaped, %d alive at t=600): %s (%d mismatches)\n",
               n[REP_DEAD], n[REP_ESCAPED], n[REP_LIVE], ok ? "ok" : "FAIL", mism);
        bad |= !ok;
    }
    printf(bad ? "SELFTEST FAILED\n" : "SELFTEST PASSED\n");
    return bad;
}

/* Header lines describing the model and the start: rule (or table), init, flip table. */
static void hdr_model(FILE *f, const char *molofsky, const char *ftab, const char *rule, double alpha, double voter,
                      const char *jflip, const Model *M, const char *init, int bg, const char *load, u64 t0)
{
    if (molofsky) fprintf(f, "# molofsky p1,p2 %s (five-cell rule, von Neumann + self, synchronous)\n", molofsky);
    else if (ftab) fprintf(f, "# ftab %s\n", ftab);
    else fprintf(f, "# rule %s alpha %.17g voter %.17g jflip %s\n", rule, alpha, voter, jflip ? jflip : "-");
    fprintf(f, "# nbhd %s\n", M->vn ? "vn" : "moore");
    fprintf(f, "# init %s bg %d load %s t0 %llu\n", init, bg, load, (unsigned long long)t0);
    int kmax = model_kmax(M);
    for (int s = 0; s < 2; s++) {
        fprintf(f, "# flip s=%d k=0..%d:", s, kmax);
        for (int k = 0; k <= kmax; k++) fprintf(f, " %.10g", model_eff(M, s, k));
        fprintf(f, "\n");
    }
}

/* --------------------------------------------------------------- --reps */

/* Spreading from a seed (--init disk:p0:R on an all-dead torus): many independent replicas,
 * each stepped on its bounding box (Rep).  Replica i uses key (seed lo, seed hi ^ 0x9E3779B9
 * (i+1)) and adds into group i % G, so errors can be taken by jackknife over groups.  Per
 * group and output time: replicas still alive (surv, escaped ones included), escaped so far
 * (esc), and sums over replicas of the interface count I, I^2, live count A, A^2, the sum of
 * r^2 over interface cells, and the box area.  Dead replicas add zeros; escaped replicas add
 * only to surv and esc after their escape.  Runs until `reps` replicas or `target` replicas
 * alive (or escaped) at T.  Resumable between replicas, also with a larger reps or target to
 * extend a finished run; the results do not depend on where it was interrupted. */

enum { C_SURV, C_ESC, C_I, C_I2, C_A, C_A2, C_IR2, C_AREA, NCOL };
static const char REPS_MAGIC[8] = "DNREPS1";

typedef struct {
    u64 Lx, Ly, T, nt, G, seed;                 /* must match on resume */
    u64 reps, target;                           /* may grow: a resumed run then continues */
    u64 next, surv_T, escaped, done;
    double updates, wall;
} RepsHead;

static void reps_acc(double *a, BoxObs o, Box b)
{
    a[C_SURV] += o.iface > 0;
    a[C_I] += (double)o.iface; a[C_I2] += (double)o.iface * o.iface;
    a[C_A] += (double)o.alive; a[C_A2] += (double)o.alive * o.alive;
    a[C_IR2] += o.r2;
    a[C_AREA] += 64.0 * (b.yb - b.ya + 1) * (b.xb - b.xa + 1);
}

static void reps_save(const char *path, const RepsHead *h, const u64 *ng, const double *acc)
{
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "wb");
    if (!f) { perror(tmp); exit(1); }
    fwrite(REPS_MAGIC, 1, 8, f); fwrite(h, sizeof *h, 1, f);
    fwrite(ng, 8, h->G, f); fwrite(acc, 8, h->G * h->nt * NCOL, f);
    if (fclose(f) || rename(tmp, path)) { perror(path); exit(1); }
}

static int reps_load(const char *path, RepsHead *h, u64 *ng, double *acc)
{
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); return 0; }
    char mg[8]; RepsHead r;
    int ok = fread(mg, 1, 8, f) == 8 && !memcmp(mg, REPS_MAGIC, 8) && fread(&r, sizeof r, 1, f) == 1 &&
             !memcmp(&r, h, 6 * sizeof(u64)) && fread(ng, 8, r.G, f) == r.G &&
             fread(acc, 8, r.G * r.nt * NCOL, f) == r.G * r.nt * NCOL;
    fclose(f);
    if (!ok) { fprintf(stderr, "%s: not a checkpoint of this run\n", path); return 0; }
    r.reps = h->reps; r.target = h->target; r.done = 0;   /* the stopping rule is re-checked */
    *h = r;
    return 1;
}

static void reps_write(const char *out, const char *hdr, const RepsHead *h, const long long *tl, const u64 *ng,
                       const double *acc)
{
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", out);
    FILE *f = fopen(tmp, "w");
    if (!f) { perror(tmp); exit(1); }
    fputs(hdr, f);
    fprintf(f, "# replicas %llu survivors_T %llu escaped %llu cell_updates %.6g wall_s %.1f\n",
            (unsigned long long)h->next, (unsigned long long)h->surv_T, (unsigned long long)h->escaped,
            h->updates, h->wall);
    fprintf(f, "# columns: group t n surv esc sumI sumI2 sumA sumA2 sumIR2 sumArea\n");
    for (u64 g = 0; g < h->G; g++)
        for (u64 k = 0; k < h->nt; k++) {
            const double *a = acc + (g * h->nt + k) * NCOL;
            fprintf(f, "%llu %lld %llu", (unsigned long long)g, tl[k], (unsigned long long)ng[g]);
            for (int c = 0; c < NCOL; c++) fprintf(f, " %.17g", a[c]);
            fputc('\n', f);
        }
    if (h->done)
        fprintf(f, "# done replicas %llu survivors_T %llu escaped %llu wall_s %.3f cell_updates_per_s %.4g\n",
                (unsigned long long)h->next, (unsigned long long)h->surv_T, (unsigned long long)h->escaped,
                h->wall, h->updates / (h->wall > 0 ? h->wall : 1e-9));
    if (fclose(f) || rename(tmp, out)) { perror(out); exit(1); }
}

static int run_reps(const Model *M, int Lx, int Ly, long long T, const char *init, int ppd, long long reps,
                    long long target, int G, const char *out, const char *ckpt, int resume, long long ckpt_reps,
                    unsigned long long seed, const char *hdr)
{
    char mode[16] = {0};
    double p0 = 1, R = 0;
    if (sscanf(init, "%15[^:]:%lf:%lf", mode, &p0, &R) != 3 || strcmp(mode, "disk") || R <= 0) {
        fprintf(stderr, "--reps needs --init disk:p0:R\n"); return 2;
    }
    long long nt = 1, cap = 1024, *tl = malloc(cap * sizeof *tl);
    tl[0] = 0;
    for (int d = 0; d <= 10; d++)
        for (int i = 0; i < ppd; i++) {
            long long t = llround(pow(10.0, d + (double)i / ppd));
            if (t <= T && t > tl[nt - 1]) {
                if (nt == cap) tl = realloc(tl, (cap *= 2) * sizeof *tl);
                tl[nt++] = t;
            }
        }
    if (tl[nt - 1] < T) {
        if (nt == cap) tl = realloc(tl, (cap *= 2) * sizeof *tl);
        tl[nt++] = T;
    }
    RepsHead h = { (u64)Lx, (u64)Ly, (u64)T, (u64)nt, (u64)G, seed, (u64)reps, (u64)target, 0, 0, 0, 0, 0, 0 };
    u64 *ng = calloc(G, sizeof *ng);
    double *acc = calloc((size_t)G * nt * NCOL, sizeof *acc);
    if (resume && access(ckpt, R_OK) == 0 && !reps_load(ckpt, &h, ng, acc)) return 2;

    Rep r = { grid_new(Lx, Ly), grid_new(Lx, Ly), {0, 0, 0, 0}, 0 };
    memset(r.g.a, 0, (size_t)r.g.W * Ly * 8);
    memset(r.nx.a, 0, (size_t)r.g.W * Ly * 8);
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    double wall0 = h.wall, last = 0;
    long long since = 0;
    while (!h.done) {
        if (h.next >= (u64)reps || (target && h.surv_T >= (u64)target)) h.done = 1;
        else {
            u64 i = h.next;
            Model Mr = *M;
            Mr.key[1] ^= 0x9E3779B9u * (u32)(i + 1);
            double *A = acc + (i % G) * nt * NCOL;
            int st = rep_seed(&r, p0, R, Mr.key);
            long long ko = 1;
            if (st == REP_LIVE) reps_acc(A, measure_box(&r.g, r.b), r.b);
            for (long long t = 1; t <= T && st == REP_LIVE; t++) {
                h.updates += 64.0 * (r.b.yb - r.b.ya + 1) * (r.b.xb - r.b.xa + 1);
                st = rep_step(&r, (u32)t, &Mr);
                if (st == REP_LIVE && tl[ko] == t) reps_acc(A + ko++ * NCOL, measure_box(&r.g, r.b), r.b);
            }
            if (st == REP_ESCAPED) {
                for (long long k = ko; k < nt; k++) { A[k * NCOL + C_SURV] += 1; A[k * NCOL + C_ESC] += 1; }
                h.escaped++;
            }
            if (st != REP_DEAD) h.surv_T++;
            rep_clear(&r);
            ng[i % G]++;
            h.next++;
            since++;
        }
        clock_gettime(CLOCK_MONOTONIC, &ts1);
        double el = (ts1.tv_sec - ts0.tv_sec) + 1e-9 * (ts1.tv_nsec - ts0.tv_nsec);
        h.wall = wall0 + el;
        if (h.done || since >= ckpt_reps || el - last >= 600) {
            if (ckpt) reps_save(ckpt, &h, ng, acc);
            reps_write(out, hdr, &h, tl, ng, acc);
            since = 0; last = el;
        }
    }
    free(r.g.a); free(r.nx.a); free(tl); free(ng); free(acc);
    return 0;
}

/* ----------------------------------------------------------------- main */

static void usage(void)
{
    fprintf(stderr,
        "usage: dnsim [options]\n"
        "  --L n | --Lx n --Ly n   lattice (Lx multiple of 64), default 1024\n"
        "  --T n                   final time (absolute), default 1000\n"
        "  --seed n                RNG key, default 1\n"
        "  --rule Bxxx/Syyy        deterministic part, default B3678/S34678\n"
        "  --alpha a               apply the rule's flips with prob a, default 1\n"
        "  --voter p               mix in voter steps with prob p, default 0\n"
        "  --jflip j:q[,j:q...]    extra flip prob q for cells with j disagreeing nbrs\n"
        "  --ftab f00,..,f08;f10,..,f18   explicit flip table (overrides the above)\n"
        "  --molofsky p1,p2        five-cell rule of Molofsky et al. 1999 (von Neumann + self):\n"
        "                          next state 1 with prob p_n, n = #1s among the 5 cells, p_0 = 0,\n"
        "                          p_5 = 1, p_3 = 1-p_2, p_4 = 1-p_1; voter point (0.2, 0.4)\n"
        "  --nbhd moore|vn         neighbourhood for --rule/--ftab tables (k = live nbrs, 0..8 or\n"
        "                          0..4; vn ignores entries k >= 5), default moore\n"
        "  --init spec             random:p0 | strip:p0:frac | disk:p0:R, default random:0.5\n"
        "  --bg 0|1                background outside strip/disk, default 0\n"
        "  --load file             start from checkpoint (continues its time)\n"
        "  --save file             write checkpoint at the end\n"
        "  --ppd n                 log-spaced outputs per decade, default 20\n"
        "  --every n               also output every n steps\n"
        "  --out file              time series (default stdout)\n"
        "  --profile file          binary per-row alive/iface counts at output times\n"
        "  --stop-absorbed         stop when no cell has a disagreeing neighbour\n"
        "  --ckpt file             checkpoint path, written atomically at the end and\n"
        "  --ckpt-every n          ... every n steps\n"
        "  --resume                if the --ckpt file exists, continue from it: outputs\n"
        "                          are truncated to t <= its time and appended to\n"
        "  --reps n                spreading mode: up to n replicas of --init disk:p0:R\n"
        "                          on an all-dead lattice, each stepped on the bounding box of\n"
        "                          its activity; --out gets per-group sums (see run_reps)\n"
        "  --target n              ... stop once n replicas are alive (or escaped) at T\n"
        "  --groups n              ... replica i adds into group i %% n, default 16\n"
        "                          (--ckpt-every then counts replicas, default 1000)\n"
        "  --selftest\n");
    exit(2);
}

int main(int argc, char **argv)
{
    int Lx = 1024, Ly = 1024, bg = 0, ppd = 20, stop_abs = 0;
    long long T = 1000, every = 0;
    unsigned long long seed = 1;
    double alpha = 1, voter = 0;
    const char *rule = "B3678/S34678", *jflip = NULL, *ftab = NULL, *init = "random:0.5";
    const char *molofsky = NULL, *nbhd = "moore";
    const char *load = NULL, *save = NULL, *out = NULL, *prof = NULL, *ckpt = NULL;
    long long ckpt_every = 0;
    int resume = 0, groups = 16;
    long long reps = 0, target = 0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
#define ARG() (i + 1 < argc ? argv[++i] : (usage(), ""))
        if (!strcmp(a, "--selftest")) return selftest();
        else if (!strcmp(a, "--L")) Lx = Ly = atoi(ARG());
        else if (!strcmp(a, "--Lx")) Lx = atoi(ARG());
        else if (!strcmp(a, "--Ly")) Ly = atoi(ARG());
        else if (!strcmp(a, "--T")) T = atoll(ARG());
        else if (!strcmp(a, "--seed")) seed = strtoull(ARG(), NULL, 0);
        else if (!strcmp(a, "--rule")) rule = ARG();
        else if (!strcmp(a, "--alpha")) alpha = atof(ARG());
        else if (!strcmp(a, "--voter")) voter = atof(ARG());
        else if (!strcmp(a, "--jflip")) jflip = ARG();
        else if (!strcmp(a, "--ftab")) ftab = ARG();
        else if (!strcmp(a, "--molofsky")) molofsky = ARG();
        else if (!strcmp(a, "--nbhd")) nbhd = ARG();
        else if (!strcmp(a, "--init")) init = ARG();
        else if (!strcmp(a, "--bg")) bg = atoi(ARG());
        else if (!strcmp(a, "--load")) load = ARG();
        else if (!strcmp(a, "--save")) save = ARG();
        else if (!strcmp(a, "--ppd")) ppd = atoi(ARG());
        else if (!strcmp(a, "--every")) every = atoll(ARG());
        else if (!strcmp(a, "--out")) out = ARG();
        else if (!strcmp(a, "--profile")) prof = ARG();
        else if (!strcmp(a, "--stop-absorbed")) stop_abs = 1;
        else if (!strcmp(a, "--ckpt")) ckpt = ARG();
        else if (!strcmp(a, "--ckpt-every")) ckpt_every = atoll(ARG());
        else if (!strcmp(a, "--resume")) resume = 1;
        else if (!strcmp(a, "--reps")) reps = atoll(ARG());
        else if (!strcmp(a, "--target")) target = atoll(ARG());
        else if (!strcmp(a, "--groups")) groups = atoi(ARG());
        else usage();
    }
    if (Lx < 64 || Lx % 64 || Ly < 3 || Lx / 64 >= (1 << 24) || T < 1 || T > 4294967295LL) {
        fprintf(stderr, "bad sizes: Lx must be a positive multiple of 64, Ly >= 3, 1 <= T < 2^32\n");
        return 2;
    }

    Model M = {0};
    M.key[0] = (u32)seed; M.key[1] = (u32)(seed >> 32);
    if (!strcmp(nbhd, "vn")) M.vn = 1;
    else if (strcmp(nbhd, "moore")) { fprintf(stderr, "bad --nbhd\n"); return 2; }
    if (molofsky) {
        double p1, p2;
        if (sscanf(molofsky, "%lf,%lf", &p1, &p2) != 2 || p1 < 0 || p1 > 1 || p2 < 0 || p2 > 1) {
            fprintf(stderr, "bad --molofsky (want p1,p2 in [0,1])\n"); return 2;
        }
        model_molofsky(&M, p1, p2);
    } else if (ftab) {
        if (sscanf(ftab, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf;%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &M.f[0][0], &M.f[0][1], &M.f[0][2], &M.f[0][3], &M.f[0][4], &M.f[0][5], &M.f[0][6],
                   &M.f[0][7], &M.f[0][8], &M.f[1][0], &M.f[1][1], &M.f[1][2], &M.f[1][3], &M.f[1][4],
                   &M.f[1][5], &M.f[1][6], &M.f[1][7], &M.f[1][8]) != 18) {
            fprintf(stderr, "bad --ftab\n"); return 2;
        }
    } else {
        int Bm = 0, Sm = 0;
        const char *p = rule;
        if (*p != 'B' && *p != 'b') { fprintf(stderr, "bad --rule\n"); return 2; }
        for (p++; *p >= '0' && *p <= '8'; p++) Bm |= 1 << (*p - '0');
        if (*p++ != '/' || (*p != 'S' && *p != 's')) { fprintf(stderr, "bad --rule\n"); return 2; }
        for (p++; *p >= '0' && *p <= '8'; p++) Sm |= 1 << (*p - '0');
        double extra[9] = {0};
        if (jflip) {
            const char *q = jflip;
            while (*q) {
                int j; double v; int n;
                if (sscanf(q, "%d:%lf%n", &j, &v, &n) != 2 || j < 0 || j > 8) { fprintf(stderr, "bad --jflip\n"); return 2; }
                extra[j] += v;
                q += n;
                if (*q == ',') q++;
            }
        }
        model_from_rule(&M, Bm, Sm, alpha, voter, extra);
    }
    model_quantise(&M);
    int absorbing = M.slot[0][0] == SLOT_NEVER && M.slot[1][model_kmax(&M)] == SLOT_NEVER;

    if ((resume || ckpt_every) && !ckpt) { fprintf(stderr, "--resume/--ckpt-every need --ckpt\n"); return 2; }
    int resumed = resume && access(ckpt, R_OK) == 0;

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    if (reps) {
        if (!absorbing || bg || load || save || prof || !out || groups < 1) {
            fprintf(stderr, "--reps needs an absorbing model, --bg 0, --out, and no --load/--save/--profile\n");
            return 2;
        }
        char *hdr = NULL;
        size_t hlen = 0;
        FILE *hs = open_memstream(&hdr, &hlen);
        fprintf(hs, "# dnsim v1.4 reps\n# Lx %d Ly %d T %lld seed %llu threads %d\n", Lx, Ly, T, seed, nthreads);
        hdr_model(hs, molofsky, ftab, rule, alpha, voter, jflip, &M, init, bg, "-", 0);
        fprintf(hs, "# reps_max %lld target %lld groups %d rescan %d escape: box within 2 rows / 1 word of the edge\n",
                reps, target, groups, REP_RESCAN);
        fprintf(hs, "# replica i: key (seed lo, seed hi ^ 0x9E3779B9*(i+1)), group i %% groups\n");
        fclose(hs);
        int rc = run_reps(&M, Lx, Ly, T, init, ppd, reps, target, groups, out, ckpt, resume,
                          ckpt_every ? ckpt_every : 1000, seed, hdr);
        free(hdr);
        return rc;
    }

    Grid g = grid_new(Lx, Ly), nx = grid_new(Lx, Ly);
    u64 t0 = 0;
    if (resumed) t0 = load_grid(&g, ckpt);
    else if (load) t0 = load_grid(&g, load);
    else init_grid(&g, init, bg, M.key);
    if ((long long)t0 >= T) {
        fprintf(stderr, "checkpoint time %llu >= T\n", (unsigned long long)t0);
        return resumed ? 0 : 2;   /* resuming a finished run is a no-op */
    }
    if (resumed) {
        if (out) truncate_text(out, t0);
        if (prof) truncate_profile(prof, Ly, t0);
    }

    FILE *fo = out ? fopen(out, resumed ? "a" : "w") : stdout;
    if (!fo) { perror(out); return 1; }
    FILE *fp = prof ? fopen(prof, resumed ? "ab" : "wb") : NULL;
    if (prof && !fp) { perror(prof); return 1; }

    if (resumed) fprintf(fo, "# resumed t0 %llu threads %d\n", (unsigned long long)t0, nthreads);
    else {
    fprintf(fo, "# dnsim v1.4\n# Lx %d Ly %d T %lld seed %llu threads %d\n", Lx, Ly, T, seed, nthreads);
    hdr_model(fo, molofsky, ftab, rule, alpha, voter, jflip, &M, load ? "-" : init, bg, load ? load : "-", t0);
    fprintf(fo, "# absorbing_homogeneous %d\n# columns: t alive_frac rho_bond flip_frac iface_frac\n", absorbing);
    }

    /* output schedule */
    long long nlog = 0, cap = 1024, *tl = malloc(cap * sizeof *tl);
    for (int d = 0; d <= 10; d++)
        for (int i = 0; i < ppd; i++) {
            long long t = llround(pow(10.0, d + (double)i / ppd));
            if (t > (long long)t0 && t <= T && (!nlog || t > tl[nlog - 1])) {
                if (nlog == cap) tl = realloc(tl, (cap *= 2) * sizeof *tl);
                tl[nlog++] = t;
            }
        }
    long long li = 0;

    const double N = (double)Lx * Ly;
    u32 *ra = fp ? malloc(sizeof(u32) * Ly) : NULL, *ri = fp ? malloc(sizeof(u32) * Ly) : NULL;
    Obs o;
    if (!resumed) {
        o = measure(&g, ra, ri);
        fprintf(fo, "%llu %.9f %.9f %.9f %.9f\n", (unsigned long long)t0, o.alive / N, o.bonds / (2 * N), 0.0, o.iface / N);
        if (fp) { long long tt = (long long)t0; fwrite(&tt, 8, 1, fp); fwrite(ra, 4, Ly, fp); fwrite(ri, 4, Ly, fp); }
        fflush(fo);
    }

    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    long long t;
    int absorbed = 0;
    for (t = (long long)t0 + 1; t <= T; t++) {
        u64 fl = step(&g, nx.a, (u32)t, &M, NULL);
        u64 *tmp = g.a; g.a = nx.a; nx.a = tmp;
        while (li < nlog && tl[li] < t) li++;
        int due = (li < nlog && tl[li] == t) || (every && t % every == 0) || t == T;
        if (due || (stop_abs && absorbing && fl == 0)) {
            o = measure(&g, ra, ri);
            if (due || o.iface == 0) {
                fprintf(fo, "%lld %.9f %.9f %.9f %.9f\n", t, o.alive / N, o.bonds / (2 * N), fl / N, o.iface / N);
                if (fp) { fwrite(&t, 8, 1, fp); fwrite(ra, 4, Ly, fp); fwrite(ri, 4, Ly, fp); }
                fflush(fo);
            }
            if (stop_abs && absorbing && o.iface == 0) { absorbed = 1; break; }
        }
        if (ckpt_every && t % ckpt_every == 0 && t < T) {
            /* outputs up to t must be on disk before the checkpoint that claims them */
            fflush(fo);
            if (fp) fflush(fp);
            save_atomic(&g, (u64)t, ckpt);
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    long long tend = absorbed ? t : T;
    double wall = (ts1.tv_sec - ts0.tv_sec) + 1e-9 * (ts1.tv_nsec - ts0.tv_nsec);
    double steps = (double)(tend - (long long)t0);
    fprintf(fo, "# done t_end %lld absorbed %d wall_s %.3f cell_updates_per_s %.4g\n", tend, absorbed, wall,
            steps * N / (wall > 0 ? wall : 1e-9));
    if (save) save_grid(&g, (u64)tend, save);
    if (fp) fflush(fp);
    fflush(fo);
    if (ckpt) save_atomic(&g, (u64)tend, ckpt);
    if (fo != stdout) fclose(fo);
    if (fp) fclose(fp);
    free(g.a); free(nx.a); free(tl); free(ra); free(ri);
    return 0;
}
