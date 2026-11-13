/* cwalk: coalescing random walks dual to the synchronous Moore voter model.
 *
 * Nine walkers start at the origin (walker 0) and its 8 Moore neighbours (walkers 1..8,
 * in the order of DX/DY). Every step each surviving lineage jumps to a uniform Moore
 * neighbour (synchronous voter: every site copies a uniform neighbour), and lineages
 * on the same site merge. A trial is dropped once <= 2 lineages remain: with a
 * complement-symmetric perturbation, 1- and 2-class patterns contribute exactly zero to
 * the drift (see vpt.py), and lineages never split.
 *
 * Lazy dual (env CWALK_R = r, default 1): each lineage jumps with probability r per step and
 * otherwise stays, the dual of the voter model in which each site updates with probability r
 * (dnsim with the flip table scaled by r).  r = 1 draws no extra bits, so output is unchanged.
 *
 * Output, at each checkpoint T: the histogram of patterns with >= 3 classes, where a
 * pattern is the sorted list of how many of the 8 neighbour sites each NON-origin class
 * holds (neighbours in the origin's class always agree with it, so they are not listed):
 *     T <T> trials <N>
 *     <count> <c1,c2,...>
 *
 *     cc -O3 -fopenmp -o cwalk cwalk.c
 *     ./cwalk <trials> <seed> T1 T2 ...      (checkpoints ascending)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int DX[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static const int DY[8] = {-1, -1, -1, 0, 0, 1, 1, 1};

#define MAXCK 64
#define MAXPAT 512

typedef struct { uint64_t code[MAXPAT]; long long n[MAXPAT]; int len; } Hist;

static uint64_t splitmix(uint64_t *s)
{
    uint64_t z = (*s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void hist_add(Hist *h, uint64_t code, long long n)
{
    for (int i = 0; i < h->len; i++)
        if (h->code[i] == code) { h->n[i] += n; return; }
    if (h->len == MAXPAT) { fprintf(stderr, "pattern table full\n"); exit(1); }
    h->code[h->len] = code; h->n[h->len++] = n;
}

/* sorted (descending) neighbour counts of non-origin classes, as decimal digits */
static uint64_t pattern(const int *cls, const int *roots, int k)
{
    int c[9], m = 0;
    for (int r = 0; r < k; r++) {
        if (roots[r] == cls[0]) continue;
        int n = 0;
        for (int i = 1; i < 9; i++) n += cls[i] == roots[r];
        c[m++] = n;
    }
    for (int i = 1; i < m; i++)
        for (int j = i; j > 0 && c[j] > c[j - 1]; j--) { int t = c[j]; c[j] = c[j - 1]; c[j - 1] = t; }
    uint64_t code = 0;
    for (int i = 0; i < m; i++) code = code * 10 + c[i];
    return code;
}

int main(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "usage: cwalk trials seed T1 [T2 ...]\n"); return 2; }
    long long trials = atoll(argv[1]);
    uint64_t seed = strtoull(argv[2], NULL, 10);
    int nck = argc - 3;
    if (nck > MAXCK) { fprintf(stderr, "too many checkpoints\n"); return 2; }
    long long ck[MAXCK];
    for (int i = 0; i < nck; i++) ck[i] = atoll(argv[3 + i]);
    Hist *H = calloc(nck, sizeof(Hist));
    const char *lr = getenv("CWALK_R");
    double lazy = lr ? atof(lr) : 1.0;
    if (!(lazy > 0 && lazy <= 1)) { fprintf(stderr, "CWALK_R must be in (0, 1]\n"); return 2; }
    uint32_t thr16 = (uint32_t)(lazy * 65536.0 + 0.5);   /* jump iff 16-bit coin < thr16 */
    fprintf(stderr, "cwalk: CWALK_R %g (coin threshold %u / 65536)\n", lazy, thr16);

#pragma omp parallel
    {
        Hist *h = calloc(nck, sizeof(Hist));
#pragma omp for schedule(dynamic, 1024)
        for (long long tr = 0; tr < trials; tr++) {
            uint64_t s = seed * 0xd1342543de82ef95ULL + (uint64_t)tr;
            splitmix(&s);
            int px[9], py[9], cls[9], roots[9], k = 9;
            px[0] = py[0] = 0;
            for (int i = 0; i < 8; i++) { px[i + 1] = DX[i]; py[i + 1] = DY[i]; }
            for (int i = 0; i < 9; i++) cls[i] = roots[i] = i;   /* roots index walkers */
            int ci = 0;
            uint64_t bits = 0, cbits = 0; int nb = 0, cnb = 0;
            for (long long t = 0; ci < nck; t++) {
                if (t > 0) {
                    for (int r = 0; r < k; r++) {
                        if (thr16 < 65536) {
                            if (cnb < 16) { cbits = splitmix(&s); cnb = 64; }
                            int move = (cbits & 0xffff) < thr16;
                            cbits >>= 16; cnb -= 16;
                            if (!move) continue;
                        }
                        if (nb < 3) { bits = splitmix(&s); nb = 64; }
                        int d = bits & 7; bits >>= 3; nb -= 3;
                        px[roots[r]] += DX[d]; py[roots[r]] += DY[d];
                    }
                    for (int a = 0; a < k; a++)
                        for (int b = a + 1; b < k; b++) {
                            int ra = roots[a], rb = roots[b];
                            if (px[ra] == px[rb] && py[ra] == py[rb]) {
                                for (int i = 0; i < 9; i++) if (cls[i] == rb) cls[i] = ra;
                                roots[b] = roots[--k];
                                b--;
                            }
                        }
                    if (k <= 2) break;
                }
                while (ci < nck && ck[ci] == t) hist_add(&h[ci++], pattern(cls, roots, k), 1);
            }
        }
#pragma omp critical
        for (int c = 0; c < nck; c++)
            for (int i = 0; i < h[c].len; i++) hist_add(&H[c], h[c].code[i], h[c].n[i]);
        free(h);
    }

    for (int c = 0; c < nck; c++) {
        printf("T %lld trials %lld\n", ck[c], trials);
        for (int i = 0; i < H[c].len; i++) {
            char buf[32]; int m = 0;
            uint64_t code = H[c].code[i];
            char dig[20]; int nd = 0;
            while (code) { dig[nd++] = '0' + code % 10; code /= 10; }
            for (int j = nd - 1; j >= 0; j--) { buf[m++] = dig[j]; if (j) buf[m++] = ','; }
            buf[m] = 0;
            printf("%lld %s\n", H[c].n[i], buf);
        }
    }
    free(H);
    return 0;
}
