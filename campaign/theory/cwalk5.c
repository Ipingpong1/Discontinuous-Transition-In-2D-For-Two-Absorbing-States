/* cwalk5: coalescing lazy random walks dual to the synchronous five-cell voter model of
 * Molofsky, Durrett, Dushoff, Griffeath & Levin, Theor. Popul. Biol. 55, 270 (1999).
 *
 * Their automaton updates every site at once; at the voter point (p_1, p_2) = (1/5, 2/5) the
 * probability that a site becomes 1 is n/5, with n the number of 1s among the five cells of
 * its von Neumann neighbourhood INCLUDING ITSELF.  That is exactly "copy a uniformly chosen
 * cell among the five": the site's new state is the old state of x + e, e uniform in
 * {0, +x, -x, +y, -y}.  Tracing a site's state backward therefore gives a lazy walk that
 * stays with probability 1/5 and steps to each of its four neighbours with probability 1/5;
 * the walks of several sites use independent choices while apart and the same choice once
 * they sit on the same site, so they coalesce on contact and never split.  This is the
 * synchronous analogue of the continuous-time dual used by Cox, Durrett & Perkins.
 *
 * Five walkers start at the origin (walker 0) and its 4 nearest neighbours (walkers 1..4,
 * in the order of DX/DY).  A trial is dropped once <= 2 lineages remain: with a
 * complement-symmetric perturbation, 1- and 2-class patterns contribute exactly zero to
 * the drift (see vpt.py), and lineages never split.
 *
 * Output, at each checkpoint T: the histogram of patterns with >= 3 classes, where a
 * pattern is the sorted list of how many of the 4 neighbour sites each NON-origin class
 * holds (neighbours in the origin's class always agree with it, so they are not listed).
 * Same format as cwalk.c, read by vpt.read_patterns:
 *     T <T> trials <N>
 *     <count> <c1,c2,...>
 *
 *     cc -O3 -fopenmp -o cwalk5 cwalk5.c
 *     ./cwalk5 <trials> <seed> T1 T2 ...      (checkpoints ascending)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NW 4                 /* neighbours */
#define NWALK (NW + 1)       /* lineages: origin + neighbours */
#define NMOVE 5              /* stay + 4 steps, each with probability 1/5 */

static const int DX[NW] = {1, -1, 0, 0};
static const int DY[NW] = {0, 0, 1, -1};
static const int MX[NMOVE] = {0, 1, -1, 0, 0};
static const int MY[NMOVE] = {0, 0, 0, 1, -1};

#define MAXCK 64
#define MAXPAT 64

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
    int c[NWALK], m = 0;
    for (int r = 0; r < k; r++) {
        if (roots[r] == cls[0]) continue;
        int n = 0;
        for (int i = 1; i < NWALK; i++) n += cls[i] == roots[r];
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
    if (argc < 4) { fprintf(stderr, "usage: cwalk5 trials seed T1 [T2 ...]\n"); return 2; }
    long long trials = atoll(argv[1]);
    uint64_t seed = strtoull(argv[2], NULL, 10);
    int nck = argc - 3;
    if (nck > MAXCK) { fprintf(stderr, "too many checkpoints\n"); return 2; }
    long long ck[MAXCK];
    for (int i = 0; i < nck; i++) ck[i] = atoll(argv[3 + i]);
    Hist *H = calloc(nck, sizeof(Hist));
    fprintf(stderr, "cwalk5: %d lazy walkers, stay 1/%d, %lld trials\n", NWALK, NMOVE, trials);

#pragma omp parallel
    {
        Hist *h = calloc(nck, sizeof(Hist));
#pragma omp for schedule(dynamic, 1024)
        for (long long tr = 0; tr < trials; tr++) {
            uint64_t s = seed * 0xd1342543de82ef95ULL + (uint64_t)tr;
            splitmix(&s);
            int px[NWALK], py[NWALK], cls[NWALK], roots[NWALK], k = NWALK;
            px[0] = py[0] = 0;
            for (int i = 0; i < NW; i++) { px[i + 1] = DX[i]; py[i + 1] = DY[i]; }
            for (int i = 0; i < NWALK; i++) cls[i] = roots[i] = i;   /* roots index walkers */
            int ci = 0;
            uint64_t bits = 0; int nb = 0;
            for (long long t = 0; ci < nck; t++) {
                if (t > 0) {
                    for (int r = 0; r < k; r++) {
                        /* uniform move in {0..4} from a 32-bit word: bias 5/2^32 */
                        if (nb < 32) { bits = splitmix(&s); nb = 64; }
                        int d = (int)(((bits & 0xffffffffULL) * NMOVE) >> 32);
                        bits >>= 32; nb -= 32;
                        px[roots[r]] += MX[d]; py[roots[r]] += MY[d];
                    }
                    for (int a = 0; a < k; a++)
                        for (int b = a + 1; b < k; b++) {
                            int ra = roots[a], rb = roots[b];
                            if (px[ra] == px[rb] && py[ra] == py[rb]) {
                                for (int i = 0; i < NWALK; i++) if (cls[i] == rb) cls[i] = ra;
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
