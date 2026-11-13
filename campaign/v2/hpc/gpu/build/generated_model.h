/* Generated from dnsim.c; SHA256 dd56be11b1e7865b52ec4607432eb5733f2306b4a714312690df801793f22784. Do not edit. */
#ifndef DNSIM_GENERATED_MODEL_H
#define DNSIM_GENERATED_MODEL_H
#include <stdint.h>
typedef uint64_t u64;
typedef uint32_t u32;
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
#endif
