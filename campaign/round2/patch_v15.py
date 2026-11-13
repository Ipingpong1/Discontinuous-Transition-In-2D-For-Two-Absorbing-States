"""Apply the v1.5 changes to a copy of the v1.4 source (exact-match replacements, each unique)."""
import sys
path = sys.argv[1]
src = open(path).read()
edits = []

def rep(old, new):
    edits.append((old, new))

rep(""" * Biol. 55, 270 (1999) on it.  The Moore path is unchanged.
 */""", """ * Biol. 55, 270 (1999) on it.  The Moore path is unchanged.
 *
 * v1.5 (10-31, campaign/round2): --pdir a,b bins the --profile along the normal of a tilted
 * band, s = (a y - b x) mod L (Lx = Ly = L; (1,0) is the row profile), for direct front tracking
 * of the orientation test; --snap file --snap-every n [--snap-from t] appends full-lattice
 * records (u64 t, then the W x Ly words) for roughness measurements; --init hole:p0:R is a
 * disk of the background inside a random:p0 lattice (an absorbing droplet in the active phase).
 * The dynamics and every existing output are unchanged.
 */""")

# init: hole mode
rep("""            else if (!strcmp(mode, "disk")) {
                in = 0;
                double cy = g->Ly / 2.0, cx = g->Lx / 2.0;
                for (int b = 0; b < 64; b++) {
                    double dx = 64.0 * x + b + 0.5 - cx, dy = y + 0.5 - cy;
                    if (dx * dx + dy * dy <= par * par) in |= 1ULL << b;
                }
            } else {""", """            else if (!strcmp(mode, "disk") || !strcmp(mode, "hole")) {
                /* disk: random inside, bg outside; hole: the complement (bg inside the disk) */
                in = 0;
                double cy = g->Ly / 2.0, cx = g->Lx / 2.0;
                for (int b = 0; b < 64; b++) {
                    double dx = 64.0 * x + b + 0.5 - cx, dy = y + 0.5 - cy;
                    if (dx * dx + dy * dy <= par * par) in |= 1ULL << b;
                }
                if (mode[0] == 'h') in = ~in;
            } else {""")

# tilted profile + snapshot helpers, before truncate_profile
rep("""/* Same for the binary profile: keep whole records with t <= t0. */
static void truncate_profile(const char *path, int Ly, u64 t0)
{
    FILE *f = fopen(path, "rb");
    if (!f) return;
    long rec = 8 + 8L * Ly, keep = 0;""", """/* Profile along the normal of a tilted band (--pdir a,b): bin s = (a y - b x) mod L of each
 * cell, Lx = Ly = L, so that (a, b) = (1, 0) gives the row profile.  The lines a y - b x = s are
 * 1/sqrt(a^2 + b^2) apart in normal distance and carry L cells each when gcd(a, b) = 1. */
static void measure_dir(const Grid *g, int a, int b, u32 *bin_alive, u32 *bin_iface)
{
    const int W = g->W, Ly = g->Ly, L = g->Lx;
    const u64 *A = g->a;
    memset(bin_alive, 0, sizeof(u32) * (size_t)L);
    memset(bin_iface, 0, sizeof(u32) * (size_t)L);
    for (int y = 0; y < Ly; y++) {
        const u64 *U = A + (size_t)(y ? y - 1 : Ly - 1) * W;
        const u64 *C = A + (size_t)y * W;
        const u64 *D = A + (size_t)(y + 1 < Ly ? y + 1 : 0) * W;
        long s0 = ((long)a * y) % L;
        for (int x = 0; x < W; x++) {
            int xl = x ? x - 1 : W - 1, xr = x + 1 < W ? x + 1 : 0;
            u64 c = C[x], m = iface_word(U, C, D, x, xl, xr);
            for (int bit = 0; bit < 64; bit++) {
                long s = (s0 - (long)b * (64L * x + bit)) % L;
                if (s < 0) s += L;
                bin_alive[s] += (u32)((c >> bit) & 1);
                bin_iface[s] += (u32)((m >> bit) & 1);
            }
        }
    }
}

/* Snapshot record: u64 t, then the W x Ly words of the lattice (same layout as a checkpoint). */
static void write_snap(FILE *f, const Grid *g, long long t)
{
    fwrite(&t, 8, 1, f);
    fwrite(g->a, 8, (size_t)g->W * g->Ly, f);
    fflush(f);
}

/* Binary files of fixed-size records that start with the time: keep whole records with t <= t0. */
static void truncate_records(const char *path, long rec, u64 t0)
{
    FILE *f = fopen(path, "rb");
    if (!f) return;
    long keep = 0;""")
rep("""        keep = i + 1;
    }
    fclose(f);
    if (truncate(path, keep * rec)) { perror(path); exit(1); }
}
""", """        keep = i + 1;
    }
    fclose(f);
    if (truncate(path, keep * rec)) { perror(path); exit(1); }
}

/* Same for the binary profile: keep whole records with t <= t0. */
static void truncate_profile(const char *path, int Ly, u64 t0)
{
    truncate_records(path, 8 + 8L * Ly, t0);
}
""")

# usage
rep("""        "  --init spec             random:p0 | strip:p0:frac | disk:p0:R, default random:0.5\\n"
""", """        "  --init spec             random:p0 | strip:p0:frac | disk:p0:R | hole:p0:R (random\\n"
        "                          outside a disk of the background), default random:0.5\\n"
""")
rep("""        "  --profile file          binary per-row alive/iface counts at output times\\n"
""", """        "  --profile file          binary per-row alive/iface counts at output times\\n"
        "  --pdir a,b              ... binned by s = (a y - b x) mod L instead of by row\\n"
        "                          (tilted band of tilt.py; needs Lx = Ly)\\n"
        "  --snap file             append the whole lattice (u64 t, then the words) ...\\n"
        "  --snap-every n          ... every n steps (and at T) ...\\n"
        "  --snap-from t           ... from time t on, default 0\\n"
""")

# options
rep("""    const char *load = NULL, *save = NULL, *out = NULL, *prof = NULL, *ckpt = NULL;
""", """    const char *load = NULL, *save = NULL, *out = NULL, *prof = NULL, *ckpt = NULL;
    const char *pdir = NULL, *snap = NULL;
    long long snap_every = 0, snap_from = 0;
    int pa = 1, pb = 0;
""")
rep("""        else if (!strcmp(a, "--profile")) prof = ARG();
""", """        else if (!strcmp(a, "--profile")) prof = ARG();
        else if (!strcmp(a, "--pdir")) pdir = ARG();
        else if (!strcmp(a, "--snap")) snap = ARG();
        else if (!strcmp(a, "--snap-every")) snap_every = atoll(ARG());
        else if (!strcmp(a, "--snap-from")) snap_from = atoll(ARG());
""")
rep("""    if ((resume || ckpt_every) && !ckpt) { fprintf(stderr, "--resume/--ckpt-every need --ckpt\\n"); return 2; }
""", """    if ((resume || ckpt_every) && !ckpt) { fprintf(stderr, "--resume/--ckpt-every need --ckpt\\n"); return 2; }
    if (pdir) {
        if (sscanf(pdir, "%d,%d", &pa, &pb) != 2 || pa < 0 || pb < 0 || (pa == 0 && pb == 0) || !prof || Lx != Ly) {
            fprintf(stderr, "--pdir a,b needs a, b >= 0 (not both 0), --profile, and Lx = Ly\\n"); return 2;
        }
    }
    if ((snap && !snap_every) || (!snap && (snap_every || snap_from))) {
        fprintf(stderr, "--snap needs --snap-every n (and vice versa)\\n"); return 2;
    }
""")
rep("""        if (prof) truncate_profile(prof, Ly, t0);
""", """        if (prof) truncate_profile(prof, Ly, t0);
        if (snap) truncate_records(snap, 8 + 8L * (Lx / 64) * Ly, t0);
""")
rep("""    FILE *fp = prof ? fopen(prof, resumed ? "ab" : "wb") : NULL;
    if (prof && !fp) { perror(prof); return 1; }
""", """    FILE *fp = prof ? fopen(prof, resumed ? "ab" : "wb") : NULL;
    if (prof && !fp) { perror(prof); return 1; }
    FILE *fs = snap ? fopen(snap, resumed ? "ab" : "wb") : NULL;
    if (snap && !fs) { perror(snap); return 1; }
""")
rep("""    fprintf(fo, "# dnsim v1.4\\n# Lx %d Ly %d T %lld seed %llu threads %d\\n", Lx, Ly, T, seed, nthreads);
    hdr_model(fo, molofsky, ftab, rule, alpha, voter, jflip, &M, load ? "-" : init, bg, load ? load : "-", t0);
""", """    fprintf(fo, "# dnsim v1.5\\n# Lx %d Ly %d T %lld seed %llu threads %d\\n", Lx, Ly, T, seed, nthreads);
    hdr_model(fo, molofsky, ftab, rule, alpha, voter, jflip, &M, load ? "-" : init, bg, load ? load : "-", t0);
    if (pdir) fprintf(fo, "# profile_dir %d,%d (bins s = (a y - b x) mod L)\\n", pa, pb);
    if (snap) fprintf(fo, "# snap every %lld from %lld\\n", snap_every, snap_from);
""")
rep("""        fprintf(hs, "# dnsim v1.4 reps\\n""", """        fprintf(hs, "# dnsim v1.5 reps\\n""")
# initial output
rep("""        if (fp) { long long tt = (long long)t0; fwrite(&tt, 8, 1, fp); fwrite(ra, 4, Ly, fp); fwrite(ri, 4, Ly, fp); }
        fflush(fo);
    }
""", """        if (fp) {
            long long tt = (long long)t0;
            if (pdir) measure_dir(&g, pa, pb, ra, ri);
            fwrite(&tt, 8, 1, fp); fwrite(ra, 4, Ly, fp); fwrite(ri, 4, Ly, fp);
        }
        if (fs && (long long)t0 >= snap_from) write_snap(fs, &g, (long long)t0);
        fflush(fo);
    }
""")
rep("""                if (fp) { fwrite(&t, 8, 1, fp); fwrite(ra, 4, Ly, fp); fwrite(ri, 4, Ly, fp); }
                fflush(fo);
            }
            if (stop_abs && absorbing && o.iface == 0) { absorbed = 1; break; }
        }
""", """                if (fp) {
                    if (pdir) measure_dir(&g, pa, pb, ra, ri);
                    fwrite(&t, 8, 1, fp); fwrite(ra, 4, Ly, fp); fwrite(ri, 4, Ly, fp);
                }
                fflush(fo);
            }
            if (stop_abs && absorbing && o.iface == 0) { absorbed = 1; break; }
        }
        if (fs && t >= snap_from && (t % snap_every == 0 || t == T)) write_snap(fs, &g, t);
""")
rep("""            if (fp) fflush(fp);
            save_atomic(&g, (u64)t, ckpt);
""", """            if (fp) fflush(fp);
            if (fs) fflush(fs);
            save_atomic(&g, (u64)t, ckpt);
""")
rep("""    if (fp) fclose(fp);
    free(g.a); free(nx.a); free(tl); free(ra); free(ri);
""", """    if (fp) fclose(fp);
    if (fs) fclose(fs);
    free(g.a); free(nx.a); free(tl); free(ra); free(ri);
""")

for old, new in edits:
    n = src.count(old)
    if n != 1:
        sys.exit(f"edit matched {n} times:\n{old[:120]}")
    src = src.replace(old, new)
open(path, "w").write(src)
print(f"applied {len(edits)} edits")
