#!/usr/bin/env python3
"""Generate CUDA update helpers and a C host from the unchanged CPU simulator."""
from pathlib import Path
import argparse
import hashlib
import re

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--source", type=Path, default=HERE / "../../../round2/dnsim.c")
parser.add_argument("--out", type=Path, default=HERE / "build")
args = parser.parse_args()
source = args.source.resolve()
text = source.read_text()
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)
digest = hashlib.sha256(source.read_bytes()).hexdigest()
banner = f"/* Generated from {source.name}; SHA256 {digest}. Do not edit. */\n"

def replace_once(s, old, new):
    if s.count(old) != 1:
        raise ValueError(f"Expected one source anchor, found {s.count(old)}: {old[:80]!r}")
    return s.replace(old, new, 1)

def function(name):
    match = re.search(r"static inline(?: __attribute__\(\(always_inline\)\))?\s+\w+\s+" +
                      re.escape(name) + r"\s*\(", text)
    if match is None:
        raise ValueError(f"Missing function {name}")
    begin = text.index("{", match.start())
    depth = 1
    end = begin + 1
    while depth:
        depth += (text[end] == "{") - (text[end] == "}")
        end += 1
    result = text[match.start():end]
    result = result.replace("static inline __attribute__((always_inline))",
                            "static __device__ __forceinline__")
    result = result.replace("static inline", "static __device__ __forceinline__")
    return result.replace("__builtin_ctz(", "gpu_ctz(")

model_begin = text.index("#define MAXV 18")
model_end = text.index("} Model;", model_begin) + len("} Model;")
model = text[model_begin:model_end]
types = "typedef uint64_t u64;\ntypedef uint32_t u32;"
(out / "generated_model.h").write_text(
    banner + "#ifndef DNSIM_GENERATED_MODEL_H\n#define DNSIM_GENERATED_MODEL_H\n"
    "#include <stdint.h>\n" + types + "\n" + model + "\n#endif\n")
rstream = re.search(r"typedef struct \{[^\n]+\} rstream;", text).group()
functions = ["philox4x32_10", "rs_init", "rs_next", "bern_table",
             "count_eq", "count_eq_vn", "step_word"]
device = banner + '#include "generated_model.h"\n' + rstream + "\n"
device += "static __device__ __forceinline__ int gpu_ctz(u32 v) { return __ffs((int)v) - 1; }\n"
device += "\n\n".join(function(name) for name in functions) + "\n"
(out / "generated_device.cuh").write_text(device)

host = replace_once(text, types, '#include "gpu_api.h"')
host = replace_once(host, model, "/* Model is shared verbatim through gpu_api.h. */")
host = replace_once(host, "int main(int argc, char **argv)", r"""
/* Optional monitoring sidecar; no lattice or RNG state is affected. */
static void gpu_phase(const char *phase)
{
    const char *path = getenv("DNSIM_PHASE_FILE");
    if (!path || !*path) return;
    char temporary[4096];
    if (snprintf(temporary, sizeof temporary, "%s.tmp", path) >= (int)sizeof temporary) {
        fprintf(stderr, "phase path too long\n"); exit(2);
    }
    FILE *file = fopen(temporary, "w");
    if (!file) { perror(temporary); exit(1); }
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    fprintf(file, "%s %.9f\n", phase, now.tv_sec + 1e-9 * now.tv_nsec);
    if (fclose(file) || rename(temporary, path)) { perror(path); exit(1); }
}

/* Return the next scheduled multiple, capped at the final time. */
static long long gpu_next_periodic(long long now, long long period,
                                   long long earliest, long long limit)
{
    if (period <= 0 || period > limit || earliest > limit) return limit;
    long long first = now + 1;
    if (earliest > first) first = earliest;
    long long q = first / period + (first % period != 0);
    long long result = q * period; /* bounded by limit + period < 2^33 */
    return result < limit ? result : limit;
}

int main(int argc, char **argv)""")
host = replace_once(host, 'if (!strcmp(a, "--selftest")) return selftest();', r"""if (!strcmp(a, "--selftest")) {
            fprintf(stderr, "Use the CPU binary for --selftest and validate.py for CPU/CUDA equivalence.\n");
            return 2;
        }""")
host = replace_once(host, "    Model M = {0};\n    M.key[0] = (u32)seed;", r"""    if (reps) {
        fprintf(stderr, "CUDA dense backend does not support --reps; use the CPU binary.\n");
        return 2;
    }
    if (every < 0 || ckpt_every < 0 || snap_every < 0 || ppd < 0) {
        fprintf(stderr, "output intervals and --ppd must be nonnegative\n");
        return 2;
    }
    Model M = {0};
    M.key[0] = (u32)seed;""")
host = replace_once(host, "    Grid g = grid_new(Lx, Ly), nx = grid_new(Lx, Ly);",
                     "    Grid g = grid_new(Lx, Ly);")
host = replace_once(host, "    /* output schedule */", r"""    fprintf(fo, "# backend CUDA dense; source SHA256 """ + digest + r"""\n");
    fprintf(fo, "# absorption checked only at scheduled time-series outputs; use CPU for exact lifetimes\n");
    /* output schedule */""")
start = host.index("    struct timespec ts0, ts1;", host.index("int main("))
stop_anchor = "    clock_gettime(CLOCK_MONOTONIC, &ts1);"
end = host.index(stop_anchor, start) + len(stop_anchor)
host = host[:start] + (HERE / "host_loop.inc").read_text().rstrip() + host[end:]
host = replace_once(host, "    free(g.a); free(nx.a); free(tl); free(ra); free(ri);",
                     "    free(g.a); free(tl); free(ra); free(ri);")
(out / "generated_host.c").write_text(banner + host)
(out / "source.sha256").write_text(digest + "\n")
print(f"Generated {out} from {source} (SHA256 {digest})")
