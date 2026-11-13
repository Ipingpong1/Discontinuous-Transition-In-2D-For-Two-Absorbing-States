"""Front v in physical-time windows (per-run linear fits of W/2 against t r), seed mean ± sem.

    python3 vwin.py RUNS_DIR/E11_lazy_front_p99 [RUNS_DIR/E10_lazy_front ...]
"""
import glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fronts
WINS = [(5e2, 1e5), (2.5e4, 1e5), (1e5, 3e5)]
for d in sys.argv[1:]:
    by = {}
    for txt in sorted(glob.glob(os.path.join(d, "out", "*.txt"))):
        n = os.path.basename(txt)[:-4]
        P = fronts.params(n)
        r = float(P.get("r", 1))
        t, W, cen, Ly, rho_b = fronts.track(txt, txt[:-4] + ".prof")
        tp = t * r
        row = []
        for a, b in WINS:
            sel = (tp >= a) & (tp <= b) & (W > 0.1 * Ly) & (W < 0.9 * Ly) & (cen > 0.7)
            row.append(np.polyfit(tp[sel], W[sel], 1)[0] / 2 if sel.sum() >= 4 else np.nan)
        by.setdefault((P["p"], r, float(P["q"])), []).append(row)
    print(f"== {d}: physical v (W/2 slope), seed mean ± sem, windows " + ", ".join(f"[{a:.0e},{b:.0e}]" for a, b in WINS))
    for (p, r, q), rows in sorted(by.items()):
        A = np.array(rows)
        cells = []
        for j in range(A.shape[1]):
            c = A[:, j][~np.isnan(A[:, j])]
            cells.append(f"{c.mean():+.5f} ± {c.std(ddof=1) / np.sqrt(len(c)) if len(c) > 1 else 0:.5f} ({len(c)})" if len(c) else "        -          ")
        print(f"p {p} r {r:g} q {q:<9g} | " + " | ".join(cells))
