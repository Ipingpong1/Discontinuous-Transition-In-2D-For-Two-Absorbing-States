#!/usr/bin/env python3
"""Reduce profile-based runs to small CSVs for the paper figures (runs on Torch, in the container).

    python3 extract_prof.py RUNS_DIR OUT_DIR fronts E1_front E4_voter ...
    python3 extract_prof.py RUNS_DIR OUT_DIR walls E8b_wall E16_wall ...

fronts: one row per strip run, from campaign/fronts.py front_v (v per step; fit window and
        interior checks as documented there) -> OUT_DIR/<EXP>_fronts.csv
walls:  per run, the wall.py series t, I, sig, rows, rin, alive, W, rc -> OUT_DIR/<EXP>/<run>.csv,
        plus the time-averaged interface profile across the wall (rows within Ly/4 of each wall,
        both walls folded onto one signed distance) at the last output time and averaged over
        t >= T/3 -> OUT_DIR/<EXP>/<run>_prof.csv
"""
import csv
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser("~/dn/campaign"))
sys.path.insert(0, os.path.expanduser("~/dn"))
import fronts  # noqa: E402
import wall  # noqa: E402


def do_fronts(runs, out, exp):
    rows = []
    for p in sorted(glob.glob(f"{runs}/{exp}/out/*.txt")):
        if not os.path.exists(p[:-4] + ".prof"):
            continue
        name = os.path.basename(p)[:-4]
        try:
            v, n, tmax, t_order, wend, rho_b = fronts.front_v(p)
        except Exception as e:  # unfinished or empty runs
            print("skip", name, e)
            continue
        rows.append([name, v, n, tmax, t_order, wend, rho_b])
    with open(os.path.join(out, f"{exp}_fronts.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run", "v", "npts", "t_fitmax", "t_order", "W_end_over_Ly", "rho_b"])
        w.writerows(rows)
    print(exp, len(rows), "front runs")


def do_walls(runs, out, exp):
    od = os.path.join(out, exp)
    os.makedirs(od, exist_ok=True)
    n = 0
    for p in sorted(glob.glob(f"{runs}/{exp}/out/*.prof")):
        prefix = p[:-5]
        name = os.path.basename(prefix)
        try:
            t, I, sig, rows, rin, al, W, rc = wall.run_series(prefix)
            Lx, Ly, tt, alive, iface = wall.read(prefix)
        except Exception as e:
            print("skip", name, e)
            continue
        with open(os.path.join(od, name + ".csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["t", "I", "sig", "rows", "rin", "alive", "W", "rc"])
            for k in range(len(t)):
                w.writerow([int(t[k]), f"{I[k]:.4f}", f"{sig[k]:.3f}", f"{rows[k]:.1f}", f"{rin[k]:.4f}",
                            f"{al[k]:.5f}", f"{W[k]:.1f}", f"{rc[k]:.4f}"])
        d = iface / Lx
        prof = []
        for c in (Ly // 2, 0):
            y = (np.arange(Ly) - c + Ly // 2) % Ly - Ly // 2
            sel = np.flatnonzero(np.abs(y) < Ly // 4)
            sel = sel[np.argsort(y[sel])]
            prof.append((y[sel], d[:, sel]))
        ys = prof[0][0]
        late = tt >= tt[-1] / 3
        last = 0.5 * (prof[0][1][-1] + prof[1][1][-1])
        mean_late = 0.5 * (prof[0][1][late].mean(0) + prof[1][1][late].mean(0))
        # centred: shift every late record so that its interface-weighted mean row sits at y = 0
        acc, nrec = np.zeros(len(ys)), 0
        for yy, dd in prof:
            for row in dd[late]:
                if row.sum() <= 0:
                    continue
                c = int(round((row * yy).sum() / row.sum()))
                acc += np.roll(row, -c)
                nrec += 1
        centred = acc / max(nrec, 1)
        with open(os.path.join(od, name + "_prof.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["y", "iface_last", "iface_late_mean", "iface_late_centred"])
            for k in range(len(ys)):
                w.writerow([int(ys[k]), f"{last[k]:.5f}", f"{mean_late[k]:.5f}", f"{centred[k]:.5f}"])
        n += 1
    print(exp, n, "wall runs")


def main():
    runs, out, mode = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    for exp in sys.argv[4:]:
        (do_fronts if mode == "fronts" else do_walls)(runs, out, exp)


if __name__ == "__main__":
    main()
