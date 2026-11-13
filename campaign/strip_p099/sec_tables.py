#!/usr/bin/env python3
"""Compact tables for NOTES.md (run by run, seed means, crossings, bulk), from analyze.py.

    python3 sec_tables.py > tables.md
"""
import os
import sys

import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, D)
import analyze  # noqa: E402


def fv(v):
    return "-" if v is None or np.isnan(v) else f"{v:+.5f}"


def main():
    rows = analyze.strips([f"{D}/strips", f"{D}/wide", f"{D}/tall"], out=open(os.devnull, "w"))
    print("### Run by run (strips)\n")
    print("v in rows/step per front; o = outer envelope (0.3 rho_b), c = Table I estimator (centre walk, rho_b/2), "
          "An = integrated width / current interior density; late = [t_x, T]; early = [2e3, 1e5] with the c estimator "
          "and its interior condition (the Table I protocol); n = records in the late window (outer); interior = central "
          "interface density relative to t ~ 1e3 at t = 1e5, t_x, T (or the last record); status as in analyze.py\n")
    print("| p | q | x_old | Lx x Ly | strip | s | T | t_end | v_early c | v_late o | v_late c | v_late An | n | interior 1e5 / t_x / end | status |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        tx = analyze.TX[r["p"]]
        T = r["T"]
        cen_T = r["cen"][-1] if not np.isnan(r["cen"][-1]) else r["cen_q"]
        c1 = r["cen"][1]
        cx = r["cen_x"]
        print(f"| {r['p']} | {r['q']:.6g} | {r['q'] / analyze.QS_OLD[r['p']]:.3f} | {r['Lx']}x{r['Ly']} | Ly/{round(1 / r['frac'])} "
              f"| {r['s']} | {1.6e6 if r['p'] == 0.995 else 1.2e6:.1e} | {T:.2g} | {fv(r['v'][('c', 'E')])} | {fv(r['v'][('o', 'L')])} "
              f"| {fv(r['v'][('c', 'L')])} | {fv(r['v'][('An', 'L')])} | {r['n'][('o', 'L')]} | "
              f"{'-' if np.isnan(c1) else f'{c1:.2f}'} / {cx:.2f} / {r['cen_q']:.2f} | {r['status']} |")
    print("\n### Seed means and crossings\n")
    analyze.cross(rows, out=sys.stdout, usable=("stationary", "filled", "decaying"))
    print("\n### Bulk\n")
    analyze.bulk([f"{D}/bulk"] + analyze.E6[0.995] + analyze.E6[0.99], out=sys.stdout)


if __name__ == "__main__":
    main()
