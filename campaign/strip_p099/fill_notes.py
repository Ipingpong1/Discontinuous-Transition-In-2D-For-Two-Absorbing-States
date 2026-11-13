#!/usr/bin/env python3
"""Fill the table placeholders of NOTES.md ([RUNTABLE], [SEEDTABLE], [CROSSTABLE], [TALLTABLE],
[BULKTABLE]) from analyze.py; the verdict placeholders are edited by hand afterwards.

    python3 fill_notes.py        # rewrites NOTES.md in place (only placeholders still present)
"""
import io
import os
import sys

import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, D)
import analyze  # noqa: E402
import sec_tables  # noqa: E402

QN = {0.995: 0.00064, 0.99: 0.001282}


def fv(v):
    return "-" if v is None or np.isnan(v) else f"{v:+.5f}"


def run_table(rows, geoms):
    out = io.StringIO()
    print("| p | q | x_new | Lx x Ly | s | t_end | v_early c | v_late o | v_late c | v_late An | n | interior 1e5 / t_x / end | status |", file=out)
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|", file=out)
    for r in rows:
        if (r["Lx"], r["Ly"]) not in geoms:
            continue
        c1 = r["cen"][1]
        print(f"| {r['p']} | {r['q']:.6g} | {r['q'] / QN[r['p']]:.3f} | {r['Lx']}x{r['Ly']} | {r['s']} | {r['T']:.2g} "
              f"| {fv(r['v'][('c', 'E')])} | {fv(r['v'][('o', 'L')])} | {fv(r['v'][('c', 'L')])} | {fv(r['v'][('An', 'L')])} "
              f"| {r['n'][('o', 'L')]} | {'-' if np.isnan(c1) else f'{c1:.2f}'} / {r['cen_x']:.2f} / {r['cen_q']:.2f} | {r['status']} |", file=out)
    return out.getvalue()


def seed_table(rows, geom):
    out = io.StringIO()
    print("| p | q | x_new | v_E | v_M | v_L | v_L (An) | v_L (c) | status |", file=out)
    print("|---|---|---|---|---|---|---|---|---|", file=out)
    for p in (0.995, 0.99):
        rs = [r for r in rows if r["p"] == p and (r["Lx"], r["Ly"]) == geom]
        for q in sorted(set(r["q"] for r in rs)):
            rq = [r for r in rs if r["q"] == q]
            cells = []
            for key in (("o", "E"), ("o", "M"), ("o", "L"), ("An", "L"), ("c", "L")):
                m, e, n = analyze.seed_mean([r["v"][key] for r in rq])
                cells.append("-" if n == 0 else (f"{m:+.5f}({e * 1e5:.0f}) [{n}]" if n > 1 else f"{m:+.5f} [1]"))
            st = "/".join(r["status"] for r in sorted(rq, key=lambda r: r["s"]))
            print(f"| {p} | {q:.6g} | {q / QN[p]:.3f} | " + " | ".join(cells) + f" | {st} |", file=out)
    return out.getvalue()


def cross_table(rows, geom):
    res = analyze.cross([r for r in rows if (r["Lx"], r["Ly"]) == geom], out=open(os.devnull, "w"),
                        usable=("stationary", "filled", "decaying"))
    out = io.StringIO()
    for p in (0.995, 0.99):
        for key, name in ((("o", "L"), "outer, late"), (("An", "L"), "integr./interior, late"), (("A", "L"), "integrated, late"),
                          (("c", "L"), "Table I estimator, late"), (("c", "E"), "Table I estimator, early [2e3, 1e5]")):
            e = res.get((p, geom[0], geom[1], key), {})
            if "qstar" in e:
                print(f"| {p} | {name} | {', '.join(f'{q:.6g}' for q in e['pts'])} | {e['qstar']:.6f} | {e['sd']:.6f} | {e['slope']:.1f} "
                      f"| {e['qstar'] * e['slope']:.3f} | {e['qstar'] / (0.5 * (1 - p)):.4f} |", file=out)
            else:
                tp = e.get("twopt")
                print(f"| {p} | {name} | two-point only | {'-' if tp is None else f'{tp:.6f}'} | - | - | - | "
                      f"{'-' if tp is None else f'{tp / (0.5 * (1 - p)):.4f}'} |", file=out)
    return out.getvalue()


def main():
    notes = open(f"{D}/NOTES.md").read()
    rows = analyze.strips([f"{D}/strips", f"{D}/wide", f"{D}/tall"], out=open(os.devnull, "w"))
    rep = {
        "[RUNTABLE]": run_table(rows, {(2048, 4096), (4096, 4096)}),
        "[SEEDTABLE]": seed_table(rows, (2048, 4096)),
        "[CROSSTABLE]": cross_table(rows, (2048, 4096)).rstrip("\n"),
        "[TALLTABLE]": run_table(rows, {(1024, 8192)}) + "\n" + seed_table(rows, (1024, 8192)),
    }
    out = io.StringIO()
    analyze.bulk([f"{D}/bulk"] + analyze.E6[0.995] + analyze.E6[0.99], out=out)
    rep["[BULKTABLE]"] = "\n".join(l for l in out.getvalue().splitlines() if l.startswith("|"))
    for k, v in rep.items():
        if k in notes:
            notes = notes.replace(k, v)
    open(f"{D}/NOTES.md", "w").write(notes)
    print("filled:", [k for k in rep if k not in notes], "left:", [k for k in rep if k in notes])


if __name__ == "__main__":
    main()
