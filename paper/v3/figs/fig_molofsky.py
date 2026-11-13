"""Five-cell closure and one-sided expansion estimates from v3 raw-profile reductions.

Run campaign/v3/analysis/molofsky.py before this figure script.
"""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import ROOT, PAGE, label, save

PROJECT = Path(ROOT)
RESULTS = PROJECT / "campaign/v3/results/molofsky"


def read_table(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key, value in row.items():
            if value in ("True", "False"):
                row[key] = value == "True"
            elif value == "":
                row[key] = np.nan
            else:
                try:
                    row[key] = float(value)
                except ValueError:
                    pass
    return rows


def main():
    paths = [PROJECT/"campaign/theory/molofsky_1e7.csv",
             RESULTS/"groups.csv", RESULTS/"extrapolations.csv"]
    theory = np.genfromtxt(paths[0], delimiter=",", names=True)
    theory = theory[theory["T"] >= 1]
    groups, fits = read_table(paths[1]), read_table(paths[2])
    fig, axes = plt.subplots(2, 2, figsize=(PAGE, 5.6))
    a, b, c, d = axes.ravel()
    for key, color, text in (("th_a", "C0", r"$f_T'(0)=0$"),
                             ("th_f", "C3", r"$\int_0^{1/2}f_T\,du=0$"),
                             ("th_h", "C2", r"$f_T'(1/2)=0$")):
        a.plot(theory["T"], theory[key], "o-", color=color, ms=2.5, label=text)
    angles = [row["angle"] for row in fits]
    a.axhspan(min(angles), max(angles), color="0.5", alpha=.3,
              label="onset extrapolations")
    a.axhspan(104.5, 105.0, color="C0", alpha=.17,
              label="initial seed-growth crossover")
    a.set(xscale="log", xlabel=r"Dual-walk cutoff $T$", ylabel=r"Closure angle ($^\circ$)", ylim=(99, 121))
    a.legend(fontsize=6.8, loc="upper right")
    label(a, "(a)")
    b.errorbar(theory["T"], theory["width"], yerr=theory["sd_w"], fmt="o-", color="k", ms=2.5)
    b.set(xscale="log", yscale="log", xlabel=r"Dual-walk cutoff $T$", ylabel=r"Bistable sector width ($^\circ$)")
    label(b, "(b)", x=.90)
    for ax, p2, marker in ((c, .45, "o"), (d, .50, "s")):
        cut_fits = [row for row in fits if row["p2"] == p2]
        fit = next(row for row in cut_fits if row["estimator"] == "early_adaptive")
        zero = fit["extrapolated_zero"]
        zeros = [row["extrapolated_zero"] for row in cut_fits]
        ax.axvspan((min(zeros)-zero)*1e4, (max(zeros)-zero)*1e4, color="0.6", alpha=.25)
        ax.axhline(0, color="0.5", linewidth=.6)
        for row in [row for row in groups if row["p2"] == p2]:
            wide = row["Lx"] == 4096
            ax.errorbar((row["p1"]-zero)*1e4, row["early_adaptive"]*1e3,
                        yerr=row["early_adaptive_sem"]*1e3,
                        fmt="^" if wide else marker,
                        markerfacecolor="k" if row["core_pass"] else "white",
                        color="0.35" if wide else "k", ms=4, capsize=2,
                        elinewidth=.7, markeredgewidth=.8)
        solid = np.array([fit["p1_min"], fit["p1_max"]])
        for interval, style in ((solid, "-"), (np.array([zero, solid[0]]), "--")):
            ax.plot((interval-zero)*1e4, (fit["gradient"]*interval+fit["intercept"])*1e3,
                    style, color="C3", linewidth=1)
        ax.plot([], [], marker, color="k", label=r"core decline $\leq0.03$")
        ax.plot([], [], marker, color="k", markerfacecolor="white", label="larger core decline")
        ax.plot([], [], "^", color="0.35", markerfacecolor="white", label=r"$L_x=4096$ check")
        ax.set(xlim=(-3.5, 2.6), ylim=(-2.7, 3.8),
               xlabel=r"$(p_1-p_1^{\rm ext})\times10^4$",
               ylabel=r"Adaptive half-width slope $\times10^3$")
        ax.text(.05, .83, f"$p_2={p2:.2f}$\n$p_1^{{\\rm ext}}={zero:.6f}$",
                ha="left", va="top", transform=ax.transAxes, fontsize=7)
        ax.legend(loc="lower right", fontsize=6.5)
        label(ax, "(c)" if p2 == .45 else "(d)")
    fig.tight_layout(h_pad=1.8, w_pad=1.3)
    save(fig, "fig_molofsky")
    plt.close(fig)
    metadata = {"inputs": {str(path.relative_to(PROJECT)): hashlib.sha256(path.read_bytes()).hexdigest()
                            for path in paths},
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "interpolation": "only positive adaptive branches; dashed one-sided extrapolation",
                "band": "range of fixed/adaptive/integrated zeros across two windows; sensitivity range"}
    (RESULTS/"figure_provenance.json").write_text(json.dumps(metadata, indent=2)+"\n")


if __name__ == "__main__":
    main()
