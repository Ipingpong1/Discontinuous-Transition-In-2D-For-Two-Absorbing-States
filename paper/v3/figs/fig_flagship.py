"""Compact view of bulk activity, front reversal, and complementary droplets."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import numpy as np

from common import ROOT, PAGE, QSTAR05, QSTAR05_HALF_RANGE, label, save
from campaign.v2.analysis.io import provenance, write_json

RESULTS = Path(ROOT) / "campaign/v3/results/fronts_review"


def read_table(name):
    with (RESULTS / (name + ".csv")).open() as stream:
        return [{key: value if key in ("run", "normal", "observable", "kind", "estimator", "eligible")
                 else float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def shade_reference(ax, horizontal=False):
    span = ax.axhspan if horizontal else ax.axvspan
    span(QSTAR05 - QSTAR05_HALF_RANGE, QSTAR05 + QSTAR05_HALF_RANGE,
         color=".8", alpha=.4, lw=0, zorder=-5)


def main():
    bulk = read_table("bulk_summary")
    bulk_runs = read_table("bulk_runs")
    points = read_table("front_points")
    runs = read_table("front_runs")
    crossings = read_table("front_crossings")
    droplets = read_table("droplet_summary")
    fig, axes = plt.subplots(2, 2, figsize=(PAGE, 4.75),
                             gridspec_kw=dict(wspace=.30, hspace=.48))
    ax = axes[0, 0]
    for length, color, marker in ((8192, "C0", "o"), (16384, "C3", "s")):
        rows = sorted((r for r in bulk if r["L"] == length), key=lambda r: r["q"])
        individuals = [r for r in bulk_runs if r["L"] == length]
        ax.scatter([r["q"] for r in individuals], [r["mean"] for r in individuals],
                   s=5, color=color, alpha=.3, linewidths=0)
        ax.errorbar([r["q"] for r in rows], [r["mean"] for r in rows],
                    [r["sem"] for r in rows], fmt=marker + "-", lw=.8,
                    mfc="white", color=color, ms=3, capsize=2, label=f"$L={length}$")
    shade_reference(ax)
    ax.set(xlabel="$q$", ylabel=r"bulk $\rho$", ylim=(.277, .295),
           title=r"Bulk, $2.5$--$5\times10^6$ steps")
    ax.legend(loc="lower right")

    ax = axes[0, 1]
    for observable, color, marker, name in (("envelope", "C0", "o", "profile threshold"),
                                          ("bond_activity", "C3", "s", "bond activity")):
        def select(row):
            return (row["L"] == 8192 and row["T"] == 1e7 and row["normal"] == "1-0"
                    and row["start"] == 5e6 and row["end"] == 1e7 and row["observable"] == observable)
        rows = sorted((r for r in points if select(r)), key=lambda r: r["q"])
        individuals = [r for r in runs if select(r)]
        ax.scatter([r["q"] for r in individuals], [r["v"] * 1e4 for r in individuals],
                   s=6, color=color, alpha=.3, linewidths=0)
        ax.errorbar([r["q"] for r in rows], [r["mean"] * 1e4 for r in rows],
                    [r["sem"] * 1e4 for r in rows], fmt=marker + "-", ms=3,
                    mfc="white", color=color, lw=.8, capsize=2, label=name)
    ax.axhline(0, color=".55", lw=.6)
    ax.set(xlabel="$q$", ylabel=r"normal $v$ ($10^{-4}$)",
           title=r"Axial front, $5$--$10\times10^6$ steps")
    ax.legend(loc="upper left", fontsize=5.8)

    ax = axes[1, 0]
    choices = ((8192, 2e6, 1e6, 2e6), (16384, 2e6, 1e6, 2e6),
               (8192, 5e6, 3.5e6, 5e6), (8192, 1e7, 5e6, 1e7))
    for normal, color, offset in (("1-0", "C0", -.14), ("2-1", "C1", 0), ("1-1", "C2", .14)):
        for i, choice in enumerate(choices):
            rows = [r for r in crossings if (r["L"], r["T"], r["start"], r["end"]) == choice
                    and r["normal"] == normal and r["observable"] == "envelope"]
            for row in rows:
                ax.errorbar(i + offset, row["estimate"],
                            [[row["estimate"] - row["low"]], [row["high"] - row["estimate"]]],
                            fmt="o", ms=3, color=color, capsize=2, lw=.8,
                            label=f"({normal.replace('-', ',')})" if i == 0 else None)
    shade_reference(ax, horizontal=True)
    ax.set(xticks=range(4), xticklabels=["1–2; 8192", "1–2; 16384", "3.5–5; 8192", "5–10; 8192"],
           xlabel=r"window ($10^6$ steps); $L$", ylabel=r"front crossing $q$",
           title="Time and direction sensitivity", ylim=(.045915, .045968))
    ax.tick_params(axis="x", labelsize=5.9)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.5f"))
    ax.legend(loc="upper center", ncol=3, fontsize=5.8, handlelength=1)

    ax = axes[1, 1]
    for kind, color, marker in (("disk", "C0", "o"), ("hole", "C3", "s")):
        rows = sorted((r for r in droplets if r["kind"] == kind and r["estimator"] == "trimmed_clipped"),
                      key=lambda r: r["q"])
        ax.errorbar([r["q"] for r in rows], [r["mean"] * 1000 for r in rows],
                    [r["sem"] * 1000 for r in rows], fmt=marker, color=color,
                    mfc="white", ms=3.5, capsize=2,
                    label="active disk" if kind == "disk" else "absorbing hole")
    shade_reference(ax)
    ax.axhline(0, color=".55", lw=.6)
    ax.set(xlabel="$q$", ylabel=r"mean $dR/dt$ ($10^{-3}$)",
           title=r"Droplets, $2$--$20\times10^4$ steps")
    ax.legend(loc="upper center", fontsize=6)
    for ax in (axes[0, 0], axes[0, 1], axes[1, 1]):
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.5f"))
        ax.xaxis.set_major_locator(MaxNLocator(4))
    axes[0, 1].set_xticks([.04593, .04594, .04595])
    axes[0, 0].set_xticks([.04590, .04594, .04598])
    for ax, letter in zip(axes.flat, "abcd"):
        label(ax, f"({letter})", x=.91, y=.14 if letter == "b" else .96)
    save(fig, "fig_flagship")
    sources = [Path(__file__), Path(__file__).with_name("common.py"),
               Path(ROOT) / "campaign/v2/results/coexistence_reference.json"]
    sources += [RESULTS / (name + ".csv") for name in
                ("bulk_summary", "bulk_runs", "front_points", "front_runs",
                 "front_crossings", "droplet_summary")]
    hashes = provenance(sources)
    for item in hashes:
        item["path"] = str(Path(item["path"]).relative_to(ROOT))
    write_json(RESULTS / "fig_flagship_provenance.json", dict(inputs=hashes))


if __name__ == "__main__":
    main()
