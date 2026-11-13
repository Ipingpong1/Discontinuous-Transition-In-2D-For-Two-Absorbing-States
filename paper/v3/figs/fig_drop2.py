"""Two-droplet trajectories, radius-dependent curvature, and area sensitivity.

Run campaign.v3.analysis.fronts_review before this figure. The analysis retains
its original trimmed/rectified radius and exports full signed-area sensitivity.
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import ROOT, PAGE, label, save
from campaign.v2.analysis.io import provenance, write_json

RESULTS = Path(ROOT) / "campaign/v3/results/fronts_review"


def read_table(name):
    with (RESULTS / (name + ".csv")).open() as stream:
        return [{key: value if key in ("run", "kind", "estimator") else float(value)
                 for key, value in row.items()} for row in csv.DictReader(stream)]


def main():
    curves = read_table("droplet_curves")
    predictions = read_table("droplet_predictions")
    summary = read_table("droplet_summary")
    qs = sorted({r["q"] for r in curves})
    colors = dict(zip(qs, ("C0", "C9", "C1", "C3")))
    fig, axes = plt.subplots(1, 3, figsize=(PAGE, 2.25),
                             gridspec_kw=dict(width_ratios=(1, 1, 1.1), wspace=.42))
    for ax, kind in zip(axes[:2], ("disk", "hole")):
        for q in qs:
            runs = sorted({r["run"] for r in curves if r["q"] == q and r["kind"] == kind})
            radii = []
            for run in runs:
                values = sorted((r for r in curves if r["run"] == run), key=lambda r: r["time"])
                time = np.array([r["time"] for r in values])
                radius = np.array([r["radius"] for r in values])
                radii.append(radius)
                ax.plot(time / 1000, radius, color=colors[q], lw=.4, alpha=.35)
            ax.plot(time / 1000, np.mean(radii, axis=0), color=colors[q], lw=1.2,
                    label=f"{q:.7f}")
            model = sorted((r for r in predictions if r["q"] == q and r["kind"] == kind),
                           key=lambda r: r["time"])
            ax.plot([r["time"] / 1000 for r in model], [r["radius"] for r in model],
                    color=colors[q], ls="--", lw=.9)
        ax.axhline(1000, color=".6", ls=":", lw=.5)
        ax.set(xlabel=r"$t$ ($10^3$ steps)", xlim=(0, 200), ylim=(450, 1520))
    axes[0].set(ylabel=r"equivalent radius $R$", title="Active disk")
    axes[1].set(title="Absorbing hole")
    axes[0].legend(fontsize=5, title="$q$", title_fontsize=5.8,
                   ncol=2, loc="upper left", handlelength=1, columnspacing=.7)
    for kind, marker in (("disk", "o"), ("hole", "s")):
        for estimator, fill in (("trimmed_clipped", True), ("full_signed", False)):
            values = sorted((r for r in summary if r["kind"] == kind
                             and r["estimator"] == estimator), key=lambda r: r["q"])
            offset = -4e-6 if estimator == "trimmed_clipped" else 4e-6
            x = np.array(qs) + offset
            y = np.array([r["mean"] * 1000 for r in values])
            err = np.array([r["sem"] * 1000 for r in values])
            axes[2].errorbar(x, y, err, fmt=marker, color="C0" if kind == "disk" else "C3",
                             mfc=("C0" if kind == "disk" else "C3") if fill else "white",
                             ms=3, capsize=1.4, lw=.7,
                             label=kind if fill else None)
    axes[2].axhline(0, color=".6", lw=.6)
    axes[2].set(xticks=qs, xticklabels=[f"{q:.7f}" for q in qs],
                xlabel="$q$", ylabel=r"mean $dR/dt$ ($10^{-3}$)", title="Area-estimator sensitivity")
    axes[2].tick_params(axis="x", labelrotation=45, labelsize=5.5)
    axes[2].legend(loc="upper center", fontsize=5.8, handlelength=1)
    for ax, letter in zip(axes, "abc"):
        label(ax, f"({letter})", x=.05 if letter == "b" else .86,
              y=.12 if letter == "b" else .96)
    save(fig, "fig_drop2")
    sources = [Path(__file__), Path(__file__).with_name("common.py")]
    sources += [RESULTS / (name + ".csv") for name in
                ("droplet_curves", "droplet_predictions", "droplet_summary")]
    hashes = provenance(sources)
    for item in hashes:
        item["path"] = str(Path(item["path"]).relative_to(ROOT))
    write_json(RESULTS / "fig_drop2_provenance.json", dict(inputs=hashes))


if __name__ == "__main__":
    main()
