"""Shared paths and matplotlib style for the paper figures (PRE: 3.4 in column, 7.0 in page)."""
import os
import sys

import matplotlib as mpl

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RUNS = os.path.join(ROOT, "data", "runs")          # local mirror of $SCRATCH/dn/runs (*.txt only)
FIGDATA = os.path.join(ROOT, "data", "figdata")    # profile reductions made on Torch (extract_prof.py)
THEORY = os.path.join(ROOT, "campaign", "theory")
OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "campaign"))
sys.path.insert(0, ROOT)

COL, PAGE = 3.4, 7.0

mpl.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8, "axes.labelsize": 8,
    "axes.titlesize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "lines.linewidth": 1.0, "lines.markersize": 3.5, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02, "pdf.fonttype": 42,
})


def label(ax, s, x=0.03, y=0.95):
    ax.text(x, y, s, transform=ax.transAxes, ha="left", va="top", fontsize=8)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200)
