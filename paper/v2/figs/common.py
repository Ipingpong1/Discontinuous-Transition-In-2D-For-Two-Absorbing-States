"""Shared paths and matplotlib style for the paper figures (PRE: 3.4 in column, 7.0 in page)."""
import os
import sys
import json
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(ROOT, "campaign", "v2", ".mplconfig"))
import matplotlib as mpl
RUNS = os.path.join(ROOT, "data", "runs")          # local mirror of $SCRATCH/dn/runs (*.txt only)
FIGDATA = os.path.join(ROOT, "data", "figdata")    # profile reductions made on Torch (extract_prof.py)
THEORY = os.path.join(ROOT, "campaign", "theory")
OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "campaign"))
sys.path.insert(0, ROOT)

COL, PAGE = 3.4, 7.0

REFERENCE_PATH = os.path.join(ROOT, 'campaign', 'v2', 'results', 'coexistence_reference.json')
with open(REFERENCE_PATH) as _reference_file:
    COEXISTENCE_REFERENCE = json.load(_reference_file)
QSTAR05 = float(os.environ.get('QSTAR05', COEXISTENCE_REFERENCE['qstar']))
QSTAR05_HALF_RANGE = float(COEXISTENCE_REFERENCE['systematic_half_range'])

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
    # Readers (paper builds and QA) only see a fully written artifact.
    for extension in ('pdf','png'):
        fd,temporary=tempfile.mkstemp(dir=OUT,prefix='.'+name+'.',suffix='.'+extension)
        os.close(fd)
        try:
            fig.savefig(temporary,**({'dpi':200} if extension=='png' else {}))
            os.replace(temporary,os.path.join(OUT,name+'.'+extension))
        finally:
            if os.path.exists(temporary):os.unlink(temporary)
