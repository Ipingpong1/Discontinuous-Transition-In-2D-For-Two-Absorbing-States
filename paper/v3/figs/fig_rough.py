"""Physical wetting-film edges and thickness distributions for v3.

Generate inputs with python -m campaign.v3.analysis.film_endpoint from project root.
The invalid free-band comparator is removed because its original snapshots are
not in the local archive. No amplitudes are reconstructed from its reductions.
"""
import csv
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/dn-v3-matplotlib')
import numpy as np
import matplotlib.pyplot as plt
from common import ROOT, COL, QSTAR05, label, save
from campaign.v2.analysis.io import write_json
from campaign.v3.analysis.film_endpoint import provenance

RESULT = Path(ROOT) / 'campaign/v3/results/film_endpoint'


def read(name):
    with (RESULT/(name+'.csv')).open() as stream:
        return list(csv.DictReader(stream))


def main():
    correlations = read('wall_correlations')
    histograms = read('wall_histograms')
    fig, (a, b) = plt.subplots(1, 2, figsize=(2*COL, 2.2), gridspec_kw=dict(wspace=.3))
    for lx, marker, color in ((256, 'v', 'C2'), (1024, 'o', 'C0'), (4096, 's', 'C3')):
        rows = [r for r in correlations if int(r['Lx']) == lx]
        r = np.array([int(row['r']) for row in rows])
        c = np.array([float(row['C']) for row in rows])
        a.plot(r, c, marker+'-', color=color, ms=3, lw=.8, label=f'$L_x={lx}$')
    r = np.array([1.5, 30])
    a.plot(r, 25*r, 'k--', lw=.6)
    a.text(8, 120, r'$r^1$', fontsize=7)
    a.set(xscale='log', yscale='log', xlabel=r'$r$',
          ylabel=r'$\langle[h(x+r)-h(x)]^2\rangle$')
    a.legend(fontsize=6, title='film edges', title_fontsize=6, loc='lower right')
    label(a, '(a)')
    cmap = plt.get_cmap('plasma')
    qs = (.043605, .044982, .045441, .0456705)
    for i, q in enumerate(qs):
        rows = [r for r in histograms if float(r['q']) == q]
        w = np.array([int(row['w']) for row in rows])
        p = np.array([float(row['P']) for row in rows])
        if np.any(np.diff(w) != 1):
            raise ValueError('Histogram smoothing requires consecutive integer bins')
        smoothed = np.convolve(p, np.ones(5)/5, 'same')
        b.plot(w, smoothed, color=cmap(.85*i/(len(qs)-1)), lw=.9,
               label=f'${1-q/QSTAR05:.3f}$')
    b.set(xlabel=r'local layer width $w$', ylabel=r'$P(w)$', xlim=(0, 600),
          yscale='log', ylim=(1e-5, .05))
    b.legend(fontsize=6, title=r'$\delta$', title_fontsize=6.5)
    label(b, '(b)', x=.9)
    save(fig, 'fig_rough')
    write_json(RESULT/'fig_rough_provenance.json', dict(
        inputs=provenance([RESULT/'wall_correlations.csv', RESULT/'wall_histograms.csv', Path(__file__)]),
        curves='Both physical film edges, equal-run means of archived per-run reductions.',
        histogram_display='Five-integer-bin moving average of run-averaged probabilities; raw probabilities retained.',
        removed='Free-band comparator combined physical edges with artificial cuts through active bulk; original snapshots are unavailable locally.'))


if __name__ == '__main__':
    main()
