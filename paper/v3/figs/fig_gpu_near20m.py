"""All individual endpoint trajectories and paired late changes, with run means.

Generate inputs with python -m campaign.v3.analysis.film_endpoint from project root.
"""
from collections import defaultdict
import csv
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/dn-v3-matplotlib')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import ROOT, PAGE, label, save
from campaign.v2.analysis.io import write_json
from campaign.v3.analysis.film_endpoint import provenance

RESULT = Path(ROOT) / 'campaign/v3/results/film_endpoint'
QS = (.000640, .000642, .000644, .000648)
NAME = 'fig_gpu_near20m'


def table(name):
    with (RESULT / (name+'.csv')).open() as stream:
        return list(csv.DictReader(stream))


def main():
    curves = table('endpoint_curves')
    run_rows = table('endpoint_runs')
    colors = ['#542788', '#2b8cbe', '#238b45', '#d95f0e']
    fig, axes = plt.subplots(1, 2, figsize=(PAGE, 2.8),
                            gridspec_kw=dict(width_ratios=[1.15, 1], wspace=.32))
    zoom = axes[1].inset_axes([.39, .13, .58, .45])
    for q, color in zip(QS, colors):
        runs = defaultdict(list)
        for row in curves:
            if float(row['q']) == q:
                runs[row['run']].append(row)
        if len(runs) != 4:
            raise ValueError('Figure requires all four independent runs')
        records = [sorted(runs[key], key=lambda r: int(r['time'])) for key in sorted(runs)]
        grid = np.unique([int(row['time']) for rows in records for row in rows])
        values = []
        for rows in records:
            t = np.array([int(row['time']) for row in rows])
            rho = np.array([float(row['rho']) for row in rows])
            if t[-1] != 20_000_000:
                raise ValueError('Missing completed or absorbing extension')
            values.append(np.interp(grid, t, rho))
        values = np.array(values)
        selected = grid >= 1000
        for rho in values:
            axes[0].plot(grid[selected]/1e6, rho[selected], color=color, alpha=.25, lw=.5)
        axes[0].plot(grid[selected]/1e6, values.mean(0)[selected], color=color,
                     label=f'${q:.6f}$', lw=1.4)
        rows = sorted([r for r in run_rows if float(r['q']) == q], key=lambda r: r['run'])
        changes = np.array([float(r['drift']) for r in rows])
        if len(changes) != 4 or {r['run'] for r in rows} != set(runs):
            raise ValueError('Paired-change and trajectory cohorts differ')
        x = q * 1e6
        # Fixed horizontal offsets make all four independent runs visible.
        for ax in (axes[1], zoom):
            ax.scatter(x+np.linspace(-.25, .25, 4), changes, color=color,
                       s=13, marker='o', facecolors='none', linewidths=.8, zorder=3)
            ax.plot(x, changes.mean(), color=color, marker='D', ms=4,
                    markeredgecolor='white', markeredgewidth=.5, zorder=4)
    axes[0].set(xlim=(0, 20), ylim=(0, .27), xticks=[0, 5, 10, 15, 20],
                xlabel=r'$t\ (10^6)$', ylabel=r'interface density $\rho$',
                title=r'$p=0.995,\ L=8192$')
    axes[0].legend(title='$q$', title_fontsize=7, fontsize=6.3, ncol=2, loc='lower left')
    axes[0].text(.97, .98, 'thin: individual runs\nthick: mean', ha='right', va='top',
                 transform=axes[0].transAxes, fontsize=6)
    axes[1].axhline(0, color='.5', lw=.7)
    axes[1].set(xlim=(639, 649), xticks=np.array(QS)*1e6, ylim=(-.122, .015),
                xlabel=r'$q\ (10^{-6})$', ylabel=r'paired late change $\Delta\rho$')
    axes[1].legend(handles=[Line2D([], [], color='.3', marker='o', markerfacecolor='none',
                                  ls='', ms=3.5, label='individual run'),
                            Line2D([], [], color='.3', marker='D', ls='', ms=3.5,
                                   label='four-run mean')], fontsize=5.8,
                   loc='lower center', bbox_to_anchor=(.5, 1.02), ncol=2,
                   borderaxespad=0)
    zoom.axhline(0, color='.5', lw=.6)
    zoom.set(xlim=(641, 649), xticks=[642, 644, 648], ylim=(-.011, .004),
             yticks=[-.01, -.005, 0], title='enlarged vertical scale')
    zoom.tick_params(labelsize=5.5)
    zoom.set_title('enlarged vertical scale', fontsize=6, pad=3)
    for ax, letter in zip(axes, 'ab'):
        label(ax, f'({letter})', x=.03)
    save(fig, NAME)
    write_json(RESULT/(NAME+'_provenance.json'), dict(
        inputs=provenance([RESULT/'endpoint_runs.csv', RESULT/'endpoint_curves.csv', Path(__file__)]),
        display='All 16 trajectories, four individual paired changes per q, equal-run means; no confidence intervals.',
        inset='Same q=0.000642, 0.000644 and 0.000648 points on an enlarged vertical scale.'))
    print(NAME)


if __name__ == '__main__':
    main()
