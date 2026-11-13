"""Reduce saved endpoint runs and valid wall-edge measurements for v3.

Run from the project root with ``python -m campaign.v3.analysis.film_endpoint``.
Only campaign/v3/results/film_endpoint is written. No simulations are launched.
"""
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from campaign.v2.analysis.bulk import extend_absorbed
from campaign.v2.analysis.io import read_timeseries, write_csv, write_json
from campaign.v2.analysis.stats import window_mean

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'campaign/v3/results/film_endpoint'
QS = (.000640, .000642, .000644, .000648)


def provenance(paths):
    return [dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size,
                 sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for p in sorted(set(map(Path, paths)))]


def endpoint():
    status_path = ROOT / 'campaign/v2/results/gpu/status.json'
    status = json.loads(status_path.read_text())
    cohorts = [e for e in status['ensembles'] if e['p'] == .995
               and e['mode'] == 'bulk' and e['Lx'] == e['Ly'] == 8192
               and e['T'] == 20_000_000 and e['q'] in QS]
    if len(cohorts) != 4 or {e['q'] for e in cohorts} != set(QS):
        raise ValueError('Expected the four declared endpoint cohorts')
    curves, runs, summaries, inputs = [], [], [], [status_path]
    for cohort in sorted(cohorts, key=lambda e: e['q']):
        if not (cohort['ensemble_complete'] and cohort['n_complete'] == 4
                and len(cohort['task_ids']) == 4):
            raise ValueError('Incomplete endpoint cohort')
        values = []
        for ident in sorted(cohort['task_ids']):
            path = ROOT / 'data/v2/gpu/runs/near_voter' / ident / 'out.txt'
            meta, series = read_timeseries(path)
            if not meta['done'] or any(meta[k] != cohort[k]
                    for k in ('p', 'alpha', 'q', 'Lx', 'Ly', 'T')):
                raise ValueError(f'Raw data disagree with cohort declaration: {ident}')
            if meta['absorbed']:
                series = extend_absorbed(series, meta['T'])
            t, rho = series[:, 0], series[:, 2]
            first = window_mean(t, rho, 10_000_000, 15_000_000)
            second = window_mean(t, rho, 15_000_000, 20_000_000)
            row = dict(run=ident, seed=meta['seed'], q=meta['q'],
                       absorbed=meta['absorbed'], t_end=meta['t_end'],
                       mean_first=first, mean_second=second,
                       mean_late=window_mean(t, rho, 10_000_000, 20_000_000),
                       drift=second-first)
            runs.append(row)
            values.append(row)
            curves.extend(dict(run=ident, q=meta['q'], time=int(tt), rho=rr)
                          for tt, rr in zip(t, rho))
            inputs.extend([path, path.with_name('done.json')])
        changes = np.array([v['drift'] for v in values])
        mean = float(changes.mean())
        # This is a transparent normal-reference sensitivity check, not a
        # distribution-free confidence statement about four realizations.
        half = 3.182446305 * float(changes.std(ddof=1)) / 2
        summaries.append(dict(q=cohort['q'], n=4,
            n_absorbed=sum(v['absorbed'] for v in values),
            mean_late=float(np.mean([v['mean_late'] for v in values])),
            mean_drift=mean, min_drift=float(changes.min()), max_drift=float(changes.max()),
            n_negative=int(np.sum(changes < 0)),
            normal_reference_t_low=mean-half, normal_reference_t_high=mean+half))
    write_csv(OUT / 'endpoint_runs.csv', runs)
    write_csv(OUT / 'endpoint_curves.csv', curves)
    write_csv(OUT / 'endpoint_summary.csv', summaries)
    return inputs


def wall_edges():
    inputs, curves, slopes, histograms = [], [], [], []
    r2 = ROOT / 'data/figdata/R2'
    for lx in (256, 1024, 4096):
        paths = sorted(r2.glob(f'rough_p0.5_a0.5_j2_q0.045441_Lx{lx}_Ly2048_s*.csv'))
        if len(paths) != (2 if lx == 4096 else 4):
            raise ValueError(f'Unexpected wall-edge run count for Lx={lx}')
        coordinates, samples = None, []
        for path in paths:
            with path.open() as stream:
                rows = list(csv.DictReader(stream))
            r = np.array([int(row['r']) for row in rows])
            if coordinates is not None and not np.array_equal(r, coordinates):
                raise ValueError('Inconsistent wall-edge lag grids')
            coordinates = r
            samples.append(np.array([(float(row['C_lo'])+float(row['C_hi']))/2 for row in rows]))
        mean = np.mean(samples, axis=0)
        fit = (coordinates >= 2) & (coordinates <= 50)
        slope = float(np.polyfit(np.log(coordinates[fit]), np.log(mean[fit]), 1)[0]/2)
        curves.extend(dict(Lx=lx, Ly=2048, q=.045441, n=len(paths), r=int(r), C=float(c))
                      for r, c in zip(coordinates, mean))
        slopes.append(dict(Lx=lx, Ly=2048, q=.045441, n=len(paths),
                           r_min=2, r_max=50, zeta_effective=slope))
        inputs.extend(paths)
    for q in (.043605, .044982, .045441, .0456705):
        paths = sorted(r2.glob(f'whist_p0.5_a0.5_j2_q{q:g}_Lx1024_Ly2048_s*.csv'))
        if not paths:
            raise ValueError(f'Missing thickness histograms at q={q}')
        probabilities = defaultdict(float)
        for path in paths:
            with path.open() as stream:
                for row in csv.DictReader(stream):
                    probabilities[int(row['w'])] += float(row['P']) / len(paths)
        histograms.extend(dict(q=q, Lx=1024, Ly=2048, n=len(paths), w=w, P=p)
                          for w, p in sorted(probabilities.items()))
        inputs.extend(paths)
    write_csv(OUT / 'wall_correlations.csv', curves)
    write_csv(OUT / 'wall_slopes.csv', slopes)
    write_csv(OUT / 'wall_histograms.csv', histograms)
    return inputs


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = endpoint() + wall_edges()
    inputs += [Path(__file__), ROOT/'campaign/v2/analysis/io.py',
               ROOT/'campaign/v2/analysis/stats.py', ROOT/'campaign/v2/analysis/bulk.py']
    write_json(OUT/'analysis.json', dict(
        endpoint_estimator='All four fixed-manifest runs per q; paired trapezoidal means on 15–20M minus 10–15M steps; known absorbing states extended to the declared horizon.',
        endpoint_figure='Every individual trajectory and paired change, plus equal-run means. No confidence bands or significance declarations.',
        t_sensitivity='t3 normal-reference mean intervals are descriptive sensitivity checks; distributional assumptions cannot be verified with four independent runs.',
        wall_estimator='Means of both physical film-edge height-difference correlations, averaged equally over stored run reductions; effective log–log slope on recorded lags 2–50.',
        removed_comparator='The old free-band reduction combined physical exterior edges with artificial cuts through active bulk. Original R2_orient .snap files are absent locally, so the comparator is removed. Existing wall-edge reductions remain valid.',
        snapshot_check=dict(local_pattern='data/runs/R2_orient/out/p0.5_a0.5_q0.04596_ab1-0_L4096_s*.snap',
                            available=bool(list((ROOT/'data/runs/R2_orient/out').glob('p0.5_a0.5_q0.04596_ab1-0_L4096_s*.snap')))),
        inputs=provenance(inputs)))
    print(f'Wrote v3 endpoint and wall summaries to {OUT}')


if __name__ == '__main__':
    main()
