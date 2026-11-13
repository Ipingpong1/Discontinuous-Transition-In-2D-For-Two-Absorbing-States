"""Fig. phase: (a) transition line q*(p) at alpha = 1/2 from strip fronts, with the voter-point theory
q* = c alpha (1 - p); (b) interface density of the active phase just above q*."""
import glob
import os

import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path
from common import COL, RUNS, ROOT, FIGDATA, QSTAR05, REFERENCE_PATH, label, save
from campaign.v2.analysis.io import read_timeseries, provenance, write_csv, write_json
from campaign.v2.analysis.stats import window_mean, mean_summary
from campaign.v2.analysis.bulk import extend_absorbed
import fronts_util as fu

d = fu.load(["E1_front", "E4_voter", "E9_alpha", "E12_phase"])
pts = []
for p in sorted({k[0] for k in d if k[0] >= 0}):
    q, v, se, n = fu.curve(d, p, 0.5)
    if len(q) < 3:
        continue
    qs, sl = fu.crossing(q, v)
    if np.isfinite(qs):
        pts.append((p, qs, qs * sl))
pts = [pt for pt in pts if pt[0] != 0.99]
# p >= 0.99: strips run to 4 t_x (1.2e6 and 1.6e6 steps), velocities for t >= t_x, 4 seeds (campaign/strip_p099,
# 2026-10-31); from 1e5-step runs the crossings were 0.00128 and 0.00062
pts.append((0.99, 0.001282, 0.045))
pts.append((0.995, 0.00064, 0.03))
pts = [(p, QSTAR05 if p == .5 else q, slope) for p, q, slope in pts]
pts = np.array(sorted(pts))
for p, qs, red in pts:
    print(f"p={p:<6g} q*={qs:.6f}  q* dv/dq={red:.3f}  c={qs / (0.5 * (1 - p)) if p < 1 else np.nan:.3f}")


def plateau(path, tmin):
    meta, rows = read_timeseries(path)
    if not meta['done']:
        raise ValueError(f'Incomplete broad-survey ensemble: {path}')
    end = max(meta['T'], meta.get('t_end', meta['T']))
    if meta['absorbed']:
        rows = extend_absorbed(rows, end)
    value = window_mean(rows[:,0], rows[:,2], tmin, end)
    if not np.isfinite(value):
        raise ValueError(f'Incomplete averaging window: {path}')
    return dict(path=str(Path(path).resolve()), p=meta['p'], q=meta['q'],
                L=meta['Lx'], seed=meta['seed'], start=tmin, end=end,
                absorbed=meta['absorbed'], rho_mean=value,
                last_over_first=float(rows[-1,2]/np.interp(tmin,rows[:,0],rows[:,2])))


# rho just above q*: (p, experiment, name glob, tmin); the lowest q/q* whose runs stay stationary
STRIPS = os.path.join(os.path.dirname(os.path.dirname(RUNS)), "campaign", "strip_p099", "bulk")
PLATEAU = [(0.0, "E2_bulk", "a0.5_j2_q0.0185_L*_s*.txt", 1e4),
           (0.01, "E12_phase_bulk", "p0.01_a0.5_j2_q0.021_L2048_s*.txt", 2e4),
           (0.02, "E12_phase_bulk", "p0.02_a0.5_j2_q0.02332_L2048_s*.txt", 2e4),
           (0.05, "E12_phase_bulk", "p0.05_a0.5_j2_q0.029_L2048_s*.txt", 2e4),
           (0.1, "E12_phase_bulk", "p0.1_a0.5_j2_q0.037904_L2048_s*.txt", 2e4),
           (0.25, "E12_phase_bulk", "p0.25_a0.5_j2_q0.05047_L2048_s*.txt", 2e4),
           (0.5, "E13_bulk05", "p0.5_a0.5_j2_q0.046359_L2048_s*.txt", 2e4),
           (0.75, "E4_voter_bulk", "p0.75_a0.5_j2_q0.0286_L*_s*.txt", 2e4),
           (0.9, "E4_voter_bulk", "p0.9_a0.5_j2_q0.0122_L*_s*.txt", 2e4),
           (0.95, "E4_voter_bulk", "p0.95_a0.5_j2_q0.0065_L*_s*.txt", 2e4),
           (0.98, "E6_long", "p0.98_a0.5_j2_q0.0025704_L2048_s*.txt", 2e5),
           (0.99, None, [("E6_long", "p0.99_a0.5_j2_q0.0012827_L2048_s*.txt"),
                   (STRIPS, "p0.99_a0.5_j2_q0.00128_L2048_s*.txt")], 2e5),
           (0.995, None, [("E6_long_p995", "p0.995_a0.5_j2_q0.0006417_L2048_s*.txt"),
                    (STRIPS, "p0.995_a0.5_j2_q0.00064_L2048_s*.txt")], 2e5)]
# STRIPS: the bulk runs of campaign/strip_p099 (L 2048, seeds 3-4, to 1.2e6 / 1.5e6 steps), at the Table I thresholds
rp = []
run_records = []
summary_records = []
for p, exp, pat, tmin in PLATEAU:
    pairs = pat if isinstance(pat, list) else [(exp, pat)]
    files = sorted(f for e, pa in pairs for f in glob.glob(os.path.join(RUNS, e, "out", pa)))
    # The broad-survey table declares L>=2048; size selection is fixed and
    # independent of activity or whether a trajectory reaches consensus.
    files = [f for f in files if read_timeseries(f)[0]['Lx'] >= 2048]
    if not files:
        print("missing", exp, pat)
        continue
    vals = [plateau(f, tmin) for f in files]
    stats = mean_summary([v['rho_mean'] for v in vals], boot=0)
    rp.append((p, stats['mean'], stats['sem'], len(vals)))
    run_records.extend(vals)
    summary_records.append(dict(p=p, **stats, q_values=';'.join(str(q) for q in sorted({v['q'] for v in vals})),
                                L_values=';'.join(str(L) for L in sorted({v['L'] for v in vals}))))
    print(f"p={p}: rho_+ = {rp[-1][1]:.4f} from {len(vals)} runs ({exp}); drift last/first {np.mean([v['last_over_first'] for v in vals]):.3f}")
rp = np.array(rp)
result_dir = Path(ROOT)/'campaign/v2/results/phase'
write_csv(result_dir/'runs.csv', run_records)
write_csv(result_dir/'summary.csv', summary_records)
write_json(result_dir/'analysis.json', dict(
    averaging='Trapezoidal per-run time means, equal independent-run weights; no activity filter',
    uncertainty='One standard error over runs; broad-survey descriptive summaries, not limiting jumps',
    sizes='L>=2048 as declared in the broad parameter table',
    windows='Declared per-point starts to actual complete horizons; exact absorbing state extended to planned T',
    inputs=provenance([Path(v['path']) for v in run_records]),
    plotted_front_thresholds=[dict(p=float(p),qstar=float(q)) for p,q,_ in pts],
    front_inputs=provenance([Path(FIGDATA)/(name+'_fronts.csv') for name in
                            ('E1_front','E4_voter','E9_alpha','E12_phase')]+[Path(REFERENCE_PATH)]),
    legacy_near_voter=dict(points={'.99':.001282,'.995':.00064},
                           interpretation='Earlier nominal finite-window estimates; not recomputed by this figure',
                           inputs=provenance([Path(ROOT)/'campaign/strip_p099'/name for name in
                                             ('analyze.py','tables.md','fill_notes.py','make_runs.py')]))))

fig, (a, b) = plt.subplots(2, 1, figsize=(COL, 3.6), sharex=True, gridspec_kw=dict(hspace=0.08, height_ratios=[1.4, 1]))
pp = np.linspace(0.6, 1, 50)
a.plot(pp, 0.255255 * 0.5 * (1 - pp), "k--", lw=0.8,
       label=r"finite-time closure, $T=65536$")
a.plot(pts[pts[:, 0] > 0, 0], pts[pts[:, 0] > 0, 1], "o-", color="C0", ms=3.5, label="strip fronts")
a.plot(pts[pts[:, 0] == 0, 0], pts[pts[:, 0] == 0, 1], "o", color="C0", mfc="w", ms=3.5, label=r"$p=0$ reference")
a.set(ylabel=r"$q^*$", ylim=(0, 0.056))
a.text(0.02, 0.003, "coarsening", fontsize=6.5, color="0.3")
a.text(0.36, 0.051, "active", fontsize=6.5, color="0.3")
a.legend(fontsize=6, loc="lower center", bbox_to_anchor=(0.5, 0.0))
label(a, "(a)", x=0.03, y=0.95)
b.errorbar(rp[:, 0], rp[:, 1], rp[:, 2], fmt="s", color="C3", ms=3.5, capsize=1.2)
b.set(xlabel=r"voter fraction $p$", ylabel=r"$\rho$ at $q\gtrsim q^*$", ylim=(0, 0.4), xlim=(-0.02, 1.02))
label(b, "(b)", x=0.03, y=0.92)
save(fig, "fig_phase")
