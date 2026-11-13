"""Compact, machine-regenerated status and complete-ensemble results for the paper."""
import csv
import json
from pathlib import Path
from .io import ROOT, write_json
from collections import defaultdict


def table(path):
    if not path.exists():return []
    with path.open() as stream:return list(csv.DictReader(stream))


def coexistence_audit(folder,status,fits):
    """Descriptive systematic envelopes, never a pooled confidence interval."""
    windows={2e6:{(1e6,2e6)},5e6:{(2e6,3.5e6),(3.5e6,5e6),(2e6,5e6)},1e7:{(2e6,5e6),(5e6,1e7)}}
    rows=[{k:v for k,v in r.items() if k not in ('coefficients','bootstrap_samples')}
          for r in fits if r['p']==.5 and r['model']=='bracket' and r['estimator']=='envelope'
          and (r['t0'],r['t1']) in windows.get(r['T'],set())]
    groups=defaultdict(list)
    for row in rows:groups[(row['L'],row['T'],row['t0'],row['t1'])].append(row)
    summaries=[]
    for key,group in sorted(groups.items()):
        central=[r['crossing'] for r in group if r['threshold']==.3]
        all_points=[r['crossing'] for r in group]
        finite=[r for r in group if r['ci_low'] is not None and r['ci_high'] is not None]
        summaries.append(dict(zip(('L','T','t0','t1'),key))|dict(
            n_direction_threshold_crossings=len(group),n_expected=3 if key[1]==1e7 else 9,
            all_directions_thresholds_available=len(group)==(3 if key[1]==1e7 else 9),
            threshold03_direction_points=central,
            point_estimate_envelope=[min(all_points),max(all_points)],
            conditional_CI_envelope=[min(r['ci_low'] for r in finite),max(r['ci_high'] for r in finite)] if finite else None,
            interval_valid_fractions=[r['bootstrap_valid']/r['bootstrap_requested'] for r in group]))
    required=[]
    for L,T in [(8192,2e6),(16384,2e6),(8192,5e6),(8192,1e7)]:
        ensembles=[e for e in status['ensembles'] if e['p']==.5 and e['mode']=='band' and e['Lx']==L and e['T']==T]
        required.append(dict(L=L,T=T,n_planned_groups=len(ensembles),
            n_complete_groups=sum(e['ensemble_complete'] for e in ensembles),
            all_complete=bool(ensembles) and all(e['ensemble_complete'] for e in ensembles)))
    all_points=[r['crossing'] for r in rows]
    finite=[r for r in rows if r['ci_low'] is not None and r['ci_high'] is not None]
    point_rows=table(folder/'front_points.csv');grid_status=[]
    for L,T in [(4096,2e6),(8192,2e6),(16384,2e6),(8192,5e6),(8192,1e7)]:
        for t0,t1 in sorted(windows[T]):
            for normal in (('1-0',) if T==1e7 else ('1-0','2-1','1-1')):
                a,b=map(int,normal.split('-'))
                declared=[e for e in status['ensembles'] if e['p']==.5 and e['mode']=='band' and e['Lx']==L and e['T']==T and (e['a'],e['b'])==(a,b)]
                finished=bool(declared) and all(e['ensemble_complete'] for e in declared)
                for threshold in (.2,.3,.5):
                    pp=[r for r in point_rows if float(r['p'])==.5 and float(r['L'])==L and float(r['T'])==T and r['ab']==normal
                        and float(r['threshold'])==threshold and r['estimator']=='envelope' and float(r['t0'])==t0 and float(r['t1'])==t1]
                    pp=sorted(pp,key=lambda r:float(r['q']));velocities=[float(r['v']) for r in pp]
                    fit=next((r for r in rows if (r['L'],r['T'],r['t0'],r['t1'],r['ab'],r['threshold'])==(L,T,t0,t1,normal,threshold)),None)
                    side=None
                    changes=sum(x*y<0 for x,y in zip(velocities[:-1],velocities[1:]))
                    if changes>1:state='multiple_sampled_mean_sign_changes'
                    elif fit:state='bracketed' if fit['ci_low'] is not None else 'bracketed_interval_withheld'
                    elif not finished:state='pending_declared_cohorts'
                    elif len(pp)<2:state='insufficient_geometry_valid_q_points'
                    elif all(v>0 for v in velocities):state='all_sampled_mean_velocities_positive';side='below smallest usable q under monotonicity; not a confidence bound'
                    elif all(v<0 for v in velocities):state='all_sampled_mean_velocities_negative';side='above largest usable q under monotonicity; not a confidence bound'
                    else:state='no_accepted_crossing'
                    grid_status.append(dict(L=L,T=T,t0=t0,t1=t1,normal=normal,threshold=threshold,state=state,
                        all_declared_q_groups_complete=finished,n_declared_q_groups=len(declared),n_usable_q_points=len(pp),n_strict_sign_changes=changes,
                        usable_q=[float(r['q']) for r in pp],mean_velocities=velocities,unbracketed_side=side))
    write_json(folder/'coexistence_audit.json',dict(
        purpose='Late-window operational crossings and separate conservative diagnostic envelopes; does not change shared qstar',
        required_long_cohorts=required,all_required_cohorts_complete=all(r['all_complete'] for r in required),
        point_estimate_envelope=[min(all_points),max(all_points)] if all_points else None,
        conditional_CI_envelope=[min(r['ci_low'] for r in finite),max(r['ci_high'] for r in finite)] if finite else None,
        uncertainty_scope='Envelopes across correlated analyses are sensitivity summaries, not pooled confidence intervals or probability distributions.',
        missing_crossing_warning='A numerical envelope of existing roots does not cover unbracketed or geometrically unavailable configurations. Inspect grid_status before choosing an operational reference.',
        grid_status=grid_status,summaries=summaries,crossings=rows))


def main():
    folder=ROOT/'campaign/v2/results/gpu'
    status=json.loads((folder/'status.json').read_text())
    out=[f"Completed tasks: {status['n_complete']}/{status['n_planned']}; eligible fixed groups: {status['n_complete_ensembles']}/{status['n_ensembles']}.",
         'Only groups containing every declared independent seed enter the results below.']
    for e in status['ensembles']:
        if e['ensemble_complete']:
            out.append(f"COMPLETE p={e['p']}, q={e['q']}, {e['mode']}, L={e['Lx']}x{e['Ly']}, T={e['T']}, normal=({e['a']},{e['b']}), n={e['n_planned']}")
    fits=json.loads((folder/'front_crossings.json').read_text())
    coexistence_audit(folder,status,fits)
    selected=[r for r in fits if r['p']==.5 and r['estimator']=='envelope' and r['model']=='bracket' and r['t1'] is not None]
    if selected:
        out.extend(['','Full-window local front crossings (conditional95% run-bootstrap intervals):'])
    for r in selected:
        out.append(f"L={r['L']}, T={r['T']}, normal={r['ab']}, threshold={r['threshold']}, t={r['t0']:g}--{r['t1']:g}: q*={r['crossing']:.9f}; CI=[{r['ci_low']},{r['ci_high']}], bootstrap valid={r['bootstrap_valid']}/{r['bootstrap_requested']}")
    contrast=json.loads((folder/'orientation_differences.json').read_text())
    contrast=[r for r in contrast if r['threshold']==.3 and r['t0']==1e6 and r['t1']==2e6 and r['estimator']=='envelope']
    if contrast:out.extend(['','Late-window direction contrasts, with simultaneous95% intervals over the three pairs:'])
    for r in contrast:
        out.append(f"L={r['L']}, T={r['T']}: {r['contrast']}: {r['difference']:.9g}; interval=[{r.get('simultaneous95_low')},{r.get('simultaneous95_high')}]; spread upper95={r.get('bootstrap_upper95')}")
    paired_path=folder/'front5m_paired_crossings.json'
    paired=json.loads(paired_path.read_text()) if paired_path.exists() else []
    paired=[r for r in paired if r['estimator']=='envelope']
    if paired:out.extend(['','Independent5m cohorts: paired late-minus-early crossing differences (2--3.5m versus3.5--5m):'])
    for r in paired:
        out.append(f"L={r['L']}, T={r['T']}, normal={r['ab']}, threshold={r['threshold']}: early={r['q_early']:.9f}, late={r['q_late']:.9f}, difference={r['difference']:.9g}, CI=[{r['ci95_low']},{r['ci95_high']}], valid={r['bootstrap_valid']}/{r['bootstrap_requested']}")
    bulk=table(folder/'bulk_window_summary.csv')
    if bulk:out.extend(['','Unconditional bulk means and within-window drift (each with conditional95% run-bootstrap interval):'])
    for r in bulk:
        out.append(f"p={r['p']}, q={r['q']}, L={r['Lx']}, T={r['T']}, window={r['start']}--{r['end']}: rho={r['rho_mean']} [{r['rho_ci_low']},{r['rho_ci_high']}]; drift={r['rho_drift_mean']} [{r['rho_drift_ci_low']},{r['rho_drift_ci_high']}]; chi={r['chi_mean']} [{r['chi_ci_low']},{r['chi_ci_high']}]")
    (folder/'review_summary.txt').write_text('\n'.join(out)+'\n')
    print(out[0]);print(folder/'review_summary.txt')

if __name__=='__main__':main()
