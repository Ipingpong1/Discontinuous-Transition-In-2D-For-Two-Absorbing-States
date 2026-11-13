"""Reduce fixed GPU ensembles; publication tables never use partial seed groups.

python3 -m campaign.v2.analysis.gpu_campaign --boot 1000
All declared tasks appear in eligibility/status outputs. Header, terminal record,
profile timeline and done.json must agree before a task is complete. Final
front summaries additionally require a valid slope for every declared seed at
the given parameter/window; geometry loss is reported, not silently filtered.
"""
from __future__ import annotations
import argparse
import math
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from . import bulk, fronts
from .io import (ROOT, complete_ensembles, completion_status, manifest_tasks,
                 profile_path, provenance, read_profile, read_timeseries,
                 scientific_group, write_csv, write_json)
from .stats import mean_summary, window_mean, group_rng
from .front_stationarity import paired_contrasts

MANIFESTS = [ROOT/'campaign/v2/hpc/gpu'/name for name in
    ['production_manifest_profiles50k.json','production_manifest_largeL.json',
     'production_manifest_near_large.json','production_manifest_orientation5m.json',
     'production_manifest_near_midpoint.json','production_manifest_axial10m.json']]
WINDOWS = [(1e6,2e6),(2e6,5e6),(5e6,1e7),(1e7,2e7)]
FRONT_WINDOWS = [(1e5,5e5),(5e5,1e6),(1e6,2e6),(2e6,3.5e6),(3.5e6,5e6),(2e6,5e6),(5e6,1e7)]


def write_table(path, rows):
    if rows:
        write_csv(path, rows)
    else:
        path.unlink(missing_ok=True)


def scheduled_times(task,end):
    """Required ppd2 plus manifest periodic samples; extra finer records allowed."""
    times={0}
    for decade in range(11):
        for half in range(2):
            t=math.floor(10.**(decade+half/2.)+.5)
            if t<=end:times.add(t)
    every=int(task['every'])
    if every>0:times.update(range(every,int(end)+1,every))
    if task['T']<=end:times.add(int(task['T']))
    return times


def task_completion(task, path):
    if not path.exists():
        return dict(complete=False,reasons=['not_copied_or_not_started']), None, None
    try:
        meta, ts = read_timeseries(path)
        errors = list(meta.get('completion_errors', []))
        missing=sorted(scheduled_times(task,meta.get('t_end',float(ts[-1,0])))-set(ts[:,0]))
        if missing:errors.append('missing_scheduled_text_records')
        for key, field in [('Lx','Lx'),('Ly','Ly'),('T','T'),('seed','seed'),('p','p'),('q','q'),('alpha','alpha')]:
            expected = float(task[field]) if key in ('p','q','alpha') else int(task[field])
            if meta.get(key) != expected:
                errors.append('manifest_mismatch_'+key)
        if task['mode']=='band':
            if not meta.get('profile_dir') or (meta.get('a'),meta.get('b')) != tuple(task['orientation']):
                errors.append('manifest_mismatch_profile_direction')
        else:
            expected_init={'bulk':('random',.5),'strip':('strip',.5,.5),'wall':('strip',1.,.5)}.get(task['mode'])
            fields=meta.get('init','').split(':')
            try: actual_init=(fields[0],)+tuple(float(v) for v in fields[1:])
            except ValueError: actual_init=None
            if expected_init is not None and actual_init != expected_init:
                errors.append('manifest_mismatch_initial_protocol')
            if meta.get('profile_dir'):
                errors.append('unexpected_directed_profile')
        if task['profile']:
            bins = meta['Lx'] if meta.get('profile_dir') else meta['Ly']
            profile, trailing = read_profile(profile_path(path), bins)
            status = completion_status(path, meta, ts, profile, trailing)
            errors += status['reasons']
        return dict(complete=not errors,reasons=sorted(set(errors)),
                    t_end=meta.get('t_end'),absorbed=meta.get('absorbed',False),
                    header_claim_complete=bool(meta.get('header_done',False)),
                    n_missing_scheduled=len(missing),missing_scheduled_times=missing,
                    latest_record=float(ts[-1,0])), meta, ts
    except (ValueError, OSError, KeyError, IndexError) as error:
        return dict(complete=False,reasons=[type(error).__name__+': '+str(error)]),None,None


def reduce_campaign(root, manifests, out, boot=1000, hashes=True):
    tasks = manifest_tasks(manifests)
    status, paths, loaded = {}, {}, {}
    for ident, task in tasks.items():
        path = root/task['family']/ident/'out.txt'
        paths[ident] = path
        status[ident], meta, ts = task_completion(task,path)
        if status[ident]['complete']:
            loaded[ident] = (meta,ts)
    ensembles = complete_ensembles(tasks,status)
    eligible_ids = {i for e in ensembles if e['ensemble_complete'] for i in e['task_ids']}
    bulk_runs, bulk_windows, bulk_curves = [], [], []
    front_runs, front_records = [], []
    inputs = list(manifests)+[Path(__file__).parent/name for name in
        ('gpu_campaign.py','front_stationarity.py','fronts.py','bulk.py','io.py','stats.py')]
    for ident,(meta,ts) in loaded.items():
        task, path = tasks[ident], paths[ident]
        inputs.append(path)
        for marker in ('done.json','import_record.json'):
            if (path.parent/marker).exists():inputs.append(path.parent/marker)
        if task['profile']:
            inputs.append(profile_path(path))
        if task['mode']=='bulk':
            row = bulk.reduce_run(path)
            row['ensemble_complete'] = ident in eligible_ids
            bulk_runs.append(row)
            series = bulk.extend_absorbed(ts,task['T']) if meta['absorbed'] else ts
            for start,end in WINDOWS:
                if end>task['T']:
                    continue
                mid = (start+end)/2
                bulk_windows.append(dict(run=ident,p=float(task['p']),q=float(task['q']),
                    Lx=task['Lx'],Ly=task['Ly'],T=task['T'],start=start,end=end,
                    ensemble_complete=ident in eligible_ids,
                    rho=window_mean(series[:,0],series[:,2],start,end),
                    rho_first=window_mean(series[:,0],series[:,2],start,mid),
                    rho_second=window_mean(series[:,0],series[:,2],mid,end),
                    rho_drift=window_mean(series[:,0],series[:,2],mid,end)-window_mean(series[:,0],series[:,2],start,mid),
                    chi=task['Lx']*task['Ly']*window_mean(series[:,0],(2*series[:,1]-1)**2,start,end)))
            if ident in eligible_ids:
                grid = np.unique(np.r_[np.geomspace(1,task['T'],180),task['T']])
                for t,rho in zip(grid,np.interp(grid,series[:,0],series[:,2])):
                    bulk_curves.append(dict(run=ident,p=float(task['p']),q=float(task['q']),
                        Lx=task['Lx'],Ly=task['Ly'],T=task['T'],time=t,rho=rho))
        else:
            windows = [w for w in FRONT_WINDOWS if w[1]<=task['T']]
            rows, records = fronts.reduce_run(path,thresholds=(.2,.3,.5),
                tmins=(1e5,5e5,1e6),windows=windows)
            for r in rows:
                r['ensemble_complete'] = ident in eligible_ids
            for r in records:
                r.update(p=float(task['p']),q=float(task['q']),L=task['Lx'],Ly=task['Ly'],
                    T=task['T'],ab=f"{(task.get('orientation') or [1,0])[0]}-{(task.get('orientation') or [1,0])[1]}",
                    ensemble_complete=ident in eligible_ids)
            if meta['absorbed'] and records[-1]['time']<task['T']:
                last = dict(records[-1]);last.update(time=task['T'],envelope=0.,component=0.,
                        hole_fraction=0.,interior_ratio=0.,touches_edge=False)
                records.append(last)
            diagnostic_t=np.array([r['time'] for r in records])
            diagnostic_i=np.array([r['interior_ratio'] for r in records])
            diagnostic_h=np.array([r['hole_fraction'] for r in records])
            for r in rows:
                r['interior_ratio_start']=float(np.interp(r['t0'],diagnostic_t,diagnostic_i))
                r['interior_ratio_end']=float(np.interp(min(r['t1'],r['T']),diagnostic_t,diagnostic_i))
                r['interior_ratio_change']=r['interior_ratio_end']-r['interior_ratio_start']
                inside=(diagnostic_t>=r['t0'])&(diagnostic_t<=min(r['t1'],r['T']))
                bounds=[r['t0'],min(r['t1'],r['T'])]
                r['hole_fraction_window_max_at03']=float(np.r_[diagnostic_h[inside],np.interp(bounds,diagnostic_t,diagnostic_h)].max())
                r['interior_ratio_window_min']=float(np.r_[diagnostic_i[inside],np.interp(bounds,diagnostic_t,diagnostic_i)].min())
            front_runs.extend(rows);front_records.extend(records)
    bulk_final = [r for r in bulk_runs if r['ensemble_complete']]
    bulk_summary, survival = bulk.summarize(bulk_final,boot=boot,stable_groups=True)
    rng = np.random.default_rng(1729)
    window_groups = defaultdict(list)
    for row in bulk_windows:
        if row['ensemble_complete']:
            key = tuple(row[k] for k in ('p','q','Lx','Ly','T','start','end'))
            window_groups[key].append(row)
    window_summary = []
    for key,rs in sorted(window_groups.items()):
        rng=group_rng(key,1729)
        item = dict(zip(('p','q','Lx','Ly','T','start','end'),key))
        for name in ('rho','rho_drift','chi'):
            item.update({name+'_'+k:v for k,v in mean_summary([r[name] for r in rs],boot,rng).items()})
        window_summary.append(item)
    # A geometric loss cannot silently remove a seed from a front fit.
    front_group = defaultdict(list)
    config = ('p','q','L','Ly','T','ab','threshold','smooth_normal','t0','t1','estimator')
    for r in front_runs:
        if r['ensemble_complete']:
            front_group[tuple(r[k] for k in config)].append(r)
    valid_front, geometry = [], []
    for key,rs in sorted(front_group.items()):
        item = dict(zip(config,key))
        good = [r for r in rs if fronts.covers_requested_window(r)]
        item.update(n_planned=len(rs),n_valid=len(good),all_valid=len(good)==len(rs),
                    valid_runs=[r['run'] for r in good],invalid_runs=[r['run'] for r in rs if r not in good])
        interior=[r['interior_ratio_mean'] for r in rs if np.isfinite(r['interior_ratio_mean'])]
        holes=[r['hole_fraction_max'] for r in rs if np.isfinite(r['hole_fraction_max'])]
        item.update(interior_ratio_mean=float(np.mean(interior)) if interior else np.nan,
                    interior_ratio_change_mean=float(np.mean([r['interior_ratio_change'] for r in rs])),
                    hole_fraction_max=float(max(holes)) if holes else np.nan,
                    hole_fraction_window_max_at03=float(max(r['hole_fraction_window_max_at03'] for r in rs)),
                    interior_ratio_window_min=float(min(r['interior_ratio_window_min'] for r in rs)))
        geometry.append(item)
        if item['all_valid']:
            valid_front.extend(rs)
    front_points, front_crossings = fronts.summarize(valid_front,boot=boot,keep_samples=True,stable_groups=True)
    out.mkdir(parents=True,exist_ok=True)
    write_json(out/'status.json',dict(tasks=status,ensembles=ensembles,
        n_planned=len(tasks),n_complete=len(loaded),n_complete_ensembles=sum(e['ensemble_complete'] for e in ensembles),
        n_ensembles=len(ensembles),settings=dict(boot=boot,
        eligibility='all manifest tasks complete; footer, text, profile and GPU marker agree',
        bootstrap_rng='SHA256 scientific group/window streams; unrelated completion cannot alter an existing cohort interval',
        mean_bootstrap='All256 ordered empirical resamples for n=4 means and paired means; percentile coverage remains approximate. Other n use requested Monte Carlo draws. Joint roots and KM estimators unchanged.',
        front_eligibility='every declared seed has a geometrically valid slope covering the full chosen window',
        unfinished_policy='diagnostic rows retained; no final means or curves from partial seed groups',
        absorption_times='GPU consensus is observed between the last positive-activity output and the first zero-activity output; tables retain both bounds')))
    schedule_audited=[dict(id=ident,**row) for ident,row in status.items() if row.get('header_claim_complete')]
    write_json(out/'scheduled_record_audit.json',dict(utc=datetime.now(timezone.utc).isoformat(),
        n_header_complete=len(schedule_audited),n_schedule_gaps=sum(r['n_missing_scheduled']>0 for r in schedule_audited),
        policy='Every manifest regular and ppd2 log time is required; extra records from an older finer cadence are allowed; no gaps are interpolated',
        failures=[r for r in schedule_audited if r['n_missing_scheduled']>0]))
    write_json(out/'front_geometry.json',geometry)
    write_table(out/'front_diagnostics.csv',[{k:v for k,v in r.items() if k not in ('valid_runs','invalid_runs')} for r in geometry])
    write_json(out/'front_crossings.json',front_crossings)
    write_json(out/'orientation_differences.json',orientation_differences(front_crossings))
    paired_drift,paired_crossings=paired_contrasts(valid_front,boot)
    write_table(out/'front5m_paired_drift.csv',paired_drift)
    write_json(out/'front5m_paired_crossings.json',paired_crossings)
    write_json(out/'provenance.json',provenance(inputs,hashes))
    simple_ensembles = [{k:v for k,v in r.items() if k not in ('task_ids','incomplete_ids')} for r in ensembles]
    write_table(out/'ensemble_eligibility.csv',simple_ensembles)
    for name,rows in [('bulk_runs',bulk_runs),('bulk_windows',bulk_windows),('bulk_summary',bulk_summary),
                      ('bulk_window_summary',window_summary),('bulk_survival',survival),('bulk_curves',bulk_curves),
                      ('front_runs',front_runs),('front_records',front_records),('front_points',front_points)]:
        write_table(out/(name+'.csv'),rows)
    write_table(out/'front_crossings.csv',[{k:v for k,v in r.items() if k not in ('coefficients','bootstrap_samples')} for r in front_crossings])
    print(f"{len(loaded)}/{len(tasks)} runs complete; {sum(e['ensemble_complete'] for e in ensembles)}/{len(ensembles)} fixed ensembles eligible; {len(front_crossings)} crossings")
    return dict(status=status,ensembles=ensembles,bulk_summary=bulk_summary,front_crossings=front_crossings)


def orientation_differences(crossings):
    """Independent-direction bootstrap contrasts, conditional on window/model.

    Pointwise95% and Bonferroni98.333% intervals are supplied for the three
    pairwise direction contrasts. Systematic window/threshold spread remains
    separate. Samples are resampled independently because direction seeds differ.
    """
    config=('p','L','Ly','T','threshold','smooth_normal','t0','t1','estimator')
    groups=defaultdict(dict)
    for fit in crossings:
        if fit['p']==.5 and fit['model']=='bracket' and fit.get('bootstrap_samples'):
            groups[tuple(fit[k] for k in config)][fit['ab']]=fit
    rows=[];rng=np.random.default_rng(5341)
    for key,group in sorted(groups.items()):
        rng=group_rng(key,5341)
        if set(group)!={'1-0','2-1','1-1'}:continue
        if any(f['bootstrap_valid']<.8*f['bootstrap_requested'] for f in group.values()):continue
        draws={ab:rng.choice(f['bootstrap_samples'],5000,replace=True) for ab,f in group.items()}
        base=dict(zip(config,key))
        pairs=[('1-0','2-1'),('1-0','1-1'),('2-1','1-1')]
        for a,b in pairs:
            delta=draws[a]-draws[b]
            lo,hi=np.quantile(delta,[.025,.975]);blo,bhi=np.quantile(delta,[.05/6,1-.05/6])
            rows.append(dict(**base,contrast=a+' minus '+b,
                difference=group[a]['crossing']-group[b]['crossing'],
                ci95_low=float(lo),ci95_high=float(hi),
                simultaneous95_low=float(blo),simultaneous95_high=float(bhi)))
        spread=np.ptp(np.asarray(list(draws.values())),axis=0)
        rows.append(dict(**base,contrast='maximum direction spread',
            difference=max(f['crossing'] for f in group.values())-min(f['crossing'] for f in group.values()),
            bootstrap_upper95=float(np.quantile(spread,.95))))
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=ROOT/'data/v2/gpu/runs')
    ap.add_argument('--manifest',type=Path,action='append')
    ap.add_argument('--out',type=Path,default=ROOT/'campaign/v2/results/gpu')
    ap.add_argument('--boot',type=int,default=1000)
    ap.add_argument('--no-hashes',action='store_true')
    args=ap.parse_args()
    reduce_campaign(args.root,args.manifest or [p for p in MANIFESTS if p.exists()],args.out,args.boot,not args.no_hashes)

if __name__=='__main__':main()
