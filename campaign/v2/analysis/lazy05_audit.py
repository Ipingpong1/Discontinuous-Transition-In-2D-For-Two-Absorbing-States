"""Corrected E18 lazy-front envelope audit; no additional simulations.

All declared three-seed groups must finish and cover each fixed fit window.
Time in input/reduced rows is in update steps; physical time is r*t, r=.25.
"""
from collections import defaultdict
from pathlib import Path
import json
import csv
import shlex
import re
import numpy as np
from . import fronts
from .io import ROOT, provenance, write_csv, write_json, profile_path

QS=(.038556,.039474,.039933,.040392,.040851,.04131,.042228,.044064)
WINDOWS=((4000,40000),(40000,100000),(100000,200000),(4000,200000))
R=.25


def main():
    folder=ROOT/'data/runs/E18_lazy05_front/out'
    out=ROOT/'campaign/v2/results/lazy05_fronts'
    manifest=folder.parent/'manifest.tsv'
    declared={Path(row['prefix']).name:shlex.split(row['args']) for row in csv.DictReader(manifest.open(),delimiter='\t')}
    if len(declared)!=24:raise ValueError('Legacy manifest has unexpected cohort count')
    reference=ROOT/'campaign/v2/results/coexistence_reference.json'
    all_rows=[];records=[];inputs=[manifest,reference,Path(fronts.__file__),Path(__file__).with_name('io.py'),Path(__file__).with_name('stats.py')]
    for q in QS:
        for seed in (1,2,3):
            name=f'p0.5_a0.5_r0.25_j2_q{q:g}_Lx2048_Ly4096_s{seed}.txt'
            path=folder/name
            if not path.exists() or not profile_path(path).exists():raise ValueError('Missing declared run '+name)
            rows,rec=fronts.reduce_run(path,thresholds=(.2,.3,.5),tmins=(),windows=WINDOWS)
            args=declared[path.stem]
            if rows[0]['seed']!=int(args[args.index('--seed')+1]):raise ValueError('RNG seed differs from manifest '+name)
            if not rows[0]['done']:raise ValueError('Incomplete declared run '+name+': '+rows[0]['completion_errors'])
            table_line=next(line for line in path.read_text().splitlines() if line.startswith('# ftab '))
            actual_table=np.array([[float(x) for x in part.split(',')] for part in table_line[7:].split(';')])
            declared_table=np.array([[float(x) for x in part.split(',')] for part in args[args.index('--ftab')+1].split(';')])
            if not np.array_equal(actual_table,declared_table):raise ValueError('Actual flip table differs from manifest '+name)
            inferred_q=float(actual_table[0,2]/R-.125)
            if abs(inferred_q-q)>1e-10:raise ValueError('Filename/table probability mismatch '+name)
            for row in rows:
                row.update(update_probability=R,physical_t0=R*row['t0'],physical_t1=R*row['t1'],
                    physical_velocity=row['v']/R,manifest_seed_label=seed)
                row['full_window']=fronts.covers_requested_window(row)
            all_rows.extend(rows);records.extend(rec);inputs.extend([path,profile_path(path)])
    groups=defaultdict(list)
    for row in all_rows:groups[tuple(row[k] for k in ('q','threshold','t0','t1','estimator'))].append(row)
    valid=[];eligibility=[]
    for key,rows in sorted(groups.items()):
        good=len(rows)==3 and all(row['full_window'] for row in rows)
        item=dict(zip(('q','threshold','t0','t1','estimator'),key))
        item.update(n_planned=3,n_valid=sum(r['full_window'] for r in rows),eligible=good)
        eligibility.append(item)
        if good:valid.extend(rows)
    points,fits=fronts.summarize(valid,boot=2000,qmin_relative=0,stable_groups=True)
    for row in fits:row.update(physical_t0=R*row['t0'],physical_t1=R*row['t1'])
    ref=json.loads(reference.read_text())
    for row in fits:row['ratio_to_shared_reference']=row['crossing']/ref['qstar']
    write_csv(out/'runs.csv',all_rows);write_csv(out/'records.csv',records)
    write_csv(out/'points.csv',points);write_csv(out/'eligibility.csv',eligibility)
    write_json(out/'analysis.json',dict(update_probability=R,time_conversion='physical time = .25 * update step; physical velocity = per-step velocity/.25',
        n_runs=24,n_seeds_per_q=3,planned_q=QS,reference=ref,
        policy='All three declared seeds completed; all three must retain valid geometry throughout each selected window. No truncation-selected means.',
        fits=fits,provenance=provenance(inputs+[Path(__file__)])))
    for row in fits:
        if row['threshold']==.3 and row['estimator']=='envelope' and row['model']=='bracket':
            print({k:row[k] for k in ('physical_t0','physical_t1','crossing','ci_low','ci_high','ratio_to_shared_reference','bootstrap_valid')})


if __name__=='__main__':main()
