"""Fixed-manifest CPU lifetime cohorts; partial groups are diagnostic only."""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import argparse,csv,json,math
from .io import ROOT,read_header,provenance,write_csv,write_json
from .bulk import collect,summarize


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data',type=Path,default=ROOT/'data/v2/cpu')
    ap.add_argument('--manifest',type=Path,default=ROOT/'campaign/v2/hpc/cpu/manifest.json')
    ap.add_argument('--out',type=Path,default=ROOT/'campaign/v2/results/cpu_lifetimes')
    ap.add_argument('--boot',type=int,default=2000);args=ap.parse_args()
    manifest=json.loads(args.manifest.read_text());planned=defaultdict(set);task_by_seed={}
    for task in manifest['tasks']['lifetimes']:
        a=task['args'];value=lambda x:a[a.index(x)+1]
        p=float(value('--voter'));q=float(value('--jflip').split(':')[1]);L=int(value('--L'));T=int(value('--T'))
        key=(p,q,L,L,T);planned[key].add(task['seed']);task_by_seed[task['seed']]=task
    paths=sorted(set(args.data.glob('*.txt'))|set(args.data.glob('lifetimes/*.txt'))
                 |set(args.data.glob('nearline/*.txt')));valid=[];rejected=[]
    for path in paths:
        try:
            h=read_header(path)
            if h.get('seed') not in task_by_seed:raise ValueError('unplannedseed')
            key=(h['p'],h['q'],h['Lx'],h['Ly'],h['T'])
            if h['seed'] not in planned[key]:raise ValueError('metadata mismatch')
            if h.get('done'):
                task=task_by_seed[h['seed']];a=task['args']
                every=int(a[a.index('--every')+1])
                ppd=int(a[a.index('--ppd')+1])
                recorded={int(line.split()[0]) for line in path.read_text().splitlines()
                          if line and not line.startswith('#')}
                expected=set(range(0,h['t_end']+1,every))|{h['t_end']}
                expected|={t for d in range(11) for i in range(ppd)
                           if (t:=math.floor(10.**(d+i/ppd)+.5))<=h['t_end']}
                missing=expected-recorded
                if missing:raise ValueError(f'incomplete scheduled observation grid: {len(missing)} missing times')
            valid.append(path)
        except (ValueError,KeyError) as error:rejected.append(dict(path=str(path),error=str(error)))
    runs=collect(valid);summary,survival=summarize(runs,boot=args.boot,stable_groups=True)
    seen=defaultdict(set)
    for r in runs:
        if r['done']:seen[(r['p'],r['q'],r['Lx'],r['Ly'],r['T'])].add(r['seed'])
    complete_keys=set()
    for r in summary:
        key=tuple(r[k] for k in ('p','q','Lx','Ly','T'))
        r['n_planned']=len(planned[key]);r['fixed_ensemble_complete']=seen[key]==planned[key]
        if r['fixed_ensemble_complete']:complete_keys.add(key)
    args.out.mkdir(parents=True,exist_ok=True)
    write_csv(args.out/'runs.csv',runs);write_csv(args.out/'summary_progress.csv',summary)
    with (args.out/'summary_complete.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(r for r in summary if r['fixed_ensemble_complete'])
    with (args.out/'survival_complete.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(survival[0]));w.writeheader()
        w.writerows(r for r in survival if tuple(r[k] for k in ('p','q','Lx','Ly','T')) in complete_keys)
    repair_audit=ROOT/'reports/v2/hpc/cpu_repairs_applied.json'
    write_json(args.out/'analysis.json',dict(settings=dict(boot=args.boot,complete_groups_only=True,
        half_threshold=.145,time_origin=1000,sampling_unit='independent planned run',
        bootstrap_streams='Stable independent stream for each parameter/size/horizon group',
        required_observations='All manifest regular and logarithmic scheduled times, plus terminal observation'),
        complete_keys=sorted(complete_keys),rejected_inputs=rejected,
        repair_provenance=provenance([repair_audit]) if repair_audit.exists() else [],
        manifest=provenance([args.manifest]),inputs=provenance(valid),summary=summary))
    print(len(runs),'started runs;',sum(r['done'] for r in runs),'finished;',len(complete_keys),'complete fixed ensembles')
    for r in summary:
        print(r['q'],r['Lx'],r['n'],'/',r['n_planned'],r['half_events'],r['absorption_events'])
if __name__=='__main__':main()
