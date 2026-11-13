#!/usr/bin/env python3
"""Declare one independent midpoint cohort before any of its outcomes exist."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json

HERE=Path(__file__).resolve().parent
NAMESPACE='dn-v2-20261103-near-midpoint:'
def seed(name):
    return int.from_bytes(hashlib.blake2b((NAMESPACE+name).encode(),digest_size=8).digest(),'little')
base_path=HERE/'production_manifest_near_large.json'
base=json.loads(base_path.read_text());template=base['tasks'][0]
tasks=[]
for replica in range(1,5):
    ident=f'near_midpoint20m_L8192x8192_q.000642_s{replica}'
    task=dict(template,id=ident,q='.000642',seed=seed(ident+':dynamics'),
        init_seed=seed(ident+':initial'),replica=replica)
    tasks.append(task)
prior=[]
for path in HERE.glob('*manifest*.json'):
    if path.name=='production_manifest_near_midpoint.json':continue
    data=json.loads(path.read_text())
    prior.extend(value for task in data.get('tasks',[]) for key in ('seed','init_seed') if (value:=task.get(key)) is not None)
new=[task[key] for task in tasks for key in ('seed','init_seed')]
assert len(set(new))==8 and not set(new)&set(prior)
assert len({v&0xffffffff for v in new})==8 and not {v&0xffffffff for v in new}&{v&0xffffffff for v in prior}
manifest=dict(schema=1,campaign='dn-v2-20261103-near-midpoint',
    declared_utc=datetime.now(timezone.utc).isoformat(),deadline_utc='2026-11-03T11:00:00+00:00',
    source_sha256=base['source_sha256'],seed_algorithm='BLAKE2b-64; namespace '+NAMESPACE[:-1],
    source_manifest=base_path.name,source_manifest_sha256=hashlib.sha256(base_path.read_bytes()).hexdigest(),
    families={'near_voter':dict(workers=4,tasks=4)},
    motivation='Independent midpoint q=.000642 between existing .000640 late-collapse trajectories and .000644 persistent trajectories; a finite-size/time crossover refinement, not an asymptotic phase-boundary determination.',
    eligibility='All four predeclared independent runs must complete. Retain absorptions; no pooling with selected survivors from prior q values.',
    analysis_windows=[[5_000_000,10_000_000],[10_000_000,20_000_000]],
    submission_policy='Four new workers behind existing submitted work, no priorities/cancellations changed; hard 11 UTC deadline.',
    tasks=tasks)
path=HERE/'production_manifest_near_midpoint.json'
path.write_text(json.dumps(manifest,indent=2)+'\n')
print(path);print(hashlib.sha256(path.read_bytes()).hexdigest())
