#!/usr/bin/env python3
"""Declare independent bounded axial10M confirmation; this script never submits."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
HERE=Path(__file__).resolve().parent
NAMESPACE='dn-v2-20261103-axial10m:'
def seed(name):
    return int.from_bytes(hashlib.blake2b((NAMESPACE+name).encode(),digest_size=8).digest(),'little')
base_path=HERE/'production_manifest_orientation5m.json'
base=json.loads(base_path.read_text());template=base['tasks'][0]
tasks=[]
for q in ['.04594','.04593','.04595']:
    for replica in range(1,5):
        ident=f'axial10m_L8192x8192_ab1-0_q{q}_s{replica}'
        tasks.append(dict(template,id=ident,observable='axial10m',q=q,orientation=[1,0],
            seed=seed(ident+':dynamics'),init_seed=seed(ident+':initial'),replica=replica,
            T=10_000_000,work=8192*8192*10_000_000))
prior=[]
for path in HERE.glob('*manifest*.json'):
    if path.name=='production_manifest_axial10m.json':continue
    prior.extend(task[key] for task in json.loads(path.read_text()).get('tasks',[]) for key in ('seed','init_seed') if key in task)
new=[task[key] for task in tasks for key in ('seed','init_seed')]
assert len(set(new))==24 and not set(new)&set(prior)
assert len({v&0xffffffff for v in new})==24 and not {v&0xffffffff for v in new}&{v&0xffffffff for v in prior}
manifest=dict(schema=1,campaign='dn-v2-20261103-axial10m',
    declared_utc=datetime.now(timezone.utc).isoformat(),deadline_utc='2026-11-03T11:05:00+00:00',
    source_sha256=base['source_sha256'],seed_algorithm='BLAKE2b-64; namespace '+NAMESPACE[:-1],
    source_manifest=base_path.name,source_manifest_sha256=hashlib.sha256(base_path.read_bytes()).hexdigest(),
    families={'orientation':dict(workers=12,tasks=12)},
    motivation='Independent extension of the axial front time horizon after the5M study resolved temporal drift; no outcomes from these new trajectories existed at declaration.',
    eligibility='All four new independent runs perq must complete and cover each fixed geometry window. No pooling withT2M/T5M cohorts; no survivor selection; unbracketed roots withheld.',
    analysis_windows=[[2_000_000,5_000_000],[5_000_000,10_000_000]],
    primary_paired_contrast='Same-run5–10M minus2–5M velocity; same within-q resample indices for paired root differences; allthreeq required.',
    submission_policy='Ordinary queue behind all existing work; unchanged validated backend and >80% utilization monitor; new campaign11:05UTC ceiling only.',tasks=tasks)
out=HERE/'production_manifest_axial10m.json'
out.write_text(json.dumps(manifest,indent=2)+'\n')
print(out);print(hashlib.sha256(out.read_bytes()).hexdigest())
