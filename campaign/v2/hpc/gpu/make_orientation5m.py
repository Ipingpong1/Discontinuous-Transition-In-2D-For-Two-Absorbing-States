#!/usr/bin/env python3
"""Predeclare independent 5m front-stationarity cohorts; never submit jobs."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

HERE=Path(__file__).resolve().parent
NAMESPACE='dn-v2-20261103-orientation5m:'
def seed(name):
    return int.from_bytes(hashlib.blake2b((NAMESPACE+name).encode(),digest_size=8).digest(),'little')

base=json.loads((HERE/'production_manifest_profiles50k.json').read_text())
tasks=[]
# Whole four-seed parameter groups stay adjacent. The fixed order improves
# complete-group turnaround without depending on any outcomes of these runs.
for ab in [(1,0),(2,1),(1,1)]:
    for q in ['.04594','.04593','.04595']:
        for replica in range(1,5):
            ident=f'orientation5m_L8192x8192_ab{ab[0]}-{ab[1]}_q{q}_s{replica}'
            tasks.append(dict(id=ident,family='orientation',observable='orientation5m',
                Lx=8192,Ly=8192,q=q,p='0.5',alpha='0.5',seed=seed(ident+':dynamics'),
                init_seed=seed(ident+':initial'),replica=replica,T=5_000_000,
                every=50_000,profile=True,orientation=list(ab),mode='band',
                checkpoint_every=1_000_000,work=8192*8192*5_000_000))
old=[]
for name in ['production_manifest_profiles50k.json','production_manifest_largeL.json','production_manifest_near_large.json','wetting_manifest.json']:
    old.extend(json.loads((HERE/name).read_text())['tasks'])
seeds=[t[k] for t in old+tasks for k in ['seed','init_seed']]
assert len(seeds)==len(set(seeds)),'Seed collision with another declared campaign'
assert len(tasks)==36
manifest=dict(schema=1,campaign='dn-v2-20261103-orientation5m',
    declared_utc=datetime.now(timezone.utc).isoformat(),deadline_utc='2026-11-03T11:00:00+00:00',
    source_sha256=base['source_sha256'],seed_algorithm='BLAKE2b-64; namespace '+NAMESPACE[:-1],
    families={'orientation':dict(workers=4,tasks=36)},
    purpose='Independent fixed-seed test of front stationarity beyond the original2m cohorts; no trajectories reused or selected by outcome.',
    analysis_windows=[[2_000_000,5_000_000],[2_000_000,3_500_000],[3_500_000,5_000_000]],
    group_order='normal(1,0),(2,1),(1,1); q=.04594,.04593,.04595; replicas1..4; every declared group retained',
    tasks=tasks)
manifest['design_amendment']=json.loads((HERE/'orientation5m_amendment_metadata.json').read_text())
out=HERE/'production_manifest_orientation5m.json'
out.write_text(json.dumps(manifest,indent=2)+'\n')
print(out);print(hashlib.sha256(out.read_bytes()).hexdigest())
print('36 tasks; 9 fixed groups; all72 initial/dynamical seeds distinct from all earlier declared seeds')
