#!/usr/bin/env python3
"""Predeclare the longer, larger near-voter bulk replication.

Motivation recorded before submission: the first L=4096, q=.000640 run
decreased from the earlier apparent plateau by t=10^7. Test finite size and
late-time drift with four new independent seeds per q, not selected continuations.
"""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
source = HERE / '../../../round2/dnsim.c'

def seed(name):
    return int.from_bytes(hashlib.blake2b(('dn-v2-20261103:' + name).encode(),
                                        digest_size=8).digest(), 'little')

tasks = []
for q in ['.000640', '.000644', '.000648']:
    for replica in range(1, 5):
        name = f'near_bulk20m_L8192x8192_q{q}_s{replica}'
        tasks.append(dict(id=name, family='near_voter', observable='near_bulk',
            Lx=8192, Ly=8192, q=q, p='0.995', alpha='0.5',
            seed=seed(name+':dynamics'), init_seed=seed(name+':initial'),
            replica=replica, T=20_000_000, every=50_000, profile=False,
            orientation=None, mode='bulk', checkpoint_every=1_000_000,
            work=8192**2*20_000_000))
tasks.sort(key=lambda t: seed(t['id']+':order'))
all_seeds = [t[k] for t in tasks for k in ('seed', 'init_seed')]
assert len(all_seeds) == len(set(all_seeds))
for path in HERE.glob('production_manifest*.json'):
    if path.name == 'production_manifest_near_large.json':
        continue
    previous = json.loads(path.read_text())
    used = {t[k] for t in previous['tasks'] for k in ('seed', 'init_seed')}
    assert not used.intersection(all_seeds), f'seed collision with {path}'
manifest = dict(schema=1, campaign='dn-v2-20261103-near-large',
    deadline_utc='2026-11-03T11:00:00+00:00',
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    seed_algorithm='BLAKE2b-64; namespace dn-v2-20261103',
    motivation='Independent longer/larger replication of late density drift at p=.995; fixed four-seed ensembles.',
    families={'near_voter': {'workers': 3, 'tasks': len(tasks)}}, tasks=tasks)
output = HERE / 'production_manifest_near_large.json'
output.write_text(json.dumps(manifest, indent=2)+'\n')
print(output)
print(len(tasks), 'runs;', sum(t['work'] for t in tasks), 'site updates')
