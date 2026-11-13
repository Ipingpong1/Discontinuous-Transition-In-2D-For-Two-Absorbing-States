"""Verify an exact-seed GPU repair against preserved quota-damaged evidence."""
import hashlib
import json
from pathlib import Path
import numpy as np
from .io import read_timeseries,read_profile,provenance


def sha256(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):value.update(block)
    return value.hexdigest()


def verify_replay(stage,evidence,task,validation):
    stage,evidence=Path(stage),Path(evidence)
    request=json.loads((evidence/'repair_request.json').read_text())
    for record in request.get('evidence',[]):
        if sha256(record['path'])!=record['sha256']:
            raise ValueError('Preserved repair evidence digest changed: '+record['path'])
    expected=request['expected_final_sha256']
    if validation['final_sha256']!=expected:
        raise ValueError('Repair final state differs from the original completed-state digest')
    _,old=read_timeseries(evidence/'out.txt');_,new=read_timeseries(stage/'out.txt')
    index={int(row[0]):i for i,row in enumerate(new)}
    if any(int(row[0]) not in index for row in old):
        raise ValueError('Repair omits a preserved text observation')
    matched=new[[index[int(t)] for t in old[:,0]]]
    if not np.array_equal(old,matched):
        bad=np.flatnonzero(np.any(old!=matched,axis=1))
        raise ValueError(f'Repair disagrees with {len(bad)} preserved text records; first t={old[bad[0],0]:g}')
    report=dict(task=task['id'],verified=True,final_sha256=expected,final_state_digest_matches=True,
        overlapping_text_records=len(old),text_columns_exact=True,
        new_text_records=len(new),filled_text_times=[int(t) for t in new[:,0] if t not in set(old[:,0])])
    files=[evidence/'out.txt',evidence/'done.json',stage/'out.txt',stage/'done.json']
    if task['profile']:
        bins=task['Lx'] if task['mode']=='band' else task['Ly']
        old_profile,old_trailing=read_profile(evidence/'profile.bin',bins)
        new_profile,new_trailing=read_profile(stage/'profile.bin',bins)
        if new_trailing:raise ValueError('Repair has a partial profile record')
        positions={int(t):i for i,t in enumerate(new_profile['t'])}
        if any(int(t) not in positions for t in old_profile['t']):
            raise ValueError('Repair omits a preserved profile observation')
        matched_profile=new_profile[[positions[int(t)] for t in old_profile['t']]]
        if not np.array_equal(old_profile,matched_profile):
            raise ValueError('Repair disagrees with preserved profile counts')
        report.update(overlapping_profile_records=len(old_profile),profile_counts_exact=True,
            old_partial_profile_bytes=old_trailing,new_profile_records=len(new_profile),
            filled_profile_times=[int(t) for t in new_profile["t"] if t not in set(old_profile["t"])])
        files += [evidence/'profile.bin',stage/'profile.bin']
    report['inputs']=provenance(files)
    return report
