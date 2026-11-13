#!/usr/bin/env python3
"""Audit v2 grid files and release only claims of confirmed failed workers.

Run on Torch after storage space is restored. Damaged files and metadata are
renamed for the audit trail, never silently accepted or erased.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import struct
import subprocess

root = Path('/scratch/mo2273/dn/v2/gpu')
code = Path('/scratch/mo2273/dn/v2/campaign/v2/hpc/gpu')
failed = {21718939,21718940,21718941,21718942,21719508,21719509,
          21720047,21720048,21720050,21720051,21720052,21720053}
accounting = subprocess.check_output(['sacct','-X','-j',','.join(map(str,sorted(failed))),
    '--format=JobIDRaw,State','-n','-P'], text=True)
states = dict(line.split('|')[:2] for line in accounting.splitlines() if line.strip())
assert all(states.get(str(j),'').split()[0] in {'FAILED','CANCELLED','COMPLETED','TIMEOUT'} for j in failed), states
tasks = {}
for name in ['production_manifest_profiles50k.json','production_manifest_largeL.json']:
    tasks.update({t['id']:t for t in json.loads((code/name).read_text())['tasks']})
audit = dict(utc=datetime.now(timezone.utc).isoformat(), states=states, files=[], released=[])

def quarantine(path):
    target = path.with_name(path.name+'.quota_invalid_20261103')
    assert not target.exists(), str(target)
    path.rename(target)
    return str(target)

for family in ['orientation','flagship_bulk','near_voter']:
    for directory in sorted((root/'runs'/family).glob('*')):
        if directory.name not in tasks:
            continue
        task = tasks[directory.name]
        validity = {}
        for name in ['initial.bin','checkpoint.bin','final.bin']:
            path = directory/name
            if not path.exists():
                continue
            with path.open('rb') as stream:
                header = stream.read(32)
            good = len(header)==32 and header[:8]==b'DNSIM1\0\0'
            dims = struct.unpack('<QQQ',header[8:]) if good else (None,None,None)
            good = good and dims[:2]==(task['Lx'],task['Ly']) and dims[2]<=task['T']
            good = good and path.stat().st_size==32+task['Lx']*task['Ly']//8
            if name=='initial.bin':
                good = good and dims[2]==0
            item = dict(path=str(path), bytes=path.stat().st_size, header=list(dims), valid=bool(good))
            if not good:
                item['quarantined_as'] = quarantine(path)
            audit['files'].append(item)
            validity[name] = (bool(good),dims[2])
        done = directory/'done.json'
        completed = False
        if done.exists():
            try:
                marker = json.loads(done.read_text())
                valid_final, final_t = validity.get('final.bin',(False,None))
                completed = (valid_final and marker.get('id')==task['id']
                             and marker.get('t_end')==final_t
                             and (final_t==task['T'] or marker.get('absorbed') is True))
                if completed:
                    completed = hashlib.sha256((directory/'final.bin').read_bytes()).hexdigest()==marker['final_sha256']
            except (ValueError,KeyError,OSError):
                completed = False
            if not completed:
                audit['files'].append(dict(path=str(done), valid=False,quarantined_as=quarantine(done)))
        if completed:
            continue
        claim = root/'claims'/family/task['id']
        owner = claim/'owner.json'
        if not owner.exists():
            continue
        info = json.loads(owner.read_text())
        if int(info['job']) not in failed:
            continue
        audit['released'].append(dict(claim=str(claim),owner=info))
        owner.unlink()
        claim.rmdir()

out = root/'quota_recovery_20261103.json'
out.write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(dict(audited_files=len(audit['files']), released_claims=len(audit['released']),
                     invalid_files=sum(not x['valid'] for x in audit['files']), report=str(out))))
