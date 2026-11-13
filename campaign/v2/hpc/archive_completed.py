#!/usr/bin/env python3
"""Verify an already-copied legacy directory, then preserve its path as a symlink.

Used only for the closed R2_orient campaign after the account's scratch quota
filled. No data content is deleted: the byte-verified archive becomes canonical.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import shutil
import sys

name = sys.argv[1] if len(sys.argv)>1 else 'R2_orient'
allowed = {'R2_orient','R2_orient2k','R2_wall','E12_phase','E9_alpha','E4_voter','E10_lazy_front'}
plan = Path('/archive/mo2273/dn-v2-storage-recovery-20261103/remaining_closed_campaigns_plan.json')
if plan.exists():
    allowed.update(x['name'] for x in json.loads(plan.read_text())['campaigns'])
assert name in allowed and '/' not in name and '..' not in name
source = Path('/scratch/mo2273/dn/runs') / name
destination = Path('/archive/mo2273/dn-v2-storage-recovery-20261103') / name
record = destination.parent / (name+'_verified_transfer.json')
assert not source.is_symlink(), 'Already archived; refusing to repeat'
assert source.is_dir() and destination.is_dir()

def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()

files = []
for path in sorted(source.rglob('*')):
    relative = path.relative_to(source)
    target = destination / relative
    assert not path.is_symlink(), f'Unexpected symlink: {path}'
    if path.is_dir():
        assert target.is_dir(), str(target)
        continue
    assert path.is_file() and target.is_file(), str(target)
    assert path.stat().st_size == target.stat().st_size, str(target)
    sha = digest(path)
    assert sha == digest(target), f'Checksum mismatch: {target}'
    files.append(dict(path=str(relative), bytes=path.stat().st_size, sha256=sha))
assert len(files) > 0, 'Unexpectedly empty source'
record.write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),
    source=str(source), destination=str(destination), verified_files=files,
    bytes=sum(f['bytes'] for f in files),
    strategy='Byte-verified archive; original directory replaced by symlink'), indent=2)+'\n')
shutil.rmtree(source)
source.symlink_to(destination, target_is_directory=True)
assert source.resolve() == destination.resolve()
print(json.dumps(dict(files=len(files), bytes=sum(f['bytes'] for f in files),
                     original_path=str(source), archive=str(destination))))
