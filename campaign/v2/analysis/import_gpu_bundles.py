"""Import only v3 COMPLETED bundles after read-only rsync from gpu/bundles/.

Checkpoint archives belong under gpu/checkpoints/ and must never be passed here.
The flat tar layout is validated, extracted to local staging, checked against
its manifest and final-state digest, then replaces the previous task directory.
Old interrupted local files are retained under previous_runs/ for the audit.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import struct
import tarfile
import tempfile
from .gpu_campaign import MANIFESTS, task_completion
from .gpu_repair import verify_replay
from .io import ROOT, manifest_tasks, write_json


def digest(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            value.update(block)
    return value.hexdigest()


def safe_unpack(bundle,directory):
    seen=set();total=0
    with tarfile.open(bundle,'r:gz') as archive:
        for member in archive:
            name=PurePosixPath(member.name)
            if (not member.isfile() or name.is_absolute() or len(name.parts)!=1
                    or member.name != name.name or name.name in ('.','..')
                    or '\\' in member.name or member.name in seen):
                raise ValueError(f'Unsafe or duplicate bundle member: {member.name}')
            seen.add(member.name);total+=member.size
            if member.size>1_000_000_000 or total>2_000_000_000:
                raise ValueError('Unexpectedly large task bundle')
            source=archive.extractfile(member)
            if source is None:raise ValueError(f'Missing member contents: {member.name}')
            with source,(Path(directory)/member.name).open('wb') as target:
                shutil.copyfileobj(source,target)
    return seen


def validate_final(directory,task):
    directory=Path(directory)
    status,meta,series=task_completion(task,directory/'out.txt')
    # task_completion normally gets the id from its parent directory. The
    # staging directory is deliberately named exactly as the declared id.
    if not status['complete']:
        raise ValueError('Incomplete task bundle: '+', '.join(status['reasons']))
    final=directory/'final.bin'
    if not final.exists() or final.stat().st_size!=32+task['Lx']*task['Ly']//8:
        raise ValueError('Invalid final-state size')
    with final.open('rb') as stream:header=stream.read(32)
    if header[:8]!=b'DNSIM1\0\0':raise ValueError('Invalid final-state magic')
    x,y,t=struct.unpack('<QQQ',header[8:])
    if (x,y,t)!=(task['Lx'],task['Ly'],meta['t_end']):
        raise ValueError('Final-state dimensions/time mismatch')
    marker=json.loads((directory/'done.json').read_text())
    sha=digest(final)
    if marker.get('final_sha256')!=sha:raise ValueError('Final-state digest mismatch')
    return dict(t_end=t,absorbed=meta['absorbed'],final_sha256=sha)


def import_one(bundle,task,root,records):
    bundle,root=Path(bundle),Path(root)
    if 'checkpoints' in bundle.parts:
        raise ValueError('Checkpoint bundles are not publication results')
    sha=digest(bundle);ident=task['id']
    destination=root/'runs'/task['family']/ident
    if records.get(ident,{}).get('bundle_sha256')==sha and destination.exists():
        return 'unchanged'
    staging=root/'import_staging';staging.mkdir(parents=True,exist_ok=True)
    parent=Path(tempfile.mkdtemp(prefix='task-',dir=staging));stage=parent/ident;stage.mkdir()
    previous=None
    try:
        safe_unpack(bundle,stage)
        validation=validate_final(stage,task)
        evidence=root/'repair_evidence'/ident
        repair_report=None
        if (evidence/'repair_request.json').exists():
            repair_report=verify_replay(stage,evidence,task,validation)
            for record in repair_report['inputs']:
                try:relative=Path(record['path']).relative_to(stage.resolve())
                except ValueError:continue
                record['path']=str((destination/relative).resolve())
            write_json(stage/'repair_verification.json',repair_report)
        if destination.is_symlink():raise ValueError('Refusing to replace a symlink task directory')
        if (destination/'done.json').exists():
            try:old=json.loads((destination/'done.json').read_text())
            except (OSError,ValueError):old={}  # Preserve quota-damaged marker in previous_runs.
            if old.get('t_end')==validation['t_end'] and old.get('final_sha256') not in (None,validation['final_sha256']):
                raise ValueError('Conflict with an already completed task final-state digest')
        note=dict(task=ident,family=task['family'],bundle=str(bundle.resolve()),bundle_sha256=sha,
                  imported_utc=datetime.now(timezone.utc).isoformat(),repair_verified=repair_report is not None,**validation)
        write_json(stage/'import_record.json',note)
        destination.parent.mkdir(parents=True,exist_ok=True)
        if destination.exists():
            previous=root/'previous_runs'/task['family']/(ident+'-'+sha[:12])
            previous.parent.mkdir(parents=True,exist_ok=True)
            if previous.exists():raise ValueError('Previous-run archive already exists; inspect before retry')
            destination.rename(previous)
        try:stage.rename(destination)
        except Exception:
            if previous is not None and not destination.exists():previous.rename(destination)
            raise
        records[ident]=note
        if repair_report is not None:
            write_json(root/'repair_verifications'/(ident+'.json'),repair_report)
        return 'imported'
    finally:
        shutil.rmtree(parent,ignore_errors=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bundles',type=Path,default=ROOT/'data/v2/gpu/bundles')
    ap.add_argument('--root',type=Path,default=ROOT/'data/v2/gpu')
    ap.add_argument('--manifest',type=Path,action='append')
    args=ap.parse_args()
    if 'checkpoints' in args.bundles.parts:raise SystemExit('Refusing checkpoint source')
    manifests=args.manifest or [p for p in MANIFESTS+[ROOT/'campaign/v2/hpc/gpu/wetting_manifest.json'] if p.exists()]
    tasks=manifest_tasks(manifests)
    ledger=args.root/'bundle_imports.json'
    records=json.loads(ledger.read_text()) if ledger.exists() else {}
    counts={'imported':0,'unchanged':0,'rejected':0};errors=[]
    for bundle in sorted(args.bundles.glob('*/*.tar.gz')):
        ident=bundle.name[:-7]
        try:
            if ident not in tasks:raise ValueError('Unknown manifest task')
            task=tasks[ident]
            if bundle.parent.name!=task['family']:raise ValueError('Bundle family mismatch')
            state=import_one(bundle,task,args.root,records);counts[state]+=1
        except (OSError,ValueError,KeyError,tarfile.TarError) as error:
            counts['rejected']+=1;errors.append(dict(bundle=str(bundle),error=str(error)))
    write_json(ledger,records)
    write_json(args.root/'bundle_import_status.json',dict(utc=datetime.now(timezone.utc).isoformat(),counts=counts,errors=errors))
    print(json.dumps(counts))
    for error in errors:print(error)
    if errors:raise SystemExit(2)

if __name__=='__main__':main()
