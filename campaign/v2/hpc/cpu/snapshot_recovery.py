#!/usr/bin/env python3
"""Preserve original16M checkpoints after its jobs stop, outside publication paths."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import io
import json
import struct
import subprocess
import tarfile
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
REMOTE = r'''
from datetime import datetime,timezone
from pathlib import Path
import hashlib,io,json,os,struct,subprocess,sys,tarfile
root=Path('/home/mo2273/dn_v2_stage/cpu')
jobs={'21723013','21725907','21725908'}
states=subprocess.check_output(['squeue','-h','-u','mo2273','-o','%i|%T'],text=True).splitlines()
active=[line for line in states if line.split('|')[0] in jobs]
if active:raise RuntimeError('CPU recovery snapshot waits for stopped jobs: '+str(active))
paths=[root/'manifest.json']+sorted((root/'checkpoints').glob('tail_*.ckpt'))
paths+=sorted((root/'out/tails').glob('*.txt'))+sorted((root/'logs').glob('*.jsonl'))
blobs={str(path.relative_to(root)):path.read_bytes() for path in paths}
assert len([name for name in blobs if name.startswith('checkpoints/')])==8
env=dict(os.environ,TZ='UTC')
blobs['slurm_final_utc.txt']=subprocess.check_output(
 ['sacct','-j',','.join(sorted(jobs)),'--format=JobIDRaw,JobName,State,ExitCode,Elapsed,End','-X','-P','-n'],env=env)
manifest=dict(snapshot_utc=datetime.now(timezone.utc).isoformat(),source_root=str(root),
 jobs=sorted(jobs),active_jobs=active,purpose='Recovery provenance only; excluded original16M ensemble, no selected-prefix inference.',
 files={name:dict(bytes=len(blob),sha256=hashlib.sha256(blob).hexdigest()) for name,blob in blobs.items()})
blobs['snapshot_manifest.json']=(json.dumps(manifest,indent=2)+'\n').encode()
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as archive:
 for name,blob in sorted(blobs.items()):
  info=tarfile.TarInfo(name);info.size=len(blob);info.mtime=0;info.mode=0o600
  archive.addfile(info,io.BytesIO(blob))
'''


def main():
    now = datetime.now(timezone.utc)
    assert now >= datetime(2026, 11, 3, 11, 25, 15, tzinfo=timezone.utc)
    directory = ROOT/'data/v2/cpu/provenance/recovery'
    directory.mkdir(parents=True, exist_ok=True)
    path = directory/f'original16m_{now:%Y%m%dT%H%M%SZ}.tar.gz'
    with path.open('wb') as stream:
        result = subprocess.run(['ssh', 'torch', 'python3', '-'], input=REMOTE.encode(),
                                stdout=stream, stderr=subprocess.PIPE, timeout=60)
    if result.returncode:
        path.unlink(missing_ok=True)
        raise RuntimeError(result.stderr.decode())
    rows = []
    with tarfile.open(path, 'r:gz') as archive:
        manifest = json.load(archive.extractfile('snapshot_manifest.json'))
        assert set(archive.getnames()) == set(manifest['files']) | {'snapshot_manifest.json'}
        for name, expected in manifest['files'].items():
            blob = archive.extractfile(name).read()
            assert len(blob) == expected['bytes']
            assert hashlib.sha256(blob).hexdigest() == expected['sha256']
            if not name.startswith('checkpoints/'):
                continue
            assert blob[:8] == b'DNREPS1\x00'
            header = struct.unpack('<12Q2d', blob[8:120])
            assert header[6:8] == (2_000_000, 0)
            nt, groups = header[3:5]
            counts = np.frombuffer(blob[120:120+8*groups], dtype='<u8')
            accumulators = np.frombuffer(blob[120+8*groups:], dtype='<f8').reshape(groups*nt, 8)
            raw_name = 'out/tails/'+Path(name).name.replace('.ckpt', '.txt')
            raw = np.loadtxt(io.BytesIO(archive.extractfile(raw_name).read()))
            agrees = (np.array_equal(raw[:,2], np.repeat(counts, nt))
                      and np.array_equal(raw[:,3:11], accumulators))
            rows.append(dict(shard=Path(name).stem[-2:], checkpointed_replicas=header[8],
                             checkpointed_survivors=header[9], escaped=header[10],
                             fixed_shard_complete=bool(header[11] and header[8] == 2_000_000),
                             raw_rows_match_checkpoint=agrees))
    assert len(rows) == 8
    receipt = dict(snapshot_utc=manifest['snapshot_utc'], archive=str(path.relative_to(ROOT)),
                   archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                   verified_archive_files=len(manifest['files']), planned_shards=8,
                   planned_replicas=16_000_000,
                   checkpointed_replicas=sum(row['checkpointed_replicas'] for row in rows),
                   completed_shards=sum(row['fixed_shard_complete'] for row in rows),
                   included_in_publication=False,
                   caveat='Recovery/progress snapshot only. No inference from partial fixed cohorts or runtime-selected prefixes.',
                   publication_outputs_modified=False, records=rows)
    out = ROOT/'reports/v2/hpc/recovery/original16m_final_snapshot.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({key: value for key, value in receipt.items() if key != 'records'}, indent=2))


if __name__ == '__main__':
    main()
