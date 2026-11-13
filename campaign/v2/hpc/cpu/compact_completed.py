#!/usr/bin/env python3
"""Keep CPU staging within the home inode quota without dropping scientific data.

Only checkpoints of validated, completed tasks are packed. Every old and new
archive member is SHA256 checked before any original checkpoint is removed.
Completed zero-byte task logs are removed; nonempty logs and all outputs remain.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, io, json, re, struct, tarfile


def completed(task):
    out = Path(task['output'])
    if not out.is_file():
        return False
    lines = out.read_text().splitlines()
    done = [line for line in lines if line.startswith('# done ')]
    if len(done) != 1:
        return False
    data = [line for line in lines if line and not line.startswith('#')]
    if not data:
        return False
    if 't_end' in done[0]:
        match = re.search(r't_end (\d+) absorbed ([01])', done[0])
        if not match or int(data[-1].split()[0]) != int(match[1]):
            return False
        rho = float(data[-1].split()[1])
        return match[2] == '0' or rho in (0.0, 1.0)
    match = re.search(r'replicas (\d+) survivors_T (\d+) escaped (\d+)', done[0])
    if not match or int(match[1]) != task.get('replicas', 2_000_000):
        return False
    final = {}
    for line in data:
        row = line.split()
        final[int(row[0])] = row
    return (sum(int(row[2]) for row in final.values()) == int(match[1])
            and sum(int(row[3]) for row in final.values()) == int(match[2])
            and sum(int(row[4]) for row in final.values()) == int(match[3]))


def compact(root):
    root = root.resolve()
    manifest = json.loads((root/'manifest.json').read_text())
    tasks = [task for group in manifest['tasks'].values() for task in group]
    complete = [task for task in tasks if completed(task)]
    archive = root/'completed_checkpoints.tar.gz'
    metadata = root/'completed_checkpoints.json'
    members = {}
    if archive.exists():
        previous = json.loads(metadata.read_text())['members']
        with tarfile.open(archive, 'r:gz') as tar:
            for member in tar:
                assert member.isfile() and '/' not in member.name
                content = tar.extractfile(member).read()
                assert hashlib.sha256(content).hexdigest() == previous[member.name]['sha256']
                members[member.name] = content
        assert set(previous) == set(members)
    originals = []
    finalized = set()
    for task in complete:
        args = task['args']
        path = Path(args[args.index('--ckpt')+1])
        if not path.exists():
            if path.name in members:
                finalized.add(task['name'])
            continue
        content = path.read_bytes()
        assert content[:8] in (b'DNSIM1\0\0', b'DNREPS1\0'), str(path)
        # The simulator writes its completion footer before its final checkpoint.
        # Do not race that final write or remove a still-open empty task log.
        if content[:8] == b'DNSIM1\0\0':
            end = int(re.search(r'# done t_end (\d+)', Path(task['output']).read_text())[1])
            if struct.unpack('<Q', content[24:32])[0] != end:
                continue
        else:
            header = struct.unpack('<12Q2d', content[8:120])
            if header[8] != task.get('replicas', 2_000_000) or not header[11]:
                continue
        if path.name in members:
            assert members[path.name] == content, str(path)
        members[path.name] = content
        originals.append(path)
        finalized.add(task['name'])
    if originals:
        temporary = archive.with_suffix('.gz.tmp')
        with tarfile.open(temporary, 'w:gz') as tar:
            for name, content in sorted(members.items()):
                info = tarfile.TarInfo(name)
                info.size = len(content)
                tar.addfile(info, io.BytesIO(content))
        hashes = {name: dict(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
                  for name, content in members.items()}
        with tarfile.open(temporary, 'r:gz') as tar:
            checked = {}
            for member in tar:
                checked[member.name] = hashlib.sha256(tar.extractfile(member).read()).hexdigest()
            assert checked == {name: info['sha256'] for name, info in hashes.items()}
        temporary.replace(archive)
        metadata.write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),
                                           members=hashes), indent=2)+'\n')
        # Recheck source bytes after archive verification, before removing them.
        for path in originals:
            assert hashlib.sha256(path.read_bytes()).hexdigest() == hashes[path.name]['sha256']
            path.unlink()
    empty_logs = []
    for task in complete:
        if task['name'] not in finalized:
            continue
        path = Path(task['log'])
        if path.is_file() and path.stat().st_size == 0:
            path.unlink()
            empty_logs.append(path.name)
    print(json.dumps(dict(completed=len(complete), checkpoints_packed=len(originals),
                          archive_members=len(members), empty_logs_removed=empty_logs,
                          tree_inodes=1+sum(1 for _ in root.rglob('*')))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/scratch/mo2273/dn/v2/cpu'))
    args = parser.parse_args()
    compact(args.root)
