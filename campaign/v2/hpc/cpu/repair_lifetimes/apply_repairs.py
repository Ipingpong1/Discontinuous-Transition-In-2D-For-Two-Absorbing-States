#!/usr/bin/env python3
"""Validate deterministic same-seed reruns before replacing quota-gapped output."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, math, re, tarfile

ROOT = Path('/home/mo2273/dn_v2_stage/cpu/repair_lifetimes')
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--apply', action='store_true')
args = parser.parse_args()
manifest = json.loads((ROOT/'manifest.json').read_text())
preserved = json.loads((ROOT/'gapped_originals.json').read_text())
sha = lambda b: hashlib.sha256(b).hexdigest()
audit_path = ROOT/'repairs_applied.json'
if audit_path.exists():
    audit = json.loads(audit_path.read_text())
    assert all(sha(Path(row['canonical']).read_bytes()) == row['repaired_sha256'] for row in audit['repairs'])
    print('Repairs already applied and verified')
    raise SystemExit


def rows_and_footer(content):
    text = content.decode()
    footer = re.search(r'^# done t_end (\d+) absorbed ([01])', text, re.M)
    if not footer:
        return None
    rows = {}
    for line in text.splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        assert len(fields) == 5
        time = int(fields[0])
        assert time not in rows
        rows[time] = line.strip()
    end, absorbed = map(int, footer.groups())
    assert max(rows) == end
    return rows, end, absorbed


pending = []
for task in manifest['tasks']['lifetimes']:
    path = Path(task['output'])
    if not path.exists() or rows_and_footer(path.read_bytes()) is None:
        pending.append(task['name'])
if pending:
    print(json.dumps(dict(ready=False, pending=pending)))
    raise SystemExit

records = []
with tarfile.open(ROOT/'gapped_originals.tar.gz', 'r:gz') as damaged, \
        tarfile.open(ROOT.parent/'completed_checkpoints.tar.gz', 'r:gz') as checkpoints:
    for task in manifest['tasks']['lifetimes']:
        original = Path(task['original_output'])
        old = damaged.extractfile(original.name).read()
        old_record = next(row for row in preserved['originals'] if row['task'] == task['name'])
        assert sha(old) == old_record['sha256']
        assert sha(original.read_bytes()) == sha(old), 'Canonical output changed since preservation'
        new = Path(task['output']).read_bytes()
        old_rows, old_end, old_absorbed = rows_and_footer(old)
        new_rows, end, absorbed = rows_and_footer(new)
        assert (end, absorbed) == (old_end, old_absorbed), task['name']
        every = int(task['args'][task['args'].index('--every')+1])
        ppd = int(task['args'][task['args'].index('--ppd')+1])
        expected = set(range(0, end+1, every))|{end}
        expected |= {t for d in range(11) for i in range(ppd)
                     if (t := math.floor(10.**(d+i/ppd)+.5)) <= end}
        assert expected <= new_rows.keys(), task['name']
        assert all(new_rows.get(time) == line for time, line in old_rows.items()), task['name']
        expected_checkpoint = checkpoints.extractfile(Path(task['original_checkpoint']).name).read()
        actual_checkpoint = Path(task['checkpoint']).read_bytes()
        assert actual_checkpoint == expected_checkpoint, 'Final checkpoint differs: '+task['name']
        records.append(dict(task=task['name'], canonical=str(original), rerun=task['output'],
                            original_sha256=sha(old), repaired_sha256=sha(new),
                            overlapping_records_verified=len(old_rows),
                            restored_records=len(new_rows)-len(old_rows),
                            final_checkpoint_sha256=sha(actual_checkpoint),
                            end=end, absorbed=absorbed))

audit = dict(utc=datetime.now(timezone.utc).isoformat(), repairs=records,
             original_archive=preserved, applied=args.apply,
             method='Same unchanged binary and seed; every overlapping numerical row and final checkpoint byte agrees. Full regular observation grid restored.')
if args.apply:
    for record in records:
        original = Path(record['canonical'])
        temporary = original.with_suffix('.repaired.tmp')
        temporary.write_bytes(Path(record['rerun']).read_bytes())
        assert sha(temporary.read_bytes()) == record['repaired_sha256']
        temporary.replace(original)
    audit_path.write_text(json.dumps(audit, indent=2)+'\n')
print(json.dumps(dict(ready=True, applied=args.apply, repairs=records)))
