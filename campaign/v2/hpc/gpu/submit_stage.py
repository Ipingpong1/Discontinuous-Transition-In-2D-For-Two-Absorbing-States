#!/usr/bin/env python3
"""Submit the audited HOME/node-local recovery workers once, recording each ID."""
from datetime import datetime, timezone
from pathlib import Path
import json
import shlex
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RECORD = ROOT / 'reports/v2/hpc/staged_submission_20261103.json'
HOME_SOURCE = '/home/mo2273/dn_v2_stage/campaign/v2/hpc/gpu'
plan = [
    ('orientation', 'production_manifest_profiles50k.json', 'orient', 1),
    ('flagship_bulk', 'production_manifest_profiles50k.json', 'bulk', 2),
    ('near_voter', 'production_manifest_profiles50k.json', 'near', 2),
    ('orientation', 'production_manifest_largeL.json', '16korient', 3),
    ('flagship_bulk', 'production_manifest_largeL.json', '16kbulk', 3),
    ('near_voter', 'production_manifest_near_large.json', 'near8192', 3),
    ('wetting', 'wetting_manifest.json', 'wetting', 4),
]
if RECORD.exists():
    raise SystemExit(f'Already submitted or partially submitted: inspect {RECORD}')
record = dict(storage_root='/home/mo2273/dn_v2_stage/gpu',
              active_storage='node-local SLURM_TMPDIR',
              deadline_utc='2026-11-03T11:00:00+00:00',
              pilot_job=21725768, pilot_family='orientation',
              pilot_manifest='production_manifest_profiles50k.json',
              prerequisite='Pilot bitwise validation and validated t=1000000 checkpoint; 100% measured GPU utilization',
              jobs=[])
def save():
    RECORD.parent.mkdir(parents=True, exist_ok=True)
    RECORD.write_text(json.dumps(record, indent=2)+'\n')
save()
for family, manifest, label, count in plan:
    for index in range(1, count+1):
        command = ['sbatch', '--parsable', f'--job-name=dnv3_{label}_{index}',
                   f'--export=ALL,DNSIM_MANIFEST={HOME_SOURCE}/{manifest}',
                   f'{HOME_SOURCE}/production_stage.sbatch', family]
        result = subprocess.run(['ssh', 'torch', shlex.join(command)],
                                check=True, text=True, capture_output=True)
        job = int(result.stdout.strip().split(';')[0])
        record['jobs'].append(dict(job=job, family=family, manifest=manifest,
                                  worker=index, utc=datetime.now(timezone.utc).isoformat(),
                                  command=command))
        save()
        print(job, family, manifest, flush=True)
metadata_path = HERE / 'submitted.json'
metadata = json.loads(metadata_path.read_text())
metadata['staged_recovery'] = record
metadata_path.write_text(json.dumps(metadata, indent=2)+'\n')
