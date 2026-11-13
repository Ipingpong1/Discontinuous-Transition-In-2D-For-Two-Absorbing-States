#!/usr/bin/env python3
"""Explicit, recorded GPU restart after verified storage/checkpoint recovery."""
from pathlib import Path
from datetime import datetime, timezone
import json
import subprocess

base=Path('/scratch/mo2273/dn/v2/campaign/v2/hpc/gpu')
report=Path('/scratch/mo2273/dn/v2/gpu/resubmission_20261103.json')
assert not report.exists(), 'Refusing duplicate submission'
assert Path('/scratch/mo2273/dn/v2/gpu/quota_recovery_20261103.json').exists()
groups=[('orientation','production_manifest_profiles50k.json',2,'orient'),
        ('flagship_bulk','production_manifest_profiles50k.json',2,'bulk'),
        ('near_voter','production_manifest_profiles50k.json',2,'near'),
        ('orientation','production_manifest_largeL.json',3,'16korient'),
        ('flagship_bulk','production_manifest_largeL.json',3,'16kbulk'),
        ('near_voter','production_manifest_near_large.json',3,'near8192')]
records=[]
for family,manifest,workers,label in groups:
    for worker in range(1,workers+1):
        cmd=['sbatch','--parsable',f'--job-name=dnv2r_{label}_{worker}',
             '--export=ALL,DNSIM_MANIFEST='+str(base/manifest)+',DNSIM_RUNNER='+str(base/'production_runner_v2.py'),
             str(base/'production.sbatch'),family]
        result=subprocess.run(cmd,text=True,capture_output=True,check=True)
        job=int(result.stdout.strip().split(';')[0])
        record=dict(job=job,family=family,manifest=manifest,worker=worker,
                    utc=datetime.now(timezone.utc).isoformat(),command=cmd)
        records.append(record)
        report.write_text(json.dumps(records,indent=2)+'\n')
        print(json.dumps(record),flush=True)
