#!/usr/bin/env python3
"""Periodically compact closed CPU checkpoints and mirror scientific outputs."""
from datetime import datetime, timezone
from pathlib import Path
import json, subprocess, time

ROOT = Path(__file__).resolve().parents[4]
DEADLINE = datetime(2026, 11, 3, 11, 26, tzinfo=timezone.utc).timestamp()
LOG = ROOT/'reports/v2/hpc/cpu_maintenance.jsonl'
commands = [
    ['ssh', 'torch', 'python3', '/home/mo2273/dn_v2_stage/cpu/compact_completed.py',
     '--root', '/home/mo2273/dn_v2_stage/cpu/confirm8m'],
    ['ssh', 'torch', 'python3', '/home/mo2273/dn_v2_stage/cpu/compact_completed.py'],
    ['rsync', '-a', '--include=*/', '--include=*.txt', '--exclude=*',
     'torch:/home/mo2273/dn_v2_stage/cpu/out/', str(ROOT/'data/v2/cpu/')+'/'],
    ['rsync', '-a', '--include=*/', '--include=*.txt', '--exclude=*',
     'torch:/home/mo2273/dn_v2_stage/cpu/confirm8m/out/', str(ROOT/'data/v2/cpu/confirm8m/')+'/'],
    ['python3', '-m', 'campaign.v2.analysis.cpu_lifetimes'],
    ['python3', '-m', 'campaign.v2.analysis.seed_tails'],
]
while True:
    for command in commands:
        try:
            result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, timeout=120)
            event = dict(utc=datetime.now(timezone.utc).isoformat(), command=command,
                         returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
        except Exception as error:
            event = dict(utc=datetime.now(timezone.utc).isoformat(), command=command, error=str(error))
        with LOG.open('a') as stream:
            stream.write(json.dumps(event)+'\n')
        print(json.dumps(event), flush=True)
    if time.time() >= DEADLINE:
        break
    # Local scheduler; tools remain free to report progress or act on new input.
    time.sleep(min(300, DEADLINE-time.time()))
