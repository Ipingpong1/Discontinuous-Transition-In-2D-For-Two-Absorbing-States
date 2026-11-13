#!/usr/bin/env python3
"""Read Slurm accounting only for this dated v2 campaign, excluding job steps."""
from datetime import datetime, timezone
from pathlib import Path
import csv
import io
import json
import re
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = ROOT/'reports/v2/hpc'
sources = [ROOT/'campaign/v2/hpc/gpu/submitted.json',
           ROOT/'campaign/v2/hpc/cpu/submitted.json',
           OUTPUT/'staged_submission_20261103.json']
jobs = {21724297,21724298,21724302,21724306}  # initial wetting attempts
def visit(value):
    if isinstance(value, dict):
        for child in value.values(): visit(child)
    elif isinstance(value, list):
        for child in value: visit(child)
    elif isinstance(value, int) or isinstance(value,str) and value.isdigit():
        number=int(value)
        if 21700000 <= number < 21800000: jobs.add(number)
for path in sources:
    if path.exists(): visit(json.loads(path.read_text()))
fields='JobIDRaw,JobName,Partition,State,Start,End,ElapsedRaw,AllocCPUS,TotalCPU,AllocTRES'
command=['sacct','-n','-P','-j',','.join(str(j) for j in sorted(jobs)),
         '--format='+fields]
result=subprocess.run(['ssh','torch',shlex.join(command)],check=True,
                      text=True,capture_output=True)
rows=[]
for line in result.stdout.splitlines():
    parts=line.split('|')
    if not parts or not parts[0].isdigit(): continue
    row=dict(zip(fields.split(','),parts))
    if int(row['JobIDRaw']) not in jobs: raise ValueError('Unexpected job')
    gres=re.search(r'(?:^|,)gres/gpu=(\d+)(?:,|$)',row['AllocTRES'])
    if not gres:
        gres=re.search(r'(?:^|,)gres/gpu:h200=(\d+)(?:,|$)',row['AllocTRES'])
    row['allocated_gpus']=int(gres.group(1)) if gres else 0
    seconds=int(row['ElapsedRaw'] or 0)
    row['gpu_hours']=row['allocated_gpus']*seconds/3600
    row['allocated_cpu_hours']=int(row['AllocCPUS'] or 0)*seconds/3600
    rows.append(row)
OUTPUT.mkdir(parents=True,exist_ok=True)
with (OUTPUT/'slurm_accounting.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(rows[0]) if rows else fields.split(','))
    writer.writeheader();writer.writerows(rows)
summary=dict(utc=datetime.now(timezone.utc).isoformat(),
             requested_job_ids=sorted(jobs),reported_jobs=len(rows),
             gpu_hours=sum(r['gpu_hours'] for r in rows),
             cpu_only_allocated_core_hours=sum(r['allocated_cpu_hours'] for r in rows if not r['allocated_gpus']),
             gpu_host_allocated_core_hours=sum(r['allocated_cpu_hours'] for r in rows if r['allocated_gpus']),
             unfinished_jobs=[r['JobIDRaw'] for r in rows if r['State'] in ('RUNNING','PENDING','COMPLETING')],
             definition='Slurm allocation elapsed time; includes validation, I/O, failed/restarted attempts. Job steps excluded to avoid double counting.',
             command=command)
(OUTPUT/'accounting_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ('requested_job_ids','command')},indent=2))
