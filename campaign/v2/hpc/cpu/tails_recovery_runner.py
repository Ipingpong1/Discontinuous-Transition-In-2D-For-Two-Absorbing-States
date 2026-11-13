#!/usr/bin/env python3
"""Shared fixed-N tail recovery queue; retains original manifest and seeds."""
import concurrent.futures as cf
from datetime import datetime, timezone
from pathlib import Path
import hashlib,json,os,signal,subprocess,threading,time
ROOT=Path('/scratch/mo2273/dn/v2/cpu')
DEADLINE=datetime(2026,11,3,11,25,tzinfo=timezone.utc).timestamp()
BIN='/scratch/mo2273/dn/bin2/dnsim_v15'
JOB=os.environ.get('SLURM_JOB_ID','local')
WORKERS=int(os.environ.get('DNSIM_CPU_WORKERS','3'))
EXCLUDE=set(os.environ.get('DNSIM_EXCLUDE_SHARDS','02').split(','))
LOCK=threading.Lock()
EVENTS=ROOT/'logs'/f'tails_recovery_{JOB}.jsonl'
CLAIMS=ROOT/'claims'/'tails_recovery'
CLAIMS.mkdir(parents=True,exist_ok=True)

def record(**obj):
    obj.update(utc=datetime.now(timezone.utc).isoformat(),job=JOB)
    with LOCK:
        with EVENTS.open('a') as f:f.write(json.dumps(obj,sort_keys=True)+'\n')
        print(json.dumps(obj,sort_keys=True),flush=True)

def complete(task):
    out=Path(task['output'])
    if not out.exists():return False
    for line in out.read_text().splitlines()[-5:]:
        if line.startswith('# done replicas '):
            return int(line.split()[3])==2_000_000
    return False

def run_one(task):
    if task['name'][-2:] in EXCLUDE or complete(task):return
    if time.time()>=DEADLINE-60:return
    claim=CLAIMS/task['name']
    try:claim.mkdir()
    except FileExistsError:return
    owner=claim/'owner.json';owner.write_text(json.dumps(dict(job=JOB,task=task['name']))+'\n')
    started=time.time();success=False
    try:
        record(status='start',task=task['name'],command=[BIN]+task['args'])
        env=dict(os.environ,OMP_NUM_THREADS='1',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1')
        with Path(task['log']).open('a') as log:
            process=subprocess.Popen([BIN]+task['args'],stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
            try:rc=process.wait(timeout=max(1,DEADLINE-time.time()))
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGTERM)
                try:rc=process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL);rc=process.wait()
        success=rc==0 and complete(task)
        record(status='finished',task=task['name'],returncode=rc,complete=success,wall_s=time.time()-started)
    finally:
        if not success:
            owner.unlink(missing_ok=True);claim.rmdir()

manifest=json.loads((ROOT/'manifest.json').read_text())
sha=hashlib.sha256(Path(BIN).read_bytes()).hexdigest()
assert sha==manifest['binary_sha256']
record(status='job_start',workers=WORKERS,excluded_shards=sorted(EXCLUDE),binary_sha256=sha,
       runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
with cf.ThreadPoolExecutor(max_workers=WORKERS) as pool:
    for future in cf.as_completed([pool.submit(run_one,t) for t in manifest['tasks']['tails']]):future.result()
record(status='job_end')
