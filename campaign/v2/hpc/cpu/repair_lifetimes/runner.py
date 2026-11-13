#!/usr/bin/env python3
import concurrent.futures as cf
from datetime import datetime,timezone
import hashlib,json,os,signal,subprocess,sys,threading,time
from pathlib import Path

ROOT=Path("/home/mo2273/dn_v2_stage/cpu/repair_lifetimes")
DEADLINE=datetime(2026,11,3,11,25,tzinfo=timezone.utc).timestamp()
BIN="/scratch/mo2273/dn/bin2/dnsim_v15"
lock=threading.Lock()

def record(path,obj):
    with lock:
        with path.open("a") as f: f.write(json.dumps(obj,sort_keys=True)+"\n")
        print(json.dumps(obj,sort_keys=True),flush=True)

def run_one(task,events):
    out=Path(task["output"]);out.parent.mkdir(parents=True,exist_ok=True)
    if out.exists() and b"# done" in out.read_bytes()[-4096:]:
        record(events,dict(task=task["name"],status="already_complete"));return
    remaining=DEADLINE-time.time()
    if remaining<60:
        record(events,dict(task=task["name"],status="deadline_not_started"));return
    env=os.environ.copy();env.update(OMP_NUM_THREADS="1",OMP_PROC_BIND="false",OPENBLAS_NUM_THREADS="1")
    started=time.time()
    record(events,dict(task=task["name"],status="start",utc=datetime.now(timezone.utc).isoformat(),
                      command=[BIN]+task["args"]))
    with Path(task["log"]).open("a") as log:
        proc=subprocess.Popen([BIN]+task["args"],stdout=log,stderr=subprocess.STDOUT,
                              env=env,start_new_session=True)
        status="finished"
        try: rc=proc.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            status="deadline_checkpoint_preserved"
            os.killpg(proc.pid,signal.SIGTERM)
            try: rc=proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL);rc=proc.wait()
    record(events,dict(task=task["name"],status=status,returncode=rc,wall_s=time.time()-started,
                      output_bytes=out.stat().st_size if out.exists() else 0))

if __name__=="__main__":
    kind=sys.argv[1]
    manifest=json.loads((ROOT/"manifest.json").read_text())
    tasks=manifest["tasks"][kind]
    events=ROOT/"logs"/(kind+"_events.jsonl")
    events.parent.mkdir(parents=True,exist_ok=True)
    record(events,dict(status="job_start",kind=kind,job=os.environ.get("SLURM_JOB_ID"),
                      utc=datetime.now(timezone.utc).isoformat(),workers=4,
                      binary_sha256=hashlib.sha256(Path(BIN).read_bytes()).hexdigest()))
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(run_one,t,events) for t in tasks]
        for f in cf.as_completed(futures): f.result()
    record(events,dict(status="job_end",kind=kind,utc=datetime.now(timezone.utc).isoformat()))
