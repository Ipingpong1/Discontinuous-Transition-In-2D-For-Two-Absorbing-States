#!/usr/bin/env python3
"""Pause only our two surviving CPU jobs, byte-verify relocation, then resume.

Never move a tree while a simulator has an atomic checkpoint/output file open.
Old trees remain for stdout/stderr descriptors and subsequent log collection.
Run on Torch login only after the coordinator confirms home inode headroom.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,os,subprocess,time
SOURCE=Path('/scratch/mo2273/dn/v2/cpu')
DEST=Path('/home/mo2273/dn_v2_stage/cpu')
BACKUP=SOURCE.with_name('cpu_scratch_before_home_20261103')
LINK=SOURCE.with_name('cpu_home_link_20261103')
JOBS=[21723013]
HELPER=r'''
from pathlib import Path
import subprocess,json,os,signal,time
job=JOB;action=ACTION;expected=EXPECTED
selected=[]
for line in subprocess.check_output(['ps','-u',str(os.getuid()),'-o','pid=,args='],text=True).splitlines():
    pid=int(line.strip().split(None,1)[0]);p=Path('/proc')/str(pid)
    try:
        args=p.joinpath('cmdline').read_bytes().replace(b'\0',b' ').decode()
        cg=p.joinpath('cgroup').read_text()
    except OSError:continue
    if f'job_{job}/' not in cg or '/scratch/mo2273/dn/v2/cpu/' not in args:continue
    if 'dnsim_v15' not in args and 'runner.py' not in args:continue
    if expected and pid not in expected:continue
    selected.append(pid)
if action=='stop':
    for pid in selected:os.kill(pid,signal.SIGSTOP)
    time.sleep(.2)
if action=='continue':
    for pid in selected:os.kill(pid,signal.SIGCONT)
files={};states={}
for pid in selected:
    p=Path('/proc')/str(pid)
    try:
        files[pid]={f.name:os.readlink(f) for f in p.joinpath('fd').iterdir()}
        states[pid]=[s for s in p.joinpath('status').read_text().splitlines() if s.startswith('State:')][0]
    except OSError:pass
visible={pid:Path(f"/proc/{pid}/root/home/mo2273/dn_v2_stage/cpu").is_dir() for pid in selected}
print(json.dumps(dict(pids=selected,files=files,states=states,destination_visible=visible)))
'''

def node(job,action,expected=()):
    code=HELPER.replace('job=JOB;action=ACTION;expected=EXPECTED',f'job={job};action={action!r};expected={list(expected)!r}')
    r=subprocess.run(['srun',f'--jobid={job}','--overlap','--nodes=1','--ntasks=1','python3','-'],
        input=code,text=True,capture_output=True,timeout=30,check=True)
    return json.loads(r.stdout)

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--apply',action='store_true');args=ap.parse_args()
if not args.apply:
    print(json.dumps(dict(source=str(SOURCE),destination=str(DEST),jobs=JOBS,apply=False)));raise SystemExit
if SOURCE.is_symlink():
    assert SOURCE.resolve()==DEST.resolve();print('Already relocated');raise SystemExit
assert SOURCE.is_dir() and not BACKUP.exists()
DEST.mkdir(parents=True,exist_ok=True)
# Create this before pausing anything: a quota failure must not strand live work.
if not LINK.exists():LINK.symlink_to(DEST,target_is_directory=True)
assert LINK.is_symlink() and LINK.resolve()==DEST.resolve()
paused={};audit=dict(utc=datetime.now(timezone.utc).isoformat(),source=str(SOURCE),destination=str(DEST),backup=str(BACKUP),jobs=JOBS)
try:
    for job in JOBS:
        state=node(job,'stop');paused[job]=state['pids']
        assert state['pids'],f'No matching live processes in job{job}'
        assert all(state['destination_visible'].values()),state
        assert all('T (stopped)' in x for x in state['states'].values()),state
        open_science=[v for fds in state['files'].values() for v in fds.values()
            if str(SOURCE)+'/' in v and not v.endswith(('.log','.slurm'))]
        assert not open_science,f'Atomic or event file open; resume and retry: {open_science}'
    audit['paused']=paused
    subprocess.run(['rsync','-a','--delete',str(SOURCE)+'/',str(DEST)+'/'],check=True)
    records=[]
    for path in sorted(SOURCE.rglob('*')):
        if path.is_dir():continue
        assert path.is_file() and not path.is_symlink(),str(path)
        target=DEST/path.relative_to(SOURCE)
        digest=sha(path);assert target.is_file() and sha(target)==digest,str(target)
        records.append(dict(path=str(path.relative_to(SOURCE)),bytes=path.stat().st_size,sha256=digest))
    audit['verified_files']=records
    SOURCE.rename(BACKUP)
    try:LINK.rename(SOURCE)
    except BaseException:
        BACKUP.rename(SOURCE);raise
    assert SOURCE.resolve()==DEST.resolve()
    audit['relocated']=True
    (DEST/'home_relocation_20261103.json').write_text(json.dumps(audit,indent=2)+'\n')
finally:
    for job,pids in paused.items():
        try:
            state=node(job,'continue',pids)
            print(json.dumps(dict(job=job,resumed=state)),flush=True)
        except Exception as error:print(json.dumps(dict(job=job,resume_error=str(error))),flush=True)
print(json.dumps(dict(relocated=audit.get('relocated',False),verified_files=len(audit.get('verified_files',[])),backup=str(BACKUP))))
