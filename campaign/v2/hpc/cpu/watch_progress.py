#!/usr/bin/env python3
"""Read-only CPU campaign watcher; local sleep, SSH metadata queries only."""
from pathlib import Path
from datetime import datetime,timezone
import subprocess,json,time
REMOTE=r'''
from pathlib import Path
import struct,json,collections,subprocess
r=Path('/scratch/mo2273/dn/v2/cpu');life=collections.Counter();started=collections.Counter()
for p in list((r/'out/lifetimes').glob('*.txt'))+list((r/'out/nearline').glob('*.txt')):
 key=p.stem.split('_L')[0]+'_L'+p.stem.split('_L')[1].split('_')[0]
 started[key]+=1
 with p.open('rb') as f:f.seek(max(0,p.stat().st_size-4096));done=b'# done t_end' in f.read()
 if done:life[key]+=1
def tail_progress(root):
 tails=[]
 for p in sorted((root/'checkpoints').glob('*.ckpt')):
  if not p.name.startswith(('tail_', 'confirm8m_')):continue
  with p.open('rb') as f:h=f.read(120)
  if len(h)!=120:continue
  a=struct.unpack('<12Q2d',h[8:]);tails.append(dict(shard=p.stem[-2:],n=a[8],survivors=a[9],escaped=a[10],done=bool(a[11])))
 for p in sorted((root/'out/tails').glob('*.txt')):
  lines=p.read_text().splitlines()
  footer=[line for line in lines if line.startswith('# done replicas ')]
  if footer:
   row=footer[-1].split();shard=p.stem[-2:]
   tails=[item for item in tails if item['shard']!=shard]
   tails.append(dict(shard=shard,n=int(row[3]),survivors=int(row[5]),escaped=int(row[7]),done=True))
 return sorted(tails,key=lambda item:item['shard'])
tails=tail_progress(r);confirm8m=tail_progress(r/'confirm8m')
states=[line for line in subprocess.check_output(['squeue','-h','-u','mo2273','-o','%i|%j|%T|%R'],text=True).splitlines() if '|dn-v2-' in line]
print(json.dumps(dict(lifetimes=dict(life),started=dict(started),tails=tails,confirm8m=confirm8m,states=states)))
'''
end=datetime(2026,11,3,11,26,tzinfo=timezone.utc).timestamp()
log=Path(__file__).resolve().parents[4]/'reports/v2/hpc/cpu_progress.jsonl';log.parent.mkdir(parents=True,exist_ok=True)
previous=None
while time.time()<end:
 try:
  result=subprocess.run(['ssh','torch','python3','-'],input=REMOTE,capture_output=True,text=True,timeout=30,check=True)
  data=json.loads(result.stdout);key=json.dumps(data,sort_keys=True)
  if key!=previous:
   data['utc']=datetime.now(timezone.utc).isoformat();line=json.dumps(data,sort_keys=True)
   with log.open('a') as f:f.write(line+'\n')
   print(line,flush=True);previous=key
 except Exception as error:print(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),error=str(error))),flush=True)
 time.sleep(60)
