"""Read-only completed-bundle synchronization and sequential local reductions.

Does not submit, cancel, or mutate anything on Torch. Run one instance only.
The scratch source is deliberately absent: it was frozen after quota failure.
"""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from .io import ROOT,manifest_tasks,manifest_groups,write_json

stopping=False

def stop(*_):
    global stopping
    stopping=True


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--interval',type=float,default=300.)
    ap.add_argument('--deadline',default='2026-11-03T11:10:00+00:00')
    args=ap.parse_args();deadline=datetime.fromisoformat(args.deadline).timestamp()
    root=ROOT/'data/v2/gpu';root.mkdir(parents=True,exist_ok=True)
    log=root/'maintenance.jsonl';pid=root/'maintenance.pid'
    if pid.exists():
        try:os.kill(int(pid.read_text()),0)
        except (ProcessLookupError,ValueError):pass
        else:raise SystemExit('A GPU maintenance process is already active')
    pid.write_text(str(os.getpid())+'\n')
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    def record(**data):
        data=dict(utc=datetime.now(timezone.utc).isoformat(),**data)
        with log.open('a') as stream:stream.write(json.dumps(data)+'\n')
        print(json.dumps(data),flush=True)
    def run(command,expected_wait=()):
        result=subprocess.run(command,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=1200)
        if result.returncode and expected_wait and result.stdout.startswith(expected_wait):
            print(result.stdout.strip(),flush=True);return True
        if result.returncode:record(event='command_error',command=command,returncode=result.returncode,output=result.stdout[-5000:])
        elif result.stdout.strip():print(result.stdout.strip(),flush=True)
        return result.returncode==0
    try:
        record(event='start',pid=os.getpid(),deadline=args.deadline,interval=args.interval)
        first_cycle=True
        while not stopping and time.time()<deadline:
            start=time.monotonic()
            ledger=root/'bundle_imports.json'
            before=json.loads(ledger.read_text()) if ledger.exists() else {}
            ok=run(['rsync','-a',"--exclude=*.tmp.*",'torch:/home/mo2273/dn_v2_stage/gpu/bundles/',str(root/'bundles')+'/'])
            if ok:ok=run([sys.executable,'-m','campaign.v2.analysis.import_gpu_bundles'])
            after=json.loads(ledger.read_text()) if ledger.exists() else {}
            changed=sorted(k for k,v in after.items() if before.get(k,{}).get('bundle_sha256')!=v['bundle_sha256'])
            if ok and (changed or first_cycle):
                ok=run([sys.executable,'-m','campaign.v2.analysis.gpu_campaign','--boot','1000'])
                if ok:ok=run([sys.executable,'-m','campaign.v2.analysis.gpu_report'])
                if ok:ok=run([sys.executable,'-m','campaign.v2.analysis.axial10m_diagnostics'])
                if ok:ok=run([sys.executable,'-m','campaign.v2.analysis.gpu_wetting','--boot','2000'])
                if ok:ok=run([sys.executable,str(ROOT/'paper/v2/figs/fig_gpu_v2.py')])
            # Root may fetch missing verified final grids without a new bundle.
            if ok:ok=run([sys.executable,'-m','campaign.v2.analysis.spatial_correlations'])
            if ok:ok=run([sys.executable,str(ROOT/'paper/v2/figs/fig_spatial_v2.py')],expected_wait=('Spatial comparison awaits complete cohorts:','Run spatial_correlations first'))
            status_path=ROOT/'campaign/v2/results/gpu/status.json'
            status=json.loads(status_path.read_text()) if status_path.exists() else {}
            wetting=manifest_tasks([ROOT/'campaign/v2/hpc/gpu/wetting_manifest.json'])
            wetting_groups=[]
            for key,ids in manifest_groups(wetting).items():
                if all(i in after for i in ids):wetting_groups.append(dict(p=key[0],q=key[2],Lx=key[3],Ly=key[4],T=key[5],n=len(ids)))
            record(event='cycle',ok=ok,new_tasks=changed,complete=status.get('n_complete'),
                complete_ensembles=status.get('n_complete_ensembles'),complete_wetting_groups=wetting_groups,
                elapsed_s=round(time.monotonic()-start,1))
            if (ok and status.get('n_complete')==status.get('n_planned')
                    and len(wetting_groups)==len(manifest_groups(wetting))):
                record(event='all_declared_gpu_runs_imported_stop',complete=status['n_complete'],wetting_groups=len(wetting_groups))
                break
            first_cycle=False
            # Short interruptible waits; no foreground tool call blocks here.
            until=min(time.time()+max(0,args.interval-(time.monotonic()-start)),deadline)
            while not stopping and time.time()<until:time.sleep(min(1,until-time.time()))
    finally:
        pid.unlink(missing_ok=True);record(event='stop',pid=os.getpid())

if __name__=='__main__':main()
