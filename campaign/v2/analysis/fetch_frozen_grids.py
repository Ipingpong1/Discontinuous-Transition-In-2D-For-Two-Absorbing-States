"""Fetch only missing final grids for already-complete frozen flagship tasks.

No old text, profile, checkpoint, or metadata file is copied over resumed data.
Every fetched grid must match the pre-existing local completion hash/time.
"""
from datetime import datetime,timezone
from pathlib import Path
import json
import subprocess
import tempfile
from .gpu_campaign import MANIFESTS,task_completion
from .import_gpu_bundles import validate_final
from .io import ROOT,manifest_tasks,write_json

def main():
    base=ROOT/'data/v2/gpu/runs';selected=[]
    for task in manifest_tasks(MANIFESTS).values():
        if task['family']!='flagship_bulk' or task['Lx']<8192: continue
        directory=base/task['family']/task['id']
        if (directory/'final.bin').exists() or not (directory/'done.json').exists():continue
        status,_,_=task_completion(task,directory/'out.txt')
        if status['complete']:selected.append(task)
    if not selected:
        print('No missing frozen final grids');return
    with tempfile.NamedTemporaryFile(mode='w',prefix='dn-v2-final-grids-',suffix='.txt') as names:
        for task in selected:
            names.write(f"{task['family']}/{task['id']}/final.bin\n")
        names.flush()
        subprocess.run(['rsync','-a','--ignore-existing','--files-from='+names.name,
                        'torch:/scratch/mo2273/dn/v2/gpu/runs/',str(base)+'/'],check=True)
    records=[]
    for task in selected:
        checked=validate_final(base/task['family']/task['id'],task)
        records.append(dict(id=task['id'],**checked))
    path=ROOT/'reports/v2/hpc/frozen_grid_fetch.json'
    prior=json.loads(path.read_text()) if path.exists() else []
    prior.append(dict(utc=datetime.now(timezone.utc).isoformat(),
                      operation='Read-only fetch of missing final.bin only; ignore-existing; checked against pre-existing local done.json',
                      grids=records))
    write_json(path,prior)
    print(f'{len(records)} frozen final grids fetched and SHA/time verified')

if __name__=='__main__':main()
