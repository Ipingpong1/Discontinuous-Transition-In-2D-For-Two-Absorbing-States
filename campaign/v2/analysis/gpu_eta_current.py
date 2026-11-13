"""Remaining independent5M work from the current four-worker event snapshot.

Read-only scheduling audit: does not change jobs, priorities, or deadlines.
The earlier front5m_eta.json remains a historical pending-slot prediction.
"""
from datetime import datetime,timezone
import heapq,json,statistics
from .io import ROOT,write_json


def main():
    folder=ROOT/'reports/v2/hpc'
    snap=json.loads((folder/'gpu_eta_current_snapshot.json').read_text())
    now=datetime.fromisoformat(snap['snapshot_utc']).timestamp()
    tasks=json.loads((ROOT/'campaign/v2/hpc/gpu/production_manifest_orientation5m.json').read_text())['tasks']
    identifiers={t['id'] for t in tasks}
    complete=set(snap['bundles']) & identifiers
    durations=[];active={}
    for job in range(21727284,21727288):
        worker=snap['workers'].get(str(job),{})
        for ident,run in worker.get('runs',{}).items():
            if ident not in identifiers:continue
            if 'end' in run:
                if run['end'].get('t_end')==5000000 and not run['end'].get('absorbed'):
                    durations.append(run['end']['wall_s'])
            elif ident not in complete:
                active[job]=(ident,datetime.fromisoformat(run['start']['utc']).timestamp())
    if not durations:raise SystemExit('No completed five-million task duration')
    typical=statistics.median(durations)
    remaining=[t for t in tasks if t['id'] not in complete]
    reports=[]
    for factor in (1.,1.2):
        queue=[t for t in remaining if t['id'] not in {v[0] for v in active.values()}]
        heap=[];finish_times={}
        for job,(ident,start) in active.items():
            finish=max(now+5,start+(typical+1.5)*factor)
            finish_times[ident]=finish;heapq.heappush(heap,(finish,job))
        # A currently idle worker is available only if Slurm still reports RUNNING.
        for job in range(21727284,21727288):
            if job not in active and snap.get('slurm_states',{}).get(str(job))=='RUNNING':heapq.heappush(heap,(now+5,job))
        if queue and not heap:raise SystemExit('Remaining tasks but no running5m worker; no ETA extrapolation')
        while queue:
            available,job=heapq.heappop(heap);task=queue.pop(0)
            finish=available+(typical+1.5)*factor;finish_times[task['id']]=finish;heapq.heappush(heap,(finish,job))
        cohorts={}
        for task in tasks:
            key=str(task['orientation'])+' q='+task['q']
            cohorts.setdefault(key,[]).append(now if task['id'] in complete else finish_times[task['id']])
        limit=datetime.fromisoformat('2026-11-03T11:00:00+00:00').timestamp()
        reports.append(dict(duration_multiplier=factor,
            all5m_finish_estimate=datetime.fromtimestamp(max(finish_times.values(),default=now),timezone.utc).isoformat(),
            cohorts_complete_by11=sum(max(v)<=limit for v in cohorts.values()),
            cohort_finish_estimates={k:datetime.fromtimestamp(max(v),timezone.utc).isoformat() for k,v in cohorts.items()}))
    result=dict(snapshot_utc=snap['snapshot_utc'],n_complete=len(complete),n_remaining=len(remaining),
        active_tasks={str(j):{'id':v[0],'start_utc':datetime.fromtimestamp(v[1],timezone.utc).isoformat()} for j,v in active.items()},
        completed_task_durations_s=durations,median_task_s=typical,
        remaining_nominal_gpu_hours=len(remaining)*(typical+1.5)/3600,
        scenarios=reports,assumptions=['The four currently running declared workers continue unchanged.','Task duration uses the median of completed5M tasks; alternative scenario increases this by20%.','An active predicted completion already passed is lower-bounded by five seconds from the snapshot.','No job, deadline, priority, or manifest is modified.'])
    write_json(folder/'front5m_eta_current.json',result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
