"""Scheduling estimate from a read-only worker-event snapshot; never changes jobs."""
from datetime import datetime,timezone,timedelta
from pathlib import Path
import heapq,json,statistics
from .io import ROOT,write_json


def main():
    snapshot=json.loads((ROOT/'reports/v2/hpc/gpu_eta_snapshot.json').read_text())
    now=datetime.fromisoformat(snapshot['snapshot_utc']).timestamp()
    status=json.loads((ROOT/'campaign/v2/results/gpu/status.json').read_text())
    complete={i for i,s in status['tasks'].items() if s['complete']}|set(snapshot['bundles'])
    folder=ROOT/'campaign/v2/hpc/gpu'
    def tasks(name,family):
        return [t for t in json.loads((folder/name).read_text())['tasks'] if t['family']==family]
    pools=[('small_orientation',[21725768,21725964],tasks('production_manifest_profiles50k.json','orientation')),
           ('large_orientation',list(range(21725969,21725972)),tasks('production_manifest_largeL.json','orientation')),
           ('large_bulk',list(range(21725972,21725975)),tasks('production_manifest_largeL.json','flagship_bulk')),
           ('near_large',list(range(21725975,21725978)),tasks('production_manifest_near_large.json','near_voter')),
           ('wetting',list(range(21725978,21725982)),tasks('wetting_manifest.json','wetting'))]
    repair='orientation_L16384x16384_ab1-1_q.04596_s2'
    # Existing family workers skip this frozen marker; a separate repair is queued.
    complete.add(repair)
    timings={};active={}
    for job,w in snapshot['workers'].items():
        for ident,run in w['runs'].items():
            cmd=run['start']['command'];Lx=int(cmd[cmd.index('--Lx')+1]);Ly=int(cmd[cmd.index('--Ly')+1]);T=int(cmd[cmd.index('--T')+1])
            profile='--profile' in cmd
            if 'end' not in run:active[int(job)]=(ident,datetime.fromisoformat(run['start']['utc']).timestamp())
            elif run['end']['t_end']==T and not run['end']['absorbed']:
                recovery=run.get('recovery') or {}
                if recovery.get('resume')=='fresh_same_seed':
                    timings.setdefault((Lx,Ly,T,profile),[]).append(run['end']['wall_s'])
    def duration(task):
        key=tuple(task[k] for k in ('Lx','Ly','T','profile'))
        if key in timings:return statistics.median(timings[key][-20:])+1.5
        candidates=[(abs(task['work']/ (x*y*T)-1),values,task['work']/(x*y*T)) for (x,y,T,p),values in timings.items() if p==task['profile']]
        _,values,scale=min(candidates,key=lambda r:r[0]);return statistics.median(values[-20:])*scale+1.5
    five=tasks('production_manifest_orientation5m.json','orientation');by_id={t['id']:t for _,_,ts in pools for t in ts};by_id.update({t['id']:t for t in five})
    scenarios=[]
    for factor in (1.,1.2):
        freed=[];pool_reports=[]
        for name,jobs,ts in pools:
            active_ids={active[j][0] for j in jobs if j in active}
            queue=[t for t in ts if t['id'] not in complete and t['id'] not in active_ids]
            heap=[]
            for job in jobs:
                if job in active:
                    ident,start=active[job];task=by_id[ident]
                    finish=max(now+5,start+duration(task)*factor)
                else:finish=now
                heapq.heappush(heap,(finish,job))
            releases=[]
            while heap:
                finish,job=heapq.heappop(heap)
                if queue:
                    task=queue.pop(0);heapq.heappush(heap,(finish+duration(task)*factor,job))
                else:freed.append((finish,job,name));releases.append(finish)
            pool_reports.append(dict(pool=name,n_remaining=sum(t['id'] not in complete for t in ts),releases_utc=[datetime.fromtimestamp(t,timezone.utc).isoformat() for t in sorted(releases)]))
        freed.sort();starts=[t+180 for t,_,_ in freed[:3]]
        active_id,start=active[21727284];first=max(now+5,start+duration(by_id[active_id])*factor)
        done_times={active_id:first};queue=[t for t in five if t['id'] not in complete and t['id']!=active_id]
        heap=[(first,21727284)]+[(t,job) for t,job in zip(starts,[21727285,21727286,21727287])];heapq.heapify(heap)
        while queue:
            available,job=heapq.heappop(heap);task=queue.pop(0);finish=available+duration(task)*factor;done_times[task['id']]=finish;heapq.heappush(heap,(finish,job))
        cohorts={}
        for task in five:
            key=str(task['orientation'])+' q='+task['q'];cohorts.setdefault(key,[]).append(now if task['id'] in complete else done_times[task['id']])
        deadlines={}
        for clock in ('11:00','11:15'):
            cutoff=datetime.fromisoformat('2026-11-03T'+clock+':00+00:00').timestamp()
            deadlines[clock]=dict(tasks_complete=sum((now if t['id'] in complete else done_times[t['id']])<=cutoff for t in five),cohorts_complete=sum(max(times)<=cutoff for times in cohorts.values()))
        scenarios.append(dict(duration_multiplier=factor,first_three_slot_releases=[dict(utc=datetime.fromtimestamp(t,timezone.utc).isoformat(),job=job,pool=name) for t,job,name in freed[:3]],
            pending5m_run_start_estimates=[datetime.fromtimestamp(t,timezone.utc).isoformat() for t in starts],
            all5m_finish_estimate=datetime.fromtimestamp(max(done_times.values()),timezone.utc).isoformat(),
            deadlines=deadlines,cohort_finish_estimates={key:datetime.fromtimestamp(max(times),timezone.utc).isoformat() for key,times in cohorts.items()},pools=pool_reports))
    observed=timings[(8192,8192,5000000,True)]
    result=dict(snapshot_utc=snapshot['snapshot_utc'],completed5m_durations_s=observed,
        measured5m_median_s=statistics.median(observed),n5m_complete=sum(t['id'] in complete for t in five),
        remaining5m_nominal_gpu_hours=sum(duration(t) for t in five if t['id'] not in complete)/3600,
        assumptions=['Ordinary scheduling starts the three older pending5m jobs on the next three released slots.','Three minutes of validation/build startup per newly scheduled worker.','Recent full fresh-run median durations; second scenario increases all predicted task durations20%.','No additional external QOS/group congestion; active durations are lower-bounded by five seconds from snapshot.','No deadlines, priorities, jobs or manifests are changed by this calculation.'],scenarios=scenarios)
    write_json(ROOT/'reports/v2/hpc/front5m_eta.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
