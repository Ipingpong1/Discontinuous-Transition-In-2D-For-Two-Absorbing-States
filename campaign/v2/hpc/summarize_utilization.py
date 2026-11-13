#!/usr/bin/env python3
"""Collect and summarize this campaign's staged-worker logs, without job actions.

Reported GPU utilization is the NVIDIA monitor's sampled quantity. It is not
GPU occupancy or allocation-wide efficiency: build, validation, inter-task
staging, and the first20s of each run are outside the poststartup statistic.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import json
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = ROOT/'reports/v2/hpc'


def records(path):
    data = path.read_bytes()
    rows = []
    lines = data.splitlines(keepends=True)
    for index, line in enumerate(lines):
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            if index == len(lines)-1 and not line.endswith(b'\n'):
                break  # A still-running worker may be in its final append.
            raise
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sync', action='store_true')
    args = parser.parse_args()
    directory = OUTPUT/'gpu_workers'
    directory.mkdir(parents=True, exist_ok=True)
    if args.sync:
        subprocess.run(['rsync','-a','--exclude=*.tmp.*',
            'torch:/home/mo2273/dn_v2_stage/gpu/workers/',str(directory)+'/'],check=True)
    rows = []
    hashes = {}
    for path in sorted(directory.glob('*/utilization.jsonl')):
        job = path.parent.name
        values = records(path)
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        samples = [r for r in values if 'gpu_util' in r]
        scientific = [r for r in samples if r.get('scientific')]
        mature = [r for r in scientific if r.get('window_mature')]
        events_path = path.with_name('events.jsonl')
        events = records(events_path) if events_path.exists() else []
        status = [r.get('status') for r in events]
        worker_end = next((r for r in reversed(events) if r.get('status')=='worker_end'),None)
        stopped = [r for r in events if 'stop' in r.get('status','') or 'failure' in r.get('status','')]
        def mean(field, selection):
            return statistics.fmean(r[field] for r in selection) if selection else None
        cumulative = next((r.get('job_poststartup_mean_util') for r in reversed(scientific)),None)
        calculated = mean('gpu_util',scientific)
        if calculated is not None and cumulative is not None and abs(calculated-cumulative)>1e-8:
            raise ValueError(f'{job}: cumulative monitor mean disagrees with raw samples')
        rows.append(dict(job=job,first_sample=samples[0]['utc'] if samples else None,
            last_sample=samples[-1]['utc'] if samples else None,
            samples=len(samples),poststartup_samples=len(scientific),mature_samples=len(mature),
            mean_util_all_logged_samples=mean('gpu_util',samples),
            mean_util_poststartup=calculated,
            minimum_mature_rolling_mean=min((r['window_mean_util'] for r in mature),default=None),
            mature_samples_below80=sum(r['window_mean_util']<80 for r in mature),
            maximum_below80_duration_s=max((r.get('low_mean_duration_s',0) for r in scientific),default=0),
            monitoring_errors=sum('monitor_error' in r for r in values),
            worker_finished=worker_end is not None,
            worker_exit_code=worker_end.get('exit_code') if worker_end else None,
            stop_or_failure_events=json.dumps(stopped,separators=(',',':'))))
    with (OUTPUT/'gpu_utilization_summary.csv').open('w',newline='') as stream:
        if rows:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    means=[r['mean_util_poststartup'] for r in rows if r['mean_util_poststartup'] is not None]
    summary=dict(utc=datetime.now(timezone.utc).isoformat(),workers=len(rows),
        finished_workers=sum(r['worker_finished'] for r in rows),
        poststartup_samples=sum(r['poststartup_samples'] for r in rows),
        minimum_worker_poststartup_mean=min(means,default=None),
        maximum_worker_poststartup_mean=max(means,default=None),
        monitoring_errors=sum(r['monitoring_errors'] for r in rows),
        workers_with_mature_below80=[r['job'] for r in rows if r['mature_samples_below80']],
        definition='Jittered NVIDIA GPU utilization samples during simulation runs; compute and in-run I/O both included. Poststartup excludes the first20s of each run. This is not allocation-wide efficiency, occupancy, or an estimate for old pre-staging attempts.',
        raw_log_sha256=hashes)
    (OUTPUT/'gpu_utilization_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='raw_log_sha256'},indent=2))


if __name__=='__main__':
    main()
