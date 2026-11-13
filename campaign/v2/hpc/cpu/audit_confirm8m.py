#!/usr/bin/env python3
"""Audit the complete, predeclared 16-shard confirmatory seed ensemble.

Mirror its JSONL worker logs into reports/v2/hpc/confirm8m_logs first. An
incomplete audit is descriptive only and exits 2; it never estimates survival
from a selected subset. The scientific reducer performs the count, physical
header, full observation schedule, and RNG-key checks.
"""
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[4]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest_path = ROOT/'campaign/v2/hpc/cpu/confirm8m/manifest.json'
    manifest = json.loads(manifest_path.read_text())
    submission = json.loads((ROOT/'campaign/v2/hpc/cpu/submitted.json').read_text())['confirmatory_8M']
    assert sha256(manifest_path) == submission['manifest_sha256']
    assert submission['primary_independent_of_original_completion']
    subprocess.run([sys.executable, '-m', 'campaign.v2.analysis.seed_tails'], cwd=ROOT, check=True)
    result = json.loads((ROOT/'campaign/v2/results/seed_tails/summary.json').read_text())
    sources = {Path(row['path']).stem: row for row in result['created_sources']}
    event_logs = []
    starts = defaultdict(list)
    finishes = defaultdict(list)
    jobs = {}
    for job in submission['jobs']:
        path = ROOT/f'reports/v2/hpc/confirm8m_logs/confirm8m_{job}.jsonl'
        events = [json.loads(line) for line in path.read_text().splitlines()]
        job_starts = [event for event in events if event['status'] == 'job_start']
        assert len(job_starts) == 1
        assert job_starts[0]['binary_sha256'] == manifest['binary_sha256']
        assert job_starts[0]['runner_sha256'] == sha256(manifest_path.with_name('runner.py'))
        for event in events:
            assert event['job'] == job
            if event['status'] == 'start':
                starts[event['task']].append(event)
            elif event['status'] == 'finished':
                finishes[event['task']].append(event)
        jobs[job] = dict(start=job_starts[0]['utc'],
                         ended=any(event['status'] == 'job_end' for event in events))
        event_logs.append(dict(path=str(path.relative_to(ROOT)), sha256=sha256(path)))
    records = []
    for task in manifest['tasks']['tails']:
        name = task['name']
        assert len(starts[name]) == 1
        assert starts[name][0]['command'] == [manifest['binary']] + task['args']
        source = sources.get(name)
        finalized = bool(source and source['completed'] and source['reached_fixed_max'])
        if finalized:
            assert source['actual_replicas'] == 500_000
            assert len(finishes[name]) == 1
            assert finishes[name][0]['returncode'] == 0 and finishes[name][0]['complete']
        records.append(dict(name=name, seed=task['seed'], job=starts[name][0]['job'],
                            finalized=finalized, source=source,
                            completion_event=finishes[name][0] if finishes[name] else None))
    completed = sum(record['finalized'] for record in records)
    complete = completed == 16
    assert complete == result['confirmatory_campaign_complete']
    if complete:
        assert all(job['ended'] for job in jobs.values())
        diagnostic = result['batch_diagnostics']['confirmatory_8M']
        assert diagnostic['replicas'] == 8_000_000
        assert diagnostic['batches'] == 16
    else:
        diagnostic = None
    audit = dict(audited_utc=datetime.now(timezone.utc).isoformat(),
                 manifest_created_utc=manifest['created_utc'],
                 manifest_sha256=sha256(manifest_path), binary_sha256=manifest['binary_sha256'],
                 primary_population='confirmatory_8M', planned_shards=16,
                 planned_replicas_per_shard=500_000, adaptive_target=0,
                 completed_shards=completed, complete=complete,
                 inference_eligible=complete,
                 selection_rule='All 16 predeclared shards required; primary independent of original16M completion or outcome.',
                 original16M='Separate supplementary ensemble; never substituted or partially pooled.',
                 validation='Declared seed and command, unchanged binary and runner hashes, success events, complete fixed counts, full physical headers and scheduled observations, integer count/footer consistency, disjoint observed RNG low words.',
                 jobs=jobs, worker_logs=event_logs, records=records,
                 scientific_diagnostic=diagnostic)
    out = ROOT/'reports/v2/hpc/confirm8m_completion_audit.json'
    out.write_text(json.dumps(audit, indent=2)+'\n')
    print(f'Confirmatory8M audit: {completed}/16 complete; inference_eligible={complete}')
    return 0 if complete else 2


if __name__ == '__main__':
    raise SystemExit(main())
