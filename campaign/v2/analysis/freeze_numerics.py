"""Hash a coordinated, quiescent numerical snapshot for paper packaging.

Run only after automatic maintenance stops and the final reference refresh.
This command reads scientific outputs and writes one snapshot manifest; it does
not regenerate figures, change references, submit jobs, or copy remote files.
"""
from datetime import datetime,timezone
from pathlib import Path
import json,os
from .io import ROOT,provenance,write_json


def main():
    pid_path=ROOT/'data/v2/gpu/maintenance.pid'
    if pid_path.exists():
        pid=int(pid_path.read_text())
        try:os.kill(pid,0)
        except ProcessLookupError:pass
        else:raise SystemExit('Refusing to freeze while GPU maintenance is active')
    results=ROOT/'campaign/v2/results'
    primary=json.loads((results/'gpu/status.json').read_text())
    wall=json.loads((results/'gpu_wetting/status.json').read_text())
    reference=json.loads((results/'coexistence_reference.json').read_text())
    paths=[]
    for base in [results,ROOT/'paper/v2/figs',ROOT/'campaign/v2/analysis']:
        paths.extend(p for p in base.rglob('*') if p.is_file() and p.suffix in ('.json','.csv','.txt','.pdf','.png','.py','.md') and '__pycache__' not in p.parts)
    paths.extend((ROOT/'campaign/v2/hpc/gpu').glob('*manifest*.json'))
    paths.extend((ROOT/'paper/v2').glob('*caption*.tex'))
    paths.extend([ROOT/'data/v2/gpu/bundle_imports.json',ROOT/'campaign/v2/hpc/gpu/submitted.json'])
    paths.extend((ROOT/'data/v2/gpu/repair_verifications').glob('*.json'))
    paths.extend((ROOT/'reports/v2/first5m_cohort_inputs').glob('*'))
    paths.extend((ROOT/'reports/v2').glob('*audit*.json'))
    paths.extend((ROOT/'reports/v2').glob('*spotcheck*.json'))
    paths.extend((ROOT/'reports/v2').glob('front5m_diagnostics.json'))
    complete=primary['n_complete']==primary['n_planned'] and wall['n_complete']==wall['n_planned']
    out=ROOT/'reports/v2/numerical_final_snapshot.json'
    write_json(out,dict(frozen_utc=datetime.now(timezone.utc).isoformat(),
        state='all_declared_runs_complete' if complete else 'frozen_with_unfinished_declared_runs',
        primary={k:primary[k] for k in ('n_planned','n_complete','n_ensembles','n_complete_ensembles')},
        wetting={k:wall[k] for k in ('n_planned','n_complete','n_ensembles','n_complete_ensembles')},
        incomplete_primary={k:v for k,v in primary['tasks'].items() if not v['complete']},
        incomplete_wetting={k:v for k,v in wall['tasks'].items() if not v['complete']},
        reference=reference,policy='Only complete predeclared cohorts enter publication summaries; geometry/interval eligibility is retained in per-analysis reports. This hashes the coordinated final numerical state, before paper compilation.',
        files=provenance(sorted(set(paths)))))
    print(out)
    print(f"Primary {primary['n_complete']}/{primary['n_planned']}; wetting {wall['n_complete']}/{wall['n_planned']}; {len(set(paths))} files hashed")


if __name__=='__main__':main()
