"""Read-only report of every complete fixed 5M cohort's geometry and interior.

No outcome-based filtering is applied. This supplementary audit does not change
front estimates or ensemble eligibility. Run after the serial GPU reduction.
"""
import csv
import json
from pathlib import Path
import numpy as np
from .fronts import widths
from .io import ROOT, manifest_tasks, provenance, read_profile, write_json


def main():
    folder=ROOT/'campaign/v2/results/gpu'
    status=json.loads((folder/'status.json').read_text())
    manifest=ROOT/'campaign/v2/hpc/gpu/production_manifest_orientation5m.json'
    tasks=manifest_tasks([manifest])
    complete_ids=set()
    groups=[]
    for e in status['ensembles']:
        if e['p']==.5 and e['mode']=='band' and e['T']==5000000 and e['ensemble_complete']:
            complete_ids.update(e['task_ids'])
            groups.append({k:e[k] for k in ('p','q','Lx','Ly','T','a','b','n_planned','n_complete')})
    drift=list(csv.DictReader((folder/'front5m_paired_drift.csv').open()))
    drift=[r for r in drift if r['estimator']=='envelope']
    files=[manifest,folder/'status.json',folder/'front5m_paired_drift.csv',Path(__file__),ROOT/'campaign/v2/analysis/fronts.py']
    rows=[]
    for task in tasks.values():
        if task['id'] not in complete_ids:continue
        source=ROOT/'data/v2/gpu/runs/orientation'/task['id']/'profile.bin'
        profile,trailing=read_profile(source,task['Lx'])
        if trailing:raise ValueError('Incomplete profile '+str(source))
        t=np.asarray(profile['t'],float)
        required=np.array([2e6,3.5e6,5e6])
        positions=[int(np.flatnonzero(t==point)[0]) for point in required]
        chosen=(t>=2e6)&(t<=5e6)
        norm=float(np.hypot(*task['orientation']))
        smooth=max(1,round(9*norm));smooth+=not smooth%2
        details=[]
        for threshold in (.2,.3,.5):
            w=widths(t,profile['iface'],task['Lx'],threshold,smooth)
            ww=w['envelope'][chosen]/norm
            lower=.15*task['Lx']/norm;upper=.85*task['Lx']/norm
            valid=(ww>lower)&(ww<upper)
            bad=np.flatnonzero(~valid)
            details.append(dict(threshold=threshold,width_min=float(ww.min()),width_max=float(ww.max()),
                geometry_lower=lower,geometry_upper=upper,
                all_geometry_valid=bool(np.all(valid)),
                first_invalid_time=float(t[chosen][bad[0]]) if len(bad) else None,
                max_hole_fraction=float(w['hole_fraction'][chosen].max()),
                touches_edge=bool(w['touches_edge'][chosen].any())))
        ratios=[float(w['interior_ratio'][i]) for i in positions]
        rows.append(dict(run=task['id'],q=task['q'],normal=task['orientation'],
            interior_times=required.tolist(),interior_ratios=ratios,
            interior_change_2_to5m=ratios[-1]-ratios[0],interior_change_35_to5m=ratios[-1]-ratios[1],
            interior_ratio_window_min=float(w['interior_ratio'][chosen].min()),threshold_geometry=details))
        files.append(source)
    report=dict(scope='Every predeclared 5M cohort with all four independent runs complete; no survivor selection',
        complete_groups=groups,n_complete_groups=len(groups),n_planned_groups=9,
        contrasts=drift,per_run=rows,
        interpretation='Nonzero paired velocity drift is a finite-horizon systematic; a zero-containing interval does not establish stationarity. Interior/reference offset alone is not proof of continuing interior decay.',
        mean_interval_method='All 256 ordered four-run empirical bootstrap samples; percentile confidence coverage remains approximate.',
        provenance=provenance(files))
    write_json(ROOT/'reports/v2/front5m_diagnostics.json',report)
    print(f'5M diagnostics: {len(groups)}/9 complete groups, {len(rows)} runs')

if __name__=='__main__':main()
