"""Renormalize retained legacy absolute crossings to the shared p=.5 reference.

This historical comparison is superseded for the paper by lazy05_audit.py.
It preserves the cached E18 estimator; it is not a fresh profile reduction.
"""
from pathlib import Path
from collections import defaultdict
import csv
import json
import re
import numpy as np
from .io import ROOT, provenance, write_json


def main():
    source=ROOT/'data/figdata/E18_lazy05_front_fronts.csv'
    reference=ROOT/'campaign/v2/results/coexistence_reference.json'
    ref=json.loads(reference.read_text());groups=defaultdict(list)
    for row in csv.DictReader(source.open()):
        match=re.fullmatch(r'p0\.5_a0\.5_r0\.25_j2_q([\d.]+)_Lx2048_Ly4096_s(\d+)',row['run'])
        if match:
            groups[float(match[1])].append((int(match[2]),float(row['v'])))
    if not groups or any({seed for seed,_ in group}!={1,2,3} or len(group)!=3 for group in groups.values()):
        raise ValueError('E18 requires all three declared seeds at every cached q')
    q=np.array(sorted(groups));v=np.array([np.mean([v for _,v in groups[x]]) for x in q])
    if not np.isfinite(v).all():raise ValueError('Nonfinite cached velocity')
    indices=np.flatnonzero((v[:-1]<0)&(v[1:]>=0))
    if len(indices)!=1:raise ValueError('E18 crossing not uniquely bracketed')
    i=int(indices[0]);crossing=float(q[i]-v[i]*(q[i+1]-q[i])/(v[i+1]-v[i]))
    result=dict(status='historical_comparison_superseded',paper_use=False,
        publication_source='campaign/v2/results/lazy05_fronts/analysis.json',
        superseded_reason='Corrected raw-profile envelope analysis with complete fixed-time windows replaces the adaptive cached tracker.',
        p=.5,update_probability=.25,qstar_lazy=crossing,
        bracket=q[i:i+2].tolist(),bracket_velocity_means=v[i:i+2].tolist(),
        n_seeds_per_q=3,reference=ref,ratio=crossing/ref['qstar'],
        legacy_reference=.04596,legacy_ratio=crossing/.04596,
        retained_legacy_quoted_ratio_uncertainty=.007*.04596/ref['qstar'],
        uncertainty_scope='The old quoted 0.007 is rescaled only; it is not a newly validated confidence interval.',
        estimator_scope='Retained E18 cached front estimator. Its negative bracket has adaptive/truncated fit windows; this script only changes the denominator.',
        provenance=provenance([source,reference,Path(__file__)]))
    write_json(ROOT/'campaign/v2/results/legacy_reference_ratios.json',result)
    print(f"HISTORICAL CACHE ONLY: legacy q*_lazy={crossing:.11f}; ratio to shared reference={result['ratio']:.8f}")


if __name__=='__main__':main()
