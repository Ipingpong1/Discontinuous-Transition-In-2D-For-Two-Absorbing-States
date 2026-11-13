"""Axis spin correlations from every final grid in complete flagship ensembles.

The Z2-symmetric ensemble has <s>=0, so raw <s_i s_j> is the ensemble
connected correlation. Also retain per-configuration mean-subtracted values
as a finite-volume diagnostic. Sites are not independent bootstrap units.
"""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import argparse
import hashlib
import json
import struct
import numpy as np
from .io import ROOT, manifest_tasks, scientific_group, write_csv, write_json, provenance
from .gpu_campaign import MANIFESTS
from .import_gpu_bundles import validate_final

DISTANCES = sorted(set(range(0,65)) | set(range(68,129,4)) |
                   set(range(136,257,8)) | set(range(288,1025,32)))

def grid(path):
    path=Path(path)
    with path.open('rb') as stream:
        header=stream.read(32)
    if len(header)!=32 or header[:8]!=b'DNSIM1\0\0':
        raise ValueError('Invalid grid header')
    Lx,Ly,t=struct.unpack('<QQQ',header[8:])
    if not Lx or Lx%64 or not Ly or path.stat().st_size!=32+Lx*Ly//8:
        raise ValueError('Invalid grid dimensions or length')
    words=np.fromfile(path,dtype='<u8',offset=32).reshape(Ly,Lx//64)
    return words,Lx,Ly,t

def correlator(words, distances):
    words=np.asarray(words,dtype=np.uint64)
    Ly,W=words.shape;Lx=64*W;N=Lx*Ly
    alive=int(np.bitwise_count(words).sum(dtype=np.uint64))
    m=2*alive/N-1
    records=[]
    for r in distances:
        if r<0 or r>min(Lx,Ly)//2: raise ValueError('Separation outside half-box')
        shift,bits=divmod(r,64)
        shifted=np.roll(words,-shift,axis=1)
        if bits:
            shifted=(shifted>>np.uint64(bits)) | (np.roll(words,-shift-1,axis=1)<<np.uint64(64-bits))
        dx=int(np.bitwise_count(words^shifted).sum(dtype=np.uint64))
        dy=int(np.bitwise_count(words^np.roll(words,-r,axis=0)).sum(dtype=np.uint64))
        cx=1-2*dx/N;cy=1-2*dy/N
        records.append(dict(r=r,different_x=dx,different_y=dy,Cx=cx,Cy=cy,
                            C=(cx+cy)/2,C_centered=(cx+cy)/2-m*m))
    return dict(alive=alive,magnetization=m,N=N),records

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/'campaign/v2/results/spatial_correlations')
    ap.add_argument('--boot',type=int,default=2000)
    args=ap.parse_args()
    tasks=manifest_tasks(MANIFESTS)
    groups=defaultdict(list)
    for task in tasks.values():
        if task['family']=='flagship_bulk' and task['mode']=='bulk' and task['Lx']>=8192:
            groups[scientific_group(task)].append(task)
    all_records=[];summaries=[];coverage=[];inputs=[]
    args.out.mkdir(parents=True,exist_ok=True)
    for key,cohort in sorted(groups.items()):
        prepared=[];reasons=[]
        for task in cohort:
            directory=ROOT/'data/v2/gpu/runs'/task['family']/task['id']
            try:
                checked=validate_final(directory,task)
                prepared.append((task,directory,checked))
            except (OSError,ValueError,KeyError) as error:
                reasons.append(dict(task=task['id'],reason=str(error)))
        coverage.append(dict(p=key[0],q=key[2],L=key[3],T=key[5],n_planned=len(cohort),
                             n_valid=len(prepared),eligible=not reasons,reasons=reasons))
        if reasons: continue
        cohort_records=[]
        for task,directory,checked in prepared:
            cache=args.out/'runs'/(task['id']+'.json')
            signature=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            payload=json.loads(cache.read_text()) if cache.exists() else {}
            if payload.get('final_sha256')!=checked['final_sha256'] or payload.get('code_sha256')!=signature:
                words,Lx,Ly,t=grid(directory/'final.bin')
                totals,records=correlator(words,DISTANCES)
                last=np.loadtxt(directory/'out.txt')[-1]
                if abs(totals['alive']/totals['N']-last[1])>6e-10:
                    raise ValueError('Final-grid alive fraction disagrees with terminal output')
                nn=next(r for r in records if r['r']==1)
                if abs((nn['different_x']+nn['different_y'])/(2*totals['N'])-last[2])>6e-10:
                    raise ValueError('Final-grid nearest-neighbor bond density disagrees with output')
                payload=dict(task=task,final_sha256=checked['final_sha256'],code_sha256=signature,
                             t_end=t,totals=totals,records=records)
                write_json(cache,payload)
            for record in payload['records']:
                all_records.append(dict(run=task['id'],p=task['p'],q=task['q'],L=task['Lx'],
                                        T=task['T'],t_end=payload['t_end'],**record))
            cohort_records.append(payload['records'])
            inputs += [directory/'final.bin',directory/'out.txt',directory/'done.json']
        rng=np.random.default_rng(int.from_bytes(hashlib.sha256(repr(key).encode()).digest()[:8],'little'))
        indices=rng.integers(0,len(cohort_records),(args.boot,len(cohort_records)))
        for i,r in enumerate(DISTANCES):
            values=np.array([rr[i]['C'] for rr in cohort_records])
            resamples=values[indices].mean(1)
            low,high=np.quantile(resamples,[.025,.975])
            summaries.append(dict(p=key[0],q=key[2],L=key[3],T=key[5],r=r,n=len(values),
                                  C=float(values.mean()),sem=float(values.std(ddof=1)/np.sqrt(len(values))),
                                  ci_low=float(low),ci_high=float(high),
                                  C_centered=float(np.mean([rr[i]['C_centered'] for rr in cohort_records]))))
    if all_records: write_csv(args.out/'correlations_runs.csv',all_records)
    else: (args.out/'correlations_runs.csv').unlink(missing_ok=True)
    if summaries: write_csv(args.out/'correlations_summary.csv',summaries)
    else: (args.out/'correlations_summary.csv').unlink(missing_ok=True)
    write_json(args.out/'analysis.json',dict(
        definition='C(r)=mean over axes and all sites of (2sigma_i-1)(2sigma_{i+r}-1); Z2 symmetry fixes ensemble mean spin to zero',
        finite_volume_diagnostic='C_centered subtracts each configuration magnetization squared; not the primary estimator',
        sampling='One terminal configuration from every run in each complete predetermined four-seed ensemble; absorbed configurations retained',
        uncertainty='Pointwise independent-run percentile bootstrap; never a bootstrap of sites',
        grid_checks='Final SHA256/dimensions/time and manifest completion; terminal alive fraction and bond density checked independently',
        distances=DISTANCES,boot=args.boot,coverage=coverage,inputs=provenance(inputs),
        manifest_inputs=provenance(MANIFESTS),
        bootstrap='One index matrix per cohort reused across all distances; cohort-specific SHA256 RNG seed',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print(f'{sum(c["eligible"] for c in coverage)}/{len(coverage)} eligible cohorts; {len(all_records)} per-run distances')

if __name__=='__main__':main()
