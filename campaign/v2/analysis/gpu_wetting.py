"""Fixed-ensemble reduction of the extended GPU wall campaign.

python3 -m campaign.v2.analysis.gpu_wetting --boot 2000
python3 -m campaign.v2.analysis.gpu_wetting --selftest

I = Ly * iface_fraction / 2 counts all interface sites per unit length of
one of the two walls. Every declared seed contributes; geometry loss and
late drift are diagnostics, never hidden activity filters. Gamma uses the
q=0 ensemble at the same (Lx,Ly), with common-baseline covariance retained.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import json
from pathlib import Path
import numpy as np

from .bulk import extend_absorbed
from .gpu_campaign import task_completion
from .io import (ROOT, complete_ensembles, manifest_tasks, profile_path,
                 provenance, read_profile, write_csv, write_json)
from .stats import mean_summary, window_mean
from .wetting import fit_models

DEFAULT_MANIFEST = ROOT/'campaign/v2/hpc/gpu/wetting_manifest.json'
DEFAULT_REFERENCE = ROOT/'campaign/v2/results/coexistence_reference.json'
DEFAULT_OUTPUT = ROOT/'campaign/v2/results/gpu_wetting'


def write_table(path, rows):
    if rows:
        write_csv(path, rows)
    else:
        path.unlink(missing_ok=True)


def adsorption_series(ts, Ly):
    """Two equal-length walls on a periodic rectangle: total sites/(2 Lx)."""
    return .5*float(Ly)*np.asarray(ts[:,4], float)


def profile_diagnostics(profile, Lx, Ly):
    """Profile moments inside each wall's initial half-period neighborhood.

    These diagnose interaction/wandering. The primary adsorption is the global
    count and never discards interface sites outside these neighborhoods.
    Width90 is a mass-quantile span, not the legacy fixed-threshold width.
    """
    d = np.asarray(profile['iface'], float)/Lx
    records = []
    for row, t in zip(d, profile['t']):
        walls = []
        for center in (0, Ly//2):
            distance = (np.arange(Ly)-center+Ly//2)%Ly-Ly//2
            selection = np.flatnonzero(abs(distance)<Ly/4)
            selection = selection[np.argsort(distance[selection])]
            y, weight = distance[selection].astype(float), row[selection]
            amount = float(weight.sum())
            if amount<=0:
                walls.append(dict(I=0.,centroid=np.nan,sigma=np.nan,width90=np.nan,margin=np.nan))
                continue
            centroid = float(y@weight/amount)
            sigma = float(np.sqrt(max(0.,((y-centroid)**2)@weight/amount)))
            positive = weight>0
            cdf = np.cumsum(weight[positive])/amount
            lo,hi = np.interp([.05,.95],cdf,y[positive])
            walls.append(dict(I=amount,centroid=centroid,sigma=sigma,
                              width90=float(hi-lo+1),margin=Ly/4-abs(centroid)-3*sigma))
        total = float(row.sum()/2)
        valid = all(np.isfinite(w['margin']) and w['margin']>0 for w in walls)
        records.append(dict(time=int(t),I_profile=total,
            I_wall0=walls[0]['I'],I_wall1=walls[1]['I'],
            centroid0=walls[0]['centroid'],centroid1=walls[1]['centroid'],
            sigma0=walls[0]['sigma'],sigma1=walls[1]['sigma'],
            width90_0=walls[0]['width90'],width90_1=walls[1]['width90'],
            minimum_margin=min(w['margin'] for w in walls),profile_geometry_valid=valid))
    return records


def reduce_wall(task, path, meta, ts):
    profile, trailing = read_profile(profile_path(path),task['Ly'])
    assert trailing==0
    profile_iface = np.asarray(profile['iface'],np.float64).sum(1)/(task['Lx']*task['Ly'])
    profile_alive = np.asarray(profile['alive'],np.float64).sum(1)/(task['Lx']*task['Ly'])
    # Normal-mode text rounds fractions to nine decimal places.
    iface_error = float(np.max(abs(profile_iface-ts[:,4])))
    alive_error = float(np.max(abs(profile_alive-ts[:,1])))
    if iface_error>6e-10 or alive_error>6e-10:
        raise ValueError('profile totals disagree with text fractions')
    geometry = profile_diagnostics(profile,task['Lx'],task['Ly'])
    series = extend_absorbed(ts,task['T']) if meta['absorbed'] else ts
    t, I = series[:,0], adsorption_series(series,task['Ly'])
    start,end = task['T']/2, float(task['T']);middle=(start+end)/2
    late = (t>=start)&(t<=end)
    early_I = window_mean(t,I,start,middle);late_I=window_mean(t,I,middle,end)
    diagnostic_late = [r for r in geometry if start<=r['time']<=end]
    alive_max_deviation = float(np.max(abs(series[late,1]-.5)))
    geometry_ok = (not meta['absorbed'] and alive_max_deviation<.2
                   and bool(diagnostic_late)
                   and all(r['profile_geometry_valid'] for r in diagnostic_late))
    margins=[r['minimum_margin'] for r in diagnostic_late if np.isfinite(r['minimum_margin'])]
    widths=[(r['width90_0']+r['width90_1'])/2 for r in diagnostic_late
            if np.isfinite(r['width90_0']) and np.isfinite(r['width90_1'])]
    row=dict(run=task['id'],path=str(path.resolve()),p=float(task['p']),alpha=float(task['alpha']),
        q=float(task['q']),Lx=task['Lx'],Ly=task['Ly'],T=task['T'],seed=task['seed'],
        t_end=meta['t_end'],absorbed=meta['absorbed'],start=start,end=end,
        I=window_mean(t,I,start,end),I_first=early_I,I_second=late_I,
        I_drift=late_I-early_I,relative_drift=(late_I-early_I)/early_I if early_I>0 else np.nan,
        alive_max_deviation=alive_max_deviation,geometry_valid=geometry_ok,
        minimum_profile_margin=min(margins) if margins else np.nan,
        width90_late_record_mean=float(np.mean(widths)) if widths else np.nan,
        profile_iface_fraction_max_error=iface_error,profile_alive_fraction_max_error=alive_error)
    for r in geometry:
        r['run']=task['id']
    return row,series,geometry


def covariance_from_sem(points):
    """Gamma(q)=I(q)-I0: baseline sampling variance is shared at fixed size."""
    C=np.diag([r['I_sem']**2 for r in points])
    for i,a in enumerate(points):
        for j,b in enumerate(points):
            if (a['Lx'],a['Ly'])==(b['Lx'],b['Ly']):
                C[i,j]+=a['I0_sem']**2
    return C


def prediction(kind, fit, delta):
    delta=np.asarray(delta,float)
    if np.any(delta<=0):raise ValueError('Wetting extrapolation requires q<qstar')
    return fit['A']*delta**(-fit['psi']) if kind=='power' else fit['A']*delta**(-1/3)+fit['B']


def prediction_jacobian(kind, fit, delta):
    delta=np.asarray(delta,float)
    if kind=='power':
        z=delta**(-fit['psi'])
        return np.column_stack([z,-fit['A']*z*np.log(delta)])
    return np.column_stack([delta**(-1/3),np.ones(len(delta))])


def coefficient_covariance(kind, fit, delta, covariance):
    J=prediction_jacobian(kind,fit,delta)
    # Conditional on measured legacy covariance. Do not inflate by residual
    # chi2: that would turn model misspecification into sampling uncertainty.
    return np.linalg.inv(J.T@np.linalg.solve(covariance,J))


def legacy_bootstrap(points, covariance, qstar, path, boot):
    rows=list(csv.DictReader(path.open()))
    groups=defaultdict(list);allgroups=defaultdict(list)
    for row in rows:
        key=(float(row['p']),float(row['q']),int(row['Lx']),int(row['Ly']))
        allgroups[key].append(row)
        if row['eligible']=='True':groups[key].append(row)
    baselines={};values=[];used=[];excluded=0
    for point in points:
        key=(point['p'],point['q'],point['Lx'],point['Ly'])
        rs=groups[key];base_key=(point['p'],0.,1024,point['baseline_Ly'])
        base=groups[base_key]
        a=np.array([float(r['I']) for r in rs]);z=np.array([float(r['I']) for r in base])
        if len(a)!=point['n'] or not np.isclose(a.mean()-z.mean(),point['Gamma'],rtol=1e-10,atol=1e-8):
            raise ValueError('Legacy run table does not reproduce frozen fit points')
        if not np.isclose(a.std(ddof=1)/np.sqrt(len(a)),point['I_sem'],rtol=1e-8):
            raise ValueError('Legacy sampling variance does not reproduce frozen fit points')
        baselines[base_key]=z;values.append((a,base_key));used.extend(rs)
        excluded+=len(allgroups[key])-len(rs)
    rng=np.random.default_rng(1729)
    y=np.empty((boot,len(points)))
    if boot:
        baseline_draws={key:z[rng.integers(0,len(z),(boot,len(z)))].mean(1)
                        for key,z in sorted(baselines.items())}
        for j,(a,key) in enumerate(values):
            y[:,j]=a[rng.integers(0,len(a),(boot,len(a)))].mean(1)-baseline_draws[key]
    delta=1-np.array([r['q'] for r in points])/qstar
    params={kind:[] for kind in ('power','third_offset')}
    for draw in y:
        fit=fit_models(delta,draw,covariance)
        params['power'].append([fit['power']['A'],fit['power']['psi']])
        params['third_offset'].append([fit['third_offset']['A'],fit['third_offset']['B']])
    return {key:np.asarray(value).reshape(boot,2) for key,value in params.items()},dict(
        sampling_unit='Independent legacy run; one resampled q=0 mean shared by all points using that baseline',
        bootstrap_seed=1729,bootstrap=boot,n_fit_points=len(points),
        n_fit_runs=len(used),n_geometry_truncated_fit_runs=sum(r['truncated']=='True' for r in used),
        n_ineligible_runs_in_fit_groups=excluded,
        n_points_with_fallback_transverse_baseline=sum(r['baseline_Ly']!=r['Ly'] for r in points),
        caveat='Legacy eligibility/windows are frozen, including geometry truncation and fallback baselines. New wall trajectories are unfiltered and use matched-size baselines.')


def legacy_comparison(qstar, reference_range, legacy_path, boot=2000):
    data=json.loads(legacy_path.read_text())['result']
    points=data['fit_points'];q=np.array([r['q'] for r in points])
    gamma=np.array([r['Gamma'] for r in points]);C=np.array(data['gamma_covariance'])
    if np.any(q>=qstar):
        return dict(status='legacy_fit_points_not_below_reference')
    delta=1-q/qstar;fit=fit_models(delta,gamma,C)
    cov={kind:coefficient_covariance(kind,fit[kind],delta,C) for kind in fit}
    params,bootstrap_settings=legacy_bootstrap(points,C,qstar,legacy_path.with_name('runs.csv'),boot)
    for kind,draws in params.items():
        if len(draws):
            low,high=np.quantile(draws,[.025,.975],axis=0)
            for j,name in enumerate(('A','psi' if kind=='power' else 'B')):
                fit[kind][name+'_ci_low']=low[j];fit[kind][name+'_ci_high']=high[j]
    sensitivity=[]
    for value in np.unique(np.r_[np.linspace(qstar-reference_range,qstar+reference_range,11),qstar]):
        if np.all(q<value):
            dd=1-q/value;ff=fit_models(dd,gamma,C)
            cc={kind:coefficient_covariance(kind,ff[kind],dd,C).tolist() for kind in ff}
            sensitivity.append(dict(qstar=float(value),**ff['third_offset'],models=ff,
                                    coefficient_covariance=cc))
    return dict(status='ok',qstar=qstar,source_qstar=data['qstar'],central=fit,
        third_coefficient_covariance=cov['third_offset'].tolist(),
        coefficient_covariance={key:value.tolist() for key,value in cov.items()},
        parameter_bootstrap={key:value.tolist() for key,value in params.items()},
        bootstrap_settings=bootstrap_settings,qstar_sensitivity=sensitivity,fit_points=points,
        fit_delta_min=float(np.min(delta)),fit_delta_max=float(np.max(delta)),
        interpretation='Both two-parameter models are fitted only to the same seven legacy points. New observations are held out; no new-data refit, exponent estimate, or automatic asymptotic interpretation.')


def bootstrap_predictions(kind, parameters, delta):
    parameters=np.asarray(parameters,float);delta=np.asarray(delta,float)
    if not len(parameters):return np.empty((0,len(delta)))
    if kind=='power':return parameters[:,0,None]*delta[None,:]**(-parameters[:,1,None])
    return parameters[:,0,None]*delta[None,:]**(-1/3)+parameters[:,1,None]


def quadratic_residual(observed, predicted, covariance):
    residual=np.asarray(observed)-predicted
    rank=int(np.linalg.matrix_rank(covariance))
    return dict(quadratic=float(residual@np.linalg.solve(covariance,residual)) if rank==len(residual) else np.nan,
                n_observations=len(residual),covariance_rank=rank,
                fitted_parameters_from_new_data=0)


def heldout_predictions(points, covariance, model, gamma_draws, size_comparisons):
    """Separate size-specific prediction tests, never a refit to new data."""
    if not points or model['status']!='ok':return [],dict(status='awaiting_complete_wall_and_baseline_ensembles')
    records=[];diagnostics={};q=np.array([r['q'] for r in points]);y=np.array([r['Gamma'] for r in points])
    qs=model['qstar'];central_delta=1-q/qs
    for kind in ('power','third_offset'):
        fit=model['central'][kind];valid=central_delta>1e-12
        ids=np.flatnonzero(valid);delta=central_delta[valid]
        pred=prediction(kind,fit,delta)
        J=prediction_jacobian(kind,fit,delta)
        Cpred=J@np.array(model['coefficient_covariance'][kind])@J.T
        boots=bootstrap_predictions(kind,model['parameter_bootstrap'][kind],delta)
        low,high=np.quantile(boots,[.025,.975],axis=0) if len(boots) else (np.full(len(ids),np.nan),)*2
        for j,index in enumerate(ids):
            point=points[index];key=(point['q'],point['Lx'],point['Ly'])
            residual=point['Gamma']-pred[j]
            row=dict(model=kind,q=point['q'],Lx=point['Lx'],Ly=point['Ly'],delta=point['delta'],
                observed=point['Gamma'],prediction=float(pred[j]),prediction_sampling_se=float(np.sqrt(Cpred[j,j])),
                prediction_ci_low=float(low[j]),prediction_ci_high=float(high[j]),
                residual=float(residual),residual_ci_low=np.nan,residual_ci_high=np.nan,
                qstar_sensitivity_min=np.nan,qstar_sensitivity_max=np.nan,
                sensitivity_domain_complete=True,geometry_valid=point['geometry_valid'],
                wall_drift_ci_excludes_zero=point['I_drift_ci_low']>0 or point['I_drift_ci_high']<0,
                baseline_drift_ci_excludes_zero=point['baseline_drift_ci_low']>0 or point['baseline_drift_ci_high']<0)
            if len(boots) and key in gamma_draws:
                rboot=gamma_draws[key]-boots[:,j]
                row['residual_ci_low'],row['residual_ci_high']=np.quantile(rboot,[.025,.975])
            sensitive=[prediction(kind,c['models'][kind],[1-point['q']/c['qstar']])[0]
                       for c in model['qstar_sensitivity'] if point['q']<c['qstar']-1e-12]
            row['sensitivity_domain_complete']=len(sensitive)==len(model['qstar_sensitivity'])
            if row['sensitivity_domain_complete']:
                row['qstar_sensitivity_min'],row['qstar_sensitivity_max']=min(sensitive),max(sensitive)
            records.append(row)
        tests=[]
        for Lx in sorted(set(r['Lx'] for r in points)):
            selection=np.array([i for i,r in enumerate(points) if r['Lx']==Lx])
            for c in model['qstar_sensitivity']:
                dd=1-q[selection]/c['qstar']
                item=dict(Lx=Lx,Ly=points[selection[0]]['Ly'],qstar=c['qstar'],
                          q_values=q[selection].tolist(),domain_valid=bool(np.all(dd>1e-12)))
                if item['domain_valid']:
                    pp=prediction(kind,c['models'][kind],dd)
                    jj=prediction_jacobian(kind,c['models'][kind],dd)
                    cp=jj@np.array(c['coefficient_covariance'][kind])@jj.T
                    total=covariance[np.ix_(selection,selection)]+cp
                    item.update(quadratic_residual(y[selection],pp,total),
                                residual_covariance=total.tolist(),prediction_covariance=cp.tolist())
                tests.append(item)
        diagnostics[kind]=dict(training_chi2=fit['chi2'],training_dof=fit['dof'],
            prediction_covariance=Cpred.tolist(),prediction_covariance_order=ids.tolist(),
            size_specific_tests=tests,n_new_points_outside_central_domain=int((~valid).sum()))
    gate=dict(n_complete_new_q_size_points=len(points),expected_points=6,
        all_fixed_ensembles_available=len(points)==6,
        all_geometry_valid=all(r['geometry_valid'] for r in points),
        drift_intervals_available=all(np.isfinite(r[k]) for r in points for k in
                                     ('I_drift_ci_low','I_drift_ci_high','baseline_drift_ci_low','baseline_drift_ci_high')),
        wall_drift_flagged=[dict(q=r['q'],Lx=r['Lx']) for r in points
                           if r['I_drift_ci_low']>0 or r['I_drift_ci_high']<0],
        baseline_drift_flagged=[dict(Lx=r['Lx']) for r in points
                               if r['baseline_drift_ci_low']>0 or r['baseline_drift_ci_high']<0],
        n_size_comparisons=len(size_comparisons),size_comparisons=size_comparisons,
        asymptotic_exponent_inference_permitted=False)
    return records,dict(status='descriptive_heldout_prediction_diagnostics',models=diagnostics,
        interpretation_checks=gate,
        uncertainty='Conditional parameter covariance and independent legacy-run bootstrap; independent new-run bootstrap with shared q=0 baseline. qstar sensitivity is separate, not a random distribution.',
        quadratic_reference='Covariance-weighted residual quadratic with no parameters fitted to new data, evaluated separately at each Lx. Approximate Gaussian sampling diagnostic, not a calibrated p-value with two/four seeds.',
        limitations='Interpret only after reviewing all geometry, drift and size checks. Failure to detect drift or size dependence does not prove stationarity/equivalence. Legacy conditioning and model misspecification remain outside sampling intervals.')


def summarize_walls(rows, series_by_id, boot=2000):
    """All rows supplied here belong to complete predetermined ensembles."""
    groups=defaultdict(list)
    for r in rows:groups[(r['q'],r['Lx'],r['Ly'],r['T'])].append(r)
    rng=np.random.default_rng(20261103)
    summaries=[];curves=[];run_samples={}
    for key,rs in sorted(groups.items()):
        q,Lx,Ly,T=key
        base=dict(q=q,Lx=Lx,Ly=Ly,T=T,n=len(rs),start=T/2,end=T,
                  n_absorbed=sum(r['absorbed'] for r in rs),
                  geometry_valid=all(r['geometry_valid'] for r in rs),
                  n_geometry_invalid=sum(not r['geometry_valid'] for r in rs),
                  maximum_alive_deviation=max(r['alive_max_deviation'] for r in rs))
        index=rng.integers(0,len(rs),(boot,len(rs))) if boot else None
        samples={}
        for name in ('I','I_drift','relative_drift'):
            values=np.array([r[name] for r in rs])
            summary=mean_summary(values,boot=0)
            if boot:
                draws=values[index].mean(1);samples[name]=draws
                summary['ci_low'],summary['ci_high']=np.quantile(draws,[.025,.975])
            base.update({name+'_'+k:v for k,v in summary.items()})
        base['drift_ci_excludes_zero']=bool(base['I_drift_ci_low']>0 or base['I_drift_ci_high']<0)
        summaries.append(base);run_samples[(q,Lx,Ly)]=samples
        grid=np.unique(np.concatenate([series_by_id[r['run']][:,0] for r in rs]))
        values=np.array([np.interp(grid,series_by_id[r['run']][:,0],
                                   adsorption_series(series_by_id[r['run']],Ly)) for r in rs])
        means=values.mean(0)
        if boot:
            draws=values[index].mean(1);lo,hi=np.quantile(draws,[.025,.975],axis=0)
        else:lo=hi=np.full(len(grid),np.nan)
        for t,mean,low,high in zip(grid,means,lo,hi):
            curves.append(dict(q=q,Lx=Lx,Ly=Ly,T=T,n=len(rs),time=float(t),
                               I=float(mean),ci_low=float(low),ci_high=float(high)))
    baseline={(r['Lx'],r['Ly']):r for r in summaries if r['q']==0}
    points=[];gamma_draws={}
    for r in summaries:
        if r['q']==0 or (r['Lx'],r['Ly']) not in baseline:continue
        b=baseline[(r['Lx'],r['Ly'])]
        point=dict(q=r['q'],Lx=r['Lx'],Ly=r['Ly'],T=r['T'],n=r['n'],
            I=r['I_mean'],I_sem=r['I_sem'],I_ci_low=r['I_ci_low'],I_ci_high=r['I_ci_high'],
            I0=b['I_mean'],I0_sem=b['I_sem'],baseline_n=b['n'],baseline_T=b['T'],
            Gamma=r['I_mean']-b['I_mean'],Gamma_sem=float(np.hypot(r['I_sem'],b['I_sem'])),
            Gamma_ci_low=np.nan,Gamma_ci_high=np.nan,
            geometry_valid=bool(r['geometry_valid'] and b['geometry_valid']),
            n_geometry_invalid=r['n_geometry_invalid'],baseline_geometry_valid=b['geometry_valid'],
            I_drift=r['I_drift_mean'],I_drift_ci_low=r['I_drift_ci_low'],I_drift_ci_high=r['I_drift_ci_high'],
            relative_drift=r['relative_drift_mean'],relative_drift_ci_low=r['relative_drift_ci_low'],
            relative_drift_ci_high=r['relative_drift_ci_high'],baseline_drift=b['I_drift_mean'],baseline_drift_ci_low=b['I_drift_ci_low'],
            baseline_drift_ci_high=b['I_drift_ci_high'],
            Gamma_first=r['I_mean']-r['I_drift_mean']/2-b['I_mean']+b['I_drift_mean']/2,
            Gamma_second=r['I_mean']+r['I_drift_mean']/2-b['I_mean']-b['I_drift_mean']/2,
            maximum_alive_deviation=r['maximum_alive_deviation'])
        if boot:
            draws=run_samples[(r['q'],r['Lx'],r['Ly'])]['I']-run_samples[(0.,r['Lx'],r['Ly'])]['I']
            point['Gamma_ci_low'],point['Gamma_ci_high']=np.quantile(draws,[.025,.975])
            gamma_draws[(r['q'],r['Lx'],r['Ly'])]=draws
        points.append(point)
    comparisons=[]
    qs=sorted(set(p['q'] for p in points));lookup={(r['q'],r['Lx'],r['Ly']):r for r in points}
    for q in qs:
        a,b=lookup.get((q,4096,8192)),lookup.get((q,16384,8192))
        if a is None or b is None:continue
        item=dict(q=q,Ly=8192,small_Lx=4096,large_Lx=16384,
            difference=b['Gamma']-a['Gamma'],ratio=b['Gamma']/a['Gamma'] if a['Gamma']!=0 else np.nan,
            difference_ci_low=np.nan,difference_ci_high=np.nan,ratio_ci_low=np.nan,ratio_ci_high=np.nan,
            geometry_valid=bool(a['geometry_valid'] and b['geometry_valid']))
        if boot:
            aa,bb=gamma_draws[(q,4096,8192)],gamma_draws[(q,16384,8192)]
            item['difference_ci_low'],item['difference_ci_high']=np.quantile(bb-aa,[.025,.975])
            good=aa>0
            if good.all():item['ratio_ci_low'],item['ratio_ci_high']=np.quantile(bb/aa,[.025,.975])
        comparisons.append(item)
    return summaries,curves,points,comparisons,gamma_draws


def analyze(root,manifest,reference,out,legacy,boot=2000,hashes=True):
    tasks=manifest_tasks([manifest]);status={};loaded={};inputs=[manifest,reference,legacy]
    for ident,task in tasks.items():
        assert task['mode']=='wall' and task['family']=='wetting'
        path=root/task['family']/ident/'out.txt'
        status[ident],meta,ts=task_completion(task,path)
        if status[ident]['complete']:
            try:
                loaded[ident]=reduce_wall(task,path,meta,ts)
            except (ValueError,AssertionError,OSError) as exc:
                status[ident]['complete']=False;status[ident]['reasons'].append(str(exc))
    ensembles=complete_ensembles(tasks,status)
    eligible={ident for e in ensembles if e['ensemble_complete'] for ident in e['task_ids']}
    rows=[];series={};geometry=[]
    for ident,(row,ts,diagnostics) in loaded.items():
        row['ensemble_complete']=ident in eligible
        rows.append(row);series[ident]=ts;geometry.extend(diagnostics)
        inputs += [Path(row['path']),profile_path(row['path']),Path(row['path']).parent/'done.json']
    final_rows=[r for r in rows if r['ensemble_complete']]
    summary,curves,points,size_comparison,gamma_draws=summarize_walls(final_rows,series,boot)
    ref=json.loads(reference.read_text());qstar=float(ref['qstar']);dq=float(ref['systematic_half_range'])
    model=legacy_comparison(qstar,dq,legacy,boot)
    inputs.append(legacy.with_name('runs.csv'))
    for p in points:
        p.update(delta=1-p['q']/qstar,delta_low=1-p['q']/(qstar-dq),delta_high=1-p['q']/(qstar+dq),
                 qstar_range_crosses_q=p['q']>=qstar-dq-1e-12,
                 legacy_third_prediction=np.nan,legacy_prediction_sampling_se=np.nan,
                 legacy_prediction_sensitivity_min=np.nan,legacy_prediction_sensitivity_max=np.nan,
                 residual_from_legacy_fit=np.nan)
        if model['status']=='ok' and p['delta']>0:
            c=model['central']['third_offset'];z=np.array([p['delta']**(-1/3),1.])
            prediction=c['A']*z[0]+c['B']
            se=np.sqrt(max(0.,z@np.array(model['third_coefficient_covariance'])@z))
            predictions=[c['A']*(1-p['q']/c['qstar'])**(-1/3)+c['B']
                         for c in model['qstar_sensitivity'] if p['q']<c['qstar']-1e-12]
            p.update(legacy_third_prediction=prediction,legacy_prediction_sampling_se=float(se),
                legacy_prediction_sensitivity_min=min(predictions) if len(predictions)==len(model['qstar_sensitivity']) else np.nan,
                legacy_prediction_sensitivity_max=max(predictions) if len(predictions)==len(model['qstar_sensitivity']) else np.nan,
                residual_from_legacy_fit=p['Gamma']-prediction)
    C=covariance_from_sem(points) if points else np.zeros((0,0))
    predictions,prediction_tests=heldout_predictions(points,C,model,gamma_draws,size_comparison)
    out.mkdir(parents=True,exist_ok=True)
    write_json(out/'status.json',dict(tasks=status,ensembles=ensembles,n_planned=len(tasks),
        n_complete=len(loaded),n_complete_ensembles=sum(e['ensemble_complete'] for e in ensembles),
        n_ensembles=len(ensembles),n_adsorption_points=len(points)))
    write_json(out/'analysis.json',dict(reference=ref,settings=dict(bootstrap=boot,
        sampling_unit='independent trajectory; all predetermined seeds required',
        adsorption='I=Ly*iface_fraction/2, all interface sites; two periodic walls',
        averaging='trapezoidal time integration on [T/2,T], then equal run weights',
        geometry='No filtering of primary I; flag absorption, max|alive-.5|>=.2, or a wall-profile mean+3sigma reaching its half-period boundary',
        baseline='Independent q=0 ensemble at the same Lx,Ly; [1e6,2e6] for baseline versus [5e6,1e7] for nonzero q',
        covariance='Independent nonzero-q variance plus a common baseline variance within each size',
        intervals='95% independent-run percentile bootstrap, approximate with 2 baseline or 4 wall seeds',
        drift='Difference of the second and first halves of the declared late window; reported rather than selected away',
        qstar_uncertainty='Separate supplied systematic sensitivity, not a probability distribution',
        inference='Three new q values test extrapolation and finite wall length; no new asymptotic exponent fit'),
        legacy_comparison=model,gamma_covariance=C.tolist(),
        gamma_covariance_order=[dict(q=p['q'],Lx=p['Lx'],Ly=p['Ly']) for p in points],
        size_comparisons=size_comparison,heldout_predictions=prediction_tests))
    write_json(out/'provenance.json',provenance(inputs,hashes))
    for name,data in [('runs',rows),('summary',summary),('curves',curves),('points',points),
                      ('size_comparison',size_comparison),('profile_geometry',geometry),('heldout_predictions',predictions)]:
        write_table(out/(name+'.csv'),data)
    print(f"{len(loaded)}/{len(tasks)} wall runs complete; {sum(e['ensemble_complete'] for e in ensembles)}/{len(ensembles)} ensembles; {len(points)} Gamma points")
    return dict(status=status,points=points,size_comparisons=size_comparison)


def selftest():
    ts=np.array([[0,.5,0,0,.1],[5,.5,0,0,.2],[10,.5,0,0,.3]],float)
    I=adsorption_series(ts,256)
    assert np.allclose(I,[12.8,25.6,38.4])
    assert abs(window_mean(ts[:,0],I,5,10)-32)<1e-12
    dtype=np.dtype([('t','<i8'),('alive','<u4',(64,)),('iface','<u4',(64,))])
    profile=np.zeros(3,dtype=dtype);profile['t']=[0,5,10]
    profile['iface'][:,0]=64;profile['iface'][:,32]=64
    diagnostics=profile_diagnostics(profile,64,64)
    assert all(abs(r['I_profile']-1)<1e-12 and r['profile_geometry_valid'] for r in diagnostics)
    points=[dict(Lx=1,Ly=2,I_sem=2.,I0_sem=1.),dict(Lx=1,Ly=2,I_sem=3.,I0_sem=1.),
            dict(Lx=4,Ly=2,I_sem=2.,I0_sem=.5)]
    C=covariance_from_sem(points)
    assert np.allclose(C,[[5,1,0],[1,10,0],[0,0,4.25]])
    # Completion depends on all planned seeds, even when the currently available
    # runs are all complete. This prevents early-extinction completion bias.
    tasks={str(i):dict(id=str(i),family='wetting',mode='wall',p='.5',alpha='.5',q='.0459',
                      Lx=64,Ly=64,T=10,orientation=None) for i in range(4)}
    state={str(i):dict(complete=True) for i in range(3)}
    assert not complete_ensembles(tasks,state)[0]['ensemble_complete']
    # Analytic parameter propagation agrees with finite differences for the
    # nonlinear power fit; a coherent residual must use its full covariance.
    delta=np.array([.002,.01,.04]);fit=dict(A=2.,psi=.45)
    J=prediction_jacobian('power',fit,delta)
    for j,name in enumerate(('A','psi')):
        lo=dict(fit);hi=dict(fit);lo[name]-=1e-6;hi[name]+=1e-6
        derivative=(prediction('power',hi,delta)-prediction('power',lo,delta))/2e-6
        assert np.allclose(J[:,j],derivative,rtol=1e-8)
    statistic=quadratic_residual([1.,1.],np.zeros(2),np.array([[2.,1.],[1.,2.]]))
    assert abs(statistic['quadratic']-2/3)<1e-12
    # A held-out point equal to its frozen legacy prediction has zero residual;
    # evaluating it must not alter the fitted legacy model.
    model=legacy_comparison(.04592,.00002,ROOT/'campaign/v2/results/wetting/analysis.json',boot=12)
    q=.04584;delta=1-q/.04592;value=float(prediction('third_offset',model['central']['third_offset'],[delta])[0])
    point=dict(q=q,Lx=4096,Ly=8192,delta=delta,Gamma=value,geometry_valid=True,
               I_drift_ci_low=-1.,I_drift_ci_high=1.,baseline_drift_ci_low=-1.,baseline_drift_ci_high=1.)
    frozen=json.dumps(model,default=lambda x:float(x))
    rr,diagnostics=heldout_predictions([point],np.eye(1),model,{(q,4096,8192):np.full(12,value)},[])
    assert abs(next(r for r in rr if r['model']=='third_offset')['residual'])<1e-12
    assert json.dumps(model,default=lambda x:float(x))==frozen
    assert not diagnostics['interpretation_checks']['asymptotic_exponent_inference_permitted']
    print('gpu_wetting selftest PASS: normalization, geometry, shared covariance, fixed ensembles, predictive covariance and no new-data fit')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',type=Path,default=DEFAULT_MANIFEST)
    ap.add_argument('--reference',type=Path,default=DEFAULT_REFERENCE)
    ap.add_argument('--root',type=Path,default=ROOT/'data/v2/gpu/runs')
    ap.add_argument('--out',type=Path,default=DEFAULT_OUTPUT)
    ap.add_argument('--legacy',type=Path,default=ROOT/'campaign/v2/results/wetting/analysis.json')
    ap.add_argument('--boot',type=int,default=2000)
    ap.add_argument('--no-hashes',action='store_true')
    ap.add_argument('--selftest',action='store_true')
    args=ap.parse_args()
    if args.selftest:return selftest()
    analyze(args.root,args.manifest,args.reference,args.out,args.legacy,args.boot,not args.no_hashes)

if __name__=='__main__':main()
