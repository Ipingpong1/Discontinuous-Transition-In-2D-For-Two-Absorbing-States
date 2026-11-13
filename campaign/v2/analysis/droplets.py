"""Finite-horizon critical-radius interpolation and conditional resampling.

The event is boundary growth by the specified horizon, not indefinite survival.
These intervals describe independent replica/binomial and planar-run sampling;
choices of horizon, radius grid, tracker and curvature law remain systematic.
"""
from __future__ import annotations
import numpy as np


def critical_radius(radius, probability):
    radius=np.asarray(radius,float);probability=np.asarray(probability,float)
    order=np.argsort(radius);radius=radius[order];probability=probability[order]
    exact=np.flatnonzero(probability==.5)
    if len(exact)==1:return float(radius[exact[0]])
    brackets=np.flatnonzero((probability[:-1]<.5)&(probability[1:]>.5))
    if len(brackets)!=1:return np.nan
    k=int(brackets[0]);p=np.clip(probability,1e-3,1-1e-3);logit=np.log(p/(1-p))
    weight=-logit[k]/(logit[k+1]-logit[k])
    return float(np.exp(np.log(radius[k])+weight*np.log(radius[k+1]/radius[k])))


def radius_bootstrap(radius, probability, replicas, boot=5000, rng=None):
    """Independent binomial parametric bootstrap conditional on sampled radii."""
    rng=rng if rng is not None else np.random.default_rng(1729)
    probability=np.asarray(probability,float);replicas=np.asarray(replicas,int)
    estimates=np.asarray([critical_radius(radius,p) for p in rng.binomial(replicas,probability,size=(boot,len(probability)))/replicas])
    valid=estimates[np.isfinite(estimates)]
    return dict(radius=critical_radius(radius,probability),
        ci_low=float(np.quantile(valid,.025)) if len(valid) else np.nan,
        ci_high=float(np.quantile(valid,.975)) if len(valid) else np.nan,
        bootstrap_requested=boot,bootstrap_valid=len(valid),samples=estimates)
