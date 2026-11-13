"""Small statistical routines with explicit sampling units and censoring."""
from __future__ import annotations

from statistics import NormalDist
from math import comb
import numpy as np
import hashlib

_FOUR_RUN_RESAMPLES = np.indices((4,4,4,4)).reshape(4,-1).T


def window_mean(t, y, start, end):
    """Trapezoidal time mean, interpolating endpoints; never silently extrapolate."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    if end <= start or start < t[0] or end > t[-1]:
        return np.nan
    inside = (t > start) & (t < end)
    tt = np.r_[start, t[inside], end]
    yy = np.r_[np.interp(start, t, y), y[inside], np.interp(end, t, y)]
    return float(np.sum(np.diff(tt) * (yy[1:] + yy[:-1]) / 2) / (end - start))


def mean_summary(values, boot=2000, rng=None):
    """Independent-run percentile intervals; time records are not replicates.

    For n=4, enumerate all256 ordered empirical-bootstrap samples. This removes
    Monte Carlo error; percentile confidence coverage remains approximate.
    Other n retain the requested Monte Carlo resamples. boot=0 disables CIs.
    """
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return dict(n=0, mean=np.nan, sem=np.nan, ci_low=np.nan, ci_high=np.nan)
    result = dict(n=len(x), mean=float(x.mean()), sem=np.nan, ci_low=np.nan, ci_high=np.nan)
    if len(x) > 1:
        result["sem"] = float(x.std(ddof=1) / np.sqrt(len(x)))
        if boot:
            rng = np.random.default_rng(1) if rng is None else rng
            sampled_indices = rng.integers(0, len(x), (boot, len(x)))
            # Retain the previous RNG advancement for callers sharing a stream
            # with KM or other estimators; only mean-interval calculation changes.
            indices = _FOUR_RUN_RESAMPLES if len(x)==4 else sampled_indices
            m = x[indices].mean(axis=1)
            result["ci_low"], result["ci_high"] = np.quantile(m, [.025, .975])
    return result


def kaplan_meier(times, events, confidence=.95):
    """Product-limit estimate; deaths precede censors at tied times.

    Greenwood log-log pointwise bands are used away from 0 and 1. If no event has
    occurred the lower band is the exact all-survivor binomial bound for the
    subjects observed through that time, avoiding a misleading zero-width band.
    Event times for half-density are first recorded crossings, not exact crossings.
    """
    times, events = np.asarray(times, float), np.asarray(events, bool)
    if times.shape != events.shape or not len(times) or np.any(times < 0):
        raise ValueError("nonempty matching nonnegative durations/events required")
    if not np.all(np.isfinite(times)):
        raise ValueError("censored durations must be finite observation endpoints")
    alpha = 1 - confidence if confidence is not None else .05
    z = NormalDist().inv_cdf(1 - alpha / 2) if confidence is not None else 0.
    survival, greenwood = 1., 0.
    rows = [dict(time=0., at_risk=len(times), events=0, censored=0,
                 survival=1., low=1., high=1.)]
    for t in np.unique(times):
        risk = int(np.sum(times >= t))
        d = int(np.sum((times == t) & events))
        c = int(np.sum((times == t) & ~events))
        survival *= 1 - d / risk
        if d and risk > d:
            greenwood += d / (risk * (risk - d))
        if confidence is None:
            low = high = np.nan
        elif survival == 0:
            low = 0.
            # A finite sample with no remaining survivors does not establish
            # exactly zero survival probability. If no earlier censoring exists,
            # use the exact zero-success binomial bound. Otherwise retain a
            # conservative upper bound rather than invent a zero-width interval.
            earlier_censor = np.any((times < t) & ~events)
            high = 1. if earlier_censor else 1-(alpha/2)**(1/len(times))
        elif survival == 1:
            if risk == len(times):
                low=(alpha/2)**(1/risk)
            else:
                # Conservative exact lower bound: earlier censored subjects
                # count as failures; all risk-set subjects survived through t.
                left,right=0.,1.
                for _ in range(50):
                    p=(left+right)/2
                    tail=sum(comb(len(times),j)*p**j*(1-p)**(len(times)-j)
                             for j in range(risk,len(times)+1))
                    if tail < alpha/2:left=p
                    else:right=p
                low=(left+right)/2
            high=1.
        else:
            sd = np.sqrt(greenwood) / abs(np.log(survival))
            ll = np.log(-np.log(survival))
            low = np.exp(-np.exp(ll + z * sd))
            high = np.exp(-np.exp(ll - z * sd))
        rows.append(dict(time=float(t), at_risk=risk, events=d, censored=c,
                         survival=float(survival), low=float(low), high=float(high)))
    return rows


def km_median(curve):
    hits = [r["time"] for r in curve if r["survival"] <= .5]
    return float(hits[0]) if hits else np.nan


def restricted_mean(curve, horizon):
    """Area under the KM step function, only up to the supplied observed horizon."""
    if horizon < 0 or horizon > curve[-1]["time"]:
        raise ValueError("RMST horizon exceeds observations")
    total, left, surv = 0., 0., 1.
    for r in curve[1:]:
        right = min(r["time"], horizon)
        total += (right - left) * surv
        left, surv = right, r["survival"]
        if left >= horizon:
            break
    return float(total)


def survival_summary(times, events, horizon=None, boot=2000, rng=None):
    times, events = np.asarray(times, float), np.asarray(events, bool)
    horizon = float(times.max()) if horizon is None else float(horizon)
    curve = kaplan_meier(times, events)
    result = dict(n=len(times), events=int(events.sum()), censored=int((~events).sum()),
                  horizon=horizon, median=km_median(curve),
                  rmst=restricted_mean(curve, horizon), rmst_low=np.nan, rmst_high=np.nan,
                  bootstrap_valid=0, curve=curve)
    if boot and len(times) > 1 and events.any():
        rng = np.random.default_rng(1) if rng is None else rng
        samples = []
        for _ in range(boot):
            ix = rng.integers(0, len(times), len(times))
            # All-early-event resamples have S=0 afterward; all-early-censor samples
            # cannot identify the requested horizon and are omitted, explicitly counted.
            c = kaplan_meier(times[ix], events[ix], confidence=None)
            if c[-1]["time"] < horizon:
                if c[-1]["survival"] == 0:
                    c.append(dict(time=horizon, survival=0.))
                else:
                    continue
            samples.append(restricted_mean(c, horizon))
        if samples:
            result["rmst_low"], result["rmst_high"] = np.quantile(samples, [.025, .975])
        result["bootstrap_valid"] = len(samples)
    return result


def group_rng(key, seed=1729):
    """Stable independent bootstrap stream for a scientific group/window."""
    digest=hashlib.sha256(repr((seed,key)).encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8],'little'))
