"""True-envelope front tracking, geometry diagnostics and crossing sensitivity.

Run from project root, for example:
 python3 -m campaign.v2.analysis.fronts data/runs/R2_orient/out --out campaign/v2/results/fronts
No velocity is silently removed for a low interior density. Density, holes and
the legacy connected-component estimator are reported for diagnosis instead.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
from .stats import group_rng

from .io import (ROOT, completion_status, expand_text_paths, profile_path, provenance, read_header, read_timeseries,
                 read_profile, run_id, write_csv, write_json)


def widths(t, iface, cells_per_bin, threshold=.3, smooth_bins=9, reference_time=1000.):
    """Track all bins above threshold, including disconnected portions of the band.

    The initial half-lattice strip is rolled to the middle. The envelope is the
    minimum/maximum qualifying bin, NOT a walk stopped by an interior hole.
    The legacy walk is returned separately to quantify that distinction.
    """
    d = np.roll(np.asarray(iface, float), iface.shape[1] // 4, axis=1) / cells_per_bin
    n, bins = d.shape
    if smooth_bins < 1 or smooth_bins > bins:
        raise ValueError("invalid smoothing width")
    center, half = bins // 2, max(1, int(.05 * bins))
    interior = d[:, center - half:center + half].mean(axis=1)
    iref = float(np.interp(reference_time, t, interior))
    if iref <= 0:
        raise ValueError("no active reference interior")
    kernel = np.ones(smooth_bins) / smooth_bins
    smoothed = np.array([np.convolve(row, kernel, "same") for row in d])
    envelope, component, holes = np.zeros(n), np.zeros(n), np.zeros(n)
    touches = np.zeros(n, bool)
    coordinate = np.arange(bins)
    for i, row in enumerate(smoothed):
        mask = row > threshold * iref
        hit = np.flatnonzero(mask)
        if not len(hit):
            continue
        lo, hi = int(hit[0]), int(hit[-1])
        envelope[i] = hi - lo + 1
        holes[i] = 1 - len(hit) / envelope[i]
        touches[i] = lo == 0 or hi == bins - 1
        total = d[i].sum()
        c = int(np.rint(d[i] @ coordinate / total)) if total else center
        c = min(max(c, 0), bins - 1)
        if mask[c]:
            left = right = c
            while left > 0 and mask[left - 1]:
                left -= 1
            while right < bins - 1 and mask[right + 1]:
                right += 1
            component[i] = right - left + 1
    return dict(envelope=envelope, component=component, hole_fraction=holes,
                interior=interior, interior_ratio=interior / iref,
                touches_edge=touches, iref=iref)


def velocity(t, width, bins, norm=1., t0=10000., t1=np.inf):
    t, width = np.asarray(t, float), np.asarray(width, float)
    valid = (t >= t0) & (t <= t1) & (width > .15 * bins) & (width < .85 * bins)
    bad = np.flatnonzero((t >= t0) & (t <= t1) & ~valid)
    if len(bad):
        valid &= t < t[bad[0]]
    if valid.sum() < 8:
        return dict(v=np.nan, t_start=np.nan, t_end=np.nan, npts=int(valid.sum()),
                    residual_rms=np.nan), valid
    tt = t[valid]
    yy = width[valid] / norm
    slope, intercept = np.polyfit(tt, yy, 1)
    rms = np.sqrt(np.mean((yy - intercept - slope * tt) ** 2))
    return dict(v=float(slope / 2), t_start=float(tt[0]), t_end=float(tt[-1]),
                npts=int(valid.sum()), residual_rms=float(rms)), valid


def covers_requested_window(row):
    """A finite slope after geometric truncation does not cover a fixed window."""
    return (np.isfinite(row['v']) and row['t_start']==row['t0']
            and row['t_end']==min(row['t1'],row['T']))


def reduce_run(path, thresholds=(.2, .3, .5), smooth_normal=9., tmins=(10000., 20000., 40000.), windows=()):
    meta, text = read_timeseries(path)
    a, b = meta.get("a", 1), meta.get("b", 0)
    norm = float(np.hypot(a, b))
    bins = meta["Lx"] if meta.get("profile_dir") else meta["Ly"]
    cells = meta["Lx"]
    profile, trailing = read_profile(profile_path(path), bins)
    t = np.asarray(profile["t"], float)
    status = completion_status(path, meta, text, profile, trailing)
    # Equal smoothing distance in normal lattice units for every orientation.
    smooth_bins = max(1, int(round(smooth_normal * norm)))
    if not smooth_bins % 2:
        smooth_bins += 1
    rows, records = [], []
    base = dict(run=run_id(path), p=meta.get("p", np.nan), q=meta.get("q", np.nan),
                L=meta["Lx"], Ly=meta["Ly"], ab=meta.get("ab", "1-0"), seed=meta["seed"],
                T=meta["T"], done=status["complete"], completion_errors=";".join(status["reasons"]), parameter_source=meta.get("parameter_source", "unknown"),
                smooth_normal=smooth_normal, smooth_bins=smooth_bins, trailing_bytes=trailing)
    for threshold in thresholds:
        w = widths(t, profile["iface"], cells, threshold, smooth_bins)
        for t0,t1 in [(t0,np.inf) for t0 in tmins]+list(windows):
            for estimator in ("envelope", "component"):
                fit, valid = velocity(t, w[estimator], bins, norm, t0, t1)
                rows.append(dict(**base, threshold=threshold, t0=t0, t1=t1, estimator=estimator, **fit,
                                 iref=w["iref"],
                                 interior_ratio_mean=float(w["interior_ratio"][valid].mean()) if valid.any() else np.nan,
                                 hole_fraction_mean=float(w["hole_fraction"][valid].mean()) if valid.any() else np.nan,
                                 hole_fraction_max=float(w["hole_fraction"][valid].max()) if valid.any() else np.nan,
                                 component_gap_max=float((w["envelope"] - w["component"])[valid].max() / norm)
                                 if valid.any() else np.nan))
        if abs(threshold - .3) < 1e-8:
            for i, ti in enumerate(t):
                records.append(dict(run=base["run"], time=ti, envelope=w["envelope"][i] / norm,
                                    component=w["component"][i] / norm, hole_fraction=w["hole_fraction"][i],
                                    interior_ratio=w["interior_ratio"][i], touches_edge=bool(w["touches_edge"][i])))
    return rows, records


def polynomial_crossing(x, v, se, degree=1):
    """Weighted polynomial in rescaled x; no extrapolated crossing is accepted."""
    x, v, se = np.asarray(x, float), np.asarray(v, float), np.asarray(se, float)
    if len(x) < degree + 1 or np.any(~np.isfinite(se)) or np.any(se <= 0):
        return None
    center, scale = float(x.mean()), float(np.ptp(x))
    if not scale:
        return None
    z = (x - center) / scale
    A = np.polynomial.polynomial.polyvander(z, degree)
    weighted = A / se[:, None]
    coef = np.linalg.lstsq(weighted, v / se, rcond=None)[0]
    roots = np.polynomial.polynomial.polyroots(coef)
    candidates = [center + scale * r.real for r in roots if abs(r.imag) < 1e-9]
    candidates = [r for r in candidates if x.min() <= r <= x.max()]
    if len(candidates) != 1:
        return None
    root = candidates[0]
    zr = (root - center) / scale
    derivative_z = np.polynomial.polynomial.polyval(zr, np.arange(1, len(coef)) * coef[1:])
    if derivative_z == 0:
        return None
    cov = np.linalg.inv(weighted.T @ weighted)
    jac = -scale * np.array([zr ** i for i in range(degree + 1)]) / derivative_z
    chi2 = float(np.sum(((v - A @ coef) / se) ** 2))
    dof = len(x) - degree - 1
    return dict(crossing=float(root), se_formal=float(np.sqrt(max(jac @ cov @ jac, 0))),
                slope=float(derivative_z / scale), chi2=chi2, dof=dof,
                reduced_chi2=chi2 / dof if dof > 0 else np.nan,
                center=center, scale=scale, coefficients=coef.tolist())


def crossing_models(q, values, boot=1000, rng=None, keep_samples=False):
    """Local interpolation, local line, all-point line and curvature sensitivity.

    Independent runs are resampled within each q. The original SE weights are
    held fixed in bootstrap fits, avoiding singular weights in tiny resamples.
    Statistical CIs remain conditional on the chosen model; model spread is
    reported separately and must not be conflated with those intervals.
    """
    rng = np.random.default_rng(1729) if rng is None else rng
    q = np.asarray(q, float)
    values = [np.asarray(v, float) for v in values]
    vm = np.array([v.mean() for v in values])
    se = np.array([v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan for v in values])
    brackets = np.flatnonzero(vm[:-1] * vm[1:] <= 0)
    if not len(brackets):
        return []
    k = int(brackets[np.argmin(abs((q[brackets] + q[brackets + 1]) / 2 - q.mean()))])
    specs = [("bracket", np.arange(k, k + 2), 1),
             ("local_linear", np.arange(max(0, k - 1), min(len(q), k + 3)), 1),
             ("all_linear", np.arange(len(q)), 1),
             ("all_quadratic", np.arange(len(q)), 2)]
    out = []
    for name, ix, degree in specs:
        fit = polynomial_crossing(q[ix], vm[ix], se[ix], degree)
        if fit is None:
            continue
        samples = []
        for _ in range(boot):
            if name == "bracket":
                # The sign of a nearly zero point fluctuates across resamples. Let
                # the bracketing pair change rather than clipping the confidence
                # interval to the original pair's endpoints.
                means = np.array([rng.choice(v, len(v), replace=True).mean() for v in values])
                pairs = np.flatnonzero(means[:-1] * means[1:] <= 0)
                if not len(pairs):
                    continue
                kk = int(pairs[np.argmin(abs((q[pairs]+q[pairs+1])/2-fit["crossing"]))])
                bx = np.arange(kk,kk+2)
                b = polynomial_crossing(q[bx],means[bx],se[bx],1)
            else:
                means = np.array([rng.choice(values[j], len(values[j]), replace=True).mean() for j in ix])
                b = polynomial_crossing(q[ix], means, se[ix], degree)
            if b is not None:
                samples.append(b["crossing"])
        low, high = np.quantile(samples, [.025, .975]) if len(samples) >= .8 * max(boot, 1) else (np.nan, np.nan)
        out.append(dict(model=name, degree=degree, n_q=len(ix), q_min=float(q[ix].min()),
                        q_max=float(q[ix].max()), **fit, ci_low=float(low), ci_high=float(high),
                        bootstrap_valid=len(samples), bootstrap_requested=boot,
                        fit_warning="lack_of_fit" if fit["dof"] > 0 and fit["reduced_chi2"] > 2 else ""))
        if keep_samples and name=='bracket':
            out[-1]['bootstrap_samples']=samples
    return out


def summarize(rows, boot=1000, qmin_relative=.995, qrefs=None, keep_samples=False, stable_groups=False):
    qrefs = qrefs or {.1: .03680, .5: .04596, .9: .01194}
    groups = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if not np.isfinite(r["v"]):
            continue
        key = tuple(r[k] for k in ("p", "L", "Ly", "T", "ab", "threshold", "smooth_normal", "t0", "t1", "estimator"))
        groups[key][r["q"]].append(r["v"])
    points, fits = [], []
    rng = np.random.default_rng(1729)
    names = ("p", "L", "Ly", "T", "ab", "threshold", "smooth_normal", "t0", "t1", "estimator")
    for key, vals in sorted(groups.items()):
        if stable_groups:rng=group_rng(key,1729)
        base = dict(zip(names, key))
        for q, vs in sorted(vals.items()):
            points.append(dict(**base, q=q, n=len(vs), v=float(np.mean(vs)),
                               sem=float(np.std(vs, ddof=1) / np.sqrt(len(vs))) if len(vs) > 1 else np.nan))
        lower = qmin_relative * qrefs.get(base["p"], 0)
        qs = sorted(q for q in vals if q >= lower)
        if len(qs) < 2:
            continue
        for f in crossing_models(qs, [vals[q] for q in qs], boot, rng, keep_samples):
            fits.append(dict(**base, **f))
    return points, fits


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("paths", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, default=ROOT / "campaign/v2/results/fronts")
    ap.add_argument("--thresholds", default="0.2,0.3,0.5")
    ap.add_argument("--tmins", default="10000,20000,40000")
    ap.add_argument("--windows",default="",help="Additional fixed intervals, e.g. 100000:500000,500000:1000000")
    ap.add_argument("--smooth-normal", type=float, default=9.)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--no-hashes", action="store_true")
    ap.add_argument("--include-incomplete",action="store_true",help="Diagnostic fits from unfinished runs; done flag is retained")
    args = ap.parse_args()
    paths = expand_text_paths(args.paths)
    windows=[tuple(map(float,x.split(":"))) for x in args.windows.split(",") if x]
    runs, records, inputs, skipped, unfinished = [], [], [], [], []
    for p in paths:
        prof = profile_path(p)
        if not prof.exists():
            skipped.append(str(p))
            continue
        if not args.include_incomplete and not read_header(p)["done"]:
            unfinished.append(str(p))
            continue
        r, rec = reduce_run(p, tuple(map(float, args.thresholds.split(","))),
                            args.smooth_normal, tuple(map(float, args.tmins.split(","))),windows)
        if not args.include_incomplete and not r[0]['done']:
            unfinished.append(str(p) + ': ' + r[0]['completion_errors'])
            continue
        runs.extend(r); records.extend(rec); inputs.extend([p, prof])
    if not runs:
        raise SystemExit("No profiles available: no replacement orientation result was produced.")
    points, fits = summarize(runs, args.boot)
    write_csv(args.out / "runs.csv", runs)
    write_csv(args.out / "records.csv", records)
    write_csv(args.out / "points.csv", points)
    write_json(args.out / "analysis.json", dict(settings=vars(args) | {"paths": list(map(str, args.paths)), "out": str(args.out)},
                                                fits=fits, skipped_missing_profiles=skipped, skipped_unfinished=unfinished,
                                                inputs=provenance(inputs, not args.no_hashes)))
    # Coefficients are kept as JSON arrays, never embedded in ambiguous CSV cells.
    if fits:
        write_csv(args.out / "crossings.csv", [{k: v for k, v in r.items() if k != "coefficients"} for r in fits])
    print(f"{len(inputs)//2} profiles, {len(runs)} run fits, {len(fits)} crossing fits; {args.out}")


if __name__ == "__main__":
    main()
