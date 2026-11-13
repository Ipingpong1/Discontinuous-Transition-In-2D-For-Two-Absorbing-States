"""Unconditioned bulk averages and separate censored half-density/absorption times.

Run: python3 -m campaign.v2.analysis.bulk --out campaign/v2/results/bulk
The main density estimate includes every completed run and every time in its
predeclared late window. A run absorbed early contributes zero thereafter.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np

from .io import ROOT, expand_text_paths, provenance, read_header, read_timeseries, run_id, write_csv, write_json
from .stats import mean_summary, survival_summary, window_mean, group_rng


def extend_absorbed(ts, end):
    """Consensus is absorbing: extend its known state, not a fitted plateau."""
    if ts[-1, 0] >= end:
        return ts
    if ts[-1, 2] != 0 or ts[-1, 1] not in (0., 1.):
        raise ValueError("absorbed output does not end in consensus")
    final = ts[-1].copy()
    final[0] = end
    final[3:] = 0
    return np.vstack([ts, final])


def reduce_run(path, rho_ref=.29, late_fraction=.5, min_time=1000.):
    meta, ts = read_timeseries(path)
    observed = float(meta.get("t_end", ts[-1, 0]))
    horizon = float(meta["T"] if meta["done"] else ts[-1, 0])
    absorption_lower = observed
    if meta['absorbed'] and meta.get('backend')=='CUDA':
        prior = ts[ts[:,4]>0,0]
        absorption_lower = float(prior[-1]) if len(prior) else 0.

    early_absorbed = bool(meta['absorbed'] and observed < min_time <= horizon)
    if early_absorbed:
        # The state at the observation origin is known exactly after consensus.
        # Without this record, extension straight to T would invent a late event.
        ts = extend_absorbed(ts, min_time)
    if meta["absorbed"]:
        ts = extend_absorbed(ts, horizon)
    t, n, rho = ts[:, 0], ts[:, 1], ts[:, 2]
    start = max(min_time, late_fraction * horizon)
    crossed = np.flatnonzero((t >= min_time) & (rho < rho_ref / 2))
    half_event = bool(len(crossed))
    k = int(crossed[0]) if half_event else len(t) - 1
    half_time = float(t[k]) if half_event else horizon
    # Interval contains the actual crossing; don't pretend log-grid outputs are exact.
    half_lower = max(min_time, float(t[max(k - 1, 0)])) if half_event else horizon
    left_censored = early_absorbed or bool(half_event and t[k] == min_time)
    if left_censored:
        half_lower = 0.
    rho_mean = window_mean(t, rho, start, horizon)
    chi_mean = meta["Lx"] * meta["Ly"] * window_mean(t, (2 * n - 1) ** 2, start, horizon)
    mid = (start + horizon) / 2
    early, late = window_mean(t, rho, start, mid), window_mean(t, rho, mid, horizon)
    return dict(path=str(Path(path).resolve()), run=run_id(path), p=meta.get("p", np.nan),
                q=meta.get("q", np.nan), alpha=meta.get("alpha", np.nan),
                parameter_source=meta.get("parameter_source", "unknown"),
                Lx=meta["Lx"], Ly=meta["Ly"], seed=meta["seed"], T=meta["T"],
                done=meta["done"], completion_errors=";".join(meta.get("completion_errors",[])), horizon=horizon, late_start=start,
                rho_all=rho_mean, chi_all=chi_mean,
                rho_early=early, rho_late=late, rho_drift=late - early,
                half_event=half_event, half_time=half_time, half_lower=half_lower,
                half_left_censored=left_censored,
                absorption_event=meta["absorbed"],
                absorption_time=observed if meta["absorbed"] else horizon,
                absorption_lower=absorption_lower if meta['absorbed'] else horizon,
                absorption_observation='scheduled_interval' if meta.get('backend')=='CUDA' else 'exact_step', 
                rho_survivor=rho_mean if not meta["absorbed"] else np.nan,
                chi_survivor=chi_mean if not meta["absorbed"] else np.nan,
                rho_no_half_event=rho_mean if not half_event else np.nan)


def collect(paths, **kwargs):
    return [reduce_run(p, **kwargs) for p in sorted(paths)]


def summarize(rows, boot=2000, seed=1729, stable_groups=False):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["p"], row["q"], row["Lx"], row["Ly"], row["T"])].append(row)
    rng = np.random.default_rng(seed)
    summaries, survival = [], []
    for key, all_rows in sorted(groups.items()):
        if stable_groups:rng=group_rng(key,seed)
        rs = [r for r in all_rows if r["done"]]
        p, q, Lx, Ly, T = key
        base = dict(p=p, q=q, Lx=Lx, Ly=Ly, T=T)
        result = dict(**base, n=len(rs), n_incomplete=len(all_rows) - len(rs),
                      half_left_censored=sum(r.get("half_left_censored",False) for r in rs))
        if not rs:
            continue
        for name in ("rho_all", "rho_survivor", "rho_no_half_event", "chi_all", "chi_survivor", "rho_drift"):
            s = mean_summary([r[name] for r in rs], boot=boot, rng=rng)
            result.update({f"{name}_{k}": v for k, v in s.items()})
        for endpoint in ("half", "absorption"):
            times = np.array([r[f"{endpoint}_time"] for r in rs])
            events = np.array([r[f"{endpoint}_event"] for r in rs])
            lower = np.array([r[f"{endpoint}_lower"] for r in rs])
            # Completed runs share planned T. Censoring at T is exact for the
            # truncated lifetime min(tau,T); interval-observed events retain
            # their last-above/first-below bounds separately from seed noise.
            result[f"{endpoint}_recording_lower_mean"] = float(lower.mean())
            result[f"{endpoint}_recording_upper_mean"] = float(times.mean())
            result[f"{endpoint}_recording_width_mean"] = float((times-lower).mean())
            # If everyone is observed to absorb, S(t)=0 until the planned T.
            horizon = min(float(T), float(times.max()))
            s = survival_summary(times, events, horizon=horizon, boot=boot, rng=rng)
            result.update({f"{endpoint}_{k}": v for k, v in s.items() if k != "curve"})
            for point in s["curve"]:
                survival.append(dict(**base, endpoint=endpoint, **point))
        summaries.append(result)
    return summaries, survival


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("paths", nargs="*", type=Path)
    ap.add_argument("--out", type=Path, default=ROOT / "campaign/v2/results/bulk")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--rho-ref", type=float, default=.29)
    ap.add_argument("--late-fraction", type=float, default=.5)
    args = ap.parse_args()
    paths = args.paths or sorted((ROOT / "data/runs/E13_bulk05/out").glob("*.txt"))
    paths = expand_text_paths(paths)
    # A GPU queue can contain both strips and bulk jobs. Never mix the protocols.
    skipped = [p for p in paths if read_header(p).get("init", "").startswith("strip")
               or read_header(p).get("profile_dir", False)]
    paths = [p for p in paths if p not in skipped]
    runs = collect(paths, rho_ref=args.rho_ref, late_fraction=args.late_fraction)
    summary, survival = summarize(runs, boot=args.boot)
    write_csv(args.out / "runs.csv", runs)
    write_csv(args.out / "summary.csv", summary)
    write_csv(args.out / "survival.csv", survival)
    write_json(args.out / "analysis.json", dict(
        settings=dict(rho_ref=args.rho_ref, late_fraction=args.late_fraction, boot=args.boot,
                      sampling_unit="independent run", averaging="trapezoidal time mean",
                      incomplete_runs="reported, excluded from completed-run summary",
                      half_event="first recorded rho < rho_ref/2 at/after t=1000; earlier absorption is left censored at1000",
                      ci="95% seed percentile bootstrap; survival bands Greenwood log-log"),
        summary=summary, skipped_nonbulk=list(map(str, skipped)), inputs=provenance(paths)))
    print(f"{len(runs)} runs, {len(summary)} groups; output {args.out}")


if __name__ == "__main__":
    main()
