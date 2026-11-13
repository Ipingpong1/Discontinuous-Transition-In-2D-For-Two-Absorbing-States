"""Rebuild the five-cell comparison from recorded profiles and lifetime runs.

Run from the project root with ``python3 -B campaign/v3/analysis/molofsky.py``.
The profile slopes use the original 32-row smoothing and geometric filters.
The three positive-branch groups on each cut are selected by their ensemble
mean interior decline of at most 0.03 between 30000 and 100000 steps. This
operational criterion defines the extrapolation sample, without assigning
asymptotic stability. No simulations or earlier analysis outputs are changed.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re

import numpy as np

PROFILE_NAME = re.compile(r"p2_([.0-9]+)_p1_([.0-9]+)_s([0-9]+)")
METHODS = ("fixed", "adaptive", "integrated")
WINDOWS = {"early": 2000, "late": 30000}
SLOPES = tuple(f"{window}_{method}" for window in WINDOWS for method in METHODS)


def profile_run(path: Path) -> dict:
    header = path.with_suffix(".txt").read_text().splitlines()
    dimensions = next(line.split() for line in header if line.startswith("# Lx"))
    lx, ly, horizon = int(dimensions[2]), int(dimensions[4]), int(dimensions[6])
    if not any(line.startswith("# done") for line in header):
        raise ValueError(f"Incomplete run: {path}")
    dtype = np.dtype([("t", "<i8"), ("alive", "<u4", ly), ("iface", "<u4", ly)])
    if path.stat().st_size % dtype.itemsize:
        raise ValueError(f"Partial profile record: {path}")
    records = np.fromfile(path, dtype=dtype)
    time = records["t"]
    if time[0] != 0 or time[-1] != horizon or horizon != 100000:
        raise ValueError(f"Unexpected observation horizon: {path}")
    if np.any(np.diff(time) <= 0):
        raise ValueError(f"Non-increasing profile times: {path}")
    density = np.roll(records["iface"] / lx, ly // 4, axis=1)
    extended = np.concatenate((density[:, -16:], density, density[:, :15]), axis=1)
    sums = np.pad(extended.cumsum(axis=1), ((0, 0), (1, 0)))
    smooth = (sums[:, 32:] - sums[:, :-32]) / 32
    center_index = ly // 2
    reference = density[time.searchsorted(1000),
                        center_index-ly//8:center_index+ly//8].mean()
    center = density[:, center_index-ly//16:center_index+ly//16].mean(axis=1)

    def width(threshold):
        edges = smooth < threshold[:, None]
        left, right = edges[:, :center_index][:, ::-1], edges[:, center_index:]
        return (np.where(left.any(axis=1), left.argmax(axis=1), center_index)
                + np.where(right.any(axis=1), right.argmax(axis=1), ly-center_index))

    fixed = width(np.full(len(time), 0.5*reference))
    adaptive = width(0.5*center)
    integrated = np.divide(density.sum(axis=1), center,
                           out=np.full(len(time), np.nan), where=center > 0)
    p2, p1, seed = PROFILE_NAME.fullmatch(path.stem).groups()
    row = dict(p2=float(p2), p1=float(p1), Lx=lx, Ly=ly, seed=int(seed),
               horizon=horizon, reference=reference,
               core_30k=float(np.interp(30000, time, center/reference)),
               core_100k=float(center[-1]/reference))
    row["core_drop"] = row["core_30k"] - row["core_100k"]
    for window, start in WINDOWS.items():
        select = (time >= start) & (fixed > 0.1*ly) & (fixed < 0.9*ly) & (center > 0.6*reference)
        row[f"{window}_records"] = int(select.sum())
        row[f"{window}_last_time"] = int(time[select][-1]) if select.any() else None
        for method, widths in zip(METHODS, (fixed, adaptive, integrated)):
            row[f"{window}_{method}"] = (float(np.polyfit(time[select], widths[select], 1)[0]/2)
                                          if select.sum() >= 4 else None)
    return row


def summarize_profiles(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["p2"], row["p1"], row["Lx"]].append(row)
    summaries = []
    for (p2, p1, lx), runs in sorted(grouped.items()):
        result = dict(p2=p2, p1=p1, Lx=lx, runs=len(runs))
        for field in ("core_30k", "core_100k", "core_drop", *SLOPES):
            values = np.array([run[field] for run in runs if run[field] is not None])
            result[field] = float(values.mean()) if len(values) else None
            result[field+"_sem"] = float(values.std(ddof=1)/np.sqrt(len(values))) if len(values)>1 else None
            result[field+"_n"] = len(values)
        result["core_pass"] = result["core_drop"] <= 0.03
        result["extrapolation_sample"] = bool(lx == 2048 and len(runs) == 4
            and result["core_pass"] and all(result[field+"_n"] == 4 for field in SLOPES))
        summaries.append(result)
    return summaries


def extrapolate(groups):
    fits = []
    for p2 in (0.45, 0.50):
        selected = [group for group in groups if group["p2"] == p2 and group["extrapolation_sample"]]
        if len(selected) != 3:
            raise ValueError(f"Unexpected extrapolation sample on p2={p2}")
        probabilities = np.array([group["p1"] for group in selected])
        for field in SLOPES:
            slopes = np.array([group[field] for group in selected])
            gradient, intercept = np.polyfit(probabilities, slopes, 1)
            zero = -intercept/gradient
            fits.append(dict(p2=p2, estimator=field, points=len(selected),
                             p1_min=float(probabilities.min()), p1_max=float(probabilities.max()),
                             gradient=float(gradient), intercept=float(intercept),
                             extrapolated_zero=float(zero),
                             angle=float(np.degrees(np.arctan2(p2-0.4, zero-0.2))),
                             bracketed=bool(slopes.min() <= 0 <= slopes.max())))
    return fits


def lifetimes(project):
    rows, files = [], []
    for p2, p1 in ((0.45, 0.18535), (0.50, 0.17057)):
        for size in (256, 512):
            observations = []
            pattern = f"life_p2_{p2:.2f}_p1_{p1:.5f}_L{size}_s*.txt"
            for path in sorted((project/"campaign/molofsky_sim/bulk_fine").glob(pattern)):
                lines = path.read_text().splitlines()
                dimensions = next(line.split() for line in lines if line.startswith("# Lx"))
                footer = lines[-1].split()
                if int(dimensions[6]) != 1000000 or footer[1] != "done":
                    raise ValueError(f"Unexpected lifetime record: {path}")
                observations.append((int(footer[3]), int(footer[5])))
                files.append(path)
            if len(observations) != 4:
                raise ValueError(f"Incomplete four-run lifetime cohort: {pattern}")
            times = np.array([time for time, _ in observations])
            events = sum(event for _, event in observations)
            rows.append(dict(p2=p2, p1=p1, L=size, runs=4, events=events,
                             horizon=1000000, restricted_mean=float(times.mean()),
                             restricted_mean_sem=float(times.std(ddof=1)/2),
                             observed_times=";".join(str(time) for time in times)))
    return rows, files


def lineage_fractions(path):
    fractions = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields[0] == "T":
            time = int(fields[1])
            fractions[time] = [0, 0]
        else:
            count = int(fields[0])
            fractions[time][0] += count
            if len(fields[1].split(",")) >= 4:
                fractions[time][1] += count
    return {time: high/total for time, (total, high) in fractions.items()}


def write_csv(path, rows):
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    project = args.root.resolve()
    output = project/"campaign/v3/results/molofsky"
    output.mkdir(parents=True, exist_ok=True)
    profiles = sorted((project/"campaign/molofsky_sim/scan_fine").rglob("p2_*_s*.prof"))
    rows = [dict(path=str(path.relative_to(project)), **profile_run(path)) for path in profiles]
    groups = summarize_profiles(rows)
    fits = extrapolate(groups)
    life, life_files = lifetimes(project)
    paths = [project/"campaign/theory/patterns.txt", project/"campaign/theory/patterns5_1e7.txt"]
    moore, five = map(lineage_fractions, paths)
    lineage = [dict(T=time, moore_ge5_given_ge3=moore[time], five_ge5_given_ge3=five[time],
                    ratio=moore[time]/five[time]) for time in (16, 64, 1024, 65536)]
    for name, records in (("runs", rows), ("groups", groups), ("extrapolations", fits),
                          ("lifetimes", life), ("lineages", lineage)):
        write_csv(output/(name+".csv"), records)
    inputs = [*profiles, *(path.with_suffix(".txt") for path in profiles), *life_files, *paths]
    provenance = dict(
        analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        inputs={str(path.relative_to(project)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in inputs},
        numpy_version=np.__version__,
        profile_method="32-row smoothing; first crossings from center; fixed and contemporaneous half-core thresholds; integrated width",
        profile_filters="fixed width between 0.1 and 0.9 Ly; current core above 0.6 reference; >=4 selected records",
        finite_window_criterion="ensemble mean core drop <=0.03 from 30000 to 100000 steps",
        fit="equal weights on three four-run Lx=2048 positive-branch groups, each with every estimator available",
        uncertainties="slope bars are SEM over independent runs at each fixed parameter; estimator spread is sensitivity, not a confidence interval",
        lifetime_endpoint="all-run mean min(T_absorb,1000000); 3/4 events at L=512 and 4/4 at L=256",
        extrapolations=fits,
    )
    (output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(f"Rebuilt {len(rows)} profile runs, {len(groups)} groups and {len(life)} lifetime cohorts in {output}")


if __name__ == "__main__":
    main()
