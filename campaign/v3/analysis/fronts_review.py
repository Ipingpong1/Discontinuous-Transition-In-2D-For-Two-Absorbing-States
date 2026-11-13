"""Reproduce the v3 flagship checks from completed, declared trajectories.

Run ``python -m campaign.v3.analysis.fronts_review`` from the project root.
The four front observables share the same full-window geometric eligibility.
Alternative integrated observables are diagnostics of estimator sensitivity.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import t as student_t

from campaign.v2.analysis.fronts import widths
from campaign.v2.analysis.io import (
    ROOT, completion_status, provenance, read_profile,
    read_timeseries, write_csv, write_json,
)
from campaign.v2.analysis.stats import mean_summary, window_mean

OUT = ROOT / "campaign/v3/results/fronts_review"
MANIFEST_NAMES = (
    "production_manifest_profiles50k.json", "production_manifest_largeL.json",
    "production_manifest_orientation5m.json", "production_manifest_axial10m.json",
)
WINDOWS = ((1e6, 2e6), (2e6, 3.5e6), (3.5e6, 5e6), (2e6, 5e6), (5e6, 1e7))
OBSERVABLES = ("envelope", "live_area", "bond_activity", "site_activity")
D = 0.20
R0 = 1000.0


def run_summary(values):
    result = mean_summary(values)
    if result["n"] > 1:
        half = student_t.ppf(.975, result["n"] - 1) * result["sem"]
        result.update(t95_low=result["mean"] - half, t95_high=result["mean"] + half)
    return result


def load_tasks():
    import json
    tasks = {}
    paths = [ROOT / "campaign/v2/hpc/gpu" / name for name in MANIFEST_NAMES]
    for path in paths:
        for task in json.loads(path.read_text())["tasks"]:
            if (float(task["p"]) == .5 and task["Lx"] >= 8192
                    and task["family"] in ("orientation", "flagship_bulk")):
                if task["id"] in tasks:
                    raise ValueError("Repeated manifest run: " + task["id"])
                tasks[task["id"]] = task
    return tasks, paths


def read_completed(task, inputs):
    path = ROOT / "data/v2/gpu/runs" / task["family"] / task["id"] / "out.txt"
    meta, series = read_timeseries(path)
    for name in ("Lx", "Ly", "T", "seed"):
        if meta[name] != task[name]:
            raise ValueError(f"Manifest mismatch {path}: {name}")
    if (meta["p"], meta["q"]) != (float(task["p"]), float(task["q"])):
        raise ValueError(f"Parameter mismatch: {path}")
    status = completion_status(path, meta, series)
    if not status["complete"]:
        raise ValueError(f"Incomplete declared run {path}: {status['reasons']}")
    inputs.extend([path, path.parent / "done.json"])
    return path, meta, series


def reduce_fronts(tasks, inputs):
    rows = []
    for task in tasks.values():
        if task["family"] != "orientation":
            continue
        path, meta, series = read_completed(task, inputs)
        length = meta["Lx"]
        profile, trailing = read_profile(path.parent / "profile.bin", length)
        status = completion_status(path, meta, series, profile, trailing)
        if not status["complete"]:
            raise ValueError(f"Incomplete profile {path}: {status['reasons']}")
        inputs.append(path.parent / "profile.bin")
        a, b = task["orientation"]
        normal_scale = np.hypot(a, b)
        smooth_bins = int(round(9 * normal_scale))
        smooth_bins += 1 - smooth_bins % 2
        tracked = widths(profile["t"], profile["iface"], length, .3, smooth_bins)
        time = series[:, 0]
        # Fixed positive normalizations leave a velocity zero invariant.
        measures = dict(
            envelope=tracked["envelope"] / normal_scale,
            live_area=2 * series[:, 1] * length / normal_scale,
            bond_activity=series[:, 2] * length / (.2866 * normal_scale),
            site_activity=series[:, 4] * length / (.705 * normal_scale),
        )
        for start, end in WINDOWS:
            if end > meta["T"]:
                continue
            selected = (time >= start) & (time <= end)
            if time[selected][0] != start or time[selected][-1] != end:
                raise ValueError(f"Missing window endpoint: {path}")
            envelope = tracked["envelope"][selected]
            eligible = bool(np.all((envelope > .15 * length) & (envelope < .85 * length)))
            for name, values in measures.items():
                rows.append(dict(
                    run=task["id"], L=length, T=meta["T"], q=meta["q"],
                    normal=f"{a}-{b}", start=start, end=end, observable=name,
                    v=float(np.polyfit(time[selected], values[selected], 1)[0] / 2),
                    eligible=eligible,
                ))
    return rows


def crossing_with_influences(group):
    q = np.array(sorted(group))
    means = np.array([np.mean(list(group[x].values())) for x in q])
    brackets = np.flatnonzero(means[:-1] * means[1:] < 0)
    if len(brackets) != 1:
        return None
    i = brackets[0]
    left, right = means[i:i + 2]
    dq = q[i + 1] - q[i]
    root = q[i] - left * dq / (right - left)
    gradient = (-right * dq / (right - left)**2, left * dq / (right - left)**2)
    influences = {}
    for x, derivative in zip(q[i:i + 2], gradient):
        mean = np.mean(list(group[x].values()))
        for run, value in group[x].items():
            influences[run] = derivative * (value - mean)
    return float(root), influences


def influence_interval(influences, estimate, multiplicity=1):
    # Each q/direction/cohort is an independent group. A repeated run carries
    # its paired early/late contribution before the variance is calculated.
    groups = defaultdict(list)
    for run, value in influences.items():
        groups[run.rsplit("_s", 1)[0]].append(value)
    components = np.array([np.var(v, ddof=1) / len(v) for v in groups.values()])
    degrees = np.array([len(v) - 1 for v in groups.values()])
    variance = components.sum()
    df = variance**2 / np.sum(components**2 / degrees)
    sem = float(np.sqrt(variance))
    half = student_t.ppf(1 - .025 / multiplicity, df) * sem
    return dict(estimate=estimate, sem=sem, degrees_of_freedom=float(df),
                low=float(estimate - half), high=float(estimate + half),
                multiplicity=multiplicity)


def front_summaries(rows):
    all_groups = defaultdict(list)
    for row in rows:
        key = tuple(row[k] for k in ("L", "T", "normal", "start", "end", "observable", "q"))
        all_groups[key].append(row)
    eligible, points = {}, []
    for key, members in sorted(all_groups.items()):
        if len(members) != 4:
            raise ValueError(f"Expected four declared seeds: {key}")
        if not all(row["eligible"] for row in members):
            continue
        eligible[key] = {row["run"]: row["v"] for row in members}
        points.append(dict(zip(("L", "T", "normal", "start", "end", "observable", "q"), key),
                           **run_summary([row["v"] for row in members])))
    curves = defaultdict(dict)
    for key, values in eligible.items():
        curves[key[:-1]][key[-1]] = values
    roots, estimates = [], {}
    for key, group in sorted(curves.items()):
        root = crossing_with_influences(group)
        if root is None:
            continue
        estimates[key] = root
        roots.append(dict(zip(("L", "T", "normal", "start", "end", "observable"), key),
                          **influence_interval(root[1], root[0])))
    contrasts = []
    specifications = [
        ("axial_minus_tilted_16384", (16384, 2000000, "1-0", 1e6, 2e6),
         (16384, 2000000, "2-1", 1e6, 2e6), 3),
        ("axial_5m_late_minus_early", (8192, 5000000, "1-0", 3.5e6, 5e6),
         (8192, 5000000, "1-0", 2e6, 3.5e6), 1),
        ("axial_10m_late_minus_early", (8192, 10000000, "1-0", 5e6, 1e7),
         (8192, 10000000, "1-0", 2e6, 5e6), 1),
    ]
    for name, first, second, multiplicity in specifications:
        for observable in OBSERVABLES:
            left, right = estimates[first + (observable,)], estimates[second + (observable,)]
            influence = dict(left[1])
            for run, value in right[1].items():
                influence[run] = influence.get(run, 0.) - value
            contrasts.append(dict(name=name, observable=observable,
                                  **influence_interval(influence, left[0] - right[0], multiplicity)))
    return points, roots, contrasts


def reduce_bulk(tasks, inputs):
    rows, curves = [], []
    for task in tasks.values():
        if task["family"] != "flagship_bulk":
            continue
        _, meta, series = read_completed(task, inputs)
        if meta["T"] != 5000000:
            raise ValueError("Unexpected bulk horizon")
        time, rho = series[:, 0], series[:, 2]
        if time[-1] < meta["T"]:
            time, rho = np.r_[time, meta["T"]], np.r_[rho, 0.]
        rows.append(dict(run=task["id"], L=meta["Lx"], q=meta["q"],
                         mean=window_mean(time, rho, 2.5e6, 5e6),
                         drift=window_mean(time, rho, 3.75e6, 5e6)
                         - window_mean(time, rho, 2.5e6, 3.75e6)))
        curves.extend(dict(run=task["id"], L=meta["Lx"], q=meta["q"], time=t, rho=r)
                      for t, r in zip(time, rho))
    groups = defaultdict(list)
    for row in rows:
        groups[row["L"], row["q"]].append(row)
    summary = []
    for (length, q), members in sorted(groups.items()):
        if len(members) != 4:
            raise ValueError("Incomplete bulk ensemble")
        summary.append(dict(L=length, q=q, **run_summary([r["mean"] for r in members]),
                            drift_mean=float(np.mean([r["drift"] for r in members]))))
    return rows, curves, summary


def curvature_solution(time, initial_radius, signed_velocity, diffusion=D):
    result = solve_ivp(lambda _, radius: signed_velocity - diffusion / radius,
                       (time[0], time[-1]), [initial_radius], t_eval=time,
                       rtol=1e-10, atol=1e-9)
    if not result.success or len(result.t) != len(time) or np.any(result.y <= 0):
        raise ValueError("Curvature trajectory reaches zero or integration failed")
    # The autonomous equation has an analytic first integral. This checks the
    # plotted integration independently of the numerical solver's status.
    radius = result.y[0]
    if signed_velocity == 0:
        residual = radius**2 - initial_radius**2 + 2 * diffusion * (time - time[0])
    else:
        integral = ((radius - initial_radius) / signed_velocity
                    + diffusion / signed_velocity**2
                    * np.log(np.abs((signed_velocity * radius - diffusion)
                                    / (signed_velocity * initial_radius - diffusion))))
        residual = integral - (time - time[0])
    if np.max(np.abs(residual)) > .02:
        raise ValueError("Curvature first-integral check failed")
    return radius


def reduce_droplets(inputs):
    source = ROOT / "campaign/v2/results/fronts/points.csv"
    with source.open() as stream:
        front = [r for r in csv.DictReader(stream) if float(r["p"]) == .5
                 and int(r["L"]) == 4096 and r["ab"] == "1-0"
                 and float(r["threshold"]) == .3 and float(r["t0"]) == 20000
                 and r["estimator"] == "envelope" and r["t1"] == "inf"]
    front.sort(key=lambda r: float(r["q"]))
    q_grid, v_grid = np.array([(float(r["q"]), float(r["v"])) for r in front]).T
    inputs.append(source)
    rows, records, groups = [], [], defaultdict(list)
    for path in sorted((ROOT / "data/runs/R2_drop2/out").glob("*.txt")):
        meta, series = read_timeseries(path)
        profile, trailing = read_profile(path.with_suffix(".prof"), meta["Lx"])
        # Legacy text outputs and profiles have different schedules.
        if not completion_status(path, meta, series)["complete"] or trailing:
            raise ValueError(f"Incomplete droplet: {path}")
        if profile["t"][-1] != meta["T"]:
            raise ValueError(f"Incomplete droplet profile: {path}")
        length, q = meta["Lx"], meta["q"]
        kind = "disk" if "_diskR" in path.name else "hole"
        chord = 2 * profile["alive"].astype(float)
        if kind == "hole":
            chord = length - chord
        crop = abs(np.arange(length) - length / 2) < 1.3 * R0 + 50
        areas = dict(trimmed_clipped=np.maximum(chord[:, crop], 0).sum(1),
                     full_signed=chord.sum(1), trimmed_signed=chord[:, crop].sum(1))
        time = profile["t"]
        selected = (time >= 20000) & (time <= 200000)
        if time[selected][0] != 20000 or time[selected][-1] != 200000:
            raise ValueError("Missing droplet fit endpoints")
        for estimator, area in areas.items():
            if np.any(area <= 0):
                raise ValueError("Nonpositive signed droplet area")
            radius = np.sqrt(area / np.pi)
            row = dict(run=path.stem, q=q, kind=kind, estimator=estimator,
                       slope=float(np.polyfit(time[selected], radius[selected], 1)[0]),
                       start=20000, end=200000, radius_final=float(radius[-1]))
            rows.append(row)
            if estimator == "trimmed_clipped":
                records.extend(dict(run=path.stem, q=q, kind=kind, time=t, radius=r)
                               for t, r in zip(time, radius))
                groups[q, kind].append((time, radius))
        inputs.extend([path, path.with_suffix(".prof")])
    predictions, summary = [], []
    for (q, kind), members in sorted(groups.items()):
        if len(members) != 4 or not all(np.array_equal(members[0][0], r[0]) for r in members):
            raise ValueError("Droplet ensemble/schedule mismatch")
        if not q_grid.min() <= q <= q_grid.max():
            raise ValueError("Front extrapolation forbidden")
        time = members[0][0]
        mean_radius = np.mean([r[1] for r in members], axis=0)
        selected = time >= 20000
        front_v = float(np.interp(q, q_grid, v_grid))
        signed_v = front_v if kind == "disk" else -front_v
        predicted = curvature_solution(time[selected], mean_radius[selected][0], signed_v)
        prediction_slope = float(np.polyfit(time[selected], predicted, 1)[0])
        predictions.extend(dict(q=q, kind=kind, time=t, radius=r)
                           for t, r in zip(time[selected], predicted))
        for estimator in areas:
            slopes = [r["slope"] for r in rows if (r["q"], r["kind"], r["estimator"]) == (q, kind, estimator)]
            summary.append(dict(q=q, kind=kind, estimator=estimator,
                                **run_summary(slopes), front_v=front_v,
                                ode_slope=prediction_slope,
                                mean_minus_ode=float(np.mean(slopes) - prediction_slope)))
    return rows, records, predictions, summary


def verify_preserved_measurements(fronts, bulk, droplets, inputs):
    """Catch reader, time-window, and estimator changes against published v2 rows."""
    paths = [ROOT / "campaign/v2/results/gpu/front_runs.csv",
             ROOT / "campaign/v2/results/gpu/bulk_runs.csv",
             ROOT / "campaign/v2/results/droplets/two_drop_rates.csv"]
    with paths[0].open() as stream:
        reference_fronts = {(r["run"], float(r["t0"]), float(r["t1"])): float(r["v"])
                            for r in csv.DictReader(stream)
                            if r["estimator"] == "envelope" and float(r["threshold"]) == .3}
    with paths[1].open() as stream:
        reference_bulk = {r["run"]: float(r["rho_all"]) for r in csv.DictReader(stream)}
    with paths[2].open() as stream:
        reference_drops = {r["run"]: float(r["radius_slope"]) for r in csv.DictReader(stream)}
    checks = dict(
        front_envelope_max_difference=max(abs(r["v"] - reference_fronts[r["run"], r["start"], r["end"]])
                                          for r in fronts if r["observable"] == "envelope" and r["eligible"]),
        bulk_time_mean_max_difference=max(abs(r["mean"] - reference_bulk[r["run"]]) for r in bulk),
        primary_droplet_slope_max_difference=max(abs(r["slope"] - reference_drops[r["run"]])
                                                 for r in droplets if r["estimator"] == "trimmed_clipped"),
    )
    if max(checks.values()) > 1e-12:
        raise ValueError(f"Existing measurements changed unexpectedly: {checks}")
    inputs.extend(paths)
    return checks


def main():
    tasks, inputs = load_tasks()
    fronts = reduce_fronts(tasks, inputs)
    points, roots, contrasts = front_summaries(fronts)
    bulk, bulk_curves, bulk_summary = reduce_bulk(tasks, inputs)
    droplets, drop_curves, predictions, drop_summary = reduce_droplets(inputs)
    checks = verify_preserved_measurements(fronts, bulk, droplets, inputs)
    outputs = dict(front_runs=fronts, front_points=points, front_crossings=roots,
                   front_contrasts=contrasts, bulk_runs=bulk, bulk_curves=bulk_curves,
                   bulk_summary=bulk_summary, droplet_runs=droplets,
                   droplet_curves=drop_curves, droplet_predictions=predictions,
                   droplet_summary=drop_summary)
    OUT.mkdir(parents=True, exist_ok=True)
    for name, rows in outputs.items():
        write_csv(OUT / (name + ".csv"), rows)
    sources = [Path(__file__), ROOT / "campaign/v2/analysis/fronts.py",
               ROOT / "campaign/v2/analysis/io.py", ROOT / "campaign/v2/analysis/stats.py"]
    hashes = provenance(inputs + sources)
    for item in hashes:
        item["path"] = str(Path(item["path"]).relative_to(ROOT))
    write_json(OUT / "analysis.json", dict(
        front_normalizations=dict(bond_activity=.2866, site_activity=.705),
        front_geometry="All four declared runs must have envelope width in (0.15L,0.85L) throughout the fixed window. Same gate for every proxy.",
        front_uncertainty="Delta method conditional on observed bracketing pair, with Welch-Satterthwaite Student-t intervals. Sensitivity analysis; no small-sample coverage guarantee. Paired times use matched runs. Direction contrast uses Bonferroni over three pairs.",
        mean_uncertainty="Empirical percentile bootstrap plus Student-t intervals over independent runs. Neither assures coverage at n=4.",
        droplet_estimator="R=sqrt(A/pi). Primary A=sum max(c_y,0) over |y-L/2|<1350, c_y=2n_y for disk and L-2n_y for hole. Sensitivities use full signed and trimmed signed sums.",
        droplet_model=dict(D=D, initial_radius=R0, fit_window=[20000, 200000],
                           prediction="Integrate dR/dt=+v-D/R for disk, -v-D/R for hole, initialized to the primary mean R at 20000; same parameters for full-signed comparison."),
        contrasts=contrasts, verification=checks,
        output_rows={k: len(v) for k, v in outputs.items()},
        provenance=hashes,
    ))
    print(f"Wrote {len(outputs)} tables for {len(tasks)} declared GPU runs and 32 droplet trajectories to {OUT}")


if __name__ == "__main__":
    main()
