"""Input/provenance helpers. Analysis never imports plotting code."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def run_id(path):
    p = Path(path)
    return p.parent.name if p.name == "out.txt" else p.stem


def expand_text_paths(paths):
    """Both legacy EXP/out/*.txt and nested GPU task/out.txt layouts."""
    found = []
    for path in map(Path, paths):
        if path.is_dir():
            found.extend(path.glob("*.txt"))
            found.extend(path.rglob("out.txt"))
        else:
            found.append(path)
    return sorted(set(found))


def profile_path(path):
    p = Path(path)
    return p.parent / "profile.bin" if p.name == "out.txt" else p.with_suffix(".prof")


def read_header(path):
    """Read actual requested parameters from headers, not rounded run names."""
    out = {"path": str(Path(path).resolve()), "done": False, "absorbed": False}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#"):
            continue
        w = line.split()
        if line.startswith("# Lx"):
            for key in ("Lx", "Ly", "T", "seed", "threads"):
                if key in w:
                    out[key] = int(w[w.index(key) + 1])
        elif line.startswith("# rule"):
            for key in ("alpha", "voter"):
                if key in w:
                    out["p" if key == "voter" else key] = float(w[w.index(key) + 1])
            if "jflip" in w:
                extra = w[w.index("jflip") + 1]
                out["q"] = sum(float(x.split(":")[1]) for x in extra.split(",")
                               if ":" in x and x.split(":")[0] == "2")
                out["parameter_source"] = "header"
        elif line.startswith("# done t_end"):
            out.update(done=True, t_end=int(w[w.index("t_end") + 1]),
                       absorbed=bool(int(w[w.index("absorbed") + 1])))
        elif line.startswith("# init"):
            out["init"] = w[2]
        elif line.startswith("# backend"):
            out["backend"] = w[2]
        elif line.startswith("# profile_dir"):
            a, b = w[2].split(",")
            out.update(a=int(a), b=int(b), ab=f"{a}-{b}", profile_dir=True)
    # Explicitly labelled fallback for older/header-table files; filenames may round q.
    for key in ("p", "q", "alpha"):
        if key not in out:
            tag = "a" if key == "alpha" else key
            m = re.search(r"(?:^|_)" + tag + r"([0-9.eE+-]+)(?:_|$)", Path(path).stem)
            if m:
                out[key] = float(m[1])
                out["parameter_source"] = "filename_fallback"
    m = re.search(r"_ab(\d+)-(\d+)_", Path(path).stem)
    if m and "ab" not in out:
        out["a"], out["b"] = int(m[1]), int(m[2])
        out["ab"] = f"{m[1]}-{m[2]}"
    return out


def read_timeseries(path):
    meta = read_header(path)
    rows = [[float(x) for x in line.split()] for line in Path(path).read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    ts = np.asarray(rows, dtype=float)
    if ts.ndim != 2 or ts.shape[1] != 5:
        raise ValueError(f"{path}: expected five-column normal-mode output")
    if np.any(np.diff(ts[:, 0]) < 0):
        raise ValueError(f"{path}: decreasing time records")
    # Restart output can repeat its boundary record. Keep the last record at each time.
    keep = np.r_[np.diff(ts[:, 0]) != 0, True]
    ts = ts[keep]
    status = completion_status(path, meta=meta, ts=ts)
    meta['header_done'] = meta['done']
    meta['done'] = status['complete']
    meta['completion_errors'] = status['reasons']
    return meta, ts


def read_profile(path, bins):
    """Little-endian dnsim profile, memory mapped; ignore only a partial final record."""
    path = Path(path)
    dtype = np.dtype([("t", "<i8"), ("alive", "<u4", (bins,)),
                      ("iface", "<u4", (bins,))])
    size = path.stat().st_size
    n, remainder = divmod(size, dtype.itemsize)
    if not n:
        raise ValueError(f"{path}: no complete profile records")
    data = np.memmap(path, dtype=dtype, mode="r", shape=(n,))
    if np.any(np.diff(data["t"]) <= 0):
        raise ValueError(f"{path}: profile times are not strictly increasing")
    return data, remainder


def completion_status(path, meta=None, ts=None, profile=None, trailing=0, require_marker=None):
    """Check completion across independently synced text, profile and GPU marker.

    A final text footer can arrive before the profile's last buffered record in
    a live rsync. It is not sufficient evidence of a complete profile. GPU runs
    additionally require the runner's completion marker, written after success.
    """
    path = Path(path)
    meta = read_header(path) if meta is None else meta
    errors = []
    if not meta.get('header_done', meta.get('done', False)):
        errors.append('missing_done_footer')
    end = meta.get('t_end')
    planned = meta.get('T')
    if end is None or planned is None:
        errors.append('missing_end_or_planned_time')
    elif end > planned or (end != planned and not meta.get('absorbed')):
        errors.append('inconsistent_planned_end')
    if ts is not None:
        if not np.all(np.isfinite(ts)) or np.any(ts[:,0]<0) or np.any(ts[:,0]!=np.floor(ts[:,0])):
            errors.append('nonfinite_or_invalid_numeric_record')
        if np.any((ts[:,1:]<0)|(ts[:,1:]>1)):
            errors.append('fraction_outside_unit_interval')
        if not len(ts) or end is None or ts[-1, 0] != end:
            errors.append('text_end_mismatch')
        if meta.get('absorbed') and len(ts):
            if ts[-1, 1] not in (0., 1.) or ts[-1, 2] != 0 or ts[-1, 4] != 0:
                errors.append('absorbed_state_mismatch')
    if profile is not None:
        if trailing:
            errors.append('partial_profile_record')
        if not len(profile) or end is None or profile['t'][-1] != end:
            errors.append('profile_end_mismatch')
        if ts is not None and not np.array_equal(profile['t'], ts[:, 0]):
            errors.append('profile_text_times_mismatch')
        if all(k in meta for k in ('Lx','Ly')) and all(k in profile.dtype.names for k in ('alive','iface')):
            capacity=meta['Lx'];area=meta['Lx']*meta['Ly']
            if np.any(profile['alive']>capacity) or np.any(profile['iface']>capacity):
                errors.append('profile_count_exceeds_bin_capacity')
            if ts is not None and len(profile)==len(ts):
                alive=profile['alive'].sum(axis=1)/area
                iface=profile['iface'].sum(axis=1)/area
                if np.any(np.abs(alive-ts[:,1])>6e-10) or np.any(np.abs(iface-ts[:,4])>6e-10):
                    errors.append('profile_text_fraction_mismatch')
    require_marker = path.name == 'out.txt' if require_marker is None else require_marker
    if require_marker:
        marker = path.parent / 'done.json'
        try:
            done = json.loads(marker.read_text())
            if (done.get('id') != run_id(path) or done.get('t_end') != end
                    or done.get('absorbed') != bool(meta.get('absorbed'))):
                errors.append('completion_marker_mismatch')
        except (OSError, ValueError):
            errors.append('missing_or_invalid_completion_marker')
    return dict(complete=not errors, reasons=errors, t_end=end, T=planned)


def manifest_tasks(paths):
    """Merge manifest revisions by id; reject conflicting scientific parameters."""
    tasks = {}
    science = ('family', 'Lx', 'Ly', 'p', 'q', 'alpha', 'T', 'seed', 'mode', 'orientation')
    for path in map(Path, paths):
        for task in json.loads(path.read_text())['tasks']:
            old = tasks.get(task['id'])
            if old is not None and any(old.get(k) != task.get(k) for k in science):
                raise ValueError(f"conflicting manifest task {task['id']}: {path}")
            tasks[task['id']] = task
    return tasks


def scientific_group(task):
    """Fixed ensemble identity; output cadence and worker assignment are not physics."""
    ab = task.get('orientation') or (1, 0)
    return (float(task['p']), float(task['alpha']), float(task['q']),
            int(task['Lx']), int(task['Ly']), int(task['T']), task['mode'],
            int(ab[0]), int(ab[1]))


def manifest_groups(tasks):
    groups = {}
    for task in tasks.values():
        groups.setdefault(scientific_group(task), []).append(task['id'])
    return {key: sorted(ids) for key, ids in groups.items()}


def complete_ensembles(tasks, status):
    """Eligibility follows all declared ids, not merely all currently copied runs."""
    rows = []
    names = ('p','alpha','q','Lx','Ly','T','mode','a','b')
    for key, ids in sorted(manifest_groups(tasks).items()):
        completed = [i for i in ids if status.get(i, {}).get('complete', False)]
        rows.append(dict(zip(names,key)) | dict(n_planned=len(ids), n_complete=len(completed),
                    ensemble_complete=len(completed)==len(ids), task_ids=ids,
                    incomplete_ids=[i for i in ids if i not in completed]))
    return rows


def provenance(paths, hashes=True):
    result = []
    for path in sorted({str(Path(p).resolve()) for p in paths}):
        p = Path(path)
        item = {"path": path, "bytes": p.stat().st_size}
        if hashes:
            h = hashlib.sha256()
            with p.open("rb") as f:
                for block in iter(lambda: f.read(1024 * 1024), b""):
                    h.update(block)
            item["sha256"] = h.hexdigest()
        result.append(item)
    return result


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean_json(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write_json(path, value):
    path=Path(path);path.parent.mkdir(parents=True, exist_ok=True)
    temporary=None
    try:
        with tempfile.NamedTemporaryFile('w',dir=path.parent,prefix='.'+path.name+'.',delete=False) as stream:
            temporary=Path(stream.name)
            stream.write(json.dumps(clean_json(value),indent=2,allow_nan=False)+'\n')
        os.replace(temporary,path)
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)


def write_csv(path, rows):
    if not rows:
        raise ValueError("refusing to write an empty table")
    path=Path(path);path.parent.mkdir(parents=True, exist_ok=True)
    temporary=None
    try:
        with tempfile.NamedTemporaryFile('w',newline='',dir=path.parent,prefix='.'+path.name+'.',delete=False) as stream:
            temporary=Path(stream.name)
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
        os.replace(temporary,path)
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)
