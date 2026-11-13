#!/usr/bin/env python3
"""GPU worker v2: jittered all-phase sampling, rolling 120-second utilization.

Stop only when a mature window (>=40 samples over>=110s) remains below80%
for30s. All post-startup samples, including profile/checkpoint I/O, enter the
window and job mean. This avoids aliasing periodic I/O with fixed10s polls.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
from collections import deque
import random
import hashlib
import json
import os
import signal
import struct
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--manifest", type=Path, default=HERE / "wetting_manifest.json")
parser.add_argument("--family", choices=["wetting"], required=True)
parser.add_argument("--binary", type=Path, required=True)
parser.add_argument("--root", type=Path, default=Path("/scratch/mo2273/dn/v2/gpu"))
parser.add_argument("--deadline", help="ISO UTC deadline, overrides manifest")
parser.add_argument("--gpu-index", default=os.environ.get("DNSIM_GPU_INDEX"),
                    help="optional nvidia-smi selector; otherwise require exactly one visible GPU")
parser.add_argument("--sample-seconds", type=float, default=2.5)
args = parser.parse_args()
manifest = json.loads(args.manifest.read_text())
deadline = datetime.fromisoformat(args.deadline or manifest["deadline_utc"]).timestamp()
job = os.environ.get("SLURM_JOB_ID", f"local-{os.getpid()}")
worker = args.root / "workers" / job
worker.mkdir(parents=True, exist_ok=True)
events = worker / "events.jsonl"
util_file = worker / "utilization.jsonl"
stopping = False

def utc():
    return datetime.now(timezone.utc).isoformat()

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def append(path, obj):
    with path.open("a") as stream:
        stream.write(json.dumps(obj, sort_keys=True) + "\n")
        stream.flush()

def record(status, **values):
    obj = dict(utc=utc(), status=status, job=job, family=args.family, **values)
    append(events, obj)
    print(json.dumps(obj, sort_keys=True), flush=True)

def atomic_json(path, obj):
    temporary = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)

def on_signal(signum, frame):
    global stopping
    stopping = True

signal.signal(signal.SIGTERM, on_signal)
signal.signal(signal.SIGINT, on_signal)
signal.signal(signal.SIGUSR1, on_signal)

def terminate(process):
    if process.poll() is None:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()

def make_band(path, task):
    import numpy as np
    L = task["Lx"]
    a, b = task["orientation"]
    rng = np.random.Generator(np.random.PCG64(task["init_seed"]))
    x = np.arange(L, dtype=np.int64)[None, :]
    temporary = path.with_suffix(".tmp")
    with temporary.open("wb") as stream:
        stream.write(b"DNSIM1\0\0")
        stream.write(struct.pack("=QQQ", L, L, 0))
        for y0 in range(0, L, 64):
            y = np.arange(y0, min(y0 + 64, L), dtype=np.int64)[:, None]
            inside = (a * y - b * x) % L < L // 2
            alive = inside & (rng.random(inside.shape) < .5)
            stream.write(np.packbits(alive, axis=1, bitorder="little").tobytes())
    os.replace(temporary, path)
    return dict(initial_sha256=digest(path), numpy_version=np.__version__,
                initial_rng="NumPy PCG64", initial_seed=task["init_seed"])

def completed(directory, task):
    done = directory / "done.json"
    if not done.exists():
        return False
    try:
        data = json.loads(done.read_text())
    except (ValueError, OSError):
        return False
    return data.get("id") == task["id"] and (data.get("t_end") == task["T"] or data.get("absorbed") is True)

def end_stats(path):
    if not path.exists():
        return None, False
    with path.open("rb") as stream:
        stream.seek(max(0, path.stat().st_size - 4096))
        lines = stream.read().decode(errors="replace").splitlines()
    for line in reversed(lines):
        if line.startswith("# done t_end "):
            fields = line.split()
            return int(fields[3]), fields[5] == "1"
    return None, False

def phase(path):
    try:
        fields = path.read_text().split()
        return fields[0], float(fields[1])
    except (OSError, ValueError, IndexError):
        return "unknown", 0.0

def utilization():
    query = ["nvidia-smi"] + (["-i", args.gpu_index] if args.gpu_index else []) + [
             "--query-gpu=uuid,utilization.gpu,utilization.memory,memory.used,power.draw",
             "--format=csv,noheader,nounits"]
    result = subprocess.run(query, check=True, capture_output=True, text=True, timeout=5)
    rows = [line for line in result.stdout.splitlines() if line.strip()]
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one accessible GPU, found {len(rows)}; set DNSIM_GPU_INDEX explicitly")
    values = rows[0].split(",")
    return dict(gpu_uuid=values[0].strip(), gpu_util=float(values[1]),
                memory_util=float(values[2]), memory_mib=float(values[3]), power_w=float(values[4]))

binary = args.binary.resolve()
source = (HERE / "../../../round2/dnsim.c").resolve()
if digest(source) != manifest["source_sha256"]:
    raise SystemExit("Original CPU source hash differs from manifest; regenerate and review.")
binary_hash = digest(binary)
provenance = dict(job=job, family=args.family, started_utc=utc(), binary=str(binary),
                  binary_sha256=binary_hash, manifest_sha256=digest(args.manifest),
                  source_sha256=digest(source), source_files={p.name:digest(p) for p in
                  [HERE/"backend.cu", HERE/"generate.py", HERE/"host_loop.inc",
                   HERE/"gpu_api.h", Path(__file__), HERE/"wetting_manifest.json"]},
                  deadline_utc=datetime.fromtimestamp(deadline, timezone.utc).isoformat(),
                  slurm={k:v for k,v in os.environ.items() if k.startswith("SLURM_")},
                  cuda_visible_devices=os.environ.get("CUDA_VISIBLE_DEVICES"))
atomic_json(worker / "provenance.json", provenance)
record("worker_start", binary_sha256=binary_hash)
queue = args.root / "claims" / args.family
queue.mkdir(parents=True, exist_ok=True)
exit_code = 0
monitor_rng = random.Random("dn-v2-utilization:" + job)
job_util_sum = 0.0
job_util_count = 0
for task in manifest["tasks"]:
    if task["family"] != args.family:
        continue
    if stopping or time.time() >= deadline - 60:
        record("deadline_or_signal_no_new_run")
        break
    directory = args.root / "runs" / args.family / task["id"]
    if completed(directory, task):
        continue
    claim = queue / task["id"]
    try:
        claim.mkdir()
    except FileExistsError:
        continue
    directory.mkdir(parents=True, exist_ok=True)
    atomic_json(claim / "owner.json", dict(job=job, utc=utc(), pid=os.getpid()))
    phase_file = directory / "phase.txt"
    phase_file.write_text(f"initialization {time.monotonic()}\n")
    checkpoint = directory / "checkpoint.bin"
    initial_info = {}
    command = [str(binary), "--Lx", str(task["Lx"]), "--Ly", str(task["Ly"]),
               "--T", str(task["T"]), "--seed", str(task["seed"]), "--alpha", task["alpha"],
               "--voter", task["p"], "--jflip", "2:" + task["q"], "--ppd", "2",
               "--every", str(task["every"]), "--out", str(directory/"out.txt"),
               "--ckpt", str(checkpoint), "--ckpt-every", str(task["checkpoint_every"]),
               "--resume", "--stop-absorbed", "--save", str(directory/"final.bin")]
    if task["profile"]:
        command += ["--profile", str(directory/"profile.bin")]
    if task["mode"] == "band":
        initial = directory / "initial.bin"
        if not initial.exists():
            initial_info = make_band(initial, task)
        else:
            initial_info = dict(initial_sha256=digest(initial), initial_seed=task["init_seed"])
        a, b = task["orientation"]
        command += ["--load", str(initial), "--pdir", f"{a},{b}"]
    elif task["mode"] == "wall":
        command += ["--init", "strip:1:0.5", "--bg", "0"]
    elif task["mode"] == "strip":
        command += ["--init", "strip:0.5:0.5", "--bg", "0"]
    else:
        command += ["--init", "random:0.5"]
    metadata = dict(task=task, command=command, job=job, binary_sha256=binary_hash,
                    source_sha256=manifest["source_sha256"], started_utc=utc(), **initial_info)
    atomic_json(directory / f"run_{job}.json", metadata)
    environment = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="1",
                       OMP_PROC_BIND="false", DNSIM_PHASE_FILE=str(phase_file))
    record("run_start", id=task["id"], command=command)
    started = time.monotonic()
    utilization_window = deque()
    low_since = None
    status = "finished"
    with (directory / f"process_{job}.log").open("a") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   env=environment, start_new_session=True)
        next_sample = time.monotonic() + args.sample_seconds
        while process.poll() is None:
            if stopping or time.time() >= deadline:
                status = "deadline_or_signal_checkpoint_preserved"
                terminate(process)
                break
            now = time.monotonic()
            if now >= next_sample:
                current_phase, phase_started = phase(phase_file)
                # NVIDIA utilization is a recent interval average, not an
                # instantaneous kernel-phase indicator. Include I/O samples and
                # de-alias the periodic profile schedule by jittering reads.
                scientific = now - started >= 20
                sample = dict(utc=utc(), job=job, id=task["id"], run_age_s=now-started,
                              phase=current_phase, phase_age_s=now-phase_started,
                              scientific=scientific, monitor_version=2)
                try:
                    sample.update(utilization())
                    if scientific:
                        value = sample["gpu_util"]
                        job_util_sum += value
                        job_util_count += 1
                        utilization_window.append((now, value))
                    while utilization_window and utilization_window[0][0] < now - 120:
                        utilization_window.popleft()
                    count = len(utilization_window)
                    span = now-utilization_window[0][0] if count else 0.
                    mean = sum(v for _,v in utilization_window)/count if count else None
                    mature = count >= 40 and span >= 110
                    if mature and mean < 80:
                        if low_since is None:
                            low_since = now
                    else:
                        low_since = None
                    sample.update(window_samples=count, window_span_s=span,
                                  window_mean_util=mean, window_mature=mature,
                                  low_mean_duration_s=now-low_since if low_since is not None else 0.,
                                  job_poststartup_mean_util=job_util_sum/job_util_count if job_util_count else None,
                                  job_poststartup_samples=job_util_count)
                    append(util_file, sample)
                    if low_since is not None and now-low_since >= 30:
                        status = "sustained_low_gpu_utilization_stop"
                        terminate(process)
                        stopping = True
                        exit_code = 3
                        break
                except (OSError, ValueError, subprocess.SubprocessError) as error:
                    sample["monitor_error"] = str(error)
                    append(util_file, sample)
                    record("monitor_failure", id=task["id"], error=str(error))
                    status = "monitor_failure_stop"
                    terminate(process)
                    stopping = True
                    exit_code = 4
                    break
                next_sample = now + monitor_rng.uniform(.8, 1.2) * args.sample_seconds
            time.sleep(.25)
        returncode = process.wait()
    finished_time, absorbed = end_stats(directory/"out.txt")
    success = returncode == 0 and (finished_time == task["T"] or absorbed)
    record("run_end", id=task["id"], returncode=returncode, result=status,
           wall_s=time.monotonic()-started, t_end=finished_time, absorbed=absorbed)
    if success:
        atomic_json(directory/"done.json", dict(id=task["id"], t_end=finished_time,
                    absorbed=absorbed, utc=utc(), job=job, binary_sha256=binary_hash,
                    final_sha256=digest(directory/"final.bin")))
    else:
        # Release failed/preempted tasks for a later explicit worker invocation.
        # Current workers do not rescan, preventing immediate retry loops.
        (claim/"owner.json").unlink(missing_ok=True)
        claim.rmdir()
        exit_code = exit_code or 2
        if returncode != 0 and status == "finished":
            record("simulation_failure_stop", id=task["id"])
            stopping = True
    if stopping:
        break
record("worker_end", exit_code=exit_code, monitor_version=2,
       poststartup_samples=job_util_count,
       poststartup_mean_util=job_util_sum/job_util_count if job_util_count else None)
sys.exit(exit_code)
