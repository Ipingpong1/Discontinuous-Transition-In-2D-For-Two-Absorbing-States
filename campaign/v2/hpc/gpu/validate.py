#!/usr/bin/env python3
"""Compare dense CPU and CUDA runs, including binary output and resumed runs."""
from pathlib import Path
import argparse
import os
import struct
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("--cpu", type=Path, required=True)
parser.add_argument("--gpu", type=Path, required=True)
parser.add_argument("--work", type=Path, help="keep validation artifacts in this directory")
parser.add_argument("--quick", action="store_true", help="skip the L=2048 workload")
args = parser.parse_args()
cpu, gpu = args.cpu.resolve(), args.gpu.resolve()
temporary = tempfile.TemporaryDirectory(prefix="dnsim-gpu-validation-") if args.work is None else None
work = Path(temporary.name) if temporary else args.work.resolve()
work.mkdir(parents=True, exist_ok=True)
environment = dict(os.environ, OMP_NUM_THREADS="1")

def run(binary, options, directory):
    directory.mkdir(parents=True, exist_ok=True)
    command = [str(binary)] + options + [
        "--out", str(directory / "out.txt"),
        "--profile", str(directory / "profile.bin"),
        "--snap", str(directory / "snap.bin"),
        "--save", str(directory / "save.bin"),
        "--ckpt", str(directory / "checkpoint.bin")]
    with (directory / "stderr.txt").open("a") as stderr:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=stderr, env=environment)

def rows(path):
    return [line for line in path.read_text().splitlines() if line and not line.startswith("#")]

def compare(name):
    left, right = work / name / "cpu", work / name / "gpu"
    assert rows(left / "out.txt") == rows(right / "out.txt"), f"{name}: time-series mismatch"
    for filename in ("profile.bin", "snap.bin", "save.bin", "checkpoint.bin"):
        assert (left / filename).read_bytes() == (right / filename).read_bytes(), f"{name}: {filename} mismatch"
    print("PASS", name, flush=True)

def pair(name, options):
    run(cpu, options, work / name / "cpu")
    run(gpu, options, work / name / "gpu")
    compare(name)

cases = [
    ("large_graph", ["--L", "2048", "--T", "1025", "--seed", "617",
                     "--alpha", ".5", "--voter", ".5", "--jflip", "2:.046",
                     "--ppd", "0", "--every", "257", "--snap-every", "9999", "--ckpt-every", "9999"]),
    ("lazy_table", ["--Lx", "128", "--Ly", "13", "--T", "1031", "--seed", "147",
                    "--ftab", "0,.0078125,.021375,.0546875,.03125,.0390625,.078125,.0859375,.09375;.09375,.0859375,.078125,.0390625,.03125,.0546875,.021375,.0078125,0",
                    "--ppd", "0", "--every", "515", "--snap-every", "9999", "--ckpt-every", "9999"]),
    ("moore_w1", ["--Lx", "64", "--Ly", "7", "--T", "37", "--seed", "18446744073709551615",
                  "--alpha", ".5", "--voter", ".5", "--jflip", "2:.04596", "--ppd", "4",
                  "--every", "5", "--snap-every", "7", "--snap-from", "3", "--ckpt-every", "11"]),
    ("graph_parity", ["--Lx", "128", "--Ly", "13", "--T", "1025", "--seed", "17",
                      "--alpha", ".5", "--voter", ".5", "--jflip", "2:.046",
                      "--init", "strip:.5:.5", "--ppd", "0", "--every", "257",
                      "--snap-every", "9999", "--ckpt-every", "9999"]),
    ("molofsky", ["--Lx", "192", "--Ly", "17", "--T", "769", "--seed", "4294967303",
                  "--molofsky", ".19,.41", "--init", "hole:.5:4", "--ppd", "0",
                  "--every", "513", "--snap-every", "9999", "--ckpt-every", "9999"]),
    ("near_voter", ["--Lx", "128", "--Ly", "11", "--T", "530", "--seed", "298",
                    "--alpha", ".5", "--voter", ".995", "--jflip", "2:.00064",
                    "--ppd", "3", "--every", "0", "--snap-every", "43", "--ckpt-every", "97"]),
    ("direction", ["--L", "128", "--T", "19", "--seed", "8", "--alpha", ".5",
                   "--voter", ".5", "--jflip", "2:.046", "--bg", "1", "--init", "disk:.3:21",
                   "--pdir", "1,2", "--ppd", "0", "--every", "3",
                   "--snap-every", "7", "--ckpt-every", "5"]),
    ("all_classes", ["--Lx", "128", "--Ly", "9", "--T", "521", "--seed", "49",
                     "--ftab", "0,.00000000023283064365386962890625,.13,.24,.35,.46,.57,.68,1;.91,.82,.73,.64,.55,.46,.37,.28,.19",
                     "--ppd", "0", "--every", "259", "--snap-every", "9999", "--ckpt-every", "9999"]),
    ("deterministic", ["--Lx", "64", "--Ly", "3", "--T", "520", "--seed", "61",
                       "--ftab", "0,1,0,1,0,1,0,1,0;1,0,1,0,1,0,1,0,1",
                       "--ppd", "0", "--every", "257", "--snap-every", "9999", "--ckpt-every", "9999"]),
]
for name, options in cases:
    if args.quick and name == "large_graph":
        continue
    pair(name, options)

resume_options = ["--Lx", "128", "--Ly", "13", "--alpha", ".5", "--voter", ".5",
                  "--jflip", "2:.046", "--seed", "391", "--ppd", "5", "--every", "61",
                  "--snap-every", "67", "--ckpt-every", "73"]
pair("resume", resume_options + ["--T", "333"])
cpu_start = work / "cpu_333.bin"
gpu_start = work / "gpu_333.bin"
cpu_start.write_bytes((work / "resume/cpu/checkpoint.bin").read_bytes())
gpu_start.write_bytes((work / "resume/gpu/checkpoint.bin").read_bytes())
pair("resume", resume_options + ["--T", "901", "--resume"])
run(cpu, resume_options + ["--T", "901", "--load", str(gpu_start)], work / "cross_restart/cpu")
run(gpu, resume_options + ["--T", "901", "--load", str(cpu_start)], work / "cross_restart/gpu")
compare("cross_restart")

# Exercise the last valid 32-bit RNG counters without billions of steps.
initial = bytearray((work / "graph_parity/cpu/save.bin").read_bytes())
struct.pack_into("=Q", initial, 24, 2**32 - 521)
load_path = work / "high_time.bin"
load_path.write_bytes(initial)
pair("high_time", ["--Lx", "128", "--Ly", "13", "--T", str(2**32 - 1),
                  "--load", str(load_path), "--seed", "17", "--alpha", ".5", "--voter", ".5",
                  "--jflip", "2:.046", "--ppd", "0", "--every", "0",
                  "--snap-every", "9999999999", "--ckpt-every", "0"])
result = subprocess.run([str(gpu), "--reps", "1"], env=environment, capture_output=True)
assert result.returncode == 2 and b"does not support --reps" in result.stderr, "GPU --reps rejection missing"
print("PASS --reps rejection")
print("All CPU/backend comparisons passed.", flush=True)
if temporary:
    temporary.cleanup()
