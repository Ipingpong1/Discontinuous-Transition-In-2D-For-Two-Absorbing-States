#!/usr/bin/env python3
"""Generate the six-worker production campaign; no simulation or submission."""
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--out", type=Path, default=HERE / "production_manifest.json")
parser.add_argument("--source", type=Path, default=HERE / "../../../round2/dnsim.c")
args = parser.parse_args()

def seed(text):
    return int.from_bytes(hashlib.blake2b(("dn-v2-20261103:" + text).encode(), digest_size=8).digest(), "little")

tasks = []
def add(family, Lx, Ly, q, replica, T, every, p, profile=False, orientation=None):
    suffix = f"_ab{orientation[0]}-{orientation[1]}" if orientation else ""
    name = f"{family}_L{Lx}x{Ly}{suffix}_q{q}_s{replica}"
    mode = "band" if orientation else "strip" if family == "near_strip" else "bulk"
    queue = "near_voter" if family.startswith("near_") else family
    tasks.append(dict(id=name, family=queue, observable=family, Lx=Lx, Ly=Ly,
                      q=q, p=p, alpha="0.5", seed=seed(name + ":dynamics"),
                      init_seed=seed(name + ":initial"), replica=replica,
                      T=T, every=every, profile=profile, orientation=orientation, mode=mode,
                      checkpoint_every=1_000_000, work=Lx * Ly * T))

for L, nseed in [(4096, 8), (8192, 4)]:
    for ab in [(1,0), (2,1), (1,1)]:
        for q in [".04588", ".04590", ".04592", ".04594", ".04596"]:
            for replica in range(1, nseed + 1):
                add("orientation", L, L, q, replica, 2_000_000, 5000, "0.5", True, ab)
for L in [2048, 4096]:
    for q in [".000632", ".000636", ".000640", ".000644", ".000648"]:
        for replica in range(1, 5):
            add("near_bulk", L, L, q, replica, 10_000_000, 10000, "0.995")
for q in [".000632", ".000636", ".000640", ".000644", ".000648"]:
    for replica in range(1, 5):
        add("near_strip", 2048, 8192, q, replica, 10_000_000, 20000, "0.995", True)
for L in [2048, 4096, 8192]:
    for q in [".04590", ".04592", ".04594", ".04596", ".04600"]:
        for replica in range(1, 5):
            add("flagship_bulk", L, L, q, replica, 5_000_000, 10000, "0.5")

# Largest jobs first limits end-of-campaign imbalance; deterministic hash order
# interleaves q/orientation/replica combinations within each workload size.
tasks.sort(key=lambda t: (t["family"], -t["work"], seed(t["id"] + ":order")))
seeds = [t[key] for t in tasks for key in ("seed", "init_seed")]
assert len(seeds) == len(set(seeds)), "seed collision"
assert len({t["id"] for t in tasks}) == len(tasks)
manifest = dict(schema=1, campaign="dn-v2-20261103-gpu",
                deadline_utc="2026-11-03T11:00:00+00:00",
                source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                seed_algorithm="BLAKE2b-64; namespace dn-v2-20261103",
                families={family: dict(workers=2, tasks=sum(t["family"] == family for t in tasks))
                          for family in ("orientation", "near_voter", "flagship_bulk")},
                tasks=tasks)
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(manifest, indent=2) + "\n")
print(args.out)
for family, summary in manifest["families"].items():
    print(family, summary["tasks"], "runs",
          sum(t["work"] for t in tasks if t["family"] == family), "site updates")
