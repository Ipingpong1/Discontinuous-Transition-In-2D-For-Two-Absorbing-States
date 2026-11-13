"""Regenerate the v2 theory figure from the recorded sensitivity tables."""
from pathlib import Path
import csv
import sys

ROOT = Path(__file__).resolve().parents[3]
ANALYSIS = ROOT / "campaign" / "v2" / "theory"
sys.path.insert(0, str(ANALYSIS))
from cutoff_sensitivity import make_figures


def read_table(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key, value in row.items():
            try:
                row[key] = float(value)
            except ValueError:
                pass
    return rows


if __name__ == "__main__":
    make_figures(read_table(ANALYSIS / "results" / "cutoff_thresholds.csv"),
                 read_table(ANALYSIS / "results" / "matching_sensitivity.csv"),
                 Path(__file__).resolve().parent)
