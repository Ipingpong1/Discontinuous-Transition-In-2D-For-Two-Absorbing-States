"""Plot recorded closure sensitivity, marking historical endpoint measurements."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
from common import save

ROOT = Path(__file__).resolve().parents[3]
TABLES = ROOT / "campaign/v2/theory/results"
SIMULATION = [(0.9,0.01194),(0.95,0.00619),(0.98,0.00252),
              (0.99,0.00128),(0.995,0.00064)]


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


def make_figures(rows, matched, destination):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":8, "axes.spines.top":False,
                         "axes.spines.right":False, "pdf.fonttype":42})
    sync = [r for r in rows if r["rate"] == 1 and r["T_steps"] >= 16]
    times = np.array([r["T_steps"] for r in sync])
    fig = plt.figure(figsize=(7.1, 5.0))
    grid = fig.add_gridspec(2, 2)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])]
    ax = axes[0]
    for key, label, color in [
        ("absorbing",r"$c_a$: absorbing stability","C3"),
        ("front",r"$c_f$: stationary front","C0"),
        ("active",r"$c_h$: active stability","C2")]:
        ax.errorbar(times, [r["c_"+key] for r in sync],
                    yerr=[1.96*r["se_"+key] for r in sync],
                    fmt="o-", ms=2.5, lw=.9, color=color, label=label)
    ax.set(xscale="log", xlabel=r"Dual-walk cutoff $T$",
           ylabel=r"Fixed-cutoff threshold $c=\kappa/\alpha$")
    ax.legend(fontsize=6.6)
    ax.text(0,1.02,"(a)",transform=ax.transAxes,va="bottom")
    ax = axes[1]
    for label, pretty, color in [
        ("absorbing_branch",r"$c_a-c_f$","C3"),
        ("active_branch",r"$c_f-c_h$","C2")]:
        ax.errorbar(times,[r[label] for r in sync],
                    yerr=[1.96*r["se_"+label] for r in sync],
                    fmt="o-",ms=2.5,lw=.9,color=color,label=pretty)
    ax.set(xscale="log",yscale="log",xlabel=r"Dual-walk cutoff $T$",
           ylabel="Closure branch width")
    ax.legend(fontsize=7)
    ax.text(0,1.02,"(b)",transform=ax.transAxes,va="bottom")
    ax = axes[2]
    points=[r for r in matched if r["A"]==1]
    for prescription,label,color in [
        ("inverse",r"$T=1/\varepsilon$","C1"),
        ("inverse_log3",r"$T=\log^3(1/\varepsilon)/\varepsilon$","C0")]:
        data=[r for r in points if r["prescription"]==prescription]
        x=np.array([r["epsilon"] for r in data])
        y=np.array([r["c_front"] for r in data])
        ax.plot(x,y,"o-",ms=3,color=color,label=label)
        low=[next(r for r in matched if r["p"]==d["p"] and r["prescription"]==prescription and r["A"]==.25)["c_front"] for d in data]
        high=[next(r for r in matched if r["p"]==d["p"] and r["prescription"]==prescription and r["A"]==4)["c_front"] for d in data]
        # Missing coverage stays missing; no extrapolated upper envelope is drawn.
        ax.fill_between(x,low,high,color=color,alpha=.12)
    survey = SIMULATION[:3]
    historical = SIMULATION[3:]
    ax.plot([1-p for p,q in survey], [q/(.5*(1-p)) for p,q in survey],
            "kD", mfc="white", ms=4, label=r"Survey, $p\leq0.98$")
    ax.plot([1-p for p,q in historical], [q/(.5*(1-p)) for p,q in historical],
            "k^", mfc="white", ms=5, label=r"Earlier estimates, $p=0.99,0.995$")
    # The original q spread is a finite-window sensitivity, without CI coverage.
    ax.errorbar([.005], [.256], yerr=[.002], fmt="none", color="0.35",
                capsize=4, label=r"Historical spread, $p=0.995$")
    for p, q in historical:
        ax.annotate(f"{p:.3g}", (1-p,q/(.5*(1-p))), xytext=(0,9),
                    textcoords="offset points", ha="center", fontsize=6.5)

    ax.set(xscale="log",xlabel=r"$\varepsilon=1-p$",ylabel=r"$c_f$ and measured $q^*/[\alpha(1-p)]$")
    ax.invert_xaxis()
    ax.legend(fontsize=6.5,loc="lower right", ncol=2)
    ax.text(0,1.02,"(c)",transform=ax.transAxes,va="bottom")
    fig.tight_layout()
    save(fig, "fig_theory_v2")
    plt.close(fig)


if __name__ == "__main__":
    sources = [TABLES/"cutoff_thresholds.csv", TABLES/"matching_sensitivity.csv"]
    make_figures(*(read_table(path) for path in sources), Path(__file__).resolve().parent)
    output = ROOT / "campaign/v3/results/theory"
    output.mkdir(parents=True, exist_ok=True)
    metadata = {
        "inputs": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in sources},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "simulation": SIMULATION,
        "historical": [0.99, 0.995],
        "p0995_spread": {"q_half_spread": 0.000005, "c_half_spread": 0.002,
                         "interpretation": "historical estimator spread, not a confidence interval"},
    }
    (output/"figure_provenance.json").write_text(json.dumps(metadata,indent=2)+"\n")
