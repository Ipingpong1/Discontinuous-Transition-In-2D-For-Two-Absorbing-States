#!/usr/bin/env python3
"""Finite-cutoff voter closure: pointwise Monte Carlo errors and matching sensitivity.

The input histograms contain marginal checkpoint counts from common walker trials.
Their within-checkpoint covariance is retained. Cross-checkpoint covariance is not
recoverable: no independent-time fit, extrapolated limit, or limit confidence
interval is produced. Matching curves are sensitivity diagnostics, not a derived
two-dimensional hydrodynamic limit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial import polynomial as P

ROOT = (Path.cwd() if __file__ == "<stdin>" else Path(__file__).resolve().parents[3])
THEORY = ROOT / "campaign" / "theory"
sys.path.insert(0, str(THEORY))
import vpt

SOURCES = [
    ("synchronous", 1.0, THEORY / "patterns.txt"),
    ("lazy_0.5", 0.5, THEORY / "patterns_r0.5.txt"),
    ("lazy_0.25", 0.25, THEORY / "patterns_r0.25.txt"),
    ("lazy_0.125", 0.125, THEORY / "patterns_r0.125.txt"),
]
# Published central values, not independently re-estimated by this script.
SIMULATION = [(0.9, 0.01194), (0.95, 0.00619), (0.98, 0.00252),
              (0.99, 0.00128), (0.995, 0.00064)]
KINDS = ["absorbing", "front", "active"]
CONTRASTS = {
    "absorbing_branch": np.array([1., -1., 0.]),
    "active_branch": np.array([0., 1., -1.]),
    "full_window": np.array([1., 0., -1.]),
}


def components(ns):
    """Per-pattern (R,delta_2) values for f'(0), integral, and f'(1/2)."""
    F = vpt.pattern_poly(ns)
    r, q = vpt.R @ F, vpt.D2 @ F
    return np.array([r[1], q[1],
                     P.polyval(.5, P.polyint(r)), P.polyval(.5, P.polyint(q)),
                     P.polyval(.5, P.polyder(r)), P.polyval(.5, P.polyder(q))])


def ratio_statistics(probabilities, values, trials):
    mu = probabilities @ values
    ratios = -mu[::2] / mu[1::2]
    # Per-trial influence for -A/B. The unreported <=2-class category is zero.
    influence = -(values[:, ::2] + ratios * values[:, 1::2]) / mu[1::2]
    cov = (influence.T * probabilities) @ influence / trials
    return ratios, cov, mu


def summarize(name, rate, T, trials, patterns, rng, nboot):
    counts = np.array([c for c, _ in patterns] +
                      [trials - sum(c for c, _ in patterns)], dtype=np.int64)
    if np.any(counts < 0) or counts.sum() != trials:
        raise ValueError("Invalid histogram counts")
    values = np.vstack([components(ns) for _, ns in patterns] + [np.zeros(6)])
    probs = counts / trials
    central, cov, mu = ratio_statistics(probs, values, trials)
    sampled = rng.multinomial(trials, probs, size=nboot) @ values / trials
    boot = -sampled[:, ::2] / sampled[:, 1::2]
    if not np.all(np.isfinite(boot)):
        raise ValueError("Bootstrap denominator unresolved; do not quote a ratio interval")
    n3plus = sum(c for c, ns in patterns if len(ns) >= 2)
    n5plus = sum(c for c, ns in patterns if len(ns) >= 4)
    row = dict(dataset=name, rate=rate, T_steps=T, T_physical=rate*T,
               trials=trials, n_ge3=n3plus, n_ge5=n5plus,
               p_ge3=n3plus/trials, ge5_over_ge3=n5plus/n3plus)
    for i, kind in enumerate(KINDS):
        lo, hi = np.quantile(boot[:, i], [.025, .975])
        row.update({f"c_{kind}":central[i],
                    f"se_{kind}":math.sqrt(max(0., cov[i, i])),
                    f"ci025_{kind}":lo, f"ci975_{kind}":hi})
    for label, contrast in CONTRASTS.items():
        lo, hi = np.quantile(boot @ contrast, [.025, .975])
        row.update({label:central @ contrast,
                    "se_"+label:math.sqrt(max(0., contrast @ cov @ contrast)),
                    "ci025_"+label:lo, "ci975_"+label:hi})
    row["corr_absorbing_front"] = cov[0, 1] / math.sqrt(cov[0, 0]*cov[1, 1])
    AR, A2, IR, I2 = mu[:4]
    BR, B2 = IR-AR/32, I2-A2/32
    identity = 32*(A2*BR-AR*B2)/(A2*(A2+32*B2))
    row["lineage_identity_residual"] = identity-row["absorbing_branch"]
    low = [(c, ns) for c, ns in patterns if len(ns) <= 3]
    F_low = vpt.drift_basis(low, trials)
    lowR, low2 = vpt.R @ F_low, vpt.D2 @ F_low
    if low2[1] != 0:
        row["c_le4_lineages"] = -lowR[1]/low2[1]
        row["le4_front_minus_absorbing"] = (
            -P.polyval(.5, P.polyint(lowR))/P.polyval(.5, P.polyint(low2))
            -row["c_le4_lineages"])
    else:
        row["c_le4_lineages"] = float("nan")
        row["le4_front_minus_absorbing"] = float("nan")
    full = vpt.drift_basis(patterns, trials)
    drift = vpt.R @ full + central[1]*(vpt.D2 @ full)
    row["front_slope_origin"] = drift[1]
    row["front_slope_half"] = P.polyval(.5, P.polyder(drift))
    roots = P.polyroots(P.polytrim(drift, 1e-14))
    interior = [x.real for x in roots if abs(x.imag) < 1e-7 and 1e-7 < x.real < .5-1e-7]
    row["front_interior_roots"] = len(interior)
    row["front_unstable_root"] = interior[0] if len(interior) == 1 else float("nan")
    return row


def interpolate(rows, time, key):
    """Linear interpolation in log time; no extrapolation.

    Endpoint errors may be correlated. Sum of interpolation weights times their
    SEs is a conservative upper bound on interpolant SE (Cauchy-Schwarz).
    It says nothing about closure or interpolation bias.
    """
    rows = sorted(rows, key=lambda r:r["T_physical"])
    times = np.array([r["T_physical"] for r in rows])
    if time < times[0] or time > times[-1]:
        return float("nan"), float("nan"), float("nan"), float("nan")
    k = int(np.searchsorted(times, time))
    if k == 0 or (k < len(times) and times[k] == time):
        r = rows[k]
        return r["c_"+key], r["se_"+key], times[k], times[k]
    left, right = rows[k-1], rows[k]
    weight = math.log(time/times[k-1])/math.log(times[k]/times[k-1])
    value = (1-weight)*left["c_"+key] + weight*right["c_"+key]
    bound = (1-weight)*left["se_"+key] + weight*right["se_"+key]
    return value, bound, times[k-1], times[k]


def matching_sensitivity(rows):
    synchronous = [r for r in rows if r["rate"] == 1]
    out = []
    for p, qstar in SIMULATION:
        eps = 1-p
        measured = qstar/(.5*eps)
        for prescription in ["inverse", "inverse_log3"]:
            base = 1/eps
            if prescription == "inverse_log3":
                base *= math.log(1/eps)**3
            for A in [.25, 1., 4.]:
                cutoff = A*base
                value, error, tlo, thi = interpolate(synchronous, cutoff, "front")
                out.append(dict(p=p, epsilon=eps, alpha=.5, qstar_source=qstar,
                                c_simulation=measured, prescription=prescription,
                                A=A, T_physical=cutoff, c_front=value,
                                mc_se_upper_bound=error, lower_checkpoint=tlo,
                                upper_checkpoint=thi, available=math.isfinite(value),
                                relative_residual=(value-measured)/measured))
    return out


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


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
    ax.plot([1-p for p,q in SIMULATION],[q/(.5*(1-p)) for p,q in SIMULATION],
            "kD",mfc="white",ms=4,label="Simulation, Table I")
    ax.set(xscale="log",xlabel=r"$\varepsilon=1-p$",ylabel=r"$c_f$ and measured $q^*/[\alpha(1-p)]$")
    ax.invert_xaxis()
    ax.legend(fontsize=6.5,loc="lower right")
    ax.text(0,1.02,"(c)",transform=ax.transAxes,va="bottom")
    fig.tight_layout()
    for extension in ["pdf","png"]:
        fig.savefig(destination / ("fig_theory_v2."+extension),dpi=220)
    plt.close(fig)


def selftest():
    # Fixed-cutoff identities, independently evaluated derivative/integral functionals,
    # and perfectly correlated artificial ratio contrasts.
    for ns in [(1,1),(1,2),(1,1,1),(2,3,1)]:
        value=components(ns)
        assert np.allclose(value[2:4],value[:2]/32,atol=1e-12)
    probs=np.array([.2,.3,.5])
    values=np.array([[1.,2.,2.,4.,3.,6.],[-1.,1.,-2.,2.,-3.,3.],[0.,0.,0.,0.,0.,0.]])
    ratios,cov,_=ratio_statistics(probs,values,1000)
    assert np.allclose(ratios,ratios[0])
    assert abs(np.array([1.,-1.,0.])@cov@np.array([1.,-1.,0.])) < 1e-15
    for T in [64,1024,65536]:
        N,patterns=vpt.read_patterns(THEORY/"patterns.txt")[T]
        row=summarize("test",1.,T,N,patterns,np.random.default_rng(1),100)
        assert abs(row["lineage_identity_residual"]) < 1e-12
        assert abs(row["le4_front_minus_absorbing"]) < 1e-12
        assert row["c_active"] < row["c_front"] < row["c_absorbing"]
        assert row["front_slope_origin"] < 0 and row["front_slope_half"] < 0
        assert row["front_interior_roots"] == 1
    print("cutoff_sensitivity selftest PASS")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selftest",action="store_true")
    parser.add_argument("--bootstrap",type=int,default=2000)
    parser.add_argument("--seed",type=int,default=20261102)
    parser.add_argument("--output",type=Path,default=ROOT/"campaign"/"v2"/"theory"/"results")
    parser.add_argument("--figures",type=Path,default=ROOT/"paper"/"v2"/"figs")
    parser.add_argument("--no-figures",action="store_true")
    args=parser.parse_args()
    if args.selftest:
        selftest()
        return
    if args.bootstrap < 100:
        parser.error("Use at least 100 bootstrap draws")
    rng=np.random.default_rng(args.seed)
    rows=[]
    sources=[]
    for name,rate,path in SOURCES:
        sources.append(dict(dataset=name,rate=rate,path=str(path.relative_to(ROOT)),
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        for T,(N,patterns) in sorted(vpt.read_patterns(path).items()):
            if T == 0:
                continue  # absorbing-state ratio is undefined at T=0
            rows.append(summarize(name,rate,T,N,patterns,rng,args.bootstrap))
    matched=matching_sensitivity(rows)
    args.output.mkdir(parents=True,exist_ok=True)
    write_csv(args.output/"cutoff_thresholds.csv",rows)
    write_csv(args.output/"matching_sensitivity.csv",matched)
    metadata=dict(seed=args.seed,bootstrap_draws=args.bootstrap,sources=sources,
                  uncertainty="Pointwise multinomial sampling only; no closure/cutoff error",
                  cross_time_covariance="Unavailable: common walk paths, marginal histograms",
                  interpolation="Linear in log time, with conservative MC SE bound; no extrapolation",
                  matching="Declared A=.25,1,4 sensitivity, not an inferred or validated cutoff",
                  simulation_source="paper/main.tex Table I central values; no simulation confidence interval implied",
                  code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/"provenance.json").write_text(json.dumps(metadata,indent=2)+"\n")
    if not args.no_figures:
        args.figures.mkdir(parents=True,exist_ok=True)
        make_figures(rows,matched,args.figures)
    print(f"Wrote {len(rows)} cutoff rows and {len(matched)} matching rows to {args.output}")
    for row in rows:
        if row["rate"]==1 and row["T_steps"] in [64,1024,65536]:
            print(f"T={row['T_steps']}: c_f={row['c_front']:.8f} ± {row['se_front']:.8f}; "
                  f"c_a-c_f={row['absorbing_branch']:.8f} ± {row['se_absorbing_branch']:.8f} (MC SE)")


if __name__ == "__main__":
    main()
