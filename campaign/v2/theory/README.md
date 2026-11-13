# Finite-cutoff theory audit for v2

Run from any working directory:

    python3 -B campaign/v2/theory/cutoff_sensitivity.py --selftest
    python3 -B campaign/v2/theory/cutoff_sensitivity.py
    python3 -B paper/v2/figs/fig_theory_v2.py

Requires NumPy and Matplotlib; no SciPy or new walk simulations are needed.
The original histograms, original manuscript, and original plotting code are not changed.

Outputs:
- results/cutoff_thresholds.csv: fixed-cutoff absorbing, front, and active stability thresholds; their pointwise errors and covariance-aware branch widths; truncated-lineage identities and root checks.
- results/matching_sensitivity.csv: declared bare and logarithmic matching prescriptions with A=1/4,1,4. No extrapolation.
- results/provenance.json: source hashes, code hash, bootstrap settings, and statistical limitations.
- paper/v2/figs/fig_theory_v2.pdf and .png.
- paper/v2/theory_revision.tex: replacement section for integration into v2.

## Statistical interpretation

Each cutoff histogram is multinomial over lineage patterns plus the omitted zero-drift category (at most two surviving lineages). For a ratio c=-A/B, per-trial influence is -(a_i+c*b_i)/B. Its multinomial covariance gives the joint covariance of the three threshold estimates. Bootstrap resampling uses the same resampled counts for every threshold. This retains the strong cancellation in differences.

Error bars are pointwise Monte Carlo sampling errors only. They are neither closure error nor cutoff error. At successive cutoffs the same random-walk paths were observed; the stored marginal histograms do not retain cross-cutoff covariance. No asymptotic fit or confidence interval for a limiting amplitude is attempted.

Matching interpolation is linear in log time. Its reported Monte Carlo SE is an upper bound obtained by adding endpoint SEs with the interpolation weights, valid even when their correlation is unknown. This bound does not include interpolation bias.

The A=1/4,1,4 grid is a declared sensitivity range, not a prior, fit, or confidence interval. The log-cubed prescription is motivated by the known two-dimensional scaling in Cox–Perkins Eq. (1.26), but is not derived for the present synchronous automaton. Its improved agreement should not be interpreted as validating an asymptotic matching prescription, especially near cancellation of the leading cubic drift.

Simulation comparisons use Table I central values from paper/main.tex; this script does not re-estimate their errors. Parent v2 analyses may supersede them.

## Exact finite-time mechanism

Only partitions with at least five lineages generate departures from the cubic h(u)=u(1-u)(1-2u). If A_X=f_X'(0), I_X=integral_0^1/2 f_X, and B_X=I_X-A_X/32, then

    c_a-c_f = 32*(A_2*B_R-A_R*B_2) / [A_2*(A_2+32*B_2)].

The B terms receive contributions only from >=5-lineage partitions. This identity is exact at fixed cutoff. Its use as a physical first-order mechanism is conditional on the finite-time voter closure and deterministic front approximation.

## Further simulation priorities

1. Actual CA tests of phase coexistence and the continuous-time limit take priority over reducing already small fixed-cutoff Monte Carlo errors.
2. If CPU time remains, an independent synchronous nine-walk sample extending to T=2^20 would cover the missing high-cutoff sensitivity points.
3. Save per-trial block histograms at every checkpoint in future walk runs. This would enable a joint bootstrap over time without treating correlated checkpoints as independent.
4. A closure validation experiment could directly measure conditional local drift in perturbed bulk/strip configurations, compare its shape with the finite-T reference family, and test whether one cutoff can reproduce more than one observable.
