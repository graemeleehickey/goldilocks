# RMST implementation validation

Run on 2026-09-08 with R 4.5.3, survival 3.8.6, and survRM2 1.0.4 on
arm64 macOS. Results describe these example designs and software versions.

## Statistical and integration checks

- Hand calculations and survRM2 comparisons cover RMST differences, Greenwood
  variances, and P-values at multiple fixed restriction times, including ties,
  events at zero, unequal allocation, and censoring.
- Completed-outcome calculations agree with Kaplan-Meier integration. Tests
  cover unsupported positive tails, curves reaching zero early, zero total
  variance, nonzero margins, arm reversal, and time-unit transformations.
- Integration checks cover observed final censoring, Rubin pooling, both
  completed-data routes, current and maximum predictive completions, saved
  horizons, caller RNG preservation, and serial/PSOCK reproducibility.
- Fixed-design simulation tests cover equal-survival, equal-RMST crossing-curve,
  and nonzero-margin nulls under arm-specific independent dropout, plus
  confidence-interval coverage. Acceptance limits account for Monte Carlo
  uncertainty; passing them does not establish exact finite-sample calibration.
- The full local test suite passed. Its installed-namespace-only test was
  skipped locally and exercised by the installed-package R CMD check.

`R CMD check --no-manual` finished with **0 errors, 0 warnings, and 0 notes**,
including vignette creation and rebuilding. Its installed-package test run
reported 1,412 passing expectations, no failures or warnings, and eight tests
skipped under CRAN settings. Longer calibration tests and the RMST PSOCK test
were also run outside CRAN settings and passed.

After grouping the source into `R/analysis_*.R`, a comparison of parsed R
expressions confirmed that all 174 top-level expressions were unchanged.
The full local suite then passed 1,432 expectations with no failures or
warnings; only the installed-namespace-only test was skipped locally.

## Calibration scenarios

Reproduce the main grid from the repository root:

```sh
Rscript benchmarks/rmst-calibration.R 500 /tmp/rmst-calibration.csv
```

Each of the 20 scenarios used 500 trials, 200 maximum participants, RMST
through month 9, follow-up through month 12, a one-sided final threshold of
0.975, and 100 predictive/final imputations where applicable. Adaptive designs
looked at enrollment 100, with Sn = 0.9, Fn = 0.05, and Qn = 1 (immediate
success disabled). Seeds 1430–1449 identify the rows in the CSV. Control
hazard was 0.1 throughout; the generation cutpoint was month 3.

The crossing-curve treatment hazards were selected analytically to give equal
RMST through month 9. The nonzero-null treatment hazard was 0.12 with its true
RMST difference (-0.4309328 months) as h0. Delayed benefit changed treatment
hazard from 0.1 to 0.04 after month 3 (true RMST difference 0.6092636 months).
Dropout variants used independent dropout CDFs of 0.1 in control and 0.2 in
treatment at month 12. The misspecified variant used a constant-hazard
predictive model despite the piecewise data generation.

No analysis failed in the 10,000 grid trials. Rates below are conditional on
successful evaluation; the CSV retains requested/completed/failed counts,
exact 95% binomial intervals, Monte Carlo standard errors, and mean enrollment.

| Scenario | Variant | Success rate | 95% Monte Carlo interval | Mean N |
|---|---|---:|---:|---:|
| equal_survival | fixed | 4.4% | 2.78%–6.59% | 200.0 |
| equal_survival | adaptive | 3.6% | 2.15%–5.63% | 151.8 |
| equal_survival | adaptive_dropout | 3.8% | 2.30%–5.87% | 151.4 |
| equal_survival | adaptive_imputed | 2.6% | 1.39%–4.41% | 155.8 |
| equal_survival | adaptive_misspecified | 2.4% | 1.25%–4.15% | 149.8 |
| equal_rmst_crossing | fixed | 2.6% | 1.39%–4.41% | 200.0 |
| equal_rmst_crossing | adaptive | 3.2% | 1.84%–5.14% | 152.8 |
| equal_rmst_crossing | adaptive_dropout | 2.8% | 1.54%–4.65% | 155.4 |
| equal_rmst_crossing | adaptive_imputed | 2.0% | 0.96%–3.65% | 159.6 |
| equal_rmst_crossing | adaptive_misspecified | 2.8% | 1.54%–4.65% | 153.0 |
| nonzero_null | fixed | 2.6% | 1.39%–4.41% | 200.0 |
| nonzero_null | adaptive | 2.6% | 1.39%–4.41% | 150.6 |
| nonzero_null | adaptive_dropout | 2.6% | 1.39%–4.41% | 150.2 |
| nonzero_null | adaptive_imputed | 2.6% | 1.39%–4.41% | 149.4 |
| nonzero_null | adaptive_misspecified | 3.4% | 1.99%–5.39% | 152.8 |
| delayed_benefit | fixed | 28.4% | 24.49%–32.57% | 200.0 |
| delayed_benefit | adaptive | 28.6% | 24.68%–32.78% | 178.0 |
| delayed_benefit | adaptive_dropout | 26.2% | 22.40%–30.29% | 175.8 |
| delayed_benefit | adaptive_imputed | 25.6% | 21.83%–29.66% | 175.0 |
| delayed_benefit | adaptive_misspecified | 29.8% | 25.82%–34.02% | 176.0 |

The initial fixed equal-survival batch rejected 22/500 trials (4.4%), above
the nominal 2.5% level. An independent, larger run with the same fixed design
(seed 1490, 10,000 trials, N_impute = 1 because no imputation was needed)
rejected 275/10,000 trials: 2.75%, with an exact 95% Monte Carlo interval of
2.438–3.090%. No analysis failed. Both runs are retained; the larger run's
interval includes the nominal target. This does not establish uniform error
control for all sample sizes, censoring patterns, horizons, or adaptive rules.

The grid is a characterization of the example thresholds, not a selected and
independently validated adaptive design. Choose and validate thresholds for
the intended design, including any immediate-success rule, using sufficiently
large independent simulations. Misspecified imputation and sparse events
require particular attention. Software tests do not replace that calibration.

## Runtime

`Rscript benchmarks/rmst.R 5` gave the following illustrative medians during
this validation session:

| Calculation | Median | Memory allocated |
|---|---:|---:|
| Completed RMST, automatic | 70.6 microseconds | 31.2 KB |
| Completed RMST, Kaplan-Meier reference | 1.09 milliseconds | 341.4 KB |
| Censored RMST, Kaplan-Meier | 1.26 milliseconds | 338.4 KB |
| Adaptive RMST trial | 45.4 milliseconds | 20.1 MB |
| Adaptive Cox trial | 49.6 milliseconds | 23.6 MB |

The completed-outcome calculation was about 15 times faster in this example.
Timings are hardware- and workload-dependent; benchmark iterations are small
and other validation jobs were running. The trial rows use corresponding
benefit directions and each method's own stopping decisions.
