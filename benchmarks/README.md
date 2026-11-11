# Maintainer Benchmarks

This directory contains optional performance benchmarks for goldilocks hot
paths. They are intended for maintainers to compare branches before and after
optimization work, not as pass/fail tests.

The benchmarks are excluded from R package builds with `.Rbuildignore`, so they
do not run on CRAN or during ordinary package checks.

## Hot-Path Benchmarks

### Enrollment

The hot-path benchmark includes constant and fractional-knot enrollment
schedules. In particular, its low-rate case guards the continuous-time
generator against work that scales with empty calendar bins rather than with
the requested number of enrollments.

### Bayesian Survival Calculations

The benchmark reports piecewise-exponential posterior sampling from a
patient-level data frame and from precomputed sufficient statistics. The
sufficient-statistic pair compares the general route with the prepared-data
route used for summaries generated inside Bayesian predictive imputation.
Both receive the same
normalized prior and reset to the same seed, isolating validation overhead
without changing the Gamma posterior calculation.

The Bayesian-survival effect pair compares `haz_to_prob()`, which returns arm probabilities and effects, with direct calculation of the effects from cumulative-hazard draws. The latter reuses fixed analysis-interval widths and avoids repeating cutpoint validation or creating a probability data frame.

### Predictive Imputation

The predictive-imputation pair generates 100 expected-success and futility
completions for the same two-arm interim data. The scalar reference calls
`impute_data()` twice per posterior draw and retains each completed data frame;
the batch path stores only subject-by-draw time and event matrices for rows that
need imputation. Both expressions reset to the same seed and use the documented
draw, cohort, and arm order, so their generated outcomes can also be checked
directly in unit tests.

### Completed-Data Analyses

The survival completed-data pairs compare the retained patient-data reference
with the prepared-outcome route used at interim looks. The reference copies the
complete interim data frame for each current and maximum-sample analysis. The
prepared route records imputation positions once, copies only follow-up and
event values, and reuses treatment assignments. Separate pairs cover log-rank,
Cox, and Bayesian piecewise-exponential analyses. The script checks identical
results before timing them; for the Bayesian pair it also resets the seed so
both paths receive the same Gamma posterior draws.

The binary completed-data pair starts from one common predictive-imputation
batch. The patient-data reference constructs a completed data
frame for every current and maximum-sample analysis. The count path instead
uses arm-level event and subject totals and evaluates each distinct
deterministic count summary once. The script verifies identical classifications
before timing the pair and prints the number and percentage of reused summaries.
The separate `sim_trial_bayes_bin` row reports the complete trial runtime,
including generation, posterior prediction, interim decisions, and final
analysis.

### Risk-Difference Tests

The risk-difference rows apply the Wald and Farrington-Manning calculations to
the same 1,000 low-event-rate, non-boundary count tables with a nonzero null
difference. They make the cost of the constrained one-dimensional likelihood
calculation visible without mixing it with predictive imputation or trial
generation.

### Probability Transformations

The stable cumulative-hazard/probability transformations are benchmarked over
100,000 values spanning near-zero inputs through their mathematical
boundaries. This guards the use of `expm1()` and `log1p()` against a material
regression in simulation hot paths.

### Cox Model Fitting

The Cox benchmark exercises 1,000 subjects with tied event times through both
the package's guarded automatic path and its forced public `survival::coxph()`
fallback. The automatic path checks the installed `survival` version and
`coxph.fit()` signature once, validates each returned coefficient/variance
structure, and falls back to `coxph()` after any incompatibility. Access to the
unexported fitter is confined to that compatibility boundary.

The public formula interface can take roughly eight to twelve times as long as
the low-level fitter for one small Cox model. Because an adaptive trial may fit
hundreds or thousands of imputed Cox models, compare `cox_guarded_auto` with
`cox_public_fallback` before changing the compatibility boundary. Exact timing
depends on the installed `survival` version and hardware; the script prints the
detected fast-path status and reason with its results.

### Running

Install the suggested `bench` package, then run from the package root:

```r
source("benchmarks/hot-paths.R")
```

When `pkgload` is installed, the script loads the package source tree directly.
Otherwise, it falls back to the installed `goldilocks` package.

To save results, set `GOLDILOCKS_BENCHMARK_OUT` to a CSV path before sourcing
the script:

```r
Sys.setenv(GOLDILOCKS_BENCHMARK_OUT = "benchmarks/results.csv")
source("benchmarks/hot-paths.R")
```

## Parallel Backends

The cross-platform parallel benchmark compares automatic selection, serial,
PSOCK, and (on Unix-like platforms) the existing fork backend over
small-overhead, representative log-rank, and heavier Bayesian-survival
workloads. Its default trial counts of 2, 8, and 32 exercise the
serial/parallel crossover as well as medium and larger task sets.

The results record requested and effective workers, the backend-selection
reason, median total time, separately measured PSOCK cluster startup, estimated
compute/dispatch time, and serial-relative speedup, plus R version and
operating-system metadata. PSOCK setup remains included in total call time; the
compute estimate subtracts the median time to create and stop an equivalent
cluster. Parent-process allocation is reported for the serial baseline;
allocation profiling is unavailable for child-process rows.

```r
source("benchmarks/parallel-backends.R")
```

By default it evaluates 2, 8, and 32 trials with 1, 2, and 4 workers. Adjust
the workload with these environment variables:

```r
Sys.setenv(
  GOLDILOCKS_PARALLEL_ITERATIONS = "5",
  GOLDILOCKS_PARALLEL_TRIALS = "8,32",
  GOLDILOCKS_PARALLEL_WORKERS = "2,4",
  GOLDILOCKS_BENCHMARK_OUT = "benchmarks/parallel-results.csv"
)
source("benchmarks/parallel-backends.R")
```

The package deliberately retains the established Unix fork path. Replace it
only when this benchmark shows reproducible results, no material regression
against the fork baseline, and a consistent capability or speed advantage.
On Windows, compare `auto` with `psock`; on macOS, compare those rows with
`fork` to quantify the platform-specific PSOCK startup cost.

## Interpreting Results

Runtime and memory allocation depend on hardware, R version, operating system,
and BLAS configuration. Compare results relatively on the same machine rather
than treating any absolute timing as a release threshold.

## Bayesian binary quadrature

Author: Graeme L. Hickey.

`beta-quadrature.R` compares the current two-arm Beta probability calculation with two experimental alternatives: integration over treatment Beta quantiles, and the same transformation with intervals placed around changes in the control posterior CDF. Both alternatives use direct upper-tail probabilities and absolute and relative integration tolerances of `1e-10`. The split version uses a bounded-error midpoint calculation on negligible intervals to avoid numerical roundoff near probability zero or one. These alternatives are defined only in the benchmark; the package implementation is unchanged.

Run from the package root with `bench`, `pkgload`, and `Rmpfr` installed:

```sh
Rscript benchmarks/beta-quadrature.R 200 5 benchmarks/beta-quadrature-results.csv
```

The arguments specify the number of repetitions per probability calculation, repetitions per ten-trial batch, and output file. Sixteen count tables cover ordinary trial sizes, sparse outcomes, Jeffreys priors, concentrated posteriors, unequal arm sizes, and a nonzero risk-difference margin. Both one-sided alternatives are evaluated, giving 32 probability comparisons. References use symmetry, positive finite sums evaluated with 128-bit arithmetic, or separately partitioned density integration checked in both arm orientations. Accuracy is assessed separately from timing.

The saved run used R 4.5.3 on an Apple ARM64 Mac. Probability timings exclude iterations affected by garbage collection; full-trial timings include it. The first row below takes the median of the ten case-specific timings for the five ordinary count tables in both directions. Full-trial times are medians of five batches using the same simulated histories.

| Workload | Current calculation | Split quantile prototype | Runtime ratio |
|:--|--:|--:|--:|
| One probability calculation, ordinary count tables | 61.0 microseconds | 739.2 microseconds | 12.1 |
| Ten adaptive trials | 0.729 seconds | 4.075 seconds | 5.59 |

The trial workload has maximum enrollment of 500, interims at 100, 200, 300, and 400, a 12-month endpoint, 100 predictive imputations, equal randomization, and no dropout. Control and treatment event risks are 30% and 22%, with accrual averaging ten participants per month. It uses Gamma(0.5, 20) shape-rate hazard priors, uniform Beta priors, `Fn = 0.05`, `Sn = 0.9`, `Qn = 1`, `prob_ha = 0.95`, and seed 92502. Both the current calculation and the split prototype completed all ten trials with identical enrollment, stopping reasons, and success classifications. Their largest final-probability difference was `1.75e-12`. This is a timing example, not a calibrated design or an estimate of power with useful precision.

Across the 32 probability comparisons, the split prototype's largest absolute difference from its reference was `9.1e-14`. The current calculation had five discrepancies exceeding `1e-8`, all in concentrated stress cases; the largest was approximately 0.750012. The unsplit quantile prototype missed a tail of approximately 0.000563 in the unequal-arm stress case and failed one of ten trial calculations with an integration error. Its incomplete trial batch is therefore recorded without a runtime comparison. Failure counts and the error message are retained in the CSV.

The split prototype improves accuracy in these checks at a substantial computational cost. These results support further work on a cheaper robust calculation before adopting this prototype throughout the package; they do not establish accuracy for every possible prior, margin, or count table. The CSV retains all case settings, probabilities, reference methods, timing quartiles, failure information, software versions, and the package source commit.

## RMST

`Rscript benchmarks/rmst.R` compares the exact completed-outcome calculation
with Kaplan-Meier integration and benchmarks a full predictive trial against
Cox. An optional integer argument sets the benchmark iterations. The empirical
calculation uses the Greenwood variance with denominator `n^2`; equivalence
checks precede timings.

`Rscript benchmarks/rmst-calibration.R 500 /tmp/rmst-calibration.csv` runs
fixed and adaptive designs under equal-survival and equal-RMST crossing-curve
nulls, a nonzero null margin, and delayed benefit. Variants include arm-specific
dropout, final imputation, and a constant-hazard predictive model fitted to
piecewise-generated outcomes. Output includes failures, success rates, exact
binomial intervals, Monte Carlo standard errors, and mean enrollment. Rates
condition on successfully evaluated trials; failures must be assessed alongside
them. These runs characterize the example thresholds; they do not establish
calibration for arbitrary designs or remove imputation-model assumptions.

The recorded implementation-validation run is in `rmst-validation.md`, with
raw grid and independent fixed-null results in `rmst-calibration-results.csv`
and `rmst-fixed-null-results.csv`.

## Reusable trial paths

`Rscript benchmarks/trial-paths.R 100 20 /tmp/trial-path-results.csv` compares one complete bank plus eight threshold replays with eight separate `sim_trials()` calls. It covers null and alternative scenarios for log-rank and Monte Carlo Bayesian binary analyses. The arguments specify the trial and predictive-imputation counts; Bayesian analyses use 100 posterior draws. The script records elapsed times, bank size, operating characteristics, failures, and R/package versions. Larger counts can be supplied for a more precise comparison.

The saved `trial-path-results.csv` run used 100 trials per scenario, 20 imputations, two interim looks, and sequential execution on macOS with R 4.5.3. The times below include bank generation and all eight replays. They describe these workloads, not a guaranteed speedup for other designs or hardware.

| Analysis | Scenario | Bank and eight replays (s) | Eight direct runs (s) | Speedup |
|:--|:--|--:|--:|--:|
| Log-rank | Null | 1.742 | 8.125 | 4.66 |
| Log-rank | Alternative | 1.469 | 8.525 | 5.80 |
| Bayesian binary | Null | 1.812 | 10.363 | 5.72 |
| Bayesian binary | Alternative | 1.734 | 11.314 | 6.52 |

No bank calculations or direct trials failed in this run. Across the 32 candidate/scenario/method comparisons, the largest absolute difference in success probability was 1.43 combined Monte Carlo standard errors. The direct simulations use a different seed; these small runs are a consistency check, not evidence of calibrated type I error control. Exact comparisons against a direct early-stopping reference using matched stage streams are covered by the package tests across all seven analysis methods, including final-imputation configurations.

The bank stores scores and cohort summaries, not participant datasets or posterior draws. Its memory requirement grows with the number of trials, looks, and predictive imputations. Use different explicit seeds and separate RDS files when a whole bank would exceed available memory. Trial identifiers start at 1 in each call, so retain a batch identifier when combining results.

## Bayesian survival final analysis

Author: Graeme L. Hickey.

`survival-final.R` compares the direct observed-data posterior introduced in commit `c251b22` with the final-imputation calculation from commit `7cc2294`. The comparison uses `method = "bayes-surv"`, `empty_interval = "prior"`, and `imputed_final = TRUE`. The old defaults use 500 predictive completions and 1,000 posterior draws per completion. Direct calculations use 1,000, 10,000, 50,000, and 500,000 independent posterior draws, plus a draw count chosen to match the estimated probability precision of the old defaults. Interim calculations and package defaults are unchanged.

Run from the package root with `pkgload` and `bench` installed:

```sh
Rscript benchmarks/survival-final.R 100 500 20 benchmarks/survival-final
Rscript benchmarks/survival-final-report.R benchmarks/survival-final
```

The first three arguments are the numbers of repeated legacy estimates, repeated direct estimates at 1,000 draws, and timing repetitions. The final argument is the output prefix. The historical commit must be available in the local Git history. The report script uses saved results and requires `ggplot2`; it does not repeat the benchmark. `survival-final-results.rds` retains the simulated participant datasets, scenario settings, sufficient statistics, individual estimates and seeds, timing measurements, reference calculations, software versions, and both final-analysis function definitions. The CSV files provide the summary measurements and exact binomial decision-error calculations. The figure `survival-final-precision.png` shows numerical decision error near the two success thresholds.

All six scenarios use a 12-month endpoint and the default independent Gamma(0.1, 0.1) shape-rate hazard priors. Two-arm designs use equal allocation and test whether the treatment-minus-control event probability is below zero. The single-arm analysis tests whether the event probability is below 0.40. Both success thresholds, 0.95 and 0.975, are considered. Dropout is independent exponential censoring within each arm; the specified probability is the dropout-time CDF at month 12, not the realized proportion lost before an event. Event hazards are constant or change at months 3 and 6. These fixed datasets are computational examples, not proposed or calibrated clinical-trial designs.

| Scenario | Total N | Event risk, control / treatment | Dropout probability, control / treatment |
|:--|--:|--:|--:|
| Constant hazard, lower dropout | 400 | 30% / 22% | 10% / 10% |
| Constant hazard, higher dropout | 400 | 30% / 22% | 30% / 30% |
| Three intervals | 400 | 36.2% / 28.1% | 30% / 30% |
| Three intervals, unequal dropout | 400 | 36.2% / 28.1% | 10% / 40% |
| Single arm | 120 | Benchmark 40% / treatment 30% | Treatment 20% |
| Small sparse trial, three intervals | 80 | 13.4% / 9.2% | 50% / 50% |

The first two datasets are deliberately selected to have posterior probabilities in (0.943, 0.947) and (0.978, 0.982), respectively. The other four use fixed seeds without selection. Every calculation within a scenario uses the same participant data. This separates Monte Carlo variation from repeated-trial sampling variation; it does not estimate power or type I error.

Constant-hazard probability references use a Gamma CDF for the single-arm case and the scaled-Gamma/Beta identity for the two-arm cases. Three-interval references use 20 million independent posterior draws, separate from all benchmark estimates. Posterior mean effects and their direct-sampling variances are available analytically from the Gamma Laplace transforms. Empirical probability and effect standard deviations use 100 legacy repetitions and 500 direct repetitions. For larger direct draw counts, the probability Monte Carlo standard error is calculated as `sqrt(p * (1 - p) / N_mcmc)`, using the reference probability. Approximate 95% intervals for empirical standard deviations are retained in the summary CSV; with 100 repetitions, the relative standard error of an approximately normal estimator's sample standard deviation is about 7%.

Timing measures complete final-analysis calls, including sufficient-statistic preparation and sampling, after loading and warming the functions. It excludes data generation and reference calculations, uses sequential execution, includes garbage collection, and reports the median and quartiles of 20 repetitions. The historical implementation is checked against the retained imputation calculation with `empty_interval = "error"` when all intervals have exposure; results and random-number states must agree exactly.

### Recorded results

The saved run used R 4.5.3 on an Apple ARM64 Mac on 11 November 2026. Times below are median milliseconds per final analysis. They are specific to this machine and workload; the saved quartiles show variation between calls.

| Scenario | Old defaults: 500 imputations x 1,000 draws | Direct: 1,000 draws | Direct: 50,000 draws | Speedup at 50,000 |
|:--|--:|--:|--:|--:|
| Constant hazard, 10% dropout | 327.6 | 0.50 | 8.83 | 37.1 |
| Constant hazard, 30% dropout | 359.2 | 0.54 | 12.41 | 28.9 |
| Three intervals, 30% dropout | 741.3 | 2.34 | 22.51 | 32.9 |
| Three intervals, unequal dropout | 507.6 | 0.83 | 19.77 | 25.7 |
| Single arm, 20% dropout | 190.5 | 0.38 | 3.64 | 52.4 |
| Small sparse trial, 50% dropout | 599.7 | 0.97 | 26.30 | 22.8 |

Monte Carlo standard errors below are on the probability scale, not percentage points. Old-method values are empirical standard deviations across 100 repetitions; direct-method values use the binomial variance. The last column estimates the direct draw count needed to match the old method's probability precision, using its empirical variance.

| Scenario | Reference probability | Old defaults MCSE | Direct 1,000 MCSE | Direct 50,000 MCSE | Equivalent direct draws |
|:--|--:|--:|--:|--:|--:|
| Constant hazard, 10% dropout | 0.946368 | 0.001136 | 0.007124 | 0.001008 | 39,339 |
| Constant hazard, 30% dropout | 0.979792 | 0.001003 | 0.004450 | 0.000629 | 19,673 |
| Three intervals, 30% dropout | 0.977111 | 0.001309 | 0.004729 | 0.000669 | 13,062 |
| Three intervals, unequal dropout | 0.849850 | 0.004920 | 0.011296 | 0.001598 | 5,272 |
| Single arm, 20% dropout | 0.852932 | 0.003936 | 0.011200 | 0.001584 | 8,096 |
| Small sparse trial, 50% dropout | 0.132172 | 0.004866 | 0.010710 | 0.001515 | 4,846 |

The direct default was 317-671 times faster, with 2.2-6.3 times the probability MCSE of the old defaults. At 50,000 draws it remained 23-52 times faster, with comparable or better precision in these cases. Matching the estimated old-method precision required about 4,800-39,300 direct draws and was 40-228 times faster. Using 500,000 direct draws, matching the old total number of completed-data posterior draws, was still 1.6-5.6 times faster and more precise in every case.

The legacy calculation's 500,000 draws are not 500,000 independent draws from the observed-data posterior: draws within an imputation share that completed dataset. If $q(M)$ is its exact completed-data posterior probability, $m$ is the number of imputations, and $b$ is the number of draws per imputation, its Monte Carlo variance is

$$
\operatorname{Var}(\widehat p_{\mathrm{legacy}})
=\frac{\operatorname{Var}\{q(M)\mid D\}}{m}
+\frac{E[q(M)\{1-q(M)\}\mid D]}{mb}.
$$

This explains why substantially fewer than 500,000 independent direct draws can match its precision. The equivalent counts above are estimates, not universal replacements for `N_impute` and `N_mcmc`.

All probability and mean-effect averages were within two estimated standard errors of their references, allowing for reference Monte Carlo error. The empirical direct probability standard deviations were 0.965-1.023 times the binomial predictions. The three simulated reference probabilities had MCSEs of 0.000033-0.000080; the other probability references and all effect means were analytical. No analysis failed. These checks support agreement of the inferential targets; they do not evaluate trial operating characteristics.

### Adequacy of the default draw count

`N_mcmc = 1000` is a coarse computational setting and is not reliably sufficient for a final decision close to 0.95 or 0.975. Its approximate 95% Monte Carlo half-width is 0.0135 at probability 0.95 and 0.0097 at probability 0.975. This limitation already applied to `imputed_final = FALSE`; the simplification exposes it when replacing the more precise former imputed-final calculation.

The three near-boundary examples illustrate the practical consequence. The errors below are probabilities that the finite-draw decision disagrees with the decision based on the reference posterior probability. Direct values use the exact binomial distribution, including the package's strict `>` comparison. Legacy values are empirical proportions from 100 repetitions; zero observed errors do not imply zero error probability.

| Reference probability | Success threshold | Old defaults, observed errors | Direct 1,000 error probability | Direct 50,000 error probability |
|--:|--:|--:|--:|--:|
| 0.946368 | 0.95 | 0 / 100 | 28.57% | 0.0131% |
| 0.979792 | 0.975 | 0 / 100 | 16.64% | less than 0.0001% |
| 0.977111 | 0.975 | 7 / 100 | 35.56% | 0.0965% |

A final-stage budget of 50,000 draws is a reasonable starting point for these workloads. Its MCSE is approximately 0.000975 near 0.95 and 0.000698 near 0.975. If the underlying posterior probability is only 0.001 from the threshold, however, even 50,000 draws can give the opposite numerical decision about 15.3% or 7.7% of the time, taking the worse side of each threshold. At 500,000 draws those values fall to about 0.063% and 0.00038%, respectively. No fixed finite budget resolves probabilities arbitrarily close to a decision boundary. Prespecify numerical precision and a larger budget for close decisions, and report the Monte Carlo uncertainty.

The present `N_mcmc` argument also controls completed-data posterior sampling inside every interim prediction. Increasing it from 1,000 to 50,000 therefore multiplies that part of the interim draw workload by 50. This benchmark supports a separate final-stage draw budget, rather than increasing the shared default without evaluating interim cost and precision. It does not establish that 1,000 draws are sufficient for interim prediction or that operating characteristics are insensitive to numerical error. No defaults have been changed.
