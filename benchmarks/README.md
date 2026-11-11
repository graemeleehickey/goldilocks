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
The separate `survival_adapt_bayes_bin` row reports the complete trial runtime,
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
