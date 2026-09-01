# Maintainer Benchmarks

This directory contains optional performance benchmarks for goldilocks hot
paths. They are intended for maintainers to compare branches before and after
optimization work, not as pass/fail tests.

The hot-path benchmark includes constant and fractional-knot enrollment
schedules. In particular, its low-rate case guards the continuous-time
generator against work that scales with empty calendar bins rather than with
the requested number of enrollments.

It also reports piecewise-exponential posterior sampling from a patient-level
data frame and from precomputed sufficient statistics. The sufficient-statistic
pair compares the general checked entry point with the trusted internal kernel
used for canonical summaries generated inside Bayesian predictive imputation.
Both receive the same normalized prior and reset to the same seed, isolating
validation overhead without changing the Gamma posterior calculation.

The Bayesian-survival effect pair compares the checked `haz_to_prop()` route
with the trusted completed-data kernel. The latter reuses fixed analysis-
interval widths and returns treatment-effect draws directly, avoiding repeated
cutpoint validation and temporary probability data frames. It deliberately
retains the outer imputation loop so seeded results and bounded peak memory are
unchanged.

The stable cumulative-hazard/probability transformations are benchmarked over
100,000 values spanning near-zero inputs through their mathematical
boundaries. This guards the use of `expm1()` and `log1p()` against a material
regression in simulation hot paths.

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

The benchmarks are excluded from R package builds with `.Rbuildignore`, so they
do not run on CRAN or during ordinary package checks.

## Running

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
