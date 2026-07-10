# Maintainer Benchmarks

This directory contains optional performance benchmarks for goldilocks hot
paths. They are intended for maintainers to compare branches before and after
optimization work, not as pass/fail tests.

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

The cross-platform parallel benchmark compares serial, PSOCK, and (on Unix-like
platforms) the existing fork backend over small-overhead, representative
log-rank, and heavier Bayesian-survival workloads. It records median elapsed
time, serial-relative speedup, R version, and platform metadata. It also
records parent-process allocation for the serial baseline; allocation profiling
is unavailable for parallel rows because bench cannot profile expressions that
create child processes. PSOCK timings include cluster startup because each
simulation call creates a cluster.

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

## Interpreting Results

Runtime and memory allocation depend on hardware, R version, operating system,
and BLAS configuration. Compare results relatively on the same machine rather
than treating any absolute timing as a release threshold.
