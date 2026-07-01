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

## Interpreting Results

Runtime and memory allocation depend on hardware, R version, operating system,
and BLAS configuration. Compare results relatively on the same machine rather
than treating any absolute timing as a release threshold.
