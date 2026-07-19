# Cumulative distribution function of the PWE for a vectorized hazard rate parameter

Extends [`PWEALL::pwe()`](https://rdrr.io/pkg/PWEALL/man/pwe.html) to
allow for vectorization over the hazard rates.

## Usage

``` r
ppwe(hazard, end_of_study, cutpoints = NULL)
```

## Arguments

- hazard:

  matrix. A matrix of hazard rate parameters with number of columns one
  greater than the length of the `cutpoints` vector. The number of rows
  can be anything, and is typically dictated by the number of MCMC
  draws.

- end_of_study:

  finite positive time at which the cumulative event probability is
  evaluated. It must be greater than every cutpoint.

- cutpoints:

  finite, positive, strictly increasing vector of interior times at
  which the hazard rate changes. The number of hazard rates must be one
  greater than the number of cutpoints. Use `NULL` for a constant
  hazard.

## Value

A vector of (0, 1) probabilities from evaluation of the PWE cumulative
distribution function. Length of the vector matches the number of rows
of the `hazard` matrix parameter.
