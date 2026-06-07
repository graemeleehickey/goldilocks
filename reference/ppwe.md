# Cumulative distribution function of the PWE for a vectorized hazard rate parameter

Extends the [`pwe`](https://rdrr.io/pkg/PWEALL/man/pwe.html) function to
allow for vectorization over the hazard rates.

## Usage

``` r
ppwe(hazard, end_of_study, cutpoints)
```

## Arguments

- hazard:

  matrix. A matrix of hazard rate parameters with number of columns
  equal to the length of the `cutpoints` vector. The number of rows can
  be anything, and is typically dictated by the number of MCMC draws.

- end_of_study:

  scalar. Length of the study; i.e. time at which endpoint will be
  evaluated.

- cutpoints:

  vector. The change-point vector indicating time when the hazard rates
  change. Note the first element of `cutpoints` should always be 0.

## Value

A vector of (0, 1) probabilities from evaluation of the PWE cumulative
distribution function. Length of the vector matches the number of rows
of the `hazard` matrix parameter.
