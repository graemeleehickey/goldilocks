# Summarize an interim decision path

Creates a compact one-row summary of the final interim look and stopping
decision. Pass a `goldilocks_trial` object to include final analysis
information, or pass its `trace` element to summarize the path only.

## Usage

``` r
summarise_trial_trace(x)
```

## Arguments

- x:

  A required `goldilocks_trial` result from
  [`sim_trial()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial.md),
  a `goldilocks_interim` result from
  [`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md),
  or an interim trace data frame.

## Value

A one-row data frame.
