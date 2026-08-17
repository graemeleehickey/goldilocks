# Summarize simulations to get operating characteristics

Summarize simulations to get operating characteristics

## Usage

``` r
summarise_sims(data)
```

## Arguments

- data:

  A complete result returned by
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md),
  a simulation `data.frame`, or a list of either form. Named list
  elements identify scenarios. Existing `scenario` columns and grouping
  variables are preserved.

## Value

Data frame reporting the operating characteristics, including the power
(which will be equal to the type I error in the null case); the
proportion of trials that stopped for early expected success, futility,
or went to the maximum sample size. The average stopping sample size
(and standard deviation) are also recorded. The proportion of trials
that stopped early for expected success, yet went to ultimately fail are
also reported. When complete
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
results are supplied, the output also reports the requested, analyzed,
and failed trial counts plus the effective backend and seed. Full call,
RNG, parallel, timing, failure, and design metadata are retained in the
`simulation_metadata` attribute.
