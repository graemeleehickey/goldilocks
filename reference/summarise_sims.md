# Summarize simulations to get operating characteristics

Summarize simulations to get operating characteristics

## Usage

``` r
summarise_sims(data)
```

## Arguments

- data:

  list (of data frames) or a single data frame. If summarizing a single
  run of simulations, `data` will be a `data.frame` object returned from
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md).
  If summarizing multiple simulation scenarios, `data` will be a `list`
  object, with each element being a `data.frame` object.

## Value

Data frame reporting the operating characteristics, including the power
(which will be equal to the type I error in the null case); the
proportion of trials that stopped for early expected success, futility,
or went to the maximum sample size. The average stopping sample size
(and standard deviation) are also recorded. The proportion of trials
that stopped early for expected success, yet went to ultimately fail are
also reported.
