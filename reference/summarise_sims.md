# Summarize simulations to get operating characteristics

Summarize simulations to get operating characteristics

## Usage

``` r
summarise_sims(data, max_mcse = NULL)
```

## Arguments

- data:

  A complete result returned by
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md),
  a simulation `data.frame`, or a list of either form. Named list
  elements identify scenarios. Existing `scenario` columns and grouping
  variables are preserved.

- max_mcse:

  Optional named numeric vector giving the largest acceptable Monte
  Carlo standard error for selected estimands. Supported names are
  `power`, `stop_success`, `stop_futility`, `stop_max_N`, `mean_N`,
  `stop_and_fail`, and `failure_rate`. A warning identifies every
  scenario and estimand whose Monte Carlo standard error exceeds its
  target.

## Value

Data frame reporting the operating characteristics, including the power
(which will be equal to the type I error in the null case); the
proportion of trials that stopped for early expected success, futility,
or went to the maximum sample size. The average stopping sample size
(and standard deviation) are also recorded. The proportion of trials
that stopped early for expected success, yet went to ultimately fail are
also reported. Each probability and mean is accompanied by its Monte
Carlo standard error and 95% Monte Carlo confidence limits, with columns
ending in `_mcse`, `_mc_lower`, and `_mc_upper`. Probability intervals
use the Wilson method; the mean sample-size interval uses a t
distribution.

These intervals describe uncertainty from using a finite number of
simulated trials under fixed design and data-generating assumptions.
They are **not clinical confidence intervals**, treatment-effect
intervals, or measures of model uncertainty.

The output always reports `n_used`, the number of successfully analyzed
simulations used by the operating-characteristic estimands. When
complete
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
results are supplied, it also reports requested, analyzed, and failed
counts, the failure rate, effective backend, and seed. For a raw
simulation data frame, requested and failed counts are unknown and are
reported as `NA`. Full call, RNG, parallel, timing, failure, and design
metadata, including the evaluated simulation arguments, are retained in
the `simulation_metadata` attribute.
