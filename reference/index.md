# Package index

## Package

- [`goldilocks`](https://graemeleehickey.github.io/goldilocks/reference/goldilocks.md)
  : goldilocks

## Adaptive trial simulation

Core functions for simulating and summarizing Goldilocks adaptive trial
designs.

- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  : Simulate and execute a single adaptive clinical trial design with a
  time-to-event endpoint
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  : Simulate one or more clinical trials subject to known design
  parameters and treatment effect
- [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  : Summarize simulations to get operating characteristics

## Trial data generation

Simulate complete trial datasets, including enrollment, randomization,
and time-to-event outcomes.

- [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
  : Simulate a complete clinical trial with event data drawn from a
  piecewise exponential distribution
- [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  : Simulate enrollment times
- [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  : Randomization allocation

## Piecewise exponential utilities

Functions for simulating, imputing, and computing distributions under
the piecewise exponential model.

- [`pwe_sim()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_sim.md)
  : Simulate piecewise exponential time-to-event outcomes
- [`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md)
  : Impute piecewise exponential time-to-event outcomes
- [`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md)
  : Cumulative distribution function of the PWE for a vectorized hazard
  rate parameter
- [`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md)
  : Estimate plausible piecewise constant hazard rates from summary
  summary event proportions
