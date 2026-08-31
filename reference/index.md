# Package index

## Package

- [`goldilocks`](https://graemeleehickey.github.io/goldilocks/reference/goldilocks.md)
  : goldilocks

## Adaptive trial simulation

Core functions for simulating and summarizing Goldilocks adaptive trial
designs.

- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  : Simulate and analyze one Goldilocks adaptive trial
- [`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md)
  : Evaluate an externally observed interim data cut
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  : Simulate one or more clinical trials subject to known design
  parameters and treatment effect
- [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  : Summarize simulations to get operating characteristics
- [`summarise_calendar_time()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_calendar_time.md)
  : Summarize operating characteristics on the calendar-time scale
- [`plot_enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/plot_enrollment.md)
  : Plot an enrollment projection
- [`plot_sim_ocs()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_ocs.md)
  : Plot operating characteristics across simulation scenarios
- [`plot_sim_decisions()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_decisions.md)
  : Plot predictive-probability decision maps
- [`summarise_trial_trace()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_trial_trace.md)
  : Summarize an interim decision path
- [`plot_trial_trace()`](https://graemeleehickey.github.io/goldilocks/reference/plot_trial_trace.md)
  : Plot predictive probabilities and enrollment at interim looks
- [`plot_sim_stopping()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_stopping.md)
  : Plot stopping outcomes from trial simulations
- [`print(`*`<goldilocks_interim>`*`)`](https://graemeleehickey.github.io/goldilocks/reference/print.goldilocks_interim.md)
  : Print an externally evaluated interim analysis
- [`print(`*`<goldilocks_trial>`*`)`](https://graemeleehickey.github.io/goldilocks/reference/print.goldilocks_trial.md)
  : Print an adaptive trial trace result
- [`print(`*`<goldilocks_calendar_summary>`*`)`](https://graemeleehickey.github.io/goldilocks/reference/print.goldilocks_calendar_summary.md)
  : Print a calendar-time operating-characteristic summary

## Trial data generation

Simulate complete trial datasets, including enrollment, randomization,
and time-to-event outcomes.

- [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
  : Simulate a complete clinical trial with event data drawn from a
  piecewise exponential distribution
- [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  : Simulate exact continuous-time enrollment
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
