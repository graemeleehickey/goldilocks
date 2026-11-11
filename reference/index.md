# Package index

## Package

- [`goldilocks`](https://graemeleehickey.github.io/goldilocks/reference/goldilocks.md)
  : Goldilocks Bayesian adaptive trial designs

## Adaptive trial simulation

Core functions for simulating and summarizing Goldilocks adaptive trial
designs.

- [`sim_trial()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial.md)
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial.md)
  : Simulate and analyze one Goldilocks adaptive trial
- [`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md)
  : Evaluate an externally observed interim data cut
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  : Estimate operating characteristics by trial simulation
- [`sim_trial_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial_paths.md)
  : Simulate complete trial paths for repeated threshold calibration
- [`apply_stopping_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_stopping_rules.md)
  : Apply stopping rules to saved trial paths
- [`print(`*`<goldilocks_paths>`*`)`](https://graemeleehickey.github.io/goldilocks/reference/print.goldilocks_paths.md)
  : Print a set of complete trial paths
- [`sim_horizon_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_trials.md)
  : Simulate trials with a prespecified endpoint-horizon schedule
- [`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md)
  : Simulate complete trial paths for prespecified endpoint horizons
- [`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md)
  : Apply a prespecified endpoint schedule to saved trial paths
- [`print(`*`<goldilocks_horizon_paths>`*`)`](https://graemeleehickey.github.io/goldilocks/reference/print.goldilocks_horizon_paths.md)
  : Print trial paths with multiple endpoint horizons
- [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  : Estimate operating characteristics from trial simulations
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
  : Print a Goldilocks adaptive trial result
- [`print(`*`<goldilocks_calendar_summary>`*`)`](https://graemeleehickey.github.io/goldilocks/reference/print.goldilocks_calendar_summary.md)
  : Print a calendar-time operating-characteristic summary

## Trial data generation

Simulate complete trial datasets, including enrollment, randomization,
and time-to-event outcomes.

- [`sim_trial_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial_data.md)
  [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial_data.md)
  : Simulate complete trial data under piecewise-exponential event rates
- [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  : Simulate exact continuous-time enrollment
- [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  : Generate a block-randomized treatment sequence

## Piecewise exponential utilities

Functions for simulating, imputing, and computing distributions under
the piecewise exponential model.

- [`pwe_sim()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_sim.md)
  : Simulate piecewise exponential time-to-event outcomes
- [`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md)
  : Impute piecewise exponential time-to-event outcomes
- [`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md)
  : Calculate endpoint event probabilities from piecewise hazards
- [`prob_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prob_to_haz.md)
  [`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prob_to_haz.md)
  : Derive piecewise-constant hazard rates from cumulative event
  probabilities
