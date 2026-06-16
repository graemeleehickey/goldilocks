# Package architecture

## Overview

The `goldilocks` package implements the Goldilocks adaptive trial design
described in Broglio et al. (2014). This vignette provides a visual
overview of how the package functions are interconnected.

## Function dependency diagram

The diagram below shows the call graph from the top-level simulation
function
([`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md))
down through the core engine
([`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md))
and into the internal analysis pipeline.

**Exported functions** are shown in blue. **Internal functions** are
shown in grey.

## Function roles

The functions fall into three layers:

### Simulation layer

- **[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)**:
  Top-level entry point. Runs
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  across multiple trials (optionally in parallel) and collates results.
- **[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)**:
  Summarizes the output of
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md),
  computing operating characteristics such as power, expected sample
  size, and stopping probabilities.

### Trial engine

- **[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)**:
  Simulates a single adaptive trial. Generates data via
  [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md),
  conducts interim analyses using `posterior()` and
  `test_stop_success()`, and performs the final analysis via
  `test_final()`.

### Data generation and analysis utilities

- **[`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)**:
  Generates a complete trial dataset by calling
  [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md),
  [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md),
  and
  [`pwe_sim()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_sim.md).
- **`posterior()`**: Estimates the posterior distribution of piecewise
  exponential hazard rates using a conjugate Gamma model.
- **`analyse_data()`**: Applies the chosen analysis method (`logrank`,
  `cox`, `bayes`, or `chisq`) to an (imputed) dataset.
- **`impute_data()`**: Imputes missing event times for censored subjects
  using
  [`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md)
  or
  [`pwe_sim()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_sim.md).
- **`haz_to_prop()`**: Converts posterior hazard rate draws to
  cumulative incidence proportions via
  [`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md).
