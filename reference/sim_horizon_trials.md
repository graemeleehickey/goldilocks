# Simulate trials with a prespecified endpoint-horizon schedule

Simulates only reached interim looks under a prespecified schedule,
using separate success and futility prediction data cuts. Use
[`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md)
when comparing many schedules or stopping thresholds.

## Usage

``` r
sim_horizon_trials(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = NULL,
  N_total,
  lambda = 0.3,
  lambda_time = NULL,
  endpoint_schedule,
  final_horizon,
  dropout_horizon = final_horizon,
  horizon_settings = NULL,
  prior_surv = c(0.1, 0.1),
  prior_bin = c(1, 1),
  bin_method = "mc",
  block = 2,
  rand_ratio = c(control = 1, treatment = 1),
  prop_loss = 0,
  alternative = "greater",
  h0 = 0,
  N_impute = 500,
  N_mcmc = 1000,
  mc_conf_level = 0.95,
  N_trials = 10,
  method = "logrank",
  imputed_final = FALSE,
  empty_interval = c("prior", "propagate", "error"),
  ncores = 1L,
  backend = c("auto", "fork", "psock", "sequential"),
  seed = NULL,
  prior_surv_final = prior_surv,
  generation_cutpoints = cutpoints,
  Fn = 0.05,
  Sn = 0.9,
  Qn = 1,
  prob_ha = 0.95,
  return_trace = FALSE
)
```

## Arguments

- hazard_treatment:

  A required numeric vector of finite, non-negative event rates for the
  treatment arm. Supply one rate per interval defined by
  `generation_cutpoints`; a single value specifies a constant event
  rate.

- hazard_control:

  `NULL` (the default) for a single-arm trial, or a numeric vector of
  finite, non-negative event rates for the control arm in a two-arm
  trial. It must contain one rate per interval defined by
  `generation_cutpoints`.

- cutpoints:

  Master piecewise-exponential analysis cutpoints, strictly inside
  `final_horizon`. Each endpoint uses only cutpoints strictly before
  that endpoint and the corresponding leading intervals of `prior_surv`
  and `prior_surv_final`. If an endpoint H lies inside a master
  interval, the same hazard parameter and unscaled Gamma prior apply,
  but its success update uses events and exposure only through H.
  Available observations beyond H are excluded even within that
  interval; H is not a new hazard-change assumption in the generating
  model.

- N_total:

  A required positive integer giving the maximum total sample size.

- lambda:

  A numeric vector of finite, positive enrollment rates per unit of
  calendar time. Supply one rate for each interval defined by
  `lambda_time`. The default is `0.3`. See
  [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  for the continuous-time enrollment model and time origin.

- lambda_time:

  `NULL` (the default), or a numeric vector of finite, positive,
  strictly increasing calendar times at which the enrollment rate
  changes. Time zero is implicit, and `length(lambda)` must equal
  `length(lambda_time) + 1`.

- endpoint_schedule:

  Data frame with exactly two numeric columns, `N` and `horizon`,
  containing every planned interim enrollment size in increasing order
  and a positive, nondecreasing endpoint horizon at each look. All
  interim sizes must be less than `N_total`. The saved horizons are
  derived as `sort(unique(c(endpoint_schedule$horizon, final_horizon)))`
  and evaluated at every look, allowing subsequent comparison of
  schedules using the same looks and endpoint times. Use a zero-row data
  frame with numeric columns when there are no interims.

- final_horizon:

  Required endpoint time for the maximum-\\N\\ final analysis, all
  futility predictions, and planned outcome collection for every
  participant from enrollment. It must be at least every scheduled
  horizon. The generating hazards must describe this entire period.

- dropout_horizon:

  Reference time at which `prop_loss` specifies the independent
  exponential dropout probability. Defaults to `final_horizon`; it
  remains fixed when the analysis endpoint changes.

- horizon_settings:

  Optional data frame with one row per distinct scheduled horizon plus
  `final_horizon`, identified by column `horizon`. Optional columns are
  `h0`, `rmst_tau`, and a list column `prior_bin`. Missing columns use
  the common `h0`, the endpoint horizon itself for RMST restriction, and
  the common `prior_bin`, respectively. These settings are fixed within
  a set of trial paths.

- prior_surv:

  A numeric vector, matrix, or named list specifying the Gamma prior for
  the piecewise-exponential hazards used to generate outcomes during
  interim prediction. A length-two vector supplies shape and rate and
  applies the same prior to every arm and interval. A `2` by
  `length(cutpoints) + 1` matrix supplies interval-specific values
  shared by all arms, with shapes in row 1 and rates in row 2. For
  independent arm-specific priors, supply a list named `control` and
  `treatment` in a two-arm design, or `treatment` in a single-arm
  design. Each list element may be a length-two vector or an
  interval-specific matrix. Both arms must be supplied; no values are
  borrowed or filled from the other arm. Rates must use the same time
  unit as event times, exposure, and cutpoints. The default is
  `c(0.1, 0.1)`.

- prior_bin:

  A length-two numeric vector of finite, positive shape parameters
  `c(a, b)` for the `Beta(a, b)` event-probability prior used when
  `method = "bayes-bin"`. The same prior is applied to both arms. The
  default is `c(1, 1)`, a uniform prior.

- bin_method:

  A single character string selecting how to calculate the posterior
  probability for `method = "bayes-bin"`. It must be one of `"mc"`
  (Monte Carlo sampling), `"normal"` (normal approximation), or
  `"quadrature"` (numerical integration). The default is `"mc"`. The
  normal approximation can be inaccurate with sparse events or
  non-events and posterior event probabilities near 0 or 1. It can
  change whether `prob_ha` is exceeded. Increasing `N_mcmc` does not
  improve this approximation; use `"quadrature"` or sufficiently precise
  `"mc"` instead. In two-arm designs, quadrature can also return an
  inaccurate posterior probability without an error when the posterior
  event risks are tightly concentrated, particularly with many
  participants but very few events or non-events. Check such
  probabilities against sufficiently precise Monte Carlo Beta draws
  (`"mc"`). Increasing `N_mcmc` does not improve quadrature. This
  integration issue does not apply to the single-arm Beta CDF
  calculation.

- block:

  A positive integer vector of permitted randomization block sizes.
  Every value must be a multiple of `sum(rand_ratio)`. The default is
  `2` and the argument is ignored for a single-arm trial.

- rand_ratio:

  A length-two positive integer vector giving the control to treatment
  randomization ratio. The default is `c(control = 1, treatment = 1)`.
  Name the values `control` and `treatment`; either supplied order is
  accepted and matched by name. A legacy unnamed vector remains accepted
  in `c(control, treatment)` order. Unequal unnamed values produce a
  warning because names may be required in a future major release. See
  [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  for more details.

- prop_loss:

  Dropout probability by `dropout_horizon`, with the same arm-specific
  conventions as
  [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md).
  The corresponding dropout hazard is held fixed across all endpoint
  horizons.

- alternative:

  A single character string specifying the alternative hypothesis. It
  must be one of `"greater"` (the default), `"less"`, or `"two.sided"`.
  One-sided alternatives (`"greater"` and `"less"`) are supported for
  `method = "bayes-surv"` and `method = "bayes-bin"`. All three options
  are supported for `method = "logrank"`, `method = "cox"`,
  `method = "rmst"`, `method = "riskdiff-wald"`, and
  `method = "riskdiff-fm"`. For an adverse event, benefit is in the
  `"greater"` direction for RMST (longer event-free time) and the
  `"less"` direction for the other methods (lower hazard or event
  probability).

- h0:

  Common null value or margin, overridden by `horizon_settings$h0` when
  supplied. Bayesian and risk-difference methods use event-probability
  scales, Cox uses log hazard ratio, RMST uses treatment-control
  differences in time units, and log-rank requires zero. See
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  for method-specific directions and allowable ranges.

- N_impute:

  A positive integer giving the number of predictive imputations used at
  each interim look and, when requested, for final multiple imputation.
  The default is `500`. An imputed Cox, RMST, or risk-difference final
  analysis requires at least two.

- N_mcmc:

  A positive integer giving the number of posterior draws used for
  Bayesian survival analyses and Bayesian binomial analyses with
  `bin_method = "mc"`. The default is `1000`. Ignored for all other
  analyses, including `bin_method = "normal"` and `"quadrature"` with
  `method = "bayes-bin"`; `NULL` is allowed when ignored.

- mc_conf_level:

  A single numeric probability strictly between `0.5` and `1`, giving
  the confidence level for one-sided exact binomial bounds reported as
  diagnostics of finite Monte Carlo uncertainty. The bounds do not alter
  completed-data success classifications or interim decisions, which use
  strict point-estimate comparisons with `prob_ha`, `Qn`, `Sn`, and
  `Fn`. The default is `0.95`.

- N_trials:

  A positive integer giving the number of independent trials to
  simulate. The default is `10`.

- method:

  A single character string specifying the completed-data and final
  analysis. Available choices are a log-rank (`method = "logrank"`)
  test, Cox proportional hazards regression model Wald test
  (`method = "cox"`), a restricted mean survival time difference Wald
  test (`method = "rmst"`), a fully-Bayesian piecewise-exponential
  analysis (`method = "bayes-surv"`), a Bayesian beta-binomial analysis
  of complete binary outcomes (`method = "bayes-bin"`), a frequentist
  risk-difference Wald test (`method = "riskdiff-wald"`), or a
  Farrington-Manning score test (`method = "riskdiff-fm"`) of complete
  binary outcomes. The deprecated `method = "riskdiff"` is accepted as
  an alias for `"riskdiff-wald"` with a warning. The default is
  `"logrank"`. See Details.

- imputed_final:

  Logical; whether to impute outcomes lost before the selected endpoint
  in final analyses. Defaults to `FALSE`. Final-imputation methods and
  restrictions are the same as in
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md):
  log-rank imputation is unavailable, and Farrington-Manning imputation
  requires no missing outcomes. Cox, RMST, and risk-difference Wald
  imputation require at least two imputations and positive pooled
  variance.

- empty_interval:

  A single character string specifying how to handle empty
  piecewise-exponential intervals when updating Gamma hazard models for
  predictive imputation and Bayesian survival analysis. An empty
  interval is an interval with no exposed subjects in a treatment arm at
  the analysis time. `"prior"` (the default) leaves the interval at zero
  exposure time and zero events, so its posterior is driven only by its
  assigned survival prior. `"propagate"` is a legacy heuristic that
  copies exposure time and event counts from the nearest non-empty
  interval in the same treatment arm and emits a warning. `"error"`
  stops when any empty interval is found.

- ncores:

  A positive integer giving the maximum number of processor cores to
  use. The default is `1L`, which runs trials sequentially. The number
  actually used cannot exceed `N_trials`; with `backend = "auto"`, at
  least two trials are required per core to justify the
  parallel-processing overhead.

- backend:

  A single character string selecting the computational method. `"auto"`
  (the default) runs sequentially when `ncores = 1` or fewer than four
  trials are requested; otherwise it uses fork-based parallelization on
  Unix-like systems and a PSOCK cluster on Windows. `"fork"`, `"psock"`,
  and `"sequential"` select a method explicitly. Forking is unavailable
  on Windows.

- seed:

  `NULL` (the default), or a single integer between `0` and
  `.Machine$integer.max`. A supplied seed gives reproducible
  simulations, including when trials are run in parallel, and leaves the
  pre-existing random-number state unchanged.

- prior_surv_final:

  Gamma prior for final-stage imputation and Bayesian survival analysis
  of completed predictive and final data. Accepts the same master-grid
  forms as `prior_surv` and defaults to it. At each endpoint, only the
  corresponding leading intervals are used. Supply both priors
  explicitly to separate predictive and analysis assumptions.

- generation_cutpoints:

  `NULL`, or a numeric vector of finite, positive, strictly increasing
  interior follow-up times defining the piecewise-exponential model used
  to generate event times. `hazard_treatment` and `hazard_control` must
  each have one value per resulting interval. Defaults to `cutpoints`,
  preserving the historical behavior in which generation and analysis
  used one partition.

- Fn:

  Futility threshold: one probability repeated at every interim or one
  per scheduled look. `NULL` or zero disables futility; the default is
  `0.05`.

- Sn:

  Expected-success threshold: one probability repeated at every interim
  or one per scheduled look. The default is `0.9`; one disables
  expected-success stopping.

- Qn:

  Immediate-success threshold: one probability repeated at every interim
  or one per scheduled look. It must be at least `Sn` at every look. The
  default, one, disables immediate-success stopping.

- prob_ha:

  A single numeric probability in `[0, 1]` defining success in each
  completed-data analysis. For Bayesian methods this is compared with
  the posterior probability of the alternative; for frequentist methods
  it is compared with `1 - P`. The default is `0.95`.

- return_trace:

  Logical; whether to include a `traces` data frame containing the
  reached interim decisions in the returned list. Defaults to `FALSE`.

## Value

A simulation result as described in
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md).
The `calculation_failures` component additionally records failed
simulation stages. The complete trial paths are not retained.

## Details

At a 12-month look with an eventual 18-month final endpoint, success
imputation is fitted only to available data through month 12; futility
imputation is fitted to available data through month 18. Both respect
the look's calendar date. The schedule changes only if the trial reaches
the corresponding look. It is not selected using interim outcomes.

Later interim calculations are skipped after a stop. The selected
cohort's final analysis and the maximum-\\N\\ final analysis are
retained, including the counterfactual futility comparison. For the same
seed, set of horizons, and other simulation settings, results agree with
[`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md)
followed by
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md).
Event-time imputation, dropout interpretation, horizon-specific analysis
settings, and the endpoint lock after a stop are described in
[`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md).

## Examples

``` r
schedule <- data.frame(N = c(20, 40), horizon = c(12, 18))
result <- sim_horizon_trials(
  hazard_treatment = 0.04, hazard_control = 0.07,
  N_total = 60, endpoint_schedule = schedule, final_horizon = 18,
  lambda = 3, method = "bayes-bin", bin_method = "quadrature",
  alternative = "less", N_impute = 5, N_trials = 2, seed = 123
)
summarise_sims(result)
#> # A tibble: 1 × 50
#>   scenario backend     seed n_requested n_analyzed n_failed n_used failure_rate
#>   <chr>    <chr>      <dbl>       <int>      <int>    <int>  <int>        <dbl>
#> 1 1        sequential   123           2          2        0      2            0
#> # ℹ 42 more variables: failure_rate_mcse <dbl>, failure_rate_mc_lower <dbl>,
#> #   failure_rate_mc_upper <dbl>, power <dbl>, power_mcse <dbl>,
#> #   power_mc_lower <dbl>, power_mc_upper <dbl>, stop_immediate_success <dbl>,
#> #   stop_immediate_success_mcse <dbl>, stop_immediate_success_mc_lower <dbl>,
#> #   stop_immediate_success_mc_upper <dbl>, stop_success <dbl>,
#> #   stop_success_mcse <dbl>, stop_success_mc_lower <dbl>,
#> #   stop_success_mc_upper <dbl>, stop_any_success <dbl>, …
```
