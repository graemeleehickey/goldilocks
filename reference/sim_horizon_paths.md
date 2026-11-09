# Simulate complete trial paths for prespecified endpoint horizons

Generates one set of event and dropout histories per trial and evaluates
every planned interim at each supported endpoint horizon. The saved
evidence can be reused by
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md)
to calibrate stopping thresholds and prespecified horizon schedules.

## Usage

``` r
sim_horizon_paths(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = NULL,
  N_total,
  lambda = 0.3,
  lambda_time = NULL,
  interim_look = NULL,
  final_horizon,
  endpoint_horizons = final_horizon,
  followup_horizon = final_horizon,
  dropout_horizon = followup_horizon,
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
  trial_offset = 0L
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
  `followup_horizon`. Each endpoint uses only cutpoints strictly before
  that endpoint and the corresponding leading intervals of `prior_surv`
  and `prior_surv_final`.

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

- interim_look:

  `NULL` (the default) for no interim analyses, or a strictly increasing
  positive integer vector giving the cumulative sample size at each
  interim look. Do not include the maximum sample size. For two-arm
  designs, each interim look must be at least the (largest) block size
  (see `block`), ensuring both treatment groups are present at every
  interim analysis; a smaller look could enroll subjects from one
  treatment group only, leaving the interim posterior undefined for the
  missing group.

- final_horizon:

  Required endpoint time for the maximum-\\N\\ final analysis and all
  futility predictions. It must equal the largest `endpoint_horizons`
  value.

- endpoint_horizons:

  Positive, strictly increasing endpoint times to retain. The default is
  `final_horizon`.

- followup_horizon:

  Planned outcome-collection duration for every participant from
  enrollment, at least `final_horizon`. Defaults to `final_horizon`. The
  generating hazards must describe this entire period.

- dropout_horizon:

  Reference time at which `prop_loss` specifies the independent
  exponential dropout probability. Defaults to `followup_horizon`; it
  remains fixed when the analysis endpoint changes.

- horizon_settings:

  Optional data frame with one row per `endpoint_horizons` value,
  identified by column `horizon`. Optional columns are `h0`, `rmst_tau`,
  and a list column `prior_bin`. Missing columns use the common `h0`,
  the endpoint horizon itself for RMST restriction, and the common
  `prior_bin`, respectively. These settings are fixed within a set of
  trial paths.

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
  `"mc"` instead.

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

  A positive integer giving the number of posterior draws used within
  each `method = "bayes-surv"` and by `method = "bayes-bin"` when
  `bin_method = "mc"`. The default is `1000`.

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

- trial_offset:

  A non-negative integer. Trial identifiers start at `trial_offset + 1`.
  With an explicit seed, consecutive offsets allow a large simulation to
  be generated, saved, and evaluated in separate batches using the same
  streams as a single call. Keep all other simulation settings fixed.

## Value

A `goldilocks_horizon_paths` object with `trials`, a calculation-level
`failures` table including `horizon`, `schema_version`, and `call`. Each
trial retains per-look evidence for each horizon and final analyses for
each cohort and horizon. Separate posterior diagnostics record the event
counts, exposure, and Gamma updates used at each data cut. Attributes
retain evaluated arguments, horizon-specific priors, enrollment design,
and random-number and computation settings. Save and restore paths with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
[`readRDS()`](https://rdrr.io/r/base/readRDS.html). Memory grows with
trials, looks, horizons, and imputation count.

## Details

Success prediction at horizon H fits its predictive hazard model using
only follow-up through H that is available at the look's calendar date.
It imputes completed outcomes through H among currently enrolled
participants. Futility prediction uses a separate model fitted to
available follow-up through `final_horizon`, and predicts success at
that horizon at maximum \\N\\. Follow-up after H cannot inform the
success model, including when a hazard interval crosses H. When H equals
`final_horizon`, the two predictions share a fitted predictive model.

Event times are generated once through `followup_horizon` and
re-censored for each analysis. Predictive imputation uses event times;
independently imputed binary outcomes at different horizons are not
supported. Predictions at different data cuts use their own posterior
draws and are not paired latent event trajectories across horizons.

The final analysis after an expected-success stop uses that look's
endpoint, even if it fails. It never resumes the trial at a longer
endpoint. Futility targets the maximum-\\N\\ final analysis, rather than
success at any remaining interim. Interim sample sizes must be less than
`N_total`; maximum \\N\\ is followed by final endpoint ascertainment.

`prop_loss` describes an exponential dropout distribution independently
of the analysis endpoint. For example,
`prop_loss = 0.05, dropout_horizon = 12, followup_horizon = 18` retains
a 5% dropout probability by month 12 while generating follow-up through
month 18. Dropout after month 12 does not make the 12-month endpoint
missing.

All other fixed design settings and analysis restrictions follow
[`sim_trial_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial_paths.md).
Explicit seeds preserve the caller's random-number state and reproduce
results across supported computation methods and consecutive
`trial_offset` batches. Changing the set of retained horizons changes
the allocation of random numbers. Changing only a candidate schedule or
thresholds does not. Select designs using one set of trial paths and
assess the complete adaptive design using independent simulations.

## See also

[`sim_horizon_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_trials.md),
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md)

## Examples

``` r
paths <- sim_horizon_paths(
  hazard_treatment = c(0.04, 0.03), hazard_control = c(0.07, 0.06),
  generation_cutpoints = 12, cutpoints = 12,
  N_total = 60, interim_look = c(20, 40),
  endpoint_horizons = c(12, 18), final_horizon = 18,
  lambda = 3, method = "bayes-bin", bin_method = "quadrature",
  alternative = "less", N_impute = 5, N_trials = 2, seed = 123
)
schedule <- data.frame(N = c(20, 40), horizon = c(12, 18))
result <- apply_horizon_rules(paths, schedule, return_trace = TRUE)
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
