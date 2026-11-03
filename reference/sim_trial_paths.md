# Simulate complete trial paths for repeated threshold calibration

Generates one reusable simulation bank under fixed design and
data-generating assumptions. Every planned interim look is evaluated,
regardless of hypothetical stopping decisions. Use
[`apply_stopping_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_stopping_rules.md)
to calibrate `Fn`, `Sn`, `Qn`, and `prob_ha` without new simulations.

## Usage

``` r
sim_trial_paths(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = NULL,
  N_total,
  lambda = 0.3,
  lambda_time = NULL,
  interim_look = NULL,
  end_of_study,
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
  binary_imputation = c("event-time", "bernoulli"),
  prior_surv_final = prior_surv,
  generation_cutpoints = cutpoints,
  rmst_tau = end_of_study,
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

  `NULL` (the default), or a numeric vector of finite, positive,
  strictly increasing interior follow-up times defining the
  piecewise-exponential model used for interim posterior estimation,
  predictive imputation, and final analysis. The number of
  interval-specific prior columns must be one greater than the number of
  cutpoints. `NULL` specifies a constant-hazard analysis model.

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

- end_of_study:

  A required finite, positive numeric value giving the planned
  subject-level follow-up time. It must be greater than the final value
  in both `cutpoints` and `generation_cutpoints`, when supplied, and use
  the same time unit.

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

  A numeric vector containing one or two probabilities in `[0, 1)`. Each
  value is the dropout-time CDF at `end_of_study`: \\P(D \le \tau) =
  p\\, where \\\tau\\ is the planned follow-up duration per subject.
  Independently of event time and enrollment, each subject's dropout
  time \\D\\ is exponentially distributed with rate \\-\log(1-p)/\tau\\.
  The observed time is the minimum of event time, dropout time, and
  `end_of_study`; an event occurring before dropout is retained. Thus,
  `prop_loss` is not the expected proportion actually censored by
  dropout: that proportion can be lower because events occur first, and
  the realized number of dropouts varies between trials. A single value
  applies the same dropout distribution to every arm. For a two-arm
  design, supply a length-two vector named `control` and `treatment` for
  arm-specific probabilities; supplied order does not matter. Single-arm
  designs require one probability. The default `0` sets dropout time to
  infinity without drawing random numbers. A value of `1` is rejected
  because it requires an infinite exponential rate.

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

  A single finite numeric value specifying the null hypothesis or
  margin. The default is `0`. For Bayesian analyses, `h0` must lie in
  `[0, 1]` for a single-arm design and `[-1, 1]` for a two-arm design.

  - When `method = "bayes-surv"`, `h0` is the null value of
    \\p\_\textrm{treatment} - p\_\textrm{control}\\. In a single-arm
    design, `h0` is the external benchmark event probability, often
    referred to as a performance goal (PG) or objective performance
    criterion (OPC).

  - When `method = "bayes-bin"`, `h0` is the null value of
    \\p\_\textrm{treatment} - p\_\textrm{control}\\ for a two-arm
    design, or the null event probability for a single-arm design.

  - When `method = "cox"`, `h0` is the null log hazard ratio for
    treatment versus control. Use `h0 = 0` for the usual hazard ratio of
    1 null, or `h0 = log(margin)` for a non-inferiority margin specified
    as a hazard ratio. A Cox non-inferiority test should usually use
    `alternative = "less"`.

  - When `method = "rmst"`, `h0` is the null treatment-control RMST
    difference in time units and must lie in `[-rmst_tau, rmst_tau]`.
    For non-inferiority allowing a loss of `m` time units, use `h0 = -m`
    and `alternative = "greater"`.

  - When `method = "riskdiff-wald"` or `method = "riskdiff-fm"`, `h0` is
    the null value of \\p\_\textrm{treatment} - p\_\textrm{control}\\
    and must lie in `[-1, 1]`.

  - When `method = "logrank"`, only `h0 = 0` is supported; this denotes
    the usual equal-survival null. Nonzero values are rejected because
    the standard log-rank statistic does not implement a nonzero effect
    margin.

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

  A single logical value indicating whether the final analysis should be
  based on imputed outcomes for subjects who were LTFU (i.e.
  right-censored with time less than `end_of_study`). The default is
  `FALSE`, which uses the observed-data analysis. If no outcomes require
  imputation, the selected complete-data test is used directly with
  either flag. With missing outcomes and `method = "cox"`, `"rmst"`, or
  `"riskdiff-wald"`, setting this to `TRUE` pools the scalar treatment
  effects and variances using Rubin's rules; this requires
  `N_impute >= 2` and positive total variance. Genuine final imputation
  is unsupported for `method = "riskdiff-fm"` because no validated FM
  pooling rule is implemented. Simulations combining FM and
  `imputed_final = TRUE` therefore require `prop_loss = 0` in both arms.
  Imputed final analyses remain unavailable for `method = "logrank"`.

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

- binary_imputation:

  A single character string selecting the predictive imputation approach
  for `method = "bayes-bin"`, `method = "riskdiff-wald"`, or
  `method = "riskdiff-fm"`. `"event-time"` (the default) draws a
  conditional piecewise-exponential event time and reduces it to event
  status at `end_of_study`. `"bernoulli"` draws the endpoint status
  directly from its conditional event probability. This argument is
  ignored for time-to-event analysis methods.

- prior_surv_final:

  A numeric vector, matrix, or named list specifying the Gamma prior
  used for final-stage piecewise-exponential imputation and, for
  `method = "bayes-surv"`, both the analysis of each hypothetical
  completed trial at interim looks and the actual final analysis. It
  accepts the same shared or arm-specific forms as `prior_surv` and
  defaults to `prior_surv`. An informative `prior_surv` can therefore
  predict outstanding outcomes while a weak `prior_surv_final` defines
  the Bayesian survival success criterion. To use different priors for
  these roles, supply `prior_surv_final` explicitly; an informative
  predictive prior is otherwise also the default analysis prior. See
  **Predictive and analysis priors** below.

- generation_cutpoints:

  `NULL`, or a numeric vector of finite, positive, strictly increasing
  interior follow-up times defining the piecewise-exponential model used
  to generate event times. `hazard_treatment` and `hazard_control` must
  each have one value per resulting interval. Defaults to `cutpoints`,
  preserving the historical behavior in which generation and analysis
  used one partition.

- rmst_tau:

  A single finite positive restriction time for `method = "rmst"`, in
  the same units as `end_of_study`. Defaults to `end_of_study` and must
  not exceed it. Prespecify the same horizon for all looks, imputations,
  and simulations. It may precede analysis cutpoints and does not
  shorten the planned follow-up or imputation horizon. Ignored for other
  methods.

- trial_offset:

  A non-negative integer. Trial identifiers start at `trial_offset + 1`.
  With an explicit seed, consecutive offsets allow a large bank to be
  generated, saved, and evaluated in separate batches using the same
  streams as a single call. Keep all other simulation settings fixed.

## Value

A `goldilocks_paths` object containing `trials` (one record per
requested trial), `failures`, `schema_version`, and `call`. Each trial
has `interims` with unthresholded `scores`, observed `context`, and
diagnostics, and `finals` with one row per possible stopping cohort.
Attributes retain evaluated `arguments`, `prior_design`,
`enrollment_design`, `rng_metadata`, `parallel_metadata`,
`runtime_metadata`, and `score_storage_bytes`. Save the whole object
with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) to retain its
metadata and posterior counts.

## Details

The bank retains the success score (posterior probability or frequentist
`1 - p`) of every predictive replicate, for the current and maximum
cohorts, plus posterior Monte Carlo counts when applicable. It also
retains a final analysis and calendar metrics for each possible stopping
cohort. Final analyses use eventual observed follow-up, including the
specified dropout and final-imputation policy. Interim evidence only
uses information available at its calendar cut.

The data-generating model, priors, method, null margin, endpoint
horizon, allocation, look schedule, maximum sample size, `N_impute`, and
`N_mcmc` are fixed within a bank. Changing them requires a new bank.
Interim looks must be strictly below `N_total`; the maximum-N analysis
is after follow-up.

A bank costs more than one early-stopping simulation but can be reused
for many candidate rules. Score storage alone requires approximately
`16 * N_trials * length(interim_look) * N_impute` bytes, or 24 bytes per
trial/look/imputation for Monte Carlo Bayesian methods, plus object and
diagnostic overhead. Large estimates are reported before simulation.
Banks are held in memory; use `trial_offset` and
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) for separate batches.

Each trial and stage receives a separate L'Ecuyer-CMRG stream or
substream. An explicit seed preserves the caller's RNG state and gives
identical paths across supported backends. With `seed = NULL`, exactly
one integer seed is drawn from the caller's RNG; all backends use that
seed. Seeded bank results need not equal
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
results, whose draw scheduling is unchanged.

Failures are retained at the calculation level. A failed predictive
replicate invalidates its entire current- or maximum-cohort probability;
replicates are never silently dropped. Later looks and cohort final
analyses are still attempted. Rule application determines which failures
prevent a particular candidate from being evaluated.

## See also

[`apply_stopping_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_stopping_rules.md),
[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)

## Examples

``` r
paths <- sim_trial_paths(
  hazard_treatment = 0.05, hazard_control = 0.08,
  N_total = 60, interim_look = c(20, 40), lambda = 5,
  end_of_study = 12, alternative = "less", method = "bayes-bin",
  bin_method = "quadrature", N_impute = 5, N_trials = 2, seed = 123
)
result <- apply_stopping_rules(paths, Fn = 0.05, Sn = 0.9, prob_ha = 0.975)
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
