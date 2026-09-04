# Estimate operating characteristics by trial simulation

Repeats
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
under fixed design and data-generating assumptions, returning
trial-level results from which operating characteristics can be
estimated.

## Usage

``` r
sim_trials(
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
  Fn = 0.05,
  Sn = 0.9,
  prob_ha = 0.95,
  N_impute = 500,
  N_mcmc = 1000,
  mc_conf_level = 0.95,
  N_trials = 10,
  method = "logrank",
  imputed_final = FALSE,
  empty_interval = c("prior", "propagate", "error"),
  return_trace = FALSE,
  ncores = 1L,
  backend = c("auto", "fork", "psock", "sequential"),
  seed = NULL,
  binary_imputation = c("event-time", "bernoulli"),
  prior_surv_final = prior_surv,
  generation_cutpoints = cutpoints,
  Qn = 1
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
  the piecewise-exponential hazards used during interim prediction. A
  length-two vector supplies shape and rate and applies the same prior
  to every arm and interval. A `2` by `length(cutpoints) + 1` matrix
  supplies interval-specific values shared by all arms, with shapes in
  row 1 and rates in row 2. For independent arm-specific priors, supply
  a list named `control` and `treatment` in a two-arm design, or
  `treatment` in a single-arm design. Each list element may be a
  length-two vector or an interval-specific matrix. Both arms must be
  supplied; no values are borrowed or filled from the other arm. Rates
  must use the same time unit as event times, exposure, and cutpoints.
  The default is `c(0.1, 0.1)`.

- prior_bin:

  A length-two numeric vector of finite, positive shape parameters
  `c(a, b)` for the `Beta(a, b)` event-probability prior used when
  `method = "bayes-bin"`. The same prior is applied to both arms. The
  default is `c(1, 1)`, a uniform prior.

- bin_method:

  A single character string selecting how to calculate the posterior
  probability for `method = "bayes-bin"`. It must be one of `"mc"`
  (Monte Carlo sampling), `"normal"` (normal approximation), or
  `"quadrature"` (numerical integration). The default is `"mc"`.

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

  A numeric vector containing one or two probabilities in `[0, 1]`. A
  single value applies the same loss-to-follow-up proportion to every
  arm. For a two-arm design, differential attrition can be specified
  with a length-two vector named `control` and `treatment`; the supplied
  order does not matter. Within each arm,
  `ceiling(prop_loss * arm size)` subjects are selected at random
  regardless of event status. Each selected subject's observed time is
  drawn from a `Uniform(0, t)` distribution, where `t` is their
  potential event or censoring time. Since the LTFU time is always less
  than `t`, the event has not yet occurred at dropout and the subject is
  right-censored. Single-arm designs require one probability. The
  default is `0`, denoting no loss to follow-up.

- alternative:

  A single character string specifying the alternative hypothesis. It
  must be one of `"greater"` (the default), `"less"`, or `"two.sided"`.
  One-sided alternatives (`"greater"` and `"less"`) are supported for
  `method = "bayes-surv"` and `method = "bayes-bin"`. All three options
  are supported for `method = "logrank"`, `method = "cox"`,
  `method = "riskdiff-wald"`, and `method = "riskdiff-fm"`. For survival
  outcomes, `"less"` corresponds to the treatment arm having a lower
  cumulative incidence (i.e., treatment is beneficial), and `"greater"`
  corresponds to the treatment arm having a higher cumulative incidence.

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

  - When `method = "riskdiff-wald"` or `method = "riskdiff-fm"`, `h0` is
    the null value of \\p\_\textrm{treatment} - p\_\textrm{control}\\
    and must lie in `[-1, 1]`.

  - When `method = "logrank"`, only `h0 = 0` is supported; this denotes
    the usual equal-survival null. Nonzero values are rejected because
    the standard log-rank statistic does not implement a nonzero effect
    margin.

- Fn:

  `NULL`, or a numeric vector of probabilities in `[0, 1]`. Each value
  is the predictive-probability threshold to stop at the \\i\\-th look
  early for futility. If there are no interim looks (i.e.
  `interim_look = NULL`), then `Fn` is not used in the simulations or
  analysis. Set `Fn = 0` to disable futility monitoring; `Fn = NULL` has
  the same effect. Supply either one value, which is repeated at every
  interim look, or exactly one value per `interim_look`. Other lengths
  are rejected rather than recycled. The default is `0.05`.

- Sn:

  A numeric vector of probabilities in `[0, 1]`. Each value is the
  predictive-probability threshold to stop accrual at the \\i\\-th look
  for expected success. If there are no interim looks (i.e.
  `interim_look = NULL`), then `Sn` is not used in the simulations or
  analysis. Supply either one value, which is repeated at every interim
  look, or exactly one value per `interim_look`. Other lengths are
  rejected rather than recycled. The default is `0.9`.

- prob_ha:

  A single numeric probability in `[0, 1]` defining success in each
  completed-data analysis. For Bayesian methods this is compared with
  the posterior probability of the alternative; for frequentist methods
  it is compared with `1 - P`. The default is `0.95`.

- N_impute:

  A positive integer giving the number of predictive imputations used at
  each interim look and, when requested, for final multiple imputation.
  The default is `500`. An imputed Cox or risk-difference final analysis
  requires at least two.

- N_mcmc:

  A positive integer giving the number of posterior draws used within
  each `method = "bayes-surv"` and by `method = "bayes-bin"` when
  `bin_method = "mc"`. The default is `1000`.

- mc_conf_level:

  A single numeric probability strictly between `0.5` and `1`, giving
  the confidence level for one-sided exact binomial bounds reported as
  diagnostics of finite Monte Carlo uncertainty. The bounds do not alter
  completed-data success classifications or interim decisions, which use
  strict point- estimate comparisons with `prob_ha`, `Qn`, `Sn`, and
  `Fn`. The default is `0.95`.

- N_trials:

  A positive integer giving the number of independent trials to
  simulate. The default is `10`.

- method:

  A single character string specifying the completed-data and final
  analysis. Available choices are a log-rank (`method = "logrank"`)
  test, Cox proportional hazards regression model Wald test
  (`method = "cox"`), a fully-Bayesian piecewise-exponential analysis
  (`method = "bayes-surv"`), a Bayesian beta-binomial analysis of
  complete binary outcomes (`method = "bayes-bin"`), a frequentist
  risk-difference Wald test (`method = "riskdiff-wald"`), or a
  Farrington-Manning score test (`method = "riskdiff-fm"`) of complete
  binary outcomes. The deprecated `method = "riskdiff"` is accepted as
  an alias for `"riskdiff-wald"` with a warning. The default is
  `"logrank"`. See Details.

- imputed_final:

  A single logical value indicating whether the final analysis should be
  based on imputed outcomes for subjects who were LTFU (i.e.
  right-censored with time less than `end_of_study`). The default is
  `FALSE`, which means that the final analysis incorporates
  right-censoring. With `method = "cox"`, `method = "riskdiff-wald"`, or
  `method = "riskdiff-fm"`, setting this to `TRUE` analyzes each imputed
  dataset and pools the scalar treatment effects and variances using
  Rubin's rules; this requires `N_impute >= 2`. The pooled
  risk-difference analysis is a Wald test rather than a
  Farrington-Manning test. Imputed final analyses remain unavailable for
  `method = "logrank"`.

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

- return_trace:

  A single logical value indicating whether to retain the compact
  interim decision trace from every simulated trial. The default,
  `FALSE`, preserves the compact output. When `TRUE`, the returned list
  also contains a `traces` data frame with a `trial` column linking each
  trace row to the corresponding original simulated trial.

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
  `method = "bayes-surv"`, final analysis. It accepts the same shared or
  arm-specific forms as `prior_surv` and defaults to `prior_surv`,
  preserving the historical behavior.

- generation_cutpoints:

  `NULL`, or a numeric vector of finite, positive, strictly increasing
  interior follow-up times defining the piecewise-exponential model used
  to generate event times. `hazard_treatment` and `hazard_control` must
  each have one value per resulting interval. Defaults to `cutpoints`,
  preserving the historical behavior in which generation and analysis
  used one partition.

- Qn:

  A numeric vector of probabilities in `[0, 1]`. Each value is the upper
  predictive-probability threshold for declaring immediate trial success
  at the \\i\\-th look. If there are no interim looks (i.e.
  `interim_look = NULL`), then `Qn` is not used in the simulations or
  analysis. Supply either one value, which is repeated at every interim
  look, or exactly one value per `interim_look`; other lengths are
  rejected. `Qn` must be greater than or equal to `Sn` at every look.
  The default, `1`, disables immediate-success stopping.

## Value

A list containing `sims`, a data frame with one row per successfully
simulated trial; `failures`, a data frame with columns `trial`,
`error_class`, and `message`; and `call`. When `return_trace = TRUE`,
the list also contains `traces`, a data frame with one row per completed
interim look and a `trial` identifier. Per-trial calendar-time metrics
are always retained in `sims`; traces additionally retain calendar time
and active follow-up at each look. See
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
for details of the summary and trace columns, and
[`summarise_calendar_time()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_calendar_time.md)
for wide operating-characteristic tables. The returned object also
retains the evaluated `decision_design` and resolved `prior_design`
attributes from
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md).
An `rng_metadata` attribute records the random-number generator,
computational method, and seed policy. A `parallel_metadata` attribute
records the requested and actual computational method and number of
cores. An `arguments` attribute contains a named list of all evaluated
argument values, including defaults. Its `prop_loss` element contains a
named value for every simulated arm, and its `rand_ratio` element is
stored in `control`, `treatment` order for two-arm designs. Its
`cutpoints` and `generation_cutpoints` elements retain the analysis and
data-generation partitions, respectively. For `method = "bayes-bin"`, it
also retains the imputation priors (`prior_surv` and
`prior_surv_final`), completed-data analysis prior (`prior_bin`), and
imputation horizon (`end_of_study`). The attribute can be saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and supplied to a
later call with `do.call(sim_trials, attr(result, "arguments"))`.

## Details

This function is a wrapper for
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
that repeatedly simulates independent trials under the same design
parameters and assumed treatment effect.

To use multiple cores (where available), the argument `ncores` can be
increased from the default of 1. The default `backend = "auto"` stays
sequential for fewer than four trials and otherwise uses no more than
one core per two trials. This avoids parallel-processing overhead for
small simulation studies. On Unix-like systems parallel trials use
forked R processes; on Windows they use PSOCK processes. Set `backend`
explicitly when a particular computational method is required.

Errors raised by an individual
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
call are isolated so other trials can finish. Failed trials are excluded
from `sims`, recorded in `failures` with their trial number, error
class, and message, and reported together in one warning. If every
requested trial fails, `sim_trials()` stops and attaches the same
failure table to the error as `failures`. With a supplied `seed`, the
original call and failed trial number reproduce the same per-trial
random-number stream.

With a supplied `seed`, each trial receives an independent random-number
stream. The resulting trial-level simulations are identical whether they
are run sequentially or with a supported parallel method, and the
pre-existing R random-number state is restored afterward. With
`seed = NULL`, the current random-number state is used and advanced.

## Examples

``` r
hc <- prop_to_haz(c(0.20, 0.30), 12, 36)
ht <- prop_to_haz(c(0.05, 0.15), 12, 36)

out <- sim_trials(
  hazard_treatment = ht,
  hazard_control = hc,
  cutpoints = 12,
  N_total = 600,
  lambda = 20,
  lambda_time = NULL,
  interim_look = c(400, 500),
  end_of_study = 36,
  prior_surv = c(0.1, 0.1),
  block = 2,
  rand_ratio = c(control = 1, treatment = 1),
  prop_loss = 0.30,
  alternative = "two.sided",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.9,
  prob_ha = 0.975,
  N_impute = 5,
  N_mcmc = 5,
  method = "logrank",
  N_trials = 2,
  ncores = 1,
  backend = "auto",
  seed = 123)
```
