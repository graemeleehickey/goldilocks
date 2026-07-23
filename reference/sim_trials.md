# Simulate one or more clinical trials subject to known design parameters and treatment effect

Simulate multiple clinical trials with fixed input parameters, and
tidily extract the relevant data to generate operating characteristics.

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
  prior = c(0.1, 0.1),
  bin_prior = c(1, 1),
  bin_method = "mc",
  block = 2,
  rand_ratio = c(1, 1),
  prop_loss = 0,
  alternative = "greater",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.9,
  prob_ha = 0.95,
  N_impute = 10,
  N_mcmc = 10,
  N_trials = 10,
  method = "logrank",
  imputed_final = FALSE,
  empty_interval = c("propagate", "prior", "error"),
  return_trace = FALSE,
  ncores = 1L,
  backend = c("auto", "fork", "psock", "sequential"),
  seed = NULL,
  binary_imputation = c("event-time", "bernoulli")
)
```

## Arguments

- hazard_treatment:

  vector. Finite non-negative constant hazard rates under the treatment
  arm.

- hazard_control:

  vector. Finite non-negative constant hazard rates under the control
  arm.

- cutpoints:

  finite, positive, strictly increasing interior times at which the
  baseline hazard changes. The number of hazards for each arm must be
  one greater than the number of cutpoints. Default is `NULL`, which
  corresponds to a simple (non-piecewise) exponential model.

- N_total:

  integer. Maximum sample size allowable

- lambda:

  finite positive enrollment rates per unit time. Supply one rate for
  each interval defined by `lambda_time`. See
  [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  for the precise continuous-time process and time-origin convention.

- lambda_time:

  `NULL`, or finite, positive, strictly increasing internal times at
  which the enrollment rate changes. The initial boundary at zero is
  implicit, so `length(lambda)` must equal `length(lambda_time) + 1`.

- interim_look:

  vector. Sample size for each interim look. Note: the maximum sample
  size should not be included. For two-arm designs, each interim look
  must be at least the (largest) block size (see `block`), ensuring both
  treatment groups are present at every interim analysis; a smaller look
  could enroll subjects from one treatment group only, leaving the
  interim posterior undefined for the missing group.

- end_of_study:

  finite study endpoint, strictly greater than the last cutpoint.

- prior:

  vector. The prior distributions for the piecewise hazard rate
  parameters are each \\Gamma(a_0, b_0)\\, where \\a_0\\ is the shape
  parameter and \\b_0\\ is the rate parameter (i.e., the inverse of the
  scale). This follows R's
  [`stats::rgamma()`](https://rdrr.io/r/stats/GammaDist.html)
  parameterization. The same prior is applied to all piecewise intervals
  and to both treatment groups. The default non-informative prior
  distribution used is `Gamma(0.1, 0.1)`, which is specified by setting
  `prior = c(0.1, 0.1)`.

- bin_prior:

  vector. Prior distribution for the event probability when
  `method = "bayes-bin"`. The two values are the shape parameters of the
  `Beta(a, b)` prior. The same prior is applied to both treatment arms.

- bin_method:

  character. Method used to calculate the posterior probability for
  `method = "bayes-bin"`, must be one of `"mc"` (Monte Carlo sampling),
  `"normal"` (normal approximation), or `"quadrature"` (numerical
  integration). The default is `"mc"`.

- block:

  scalar. Block size for generating the randomization schedule.

- rand_ratio:

  vector. Randomization allocation for the ratio of control to
  treatment. Integer values mapping the size of the block. See
  [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  for more details.

- prop_loss:

  scalar. Overall proportion of subjects lost to follow-up. Subjects are
  selected at random for LTFU regardless of treatment assignment or
  event status. Each LTFU subject's observed time is drawn from a
  `Uniform(0, t)` distribution, where `t` is their potential event or
  censoring time. Since the LTFU time is always less than `t`, the event
  has not yet occurred at dropout and the subject is right-censored.
  Defaults to zero.

- alternative:

  character. The string specifying the alternative hypothesis, must be
  one of `"greater"` (default), `"less"` or `"two.sided"`. One-sided
  alternatives (`"greater"` and `"less"`) are supported for
  `method = "bayes-surv"` and `method = "bayes-bin"`. All three options
  are supported for `method = "logrank"`, `method = "cox"`, and
  `method = "riskdiff"`. For survival outcomes, `"less"` corresponds to
  the treatment arm having a lower cumulative incidence (i.e., treatment
  is beneficial), and `"greater"` corresponds to the treatment arm
  having a higher cumulative incidence.

- h0:

  single finite numeric null hypothesis value or margin. Default is
  `h0 = 0`. For Bayesian analyses, `h0` must lie in `[0, 1]` for a
  single-arm design and `[-1, 1]` for a two-arm design.

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

  - When `method = "riskdiff"`, `h0` is the null value of
    \\p\_\textrm{treatment} - p\_\textrm{control}\\ and must lie in
    `[-1, 1]`.

  - The argument is ignored for `method = "logrank"` after its
    finite-value validation; the usual equal-survival null hypothesis is
    used.

- Fn:

  vector of values between 0 and 1. Each element is the probability
  threshold to stop at the \\i\\-th look early for futility. If there
  are no interim looks (i.e. `interim_look = NULL`), then `Fn` is not
  used in the simulations or analysis. Set `Fn = 0` to disable futility
  monitoring. The length of `Fn` should be the same as `interim_look`,
  else the values are recycled.

- Sn:

  vector of values between 0 and 1. Each element is the probability
  threshold to stop at the \\i\\-th look early for expected success. If
  there are no interim looks (i.e. `interim_look = NULL`), then `Sn` is
  not used in the simulations or analysis. The length of `Sn` should be
  the same as `interim_look`, else the values are recycled.

- prob_ha:

  scalar value between 0 and 1. Probability threshold of alternative
  hypothesis.

- N_impute:

  integer. Number of imputations for Monte Carlo simulation of missing
  data. An imputed Cox or risk-difference final analysis requires at
  least two.

- N_mcmc:

  integer. Number of posterior samples used by `method = "bayes-surv"`
  and by `method = "bayes-bin"` when `bin_method = "mc"`.

- N_trials:

  integer. Number of trials to simulate.

- method:

  character. For an imputed data set (or the final data set after
  follow-up is complete), whether the analysis should be a log-rank
  (`method = "logrank"`) test, Cox proportional hazards regression model
  Wald test (`method = "cox"`), a fully-Bayesian piecewise-exponential
  analysis (`method = "bayes-surv"`), a Bayesian beta-binomial analysis
  of complete binary outcomes (`method = "bayes-bin"`), or a frequentist
  risk-difference Wald test of complete binary outcomes
  (`method = "riskdiff"`). See Details section.

- imputed_final:

  logical. Should the final analysis (after all subjects have been
  followed-up to the study end) be based on imputed outcomes for
  subjects who were LTFU (i.e. right-censored with time less than
  `end_of_study`)? Default is `FALSE`, which means that the final
  analysis incorporates right-censoring. With `method = "cox"` or
  `method = "riskdiff"`, setting this to `TRUE` analyzes each imputed
  dataset and pools the scalar treatment effects and variances using
  Rubin's rules; this requires `N_impute >= 2`. Imputed final analyses
  remain unavailable for `method = "logrank"`.

- empty_interval:

  character. Policy for empty piecewise-exponential intervals in
  `method = "bayes-surv"` posterior calculations. An empty interval is
  an interval with no exposed subjects in a treatment arm at the
  analysis time. `"propagate"` (the default, matching earlier package
  behavior) copies exposure time and event counts from the nearest
  non-empty interval in the same treatment arm and emits a warning.
  `"prior"` leaves the interval at zero exposure time and zero events,
  so its posterior is driven only by `prior`. `"error"` stops when any
  empty interval is found.

- return_trace:

  logical. Should the compact interim decision trace from every
  simulated trial be retained? The default, `FALSE`, preserves the
  compact output. When `TRUE`, the returned list also contains a
  `traces` data frame with a `trial` column linking each trace row to
  the corresponding row of `sims`.

- ncores:

  positive integer. Number of cores to use for parallel processing.
  Defaults to `1L` (serial execution).

- backend:

  character. Parallel backend. "auto" (the default) uses serial
  execution for `ncores = 1`, the existing fork backend on Unix-like
  platforms, and a PSOCK cluster on Windows. "fork", "psock", and
  "sequential" select a backend explicitly.

- seed:

  optional integer. Seed used to generate independent per-trial
  `"L'Ecuyer-CMRG"` random-number streams. The default, `NULL`, does not
  reset the global RNG state, preserving the usual unseeded simulation
  behavior.

- binary_imputation:

  character. Predictive imputation approach for `method = "bayes-bin"`
  or `method = "riskdiff"`. `"event-time"` (the default) draws a
  conditional piecewise-exponential event time and reduces it to event
  status at `end_of_study`. `"bernoulli"` draws the endpoint status
  directly from its conditional event probability. This argument is
  ignored for time-to-event analysis methods.

## Value

A list containing `sims`, a data frame with one row per simulated trial,
and `call`. When `return_trace = TRUE`, the list also contains `traces`,
a data frame with one row per completed interim look and a `trial`
identifier. See
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
for details of the summary and trace columns.

## Details

This is basically a wrapper function for
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md),
whereby we repeatedly run the function for independent trials (all with
the same input design parameters and treatment effect).

To use multiple cores (where available), the argument `ncores` can be
increased from the default of 1. The default `backend = "auto"` uses
[`pbmcapply::pbmclapply()`](https://rdrr.io/pkg/pbmcapply/man/pbmclapply.html)
on Unix-like platforms and a PSOCK cluster on Windows, where forked
processes are unavailable. Set `backend` explicitly to compare backends
or to require serial execution.

Set `seed` to make `sim_trials()` reproducible. When a seed is supplied,
`sim_trials()` first generates one independent `"L'Ecuyer-CMRG"` stream
for each simulated trial, then each call to
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
runs with its own per-trial stream. This avoids reusing the same
random-number stream across workers when `ncores > 1`, and produces
identical seeded results across supported backends. A seeded call
restores the caller's RNG state on exit. With `seed = NULL`, the
function uses and advances R's current global RNG state.

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
  prior = c(0.1, 0.1),
  block = 2,
  rand_ratio = c(1, 1),
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
