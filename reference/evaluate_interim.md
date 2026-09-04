# Evaluate an externally observed interim data cut

Applies the same posterior predictive decision calculation used by
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
to subject-level data observed at one interim look. The function does
not simulate a trial or modify `data`. Participants who could still be
enrolled before the maximum sample size are represented according to
`N_total`, the observed arm counts, and `rand_ratio`.

## Usage

``` r
evaluate_interim(
  data,
  data_cut,
  look,
  N_total,
  end_of_study,
  cutpoints = NULL,
  prior_surv = c(0.1, 0.1),
  prior_bin = c(1, 1),
  bin_method = "mc",
  rand_ratio = c(control = 1, treatment = 1),
  single_arm = FALSE,
  alternative = "greater",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.9,
  prob_ha = 0.95,
  N_impute = 500,
  N_mcmc = 1000,
  mc_conf_level = 0.95,
  empty_interval = c("prior", "propagate", "error"),
  method = "logrank",
  binary_imputation = c("event-time", "bernoulli"),
  seed = NULL
)
```

## Arguments

- data:

  A required data frame with one row per enrolled subject and columns
  `id`, `treatment`, `enrollment`, `time`, `event`, and `status`.
  Treatment is coded `1` for treatment and `0` for control; single-arm
  data use `1`. `enrollment` is measured from first participant
  randomization, which must be zero, and `time` is follow-up from that
  subject's randomization. See Details for the permitted character
  values in `status`.

- data_cut:

  A required single finite, non-negative numeric value giving the
  calendar time of the interim data cut, measured from the same origin
  and in the same units as `enrollment`, `time`, `end_of_study`, and
  `cutpoints`.

- look:

  A required positive integer identifying the prespecified interim look.

- N_total:

  A required positive integer giving the maximum total sample size.

- end_of_study:

  A required single finite, positive numeric value giving the planned
  follow-up time for each subject.

- cutpoints:

  `NULL` (the default), or a numeric vector of finite, positive,
  strictly increasing interior follow-up times defining the
  piecewise-exponential model used for interim posterior estimation,
  predictive imputation, and final analysis. The number of
  interval-specific prior columns must be one greater than the number of
  cutpoints. `NULL` specifies a constant-hazard analysis model.

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

- rand_ratio:

  A length-two positive integer vector giving the control to treatment
  allocation ratio at the maximum sample size. The default is
  `c(control = 1, treatment = 1)`. Name the values `control` and
  `treatment`; either order is accepted. A legacy unnamed vector is
  interpreted as `c(control, treatment)`. The maximum sample size must
  divide exactly according to this ratio. Ignored for single-arm
  designs.

- single_arm:

  A single logical value indicating whether the design has one treatment
  arm and no control arm. The default is `FALSE`.

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

  `NULL`, or a single numeric probability in `[0, 1]` giving the
  threshold for stopping for futility at this look. Futility is declared
  when predictive success at the maximum sample size is strictly less
  than `Fn`. Set `Fn = 0` or `NULL` to disable the maximum-sample
  calculation. The default is `0.05`.

- Sn:

  A single numeric probability in `[0, 1]` giving the threshold for
  stopping for expected success at this look. Expected success is
  declared when predictive success among the currently enrolled
  participants is strictly greater than `Sn`. The default is `0.9`.

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
  strict point- estimate comparisons with `prob_ha`, `Sn`, and `Fn`. The
  default is `0.95`.

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

- binary_imputation:

  A single character string selecting the predictive imputation approach
  for `method = "bayes-bin"`, `method = "riskdiff-wald"`, or
  `method = "riskdiff-fm"`. `"event-time"` (the default) draws a
  conditional piecewise-exponential event time and reduces it to event
  status at `end_of_study`. `"bernoulli"` draws the endpoint status
  directly from its conditional event probability. This argument is
  ignored for time-to-event analysis methods.

- seed:

  `NULL` (the default), or a single non-negative integer used for the
  predictive Monte Carlo calculation. A supplied seed makes the result
  reproducible and leaves the existing random-number state unchanged.
  With `seed = NULL`, the call uses and advances the current
  random-number state.

## Value

An object of class `goldilocks_interim`, containing:

- `decision`: a one-row decision summary;

- `probabilities`: current- and maximum-sample predictive probabilities;

- `monte_carlo`: estimates, standard errors, and diagnostic bounds;

- `diagnostics`: observed status counts, potential accruals, warnings,
  imputation diagnostics, and resolved Gamma prior and posterior
  parameters by arm and interval;

- `trace`: a one-row decision trace compatible with
  [`plot_trial_trace()`](https://graemeleehickey.github.io/goldilocks/reference/plot_trial_trace.md)
  and
  [`summarise_trial_trace()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_trial_trace.md);

- `metadata`: the evaluated design, resolved prior design, package
  version, time-origin, data-cut, and random-number policy. For
  `method = "bayes-bin"`, `metadata$design` retains the normalized
  imputation prior (`prior_surv`), completed-data analysis prior
  (`prior_bin`), and imputation horizon (`end_of_study`).

## Details

The arm-specific maximum enrollment is
`N_total * rand_ratio / sum(rand_ratio)`. Potential future accrual in
each arm is its maximum enrollment minus its observed enrollment.
Randomization block information and concealed future assignment order
are unnecessary: the maximum-sample predictive calculation requires only
the number of potential future participants in each arm.

The `data$status` values distinguish observed follow-up: `"event"` for
an observed endpoint event, `"complete"` for event-free completion of
`end_of_study`, `"pending"` for a subject still under follow-up, and
`"censored"` for permanent early censoring. Pending and censored
outcomes are predictively imputed conditional on `time`.

`Sn` and `Fn` are scalar thresholds for this look. Expected success is
declared when the estimated probability of completed-data success among
the current participants is strictly greater than `Sn`. Futility is
declared when the corresponding maximum-sample probability is strictly
less than `Fn`. Set `Fn = 0` or `NULL` to disable futility. Exact
one-sided Monte Carlo bounds are returned as diagnostics and do not
drive either decision.

The function requires treatment assignments to perform the arm-specific
posterior and completed-data analyses. In a blinded trial, an
independent unblinded statistician or service should join the treatment
assignments, run this function, and return the aggregate decision
without distributing the subject-level input.

## Examples

``` r
interim_data <- data.frame(
  id = 1:6,
  treatment = c(0, 1, 0, 1, 0, 1),
  enrollment = c(0, 1, 2, 3, 4, 5),
  time = c(6, 5, 4, 3, 2, 1),
  event = c(1, 0, 0, 1, 0, 0),
  status = c("event", "pending", "pending", "event", "pending", "pending")
)

evaluate_interim(
  data = interim_data,
  data_cut = 6,
  look = 1,
  N_total = 10,
  end_of_study = 12,
  rand_ratio = c(control = 1, treatment = 1),
  alternative = "less",
  N_impute = 5,
  seed = 2026
)
#> Goldilocks interim evaluation
#> Look: 1 (N = 6)
#> Decision: stop_futility
#> Predictive success if accrual stops now: 0
#> Predictive success at maximum sample size: 0
```
