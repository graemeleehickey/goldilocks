# Simulate and analyze one Goldilocks adaptive trial

Simulate and analyze one Goldilocks adaptive trial

## Usage

``` r
survival_adapt(
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
  empty_interval = c("prior", "propagate", "error"),
  method = "logrank",
  imputed_final = FALSE,
  return_trace = FALSE,
  binary_imputation = c("event-time", "bernoulli"),
  prior_surv_final = prior_surv,
  generation_cutpoints = cutpoints
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

  finite, positive, strictly increasing interior follow-up times
  defining the piecewise-exponential model used for interim posterior
  estimation, predictive imputation, and final analysis. The number of
  interval-specific prior columns must be one greater than the number of
  cutpoints. Default is `NULL`, which uses a constant-hazard analysis
  model.

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

  finite study endpoint, strictly greater than the final value in both
  `cutpoints` and `generation_cutpoints` when either is supplied.

- prior_surv:

  numeric vector, matrix, or named list. Gamma prior for the
  piecewise-exponential hazards used during interim prediction. A
  length-two vector supplies shape and rate and is broadcast across all
  arms and intervals. A `2` by `length(cutpoints) + 1` matrix supplies
  interval-specific values shared by all arms, with shapes in row 1 and
  rates in row 2. For independent arm-specific priors, supply a list
  named `control` and `treatment` in a two-arm design, or `treatment` in
  a single-arm design. Each list element may be a length-two vector or
  an interval-specific matrix. Both arms must be supplied; no values are
  borrowed or filled from the other arm. Rates must use the same time
  unit as event times, exposure, and cutpoints. The default is
  `c(0.1, 0.1)`.

- prior_bin:

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

  length-two positive integer randomization allocation. Name the values
  `control` and `treatment`; either supplied order is accepted and
  normalized internally. A legacy unnamed vector remains accepted in
  `c(control, treatment)` order. Unequal unnamed values produce a
  warning because names may be required in a future major release. See
  [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  for more details.

- prop_loss:

  one or two probabilities. A scalar applies the same LTFU proportion to
  every arm. For a two-arm design, differential attrition can be
  specified with a length-two vector named `control` and `treatment`;
  the supplied order does not matter. Within each arm,
  `ceiling(prop_loss * arm size)` subjects are selected at random
  regardless of event status. Each selected subject's observed time is
  drawn from a `Uniform(0, t)` distribution, where `t` is their
  potential event or censoring time. Since the LTFU time is always less
  than `t`, the event has not yet occurred at dropout and the subject is
  right-censored. Single-arm designs require one probability. Defaults
  to zero.

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

  - When `method = "logrank"`, only `h0 = 0` is supported; this denotes
    the usual equal-survival null. Nonzero values are rejected because
    the standard log-rank statistic does not implement a nonzero effect
    margin.

- Fn:

  vector of values between 0 and 1. Each element is the probability
  threshold to stop at the \\i\\-th look early for futility. If there
  are no interim looks (i.e. `interim_look = NULL`), then `Fn` is not
  used in the simulations or analysis. Set `Fn = 0` to disable futility
  monitoring; `Fn = NULL` has the same effect. Supply either one value,
  which is broadcast to every interim look, or exactly one value per
  `interim_look`. Other lengths are rejected rather than recycled.

- Sn:

  vector of values between 0 and 1. Each element is the probability
  threshold to stop at the \\i\\-th look early for expected success. If
  there are no interim looks (i.e. `interim_look = NULL`), then `Sn` is
  not used in the simulations or analysis. Supply either one value,
  which is broadcast to every interim look, or exactly one value per
  `interim_look`. Other lengths are rejected rather than recycled.

- prob_ha:

  scalar value between 0 and 1. Probability threshold of alternative
  hypothesis.

- N_impute:

  integer. Number of predictive imputations used at each interim look
  and, when requested, for final multiple imputation. The default
  is 500. An imputed Cox or risk-difference final analysis requires at
  least two.

- N_mcmc:

  integer. Number of posterior samples used within each
  `method = "bayes-surv"` and by `method = "bayes-bin"` when
  `bin_method = "mc"`. The default is 1000.

- mc_conf_level:

  probability strictly between 0.5 and 1. Confidence level for one-sided
  exact binomial bounds reported as diagnostics of finite Monte Carlo
  uncertainty. The bounds do not alter completed-data success
  classifications or interim decisions, which use strict point- estimate
  comparisons with `prob_ha`, `Sn`, and `Fn`. The default is 0.95.

- empty_interval:

  character. Policy for empty piecewise-exponential intervals in
  `method = "bayes-surv"` posterior calculations. An empty interval is
  an interval with no exposed subjects in a treatment arm at the
  analysis time. `"prior"` (the default) leaves the interval at zero
  exposure time and zero events, so its posterior is driven only by its
  assigned survival prior. `"propagate"` is a legacy heuristic that
  copies exposure time and event counts from the nearest non-empty
  interval in the same treatment arm and emits a warning. `"error"`
  stops when any empty interval is found.

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

- return_trace:

  logical. Should the interim decision path be returned in addition to
  the usual final summary? The default, FALSE, returns the historical
  one-row data frame. When TRUE, the result is a goldilocks_trial object
  with summary, trace, prior and posterior diagnostics, and call
  elements.

- binary_imputation:

  character. Predictive imputation approach for `method = "bayes-bin"`
  or `method = "riskdiff"`. `"event-time"` (the default) draws a
  conditional piecewise-exponential event time and reduces it to event
  status at `end_of_study`. `"bernoulli"` draws the endpoint status
  directly from its conditional event probability. This argument is
  ignored for time-to-event analysis methods.

- prior_surv_final:

  numeric vector, matrix, or named list. Gamma prior used for
  final-stage piecewise-exponential imputation and, for
  `method = "bayes-surv"`, final analysis. It accepts the same shared or
  arm-specific forms as `prior_surv` and defaults to `prior_surv`,
  preserving the historical behavior.

- generation_cutpoints:

  finite, positive, strictly increasing interior follow-up times
  defining the piecewise-exponential model used to generate event times.
  `hazard_treatment` and `hazard_control` must each have one value per
  resulting interval. Defaults to `cutpoints`, preserving the historical
  behavior in which generation and analysis used one partition.

## Value

With return_trace = FALSE (the default), a data frame containing some
input parameters (arguments) as well as statistics from the analysis,
including:

- `N_treatment`: Number of patients enrolled in the treatment arm.

- `N_control`: Number of patients enrolled in the control arm.

- `est_final`: Treatment effect estimated at the final analysis. The
  final analysis occurs when either the maximum sample size is reached
  and follow-up is complete, or the interim analysis triggered early
  stopping of enrollment/accrual and follow-up for those subjects is
  complete.

- `post_prob_ha`: Posterior probability from the final analysis. If a
  Bayesian method uses `imputed_final = TRUE`, this is calculated for
  each imputed final-analysis dataset and averaged over `N_impute`
  imputations. For an imputed Cox analysis it is \\1 - P\\ from the
  Rubin-pooled Wald test. The same interpretation applies to imputed
  risk-difference analyses. For non-imputed frequentist analyses it is
  \\1 - P\\ from the corresponding test.

- `stop_futility`: Logical indicator of whether the trial stopped early
  for futility.

- `stop_expected_success`: Logical indicator of whether the trial
  stopped early for expected success.

- `stopping_reason`: One of `"expected_success"`, `"futility"`, or
  `"maximum_sample_size"`.

- `accrual_stop_time`: Calendar time of the last enrollment in the
  trial.

- `analysis_ready_time`: Calendar time at which the last enrolled
  subject's observed event or censoring becomes available. This excludes
  external data-cleaning and database-lock delays.

- `planned_completion_time`: Calendar time at which the last enrolled
  subject would complete the full planned follow-up.

- `followup_person_time`: Sum of observed follow-up times across
  enrolled subjects.

- `peak_active_followup`: Largest number of enrolled subjects
  concurrently under follow-up.

Calendar time is measured from the first patient's enrollment at time
zero. Times use the same units as `lambda_time`, `cutpoints`,
`generation_cutpoints`, and `end_of_study`.

The returned object has a `decision_design` attribute containing
`interim_look`, `Fn`, `Sn`, and the Monte Carlo settings. Thresholds in
this metadata are normalized to one value per interim look (and have
length zero when no interim looks are planned). A `prior_design`
attribute contains the fully broadcast Gamma shape, rate, mean hazard,
and standard deviation for every stage, arm, and interval.

Both return forms have an `arguments` attribute containing a named list
of the evaluated argument values, including defaults. It can be saved
with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and supplied to
a later call with `do.call(survival_adapt, attr(result, "arguments"))`.
Its `prop_loss` element is normalized to a named value for every
simulated arm. Its `rand_ratio` element is normalized to `control`,
`treatment` order for two-arm designs. Its `cutpoints` and
`generation_cutpoints` elements retain the analysis and data-generation
partitions, respectively. For `method = "bayes-bin"`, this metadata
explicitly retains the imputation priors (`prior_surv` and
`prior_surv_final`), completed-data analysis prior (`prior_bin`), and
imputation horizon (`end_of_study`). The separate `prior_design`
attribute gives the resolved Gamma parameters by stage, arm, and
interval.

With return_trace = TRUE, a goldilocks_trial object is returned. Its
`summary` element is the same data frame and its `trace` element has one
row per interim look. `prior_diagnostics` contains the resolved interim
and final priors. `posterior_diagnostics` reports observed and effective
sufficient statistics and conjugate posterior parameters by completed
look, arm, and interval. The trace records calendar time, the number of
subjects actively under follow-up, enrollment and observed events by
arm, predictive probabilities, diagnostic Monte Carlo standard errors
and exact bounds, draw counts, thresholds, the decision and reason,
empty-interval fallback diagnostics, and warnings raised during that
look. It deliberately excludes imputed data sets and posterior draws to
keep the output compact.

## Details

Implements the Goldilocks design method described in Broglio et al.
(2014). At each interim analysis, two probabilities are computed:

1.  **The posterior predictive probability of eventual success.** This
    is calculated as the proportion of imputed datasets at the *current*
    sample size that would go on to be success at the specified
    threshold. At each interim analysis this proportion is compared to
    the corresponding element of `Sn`, and if it exceeds the threshold,
    accrual/enrollment is suspended and the outstanding follow-up
    allowed to complete before conducting the pre-specified final
    analysis.

2.  **The posterior predictive probability of final success**. This is
    calculated as the proportion of imputed datasets at the *maximum*
    threshold that would go on to be successful. Similar to above, it is
    compared to the corresponding element of `Fn`, and if it is below
    the threshold, accrual/enrollment is suspended and the trial
    terminated. Typically this would be a binding decision. If it is not
    a binding decision, then one should also explore the simulations
    with `Fn = 0`.

Hence, at each interim analysis look, 3 decisions are allowed:

1.  **Stop for expected success**

2.  **Stop for futility**

3.  **Continue to enroll** new subjects, or if at maximum sample size,
    proceed to final analysis.

At each interim (and final) analysis methods as:

- Log-rank test (`method = "logrank"`). Each (imputed) dataset with both
  treatment and control arms can be compared using a standard log-rank
  test. The output is a *P*-value, and there is no treatment effect
  reported. The function returns \\1 - P\\, which is reported in
  `post_prob_ha`. Whilst not a posterior probability, it can be
  contrasted in the same manner. For example, if the success threshold
  is \\P \< 0.05\\, then one requires `post_prob_ha` \\\> 0.95\\. The
  reason for this is to enable simple switching between Bayesian and
  frequentist paradigms for analysis. When `alternative = "less"` or
  `"greater"`, a one-sided *P*-value is computed from the log-rank
  z-statistic.

- Cox proportional hazards regression Wald test (`method = "cox"`).
  Similar to the log-rank test, a *P*-value is calculated and \\1 - P\\
  is reported in `post_prob_ha`. When `alternative = "two.sided"`, the
  standard two-sided Wald *P*-value is used when `h0 = 0`. For other
  values of `h0`, the Wald test is centered on the specified null log
  hazard ratio. When `alternative = "less"` or `"greater"`, a one-sided
  *P*-value is derived from the Wald z-statistic relative to `h0`. The
  treatment effect (log hazard ratio) is also reported. When
  `imputed_final = TRUE`, the Cox model is fitted separately to each of
  at least two imputed datasets. The log hazard ratios and their
  within-imputation variances are combined using Rubin's rules; the
  pooled Wald test uses Rubin's large-sample degrees of freedom. When
  `imputed_final = FALSE`, the existing single Cox model is fitted
  directly to the observed right-censored data.

- Bayesian absolute difference (`method = "bayes-surv"`). Each imputed
  dataset is used to update the conjugate Gamma prior (defined by
  `prior_surv` at interim looks and `prior_surv_final` at the final
  stage), yielding a posterior distribution for the piecewise
  exponential rate parameters. In turn, the posterior distribution of
  the cumulative incidence function (\\1 - S(t)\\, where \\S(t)\\ is the
  survival function) evaluated at time `end_of_study` is calculated. If
  a single-arm study, then this summarizes the treatment effect, else,
  if a two-armed study, the independent posteriors are used to estimate
  the posterior distribution of the difference. A posterior probability
  is calculated according to the specification of the test type
  (`alternative`) and the value of the null hypothesis (`h0`).

  For piecewise-exponential analyses, an interim or final dataset may
  contain intervals with no exposed subjects in one treatment arm,
  especially when later cutpoints occur after the available follow-up at
  early looks. The `empty_interval` argument controls this case. The
  default, `"prior"`, leaves an empty interval prior-driven, making the
  absence of interval data explicit. The legacy `"propagate"` option
  borrows sufficient statistics from the nearest non-empty interval
  within the same treatment arm. It is operationally stable but
  statistically consequential because adjacent observed data then inform
  the empty interval's posterior. `"error"` is strict and stops the
  simulation or analysis when an empty interval is encountered.

- Bayesian beta-binomial analysis (`method = "bayes-bin"`). Each
  complete or imputed dataset is reduced to binary event outcomes at
  `end_of_study`. A conjugate `Beta(a, b)` prior, specified with
  `prior_bin`, is updated with the number of events and non-events in
  each arm. In a single-arm study, inference is based on the posterior
  event probability. In a two-arm study, inference is based on
  \\p\_\textrm{treatment} - p\_\textrm{control}\\. This posterior
  probability can be calculated using Monte Carlo beta draws
  (`bin_method = "mc"`), a normal approximation (`"normal"`), or
  numerical quadrature (`"quadrature"`). Like the risk-difference test,
  this method requires complete binary outcomes: censored subjects must
  either be followed to `end_of_study`, imputed, or excluded when
  `imputed_final = FALSE`.

  Two equivalent predictive imputation approaches are available through
  `binary_imputation`. With `"event-time"`, the package samples a future
  event time conditional on the available event-free follow-up and then
  records whether it falls by `end_of_study`. With `"bernoulli"`, it
  calculates the same endpoint probability directly. If \\T\\ is the
  observed event-free follow-up, \\T^\*\\ is `end_of_study`, \\S(t)\\ is
  the survival function, and \\H(t)\\ is the cumulative hazard, that
  probability is

  \$\$\Pr(X = 1 \mid T\_\mathrm{event} \> T) = \frac{S(T) -
  S(T^\*)}{S(T)} = 1 - \exp\\-\[H(T^\*) - H(T)\]\\.\$\$

  A Bernoulli outcome is drawn with this probability. For a subject not
  yet enrolled, \\T = 0\\; observed events are retained unchanged.
  Because no precise event time is generated, the imputed `time` is set
  to `end_of_study` and only the binary `event` status is analyzed. Each
  imputation still uses a sampled posterior hazard draw, so uncertainty
  in the piecewise-exponential model is retained.

- Frequentist risk difference (`method = "riskdiff"`). Each complete or
  imputed dataset is reduced to binary event outcomes at `end_of_study`.
  The estimated treatment effect is \\p\_\textrm{treatment} -
  p\_\textrm{control}\\, with an unpooled binomial variance. A Wald test
  compares this estimate with `h0`, and \\1 - P\\ is reported in
  `post_prob_ha`. All three alternatives are supported. Because the test
  requires complete binary outcomes, lost-to-follow-up subjects are
  excluded when `imputed_final = FALSE`. When `imputed_final = TRUE`,
  estimates and within-imputation variances from at least two completed
  datasets are combined using Rubin's rules.

- Imputed final analysis (`imputed_final`). The overall final analysis
  conducted after accrual is suspended and follow-up is complete can be
  analyzed on imputed datasets for Bayesian methods (`"bayes-surv"` and
  `"bayes-bin"`), Cox regression, and the frequentist risk-difference
  analysis, or on the non-imputed dataset. Since the
  imputations/predictions used during the interim analyses assume all
  subjects are imputed (since loss to follow-up is not yet known), it
  would seem most appropriate to conduct the trial in the same manner,
  especially if loss to follow-up rates are appreciable. Note, this only
  applies to subjects who are right-censored due to loss to follow-up,
  which we assume is a non-informative process. For Cox regression the
  final estimates and variances are pooled with Rubin's rules. It cannot
  be used with `method = "logrank"`.

When imputation is involved, either at interim analyses or through
`imputed_final = TRUE`, the package uses a modular impute-then-analyze
procedure. First, the piecewise-exponential model is fitted to the
*observed* time-to-event data and used to complete pending outcomes.
Second, each completed dataset is analyzed using the model selected by
`method`.

For `method = "bayes-bin"`, these are deliberately separate models. The
imputation model has piecewise hazards with Gamma prior `prior_surv` (or
`prior_surv_final` during final imputation), whereas the completed-data
analysis has an event probability at `end_of_study` with Beta prior
`prior_bin`. The Beta prior is not derived from the Gamma prior, and the
two stages are not one joint Bayesian model. Both prior specifications
can therefore affect predictive decisions when outcomes require
imputation. The `"event-time"` and `"bernoulli"` binary-imputation
options use the same piecewise-exponential prediction model and do not
change this separation.

For `method = "bayes-surv"`, the second analysis instead forms a fresh
piecewise-exponential posterior from the completed data and the original
survival prior. For frequentist methods (`"logrank"`, `"cox"`,
`"riskdiff"`), each completed dataset uses a standard test rather than a
posterior. Imputed Cox and risk-difference final analyses pool estimates
and variances using Rubin's rules.

At each interim look, follow-up times are masked (censored) to reflect
the calendar time of the analysis. The package treats enrollment and
randomization as occurring at the same time. Subjects enrolled at the
exact interim boundary have zero follow-up time. These times are clamped
to `.Machine$double.eps` (approximately \\2.2 \times 10^{-16}\\) so that
they contribute negligible but non-zero exposure to the interim
posterior. This affects at most one subject per interim look.

## References

Broglio KR, Connor JT, Berry SM. Not too big, not too small: a
Goldilocks approach to sample size selection. *Journal of
Biopharmaceutical Statistics*, 2014; 24(3): 685–705.

## Examples

``` r
# RCT with exponential hazard (no piecewise breaks)
# Note: the number of imputations is small to enable this example to run
#       quickly on CRAN tests. In practice, much larger values are needed.
survival_adapt(
 hazard_treatment = -log(0.85) / 36,
 hazard_control = -log(0.7) / 36,
 cutpoints = NULL,
 N_total = 600,
 lambda = 20,
 lambda_time = NULL,
 interim_look = 400,
 end_of_study = 36,
 prior_surv = c(0.1, 0.1),
 block = 2,
 rand_ratio = c(control = 1, treatment = 1),
 prop_loss = 0.30,
 alternative = "less",
 h0 = 0,
 Fn = 0.05,
 Sn = 0.9,
 prob_ha = 0.975,
 N_impute = 10,
 N_mcmc = 10,
 method = "bayes-surv")
#>   prob_threshold margin alternative N_treatment N_control N_enrolled N_max
#> 1          0.975      0        less         300       300        600   600
#>   post_prob_ha   est_final ppp_success stop_futility stop_expected_success
#> 1            1 -0.08353039         0.6             0                     0
#>       stopping_reason accrual_stop_time analysis_ready_time
#> 1 maximum_sample_size          27.54104            63.52808
#>   planned_completion_time followup_person_time peak_active_followup
#> 1                63.54104             16439.45                  480
```
