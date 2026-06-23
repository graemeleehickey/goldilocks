# Simulate and execute a single adaptive clinical trial design with a time-to-event endpoint

Simulate and execute a single adaptive clinical trial design with a
time-to-event endpoint

## Usage

``` r
survival_adapt(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = 0,
  N_total,
  lambda = 0.3,
  lambda_time = 0,
  interim_look = NULL,
  end_of_study,
  prior = c(0.1, 0.1),
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
  method = "logrank",
  imputed_final = FALSE
)
```

## Arguments

- hazard_treatment:

  vector. Constant hazard rates under the treatment arm.

- hazard_control:

  vector. Constant hazard rates under the control arm.

- cutpoints:

  vector. Times at which the baseline hazard changes. Default is
  `cutpoints = 0`, which corresponds to a simple (non-piecewise)
  exponential model.

- N_total:

  integer. Maximum sample size allowable

- lambda:

  vector. Enrollment rates across simulated enrollment times. See
  [`enrollment`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  for more details.

- lambda_time:

  vector. Enrollment time(s) at which the enrollment rates change. Must
  be same length as lambda. See
  [`enrollment`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  for more details.

- interim_look:

  vector. Sample size for each interim look. Note: the maximum sample
  size should not be included. For two-arm designs, each interim look
  must be at least the (largest) block size (see `block`), ensuring both
  treatment arms are present at every interim analysis; a smaller look
  could enrol subjects from a single arm only, leaving the interim
  posterior undefined for the missing arm.

- end_of_study:

  scalar. Length of the study; i.e. time at which endpoint will be
  evaluated.

- prior:

  vector. The prior distributions for the piecewise hazard rate
  parameters are each \\Gamma(a_0, b_0)\\, where \\a_0\\ is the shape
  parameter and \\b_0\\ is the rate parameter (i.e., the inverse of the
  scale). This follows R's
  [`rgamma`](https://rdrr.io/r/stats/GammaDist.html) parameterization.
  The same prior is applied to all piecewise intervals and to both
  treatment arms. The default non-informative prior distribution used is
  `Gamma(0.1, 0.1)`, which is specified by setting
  `prior = c(0.1, 0.1)`.

- block:

  scalar. Block size for generating the randomization schedule.

- rand_ratio:

  vector. Randomization allocation for the ratio of control to
  treatment. Integer values mapping the size of the block. See
  [`randomization`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  for more details.

- prop_loss:

  scalar. Overall proportion of subjects lost to follow-up. Subjects are
  selected at random for LTFU regardless of treatment arm or event
  status. Each LTFU subject's observed time is drawn from a
  `Uniform(0, t)` distribution, where `t` is their potential event or
  censoring time. Since the LTFU time is always less than `t`, the event
  has not yet occurred at dropout and the subject is right-censored.
  Defaults to zero.

- alternative:

  character. The string specifying the alternative hypothesis, must be
  one of `"greater"` (default), `"less"` or `"two.sided"`. All three
  options are supported for `method = "bayes"`, `"logrank"`, and
  `"cox"`. The chi-square test (`method = "chisq"`) only supports
  `"two.sided"`. For survival outcomes, `"less"` corresponds to the
  treatment group having a lower cumulative incidence (i.e., treatment
  is beneficial), and `"greater"` corresponds to the treatment group
  having a higher cumulative incidence.

- h0:

  scalar. Null hypothesis value of \\p\_\textrm{treatment} -
  p\_\textrm{control}\\ when `method = "bayes"`. Default is `h0 = 0`. In
  a single-arm design, `h0` is the external benchmark event probability,
  often referred to as a performance goal (PG) or objective performance
  criterion (OPC). The argument is ignored for non-Bayesian analysis
  methods; in those cases the usual method-specific null hypothesis is
  used.

- Fn:

  vector of `[0, 1]` values. Each element is the probability threshold
  to stop at the \\i\\-th look early for futility. If there are no
  interim looks (i.e. `interim_look = NULL`), then `Fn` is not used in
  the simulations or analysis. The length of `Fn` should be the same as
  `interim_look`, else the values are recycled.

- Sn:

  vector of `[0, 1]` values. Each element is the probability threshold
  to stop at the \\i\\-th look early for expected success. If there are
  no interim looks (i.e. `interim_look = NULL`), then `Sn` is not used
  in the simulations or analysis. The length of `Sn` should be the same
  as `interim_look`, else the values are recycled.

- prob_ha:

  scalar `[0, 1]`. Probability threshold of alternative hypothesis.

- N_impute:

  integer. Number of imputations for Monte Carlo simulation of missing
  data.

- N_mcmc:

  integer. Number of samples to draw from the posterior distribution
  when using a Bayesian test (`method = "bayes"`).

- method:

  character. For an imputed data set (or the final data set after
  follow-up is complete), whether the analysis should be a log-rank
  (`method = "logrank"`) test, Cox proportional hazards regression model
  Wald test (`method = "cox"`), a fully-Bayesian analysis
  (`method = "bayes"`), or a chi-square test (`method = "chisq"`). See
  Details section.

- imputed_final:

  logical. Should the final analysis (after all subjects have been
  followed-up to the study end) be based on imputed outcomes for
  subjects who were LTFU (i.e. right-censored with time
  `<end_of_study`)? Default is `TRUE`. Setting to `FALSE` means that the
  final analysis would incorporate right-censoring.

## Value

A data frame containing some input parameters (arguments) as well as
statistics from the analysis, including:

- `N_treatment:`:

  integer. The number of patients enrolled in the treatment arm for each
  simulation.

- `N_control:`:

  integer. The number of patients enrolled in the control arm for each
  simulation.

- `est_final:`:

  scalar. The treatment effect that was estimated at the final analysis.
  Final analysis occurs when either the maximum sample size is reached
  and follow-up complete, or the interim analysis triggered an early
  stopping of enrollment/accrual and follow-up for those subjects is
  complete.

- `post_prob_ha:`:

  scalar. The corresponding posterior probability from the final
  analysis. If `imputed_final` is true, this is calculated as the
  posterior probability of efficacy (or equivalent, depending on how
  `alternative:` and `h0` were specified) for each imputed final
  analysis dataset, and then averaged over the `N_impute` imputations.
  If `method = "logrank"`, `post_prob_ha` is calculated in the same
  fashion, but value represents \\1 - P\\, where \\P\\ denotes the
  frequentist \\P\\-value.

- `stop_futility:`:

  integer. A logical indicator of whether the trial was stopped early
  for futility.

- `stop_expected_success:`:

  integer. A logical indicator of whether the trial was stopped early
  for expected success.

## Details

Implements the Goldilocks design method described in Broglio et al.
(2014). At each interim analysis, two probabilities are computed:

1.  **The posterior predictive probability of eventual success.** This
    is calculated as the proportion of imputed datasets at the *current*
    sample size that would go on to be success at the specified
    threshold. At each interim analysis it is compared to the
    corresponding element of `Sn`, and if it exceeds the threshold,
    accrual/enrollment is suspended and the outstanding follow-up
    allowed to complete before conducting the pre-specified final
    analysis.

2.  **The posterior predictive probability of final success**. This is
    calculated as the proportion of imputed datasets at the *maximum*
    threshold that would go on to be successful. Similar to above, it is
    compared to the corresponding element of `Fn`, and if it is less
    than the threshold, accrual/enrollment is suspended and the trial
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
  standard two-sided Wald *P*-value is used. When `alternative = "less"`
  or `"greater"`, a one-sided *P*-value is derived from the Wald
  z-statistic. The treatment effect (log hazard ratio) is also reported.

- Bayesian absolute difference (`method = "bayes"`). Each imputed
  dataset is used to update the conjugate Gamma prior (defined by
  `prior`), yielding a posterior distribution for the piecewise
  exponential rate parameters. In turn, the posterior distribution of
  the cumulative incidence function (\\1 - S(t)\\, where \\S(t)\\ is the
  survival function) evaluated at time `end_of_study` is calculated. If
  a single arm study, then this summarizes the treatment effect, else,
  if a two-armed study, the independent posteriors are used to estimate
  the posterior distribution of the difference. A posterior probability
  is calculated according to the specification of the test type
  (`alternative`) and the value of the null hypothesis (`h0`).

- Chi-square test (`method = "chisq"`). Each (imputed) dataset with both
  treatment and control arms can be compared using a standard chi-square
  test on the final event status, which discards the event time
  information. The output is a *P*-value, and there is no treatment
  effect reported. The function returns \\1 - P\\, which is reported in
  `post_prob_ha`. Whilst not a posterior probability, it can be
  contrasted in the same manner. For example, if the success threshold
  is \\P \< 0.05\\, then one requires `post_prob_ha` \\\> 0.95\\. The
  reason for this is to enable simple switching between Bayesian and
  frequentist paradigms for analysis. Because the chi-square test cannot
  handle right-censored observations, subjects lost to follow-up are
  excluded from the final analysis when `imputed_final = FALSE`. When
  `imputed_final = TRUE`, LTFU subjects are imputed before the test is
  applied, so all subjects are included.

- Imputed final analysis (`imputed_final`). The overall final analysis
  conducted after accrual is suspended and follow-up is complete can be
  analyzed on imputed datasets (default) or on the non-imputed dataset.
  Since the imputations/predictions used during the interim analyses
  assume all subjects are imputed (since loss to follow-up is not yet
  known), it would seem most appropriate to conduct the trial in the
  same manner, especially if loss to follow-up rates are appreciable.
  Note, this only applies to subjects who are right-censored due to loss
  to follow-up, which we assume is a non-informative process. This can
  be used with any `method`.

When `method = "bayes"` and imputation is involved (either at interim
analyses or via `imputed_final = TRUE`), a two-stage posterior procedure
is used. First, the posterior distribution of the piecewise hazard rates
is estimated from the *observed* data and used to draw imputed event
times for censored subjects. Second, a *new* posterior is estimated from
the combined observed and imputed data, and this posterior is used for
inference. This is consistent with the predictive probability framework
described in Broglio et al. (2014), but users should be aware that the
imputation model's posterior influences the analysis posterior. For
frequentist methods (`"logrank"`, `"cox"`, `"chisq"`), the second stage
uses a standard test rather than a posterior, so this feedback loop does
not arise.

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
 cutpoints = 0,
 N_total = 600,
 lambda = 20,
 lambda_time = 0,
 interim_look = 400,
 end_of_study = 36,
 prior = c(0.1, 0.1),
 block = 2,
 rand_ratio = c(1, 1),
 prop_loss = 0.30,
 alternative = "less",
 h0 = 0,
 Fn = 0.05,
 Sn = 0.9,
 prob_ha = 0.975,
 N_impute = 10,
 N_mcmc = 10,
 method = "bayes")
#>   prob_threshold margin alternative N_treatment N_control N_enrolled N_max
#> 1          0.975      0        less         300       300        600   600
#>   post_prob_ha  est_final ppp_success stop_futility stop_expected_success
#> 1            1 -0.1309505           0             0                     0
```
