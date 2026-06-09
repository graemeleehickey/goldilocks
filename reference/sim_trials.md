# Simulate one or more clinical trials subject to known design parameters and treatment effect

Simulate multiple clinical trials with fixed input parameters, and
tidily extract the relevant data to generate operating characteristics.

## Usage

``` r
sim_trials(
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
  alternative = "two.sided",
  h0 = 0,
  Fn = 0.1,
  Sn = 0.9,
  prob_ha = 0.95,
  N_impute = 10,
  N_mcmc = 10,
  N_trials = 10,
  method = "logrank",
  imputed_final = FALSE,
  ncores = 1L
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
  size should not be included.

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
  p\_\textrm{control}\\ when `method = "bayes"`. Default is `h0 = 0`.
  The argument is ignored when `method = "logrank"` or `= "cox"`; in
  those cases the usual test of non-equal hazards is assumed.

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

- N_trials:

  integer. Number of trials to simulate.

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

- ncores:

  integer. Number of cores to use for parallel processing.

## Value

Data frame with 1 row per simulated trial and columns for key summary
statistics. See
[`survival_adapt`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
for details of what is returned in each row.

## Details

This is basically a wrapper function for
[`survival_adapt`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md),
whereby we repeatedly run the function for a independent number of
trials (all with the same input design parameters and treatment effect).

To use will multiple cores (where available), the argument `ncores` can
be increased from the default of 1. Note: on Windows machines, it is not
possible to use the
[`mclapply`](https://rdrr.io/r/parallel/mclapply.html) function with
`ncores` \\\>1\\.

## Examples

``` r
hc <- prop_to_haz(c(0.20, 0.30), c(0, 12), 36)
ht <- prop_to_haz(c(0.05, 0.15), c(0, 12), 36)

out <- sim_trials(
  hazard_treatment = ht,
  hazard_control = hc,
  cutpoints = c(0, 12),
  N_total = 600,
  lambda = 20,
  lambda_time = 0,
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
  ncores = 1)
```
