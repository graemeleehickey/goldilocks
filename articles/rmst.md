# Restricted mean survival time designs

## Effect and hypothesis

`method = "rmst"` compares the area under the treatment and control
Kaplan-Meier curves through one prespecified restriction time,
`rmst_tau`. The estimated effect is **treatment minus control RMST**,
measured in the same time units as follow-up. For an adverse event,
positive values mean additional event-free time on treatment. This
unadjusted, two-arm Wald analysis does not assume proportional hazards.

The default restriction time is `end_of_study`. A shorter `rmst_tau` may
be chosen in advance, but it must be positive and cannot exceed
`end_of_study`. The same value is used for every arm, look, imputation,
and simulated trial. The planned follow-up and imputation horizon remain
`end_of_study`.

For superiority use `h0 = 0` and `alternative = "greater"`. For
non-inferiority allowing a loss of one time unit use `h0 = -1` with the
same alternative. A two-sided difference test uses
`alternative = "two.sided"`. Null values must lie in
`[-rmst_tau, rmst_tau]`.

## A trial with a delayed treatment effect

Suppose time is measured in months and treatment begins reducing the
event hazard after month three. The RMST endpoint summarizes event-free
time through month nine, while follow-up continues through month twelve.

``` r

set.seed(1410)
trial <- survival_adapt(
  hazard_control = c(0.10, 0.10),
  hazard_treatment = c(0.10, 0.04),
  cutpoints = 3,
  N_total = 160,
  lambda = 8,
  interim_look = c(80, 120),
  end_of_study = 12,
  method = "rmst",
  rmst_tau = 9,
  alternative = "greater",
  h0 = 0,
  prob_ha = 0.975,
  N_impute = 100,
  return_trace = TRUE
)
trial$summary[, c("est_final", "post_prob_ha", "N_enrolled", "trial_success")]
#>   est_final post_prob_ha N_enrolled trial_success
#> 1 0.9182006    0.9601273        160         FALSE
trial$trace[, c("planned_N", "ppp_stop_now", "decision")]
#>   planned_N ppp_stop_now decision
#> 1        80         0.50 continue
#> 2       120         0.02 continue
```

`est_final` is the estimated additional event-free time in months
through month nine. `post_prob_ha` is one minus the Wald-test P-value,
not a posterior probability. Final success requires it to be strictly
greater than `prob_ha`. An immediate-success decision instead ends the
trial without a later final analysis, leaving these final-analysis
summaries unavailable.

## Prediction and loss to follow-up

The selected RMST analysis is applied to every predictively completed
trial at both the current and maximum sample sizes.
[`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md)
accepts the same `method` and `rmst_tau` for an observed interim data
cut. Predictive imputation uses the piecewise-exponential model and the
Gamma hazard priors. `N_impute` controls its Monte Carlo resolution;
`N_mcmc` and `binary_imputation` do not change the RMST test.

With `imputed_final = FALSE`, participants lost to follow-up contribute
their observed right-censored data to Kaplan-Meier estimation. Each arm
must have follow-up through the fixed restriction time, or its survival
curve must already have reached zero. A positive survival tail ending
earlier causes an explicit non-estimability error; the package does not
shorten the horizon or extrapolate the tail. Independent censoring is
required for this inference.

With `imputed_final = TRUE`, missing event times are generated from the
final-stage hazard posterior using `prior_surv_final`. The RMST
differences and within-imputation Greenwood variances are pooled using
Rubin’s scalar rules with a t reference distribution. At least two
imputations are needed. The test requires positive total variance,
including after pooling; an individual arm or imputation may contribute
zero variance. Final pooling is distinct from the interim calculation,
which tests each completed trial separately and averages the success
indicators.

## Evaluate the design

Use
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
with the same arguments and inspect both results and `failures`. The
following small run demonstrates the interface; substantially more
trials are needed to assess operating characteristics precisely.

``` r

sims <- sim_trials(
  hazard_control = 0.10,
  hazard_treatment = 0.10,
  N_total = 160,
  lambda = 8,
  interim_look = 80,
  end_of_study = 12,
  method = "rmst",
  rmst_tau = 9,
  alternative = "greater",
  prob_ha = 0.975,
  N_impute = 50,
  N_trials = 20,
  backend = "sequential",
  seed = 1411
)
summarise_sims(sims)
#> # A tibble: 1 × 44
#>   scenario backend     seed n_requested n_analyzed n_failed n_used failure_rate
#>   <chr>    <chr>      <dbl>       <int>      <int>    <int>  <int>        <dbl>
#> 1 1        sequential  1411          20         20        0     20            0
#> # ℹ 36 more variables: failure_rate_mcse <dbl>, failure_rate_mc_lower <dbl>,
#> #   failure_rate_mc_upper <dbl>, power <dbl>, power_mcse <dbl>,
#> #   power_mc_lower <dbl>, power_mc_upper <dbl>, stop_immediate_success <dbl>,
#> #   stop_immediate_success_mcse <dbl>, stop_immediate_success_mc_lower <dbl>,
#> #   stop_immediate_success_mc_upper <dbl>, stop_success <dbl>,
#> #   stop_success_mcse <dbl>, stop_success_mc_lower <dbl>,
#> #   stop_success_mc_upper <dbl>, stop_any_success <dbl>, …
sims$failures
#> [1] trial       error_class message    
#> <0 rows> (or 0-length row.names)
```

A nominal final-test threshold alone does not establish adaptive type I
error control. Prespecify and calibrate the complete stopping rule,
including any immediate-success boundary. Assess equal-survival nulls,
equal-RMST nulls with crossing curves, nonzero margins, delayed
benefits, dropout, and discrepancies between the generating and
predictive hazard models. RMST avoids the proportional-hazards
assumption for the completed-data test; it does not remove assumptions
from prediction or final imputation.

For maintainer validation, `benchmarks/rmst-calibration.R` runs fixed
and adaptive scenarios with failure counts and Monte Carlo intervals,
and `benchmarks/rmst.R` compares completed-data and predictive runtimes.
The technical methods vignette gives the variance and pooling formulas.

## References

Uno H, Claggett B, Tian L, et al. Moving beyond the hazard ratio in
quantifying the between-group difference in survival analysis. *Journal
of Clinical Oncology*. 2014;32:2380-2385.
<https://doi.org/10.1200/JCO.2014.55.2208>.

The numerical reference tests use
[`survRM2::rmst2()`](https://search.r-project.org/CRAN/refmans/survRM2/html/rmst2.html).
Runtime calculations use the existing `survival` dependency; when every
truncated event time is known, the equivalent empirical mean and
Greenwood variance avoid rebuilding a survival fit for each predictive
completion.
