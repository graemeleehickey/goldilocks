# Applying predictive decisions to an observed interim data cut

``` r

library(goldilocks)
```

[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
simulates a complete adaptive trial.
[`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md)
instead applies the same predictive decision calculation to one observed
interim data cut. It does not generate event, censoring, or enrollment
times, and it does not retain a participant-level copy of the input in
its result.

## Preparing the data cut

The input has one row for every enrolled participant and uses a common
time unit throughout. Calendar-time `enrollment` is measured from first
participant randomization at time zero. Subject-level `time` is measured
from that participant’s randomization.

The status values make the reason for incomplete follow-up explicit:

| `status`   | Interpretation                               |
|------------|----------------------------------------------|
| `event`    | The endpoint event has been observed.        |
| `complete` | Event-free follow-up reached `end_of_study`. |
| `pending`  | The participant remains under follow-up.     |
| `censored` | Follow-up ended early without an event.      |

A blinded data-management team can prepare enrollment, follow-up, event,
and status fields without access to treatment assignments:

``` r

blinded_cut <- data.frame(
  id = 1:8,
  enrollment = 0:7,
  time = c(8, 7, 6, 5, 4, 3, 2, 1),
  event = c(1, 0, 0, 1, 0, 0, 0, 0),
  status = c(
    "event", "pending", "censored", "event",
    "pending", "pending", "pending", "pending"
  )
)

blinded_cut
#>   id enrollment time event   status
#> 1  1          0    8     1    event
#> 2  2          1    7     0  pending
#> 3  3          2    6     0 censored
#> 4  4          3    5     1    event
#> 5  5          4    4     0  pending
#> 6  6          5    3     0  pending
#> 7  7          6    2     0  pending
#> 8  8          7    1     0  pending
```

The calculation itself is arm-specific. An independent unblinded
statistician or service therefore joins the treatment assignment before
running it. The participant-level input need not be distributed with the
aggregate result.

``` r

randomization_export <- data.frame(
  id = 1:8,
  treatment = c(0, 1, 0, 1, 0, 1, 0, 1)
)

interim_cut <- blinded_cut
interim_cut$treatment <- randomization_export$treatment[
  match(interim_cut$id, randomization_export$id)
]
interim_cut <- interim_cut[
  c("id", "treatment", "enrollment", "time", "event", "status")
]
```

## Evaluating the look

Suppose this is the second planned look of a trial with maximum sample
size 12, equal randomization, and ten time units of follow-up per
participant:

``` r

interim_result <- evaluate_interim(
  data = interim_cut,
  data_cut = 8,
  look = 2,
  N_total = 12,
  end_of_study = 10,
  rand_ratio = c(control = 1, treatment = 1),
  method = "logrank",
  alternative = "less",
  Fn = 0.05,
  Sn = 0.90,
  Qn = 1,
  prob_ha = 0.95,
  N_impute = 20,
  seed = 20260831
)

interim_result
#> Goldilocks interim evaluation
#> Look: 2 (N = 8)
#> Decision: continue
#> Predictive success if accrual stops now: 0
#> Predictive success at maximum sample size: 0.2
interim_result$decision
#>   look planned_N calendar_time decision                 decision_reason
#> 1    2         8             8 continue continue_thresholds_not_crossed
#>   ppp_stop_now success_threshold immediate_success_threshold
#> 1            0               0.9                           1
#>   immediate_success_crossed expected_success_crossed ppp_success_at_max
#> 1                     FALSE                    FALSE                0.2
#>   futility_threshold futility_crossed
#> 1               0.05            FALSE
interim_result$monte_carlo
#>              estimand successes draws estimate       mcse      lower     upper
#> 1 success_if_stop_now         0    20      0.0 0.00000000 0.00000000 0.1391083
#> 2  success_at_maximum         4    20      0.2 0.08944272 0.07135388 0.4010281
#>   threshold direction point_crossed bound_crossed
#> 1      0.90   greater         FALSE         FALSE
#> 2      0.05      less         FALSE         FALSE
```

Here `Qn = 1` disables an immediate success claim. Crossing `Sn` stops
accrual for expected success while enrolled participants continue
follow-up. The method, alternative, and thresholds must match the
prespecified trial design.

At the maximum sample size, equal randomization implies six participants
per arm. Four have already accrued in each arm, so the predictive
calculation includes two potential future participants per arm:

``` r

interim_result$diagnostics$target_allocation
#>   control treatment 
#>         6         6
interim_result$diagnostics$current_allocation
#>   control treatment 
#>         4         4
interim_result$diagnostics$potential_accruals
#>   control treatment 
#>         2         2
```

No block size or future randomization sequence is needed. The design
must have an exactly allocable maximum sample size: for example, a 1:2
ratio requires `N_total` to be divisible by three. Observed enrollment
in either arm cannot already exceed its implied maximum.

## Using separate predictive and analysis priors

Suppose a sponsor wants external evidence to inform the prediction of
future events, while the Bayesian survival success test uses a diffuse
prior. Supply both `prior_surv` and `prior_surv_final`. **The latter is
used now, inside each interim predictive replicate, as well as at the
actual final analysis.**

| Step in one interim predictive replicate | Prior used |
|----|----|
| Update observed events and exposure, draw hazards, and generate outstanding outcomes | `prior_surv` |
| Analyze that hypothetical completed trial and compare its posterior probability with `prob_ha` | `prior_surv_final` |

Repeat these steps and calculate the fraction of completed trials that
pass. The same pair of priors is used for current-sample and
maximum-sample predictions. Each completed-data analysis starts from the
specified analysis prior, rather than using the prediction posterior as
a new prior.

For this illustration, interpret the data’s time units as months. The
predictive Gamma priors below have mean hazards of 0.05 and 0.03 per
month for control and treatment. Their rate parameters use
patient-months. These values are illustrative; the actual predictive
prior needs justification from the external evidence and a
prior-predictive assessment.

``` r

predictive_prior <- list(
  control = c(shape = 10, rate = 200),
  treatment = c(shape = 6, rate = 200)
)
analysis_prior <- c(shape = 0.1, rate = 0.1)

bayes_prior_result <- evaluate_interim(
  data = interim_cut,
  data_cut = 8,
  look = 2,
  N_total = 12,
  end_of_study = 10,
  method = "bayes-surv",
  alternative = "less",
  h0 = 0,
  prior_surv = predictive_prior,    # Generates outstanding outcomes
  prior_surv_final = analysis_prior, # Tests each hypothetical completed trial
  Fn = 0.05,
  Sn = 0.90,
  Qn = 1,
  prob_ha = 0.975,
  N_impute = 100,
  N_mcmc = 1000,
  seed = 20260909
)

knitr::kable(
  bayes_prior_result$diagnostics$prior[
    c("stage", "arm", "shape", "rate", "mean_hazard")
  ],
  digits = 3,
  caption = "Gamma priors used in the two parts of this interim calculation."
)
```

| stage   | arm       | shape |  rate | mean_hazard |
|:--------|:----------|------:|------:|------------:|
| interim | control   |  10.0 | 200.0 |        0.05 |
| interim | treatment |   6.0 | 200.0 |        0.03 |
| final   | control   |   0.1 |   0.1 |        1.00 |
| final   | treatment |   0.1 |   0.1 |        1.00 |

Gamma priors used in the two parts of this interim calculation. {.table}

``` r

bayes_prior_result$probabilities
#>              estimand probability threshold direction threshold_crossed
#> 1 success_if_stop_now        0.00      0.90   greater             FALSE
#> 2  success_at_maximum        0.03      0.05      less              TRUE
```

Here the rows labelled `final` describe the analysis prior already used
to test the hypothetical completed datasets. `diagnostics$posterior`
describes the predictive posterior based on the observed interim data.

**Omitting `prior_surv_final` would also use `predictive_prior` in the
success tests. The package does not automatically replace an informative
predictive prior with a weak analysis prior.** The two priors need not
be distinct; the default intentionally uses one prior for both roles.

Use these same two arguments in
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
or
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
when calibrating the design. The small Monte Carlo sizes here
demonstrate the API; precision and operating characteristics must be
assessed for the deployed design. At the actual final analysis,
`prior_surv_final` is also used for any optional imputation of missing
outcomes. Even with a separate analysis prior, predictive borrowing can
affect the stopping decision and sample size.

This example concerns `method = "bayes-surv"`. With
`method = "bayes-bin"`, the completed-data success tests use `prior_bin`
at both interim and final; `prior_surv_final` controls only optional
final imputation.

## Using a prespecified RMST endpoint

For a design that prespecifies RMST as its analysis, supply
`method = "rmst"` and the fixed restriction time. The same illustrative
data cut can demonstrate this separate design, targeting
treatment-minus-control RMST through six time units while follow-up
continues through ten:

``` r

rmst_result <- evaluate_interim(
  data = interim_cut,
  data_cut = 8,
  look = 2,
  N_total = 12,
  end_of_study = 10,
  rand_ratio = c(control = 1, treatment = 1),
  method = "rmst",
  rmst_tau = 6,
  alternative = "greater",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.90,
  Qn = 1,
  prob_ha = 0.95,
  N_impute = 20,
  seed = 20260908
)

rmst_result$decision
#>   look planned_N calendar_time      decision                   decision_reason
#> 1    2         8             8 stop_futility futility_estimate_below_threshold
#>   ppp_stop_now success_threshold immediate_success_threshold
#> 1            0               0.9                           1
#>   immediate_success_crossed expected_success_crossed ppp_success_at_max
#> 1                     FALSE                    FALSE                  0
#>   futility_threshold futility_crossed
#> 1               0.05             TRUE
rmst_result$monte_carlo
#>              estimand successes draws estimate mcse lower     upper threshold
#> 1 success_if_stop_now         0    20        0    0     0 0.1391083      0.90
#> 2  success_at_maximum         0    20        0    0     0 0.1391083      0.05
#>   direction point_crossed bound_crossed
#> 1   greater         FALSE         FALSE
#> 2      less          TRUE         FALSE
rmst_result$metadata$design[c("method", "rmst_tau", "alternative", "h0")]
#> $method
#> [1] "rmst"
#> 
#> $rmst_tau
#> [1] 6
#> 
#> $alternative
#> [1] "greater"
#> 
#> $h0
#> [1] 0
```

Positive effects mean longer event-free time on treatment, so benefit
uses `alternative = "greater"`. Both `h0` and `rmst_tau` use the data’s
time units. The horizon remains fixed across looks and predictive
imputations. Input `status = "complete"` still requires event-free
follow-up to `end_of_study`, even when a participant has already reached
`rmst_tau`. The [RMST
vignette](https://graemeleehickey.github.io/goldilocks/articles/rmst.md)
explains the final analysis and its support and variance requirements.
These small Monte Carlo examples demonstrate the API; the deployed
design needs its own calibration and precision assessment.

## Reviewing and retaining the audit trail

The decision uses the Monte Carlo point estimates. The standard errors
and exact one-sided bounds describe finite-draw uncertainty but do not
change the decision. The interim record contains the same quantities
retained for a simulated adaptive trial and can be summarized or plotted
in the same way:

``` r

summarise_trial_trace(interim_result)
#>   interim_looks_completed last_look last_decision final_N final_post_prob_ha
#> 1                       1         2      continue      NA                 NA
#>   ppp_stop_now ppp_success_at_max warning_count
#> 1            0                0.2             0
interim_result$trace
#>   look planned_N calendar_time active_followup N_enrolled N_treatment N_control
#> 1    2         8             8               5          8           4         4
#>   events_treatment events_control N_pending N_not_enrolled ppp_stop_now
#> 1                1              1         6              4            0
#>   ppp_stop_now_mcse ppp_stop_now_lower ppp_stop_now_upper ppp_stop_now_draws
#> 1                 0                  0          0.1391083                 20
#>   success_threshold immediate_success_threshold immediate_success_crossed
#> 1               0.9                           1                     FALSE
#>   expected_success_crossed ppp_success_at_max ppp_success_at_max_mcse
#> 1                    FALSE                0.2              0.08944272
#>   ppp_success_at_max_lower ppp_success_at_max_upper ppp_success_at_max_draws
#> 1               0.07135388                0.4010281                       20
#>   futility_threshold futility_crossed inner_mc_uncertain_stop_now
#> 1               0.05            FALSE                           0
#>   inner_mc_uncertain_success_at_max decision                 decision_reason
#> 1                                 0 continue continue_thresholds_not_crossed
#>   empty_interval_fallback_count empty_interval_fallbacks warning_count
#> 1                             0                                      0
#>   warning_messages
#> 1
```

The saved audit information records the evaluated design, package
version, data cut, time-origin convention, and simulation seed.
Supplying `seed` makes the predictive calculation exactly reproducible.
