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
  alternative = "less",
  Fn = 0.05,
  Sn = 0.90,
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
#>   ppp_stop_now success_threshold ppp_success_at_max futility_threshold
#> 1            0               0.9                0.2               0.05
interim_result$monte_carlo
#>              estimand successes draws estimate       mcse      lower     upper
#> 1 success_if_stop_now         0    20      0.0 0.00000000 0.00000000 0.1391083
#> 2  success_at_maximum         4    20      0.2 0.08944272 0.07135388 0.4010281
#>   threshold direction point_crossed bound_crossed
#> 1      0.90   greater         FALSE         FALSE
#> 2      0.05      less         FALSE         FALSE
```

At the maximum sample size, equal randomization implies six participants
per arm. Four have already accrued in each arm, so the function
internally adds two potential future participants per arm for the
maximum-sample calculation:

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

## Reviewing and retaining the audit trail

The decision uses the Monte Carlo point estimates. The standard errors
and exact one-sided bounds describe finite-draw uncertainty but do not
change the decision. The one-row trace uses the same schema as a
simulated adaptive-trial trace and can be summarized or plotted with the
existing helpers:

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
#>   success_threshold ppp_success_at_max ppp_success_at_max_mcse
#> 1               0.9                0.2              0.08944272
#>   ppp_success_at_max_lower ppp_success_at_max_upper ppp_success_at_max_draws
#> 1               0.07135388                0.4010281                       20
#>   futility_threshold inner_mc_uncertain_stop_now
#> 1               0.05                           0
#>   inner_mc_uncertain_success_at_max decision                 decision_reason
#> 1                                 0 continue continue_thresholds_not_crossed
#>   empty_interval_fallback_count empty_interval_fallbacks warning_count
#> 1                             0                                      0
#>   warning_messages
#> 1
```

The metadata records the evaluated design, package version, data cut,
fixed time-origin convention, and random-number policy. Supplying `seed`
makes the calculation reproducible while preserving the caller’s
random-number state. With `seed = NULL`, the calculation instead uses
and advances the current R random-number state.
