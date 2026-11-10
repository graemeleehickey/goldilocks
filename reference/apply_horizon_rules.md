# Apply a prespecified endpoint schedule to saved trial paths

Applies endpoint horizons and stopping thresholds without further
imputation, fitting, or random-number generation.

## Usage

``` r
apply_horizon_rules(
  paths,
  endpoint_schedule,
  Fn = 0.05,
  Sn = 0.9,
  Qn = 1,
  prob_ha = 0.95,
  return_trace = FALSE
)
```

## Arguments

- paths:

  A complete `goldilocks_horizon_paths` object from
  [`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md).

- endpoint_schedule:

  Data frame with exactly two columns, `N` and `horizon`, containing
  every planned interim enrollment size in increasing order and a
  nondecreasing endpoint horizon at each look. Horizons must be
  available in the saved paths. Use a zero-row data frame with numeric
  columns when there are no interims.

- Fn:

  Futility threshold: one probability repeated at every interim or one
  per scheduled look. `NULL` or zero disables futility; the default is
  `0.05`.

- Sn:

  Expected-success threshold: one probability repeated at every interim
  or one per scheduled look. The default is `0.9`; one disables
  expected-success stopping.

- Qn:

  Immediate-success threshold: one probability repeated at every interim
  or one per scheduled look. It must be at least `Sn` at every look. The
  default, one, disables immediate-success stopping.

- prob_ha:

  A single numeric probability in `[0, 1]` defining success in each
  completed-data analysis. For Bayesian methods this is compared with
  the posterior probability of the alternative; for frequentist methods
  it is compared with `1 - P`. The default is `0.95`.

- return_trace:

  Logical; whether to include a `traces` data frame containing the
  reached interim decisions in the returned list. Defaults to `FALSE`.

## Value

A list compatible with
[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
and
[`summarise_calendar_time()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_calendar_time.md),
containing `sims`, `failures`, and `call`, and optionally `traces`.
Simulation rows additionally identify `endpoint_horizon`,
`final_horizon`, and `followup_horizon`, and the null `margin` used for
the selected endpoint. Treatment-effect estimates from different
endpoint horizons should be reported separately. Overall success is
success under the complete prespecified adaptive design, which may
declare success at different endpoints.

Calendar columns `analysis_ready_time`, `planned_completion_time`, and
`followup_person_time` refer to the selected endpoint. `decision_time`
equals the interim date for immediate success or futility and
`analysis_ready_time` for an expected-success or maximum-N final
analysis. An expected-success stop therefore ends recruitment at
`accrual_stop_time` but determines final success or failure later. Event
and dropout times can precede planned completion. The retained result
column `followup_horizon` equals `final_horizon` for newly simulated
paths; older saved paths retain their original collection duration.
Trace columns identify both prediction data horizons and futility event
and pending-outcome counts; the unprefixed counts refer to success
prediction. Saved posterior diagnostics remain separate by horizon.

## Details

Success uses the scheduled endpoint and its own prediction data cut.
Futility always uses the final endpoint specified when simulating the
paths and its separate data cut, regardless of the scheduled success
endpoint. Immediate success, expected success, and binding futility
retain their existing priority and strict threshold comparisons. The
first stop is terminal. An expected-success stop fixes the endpoint for
that cohort's final analysis; failure there does not trigger extension.

Only calculations needed by a candidate determine whether its trial is
evaluable. A failed maximum-\\N\\ counterfactual makes that comparison
unavailable without discarding a known adaptive outcome. Both the
ordinary futility probability and joint futility with counterfactual
success are available through
[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md).
The latter always refers to maximum \\N\\ at `final_horizon`.

Under a generating scenario where every endpoint null hypothesis is
true, the fraction of trials declaring success at any selected endpoint
estimates the overall type I error. Calibrate the schedule and
thresholds together. If some endpoint hypotheses are false, overall
success includes legitimate findings; assessing false positives then
requires identifying success at a truly null endpoint. A success at an
earlier endpoint does not by itself reject the final-horizon null.

## See also

[`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md),
[`sim_horizon_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_trials.md)
