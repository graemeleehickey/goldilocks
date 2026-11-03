# Apply stopping rules to a saved trial-path bank

Reconstructs adaptive trial outcomes from
[`sim_trial_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial_paths.md)
without predictive imputation, model fitting, or random-number
generation.

## Usage

``` r
apply_stopping_rules(
  paths,
  Fn = 0.05,
  Sn = 0.9,
  Qn = 1,
  prob_ha = 0.95,
  return_trace = FALSE
)
```

## Arguments

- paths:

  A complete `goldilocks_paths` object from
  [`sim_trial_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial_paths.md).

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

- Qn:

  A numeric vector of probabilities in `[0, 1]`. Each value is the upper
  predictive-probability threshold for declaring immediate trial success
  at the \\i\\-th look. If there are no interim looks (i.e.
  `interim_look = NULL`), then `Qn` is not used in the simulations or
  analysis. Supply either one value, which is repeated at every interim
  look, or exactly one value per `interim_look`; other lengths are
  rejected. `Qn` must be greater than or equal to `Sn` at every look.
  The default, `1`, disables immediate-success stopping.

- prob_ha:

  A single numeric probability in `[0, 1]` defining success in each
  completed-data analysis. For Bayesian methods this is compared with
  the posterior probability of the alternative; for frequentist methods
  it is compared with `1 - P`. The default is `0.95`.

- return_trace:

  A single logical value indicating whether the interim decision path
  should be returned in addition to the usual final summary. The
  default, `FALSE`, returns the historical one-row data frame. When
  `TRUE`, the result is a `goldilocks_trial` object with summary, trace,
  prior and posterior diagnostics, and call elements.

## Value

A list compatible with
[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md),
[`summarise_calendar_time()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_calendar_time.md),
and simulation plots, containing `sims`, `failures`, and `call`, plus
`traces` when requested. Each simulation row retains its original
`trial` identifier and additionally reports `success_at_max` (the final
outcome if enrollment continued to maximum N without interim stopping)
and `futility_and_success_at_max` (binding futility together with that
counterfactual success). The latter is `FALSE` for non-futility trials
and `NA` for futility trials whose maximum-N final analysis failed.
These quantities use the candidate's `prob_ha` and the bank's final
analysis.

[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
reports the joint probability as `futility_and_success_at_max`, with its
available and unavailable counts. Its ordinary `stop_futility` estimate
is the futility probability under the simulated scenario; under an
alternative this is the requested alternative-scenario futility rate.
The joint probability is not the net power loss relative to a fixed
maximum-N design.

Attributes include the evaluated design, bank provenance
(`path_metadata`), and `stopping_rules`. The bank itself is not copied
into the result. If no trial is evaluable, the function stops with a
`goldilocks_all_trials_failed` error carrying the candidate failure
table.

## Details

At each look, predictive success is the fraction of saved scores
strictly greater than `prob_ha`. Decisions retain the ordered immediate
success (`Qn`), expected success (`Sn`), and binding futility (`Fn`)
rules. The first stopping decision is terminal. Expected success
requires a final analysis of the stopping cohort; immediate success and
futility do not.

Only failures needed by the candidate invalidate a trial. Failures at
unreached looks, unused maximum-cohort predictions, or optional final
analyses remain in the bank but do not remove an evaluable trial. A
missing maximum-N final analysis makes the counterfactual comparison
unavailable, without invalidating an otherwise known adaptive outcome.

Candidate results from the same bank are paired by `trial`. Use paired
uncertainty calculations for differences; marginal Monte Carlo intervals
do not describe uncertainty in a difference. Screen candidate rules on
one bank and validate selected rules on a fresh independent bank.

## See also

[`sim_trial_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial_paths.md),
[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)

## Examples

``` r
paths <- sim_trial_paths(
  hazard_treatment = 0.05, hazard_control = 0.08,
  N_total = 60, interim_look = c(20, 40), end_of_study = 12,
  lambda = 5, method = "bayes-bin", bin_method = "quadrature",
  alternative = "less", N_impute = 5, N_trials = 2, seed = 123
)
candidates <- list(
  original = list(Fn = 0.05, Sn = 0.90, Qn = 1, prob_ha = 0.975),
  stricter = list(Fn = 0.10, Sn = 0.95, Qn = 1, prob_ha = 0.99)
)
results <- lapply(candidates, function(rule) {
  do.call(apply_stopping_rules, c(list(paths = paths), rule))
})
summarise_sims(results)
#> # A tibble: 2 × 50
#>   scenario backend     seed n_requested n_analyzed n_failed n_used failure_rate
#>   <chr>    <chr>      <dbl>       <int>      <int>    <int>  <int>        <dbl>
#> 1 original sequential   123           2          2        0      2            0
#> 2 stricter sequential   123           2          2        0      2            0
#> # ℹ 42 more variables: failure_rate_mcse <dbl>, failure_rate_mc_lower <dbl>,
#> #   failure_rate_mc_upper <dbl>, power <dbl>, power_mcse <dbl>,
#> #   power_mc_lower <dbl>, power_mc_upper <dbl>, stop_immediate_success <dbl>,
#> #   stop_immediate_success_mcse <dbl>, stop_immediate_success_mc_lower <dbl>,
#> #   stop_immediate_success_mc_upper <dbl>, stop_success <dbl>,
#> #   stop_success_mcse <dbl>, stop_success_mc_lower <dbl>,
#> #   stop_success_mc_upper <dbl>, stop_any_success <dbl>, …
```
