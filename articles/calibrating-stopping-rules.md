# Calibrating stopping rules from complete trial paths

``` r

library(goldilocks)
```

Repeated threshold calibration can spend most of its time generating and
analysing the same kinds of predictive completions.
[`sim_trial_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial_paths.md)
performs these calculations once at every planned look.
[`apply_stopping_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_stopping_rules.md)
then evaluates different `Fn`, `Sn`, `Qn`, and `prob_ha` values using
the saved evidence, without further imputation or model fitting.

A bank represents one fixed data-generating scenario and analysis
specification. Separate banks are required for the null, alternatives,
or changes to the priors, analysis method, endpoint horizon, allocation,
maximum sample size, look schedule, or Monte Carlo settings. This
workflow supplements
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md),
which remains useful for evaluating individual designs.

## What is saved

At each interim look, the bank retains the success score from every
predictive completion of the currently enrolled cohort and the maximum
cohort. A score is a posterior probability for a Bayesian analysis or
1-p for a frequentist analysis. For a candidate completed-data threshold
c, the predictive probability is reconstructed as

\widehat P(c) = \frac{1}{M}\sum\_{m=1}^{M} I(s_m \> c),

where M is `N_impute`. Retaining these scores is necessary because
changing `prob_ha` changes the predictive probabilities themselves. A
trace containing only the predictive probabilities for one threshold
cannot support this recalibration.

The bank also retains the final analysis and calendar-time measures for
each possible stopping cohort. If accrual stops at 60 participants for
expected success, the final outcome comes from those 60 participants
after follow-up. The final result from the maximum cohort is reserved
for trials that continue to maximum N and for the separate
counterfactual comparison.

Interim calculations use calendar-masked data. The eventual outcomes
used for cohort final analyses are not made available to the interim
predictions. Dropout and final imputation follow the same rules as
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md).

## Create null and alternative banks

The example has four analyses: interim enrollment looks at 40, 60, and
80 participants, followed by the final analysis of 100 participants
after follow-up. Accrual is six participants per month, follow-up is 12
months, allocation is equal, and there is no dropout. The control event
risk is 30% at 12 months; the alternative treatment risk is 20%. The
primary analysis is a Bayesian beta-binomial comparison of the
treatment-minus-control event risk, using a uniform Beta(1, 1) prior in
each arm and numerical quadrature.

The small simulation counts below demonstrate the workflow. They do not
provide enough precision to establish type I error control or select a
design.

``` r

fixed_design <- list(
  hazard_control = prop_to_haz(0.30, endtime = 12),
  N_total = 100,
  interim_look = c(40, 60, 80),
  lambda = 6,
  lambda_time = NULL,
  end_of_study = 12,
  cutpoints = NULL,
  generation_cutpoints = NULL,
  prior_surv = c(0.1, 0.1),
  prior_surv_final = c(0.1, 0.1),
  block = 2,
  rand_ratio = c(control = 1, treatment = 1),
  prop_loss = c(control = 0, treatment = 0),
  alternative = "less",
  h0 = 0,
  method = "bayes-bin",
  prior_bin = c(1, 1),
  bin_method = "quadrature",
  imputed_final = FALSE,
  N_impute = 10,
  N_mcmc = 1000,
  mc_conf_level = 0.95,
  ncores = 1,
  backend = "sequential"
)

null_args <- c(fixed_design, list(
  hazard_treatment = prop_to_haz(0.30, endtime = 12),
  N_trials = 40, seed = 8201
))
alternative_args <- c(fixed_design, list(
  hazard_treatment = prop_to_haz(0.20, endtime = 12),
  N_trials = 40, seed = 8202
))

banks <- list(
  null = do.call(sim_trial_paths, null_args),
  alternative = do.call(sim_trial_paths, alternative_args)
)
banks$null
#> Goldilocks trial paths
#> 40 trials; 3 interim looks; maximum N = 100 
#> Failed calculations: 0 
#> Memory: 1.2 Mb
```

`N_impute` is part of the implemented stopping algorithm. It is fixed
within each bank and should remain fixed when validating a selected
candidate. `N_mcmc` likewise remains fixed for Monte Carlo Bayesian
analyses; it is unused for the quadrature calculation in this example.
Increasing `N_trials` improves operating-characteristic precision
without changing the decision algorithm.

## Evaluate candidate rules

This grid varies futility, expected-success, and completed-data
thresholds. Immediate success remains disabled with `Qn = 1`. Scalars
are applied at every look; vectors of length three would specify
separate interim boundaries for this design.

``` r

grid <- expand.grid(
  Fn = c(0.05, 0.10),
  Sn = c(0.85, 0.95),
  Qn = 1,
  prob_ha = c(0.975, 0.99)
)
rules <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))
names(rules) <- paste0("rule_", seq_along(rules))

results <- lapply(banks, function(bank) {
  lapply(rules, function(rule) {
    do.call(apply_stopping_rules, c(list(paths = bank), rule))
  })
})
null_summary <- summarise_sims(results$null)
alternative_summary <- summarise_sims(results$alternative)
knitr::kable(cbind(candidate = names(rules), grid), digits = 3)
```

| candidate |   Fn |   Sn |  Qn | prob_ha |
|:----------|-----:|-----:|----:|--------:|
| rule_1    | 0.05 | 0.85 |   1 |   0.975 |
| rule_2    | 0.10 | 0.85 |   1 |   0.975 |
| rule_3    | 0.05 | 0.95 |   1 |   0.975 |
| rule_4    | 0.10 | 0.95 |   1 |   0.975 |
| rule_5    | 0.05 | 0.85 |   1 |   0.990 |
| rule_6    | 0.10 | 0.85 |   1 |   0.990 |
| rule_7    | 0.05 | 0.95 |   1 |   0.990 |
| rule_8    | 0.10 | 0.95 |   1 |   0.990 |

The first stopping decision is terminal. Immediate success takes
precedence over expected success, and either success decision takes
precedence over futility. Expected success stops accrual and requires
the corresponding cohort final analysis. Binding futility is an official
failure even if the maximum-N counterfactual would succeed. All boundary
comparisons are strict; equality does not cross a boundary.

Under the null, the `power` column is the type I error estimate. Its
Monte Carlo interval reflects uncertainty from the finite number of
simulated trials.

``` r

display <- null_summary[c(
  "scenario", "n_used", "power", "power_mc_lower", "power_mc_upper", "mean_N"
)]
names(display)[1:5] <- c("Candidate", "Trials", "Type I error", "Lower", "Upper")
knitr::kable(display, digits = 3)
```

| Candidate | Trials | Type I error | Lower | Upper | mean_N |
|:----------|-------:|-------------:|------:|------:|-------:|
| rule_1    |     40 |            0 |     0 | 0.088 |   61.5 |
| rule_2    |     40 |            0 |     0 | 0.088 |   61.5 |
| rule_3    |     40 |            0 |     0 | 0.088 |   63.5 |
| rule_4    |     40 |            0 |     0 | 0.088 |   63.5 |
| rule_5    |     40 |            0 |     0 | 0.088 |   59.5 |
| rule_6    |     40 |            0 |     0 | 0.088 |   59.5 |
| rule_7    |     40 |            0 |     0 | 0.088 |   61.0 |
| rule_8    |     40 |            0 |     0 | 0.088 |   61.0 |

## Two distinct futility measures

`stop_futility` estimates the probability of binding futility under the
specified data-generating scenario. Under the alternative, it answers
how often that design stops for futility despite the assumed treatment
effect.

`futility_and_success_at_max` estimates the joint probability of binding
futility and success had the same simulated trial continued to maximum N
without interim stopping. The counterfactual uses the same final
analysis and the candidate’s `prob_ha`. Its denominator is all available
joint-outcome comparisons, not only trials stopping for futility. This
is not the net power loss relative to a fixed maximum-N design: some
adaptive trials can succeed when their maximum-N counterfactual would
fail.

``` r

display <- alternative_summary[c(
  "scenario", "power", "stop_futility", "futility_and_success_at_max", "mean_N"
)]
names(display) <- c(
  "Candidate", "Power", "Futility", "Futility and max-N success", "Mean N"
)
knitr::kable(display, digits = 3)
```

| Candidate | Power | Futility | Futility and max-N success | Mean N |
|:----------|------:|---------:|---------------------------:|-------:|
| rule_1    | 0.250 |     0.25 |                      0.050 |   74.0 |
| rule_2    | 0.250 |     0.25 |                      0.050 |   74.0 |
| rule_3    | 0.225 |     0.25 |                      0.050 |   77.5 |
| rule_4    | 0.225 |     0.25 |                      0.050 |   77.5 |
| rule_5    | 0.100 |     0.30 |                      0.025 |   78.5 |
| rule_6    | 0.100 |     0.30 |                      0.025 |   78.5 |
| rule_7    | 0.100 |     0.30 |                      0.025 |   82.0 |
| rule_8    | 0.100 |     0.30 |                      0.025 |   82.0 |

## Audit decisions and failures

Request a candidate trace when inspecting individual decisions. It
includes only reached looks. The original bank retains all hypothetical
looks.

``` r

selected_example <- do.call(apply_stopping_rules, c(
  list(paths = banks$alternative, return_trace = TRUE), rules[[1]]
))
selected_example$traces[selected_example$traces$trial == 1,
  c("trial", "planned_N", "ppp_stop_now", "ppp_success_at_max", "decision")]
#>   trial planned_N ppp_stop_now ppp_success_at_max      decision
#> 1     1        40            0                  0 stop_futility
banks$alternative$trials[[1]]$finals[
  c("planned_N", "post_prob_ha", "analysis_ready_time")
]
#>   planned_N post_prob_ha analysis_ready_time
#> 1        40    0.6331412            18.45936
#> 2        60    0.8647990            21.70470
#> 3        80    0.7664969            24.49464
#> 4       100    0.8664249            27.91432
```

Inspect `bank$failures` for failures by trial, look, and calculation.
One failed predictive replicate invalidates its entire cohort predictive
probability; the implementation does not change the denominator by
dropping that replicate. Later calculations are still attempted. A
candidate only loses a trial when a calculation required by its decision
path is unavailable. Its own `result$failures` records those exclusions.

A failed maximum-N final analysis does not invalidate an otherwise known
early adaptive outcome. For a futility stop it makes
`futility_and_success_at_max` unavailable. The summary reports
`futility_and_success_at_max_n_used` and
`futility_and_success_at_max_n_missing`; its interval uses available
comparisons only. Missing comparisons may be informative and should be
reviewed alongside the estimated rates.

## Compare candidates and validate independently

All candidates within a scenario share the same simulated trials.
Marginal Monte Carlo intervals are suitable for individual operating
characteristics, but differences between candidates require paired
calculations. Align results using the original `trial` identifiers,
particularly if candidates have different failures.

``` r

a <- results$alternative[[1]]$sims
b <- results$alternative[[2]]$sims
paired <- merge(a[c("trial", "trial_success")], b[c("trial", "trial_success")],
  by = "trial", suffixes = c("_a", "_b"))
difference <- as.numeric(paired$trial_success_b) - as.numeric(paired$trial_success_a)
c(difference = mean(difference), mcse = sd(difference) / sqrt(length(difference)))
#> difference       mcse 
#>          0          0
```

Prespecify calibration targets and candidate-selection criteria. For
example, require the upper null Monte Carlo bound to be below the type I
error target, then assess power, futility, enrollment, and duration
among qualifying candidates. A composite null may require several null
banks. Screening many candidates on one bank does not provide
independent confirmation of the selected rule.

The following template generates an independent validation bank after a
candidate has been selected using an adequately sized screening study.
`selected_rule` must contain that candidate’s four threshold settings.
The Monte Carlo settings governing interim decisions remain unchanged.

``` r

validation_args <- null_args
validation_args$N_trials <- 10000
validation_args$seed <- 9201
validation_bank <- do.call(sim_trial_paths, validation_args)
validation_result <- do.call(apply_stopping_rules, c(
  list(paths = validation_bank), selected_rule
))
summarise_sims(validation_result)
```

## Save and process batches

Save the complete bank with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html); exporting only
scores to CSV loses metadata and posterior Monte Carlo counts. For large
studies, generate batches with one explicit seed and consecutive
`trial_offset` values. Each batch is held in memory and can be saved
separately. Keeping all other settings fixed produces the same paths and
trial identifiers as a single large call.

``` r

for (batch in 0:9) {
  batch_args <- null_args
  batch_args$N_trials <- 1000
  batch_args$trial_offset <- batch * 1000L
  bank <- do.call(sim_trial_paths, batch_args)
  saveRDS(bank, sprintf("null-paths-%02d.rds", batch + 1L))
}
```

Evaluate each saved batch separately when memory is limited. Combine
trial outcomes and failure records before calculating overall operating
characteristics; do not average batch confidence limits. Retain the
requested trial count, scenario, seed, candidate rules, and original
trial identifiers when combining results.

The same explicit seed produces identical bank paths across sequential,
fork, and PSOCK backends where supported, and preserves the caller’s
random-number state. The bank uses separate streams for each trial and
calculation stage. Its seeded results therefore need not match
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md),
whose existing draw scheduling is preserved. Rule application itself
consumes no random numbers.
