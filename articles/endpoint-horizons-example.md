# Extending an endpoint from 12 to 18 months

``` r

library(goldilocks)
```

This example starts with a 12-month endpoint and extends it to 18 months
if the trial reaches enrollment of 400 participants. The change is
prespecified; it is not chosen from the interim results. All
participants have planned outcome collection through 18 months from
enrollment, so earlier recruits can contribute their observed
longer-term outcomes after the switch.

At enrollment looks of 100, 200, and 300, success imputation uses only
available follow-up through 12 months. Futility imputation uses a
separate model fitted to available follow-up through 18 months and
predicts success at maximum enrollment of 500 with the 18-month
endpoint. At the 400-participant look, both predictions use the 18-month
endpoint. If the trial continues, enrollment ends at 500 and the final
analysis follows endpoint ascertainment. There is no
enrollment-triggered interim at 500.

The [technical
companion](https://graemeleehickey.github.io/goldilocks/articles/endpoint-horizons-technical.md)
explains the separate data cuts, predictive models, and stopping
sequence with generic flow charts. The settings below illustrate the
workflow; they have not been selected to satisfy a type I error or power
requirement.

## Specify the endpoint schedule

``` r

endpoint_schedule <- data.frame(
  N = c(100, 200, 300, 400),
  horizon = c(12, 12, 12, 18)
)

knitr::kable(data.frame(
  Analysis = c(paste("Interim at N =", endpoint_schedule$N), "Maximum-N final"),
  Success_endpoint = c(endpoint_schedule$horizon, 18),
  Futility_target = c(rep("N = 500, month 18", 4), "Not applicable")
), col.names = c("Analysis", "Success endpoint (months)", "Futility target"))
```

| Analysis           | Success endpoint (months) | Futility target   |
|:-------------------|--------------------------:|:------------------|
| Interim at N = 100 |                        12 | N = 500, month 18 |
| Interim at N = 200 |                        12 | N = 500, month 18 |
| Interim at N = 300 |                        12 | N = 500, month 18 |
| Interim at N = 400 |                        18 | N = 500, month 18 |
| Maximum-N final    |                        18 | Not applicable    |

If accrual stops for expected success at N=300, that cohort’s final
analysis uses the 12-month endpoint. A subsequent failure at 12 months
does not reopen the trial at 18 months. Futility nevertheless remains a
prediction about the eventual 18-month maximum-N analysis.

## Define the generating and analysis assumptions

The outcome is an adverse event. Under the null scenario, both arms have
a 30% event risk by month 12 and a further 15% conditional event risk
between months 12 and 18 among those event-free at month 12. Under the
alternative, the treatment risks are 22% and 12%, respectively. These
assumptions define piecewise-constant hazards over both intervals from
the outset.

``` r

control_hazard <- c(
  -log1p(-0.30) / 12,
  -log1p(-0.15) / 6
)
treatment_hazard <- c(
  -log1p(-0.22) / 12,
  -log1p(-0.12) / 6
)

risk_table <- data.frame(
  Arm = c("Control / null treatment", "Alternative treatment"),
  Hazard_0_12 = c(control_hazard[1], treatment_hazard[1]),
  Hazard_12_18 = c(control_hazard[2], treatment_hazard[2]),
  Risk_12 = 1 - exp(-12 * c(control_hazard[1], treatment_hazard[1])),
  Risk_18 = 1 - exp(-12 * c(control_hazard[1], treatment_hazard[1]) -
    6 * c(control_hazard[2], treatment_hazard[2]))
)
knitr::kable(risk_table, digits = 4, col.names = c(
  "Arm", "Hazard 0-12", "Hazard 12-18", "12-month risk", "18-month risk"
))
```

| Arm | Hazard 0-12 | Hazard 12-18 | 12-month risk | 18-month risk |
|:---|---:|---:|---:|---:|
| Control / null treatment | 0.0297 | 0.0271 | 0.30 | 0.4050 |
| Alternative treatment | 0.0207 | 0.0213 | 0.22 | 0.3136 |

The 18-month risks are cumulative risks from enrollment, not the
conditional risks in the second interval. All times are measured in
months and hazards in events per participant-month.

Participants are randomized equally in blocks of two and accrue at an
average rate of 10 per month. There is no dropout in the main example.
Predictive imputation uses independent Gamma(0.5, 20) shape-rate priors
for each arm and hazard interval. These have mean hazard 0.025 per
month, with prior contributions of 0.5 events and 20 participant-months
in each interval. The second interval can be entirely prior-driven at
early looks, so its prior is part of the design specification even while
success is assessed at 12 months.

The completed-data analysis is Bayesian binary, comparing
treatment-minus-control event risk against zero in the `"less"`
direction. It uses independent Beta(1, 1) priors and numerical
quadrature. The predictive Gamma model generates missing event times; it
is distinct from this beta-binomial analysis. With no missing final
endpoint outcomes, `prior_surv_final` does not enter the completed
binary analysis, but it is specified explicitly alongside the other
simulation settings.

``` r

fixed_design <- list(
  hazard_control = control_hazard,
  generation_cutpoints = 12,
  cutpoints = 12,
  N_total = 500,
  interim_look = endpoint_schedule$N,
  endpoint_horizons = c(12, 18),
  final_horizon = 18,
  followup_horizon = 18,
  dropout_horizon = 12,
  lambda = 10,
  lambda_time = NULL,
  block = 2,
  rand_ratio = c(control = 1, treatment = 1),
  prop_loss = c(control = 0, treatment = 0),
  prior_surv = c(0.5, 20),
  prior_surv_final = c(0.5, 20),
  empty_interval = "prior",
  method = "bayes-bin",
  prior_bin = c(1, 1),
  bin_method = "quadrature",
  alternative = "less",
  h0 = 0,
  imputed_final = FALSE,
  N_impute = 100,
  N_trials = 500,
  ncores = 1,
  backend = "sequential"
)

rule <- list(Fn = 0.05, Sn = 0.90, Qn = 1, prob_ha = 0.975)
```

`generation_cutpoints` describes the true hazard changes; `cutpoints`
describes the predictive hazard model. They coincide here but need not
coincide in a sensitivity analysis. The null value and beta priors are
the same at both endpoints, so `horizon_settings` is not needed. That
argument can specify different null values or analysis settings by
horizon when the design requires them.

There are 500 trial replicates per scenario and 100 predictive
imputations per look and endpoint. Quadrature does not use `N_mcmc`.
These counts keep the example manageable; the resulting Monte Carlo
intervals show the precision available. Increasing `N_trials` improves
the operating-characteristic estimates. Changing `N_impute` changes the
precision of the stopping calculation itself, so it should be specified
when calibrating and validating the design.

## Simulate all looks and both horizons

We first simulate complete trial paths. Both endpoint horizons are
evaluated at every planned look, including combinations not used by the
current schedule. This permits later comparisons of when to switch
endpoints without repeating model fitting or imputation.

``` r

null_args <- c(fixed_design, list(
  hazard_treatment = control_hazard,
  seed = 18011
))
alternative_args <- c(fixed_design, list(
  hazard_treatment = treatment_hazard,
  seed = 18012
))

paths <- list(
  null = do.call(sim_horizon_paths, null_args),
  alternative = do.call(sim_horizon_paths, alternative_args)
)
paths$alternative
#> Goldilocks endpoint-horizon paths
#> 500 trials; 4 interim looks; horizons: 12, 18 
#> Maximum N: 500 ; final endpoint: 18 
#> Failed calculations: 0 
#> Memory: 42.8 Mb
```

Separate seeds generate independent null and alternative simulations.
Within a scenario, candidate schedules will share the same trial
histories and saved predictive evidence.

``` r

results <- lapply(paths, function(scenario_paths) {
  do.call(apply_horizon_rules, c(list(
    paths = scenario_paths,
    endpoint_schedule = endpoint_schedule,
    return_trace = TRUE
  ), rule))
})

failure_counts <- data.frame(
  Scenario = names(paths),
  Failed_calculations = vapply(paths, function(x) nrow(x$failures), integer(1)),
  Excluded_trials = vapply(results, function(x) nrow(x$failures), integer(1))
)
knitr::kable(failure_counts, row.names = FALSE, col.names = c(
  "Scenario", "Failed path calculations", "Trials excluded by this rule"
))
```

| Scenario    | Failed path calculations | Trials excluded by this rule |
|:------------|-------------------------:|-----------------------------:|
| null        |                        0 |                            0 |
| alternative |                        0 |                            0 |

These are different failure counts. A calculation at an unreached look
can fail without invalidating the candidate’s known trial outcome.
Required failures exclude that trial and should be investigated before
interpreting its operating characteristics. The complete path object
retains the calculation details.

## Examine operating characteristics

[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
reports the probability of success as `power`. Under the joint null used
here, this is the type I error; under the alternative, it is power for
the complete adaptive design, which may establish success at either
endpoint.

``` r

operating <- summarise_sims(results)
operating <- operating[match(names(results), operating$scenario), ]
knitr::kable(operating[c(
  "scenario", "n_used", "power", "power_mc_lower", "power_mc_upper", "mean_N"
)], digits = 3, col.names = c(
  "Scenario", "Trials used", "Success", "MC lower", "MC upper", "Mean N"
))
```

| Scenario    | Trials used | Success | MC lower | MC upper | Mean N |
|:------------|------------:|--------:|---------:|---------:|-------:|
| null        |         500 |   0.036 |    0.023 |    0.056 |  338.4 |
| alternative |         500 |   0.554 |    0.510 |    0.597 |  403.0 |

The limits are 95% Monte Carlo intervals, describing uncertainty from a
finite number of trial replicates. They are not clinical confidence
intervals. This example does not establish type I error control at
0.025: the thresholds and endpoint schedule must be evaluated together,
under relevant joint null scenarios, and selected settings checked using
independent simulations with adequate precision.

Futility is reported in two ways. `stop_futility` is the ordinary
probability of stopping for binding futility; under the alternative, it
is the probability of a futility stop despite the assumed treatment
benefit. `futility_and_success_at_max` is the joint probability of a
futility stop and success had the same simulated trial continued to
N=500 with the 18-month endpoint. It is not conditional on having
stopped for futility.

``` r

knitr::kable(operating[c(
  "scenario", "stop_futility", "futility_and_success_at_max",
  "futility_and_success_at_max_mc_lower", "futility_and_success_at_max_mc_upper",
  "futility_and_success_at_max_n_missing"
)], digits = 3, col.names = c(
  "Scenario", "Futility", "Futility and max-N success", "MC lower", "MC upper",
  "Missing comparisons"
))
```

| Scenario | Futility | Futility and max-N success | MC lower | MC upper | Missing comparisons |
|:---|---:|---:|---:|---:|---:|
| null | 0.724 | 0.000 | 0.000 | 0.008 | 0 |
| alternative | 0.162 | 0.028 | 0.017 | 0.046 | 0 |

The joint outcome is known to be false for trials without a futility
stop. For futility trials, an unavailable maximum-N final analysis makes
it missing. Its interval uses available comparisons only, and the
missing-comparison count should be retained when reporting it.
[Calibrating stopping rules from complete trial
paths](https://graemeleehickey.github.io/goldilocks/articles/calibrating-stopping-rules.md)
explains this distinction further.

## Inspect a trial and its two predictive models

The trace contains only looks reached under the applied rule. We select
an alternative-scenario trial that reached the latest interim observed
in these simulations, so the endpoint schedule is visible. This selected
trace is illustrative; operating characteristics above use all evaluable
trials.

``` r

trace <- results$alternative$traces
trial_id <- trace$trial[which.max(trace$planned_N)]
trial_trace <- trace[trace$trial == trial_id, ]
knitr::kable(trial_trace[c(
  "planned_N", "calendar_time", "success_data_horizon", "futility_data_horizon",
  "ppp_stop_now", "ppp_success_at_max", "decision"
)], digits = 3, col.names = c(
  "N", "Calendar month", "Success cut", "Futility cut", "P current", "P max", "Decision"
))
```

|   N | Calendar month | Success cut | Futility cut | P current | P max | Decision |
|----:|---------------:|------------:|-------------:|----------:|------:|:---------|
| 100 |         11.820 |          12 |           18 |      0.00 |  0.26 | continue |
| 200 |         22.258 |          12 |           18 |      0.01 |  0.12 | continue |
| 300 |         31.316 |          12 |           18 |      0.01 |  0.32 | continue |
| 400 |         39.159 |          18 |           18 |      0.02 |  0.05 | continue |

`P current` predicts success in the enrolled cohort at that look’s
endpoint; `P max` predicts success at N=500 and month 18. The
current-cohort probability is compared with `Qn` and then `Sn`. Futility
is checked only if neither success rule applies. Equality with a
boundary does not trigger stopping.

The saved paths also contain the separate posterior diagnostics. The
following table shows the Gamma updates for this trial at N=300. With a
master cutpoint at month 12, the success model uses one interval per arm
and the futility model uses two. Only follow-up already available on the
look date enters either calculation.

``` r

trial_position <- match(trial_id, vapply(
  paths$alternative$trials, function(x) x$trial, integer(1)
))
look_position <- match(300, fixed_design$interim_look)
evidence <- paths$alternative$trials[[trial_position]]$interims[[look_position]]
posterior_table <- do.call(rbind, lapply(seq_along(evidence), function(j) {
  posterior <- evidence[[j]]$diagnostics$posterior
  data.frame(
    Horizon = fixed_design$endpoint_horizons[j],
    posterior[c("arm", "interval_start", "interval_end", "observed_events",
      "observed_exposure", "posterior_mean_hazard")]
  )
}))
knitr::kable(posterior_table, digits = 3, row.names = FALSE, col.names = c(
  "Cut", "Arm", "From", "To", "Events", "Exposure", "Mean hazard"
))
```

| Cut | Arm       | From |  To | Events | Exposure | Mean hazard |
|----:|:----------|-----:|----:|-------:|---------:|------------:|
|  12 | control   |    0 |  12 |     34 | 1202.142 |       0.028 |
|  12 | treatment |    0 |  12 |     30 | 1250.958 |       0.024 |
|  18 | control   |    0 |  12 |     34 | 1202.142 |       0.028 |
|  18 | control   |   12 |  18 |      7 |  272.790 |       0.026 |
|  18 | treatment |    0 |  12 |     30 | 1250.958 |       0.024 |
|  18 | treatment |   12 |  18 |      5 |  280.367 |       0.018 |

The 12-month success model has no contribution from exposure or events
after month 12. The 18-month futility model includes available
observations from the second interval. [Technical details of
prespecified endpoint
horizons](https://graemeleehickey.github.io/goldilocks/articles/endpoint-horizons-technical.html#two-data-cuts-at-the-same-calendar-date)
gives the general definition, including when a hazard interval spans the
shorter endpoint.

## Compare an earlier switch using the same paths

An alternative schedule switches at N=300. It uses the same interim
enrollment sizes, supported horizons, final endpoint, and analysis
assumptions, so the saved paths are sufficient for this comparison.

``` r

earlier_schedule <- data.frame(
  N = c(100, 200, 300, 400),
  horizon = c(12, 12, 18, 18)
)
earlier_results <- lapply(paths, function(scenario_paths) {
  do.call(apply_horizon_rules, c(list(
    paths = scenario_paths,
    endpoint_schedule = earlier_schedule
  ), rule))
})
comparison <- rbind(
  data.frame(Switch = "N = 400", summarise_sims(results)),
  data.frame(Switch = "N = 300", summarise_sims(earlier_results))
)
comparison <- comparison[order(
  match(comparison$scenario, names(paths)),
  match(comparison$Switch, c("N = 400", "N = 300"))
), ]
knitr::kable(comparison[c(
  "scenario", "Switch", "power", "stop_futility", "mean_N"
)], digits = 3, row.names = FALSE, col.names = c(
  "Scenario", "Switch to 18 months", "Success", "Futility", "Mean N"
))
```

| Scenario    | Switch to 18 months | Success | Futility | Mean N |
|:------------|:--------------------|--------:|---------:|-------:|
| null        | N = 400             |   0.036 |    0.724 |  338.4 |
| null        | N = 300             |   0.032 |    0.726 |  338.0 |
| alternative | N = 400             |   0.554 |    0.162 |  403.0 |
| alternative | N = 300             |   0.554 |    0.162 |  401.8 |

Different `Fn`, `Sn`, `Qn`, and `prob_ha` values can also be passed to
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md)
without repeating simulation. A new generating scenario, prior, null
margin, final horizon, or imputation count requires new paths. Candidate
outcomes are paired by trial, so uncertainty in a difference should use
paired outcomes rather than combine the two marginal Monte Carlo
intervals as if they were independent. Selecting a schedule from these
results calls for independent validation; this comparison does not
select an optimal switch point.

## Direct simulation and continued outcome collection

When evaluating one fixed schedule,
[`sim_horizon_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_trials.md)
skips interim calculations after a terminal stop. The short run below
reproduces the first three alternative-scenario trials from the
complete-path calculation. It retains the same set of horizons,
settings, and seed; only the number of requested trials is reduced.

``` r

direct_args <- alternative_args[setdiff(names(alternative_args), "interim_look")]
direct_args$endpoint_schedule <- endpoint_schedule
direct_args$N_trials <- 3
direct <- do.call(sim_horizon_trials, c(direct_args, rule))

reference <- results$alternative$sims[results$alternative$sims$trial <= 3, ]
stopifnot(identical(direct$sims, reference))
knitr::kable(direct$sims[c(
  "trial", "N_enrolled", "endpoint_horizon", "stopping_reason", "trial_success"
)], col.names = c("Trial", "N enrolled", "Endpoint", "Stopping reason", "Success"))
```

| Trial | N enrolled | Endpoint | Stopping reason     | Success |
|------:|-----------:|---------:|:--------------------|:--------|
|     1 |        500 |       18 | maximum_sample_size | FALSE   |
|     2 |        200 |       12 | expected_success    | TRUE    |
|     3 |        500 |       18 | maximum_sample_size | FALSE   |

The endpoint used at stopping and the planned collection duration are
recorded separately. An expected-success stop at month 12 uses the
12-month endpoint for the primary analysis, but outcome collection is
still planned through month 18 for the enrolled cohort. The table below
summarizes timing in the alternative scenario, conditional on the
endpoint attached to the stopping cohort. It does not average
treatment-effect estimates across endpoints.

``` r

timing <- results$alternative$sims |>
  dplyr::group_by(endpoint_horizon) |>
  dplyr::summarise(
    Trials = dplyr::n(),
    Mean_analysis_ready = mean(analysis_ready_time),
    Mean_endpoint_completion = mean(planned_completion_time),
    Mean_collection_completion = mean(planned_collection_time),
    .groups = "drop"
  )
knitr::kable(timing, digits = 1, col.names = c(
  "Endpoint", "Trials", "Analysis ready", "Planned endpoint completion",
  "Planned collection completion"
))
```

| Endpoint | Trials | Analysis ready | Planned endpoint completion | Planned collection completion |
|---:|---:|---:|---:|---:|
| 12 | 148 | 35.5 | 35.6 | 41.6 |
| 18 | 352 | 65.2 | 65.3 | 65.3 |

These are calendar months from first enrollment. Observed events can
make the endpoint known before its planned completion time. For futility
stops, `analysis_ready_time` describes the stopping cohort’s endpoint
availability, not the date on which futility was declared; that date is
`decision_time`. The [technical
companion](https://graemeleehickey.github.io/goldilocks/articles/endpoint-horizons-technical.html#shared-histories-dropout-and-timing)
defines the endpoint and collection measures.

## Add dropout or save results

For a sensitivity analysis with a 5% independent dropout probability by
month 12, retain `dropout_horizon = 12` and set `prop_loss = 0.05`. This
corresponds to a dropout probability of approximately 7.4% by month 18,
before accounting for events that occur first. The snippet also requests
imputation of missing final outcomes. It requires new simulation paths
because the generating and final-analysis assumptions have changed.

``` r

dropout_args <- alternative_args
dropout_args$prop_loss <- c(control = 0.05, treatment = 0.05)
dropout_args$dropout_horizon <- 12
dropout_args$imputed_final <- TRUE
dropout_args$seed <- 18013
dropout_paths <- do.call(sim_horizon_paths, dropout_args)
dropout_result <- do.call(apply_horizon_rules, c(list(
  paths = dropout_paths, endpoint_schedule = endpoint_schedule
), rule))
```

Saving the complete paths retains the fixed design settings, posterior
diagnostics, predictive scores, and failure records needed for
subsequent calibration. Save the applied schedule and results as well
when retaining evidence for a particular design comparison.

``` r

saveRDS(list(
  paths = paths,
  endpoint_schedule = endpoint_schedule,
  rule = rule,
  results = results
), "endpoint-horizon-calibration.rds")
```

For larger simulation studies, use consecutive `trial_offset` batches
with a fixed seed and otherwise unchanged settings. See [Calibrating
stopping rules from complete trial
paths](https://graemeleehickey.github.io/goldilocks/articles/calibrating-stopping-rules.md)
and [Calibrating prob_ha for type I error
control](https://graemeleehickey.github.io/goldilocks/articles/calibrating-prob-ha.md)
for further calibration guidance.
