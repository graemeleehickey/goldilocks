# goldilocks 0.6.0.9000

## Improvements

* For fixed-time binary endpoints, each completed predictive replicate is now
  summarized by the event count and sample size in each arm. These are the
  sufficient statistics for the risk-difference and beta-binomial analyses.
  Identical sufficient statistics are evaluated once for deterministic
  analyses. With `bin_method = "mc"`, every replicate continues to receive
  fresh beta-posterior draws, even when its counts match an earlier replicate.
  Predictive probabilities and stopping rules are unchanged. Interim results
  report how often identical completed-data summaries occurred and were reused
  (#91).
* Additional statistical validation now covers piecewise-exponential event and
  conditional-imputation distributions, conjugate Gamma posterior moments,
  frequentist type I error and confidence-interval coverage, and invariance to
  the choice of time unit. Longer non-CRAN reference calculations also check
  predictive probabilities close to the futility and expected-success
  thresholds.
  Acceptance limits reflect Monte Carlo uncertainty and identify the estimand,
  target, estimate, and uncertainty when a check fails.
* Function documentation now states argument types, accepted lengths and
  ranges, defaults, and required inputs more consistently. Titles,
  descriptions, analysis-method explanations, and return values have been
  revised to emphasize their statistical interpretation.
* All `N_impute` completions at an interim look are now generated together from
  the corresponding posterior hazard draws. The statistical imputation model
  and completed-data analyses are unchanged, but repeated preliminary work is
  avoided. Imputations are generated in posterior-draw order: current treatment,
  current control, future treatment, then future control. Because all
  imputations now precede the completed-data analyses, a given seed can produce
  different Bayesian results than earlier package versions. Repeated runs and
  supported serial or parallel calculations remain reproducible (#37).
* Bayesian survival calculations now use analysis-interval widths that always
  extend from time zero to `end_of_study`. They are not shortened when no
  participant has yet reached the maximum follow-up time at an interim or final
  analysis. This preserves the specified cumulative-hazard estimand while
  avoiding repeated calculation of the same interval widths (#90).
* Repeated Bayesian survival analyses of completed predictive datasets now work
  directly from arm-by-interval event counts and exposure times. These are the
  sufficient statistics for the conjugate Gamma posterior. The posterior
  distribution, treatment-effect calculation, empty-interval policy, and
  reported results are unchanged (#89).
* Documentation now states explicitly that Bayesian binary prediction uses a
  piecewise-exponential Gamma model to impute pending outcomes and a separate
  beta-binomial model to analyze completed endpoint statuses. These are
  separate models rather than one joint Bayesian model. Retained analysis
  information identifies both priors and the imputation horizon (#64).
* `prior_surv` and `prior_surv_final` now accept independent arm-specific Gamma
  priors supplied as lists named `control` and `treatment` (or `treatment` for a
  single-arm design). Each arm may use one shape-rate pair across intervals or
  an interval-specific matrix. Existing shared priors remain supported and
  preserve seeded results. Individual trial traces and externally evaluated
  interim looks report prior and posterior parameters by arm and interval.
  Priors remain independent; there is no hierarchical or implicit borrowing
  between arms (#88).
* New `evaluate_interim()` applies the same posterior predictive decision rule
  used by `survival_adapt()` to an observed interim data cut. Follow-up status
  distinguishes events, completed follow-up, pending outcomes, and early
  censoring. The potential number of future participants in each arm is derived
  from `N_total`, observed enrollment, and `rand_ratio`; block details and
  unobserved future assignments are not required. Results include predictive
  probabilities, Monte Carlo uncertainty, allocation and imputation summaries,
  a one-look decision trace, and the design information needed for audit (#81).
* `survival_adapt()` and `sim_trials()` now allow event-time generation and the
  predictive model to use different piecewise-exponential partitions.
  `generation_cutpoints` defines the intervals for `hazard_treatment` and
  `hazard_control`; `cutpoints` defines the intervals for posterior estimation,
  predictive imputation, Bayesian survival analysis, and interval-specific
  priors. `generation_cutpoints` defaults to `cutpoints`, preserving existing
  calls and seeded results. In direct calls to `sim_comp_data()`, the previous
  named argument `cutpoints` must now be written as `generation_cutpoints`
  (#95).
* `rand_ratio` and `randomization(allocation =)` now accept vectors named
  `control` and `treatment` in either order. Unnamed vectors continue to mean
  `c(control, treatment)`; unequal unnamed allocations warn that names may be
  required in a future major release. Retained analysis information records
  the control and treatment allocation weights in that order.
* `prop_loss` now accepts separate loss-to-follow-up proportions for the
  control and treatment arms through a named vector. A single value continues
  to apply the same proportion to both arms. Attrition is sampled separately
  within each randomized arm, and single-arm designs continue to require one
  value (#69).
* Piecewise event and enrollment times now use the survival counting-process
  convention `(start, stop]`: a time exactly at a cutpoint belongs to the
  interval ending at that cutpoint. Posterior sufficient statistics and
  piecewise enrollment calculations use this convention consistently. PWEALL
  continues to generate continuous event times using `[start, stop)` hazard
  intervals; the two conventions differ only at individual boundary points,
  which have probability zero under a continuous distribution (#70).
* Results from `survival_adapt()` and `sim_trials()` now retain all evaluated
  design and analysis settings, including defaults. These values can be saved
  and supplied to the same function with `do.call()` to reproduce or modify an
  analysis without reconstructing the original call (#86).
* `summarise_sims()` now reports the numbers of requested, successfully
  analyzed, failed, and included trials. Every estimated probability and mean
  sample size is accompanied by a 95% Monte Carlo interval. Wilson intervals
  remain informative when no trials, or all trials, have the outcome. The
  optional `max_mcse` criterion warns when a scenario has insufficient Monte
  Carlo precision and distinguishes simulation uncertainty from clinical or
  model uncertainty (#68).
* New `summarise_calendar_time()` reports operating characteristics for trial
  duration, accrual, analysis readiness, person-time, peak concurrent follow-up,
  and interim-look timing. Calendar time begins at first enrollment, and the
  requested, analyzed, and failed trial counts remain explicit (#87).

## Bug fixes

* Parallel simulations using PSOCK now load `goldilocks` and its compiled
  log-rank calculation correctly on Windows, including when vignettes are built
  from the source package.
* `summarise_sims()` now calculates the probability of reaching the maximum
  sample size from the trial-level stopping indicators. Previously this
  probability could be incorrect because the futility probability and the
  trial-level stopping indicator used the same name during the calculation
  (#68).
* `sim_trials()` now limits the number of workers to the number useful for the
  requested trials and starts parallel calculation only when there are at least
  two trials per worker. This avoids parallel setup for fewer than four trials.
  Worker processes are closed after an error, and timing summaries separate
  setup time from simulation time (#93).
* A failure in one simulated trial no longer discards the other trials in the
  same `sim_trials()` call. Successful trials are retained, failed trials are
  listed separately, and one warning gives the number affected. If every trial
  fails, the calculation still stops and provides the failure details. Summary
  denominators report requested, analyzed, and failed trials explicitly (#72).
* Repeated Cox analyses use a faster calculation when it is compatible
  with the installed `survival` package and otherwise use
  `survival::coxph()`. Both calculations target the same Wald test. Singular
  fits, analyses with no events, and convergence warnings now produce specific
  non-estimability errors (#77).
* `summarise_sims()` now accepts a complete `sim_trials()` result as well as a
  trial-level results table. Named or grouped scenario identifiers are retained,
  together with the seed, parallel method, failures, and design settings.
  Requested, analyzed, and failed trial counts remain explicit so failures
  cannot silently change the denominator (#78).
* Interim Monte Carlo decisions use strict point-estimate comparisons at both
  the completed-data posterior layer and the outer predictive-imputation layer,
  preserving the Goldilocks stopping rule. The defaults are now `N_impute =
  500` and `N_mcmc = 1000`, and decision traces report estimates, diagnostic
  Monte Carlo standard errors and exact bounds, draw counts, uncertainty
  counts, and the decision reason (#60).
* Empty piecewise-exponential intervals are now prior-driven by default. The
  former `empty_interval = "propagate"` behavior remains available only as an
  explicit legacy heuristic. To reproduce earlier results, set it explicitly;
  otherwise review the interval-specific prior because it now supplies all
  information for an empty interval. Decision traces identify every interval
  handled by either policy (#62).
* Unseeded PSOCK simulations now derive independent random-number streams for
  individual trials from the caller's random-number state. Resetting that state
  reproduces the calculation, whereas consecutive calls advance rather than
  reuse the same streams. Results record the parallel method and seed policy
  (#71).
* Log-rank analyses now reject nonzero `h0` values instead of silently testing
  the equal-survival null. Cox margins remain on the log-hazard-ratio scale,
  while Bayesian and risk-difference margins remain on their documented event-
  probability scales (#73).
* Futility and expected-success thresholds may now be supplied either as one
  value applied to every interim look or as one value per look. Incorrect
  lengths produce an error that states the expected and observed numbers, and
  the resolved thresholds are retained with the design information (#74).
* Exported piecewise-exponential utilities now validate draw counts and
  observed times consistently, with explicit zero-draw and zero-length
  behavior. Related binary-outcome validation now rejects missing and malformed
  event indicators before analysis (#75).
* Cumulative-hazard/event-probability conversions now use stable `expm1()` and
  `log1p()` transformations, retaining precision near zero and defining exact
  behavior at zero and infinite boundaries (#76).
* Randomization errors now state correctly that every block size must be a
  multiple of the sum of the two-arm allocation weights and include the
  observed values (#79).

# goldilocks 0.6.0

## Improvements

* New `plot_enrollment()` displays the expected enrollment curve, simulated
  trial trajectories, interim looks, and the maximum sample size for constant
  or piecewise enrollment rates. Design settings may be supplied directly or
  taken from a `survival_adapt()` or `sim_trials()` result.
* Survival priors are now specified through `prior_surv` and
  `prior_surv_final`, and the Bayesian binary prior through `prior_bin`.
  `prior_surv_final` defaults to `prior_surv` but may specify a different Gamma
  prior for final imputation and analysis. A survival prior may be one
  shape-rate pair shared across intervals or an interval-specific matrix. The
  ADVENT vignette illustrates interval- and stage-specific priors (#59).
* **Major change from goldilocks 0.5.0 and earlier:** `enrollment()` now
  simulates exact continuous arrival times from a piecewise-constant Poisson
  process. This handles fractional changes in enrollment rate without dividing
  follow-up into unit-time bins. `lambda_time = NULL` represents a constant
  rate, and every positive change time requires one additional value of
  `lambda`. Time zero is the first-patient-in origin (#55). Existing calls must
  replace `lambda_time = 0` with `NULL` and remove the leading zero from
  piecewise schedules. Results obtained with a fixed seed will differ from
  earlier versions, and estimated operating characteristics may change,
  particularly for low rates or fractional change times.
* The frequentist binary `method = "chisq"` analysis has been replaced by
  `method = "riskdiff"`. It estimates the treatment-minus-control event-risk
  difference and supports one- or two-sided Wald tests against `h0`. With
  `imputed_final = TRUE`, estimates and within-imputation variances are combined
  using Rubin's rules. The previous mean-of-`1 - P` rule is no longer used.
  Imputed log-rank final analyses remain unsupported (#56).
* Cox regression now supports `imputed_final = TRUE`. Log hazard ratios and
  within-imputation variances are combined using Rubin's rules, and
  `post_prob_ha` reports `1 - P` from the pooled Wald test. Analyses without
  final imputation are unchanged (#56).
* Bayesian survival and Monte Carlo beta-binomial analyses now use `N_mcmc` for
  the number of posterior draws; the separate `bin_N` argument has been
  removed. Results retain the posterior mean treatment effect rather than all
  individual posterior draws.
* `survival_adapt()` can now retain an interim decision trace containing
  predictive probabilities, decision thresholds, stopping outcomes, arm-level
  sample sizes, and relevant warnings. New functions summarize and plot traces
  from individual trials and simulation studies (#57).
* `plot_sim_stopping()` now provides marginal, conditional, cumulative, and
  flowchart summaries of stopping outcomes by enrolled sample size. Plot
  denominators are stated explicitly, and reached looks with no stopping events
  are retained when simulation traces are available.
* New `plot_sim_ocs()` displays success probability, stopping probability, and
  expected sample size across treatment-effect scenarios.
* `sim_trials(return_trace = TRUE)` retains interim results across simulated
  trials, and new `plot_sim_decisions()` displays predictive probabilities and
  decision regions at each interim look.
* `sim_trials()` now supports reproducible parallel simulation on Windows,
  macOS, and Linux. The `backend` argument selects the parallel method, with
  `backend = "auto"` choosing an appropriate method for the operating system.
  Setting `seed` leaves the caller's random-number state unchanged (#40).
* Bayesian binary calculations now use a single consistent beta-binomial
  analysis for the Monte Carlo, normal-approximation, and quadrature methods.
* The piecewise-exponential Bayesian survival method is now named
  `method = "bayes-surv"`, distinguishing it from the new beta-binomial
  `method = "bayes-bin"` for complete binary outcomes.
* Binary endpoint analyses may set `binary_imputation = "bernoulli"` to draw
  event status from the conditional piecewise-exponential event probability.
  Conditional event-time imputation remains the default (#21).
* Repeated Cox model analyses are faster while retaining the same estimated log
  hazard ratio and Wald test.
* For Cox analyses, `h0` is now the null log hazard ratio. Non-inferiority may
  therefore be specified as `h0 = log(margin)`.
* Performance checks now cover posterior probabilities, posterior sampling,
  imputation, and representative `survival_adapt()` simulations (#42).
* Treatment assignment is described consistently as `treatment = 1` for the
  treatment arm and `treatment = 0` for the control arm (#35).
* The `seed` argument to `sim_trials()` creates independent random-number
  streams for individual trials, giving reproducible results with one or
  several cores (#41).
* Piecewise-exponential event probabilities and posterior sufficient
  statistics are now calculated directly from cumulative hazards, events, and
  exposure times, reducing repeated work in simulation studies (#34).
* After predictive imputation, Bayesian survival analyses update a fresh
  completed-data posterior from arm-by-interval events and exposure. This
  retains the documented two-stage posterior procedure while reducing
  simulation time (#38).
* Enrollment, randomization, and assembly of simulation summaries have been
  made more efficient (#43, #44).

## Bug fixes

* Piecewise Bayesian analyses now work when `N_mcmc = 1`.
* `sim_trials()` now uses the same default alternative hypothesis and futility
  threshold as `survival_adapt()`, so omitted arguments define the same design
  (#47).
* Interim looks must be strictly increasing, preventing duplicate or
  non-chronological analyses (#48).
* Piecewise hazards must be finite and non-negative. A finite `maxtime` is
  required when the final hazard is zero (#49).
* `h0` must be one finite value and must satisfy the relevant probability bounds
  for Bayesian analyses (#50).
* `cutpoints` now contains only the interior hazard-change times. `NULL`
  represents a constant hazard, and every cutpoint requires one additional
  hazard rate (#58). Cutpoints must be finite, positive, and strictly
  increasing, with the study endpoint after the final cutpoint (#52).
* `enrollment()` now checks the full schedule, including the sample size,
  positive enrollment rates, and strictly increasing change times, before
  generating data (#53).
* `sim_trials()` now defaults to `ncores = 1L`, avoiding unexpected parallel
  calculation when the argument is omitted (#54).
* Designs without interim looks now use the documented final success rule and
  ignore the unused interim success and futility thresholds.
* Invalid probabilities, priors, and sample-size arguments are now identified
  before simulation begins.
* `sim_trials()` uses one core when the available number of cores cannot be
  determined.
* Log-rank and Cox analyses now give a specific non-estimability error when the
  test statistic cannot be calculated reliably.
* Two-arm randomization now requires two positive integer allocation weights;
  invalid schedules such as `c(0, 1)` are rejected.
* `prop_to_haz()` now rejects cumulative event probabilities that are
  non-finite, outside `[0, 1)`, or decreasing over time.
* `pwe_impute()` now rejects a `maxtime` earlier than the participant's observed
  follow-up time.
* When `Fn = 0`, futility monitoring remains disabled throughout the trial
  (#33).
* `empty_interval` now determines whether an interval with no observed exposure
  borrows information from an adjacent interval, uses its prior alone, or
  produces an error.

## Documentation

* Updated the package title, description, help landing page, README, and
  architecture overview to cover survival and fixed-time binary endpoints, all
  five final-analysis methods, and the available plotting functions.
* Corrected the technical-methods vignette's section numbering and updated its overview to cover the package's binary endpoint methods.
* Expanded the simulation, design, and technical vignettes with examples of `plot_enrollment()`, `plot_sim_ocs()`, `plot_sim_stopping()`, and `plot_sim_decisions()`.
* Revised function help pages to use clearer tables, links, mathematical
  notation, and formatting.
* Added a new vignette, "Technical details of the Goldilocks design", documenting the design notation, piecewise-exponential event-time model, Gamma posterior updating, posterior predictive probabilities, interim decision rules, final analysis options, and simulation-based calibration.
* Added a new vignette, "Bayesian binary outcome designs", documenting `method = "bayes-bin"` for two-arm and single-arm complete binary endpoint analyses.
* Added a new vignette, "ADVENT: a published Goldilocks design", showing how the published ADVENT pulsed field ablation trial maps to `method = "bayes-bin"` with beta-binomial non-inferiority endpoints, adaptive sample-size thresholds, and pre-computed simulation examples.
* Renamed vignettes for consistency: "Two-arm randomized trials", "Bayesian piecewise-exponential designs", "Single-arm designs with a performance goal", and "Package architecture".
* Clarified that `goldilocks` treats enrollment time and randomization time as the same time point in its time-to-event simulations.
* Clarified that the single-arm external benchmark `h0` is often referred to as a performance goal (PG) or objective performance criterion (OPC).
* Expanded technical documentation for the Goldilocks predictive-probability algorithm, including notation for final analysis quantities, operating characteristics, and method-specific decision rules.
* Clarified the `Fn` documentation to state that `Fn = 0` disables futility monitoring.
* Added light, dark, and automatic display modes to the documentation website.
* Updated the README with CRAN status, download, and license information.

# goldilocks 0.5.0

## Improvements

* Cox and log-rank analyses now support one-sided alternatives through
  `alternative = "greater"` or `"less"`. The chi-square analysis remains
  two-sided (#20).
* For a non-imputed chi-square final analysis, participants lost to follow-up
  are now excluded. Previously they were counted as non-events, diluting the
  event rate and biasing the comparison (#22).
* `posterior()` now warns when a piecewise interval contains no participants and
  information is carried forward from an adjacent interval.
* In a two-arm design, every interim look must be at least as large as the
  randomization block. This prevents an interim analysis before both arms can
  be represented.
* Package maintenance removed calculations inherited from the deprecated
  `fastlogranktest` package that were no longer used and adopted the base R pipe
  in package calculations.

## Bug fixes

* For a piecewise interval with no exposure, information is now carried from an
  adjacent interval within the same treatment arm. Previously, an empty first
  interval in the treatment arm could incorrectly use control-arm information,
  contaminating the treatment posterior. The expected arms must also be present
  in the supplied data; an absent arm now produces an informative error.
* Designs with `interim_look = NULL` now proceed directly to the final analysis.
* `impute_data()` now continues to work when the supplied dataset contains
  additional columns (#26).
* Multi-element randomization block schedules now cycle correctly rather than
  producing a missing next block (#31).
* Cox model estimates are now extracted reliably from the fitted model results
  (#29).
* One-sided Cox and log-rank P values no longer depend on which other packages
  the user has attached.
* Enrollment exactly at a piecewise rate-change time now uses the intended
  enrollment rate (#28).
* Interim analyses no longer add a small value to every survival time. Only a
  participant with exactly zero follow-up is adjusted by machine precision when
  required for the piecewise calculation (#24).
* Removed an unused loss-to-follow-up calculation.

## Documentation

* Removed the `est_interim` element from the `survival_adapt()` return-value documentation. This field was documented but never computed or returned.
* Documented the two-stage posterior procedure used when `method = "bayes-surv"` with imputation, clarifying that the imputation model's posterior influences the analysis posterior (#27).
* Clarified `prop_loss` parameter documentation, explaining that LTFU times are drawn from `Uniform(0, t)` and that the event has not yet occurred at the dropout time (#25).
* Documented the minimum `interim_look` requirement (at least the block size for two-arm designs) in the `survival_adapt()` `interim_look` parameter.
* Improved the "Example: Two-armed RCT" vignette: the `summarise_sims()` operating characteristics are now rendered as captioned tables, a section documents one-sided tests (including that `method = "bayes-surv"` requires a one-sided alternative and measures the effect on the cumulative-failure-probability scale `p_treatment - p_control` against `h0`), and the `cutpoint` argument name was corrected to `cutpoints`.
* Added a new vignette, "Bayesian decisions with piecewise-exponential hazards", demonstrating `method = "bayes-surv"` with a piecewise hazard via `cutpoints` and `prop_to_haz()`, the Gamma-prior / posterior decision rule on the cumulative-failure-probability scale, and a worked single-trial example.
* Added a new vignette, "Single-arm trials", documenting the `hazard_control = NULL` mode (Bayesian-only), the role of `h0` as a benchmark failure rate, the success rule `Pr(p_treatment < h0) > prob_ha`, and a worked single-trial example with operating-characteristics templates.

## Package maintenance

* Expanded statistical checks for enrollment, randomization, event-time
  simulation and imputation, cumulative event probabilities, complete-data
  simulation, operating-characteristic summaries, completed-data analyses, and
  posterior calculations.
* Added a package architecture vignette with a function-dependency diagram.
* Added automated building and publication of the documentation website and
  updated the associated GitHub Actions.
* Improved accessibility text for the documentation website and README images.
* Removed obsolete AppVeyor configuration and excluded generated website files
  from the source package.
* Clarified that the C++ log-rank calculation was adapted from the deprecated
  `fastlogranktest` package.
* Clarified that `prior` uses a Gamma shape-rate parameterization and is shared
  across piecewise intervals and treatment arms.

# goldilocks 0.4.0

## Main updates

* The log-rank calculation previously obtained from `fastlogranktest` is now
  included in `goldilocks`, because `fastlogranktest` is no longer available on
  CRAN.

## Package maintenance

* Updated automated package checks and README status information.
* Corrected a condition used when summarizing simulated trials.

# goldilocks 0.3.1

## Features

* Added a chi-square test for binary outcomes.

## Bugs

* Corrected an error in `summarise_sims()`.

## Package maintenance

* Assigned complete-trial simulation, interim stopping decisions, and final
  analyses to `sim_comp_data()`, `test_stop_success()`, and `test_final()`,
  respectively.
* Updated package links and documentation and removed the obsolete AppVeyor
  check.

# goldilocks 0.3.0

* Added a `NEWS.md` file to track changes to the package.
* First CRAN submission.
