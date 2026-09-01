# Changelog

## goldilocks 0.6.0.9000

### Improvements

- Repeated Bayesian-survival completed-data analyses now use a guarded
  internal Gamma-posterior kernel when both the sufficient statistics
  and normalized priors were created by `goldilocks`. General entry
  points retain their full validation and flexible input handling, while
  exact seeded results, random number consumption, empty-interval
  behavior, and public outputs are unchanged
  ([\#89](https://github.com/graemeleehickey/goldilocks/issues/89)).
- Documentation now states explicitly that Bayesian binary prediction
  uses a piecewise-exponential Gamma model to impute pending outcomes
  and a separate beta-binomial model to analyze completed endpoint
  statuses. These stages are modular rather than one joint Bayesian
  model. The documented retained metadata identifies both prior
  specifications and the imputation horizon
  ([\#64](https://github.com/graemeleehickey/goldilocks/issues/64)).
- `prior_surv` and `prior_surv_final` now accept independent
  arm-specific Gamma priors as lists named `control` and `treatment` (or
  `treatment` for a single-arm design). Each arm may use a shared
  shape-rate vector or its own interval-specific matrix. Existing
  vectors and matrices still broadcast across arms and preserve seeded
  results. Resolved prior parameters are retained in metadata, while
  individual traced trials and externally evaluated interim looks report
  conjugate posterior diagnostics by arm and interval. Priors remain
  independent: no hierarchical or implicit borrowing is performed
  ([\#88](https://github.com/graemeleehickey/goldilocks/issues/88)).
- New
  [`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md)
  applies the same posterior-predictive decision calculation used by
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  to an externally observed interim data cut. Explicit follow-up
  statuses distinguish events, completed follow-up, pending outcomes,
  and early censoring. Potential future accrual is derived from
  `N_total`, observed arm counts, and `rand_ratio`, without requiring
  randomization blocks or concealed future assignments. The returned
  `goldilocks_interim` object includes predictive probabilities,
  diagnostic Monte Carlo uncertainty, imputation and allocation
  diagnostics, a compatible one-look decision trace, and auditable
  design, data-cut, package-version, and RNG metadata
  ([\#81](https://github.com/graemeleehickey/goldilocks/issues/81)).
- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  and
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now allow the event-time generator and predictive analysis to use
  different piecewise-exponential partitions. `generation_cutpoints`
  controls `hazard_treatment` and `hazard_control`, while `cutpoints`
  continues to control posterior estimation, predictive imputation,
  final analysis, and interval-specific priors. The new argument
  defaults to `cutpoints`, preserving existing calls and seeded results.
  The generation-only
  [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
  interface now names its partition `generation_cutpoints` explicitly;
  direct named calls to `sim_comp_data(cutpoints =)` must use
  `sim_comp_data(generation_cutpoints =)` instead
  ([\#95](https://github.com/graemeleehickey/goldilocks/issues/95)).
- `rand_ratio` and `randomization(allocation =)` now accept vectors
  named `control` and `treatment` in either order and normalize them
  internally. Legacy unnamed vectors remain interpreted as
  `c(control, treatment)`; unequal unnamed allocations warn that names
  may be required in a future major release. Allocation and arm-specific
  attrition now share one generic arm-vector normalizer, and retained
  argument metadata stores allocation in canonical arm order.
- `prop_loss` now accepts a named `control`/`treatment` vector for
  arm-specific attrition in two-arm trials, while a scalar applies the
  same proportion to every arm. Loss is selected separately within each
  randomized arm, retained argument metadata stores the normalized named
  proportions, and single-arm designs continue to require one value
  ([\#69](https://github.com/graemeleehickey/goldilocks/issues/69)).
- Piecewise event and enrollment times now use one explicit survival
  counting-process convention when assigned to intervals:
  `(start, stop]`, so a realized time exactly at a cutpoint belongs to
  the interval ending there. Posterior sufficient statistics and
  enrollment intensity inversion share the same boundary helper and
  regression tests cover times immediately before, at, and after
  multiple cutpoints. PWEALL remains the continuous event-time generator
  because its `[start, stop)` hazard representation differs only at
  zero-probability singleton boundaries; documentation and vignettes now
  make that distinction explicit
  ([\#70](https://github.com/graemeleehickey/goldilocks/issues/70)).
- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  and
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  results now retain an `arguments` attribute containing all evaluated
  argument values, including defaults. The plain named list can be
  serialized with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
  passed back to the corresponding function with
  [`do.call()`](https://rdrr.io/r/base/do.call.html) without relying on
  symbols retained in the original call
  ([\#86](https://github.com/graemeleehickey/goldilocks/issues/86)).
- [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  now reports explicit requested, analyzed, failed, and used simulation
  counts together with 95% Monte Carlo intervals for every probability
  and mean sample size. Wilson intervals remain informative for zero or
  one observed outcome, and an optional named `max_mcse` target warns
  when a scenario has insufficient Monte Carlo precision. Documentation
  and column names distinguish finite-simulation uncertainty from
  clinical or model uncertainty
  ([\#68](https://github.com/graemeleehickey/goldilocks/issues/68)).
- New
  [`summarise_calendar_time()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_calendar_time.md)
  reports wide operating-characteristic tables for trial duration,
  accrual, analysis readiness, person-time, peak concurrent follow-up,
  and interim-look timing. Calendar time starts at first enrollment;
  requested, analyzed, and failed denominators remain explicit, and no
  new simulation arguments are required. Per-trial results and retained
  traces now carry the compact timing values needed for these summaries
  ([\#87](https://github.com/graemeleehickey/goldilocks/issues/87)).

### Bug fixes

- Package-owned PSOCK workers now inherit the parent package library and
  load the `goldilocks` namespace before trial functions are dispatched.
  This keeps registered compiled routines available for multi-core
  log-rank simulations on Windows, including during source-package
  vignette builds.
- [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  now calculates the probability of reaching the maximum sample size
  from trial-level stopping indicators. Previously, the newly summarized
  futility probability shadowed the input indicator inside
  [`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html),
  producing incorrect `stop_max_N` values
  ([\#68](https://github.com/graemeleehickey/goldilocks/issues/68)).
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now caps package-owned workers at useful trial tasks and uses a
  documented automatic crossover of at least two trials per worker,
  avoiding parallel startup for workloads smaller than four trials.
  Package-owned PSOCK clusters are cleaned up on errors, returned
  parallel metadata records the deterministic execution plan, and
  benchmarks report setup separately from estimated compute and dispatch
  time
  ([\#93](https://github.com/graemeleehickey/goldilocks/issues/93)).
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now isolates trial-level errors across sequential and parallel
  execution. Successful trials are retained, failures are excluded into
  a compact `failures` table, and one warning reports the affected
  count; an all-failed batch still stops with the failure table attached
  to its error. Complete-result summaries retain the requested count and
  report the analyzed and failed counts without adding a new simulation
  argument
  ([\#72](https://github.com/graemeleehickey/goldilocks/issues/72)).
- Cox analyses now isolate the low-overhead
  [`survival::coxph.fit()`](https://rdrr.io/pkg/survival/man/agreg.fit.html)
  path behind a cached version/signature compatibility check and
  validate its returned coefficient and variance structure. Incompatible
  installations fall back automatically to exported
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html),
  without a new public argument. Singular fits, zero-event data, and
  convergence warnings produce targeted non-estimability errors; tests
  exercise both engines and benchmarks document their material
  performance tradeoff on realistic tied data
  ([\#77](https://github.com/graemeleehickey/goldilocks/issues/77)).
- [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  now accepts complete
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  results directly as well as the existing data-frame form. Named and
  grouped scenario identifiers are retained; summaries of complete
  results preserve backend, seed, failure, RNG, parallel, call, and
  design metadata and expose requested, analyzed, and failed trial
  counts so failed runs cannot silently leave the denominator unchanged
  ([\#78](https://github.com/graemeleehickey/goldilocks/issues/78)).
- Interim Monte Carlo decisions use strict point-estimate comparisons at
  both the completed-data posterior layer and the outer
  predictive-imputation layer, preserving the Goldilocks stopping rule.
  The defaults are now `N_impute = 500` and `N_mcmc = 1000`, and
  decision traces report estimates, diagnostic Monte Carlo standard
  errors and exact bounds, draw counts, uncertainty counts, and the
  decision reason
  ([\#60](https://github.com/graemeleehickey/goldilocks/issues/60)).
- Empty piecewise-exponential intervals are now prior-driven by default.
  The former `empty_interval = "propagate"` behavior remains available
  only as an explicit legacy heuristic. To reproduce earlier results,
  set it explicitly; otherwise review the interval-specific prior
  because it now supplies all information for an empty interval.
  Decision traces identify every interval handled by either policy
  ([\#62](https://github.com/graemeleehickey/goldilocks/issues/62)).
- Unseeded PSOCK simulations now derive independent per-trial streams
  from one draw from the caller’s RNG state. Resetting the caller state
  reproduces the call, consecutive calls advance rather than reuse
  streams, and returned RNG metadata records the effective backend and
  seed policy
  ([\#71](https://github.com/graemeleehickey/goldilocks/issues/71)).
- Log-rank analyses now reject nonzero `h0` values instead of silently
  testing the equal-survival null. Cox margins remain on the
  log-hazard-ratio scale, while Bayesian and risk-difference margins
  remain on their documented event- probability scales
  ([\#73](https://github.com/graemeleehickey/goldilocks/issues/73)).
- Interim futility and expected-success thresholds now share a
  documented scalar-or-exact-length contract, report expected and
  observed lengths in errors, and are retained in normalized
  decision-design metadata
  ([\#74](https://github.com/graemeleehickey/goldilocks/issues/74)).
- Exported piecewise-exponential utilities now validate draw counts and
  observed times consistently, with explicit zero-draw and zero-length
  behavior. Related binary-outcome validation now rejects missing and
  malformed event indicators before analysis
  ([\#75](https://github.com/graemeleehickey/goldilocks/issues/75)).
- Cumulative-hazard/event-probability conversions now use stable
  [`expm1()`](https://rdrr.io/r/base/Log.html) and
  [`log1p()`](https://rdrr.io/r/base/Log.html) transformations,
  retaining precision near zero and defining exact behavior at zero and
  infinite boundaries
  ([\#76](https://github.com/graemeleehickey/goldilocks/issues/76)).
- Randomization errors now state correctly that every block size must be
  a multiple of the sum of the two-arm allocation weights and include
  the observed values
  ([\#79](https://github.com/graemeleehickey/goldilocks/issues/79)).

## goldilocks 0.6.0

CRAN release: 2026-07-29

### Improvements

- New
  [`plot_enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/plot_enrollment.md)
  draws the expected enrollment projection, simulated trial
  trajectories, and interim and maximum-sample-size milestones for
  constant or piecewise enrollment rates. It accepts explicit design
  arguments or extracts the evaluated design retained on
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  and
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  results, with optional follow-up and time-unit annotations.
- Survival-prior arguments are now named `prior_surv` and
  `prior_surv_final`, while the Bayesian binary prior is named
  `prior_bin`. `prior_surv_final` defaults to `prior_surv`, but can
  specify a different Gamma prior for final imputation and final
  piecewise-exponential analysis. Survival priors may be length-two
  shape-rate vectors shared across intervals or two-row matrices with
  one shape-rate column per piecewise interval. The ADVENT vignette now
  expresses its interval- and stage-specific imputation priors directly
  ([\#59](https://github.com/graemeleehickey/goldilocks/issues/59)).
- **Major behavior change from goldilocks 0.5.0 and earlier:**
  [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  now simulates exact continuous-time arrivals from a piecewise-constant
  Poisson process by inverting its cumulative intensity. Fractional
  enrollment-rate changes are handled exactly, and runtime depends on
  the requested sample size rather than the number of empty unit-time
  bins. Enrollment schedules now follow the internal-knot convention
  used by hazard schedules: `lambda_time = NULL` represents a constant
  rate, each supplied positive internal knot requires one additional
  `lambda` value, and time zero remains the implicit first-patient-in
  origin.
  [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
  now consumes these continuous enrollment times directly without
  post-hoc uniform jitter
  ([\#55](https://github.com/graemeleehickey/goldilocks/issues/55)).
  Previously,
  [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  generated Poisson counts in unit-time bins, returned rebased integer
  bin times, and
  [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
  added uniform jitter. Existing calls must replace `lambda_time = 0`
  with `NULL` and remove the leading zero from piecewise schedules.
  Seeded simulations will not reproduce results from version 0.5.0 or
  earlier, and operating-characteristic estimates may change, especially
  for fractional knots or low enrollment rates.
- The frequentist binary `method = "chisq"` analysis has been replaced
  by `method = "riskdiff"`. It estimates the treatment-minus-control
  event-risk difference, supports one- and two-sided Wald tests against
  `h0`, and reports that difference in `est_final`. When
  `imputed_final = TRUE`, estimates and within-imputation variances are
  combined using Rubin’s rules. The historical mean-of-`1 - P` rule is
  no longer used. Imputed log-rank final analyses remain unsupported
  ([\#56](https://github.com/graemeleehickey/goldilocks/issues/56)).
- Cox regression now supports `imputed_final = TRUE`. The final log
  hazard ratios and within-imputation variances are combined using
  Rubin’s rules, and `post_prob_ha` reports `1 - P` from the pooled Wald
  test. Non-imputed Cox analyses are unchanged
  ([\#56](https://github.com/graemeleehickey/goldilocks/issues/56)).
- Bayesian survival and Monte Carlo beta-binomial analyses now use the
  shared `N_mcmc` argument for posterior sampling; the separate `bin_N`
  argument has been removed. Their internal analysis results return the
  posterior mean effect instead of all posterior draws.
- survival_adapt() can now return an optional, tidy interim decision
  trace with predictive probabilities, thresholds, stopping decisions,
  arm-level counts, and relevant warnings. New helpers summarize and
  plot individual traces and simulation stopping outcomes
  ([\#57](https://github.com/graemeleehickey/goldilocks/issues/57)).
- [`plot_sim_stopping()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_stopping.md)
  now offers marginal, conditional, cumulative, and flowchart views of
  stopping outcomes by enrolled sample size. Plot subtitles state
  bar-chart denominators explicitly, while the flowchart shows trial
  counts branching through futility, continued enrollment, and early
  success at successive looks. Bar-chart legends are placed beyond the
  plotting region with margins sized to their labels, avoiding overlap
  across device themes and text sizes. Subtitles are left-aligned so
  long denominator descriptions stay visible, and the cumulative view
  uses a compact y-axis label that remains clear of its title.
  Percentage labels use a compact, consistent size across bar-chart
  views to avoid collisions between adjacent looks. When retained
  simulation traces are supplied, conditional, cumulative, and flowchart
  views include reached looks at which no trial stopped.
- New
  [`plot_sim_ocs()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_ocs.md)
  plots success and stopping probabilities together with expected sample
  size across treatment-effect simulation scenarios.
- `sim_trials(return_trace = TRUE)` now retains compact interim traces
  across simulations, and new
  [`plot_sim_decisions()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_decisions.md)
  visualizes their predictive-probability decision regions by interim
  look.
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now supports reproducible PSOCK parallel execution on Windows and an
  explicit `backend` argument. The default Unix fork path is retained;
  `backend = "auto"` selects the appropriate implementation for the
  platform. Seeded simulations now also restore the caller’s RNG state
  ([\#40](https://github.com/graemeleehickey/goldilocks/issues/40)).
- Consolidated the Bayesian binomial analysis implementation so its
  posterior calculations and method-specific branches are co-located in
  `bayes_binomial_test()`.
- Renamed the piecewise-exponential Bayesian survival analysis method
  from `method = "bayes"` to `method = "bayes-surv"` to distinguish it
  from `method = "bayes-bin"`.
- Added `method = "bayes-bin"` for Bayesian beta-binomial analysis of
  complete binary outcomes, with Monte Carlo, normal approximation, and
  quadrature options for treatment-control differences.
- Binary endpoint analyses can now set `binary_imputation = "bernoulli"`
  to draw event status directly from the conditional
  piecewise-exponential event probability. The existing conditional
  event-time approach remains the default
  ([\#21](https://github.com/graemeleehickey/goldilocks/issues/21)).
- Cox model analyses now use a lower-level survival fit for repeated
  Wald tests, avoiding formula and summary overhead in simulation hot
  paths.
- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  now uses `h0` as the null log hazard ratio for Cox model tests,
  allowing non-inferiority testing with `h0 = log(margin)`.
- Added maintainer performance benchmarks for simulation hot paths,
  including posterior probability conversion, posterior sampling,
  imputation, and representative
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  runs
  ([\#42](https://github.com/graemeleehickey/goldilocks/issues/42)).
- Harmonized treatment-assignment terminology in internal simulation
  helpers and documentation: data use `treatment = 1` for treatment and
  `treatment = 0` for control, while posterior/imputation array indexing
  is now described separately as hazard slices
  ([\#35](https://github.com/graemeleehickey/goldilocks/issues/35)).
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now accepts a `seed` argument that creates independent per-trial
  `"L'Ecuyer-CMRG"` random-number streams for reproducible simulations,
  including when using multiple cores
  ([\#41](https://github.com/graemeleehickey/goldilocks/issues/41)).
- [`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md)
  now computes piecewise-exponential cumulative event probabilities
  directly from the cumulative hazard, avoiding row-wise calls to
  [`PWEALL::pwe()`](https://rdrr.io/pkg/PWEALL/man/pwe.html) in Bayesian
  posterior summaries
  ([\#34](https://github.com/graemeleehickey/goldilocks/issues/34)).
- `posterior()` now computes piecewise-exponential sufficient statistics
  directly, avoiding
  [`survSplit()`](https://rdrr.io/pkg/survival/man/survSplit.html) and
  grouped `dplyr` summarization in a simulation hot path.
- Bayesian survival analyses after predictive imputation now draw their
  fresh completed-data posterior directly from arm-by-interval exposure
  and event sufficient statistics. This preserves the documented
  two-stage posterior procedure while avoiding repeated patient-level
  posterior setup
  ([\#38](https://github.com/graemeleehickey/goldilocks/issues/38)).
- [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  and
  [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  no longer grow vectors repeatedly inside their simulation loops
  ([\#44](https://github.com/graemeleehickey/goldilocks/issues/44)).
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  and
  [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  now use
  [`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
  for faster binding of simulation result data frames
  ([\#43](https://github.com/graemeleehickey/goldilocks/issues/43)).

### Bug fixes

- Piecewise Bayesian analyses now retain the posterior-draw dimension
  when `N_mcmc = 1`, avoiding a matrix-validation error in probability
  conversion.
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now uses the same default alternative hypothesis and futility
  threshold as
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md),
  so omitted arguments define the same adaptive design
  ([\#47](https://github.com/graemeleehickey/goldilocks/issues/47)).
- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  now requires interim looks to be strictly increasing, preventing
  non-chronological or duplicated analyses
  ([\#48](https://github.com/graemeleehickey/goldilocks/issues/48)).
- Piecewise hazard inputs are now validated as finite and non-negative
  across simulation, imputation, and probability helpers. A finite
  `maxtime` is now required when the final hazard is zero
  ([\#49](https://github.com/graemeleehickey/goldilocks/issues/49)).
- `h0` must now be a single finite value, with probability-scale bounds
  enforced for Bayesian analyses
  ([\#50](https://github.com/graemeleehickey/goldilocks/issues/50)).
- The `cutpoints` API now accepts only interior hazard-change times:
  `NULL` represents a constant hazard, and each supplied cutpoint
  requires one additional hazard rate. The implicit interval start at
  zero is handled internally
  ([\#58](https://github.com/graemeleehickey/goldilocks/issues/58)).
  Piecewise model inputs require finite, positive, strictly increasing
  cutpoints and a finite study endpoint after the final cutpoint; this
  validation is shared by simulation and probability helpers
  ([\#52](https://github.com/graemeleehickey/goldilocks/issues/52)).
- [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  now validates its complete schedule before generating data, including
  integer sample size, finite positive rates, and finite strictly
  increasing knots
  ([\#53](https://github.com/graemeleehickey/goldilocks/issues/53)).
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now defaults to serial execution with `ncores = 1L`, matching its
  documentation and avoiding unexpected use of available cores
  ([\#54](https://github.com/graemeleehickey/goldilocks/issues/54)).
- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  now works with the documented default success and futility thresholds
  when `interim_look = NULL`; thresholds are ignored when there are no
  interim looks.
- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md),
  [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md),
  and
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now validate probability, prior, and positive-integer count arguments
  up front, avoiding invalid simulations or low-level downstream errors.
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now falls back to `ncores = 1L` when the available core count cannot
  be detected.
- Log-rank and Cox analyses now fail with clear non-estimable-analysis
  errors when the test statistic cannot be computed reliably, rather
  than allowing `NA` or `NaN` values to break adaptive decision logic.
- [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  now requires a two-element positive integer allocation vector,
  preventing invalid two-arm allocation schedules such as `c(0, 1)`.
- [`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md)
  now rejects invalid cumulative event probabilities, including
  non-finite values, values outside `[0, 1)`, and values that decrease
  over time.
- [`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md)
  now errors if `maxtime` is earlier than the observed conditioning
  time, preventing imputed records from moving backward in time.
- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  no longer dereferences the disabled futility result when `Fn = 0`, so
  futility counters remain inert when futility monitoring is turned off
  ([\#33](https://github.com/graemeleehickey/goldilocks/issues/33)).
- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  and
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  now expose an `empty_interval` policy for Bayesian
  piecewise-exponential posterior calculations, allowing empty intervals
  to use legacy propagation, prior-only updating, or strict errors.

### Documentation

- Updated the package title, description, help landing page, README, and
  architecture overview to reflect the current support for both survival
  and fixed-time binary endpoints, all five final-analysis methods, and
  the complete plotting API.
- Corrected the technical-methods vignette’s section numbering and
  updated its overview to cover the package’s binary endpoint methods.
- Expanded the simulation, design, and technical vignettes with examples
  of
  [`plot_enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/plot_enrollment.md),
  [`plot_sim_ocs()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_ocs.md),
  [`plot_sim_stopping()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_stopping.md),
  and
  [`plot_sim_decisions()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_decisions.md).
- Modernized roxygen2 source comments to use markdown tables, links,
  code spans, and emphasis in place of older Rd markup where
  appropriate.
- Added a new vignette, “Technical details of the Goldilocks design”,
  documenting the design notation, piecewise-exponential event-time
  model, Gamma posterior updating, posterior predictive probabilities,
  interim decision rules, final analysis options, and simulation-based
  calibration.
- Added a new vignette, “Bayesian binary outcome designs”, documenting
  `method = "bayes-bin"` for two-arm and single-arm complete binary
  endpoint analyses.
- Added a new vignette, “ADVENT: a published Goldilocks design”, showing
  how the published ADVENT pulsed field ablation trial maps to
  `method = "bayes-bin"` with beta-binomial non-inferiority endpoints,
  adaptive sample-size thresholds, and cached simulation examples.
- Renamed vignettes for consistency: “Two-arm randomized trials”,
  “Bayesian piecewise-exponential designs”, “Single-arm designs with a
  performance goal”, and “Package architecture”.
- Clarified that `goldilocks` treats enrollment time and randomization
  time as the same time point in its time-to-event simulations.
- Clarified that the single-arm external benchmark `h0` is often
  referred to as a performance goal (PG) or objective performance
  criterion (OPC).
- Expanded technical documentation for the Goldilocks
  predictive-probability algorithm, including notation for final
  analysis quantities, operating characteristics, and method-specific
  decision rules.
- Clarified the `Fn` documentation to state that `Fn = 0` disables
  futility monitoring.
- Added a pkgdown light switch so the documentation site supports light,
  dark, and automatic color modes.
- Added CRAN checks, CRAN monthly downloads, and GPL-3 license badges to
  the README, and reordered the badge block.

## goldilocks 0.5.0

CRAN release: 2026-06-10

### Improvements

- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  now supports one-sided tests (`alternative = "greater"` or `"less"`)
  for `method = "cox"` and `method = "logrank"`. The chi-square test
  remains two-sided only
  ([\#20](https://github.com/graemeleehickey/goldilocks/issues/20)).
- When `method = "chisq"` and `imputed_final = FALSE`, subjects lost to
  follow-up are now excluded from the final analysis. Previously, LTFU
  subjects were counted as non-events, which diluted the event rate and
  biased the chi-square test
  ([\#22](https://github.com/graemeleehickey/goldilocks/issues/22)).
- Removed unused C++ global variables and dead threading code inherited
  from the deprecated `fastlogranktest` package.
- Replaced all uses of the magrittr pipe (`%>%`) with the base pipe
  (`|>`) in
  [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  and `posterior()`, and removed the `dplyr::%>%` re-export from the
  NAMESPACE.
- `posterior()` now warns when a piecewise interval has zero subjects
  and data is propagated from an adjacent interval.
- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  now validates that each `interim_look` in a two-arm design is at least
  the block size, so an interim look that could enrol a single treatment
  arm only is rejected as an input error rather than producing an
  undefined interim posterior.

### Bug fixes

- `posterior()` now propagates data for zero-exposure piecewise
  intervals *within* each treatment arm. Previously, propagation walked
  the flat row order of the per-arm summary, so an empty leading
  interval in the treatment arm (e.g. a look where the treatment arm has
  no subjects) could copy the control arm’s data into the treatment
  posterior, contaminating the estimate across arms. `posterior()` also
  now checks up front that the expected treatment arms are present in
  the supplied data, erroring with an informative message instead of
  silently returning an all-`NA` posterior slice for an absent arm
  (which could occur at a small interim look where one arm has no
  enrolled subjects yet).

- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  no longer errors when called without interim looks
  (`interim_look = NULL`). The final analysis previously relied on an
  undefined loop index variable, which has been replaced with
  `stage_trial_stopped`.

- `impute_data()` no longer uses hard-coded positional column subsetting
  (`[, 1:10]`). Temporary columns are now dropped by name, making the
  function robust to upstream changes in the data frame structure
  ([\#26](https://github.com/graemeleehickey/goldilocks/issues/26)).

- [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  no longer produces `NA` for `next_block` when the loop exhausts all
  elements of a multi-element `block` vector. The index now wraps around
  cyclically
  ([\#31](https://github.com/graemeleehickey/goldilocks/issues/31)).

- `analyse_data()` now uses explicit row/column indexing when extracting
  Cox model results, preventing silent errors if the summary matrix
  structure changes
  ([\#29](https://github.com/graemeleehickey/goldilocks/issues/29)).

- `analyse_data()` now imports
  [`stats::pnorm`](https://rdrr.io/r/stats/Normal.html), which is used
  to compute one-sided `p`-values for the Cox and log-rank tests.
  Previously this relied on `stats` being attached.

- Removed a dead `loss_to_fu <- NA` assignment in
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  that was shadowed by the `loss_to_fu` column inside
  [`within()`](https://rdrr.io/r/base/with.html) and never used.

- [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  now correctly selects the enrollment rate at piecewise changepoints.
  Previously, the rate at exact changepoint boundaries could use the
  rate from the prior interval
  ([\#28](https://github.com/graemeleehickey/goldilocks/issues/28)).

- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  no longer adds a systematic perturbation (`sd(time) / 1e4`) to all
  survival times at interim looks. Instead, only the boundary subject
  with zero follow-up time is clamped to `.Machine$double.eps` to
  satisfy
  [`survSplit()`](https://rdrr.io/pkg/survival/man/survSplit.html)
  requirements
  ([\#24](https://github.com/graemeleehickey/goldilocks/issues/24)).

### Documentation

- Removed the `est_interim` element from the
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  return-value documentation. This field was documented but never
  computed or returned.
- Documented the two-stage posterior procedure used when
  `method = "bayes-surv"` with imputation, clarifying that the
  imputation model’s posterior influences the analysis posterior
  ([\#27](https://github.com/graemeleehickey/goldilocks/issues/27)).
- Clarified `prop_loss` parameter documentation, explaining that LTFU
  times are drawn from `Uniform(0, t)` and that the event has not yet
  occurred at the dropout time
  ([\#25](https://github.com/graemeleehickey/goldilocks/issues/25)).
- Documented the minimum `interim_look` requirement (at least the block
  size for two-arm designs) in the
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  `interim_look` parameter.
- Improved the “Example: Two-armed RCT” vignette: the
  [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  operating characteristics are now rendered as captioned tables, a
  section documents one-sided tests (including that
  `method = "bayes-surv"` requires a one-sided alternative and measures
  the effect on the cumulative-failure-probability scale
  `p_treatment - p_control` against `h0`), and the `cutpoint` argument
  name was corrected to `cutpoints`.
- Added a new vignette, “Bayesian decisions with piecewise-exponential
  hazards”, demonstrating `method = "bayes-surv"` with a piecewise
  hazard via `cutpoints` and
  [`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md),
  the Gamma-prior / posterior decision rule on the
  cumulative-failure-probability scale, and a worked single-trial
  example.
- Added a new vignette, “Single-arm trials”, documenting the
  `hazard_control = NULL` mode (Bayesian-only), the role of `h0` as a
  benchmark failure rate, the success rule
  `Pr(p_treatment < h0) > prob_ha`, and a worked single-trial example
  with operating-characteristics templates.

### Housekeeping

- Added unit tests for
  [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md),
  [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md),
  [`pwe_sim()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_sim.md),
  [`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md),
  [`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md),
  [`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md),
  [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md),
  [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md),
  `analyse_data()`, and `posterior()`.
- Added a package architecture vignette with a function dependency
  diagram.
- Added `_pkgdown.yml` configuration for a documentation website.
- Added GitHub Actions workflow for pkgdown site deployment.
- Updated GitHub Actions (`actions/checkout`, `actions/upload-artifact`)
  from v4 to v5 for Node.js 24 compatibility.
- Added `.positai`, `_pkgdown.yml`, and `docs` to `.Rbuildignore` to
  suppress `R CMD check` NOTEs.
- Removed the unused `appveyor.yml` CI configuration file and its stale
  `.Rbuildignore` entry.
- Added `docs/` to `.gitignore`, since the pkgdown site is built and
  deployed to the `gh-pages` branch by GitHub Actions rather than
  committed to the main branch.
- Added an `aria-label` to the pkgdown navbar GitHub icon and alt-text
  to the README hex logo and downloads badge to address pkgdown
  accessibility warnings.
- Updated README to clarify that the C++ log-rank code was ported from
  the now-deprecated `fastlogranktest` package.
- Clarified `prior` parameter documentation to explicitly state the
  Gamma rate parameterization and that the same prior is shared across
  all piecewise intervals and treatment arms.

## goldilocks 0.4.0

CRAN release: 2025-01-08

### Main updates

- Because `fastlogranktest` is no longer available on CRAN, a copy of
  the C++ code and wrapper from the CRAN archive have been included
  directly into the source code of this package.

### Housekeeping

- Updated GitHub actions workflows
- Updated README badges
- Changed logic check in `summarise_sims.R`

## goldilocks 0.3.1

### Features

- Added a chi-square test for binomial data

### Bugs

- Fixed a small bug in
  [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)

### Housekeeping

- Restructured code so complete trial data is generated by
  [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
- Restructured code so interim stopping tests are done by
  `test_stop_success()`
- Restructured code so final analysis is done by `test_final()`
- Added URL links to DESCRIPTION
- Updates to documentation
- Removed appveyor check and badge

## goldilocks 0.3.0

CRAN release: 2021-05-10

- Added a `NEWS.md` file to track changes to the package.
- First CRAN submission.
