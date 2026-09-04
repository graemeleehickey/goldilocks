#' @title Simulate and analyze one Goldilocks adaptive trial
#'
#' @description Simulates one single-arm or randomized two-arm trial under a
#'   Goldilocks sample-size design. At each planned interim look, posterior
#'   predictive probabilities determine whether to stop accrual for expected
#'   success, stop for futility, or continue toward the maximum sample size.
#'
#' @inheritParams sim_comp_data
#' @param cutpoints `NULL` (the default), or a numeric vector of finite,
#'   positive, strictly increasing interior follow-up
#'   times defining the piecewise-exponential model used for interim posterior
#'   estimation, predictive imputation, and final analysis. The number of
#'   interval-specific prior columns must be one greater than the number of
#'   cutpoints. `NULL` specifies a constant-hazard analysis model.
#' @param generation_cutpoints `NULL`, or a numeric vector of finite, positive,
#'   strictly increasing interior
#'   follow-up times defining the piecewise-exponential model used to generate
#'   event times. `hazard_treatment` and `hazard_control` must each have one
#'   value per resulting interval. Defaults to `cutpoints`, preserving the
#'   historical behavior in which generation and analysis used one partition.
#' @param end_of_study A required finite, positive numeric value giving the
#'   planned subject-level follow-up time. It must be greater than the final
#'   value in both `cutpoints` and `generation_cutpoints`, when supplied, and
#'   use the same time unit.
#' @param interim_look `NULL` (the default) for no interim analyses, or a
#'   strictly increasing positive integer vector giving the cumulative sample
#'   size at each interim look. Do not include the maximum sample size. For
#'   two-arm designs, each
#'   interim look must be at least the (largest) block size (see `block`),
#'   ensuring both treatment groups are present at every interim analysis; a
#'   smaller look could enroll subjects from one treatment group only, leaving
#'   the interim posterior undefined for the missing group.
#' @param prior_surv A numeric vector, matrix, or named list specifying the
#'   Gamma prior for the
#'   piecewise-exponential hazards used during interim prediction. A length-two
#'   vector supplies shape and rate and applies the same prior to every arm and
#'   interval. A `2` by `length(cutpoints) + 1` matrix supplies
#'   interval-specific values shared by all arms, with shapes in row 1 and rates
#'   in row 2. For independent arm-specific priors, supply a list named
#'   `control` and `treatment` in a two-arm design, or `treatment` in a
#'   single-arm design. Each list element may be a length-two vector or an
#'   interval-specific matrix. Both arms must be supplied; no values are
#'   borrowed or filled from the other arm. Rates must use the same time unit as
#'   event times, exposure, and cutpoints. The default is `c(0.1, 0.1)`.
#' @param prior_surv_final A numeric vector, matrix, or named list specifying
#'   the Gamma prior
#'   used for final-stage piecewise-exponential imputation and, for `method =
#'   "bayes-surv"`, final analysis. It accepts the same shared or arm-specific
#'   forms as `prior_surv` and defaults to `prior_surv`, preserving the
#'   historical behavior.
#' @param prior_bin A length-two numeric vector of finite, positive shape
#'   parameters `c(a, b)` for the `Beta(a, b)` event-probability prior used when
#'   `method = "bayes-bin"`. The same prior is applied to both arms. The default
#'   is `c(1, 1)`, a uniform prior.
#' @param bin_method A single character string selecting how to calculate the
#'   posterior probability for `method = "bayes-bin"`. It must be one of
#'   `"mc"` (Monte Carlo
#'   sampling), `"normal"` (normal approximation), or `"quadrature"` (numerical
#'   integration). The default is `"mc"`.
#' @param binary_imputation A single character string selecting the predictive
#'   imputation approach for
#'   `method = "bayes-bin"`, `method = "riskdiff-wald"`, or `method =
#'   "riskdiff-fm"`. `"event-time"` (the default) draws a conditional
#'   piecewise-exponential event time and reduces it to event status at
#'   `end_of_study`. `"bernoulli"` draws the endpoint status directly from its
#'   conditional event probability. This argument is ignored for time-to-event
#'   analysis methods.
#' @param alternative A single character string specifying the alternative
#'   hypothesis. It must be one of `"greater"` (the default), `"less"`, or
#'   `"two.sided"`. One-sided alternatives (`"greater"` and `"less"`) are
#'   supported for `method = "bayes-surv"` and `method = "bayes-bin"`. All three
#'   options are supported for `method = "logrank"`, `method = "cox"`,
#'   `method = "riskdiff-wald"`, and `method = "riskdiff-fm"`. For
#'   survival outcomes, `"less"` corresponds to the treatment arm having a lower
#'   cumulative incidence (i.e., treatment is beneficial), and `"greater"`
#'   corresponds to the treatment arm having a higher cumulative incidence.
#' @param h0 A single finite numeric value specifying the null hypothesis or
#'   margin. The default is `0`. For Bayesian analyses, `h0` must lie in
#'   `[0, 1]` for a
#'   single-arm design and `[-1, 1]` for a two-arm design.
#'   * When `method = "bayes-surv"`, `h0` is the null value of
#'     \eqn{p_\textrm{treatment} - p_\textrm{control}}. In a single-arm design,
#'     `h0` is the external benchmark event probability, often referred to
#'     as a performance goal (PG) or objective performance criterion (OPC).
#'   * When `method = "bayes-bin"`, `h0` is the null value of
#'     \eqn{p_\textrm{treatment} - p_\textrm{control}} for a two-arm design, or
#'     the null event probability for a single-arm design.
#'   * When `method = "cox"`, `h0` is the null log hazard ratio for
#'     treatment versus control. Use `h0 = 0` for the usual hazard ratio
#'     of 1 null, or `h0 = log(margin)` for a non-inferiority margin
#'     specified as a hazard ratio. A Cox non-inferiority test should usually
#'     use `alternative = "less"`.
#'   * When `method = "riskdiff-wald"` or `method = "riskdiff-fm"`, `h0` is the
#'     null value of
#'     \eqn{p_\textrm{treatment} - p_\textrm{control}} and must lie in
#'     `[-1, 1]`.
#'   * When `method = "logrank"`, only `h0 = 0` is supported; this denotes the
#'     usual equal-survival null. Nonzero values are rejected because the
#'     standard log-rank statistic does not implement a nonzero effect margin.
#' @param Fn `NULL`, or a numeric vector of probabilities in `[0, 1]`. Each
#'   value is the predictive-probability
#'   threshold to stop at the \eqn{i}-th look early for futility. If there are
#'   no interim looks (i.e. `interim_look = NULL`), then `Fn` is not
#'   used in the simulations or analysis. Set `Fn = 0` to disable futility
#'   monitoring; `Fn = NULL` has the same effect. Supply either one value,
#'   which is repeated at every interim look, or exactly one value per
#'   `interim_look`. Other lengths are rejected rather than recycled. The
#'   default is `0.05`.
#' @param Sn A numeric vector of probabilities in `[0, 1]`. Each value is the
#'   predictive-probability
#'   threshold to stop at the \eqn{i}-th look early for expected success. If
#'   there are no interim looks (i.e. `interim_look = NULL`), then
#'   `Sn` is not used in the simulations or analysis. Supply either one value,
#'   which is repeated at every interim look, or exactly one value per
#'   `interim_look`. Other lengths are rejected rather than recycled. The
#'   default is `0.9`.
#' @param prob_ha A single numeric probability in `[0, 1]` defining success in
#'   each completed-data analysis. For Bayesian methods this is compared with
#'   the posterior probability of the alternative; for frequentist methods it
#'   is compared with `1 - P`. The default is `0.95`.
#' @param N_impute A positive integer giving the number of predictive
#'   imputations used at each
#'   interim look and, when requested, for final multiple imputation. The
#'   default is `500`. An imputed Cox or risk-difference final analysis requires
#'   at least two.
#' @param N_mcmc A positive integer giving the number of posterior draws used
#'   within each
#'   `method = "bayes-surv"` and by `method = "bayes-bin"` when
#'   `bin_method = "mc"`. The default is `1000`.
#' @param mc_conf_level A single numeric probability strictly between `0.5` and
#'   `1`, giving the confidence
#'   level for one-sided exact binomial bounds reported as diagnostics of
#'   finite Monte Carlo uncertainty. The bounds do not alter completed-data
#'   success classifications or interim decisions, which use strict point-
#'   estimate comparisons with `prob_ha`, `Sn`, and `Fn`. The default is `0.95`.
#' @param empty_interval A single character string specifying how to handle
#'   empty piecewise-exponential intervals when updating Gamma hazard models for
#'   predictive imputation and Bayesian survival analysis. An empty
#'   interval is an interval with no exposed subjects in a treatment arm at the
#'   analysis time. `"prior"` (the default) leaves the interval at zero exposure
#'   time and zero events, so its posterior is driven only by its assigned
#'   survival prior. `"propagate"` is a legacy heuristic that copies exposure
#'   time and event counts from the nearest non-empty interval in the same
#'   treatment arm and emits a warning. `"error"` stops when any empty interval
#'   is found.
#' @param method A single character string specifying the completed-data and
#'   final analysis. Available choices are a log-rank
#'   (`method = "logrank"`) test, Cox proportional hazards regression model
#'   Wald test (`method = "cox"`), a fully-Bayesian piecewise-exponential
#'   analysis (`method = "bayes-surv"`), a Bayesian beta-binomial analysis of
#'   complete binary outcomes (`method = "bayes-bin"`), a frequentist
#'   risk-difference Wald test (`method = "riskdiff-wald"`), or a
#'   Farrington-Manning score test (`method = "riskdiff-fm"`) of complete
#'   binary outcomes. The deprecated `method = "riskdiff"` is accepted as an
#'   alias for `"riskdiff-wald"` with a warning. The default is `"logrank"`.
#'   See Details.
#' @param imputed_final A single logical value indicating whether the final
#'   analysis should be based on imputed outcomes for
#'   subjects who were LTFU (i.e. right-censored with time less than
#'   `end_of_study`). The default is `FALSE`, which means that the final analysis
#'   incorporates right-censoring. With `method = "cox"`, `method =
#'   "riskdiff-wald"`, or `method = "riskdiff-fm"`, setting this to `TRUE`
#'   analyzes each imputed dataset and pools the scalar treatment effects and
#'   variances using Rubin's rules; this requires `N_impute >= 2`. The pooled
#'   risk-difference analysis is a Wald test rather than a Farrington-Manning
#'   test. Imputed final analyses remain unavailable for `method = "logrank"`.
#' @param return_trace A single logical value indicating whether the interim
#'   decision path should be returned in addition to the usual final summary.
#'   The default, `FALSE`, returns the historical one-row data frame. When
#'   `TRUE`, the result is a `goldilocks_trial` object with summary, trace, prior and posterior
#'   diagnostics, and call elements.
#'
#' @details Implements the Goldilocks design method described in Broglio et al.
#'   (2014). At each interim analysis, two probabilities are computed:
#'
#'   1. **The posterior predictive probability of eventual success.** This is
#'      calculated as the proportion of imputed datasets at the *current* sample
#'      size that satisfy the completed-data success criterion. At each
#'      interim analysis this proportion is compared to the corresponding
#'      element of `Sn`, and if it exceeds the threshold,
#'      accrual/enrollment is suspended and the outstanding follow-up allowed to
#'      complete before conducting the pre-specified final analysis.
#'
#'   2. **The posterior predictive probability of success at the maximum sample
#'      size.** This is
#'      calculated as the proportion of imputed datasets at the *maximum*
#'      sample size that satisfy the completed-data success criterion. It is
#'      compared to the corresponding element of `Fn`, and if it is below the
#'      threshold, accrual/enrollment is suspended and the
#'      trial terminated. Typically this would be a binding decision. If it is
#'      not a binding decision, then one should also explore the simulations
#'      with `Fn = 0`.
#'
#'   Hence, each interim look has three possible decisions:
#'
#'   1. **Stop for expected success**
#'   2. **Stop for futility**
#'   3. **Continue to enroll** new subjects, or if at maximum sample size,
#'      proceed to final analysis.
#'
#'   The following completed-data analysis methods are available at interim and
#'   final analyses:
#'
#'  * Log-rank test (`method = "logrank"`).
#'      Each (imputed) dataset with both treatment and control arms can be
#'      compared using a standard log-rank test. The output is a *P*-value,
#'      and there is no treatment effect reported. The function returns \eqn{1 -
#'      P}, which is reported in `post_prob_ha`. Whilst not a posterior
#'      probability, it can be contrasted in the same manner. For example, if
#'      the success threshold is \eqn{P < 0.05}, then one requires
#'      `post_prob_ha` \eqn{> 0.95}. The reason for this is to enable
#'      simple switching between Bayesian and frequentist paradigms for
#'      analysis. When `alternative = "less"` or `"greater"`, a
#'      one-sided *P*-value is computed from the log-rank z-statistic.
#'
#'   * Cox proportional hazards regression Wald test (`method = "cox"`).
#'      Similar to the log-rank test, a *P*-value is calculated and
#'      \eqn{1 - P} is reported in `post_prob_ha`. When
#'      `alternative = "two.sided"`, the standard two-sided Wald
#'      *P*-value is used when `h0 = 0`. For other values of
#'      `h0`, the Wald test is centered on the specified null log hazard
#'      ratio. When `alternative = "less"` or
#'      `"greater"`, a one-sided *P*-value is derived from the Wald
#'      z-statistic relative to `h0`. The treatment effect (log hazard
#'      ratio) is also reported. When `imputed_final = TRUE`, the Cox model is
#'      fitted separately to each of at least two imputed datasets. The log
#'      hazard ratios and their within-imputation variances are combined using
#'      Rubin's rules; the pooled Wald test uses Rubin's large-sample degrees of
#'      freedom. When `imputed_final = FALSE`, the existing single Cox model is
#'      fitted directly to the observed right-censored data.
#'
#'   * Bayesian difference in cumulative event probability
#'     (`method = "bayes-surv"`).
#'      Each imputed dataset is used to update the conjugate Gamma prior
#'      (defined by `prior_surv` at interim looks and `prior_surv_final` at the
#'      final stage), yielding a posterior distribution for the
#'      piecewise exponential rate parameters. In turn, the posterior
#'      distribution of the cumulative incidence function (\eqn{1 - S(t)}, where
#'      \eqn{S(t)} is the survival function) evaluated at time
#'      `end_of_study` is calculated. In a single-arm study, inference concerns
#'      the treatment-arm event probability. In a two-arm study, the independent
#'      arm-specific posteriors define the posterior distribution of the
#'      treatment-minus-control difference. The reported posterior probability
#'      is determined by `alternative` and `h0`.
#'
#'      For piecewise-exponential analyses, an interim or final dataset may
#'      contain intervals with no exposed subjects in one treatment arm,
#'      especially when later cutpoints occur after the available follow-up at
#'      early looks. The `empty_interval` argument controls this case.
#'      The default, `"prior"`, leaves an empty interval prior-driven, making
#'      the absence of interval data explicit. The legacy `"propagate"` option
#'      borrows sufficient statistics from the nearest non-empty interval
#'      within the same treatment arm. It is operationally stable but
#'      statistically consequential because adjacent observed data then inform
#'      the empty interval's posterior. `"error"` is strict and stops the
#'      simulation or analysis when an empty interval is encountered.
#'
#'   * Bayesian beta-binomial analysis (`method = "bayes-bin"`).
#'      Each complete or imputed dataset is reduced to binary event outcomes at
#'      `end_of_study`. A conjugate `Beta(a, b)` prior, specified with
#'      `prior_bin`, is updated with the number of events and non-events in
#'      each arm. In a single-arm study, inference is based on the posterior
#'      event probability. In a two-arm study, inference is based on
#'      \eqn{p_\textrm{treatment} - p_\textrm{control}}. This posterior
#'      probability can be calculated using Monte Carlo beta draws
#'      (`bin_method = "mc"`), a normal approximation (`"normal"`), or numerical
#'      quadrature (`"quadrature"`). Like the risk-difference test, this method
#'      requires complete binary outcomes: censored subjects must either be
#'      followed to `end_of_study`, imputed, or excluded when
#'      `imputed_final = FALSE`.
#'
#'      Two equivalent predictive imputation approaches are available through
#'      `binary_imputation`. With `"event-time"`, the package samples a future
#'      event time conditional on the available event-free follow-up and then
#'      records whether it falls by `end_of_study`. With `"bernoulli"`, it
#'      calculates the same endpoint probability directly. If \eqn{T} is the
#'      observed event-free follow-up, \eqn{T^*} is `end_of_study`,
#'      \eqn{S(t)} is the survival function, and \eqn{H(t)} is the cumulative
#'      hazard, that probability is
#'
#'      \deqn{\Pr(X = 1 \mid T_\mathrm{event} > T)
#'        = \frac{S(T) - S(T^*)}{S(T)}
#'        = 1 - \exp\{-[H(T^*) - H(T)]\}.}
#'
#'      A Bernoulli outcome is drawn with this probability. For a subject not
#'      yet enrolled, \eqn{T = 0}; observed events are retained unchanged.
#'      Because no precise event time is generated, the imputed `time` is set
#'      to `end_of_study` and only the binary `event` status is analyzed. Each
#'      imputation still uses a sampled posterior hazard draw, so uncertainty
#'      in the piecewise-exponential model is retained.
#'
#'  * Frequentist risk difference (`method = "riskdiff-wald"` or
#'    `"riskdiff-fm"`).
#'      Each complete or imputed dataset is reduced to binary event outcomes at
#'      `end_of_study`. The estimated treatment effect is
#'      \eqn{p_\textrm{treatment} - p_\textrm{control}}. `"riskdiff-wald"`
#'      uses the observed arm risks in an unpooled Wald variance.
#'      `"riskdiff-fm"` instead uses maximum likelihood arm risks constrained
#'      by the null difference `h0` in a Farrington-Manning score variance. The
#'      latter remains defined for common sparse tables, including equal-arm
#'      all-zero and all-one outcomes. Both methods report \eqn{1 - P} in
#'      `post_prob_ha` and support all three alternatives. Because they require
#'      complete binary outcomes, lost-to-follow-up subjects are excluded when
#'      `imputed_final = FALSE`. When `imputed_final = TRUE`, estimates and
#'      within-imputation variances from at least two completed datasets are
#'      combined using Rubin's rules, producing a pooled Wald analysis for
#'      either method setting.
#'
#'  * Imputed final analysis (`imputed_final`).
#'      The overall final analysis conducted after accrual is suspended and
#'      follow-up is complete can be analyzed on imputed datasets for Bayesian
#'      methods (`"bayes-surv"` and `"bayes-bin"`), Cox regression, and the
#'      frequentist risk-difference analysis, or on the non-imputed dataset.
#'      Interim prediction completes outcomes that are not yet observed,
#'      whereas final imputation applies only to subjects right-censored because
#'      of loss to follow-up before `end_of_study`. Design evaluations should
#'      prespecify whether the final analysis imputes these outcomes and assess
#'      sensitivity to that choice, particularly when appreciable attrition is
#'      expected. Loss to follow-up is assumed to be non-informative. For Cox
#'      regression the final estimates and
#'      variances are pooled with Rubin's rules. It cannot be used with
#'      `method = "logrank"`.
#'
#'   When imputation is involved, either at interim analyses or through
#'   `imputed_final = TRUE`, the package uses a two-stage impute-then-analyze
#'   procedure. First, the piecewise-exponential model is fitted to the
#'   *observed* time-to-event data and used to complete pending outcomes.
#'   Second, each completed dataset is analyzed using the model selected by
#'   `method`.
#'
#'   For `method = "bayes-bin"`, these are deliberately separate models. The
#'   imputation model has piecewise hazards with Gamma prior `prior_surv` (or
#'   `prior_surv_final` during final imputation), whereas the completed-data
#'   analysis has an event probability at `end_of_study` with Beta prior
#'   `prior_bin`. The Beta prior is not derived from the Gamma prior, and the two
#'   stages are not one joint Bayesian model. Both prior specifications can
#'   therefore affect predictive decisions when outcomes require imputation.
#'   The `"event-time"` and `"bernoulli"` binary-imputation options use the same
#'   piecewise-exponential prediction model and do not change this separation.
#'
#'   For `method = "bayes-surv"`, the second analysis instead forms a fresh
#'   piecewise-exponential posterior from the completed data and the original
#'   survival prior. For frequentist methods (`"logrank"`, `"cox"`,
#'   `"riskdiff-wald"`, and `"riskdiff-fm"`), each completed dataset uses a
#'   standard test rather than a posterior. Imputed Cox and risk-difference
#'   final analyses pool estimates and variances using Rubin's rules.
#'
#'   At each interim look, follow-up times are masked (censored) to reflect
#'   the calendar time of the analysis. The package treats enrollment and
#'   randomization as occurring at the same time. Subjects enrolled at the exact
#'   interim boundary have zero follow-up time. These times are clamped to
#'   `.Machine$double.eps` (approximately \eqn{2.2 \times 10^{-16}}) so
#'   that they contribute negligible but non-zero exposure to the interim
#'   posterior. This affects at most one subject per interim look.
#'
#' @return With `return_trace = FALSE` (the default), a one-row data frame
#'   containing the evaluated design and final trial results, including:
#'
#'   - `N_treatment`: Number of patients enrolled in the treatment arm.
#'   - `N_control`: Number of patients enrolled in the control arm.
#'   - `est_final`: Treatment effect estimated at the final analysis. The final
#'     analysis occurs when either the maximum sample size is reached and
#'     follow-up is complete, or the interim analysis triggered early stopping
#'     of enrollment/accrual and follow-up for those subjects is complete.
#'   - `post_prob_ha`: Posterior probability from the final analysis. If
#'     a Bayesian method uses `imputed_final = TRUE`, this is calculated for
#'     each imputed final-analysis dataset and averaged over `N_impute`
#'     imputations. For an imputed Cox analysis it is \eqn{1 - P} from the
#'     Rubin-pooled Wald test. The same interpretation applies to imputed
#'     risk-difference analyses. For non-imputed frequentist analyses it is
#'     \eqn{1 - P} from the corresponding test.
#'   - `stop_futility`: Logical indicator of whether the trial stopped early for
#'     futility.
#'   - `stop_expected_success`: Logical indicator of whether the trial stopped
#'     early for expected success.
#'   - `stopping_reason`: One of `"expected_success"`, `"futility"`, or
#'     `"maximum_sample_size"`.
#'   - `accrual_stop_time`: Calendar time of the last enrollment in the trial.
#'   - `analysis_ready_time`: Calendar time at which the last enrolled
#'     subject's observed event or censoring becomes available. This excludes
#'     external data-cleaning and database-lock delays.
#'   - `planned_completion_time`: Calendar time at which the last enrolled
#'     subject would complete the full planned follow-up.
#'   - `followup_person_time`: Sum of observed follow-up times across enrolled
#'     subjects.
#'   - `peak_active_followup`: Largest number of enrolled subjects concurrently
#'     under follow-up.
#'
#'   Calendar time is measured from the first patient's enrollment at time
#'   zero. Times use the same units as `lambda_time`, `cutpoints`,
#'   `generation_cutpoints`, and `end_of_study`.
#'
#'   The returned object has a `decision_design` attribute containing
#'   `interim_look`, `Fn`, `Sn`, and the Monte Carlo settings. Thresholds in
#'   this information are stored as one value per interim look (and have
#'   length zero when no interim looks are planned).
#'   A `prior_design` attribute contains the resolved Gamma shape, rate,
#'   mean hazard, and standard deviation for every stage, arm, and interval.
#'
#'   Both return forms have an `arguments` attribute containing a named list of
#'   the evaluated argument values, including defaults. It can be saved with
#'   [saveRDS()] and supplied to a later call with
#'   `do.call(survival_adapt, attr(result, "arguments"))`. Its `prop_loss`
#'   element contains a named value for every simulated arm.
#'   Its `rand_ratio` element is stored in `control`, `treatment` order for
#'   two-arm designs. Its `cutpoints` and `generation_cutpoints` elements retain
#'   the analysis and data-generation partitions, respectively.
#'   For `method = "bayes-bin"`, this metadata explicitly retains the
#'   imputation priors (`prior_surv` and `prior_surv_final`), completed-data
#'   analysis prior (`prior_bin`), and imputation horizon (`end_of_study`). The
#'   separate `prior_design` attribute gives the resolved Gamma parameters by
#'   stage, arm, and interval.
#'
#'   With `return_trace = TRUE`, a `goldilocks_trial` object is returned. Its
#'   `summary` element is the same data frame and its `trace` element has one row
#'   per interim look. `prior_diagnostics` contains the resolved interim and
#'   final priors. `posterior_diagnostics` reports observed and effective
#'   sufficient statistics and conjugate posterior parameters by completed
#'   look, arm, and interval. The trace records calendar time, the number of subjects
#'   actively under follow-up, enrollment and observed events by arm,
#'   predictive probabilities, diagnostic Monte Carlo standard errors and exact
#'   bounds, draw counts, thresholds, the decision and reason, empty-interval
#'   fallback diagnostics, and warnings raised during that look. It deliberately
#'   excludes imputed data sets and posterior draws to keep the output compact.
#'
#' @references
#' Broglio KR, Connor JT, Berry SM. Not too big, not too small: a Goldilocks
#' approach to sample size selection. *Journal of Biopharmaceutical Statistics*,
#' 2014; 24(3): 685–705.
#'
#' @importFrom stats pexp coef
#' @export
#'
#' @examples
#' # RCT with exponential hazard (no piecewise breaks)
#' # Note: the number of imputations is small to enable this example to run
#' #       quickly on CRAN tests. In practice, much larger values are needed.
#' survival_adapt(
#'  hazard_treatment = -log(0.85) / 36,
#'  hazard_control = -log(0.7) / 36,
#'  cutpoints = NULL,
#'  N_total = 600,
#'  lambda = 20,
#'  lambda_time = NULL,
#'  interim_look = 400,
#'  end_of_study = 36,
#'  prior_surv = c(0.1, 0.1),
#'  block = 2,
#'  rand_ratio = c(control = 1, treatment = 1),
#'  prop_loss = 0.30,
#'  alternative = "less",
#'  h0 = 0,
#'  Fn = 0.05,
#'  Sn = 0.9,
#'  prob_ha = 0.975,
#'  N_impute = 10,
#'  N_mcmc = 10,
#'  method = "bayes-surv")
survival_adapt <- function(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = NULL,
  N_total,
  lambda = 0.3,
  lambda_time = NULL,
  interim_look = NULL,
  end_of_study,
  prior_surv = c(0.1, 0.1),
  prior_bin = c(1, 1),
  bin_method = "mc",
  block = 2,
  rand_ratio = c(control = 1, treatment = 1),
  prop_loss = 0,
  alternative = "greater",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.9,
  prob_ha = 0.95,
  N_impute = 500,
  N_mcmc = 1000,
  mc_conf_level = 0.95,
  empty_interval = c("prior", "propagate", "error"),
  method = "logrank",
  imputed_final = FALSE,
  return_trace = FALSE,
  binary_imputation = c("event-time", "bernoulli"),
  prior_surv_final = prior_surv,
  generation_cutpoints = cutpoints
) {
  Call <- match.call()
  Arguments <- capture_arguments(survival_adapt, environment())
  method <- normalize_analysis_method(method)
  Arguments$method <- method
  ##############################################################################
  ### Derive variables
  ##############################################################################

  # Indicator of whether single-arm study
  single_arm <- is.null(hazard_control)

  # Interim look and final look
  analysis_at_enrollnumber <- c(interim_look, N_total)

  # Number of looks
  N_looks <- length(analysis_at_enrollnumber)

  ##############################################################################
  ### Run checks on arguments
  ##############################################################################

  validate_positive_integer_scalar(N_total, "N_total")
  if (!single_arm) {
    rand_ratio <- validate_randomization_args(
      N_total,
      block,
      rand_ratio,
      allocation_name = "rand_ratio"
    )
    Arguments$rand_ratio <- rand_ratio
  }
  prop_loss <- normalize_prop_loss(prop_loss, single_arm)
  Arguments$prop_loss <- prop_loss
  validate_single_probability(prob_ha, "prob_ha")
  validate_positive_integer_scalar(N_impute, "N_impute")
  validate_positive_integer_scalar(N_mcmc, "N_mcmc")
  validate_single_probability(
    mc_conf_level,
    "mc_conf_level",
    upper_open = TRUE
  )
  if (mc_conf_level <= 0.5) {
    stop("'mc_conf_level' must be greater than 0.5 and less than 1")
  }
  validate_cutpoints(cutpoints)
  validate_endpoint_time(end_of_study, cutpoints, "end_of_study")
  n_intervals <- length(cutpoints) + 1L
  prior_surv <- normalize_gamma_prior(
    prior_surv,
    n_intervals = n_intervals,
    single_arm = single_arm,
    name = "prior_surv"
  )
  prior_surv_final <- normalize_gamma_prior(
    prior_surv_final,
    n_intervals = n_intervals,
    single_arm = single_arm,
    name = "prior_surv_final"
  )
  prior_design <- rbind(
    gamma_prior_diagnostics(
      prior_surv = prior_surv,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      single_arm = single_arm,
      stage = "interim"
    ),
    gamma_prior_diagnostics(
      prior_surv = prior_surv_final,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      single_arm = single_arm,
      stage = "final"
    )
  )
  rownames(prior_design) <- NULL
  empty_interval <- match.arg(empty_interval)
  binary_imputation <- match.arg(binary_imputation)
  Arguments$empty_interval <- empty_interval
  Arguments$binary_imputation <- binary_imputation
  validate_logical_scalar(return_trace, "return_trace")
  validate_analysis_configuration(
    method,
    alternative,
    single_arm,
    imputed_final
  )
  if (
    imputed_final &&
      method %in% c("cox", "riskdiff-wald", "riskdiff-fm") &&
      N_impute < 2
  ) {
    stop(
      "Frequentist final-analysis imputation requires at least two imputations ",
      "to apply Rubin's rules"
    )
  }

  if (!is.null(interim_look)) {
    validate_interim_looks(
      interim_look,
      N_total,
      min_look = if (single_arm) NULL else max(block)
    )
  }

  validate_h0(h0, method, single_arm)

  # Check: Bayesian binomial test arguments
  if (method == "bayes-bin") {
    validate_bayes_binomial_args(prior_bin, bin_method, N_mcmc)
  }

  # A scalar threshold is broadcast; any non-scalar must supply exactly one
  # value per interim look. This normalization is shared by both decision rules
  # and retained in the returned design metadata.
  N_interims <- N_looks - 1L
  Sn <- normalize_interim_threshold(Sn, N_interims, "Sn")
  Fn <- normalize_interim_threshold(
    Fn,
    N_interims,
    "Fn",
    null_disables = TRUE
  )
  check_futility <- N_interims > 0L && any(Fn != 0)
  decision_design <- list(
    interim_look = interim_look,
    Fn = Fn,
    Sn = Sn,
    N_impute = N_impute,
    N_mcmc = N_mcmc,
    mc_conf_level = mc_conf_level
  )

  # Posterior samples are not used by frequentist methods
  if (!method %in% c("bayes-surv", "bayes-bin")) {
    N_mcmc <- 1
  }

  ##############################################################################
  ### Simulate complete data for trial
  ##############################################################################

  data_total <- sim_comp_data(
    hazard_treatment = hazard_treatment,
    hazard_control = hazard_control,
    generation_cutpoints = generation_cutpoints,
    N_total = N_total,
    lambda = lambda,
    lambda_time = lambda_time,
    end_of_study = end_of_study,
    block = block,
    rand_ratio = rand_ratio,
    prop_loss = prop_loss
  )

  ##############################################################################
  ### Evaluate trial at each interim analysis
  ##############################################################################

  # Assigning stop_futility and stop_expected_success
  stop_futility <- 0
  stop_expected_success <- 0
  trace_rows <- if (return_trace) {
    vector("list", max(N_looks - 1L, 0L))
  } else {
    NULL
  }
  posterior_diagnostic_rows <- if (return_trace) {
    vector("list", max(N_looks - 1L, 0L))
  } else {
    NULL
  }

  if (N_looks > 1) {
    for (i in 1:(N_looks - 1)) {
      # Analysis at the `analysis_at_enrollnumber` look (not incl. final look)
      # Indicators for subject type:
      # - subject_enrolled:        subject has data present in the current look
      # - subject_impute_success:  subject has data present in the current look but has not
      #                            reached end of study due to staged entry or LTFU;
      #                            needs imputation
      # - subject_impute_futility: subject has no data present in the current look;
      #                            needs imputation
      # - time_from_rand_at_look:  time from enrollment/randomization to sample size look
      #                            e.g. if patient enrolled month 3, but look occurs month 7,
      #                            then patient could potentially be observed for 4 months

      data_interim <- within(data_total, {
        subject_enrolled <- (id <= analysis_at_enrollnumber[i])
        subject_impute_futility <- !subject_enrolled
        time_from_rand_at_look <- enrollment[analysis_at_enrollnumber[i]] -
          enrollment
        subject_impute_success <-
          # Had event, but has not occurred yet (based on interim look)
          ((event == 1) * (time_from_rand_at_look < time) & subject_enrolled) |
          # Event-free and not had opportunity to complete full follow
          ((event == 0) *
            (time_from_rand_at_look < end_of_study) &
            subject_enrolled) |
          (loss_to_fu & subject_enrolled)
      })

      # Mask the data at time of look
      # Note: subjects at the exact interim boundary have
      # time_from_rand_at_look = 0, yielding time = 0 after masking.
      # Clamp to .Machine$double.eps so the boundary subject contributes
      # negligible but non-zero exposure to the interim posterior.
      data_interim <- within(data_interim, {
        time <- pmax(pmin(time, time_from_rand_at_look), .Machine$double.eps)
        event <- ifelse(subject_impute_success, 0, event)
      })

      look_time <- data_total$enrollment[analysis_at_enrollnumber[i]]
      interim_result <- evaluate_interim_decision(
        data_interim = data_interim,
        look = i,
        planned_N = analysis_at_enrollnumber[i],
        calendar_time = look_time,
        active_followup = active_followup_at(data_total, look_time),
        end_of_study = end_of_study,
        cutpoints = cutpoints,
        single_arm = single_arm,
        prior_surv = prior_surv,
        prior_bin = prior_bin,
        bin_method = bin_method,
        alternative = alternative,
        h0 = h0,
        Fn = Fn[i],
        Sn = Sn[i],
        prob_ha = prob_ha,
        N_impute = N_impute,
        N_mcmc = N_mcmc,
        mc_conf_level = mc_conf_level,
        empty_interval = empty_interval,
        method = method,
        binary_imputation = binary_imputation,
        check_futility = check_futility
      )
      ppp_success <- interim_result$ppp_success
      ppp_success_at_max <- interim_result$ppp_success_at_max
      decision <- interim_result$decision
      decision_reason <- interim_result$decision_reason

      if (return_trace) {
        trace_rows[[i]] <- interim_result$trace
        posterior_diagnostics <- interim_result$diagnostics$posterior
        posterior_diagnostics$look <- as.integer(i)
        posterior_diagnostics$planned_N <- as.integer(
          analysis_at_enrollnumber[i]
        )
        posterior_diagnostics$calendar_time <- look_time
        posterior_diagnostic_rows[[i]] <- posterior_diagnostics[c(
          "look",
          "planned_N",
          "calendar_time",
          setdiff(
            names(posterior_diagnostics),
            c("look", "planned_N", "calendar_time")
          )
        )]
      }

      if (decision == "stop_expected_success") {
        stop_expected_success <- 1
        stage_trial_stopped <- analysis_at_enrollnumber[i]
        break # No further SS looks
      }

      # Test if futility success criteria is met
      if (decision == "stop_futility") {
        stop_futility <- 1
        stage_trial_stopped <- analysis_at_enrollnumber[i]
        break # No further SS looks
      }

      # Stop study if at last interim look
      if (analysis_at_enrollnumber[i + 1] == N_total) {
        stage_trial_stopped <- analysis_at_enrollnumber[i + 1]
        break # No further SS looks
      }
    }

    # Number of patients enrolled at trial stop
    N_enrolled <- nrow(data_interim[data_interim$id <= stage_trial_stopped, ])
  } else {
    # Assigning stage trial stopped given no interim look
    N_enrolled <- N_total
    stage_trial_stopped <- N_total
    stop_futility <- 0
    stop_expected_success <- 0
    ppp_success <- NA
  }

  ##############################################################################
  ### Final analysis (after enrollment complete)
  ##############################################################################

  # All patients that have made it to the end of study
  # - complete follow-up (except any censoring)
  data_final <- subset(data_total, id <= stage_trial_stopped)
  data_final <- within(data_final, {
    time_from_rand_at_look <- enrollment[stage_trial_stopped] - enrollment
    subject_impute_success <- ((event == 0) & (time < end_of_study))
  })

  results_final <- analyse_final(
    data_in = data_final,
    cutpoints = cutpoints,
    prior_surv_final = prior_surv_final,
    N_mcmc = N_mcmc,
    single_arm = single_arm,
    imputed_final = imputed_final,
    method = method,
    N_impute = N_impute,
    alternative = alternative,
    h0 = h0,
    prior_bin = prior_bin,
    bin_method = bin_method,
    binary_imputation = if (
      method %in% c("bayes-bin", "riskdiff-wald", "riskdiff-fm")
    ) {
      binary_imputation
    } else {
      "event-time"
    },
    empty_interval = empty_interval,
    end_of_study = end_of_study
  )

  post_paa <- results_final[1]
  est_final <- results_final[2]

  N_treatment <- sum(data_final$treatment == 1) # Total analyzed: treatment
  N_control <- sum(data_final$treatment == 0) # Total analyzed: control
  calendar_metrics <- trial_calendar_metrics(data_final, end_of_study)
  stopping_reason <- calendar_stopping_reason(
    stop_futility,
    stop_expected_success
  )

  ##############################################################################
  ### Output
  ##############################################################################

  results <- data.frame(
    prob_threshold = prob_ha,
    margin = h0,
    alternative = alternative,
    N_treatment = N_treatment,
    N_control = N_control,
    N_enrolled = N_treatment + N_control,
    N_max = N_total,
    post_prob_ha = post_paa,
    est_final = est_final,
    ppp_success = ppp_success,
    stop_futility = stop_futility,
    stop_expected_success = stop_expected_success,
    stopping_reason = stopping_reason,
    accrual_stop_time = calendar_metrics$accrual_stop_time,
    analysis_ready_time = calendar_metrics$analysis_ready_time,
    planned_completion_time = calendar_metrics$planned_completion_time,
    followup_person_time = calendar_metrics$followup_person_time,
    peak_active_followup = calendar_metrics$peak_active_followup
  )
  enrollment_design <- new_enrollment_design(
    lambda = lambda,
    N_total = N_total,
    lambda_time = lambda_time,
    interim_look = interim_look,
    end_of_study = end_of_study
  )
  attr(results, "enrollment_design") <- enrollment_design
  attr(results, "decision_design") <- decision_design
  attr(results, "prior_design") <- prior_design
  attr(results, "arguments") <- Arguments

  if (return_trace) {
    out <- list(
      summary = results,
      trace = new_trial_trace(trace_rows),
      prior_diagnostics = prior_design,
      posterior_diagnostics = new_interim_posterior_diagnostics(
        posterior_diagnostic_rows
      ),
      call = Call
    )
    class(out) <- "goldilocks_trial"
    attr(out, "enrollment_design") <- enrollment_design
    attr(out, "decision_design") <- decision_design
    attr(out, "prior_design") <- prior_design
    attr(out, "arguments") <- Arguments
    return(out)
  }

  return(results)
}
