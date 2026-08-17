#' @title Simulate and analyze one Goldilocks adaptive trial
#'
#' @inheritParams sim_comp_data
#' @param interim_look vector. Sample size for each interim look. Note: the
#'   maximum sample size should not be included. For two-arm designs, each
#'   interim look must be at least the (largest) block size (see `block`),
#'   ensuring both treatment groups are present at every interim analysis; a
#'   smaller look could enroll subjects from one treatment group only, leaving
#'   the interim posterior undefined for the missing group.
#' @param prior_surv numeric vector or matrix. Gamma prior for the
#'   piecewise-exponential hazards used during interim prediction. A length-two
#'   vector supplies shape and rate and is broadcast across all intervals. A
#'   `2` by `length(cutpoints) + 1` matrix supplies interval-specific values,
#'   with shapes in row 1, rates in row 2, and columns ordered from the earliest
#'   to the latest interval. The same interval prior is applied to both
#'   treatment groups. Rates must use the same time unit as event times,
#'   exposure, and cutpoints. The default is `c(0.1, 0.1)`.
#' @param prior_surv_final numeric vector or matrix. Gamma prior used for
#'   final-stage piecewise-exponential imputation and, for `method =
#'   "bayes-surv"`, final analysis. It accepts the same forms as `prior_surv`
#'   and defaults to `prior_surv`, preserving the historical behavior.
#' @param prior_bin vector. Prior distribution for the event probability when
#'   `method = "bayes-bin"`. The two values are the shape parameters of the
#'   `Beta(a, b)` prior. The same prior is applied to both treatment arms.
#' @param bin_method character. Method used to calculate the posterior
#'   probability for `method = "bayes-bin"`, must be one of `"mc"` (Monte Carlo
#'   sampling), `"normal"` (normal approximation), or `"quadrature"` (numerical
#'   integration). The default is `"mc"`.
#' @param binary_imputation character. Predictive imputation approach for
#'   `method = "bayes-bin"` or `method = "riskdiff"`. `"event-time"` (the default)
#'   draws a conditional piecewise-exponential event time and reduces it to
#'   event status at `end_of_study`. `"bernoulli"` draws the endpoint status
#'   directly from its conditional event probability. This argument is ignored
#'   for time-to-event analysis methods.
#' @param alternative character. The string specifying the alternative
#'   hypothesis, must be one of `"greater"` (default), `"less"` or
#'   `"two.sided"`. One-sided alternatives (`"greater"` and `"less"`) are
#'   supported for `method = "bayes-surv"` and `method = "bayes-bin"`. All three
#'   options are supported for `method = "logrank"`, `method = "cox"`, and
#'   `method = "riskdiff"`. For
#'   survival outcomes, `"less"` corresponds to the treatment arm having a lower
#'   cumulative incidence (i.e., treatment is beneficial), and `"greater"`
#'   corresponds to the treatment arm having a higher cumulative incidence.
#' @param h0 single finite numeric null hypothesis value or margin. Default is
#'   `h0 = 0`. For Bayesian analyses, `h0` must lie in `[0, 1]` for a
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
#'   * When `method = "riskdiff"`, `h0` is the null value of
#'     \eqn{p_\textrm{treatment} - p_\textrm{control}} and must lie in
#'     `[-1, 1]`.
#'   * When `method = "logrank"`, only `h0 = 0` is supported; this denotes the
#'     usual equal-survival null. Nonzero values are rejected because the
#'     standard log-rank statistic does not implement a nonzero effect margin.
#' @param Fn vector of values between 0 and 1. Each element is the probability
#'   threshold to stop at the \eqn{i}-th look early for futility. If there are
#'   no interim looks (i.e. `interim_look = NULL`), then `Fn` is not
#'   used in the simulations or analysis. Set `Fn = 0` to disable futility
#'   monitoring; `Fn = NULL` has the same effect. Supply either one value,
#'   which is broadcast to every interim look, or exactly one value per
#'   `interim_look`. Other lengths are rejected rather than recycled.
#' @param Sn vector of values between 0 and 1. Each element is the probability
#'   threshold to stop at the \eqn{i}-th look early for expected success. If
#'   there are no interim looks (i.e. `interim_look = NULL`), then
#'   `Sn` is not used in the simulations or analysis. Supply either one value,
#'   which is broadcast to every interim look, or exactly one value per
#'   `interim_look`. Other lengths are rejected rather than recycled.
#' @param prob_ha scalar value between 0 and 1. Probability threshold of alternative
#'   hypothesis.
#' @param N_impute integer. Number of predictive imputations used at each
#'   interim look and, when requested, for final multiple imputation. The
#'   default is 500. An imputed Cox or risk-difference final analysis requires
#'   at least two.
#' @param N_mcmc integer. Number of posterior samples used within each
#'   `method = "bayes-surv"` and by `method = "bayes-bin"` when
#'   `bin_method = "mc"`. The default is 1000.
#' @param mc_conf_level probability strictly between 0.5 and 1. Confidence
#'   level for one-sided exact binomial bounds that guard finite Monte Carlo
#'   interim decisions. A completed-data Bayesian analysis is counted as
#'   successful only when its lower Monte Carlo bound exceeds `prob_ha`.
#'   Expected-success stopping requires the lower bound for the outer
#'   predictive estimate to exceed `Sn`; futility stopping requires its upper
#'   bound to be below `Fn`. Equality or an unresolved bound continues
#'   enrollment. The default is 0.95.
#' @param empty_interval character. Policy for empty piecewise-exponential
#'   intervals in `method = "bayes-surv"` posterior calculations. An empty
#'   interval is an interval with no exposed subjects in a treatment arm at the
#'   analysis time. `"prior"` (the default) leaves the interval at zero exposure
#'   time and zero events, so its posterior is driven only by its assigned
#'   survival prior. `"propagate"` is a legacy heuristic that copies exposure
#'   time and event counts from the nearest non-empty interval in the same
#'   treatment arm and emits a warning. `"error"` stops when any empty interval
#'   is found.
#' @param method character. For an imputed data set (or the final data set after
#'   follow-up is complete), whether the analysis should be a log-rank
#'   (`method = "logrank"`) test, Cox proportional hazards regression model
#'   Wald test (`method = "cox"`), a fully-Bayesian piecewise-exponential
#'   analysis (`method = "bayes-surv"`), a Bayesian beta-binomial analysis of
#'   complete binary outcomes (`method = "bayes-bin"`), or a frequentist
#'   risk-difference Wald test of complete binary outcomes (`method =
#'   "riskdiff"`). See Details section.
#' @param imputed_final logical. Should the final analysis (after all subjects
#'   have been followed-up to the study end) be based on imputed outcomes for
#'   subjects who were LTFU (i.e. right-censored with time less than
#'   `end_of_study`)? Default is `FALSE`, which means that the final analysis
#'   incorporates right-censoring. With `method = "cox"` or `method =
#'   "riskdiff"`, setting this to `TRUE` analyzes each imputed dataset and pools
#'   the scalar treatment effects and variances using Rubin's rules; this
#'   requires `N_impute >= 2`. Imputed final analyses remain unavailable for
#'   `method = "logrank"`.
#' @param return_trace logical. Should the interim decision path be returned in
#'   addition to the usual final summary? The default, FALSE, returns the
#'   historical one-row data frame. When TRUE, the result is a
#'   goldilocks_trial object with summary, trace, and call elements.
#'
#' @details Implements the Goldilocks design method described in Broglio et al.
#'   (2014). At each interim analysis, two probabilities are computed:
#'
#'   1. **The posterior predictive probability of eventual success.** This is
#'      calculated as the proportion of imputed datasets at the *current* sample
#'      size that would go on to be success at the specified threshold. At each
#'      interim analysis its one-sided exact lower Monte Carlo confidence bound
#'      is compared to the corresponding element of `Sn`, and if the bound
#'      exceeds the threshold,
#'      accrual/enrollment is suspended and the outstanding follow-up allowed to
#'      complete before conducting the pre-specified final analysis.
#'
#'   2. **The posterior predictive probability of final success**. This is
#'      calculated as the proportion of imputed datasets at the *maximum*
#'      threshold that would go on to be successful. Similar to above, it is
#'      compared to the corresponding element of `Fn`, and if its one-sided
#'      exact upper Monte Carlo confidence bound is below the threshold,
#'      accrual/enrollment is suspended and the
#'      trial terminated. Typically this would be a binding decision. If it is
#'      not a binding decision, then one should also explore the simulations
#'      with `Fn = 0`.
#'
#'   Hence, at each interim analysis look, 3 decisions are allowed:
#'
#'   1. **Stop for expected success**
#'   2. **Stop for futility**
#'   3. **Continue to enroll** new subjects, or if at maximum sample size,
#'      proceed to final analysis.
#'
#'   At each interim (and final) analysis methods as:
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
#'   * Bayesian absolute difference (`method = "bayes-surv"`).
#'      Each imputed dataset is used to update the conjugate Gamma prior
#'      (defined by `prior_surv` at interim looks and `prior_surv_final` at the
#'      final stage), yielding a posterior distribution for the
#'      piecewise exponential rate parameters. In turn, the posterior
#'      distribution of the cumulative incidence function (\eqn{1 - S(t)}, where
#'      \eqn{S(t)} is the survival function) evaluated at time
#'      `end_of_study` is calculated. If a single-arm study, then this
#'      summarizes the treatment effect, else, if a two-armed study, the
#'      independent posteriors are used to estimate the posterior distribution
#'      of the difference. A posterior probability is calculated according to the
#'      specification of the test type (`alternative`) and the value of the null
#'      hypothesis (`h0`).
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
#'  * Frequentist risk difference (`method = "riskdiff"`).
#'      Each complete or imputed dataset is reduced to binary event outcomes at
#'      `end_of_study`. The estimated treatment effect is
#'      \eqn{p_\textrm{treatment} - p_\textrm{control}}, with an unpooled
#'      binomial variance. A Wald test compares this estimate with `h0`, and
#'      \eqn{1 - P} is reported in `post_prob_ha`. All three alternatives are
#'      supported. Because the test requires complete binary outcomes,
#'      lost-to-follow-up subjects are excluded when `imputed_final = FALSE`.
#'      When `imputed_final = TRUE`, estimates and within-imputation variances
#'      from at least two completed datasets are combined using Rubin's rules.
#'
#'  * Imputed final analysis (`imputed_final`).
#'      The overall final analysis conducted after accrual is suspended and
#'      follow-up is complete can be analyzed on imputed datasets for Bayesian
#'      methods (`"bayes-surv"` and `"bayes-bin"`), Cox regression, and the
#'      frequentist risk-difference analysis, or on the non-imputed dataset.
#'      Since the imputations/predictions used
#'      during the interim analyses assume all subjects are imputed (since loss
#'      to follow-up is not yet known), it would seem most appropriate to
#'      conduct the trial in the same manner, especially if loss to follow-up
#'      rates are appreciable. Note, this only applies to subjects who are
#'      right-censored due to loss to follow-up, which we assume is a
#'      non-informative process. For Cox regression the final estimates and
#'      variances are pooled with Rubin's rules. It cannot be used with
#'      `method = "logrank"`.
#'
#'   When `method = "bayes-surv"` or `method = "bayes-bin"` and imputation is
#'   involved (either at interim
#'   analyses or via `imputed_final = TRUE`), a two-stage posterior
#'   procedure is used. First, the posterior distribution of the piecewise
#'   hazard rates is estimated from the *observed* data and used to draw
#'   imputed event times for censored subjects. Second, a *new* posterior is
#'   estimated from the combined observed and imputed data: the
#'   piecewise-exponential posterior for `method = "bayes-surv"` or the beta
#'   posterior for `method = "bayes-bin"`. This posterior is used for
#'   inference. This is consistent with the predictive probability framework
#'   described in Broglio et al. (2014), but users should be aware that the
#'   imputation model's posterior influences the analysis posterior. For
#'   frequentist methods (`"logrank"`, `"cox"`, `"riskdiff"`), each completed
#'   dataset uses a standard test rather than a posterior, so this feedback loop
#'   does not arise. Imputed Cox final analyses then pool the completed-data
#'   estimates and variances using Rubin's rules. Imputed risk-difference final
#'   analyses use the same scalar combining rule.
#'
#'   At each interim look, follow-up times are masked (censored) to reflect
#'   the calendar time of the analysis. The package treats enrollment and
#'   randomization as occurring at the same time. Subjects enrolled at the exact
#'   interim boundary have zero follow-up time. These times are clamped to
#'   `.Machine$double.eps` (approximately \eqn{2.2 \times 10^{-16}}) so
#'   that they contribute negligible but non-zero exposure to the interim
#'   posterior. This affects at most one subject per interim look.
#'
#' @return With return_trace = FALSE (the default), a data frame containing
#'   some input parameters (arguments) as well as statistics from the analysis,
#'   including:
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
#'
#'   The returned object has a `decision_design` attribute containing
#'   `interim_look`, `Fn`, `Sn`, and the Monte Carlo settings. Thresholds in
#'   this metadata are normalized to one value per interim look (and have
#'   length zero when no interim looks are planned).
#'
#'   With return_trace = TRUE, a goldilocks_trial object is returned. Its
#'   summary element is the same data frame and its trace element has one row
#'   per interim look. The trace records enrollment and observed events by arm,
#'   predictive probabilities, Monte Carlo standard errors and exact bounds,
#'   draw counts, thresholds, the decision and reason, empty-interval fallback
#'   diagnostics, and warnings raised during that look. It deliberately excludes imputed data
#'   sets and posterior draws to keep the output compact.
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
#'  rand_ratio = c(1, 1),
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
  rand_ratio = c(1, 1),
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
  prior_surv_final = prior_surv
) {
  Call <- match.call()
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
  validate_single_probability(prop_loss, "prop_loss")
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
  n_intervals <- length(cutpoints) + 1L
  prior_surv <- normalize_gamma_prior(
    prior_surv,
    n_intervals = n_intervals,
    name = "prior_surv"
  )
  prior_surv_final <- normalize_gamma_prior(
    prior_surv_final,
    n_intervals = n_intervals,
    name = "prior_surv_final"
  )
  empty_interval <- match.arg(empty_interval)
  binary_imputation <- match.arg(binary_imputation)
  validate_logical_scalar(return_trace, "return_trace")
  validate_analysis_configuration(
    method,
    alternative,
    single_arm,
    imputed_final
  )
  if (imputed_final && method %in% c("cox", "riskdiff") && N_impute < 2) {
    stop(
      "Frequentist final-analysis imputation requires at least two imputations ",
      "to apply Rubin's rules"
    )
  }

  if (!single_arm) {
    validate_randomization_args(N_total, block, rand_ratio)
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
    cutpoints = cutpoints,
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

      # Carry out interim analysis on patients with complete data only
      # - Set-up new 'data' data frame
      data <- subset(
        data_interim,
        subset = subject_enrolled,
        select = c(time, event, treatment)
      )

      # Capture warnings for this look while preserving their usual output.
      warning_state <- new.env(parent = emptyenv())
      warning_state$messages <- character()
      warning_state$empty_interval_fallbacks <- character()
      capture_warning <- function(warning) {
        if (return_trace) {
          warning_state$messages <- unique(c(
            warning_state$messages,
            conditionMessage(warning)
          ))
        }
      }
      capture_empty_interval <- function(condition) {
        if (return_trace) {
          detail <- condition$details
          entries <- paste0(
            condition$policy,
            ": treatment=",
            detail$treatment,
            ", interval=",
            detail$interval
          )
          warning_state$empty_interval_fallbacks <- unique(c(
            warning_state$empty_interval_fallbacks,
            entries
          ))
        }
      }

      # Posterior distribution of lambdas: current data
      post_lambda <- withCallingHandlers(
        posterior(
          data = data,
          cutpoints = cutpoints,
          prior_surv = prior_surv,
          N_mcmc = N_impute,
          single_arm = single_arm,
          empty_interval = empty_interval
        ),
        warning = capture_warning,
        goldilocks_empty_interval = capture_empty_interval
      )

      ##########################################################################
      ### Loop over multiple imputations
      ##########################################################################

      futility_test <- 0
      expected_success_test <- 0
      inner_mc_uncertain_now <- 0L
      inner_mc_uncertain_max <- 0L
      for (j in 1:N_impute) {
        h <- post_lambda[j, , , drop = FALSE]

        stop_check <- withCallingHandlers(
          test_stop_success(
            data = data_interim,
            hazard = h,
            end_of_study = end_of_study,
            cutpoints = cutpoints,
            single_arm = single_arm,
            prior_surv = prior_surv,
            N_mcmc = N_mcmc,
            method = method,
            alternative = alternative,
            h0 = h0,
            prior_bin = prior_bin,
            bin_method = bin_method,
            binary_imputation = if (method %in% c("bayes-bin", "riskdiff")) {
              binary_imputation
            } else {
              "event-time"
            },
            empty_interval = empty_interval,
            check_futility = check_futility,
            mc_conf_level = mc_conf_level,
            prob_ha = prob_ha
          ),
          warning = capture_warning,
          goldilocks_empty_interval = capture_empty_interval
        )

        # Bayesian completed-data analyses use an inner exact Monte Carlo
        # bound; deterministic frequentist and quadrature analyses retain their
        # strict point comparison with prob_ha.
        analysis_now <- stop_check$classification_now
        if (analysis_now$crossed) {
          expected_success_test <- expected_success_test + 1
        }
        inner_mc_uncertain_now <- inner_mc_uncertain_now +
          as.integer(analysis_now$uncertain)

        if (check_futility) {
          analysis_max <- stop_check$classification_max
          if (analysis_max$crossed) {
            futility_test <- futility_test + 1
          }
          inner_mc_uncertain_max <- inner_mc_uncertain_max +
            as.integer(analysis_max$uncertain)
        }
      }

      # The outer Monte Carlo estimands are probabilities that a completed
      # dataset meets its final success rule. Exact one-sided bounds prevent a
      # coarse point estimate from triggering an unsupported trial decision.
      ppp_now_mc <- monte_carlo_probability_summary(
        successes = expected_success_test,
        draws = N_impute,
        threshold = Sn[i],
        direction = "greater",
        confidence = mc_conf_level
      )
      ppp_success <- ppp_now_mc$estimate
      ppp_max_mc <- if (check_futility) {
        monte_carlo_probability_summary(
          successes = futility_test,
          draws = N_impute,
          threshold = Fn[i],
          direction = "less",
          confidence = mc_conf_level
        )
      } else {
        NULL
      }
      ppp_success_at_max <- if (check_futility) {
        ppp_max_mc$estimate
      } else {
        NA_real_
      }

      decision <- if (ppp_now_mc$crossed) {
        "stop_expected_success"
      } else if (check_futility && ppp_max_mc$crossed) {
        "stop_futility"
      } else {
        "continue"
      }
      decision_reason <- if (decision == "stop_expected_success") {
        "expected_success_lower_bound_above_threshold"
      } else if (decision == "stop_futility") {
        "futility_upper_bound_below_threshold"
      } else if (
        inner_mc_uncertain_now > 0L || inner_mc_uncertain_max > 0L
      ) {
        "continue_inner_monte_carlo_uncertain"
      } else if (ppp_now_mc$reason == "monte_carlo_uncertain") {
        "continue_expected_success_monte_carlo_uncertain"
      } else if (
        check_futility && ppp_max_mc$reason == "monte_carlo_uncertain"
      ) {
        "continue_futility_monte_carlo_uncertain"
      } else {
        "continue_thresholds_not_crossed"
      }

      if (return_trace) {
        trace_rows[[i]] <- data.frame(
          look = i,
          planned_N = analysis_at_enrollnumber[i],
          calendar_time = data_total$enrollment[analysis_at_enrollnumber[i]],
          N_enrolled = nrow(data),
          N_treatment = sum(data$treatment == 1),
          N_control = sum(data$treatment == 0),
          events_treatment = sum(data$event[data$treatment == 1]),
          events_control = sum(data$event[data$treatment == 0]),
          N_pending = sum(
            data_interim$subject_enrolled &
              data_interim$subject_impute_success
          ),
          N_not_enrolled = sum(data_interim$subject_impute_futility),
          ppp_stop_now = ppp_success,
          ppp_stop_now_mcse = ppp_now_mc$mcse,
          ppp_stop_now_lower = ppp_now_mc$lower,
          ppp_stop_now_upper = ppp_now_mc$upper,
          ppp_stop_now_draws = ppp_now_mc$draws,
          success_threshold = Sn[i],
          ppp_success_at_max = ppp_success_at_max,
          ppp_success_at_max_mcse = if (check_futility) {
            ppp_max_mc$mcse
          } else {
            NA_real_
          },
          ppp_success_at_max_lower = if (check_futility) {
            ppp_max_mc$lower
          } else {
            NA_real_
          },
          ppp_success_at_max_upper = if (check_futility) {
            ppp_max_mc$upper
          } else {
            NA_real_
          },
          ppp_success_at_max_draws = if (check_futility) {
            ppp_max_mc$draws
          } else {
            NA_integer_
          },
          futility_threshold = if (check_futility) Fn[i] else NA_real_,
          inner_mc_uncertain_stop_now = inner_mc_uncertain_now,
          inner_mc_uncertain_success_at_max = if (check_futility) {
            inner_mc_uncertain_max
          } else {
            NA_integer_
          },
          decision = decision,
          decision_reason = decision_reason,
          empty_interval_fallback_count = length(
            warning_state$empty_interval_fallbacks
          ),
          empty_interval_fallbacks = paste(
            warning_state$empty_interval_fallbacks,
            collapse = " | "
          ),
          warning_count = length(warning_state$messages),
          warning_messages = paste(warning_state$messages, collapse = " | "),
          stringsAsFactors = FALSE
        )
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

  results_final <- test_final(
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
    binary_imputation = if (method %in% c("bayes-bin", "riskdiff")) {
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
    stop_expected_success = stop_expected_success
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

  if (return_trace) {
    out <- list(
      summary = results,
      trace = new_trial_trace(trace_rows),
      call = Call
    )
    class(out) <- "goldilocks_trial"
    attr(out, "enrollment_design") <- enrollment_design
    attr(out, "decision_design") <- decision_design
    return(out)
  }

  return(results)
}
