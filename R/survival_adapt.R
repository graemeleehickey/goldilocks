#' @title Simulate and execute a single adaptive clinical trial design with a
#'   time-to-event endpoint
#'
#' @inheritParams sim_comp_data
#' @param interim_look vector. Sample size for each interim look. Note: the
#'   maximum sample size should not be included. For two-arm designs, each
#'   interim look must be at least the (largest) block size (see `block`),
#'   ensuring both treatment groups are present at every interim analysis; a
#'   smaller look could enroll subjects from one treatment group only, leaving
#'   the interim posterior undefined for the missing group.
#' @param prior vector. The prior distributions for the piecewise hazard rate
#'   parameters are each \eqn{Gamma(a_0, b_0)}, where \eqn{a_0} is the shape
#'   parameter and \eqn{b_0} is the rate parameter (i.e., the inverse of the
#'   scale). This follows R's [stats::rgamma()] parameterization. The
#'   same prior is applied to all piecewise intervals and to both treatment
#'   groups. The default non-informative prior distribution used is
#'   `Gamma(0.1, 0.1)`, which is specified by setting `prior = c(0.1, 0.1)`.
#' @param bin_prior vector. Prior distribution for the event probability when
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
#'   * The argument is ignored for `method = "logrank"` after its finite-value
#'     validation; the usual equal-survival null hypothesis is used.
#' @param Fn vector of values between 0 and 1. Each element is the probability
#'   threshold to stop at the \eqn{i}-th look early for futility. If there are
#'   no interim looks (i.e. `interim_look = NULL`), then `Fn` is not
#'   used in the simulations or analysis. Set `Fn = 0` to disable futility
#'   monitoring. The length of `Fn` should be the same as
#'   `interim_look`, else the values are recycled.
#' @param Sn vector of values between 0 and 1. Each element is the probability
#'   threshold to stop at the \eqn{i}-th look early for expected success. If
#'   there are no interim looks (i.e. `interim_look = NULL`), then
#'   `Sn` is not used in the simulations or analysis. The length of
#'   `Sn` should be the same as `interim_look`, else the values are
#'   recycled.
#' @param prob_ha scalar value between 0 and 1. Probability threshold of alternative
#'   hypothesis.
#' @param N_impute integer. Number of imputations for Monte Carlo simulation of
#'   missing data. An imputed Cox or risk-difference final analysis requires at
#'   least two.
#' @param N_mcmc integer. Number of posterior samples used by
#'   `method = "bayes-surv"` and by `method = "bayes-bin"` when
#'   `bin_method = "mc"`.
#' @param empty_interval character. Policy for empty piecewise-exponential
#'   intervals in `method = "bayes-surv"` posterior calculations. An empty
#'   interval is an interval with no exposed subjects in a treatment arm at the
#'   analysis time. `"propagate"` (the default, matching earlier package
#'   behavior) copies exposure time and event counts from the nearest non-empty
#'   interval in the same treatment arm and emits a warning. `"prior"` leaves
#'   the interval at zero exposure time and zero events, so its posterior is
#'   driven only by `prior`. `"error"` stops when any empty interval is found.
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
#'      interim analysis it is compared to the corresponding element of
#'      `Sn`, and if it exceeds the threshold,
#'      accrual/enrollment is suspended and the outstanding follow-up allowed to
#'      complete before conducting the pre-specified final analysis.
#'
#'   2. **The posterior predictive probability of final success**. This is
#'      calculated as the proportion of imputed datasets at the *maximum*
#'      threshold that would go on to be successful. Similar to above, it is
#'      compared to the corresponding element of `Fn`, and if it
#'      is less than the threshold, accrual/enrollment is suspended and the
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
#'      (defined by `prior`), yielding a posterior distribution for the
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
#'      The default, `"propagate"`, preserves historical package behavior by
#'      borrowing sufficient statistics from the nearest non-empty interval
#'      within the same treatment arm. This is operationally stable but
#'      statistically consequential because the empty interval's posterior is
#'      informed by adjacent observed data. `"prior"` instead leaves the empty
#'      interval prior-driven, making the absence of interval data explicit.
#'      `"error"` is strict and stops the simulation or analysis when an empty
#'      interval is encountered.
#'
#'   * Bayesian beta-binomial analysis (`method = "bayes-bin"`).
#'      Each complete or imputed dataset is reduced to binary event outcomes at
#'      `end_of_study`. A conjugate `Beta(a, b)` prior, specified with
#'      `bin_prior`, is updated with the number of events and non-events in
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
#'   With return_trace = TRUE, a goldilocks_trial object is returned. Its
#'   summary element is the same data frame and its trace element has one row
#'   per interim look. The trace records enrollment and observed events by arm,
#'   predictive probabilities and their thresholds, the decision taken, and
#'   warnings raised during that look. It deliberately excludes imputed data
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
#'  prior = c(0.1, 0.1),
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
  prior = c(0.1, 0.1),
  bin_prior = c(1, 1),
  bin_method = "mc",
  block = 2,
  rand_ratio = c(1, 1),
  prop_loss = 0,
  alternative = "greater",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.9,
  prob_ha = 0.95,
  N_impute = 10,
  N_mcmc = 10,
  empty_interval = c("propagate", "prior", "error"),
  method = "logrank",
  imputed_final = FALSE,
  return_trace = FALSE,
  binary_imputation = c("event-time", "bernoulli")
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
  validate_gamma_prior(prior)
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
    validate_bayes_binomial_args(bin_prior, bin_method, N_mcmc)
  }

  # Assign: if no interim looks, set thresholds to 0, as they are not needed
  if (N_looks == 1) {
    Sn <- 0
    Fn <- 0
    check_futility <- FALSE
  } else {
    validate_probability_vector(Sn, "Sn")
    if (!is.null(Fn)) {
      validate_probability_vector(Fn, "Fn")
    }

    N_interims <- N_looks - 1

    if (length(Sn) == 1) {
      Sn <- rep(Sn, N_interims)
    } else if (length(Sn) != N_interims) {
      stop("More thresholds specified than actual interim looks")
    }

    if (is.null(Fn)) {
      Fn <- rep(0, N_interims)
    } else if (length(Fn) == 1) {
      Fn <- rep(Fn, N_interims)
    } else if (length(Fn) != N_interims) {
      stop("More thresholds specified than actual interim looks")
    }

    check_futility <- any(Fn != 0)
  }

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
      capture_warning <- function(warning) {
        if (return_trace) {
          warning_state$messages <- unique(c(
            warning_state$messages,
            conditionMessage(warning)
          ))
        }
      }

      # Posterior distribution of lambdas: current data
      post_lambda <- withCallingHandlers(
        posterior(
          data = data,
          cutpoints = cutpoints,
          prior = prior,
          N_mcmc = N_impute,
          single_arm = single_arm,
          empty_interval = empty_interval
        ),
        warning = capture_warning
      )

      ##########################################################################
      ### Loop over multiple imputations
      ##########################################################################

      futility_test <- 0
      expected_success_test <- 0
      for (j in 1:N_impute) {
        h <- post_lambda[j, , , drop = FALSE]

        stop_check <- withCallingHandlers(
          test_stop_success(
            data = data_interim,
            hazard = h,
            end_of_study = end_of_study,
            cutpoints = cutpoints,
            single_arm = single_arm,
            prior = prior,
            N_mcmc = N_mcmc,
            method = method,
            alternative = alternative,
            h0 = h0,
            bin_prior = bin_prior,
            bin_method = bin_method,
            binary_imputation = if (method %in% c("bayes-bin", "riskdiff")) {
              binary_imputation
            } else {
              "event-time"
            },
            empty_interval = empty_interval,
            check_futility = check_futility
          ),
          warning = capture_warning
        )

        # Increment counter if P(efficacy | data) > prob_ha
        prob_now <- stop_check$success_now$success
        if (prob_now > prob_ha) {
          expected_success_test <- expected_success_test + 1
        }

        if (check_futility) {
          # Increase futility counter by 1 if P(efficacy | data) > prob_ha
          prob_max <- stop_check$success_max$success
          if (prob_max > prob_ha) {
            futility_test <- futility_test + 1
          }
        }
      }

      # Test if expected success criteria met
      # Note: ppp_success = posterior predictive probability of eventual success
      ppp_success <- expected_success_test / N_impute
      ppp_success_at_max <- if (check_futility) {
        futility_test / N_impute
      } else {
        NA_real_
      }

      decision <- if (ppp_success > Sn[i]) {
        "stop_expected_success"
      } else if (check_futility && ppp_success_at_max < Fn[i]) {
        "stop_futility"
      } else {
        "continue"
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
          success_threshold = Sn[i],
          ppp_success_at_max = ppp_success_at_max,
          futility_threshold = if (check_futility) Fn[i] else NA_real_,
          decision = decision,
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
    prior = prior,
    N_mcmc = N_mcmc,
    single_arm = single_arm,
    imputed_final = imputed_final,
    method = method,
    N_impute = N_impute,
    alternative = alternative,
    h0 = h0,
    bin_prior = bin_prior,
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

  if (return_trace) {
    out <- list(
      summary = results,
      trace = new_trial_trace(trace_rows),
      call = Call
    )
    class(out) <- "goldilocks_trial"
    return(out)
  }

  return(results)
}
