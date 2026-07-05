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
#' @param alternative character. The string specifying the alternative
#'   hypothesis, must be one of `"greater"` (default), `"less"` or
#'   `"two.sided"`. All three options are supported for `method = "bayes"`,
#'   `"logrank"`, and `"cox"`. The chi-square test (`method = "chisq"`) only
#'   supports `"two.sided"`. For survival outcomes, `"less"` corresponds to the
#'   treatment arm having a lower cumulative incidence (i.e., treatment is
#'   beneficial), and `"greater"`
#'   corresponds to the treatment arm having a higher cumulative incidence.
#' @param h0 scalar. Null hypothesis value or margin. Default is `h0 = 0`.
#'   * When `method = "bayes"`, `h0` is the null value of
#'     \eqn{p_\textrm{treatment} - p_\textrm{control}}. In a single-arm design,
#'     `h0` is the external benchmark event probability, often referred to
#'     as a performance goal (PG) or objective performance criterion (OPC).
#'   * When `method = "cox"`, `h0` is the null log hazard ratio for
#'     treatment versus control. Use `h0 = 0` for the usual hazard ratio
#'     of 1 null, or `h0 = log(margin)` for a non-inferiority margin
#'     specified as a hazard ratio. A Cox non-inferiority test should usually
#'     use `alternative = "less"`.
#'   * The argument is ignored for `method = "logrank"` and
#'     `method = "chisq"`; in those cases the usual method-specific null
#'     hypothesis is used.
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
#'   missing data.
#' @param N_mcmc integer. Number of samples to draw from the posterior
#'   distribution when using a Bayesian test (`method = "bayes"`).
#' @param method character. For an imputed data set (or the final data set after
#'   follow-up is complete), whether the analysis should be a log-rank
#'   (`method = "logrank"`) test, Cox proportional hazards regression model
#'   Wald test (`method = "cox"`), a fully-Bayesian analysis (`method
#'   = "bayes"`), or a chi-square test (`method = "chisq"`). See Details
#'   section.
#' @param imputed_final logical. Should the final analysis (after all subjects
#'   have been followed-up to the study end) be based on imputed outcomes for
#'   subjects who were LTFU (i.e. right-censored with time less than
#'   `end_of_study`)? Default is `TRUE`. Setting to `FALSE`
#'   means that the final analysis would incorporate right-censoring.
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
#'      ratio) is also reported.
#'
#'   * Bayesian absolute difference (`method = "bayes"`).
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
#'  * Chi-square test (`method = "chisq"`).
#'      Each (imputed) dataset with both treatment and control arms can be
#'      compared using a standard chi-square test on the final event status,
#'      which discards the event time information. The output is a *P*-value,
#'      and there is no treatment effect reported. The function returns \eqn{1 -
#'      P}, which is reported in `post_prob_ha`. Whilst
#'      not a posterior probability, it can be contrasted in the same manner.
#'      For example, if the success threshold is \eqn{P < 0.05}, then one
#'      requires `post_prob_ha` \eqn{> 0.95}. The reason for this is to
#'      enable simple switching between Bayesian and frequentist paradigms for
#'      analysis. Because the chi-square test cannot handle right-censored
#'      observations, subjects lost to follow-up are excluded from the final
#'      analysis when `imputed_final = FALSE`. When
#'      `imputed_final = TRUE`, LTFU subjects are imputed before the
#'      test is applied, so all subjects are included.
#'
#'  * Imputed final analysis (`imputed_final`).
#'      The overall final analysis conducted after accrual is suspended and
#'      follow-up is complete can be analyzed on imputed datasets (default) or
#'      on the non-imputed dataset. Since the imputations/predictions used
#'      during the interim analyses assume all subjects are imputed (since loss
#'      to follow-up is not yet known), it would seem most appropriate to
#'      conduct the trial in the same manner, especially if loss to follow-up
#'      rates are appreciable. Note, this only applies to subjects who are
#'      right-censored due to loss to follow-up, which we assume is a
#'      non-informative process. This can be used with any `method`.
#'
#'   When `method = "bayes"` and imputation is involved (either at interim
#'   analyses or via `imputed_final = TRUE`), a two-stage posterior
#'   procedure is used. First, the posterior distribution of the piecewise
#'   hazard rates is estimated from the *observed* data and used to draw
#'   imputed event times for censored subjects. Second, a *new* posterior
#'   is estimated from the combined observed and imputed data, and this
#'   posterior is used for inference. This is consistent with the predictive
#'   probability framework described in Broglio et al. (2014), but users
#'   should be aware that the imputation model's posterior influences the
#'   analysis posterior. For frequentist methods (`"logrank"`, `"cox"`,
#'   `"chisq"`), the second stage uses a standard test
#'   rather than a posterior, so this feedback loop does not arise.
#'
#'   At each interim look, follow-up times are masked (censored) to reflect
#'   the calendar time of the analysis. The package treats enrollment and
#'   randomization as occurring at the same time. Subjects enrolled at the exact
#'   interim boundary have zero follow-up time. These times are clamped to
#'   `.Machine$double.eps` (approximately \eqn{2.2 \times 10^{-16}}) so
#'   that they contribute negligible but non-zero exposure to the interim
#'   posterior. This affects at most one subject per interim look.
#'
#' @return A data frame containing some input parameters (arguments) as well as
#'   statistics from the analysis, including:
#'
#' | Column | Description |
#' | --- | --- |
#' | `N_treatment` | Number of patients enrolled in the treatment arm. |
#' | `N_control` | Number of patients enrolled in the control arm. |
#' | `est_final` | Treatment effect estimated at the final analysis. The final analysis occurs when either the maximum sample size is reached and follow-up is complete, or the interim analysis triggered early stopping of enrollment/accrual and follow-up for those subjects is complete. |
#' | `post_prob_ha` | Posterior probability from the final analysis. If `imputed_final` is `TRUE`, this is calculated for each imputed final-analysis dataset and averaged over `N_impute` imputations. If `method = "logrank"`, the value represents \eqn{1 - P}, where \eqn{P} is the frequentist *P*-value. |
#' | `stop_futility` | Logical indicator of whether the trial stopped early for futility. |
#' | `stop_expected_success` | Logical indicator of whether the trial stopped early for expected success. |
#'
#' @references
#' Broglio KR, Connor JT, Berry SM. Not too big, not too small: a Goldilocks
#' approach to sample size selection. *Journal of Biopharmaceutical Statistics*,
#' 2014; 24(3): 685–705.
#'
#' @importFrom stats pexp coef chisq.test
#' @export
#'
#' @examples
#' # RCT with exponential hazard (no piecewise breaks)
#' # Note: the number of imputations is small to enable this example to run
#' #       quickly on CRAN tests. In practice, much larger values are needed.
#' survival_adapt(
#'  hazard_treatment = -log(0.85) / 36,
#'  hazard_control = -log(0.7) / 36,
#'  cutpoints = 0,
#'  N_total = 600,
#'  lambda = 20,
#'  lambda_time = 0,
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
#'  method = "bayes")
survival_adapt <- function(
  hazard_treatment,
  hazard_control    = NULL,
  cutpoints         = 0,
  N_total,
  lambda            = 0.3,
  lambda_time       = 0,
  interim_look      = NULL,
  end_of_study,
  prior             = c(0.1, 0.1),
  block             = 2,
  rand_ratio        = c(1, 1),
  prop_loss         = 0,
  alternative       = "greater",
  h0                = 0,
  Fn                = 0.05,
  Sn                = 0.9,
  prob_ha           = 0.95,
  N_impute          = 10,
  N_mcmc            = 10,
  method            = "logrank",
  imputed_final     = FALSE) {

  ##############################################################################
  ### Derive variables
  ##############################################################################

  # Indicator of whether single-arm study
  single_arm <- is.null(hazard_control)

  # Interim look and final look
  analysis_at_enrollnumber <- c(interim_look, N_total)

  # Futility assessment required?
  if (is.null(Fn) | all(Fn == 0)) {
    check_futility <- FALSE
  } else {
    check_futility <- TRUE
  }

  # Number of looks
  N_looks <- length(analysis_at_enrollnumber)

  # If not using Bayesian test, then set N_mcmc = 1
  if (method != "bayes") {
    N_mcmc <- 1
  }

  ##############################################################################
  ### Run checks on arguments
  ##############################################################################

  # Check: 'interim_look' bounded by maximum sample size
  if (!is.null(interim_look)) {
    stopifnot(all(N_total > interim_look))

    # Check: each interim look is large enough to (with block randomization)
    # guarantee both treatment groups are represented. A look smaller than one
    # full block can enroll subjects from one treatment group only, which would
    # make the interim posterior undefined for the missing group. For two-arm
    # designs we require each look to be at least the (largest) block size.
    if (!single_arm) {
      min_look <- max(block)
      if (any(interim_look < min_look)) {
        stop(
          "Each 'interim_look' must be at least the block size (",
          min_look, ") so that both treatment groups are present at every ",
          "interim analysis. Smallest 'interim_look' given: ",
          min(interim_look), "."
        )
      }
    }
  }

  # Check: 'alternative' is correctly specified
  if (!alternative %in% c("two.sided", "greater", "less")) {
    stop("The input for alternative is wrong")
  }

  # Check: Bayesian test only available as a one-sided test
  if (alternative == "two.sided" & method == "bayes") {
    stop("The Bayes test can only be used with alternative equal to 'greater' or 'less'")
  }

  # Check: chi-square test only available as two-sided
  if (alternative != "two.sided" & method == "chisq") {
    stop("The chi-square test can only be applied as a two-sided test")
  }

  # Check: frequentist tests only available for two-armed trials
  if (single_arm & method %in% c("logrank", "cox", "chisq")) {
    stop("The selected method can only be used for two-armed trials")
  }

  # Check: thresholds of consistent dimension
  if (length(Sn) != length(Fn)) {
    stop("Probability threshold vectors need to be the same length")
  }

  # Check: thresholds available for each interim look
  if (N_looks <= length(Fn)) {
    stop("More thresholds specified than actual interim looks")
  }
  if (N_looks > 1 & length(Fn) == 1) {
    # Recycle thresholds if needed
    # Note - we already enforce that Sn = Fn
    Sn <- rep(Sn, N_looks - 1)
    Fn <- rep(Fn, N_looks - 1)
  }

  # Assign: if no interim looks, set thresholds to 0, as they are not needed
  if (is.null(Fn) & N_looks == 1) {
    Sn <- 0
    Fn <- 0
    check_futility <- FALSE
  }

  ##############################################################################
  ### Simulate complete data for trial
  ##############################################################################

  data_total <- sim_comp_data(
    hazard_treatment = hazard_treatment,
    hazard_control   = hazard_control,
    cutpoints        = cutpoints,
    N_total          = N_total,
    lambda           = lambda,
    lambda_time      = lambda_time,
    end_of_study     = end_of_study,
    block            = block,
    rand_ratio       = rand_ratio,
    prop_loss        = prop_loss
  )

  ##############################################################################
  ### Evaluate trial at each interim analysis
  ##############################################################################

  # Assigning stop_futility and stop_expected_success
  stop_futility         <- 0
  stop_expected_success <- 0

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
        subject_enrolled = (id <= analysis_at_enrollnumber[i])
        subject_impute_futility = !subject_enrolled
        time_from_rand_at_look = enrollment[analysis_at_enrollnumber[i]] - enrollment
        subject_impute_success =
          # Had event, but has not occurred yet (based on interim look)
          ((event == 1) * (time_from_rand_at_look < time) & subject_enrolled) |
          # Event-free and not had opportunity to complete full follow
          ((event == 0) * (time_from_rand_at_look < end_of_study) & subject_enrolled) |
          (loss_to_fu & subject_enrolled)
      })

      # Mask the data at time of look
      # Note: subjects at the exact interim boundary have
      # time_from_rand_at_look = 0, yielding time = 0 after masking.
      # Clamp to .Machine$double.eps so the boundary subject contributes
      # negligible but non-zero exposure to the interim posterior.
      data_interim <- within(data_interim, {
        time  = pmax(pmin(time, time_from_rand_at_look), .Machine$double.eps)
        event = ifelse(subject_impute_success, 0, event)
      })

      # Carry out interim analysis on patients with complete data only
      # - Set-up new 'data' data frame
      data <- subset(data_interim,
                     subset = subject_enrolled,
                     select = c(time, event, treatment))

      # Posterior distribution of lambdas: current data
      post_lambda <- posterior(data       = data,
                               cutpoints  = cutpoints,
                               prior      = prior,
                               N_mcmc     = N_impute,
                               single_arm = single_arm)

      ##########################################################################
      ### Loop over multiple imputations
      ##########################################################################

      futility_test         <- 0
      expected_success_test <- 0
      for (j in 1:N_impute) {

        h <- post_lambda[j, , , drop = FALSE]

        stop_check <- test_stop_success(
          data           = data_interim,
          hazard         = h,
          end_of_study   = end_of_study,
          cutpoints      = cutpoints,
          single_arm     = single_arm,
          prior          = prior,
          N_mcmc         = N_mcmc,
          method         = method,
          alternative    = alternative,
          h0             = h0,
          check_futility = check_futility
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
      if (ppp_success > Sn[i]) {
        stop_expected_success <- 1
        stage_trial_stopped   <- analysis_at_enrollnumber[i]
        break # No further SS looks
      }

      # Test if futility success criteria is met
      if (futility_test / N_impute < Fn[i]) {
        stop_futility       <- 1
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
    N_enrolled            <- N_total
    stage_trial_stopped   <- N_total
    stop_futility         <- 0
    stop_expected_success <- 0
    ppp_success           <- NA
  }

  ##############################################################################
  ### Final analysis (after enrollment complete)
  ##############################################################################

  # All patients that have made it to the end of study
  # - complete follow-up (except any censoring)
  data_final <- subset(data_total, id <= stage_trial_stopped)
  data_final <- within(data_final, {
    time_from_rand_at_look = enrollment[stage_trial_stopped] - enrollment
    subject_impute_success = ((event == 0) & (time < end_of_study))
  })

  results_final <- test_final(data_in       = data_final,
                              cutpoints     = cutpoints,
                              prior         = prior,
                              N_mcmc        = N_mcmc,
                              single_arm    = single_arm,
                              imputed_final = imputed_final,
                              method        = method,
                              N_impute      = N_impute,
                              alternative   = alternative,
                              h0            = h0,
                              end_of_study  = end_of_study)

  post_paa  <- results_final[1]
  est_final <- results_final[2]

  N_treatment  <- sum(data_final$treatment == 1) # Total analyzed: treatment
  N_control    <- sum(data_final$treatment == 0) # Total analyzed: control

  ##############################################################################
  ### Output
  ##############################################################################

    results <- data.frame(
    prob_threshold        = prob_ha,
    margin                = h0,                       # Margin
    alternative           = alternative,              # Alternative hypothesis
    N_treatment           = N_treatment,
    N_control             = N_control,
    N_enrolled            = N_treatment + N_control,
    N_max                 = N_total, 				          # Total potential sample size
    post_prob_ha          = post_paa,                 # Posterior probability that alternative hypothesis is true
    est_final             = est_final,                # Posterior treatment effect at final analysis
    ppp_success           = ppp_success,              # Posterior predictive probability of eventual success when trial stopped
    stop_futility         = stop_futility,            # Did the trial stop for futility
    stop_expected_success = stop_expected_success     # Did the trial stop for expected success
  )

  return(results)

}
