#' @title Simulate a single adaptive clinical trial with a time-to-event
#'   endpoint
#'
#' @param hazard_treatment vector. Constant hazard rates under the treatment arm.
#' @param hazard_control vector. Constant hazard rates under the control arm.
#' @param cutpoints  vector. Times at which the baseline hazard changes. Default
#'   is \code{cutpoints = 0}, which corresponds to a simple (non-piecewise)
#'   exponential model.
#' @param N_total integer. Maximum sample size allowable
#' @param lambda vector. Enrollment rates across simulated enrollment times. See
#'   \code{\link{enrollment}} for more details.
#' @param lambda_time vector. Enrollment time(s) at which the enrollment rates
#'   change. Must be same length as lambda. See
#'   \code{\link{enrollment}} for more details.
#' @param interim_look vector. Sample size for each interim look. Note: the
#'   maximum sample size should not be included.
#' @param end_of_study scalar. Length of the study; i.e. time at which endpoint
#'   will be evaluated.
#' @param prior vector. The prior distributions for the piecewise hazard rate
#'   parameters are each \eqn{Gamma(a_0, b_0)}, with specified (known)
#'   hyper-parameters \eqn{a_0} and \eqn{b_0}. The default non-informative prior
#'   distribution used is Gamma(0.1, 0.1), which is specified by setting
#'   \code{prior = c(0.1, 0.1)}.
#' @param block scalar. Block size for generating the randomization schedule.
#' @param rand_ratio vector. Randomization allocation for the ratio of control
#'   to treatment. Integer values mapping the size of the block. See
#'   \code{\link{randomization}} for more details.
#' @param prop_loss scalar. Overall proportion of subjects lost to
#'   follow-up. Defaults to zero.
#' @param alternative character. The string specifying the alternative
#'   hypothesis, must be one of \code{"greater"} (default), \code{"less"} or
#'   \code{"two.sided"}.
#' @param h0 scalar. Null hypothesis value of \eqn{p_\textrm{treatment} -
#'   p_\textrm{control}} when \code{method = "bayes"}. Default is \code{h0 = 0}.
#'   The argument is ignored when \code{method = "logrank"} or \code{= "cox"};
#'   in those cases the usual test of non-equal hazards is assumed.
#' @param Fn vector of \code{[0, 1]} values. Each element is the probability
#'   threshold to stop at the \eqn{i}-th look early for futility. If there are
#'   no interim looks (i.e. \code{interim_look = NULL}), then \code{Fn} is not
#'   used in the simulations or analysis. The length of \code{Fn} should be the
#'   same as \code{interim_look}, else the values are recycled.
#' @param Sn vector of \code{[0, 1]} values. Each element is the probability
#'   threshold to stop at the \eqn{i}-th look early for expected success. If
#'   there are no interim looks (i.e. \code{interim_look = NULL}), then
#'   \code{Sn} is not used in the simulations or analysis. The length of
#'   \code{Sn} should be the same as \code{interim_look}, else the values are
#'   recycled.
#' @param prob_ha scalar \code{[0, 1]}. Probability threshold of alternative
#'   hypothesis.
#' @param N_impute integer. Number of imputations for Monte Carlo simulation of
#'   missing data.
#' @param N_mcmc integer. Number of samples to draw from the posterior
#'   distribution when using a Bayesian test (\code{method = "bayes"}).
#' @param method character. For an imputed data set (or the final data set after
#'   follow-up is complete), whether the analysis should be a log-rank
#'   (\code{method = "logrank"}) test, Cox proportional hazards regression model
#'   Wald test (\code{method = "cox"}), or a fully-Bayesian analysis
#'   (\code{method = "bayes"}). See Details section.
#' @param imputed_final logical. Should the final analysis (after all subjects
#'   have been followed-up to the study end) be based on imputed outcomes for
#'   subjects who were LTFU (i.e. right-censored with time
#'   \code{<end_of_study})? Default is \code{TRUE}. Setting to \code{FALSE}
#'   means that the final analysis would incorporate right-censoring.
#' @param debug logical. If \code{TRUE} can be used to debug aspects of the
#'   code, including producing Kaplan-Meier graphs at each step of the
#'   algorithm. Default is \code{debug = FALSE}.
#'
#' @details Implements the Goldilocks design method described in Broglio et al.
#'   (2014). At each interim analysis, two probabilities are computed:
#'
#'   1. **The posterior predictive probability of eventual success.** This is
#'      calculated as the proportion of imputed datasets at the *current* sample
#'      size that would go on to be success at the specified threshold. At each
#'      interim analysis it is compared to the corresponding element of
#'      \code{Sn}, and if it exceeds the threshold,
#'      accrual/enrollment is suspended and the outstanding follow-up allowed to
#'      complete before conducting the pre-specified final analysis.
#'
#'   2. **The posterior predictive probability of final success**. This is
#'      calculated as the proportion of imputed datasets at the *maximum*
#'      threshold that would go on to be successful. Similar to above, it is
#'      compared to the corresponding element of \code{Fn}, and if it
#'      is less than the threshold, accrual/enrollment is suspended and the
#'      trial terminated. Typically this would be a binding decision. If it is
#'      not a binding decision, then one should also explore the simulations
#'      with \code{Fn = 0}.
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
#'  * Log-rank test (\code{method = "logrank"}).
#'      Each (imputed) dataset with both treatment and control arms can be
#'      compared using a standard log-rank test. The output is a \emph{P}-value,
#'      and there is no treatment effect reported. The function returns \eqn{1 -
#'      P}, which is reported in \code{post_prob_ha}. Whilst not a posterior
#'      probability, it can be contrasted in the same manner. For example, if
#'      the success threshold is \eqn{P < 0.05}, then one requires
#'      \code{post_prob_ha} \eqn{> 0.95}. The reason for this is to enable
#'      simple switching between Bayesian and frequentist paradigms for
#'      analysis.
#'
#'   * Cox proportional hazards regression Wald test (\code{method = "cox"}).
#'      Similar to the log-rank test, a \emph{P}-value is calculated based on a
#'      two-sided test. However, for consistency, \eqn{1 - P}, which is
#'      reported in \code{post_prob_ha}. Whilst not a posterior probability, it
#'      can be contrasted in the same manner. For example, if the success
#'      threshold is \eqn{P < 0.05}, then one requires \code{post_prob_ha}
#'      \eqn{> 0.95}.
#'
#'   * Bayesian absolute difference (\code{method = "bayes"}).
#'      Each imputed dataset is used to update the conjugate Gamma prior
#'      (defined by \code{prior}), yielding a posterior distribution for the
#'      piecewise exponential rate parameters. In turn, the posterior
#'      distribution of the cumulative incidence function (\eqn{1 - S(t)}, where
#'      \eqn{S(t)} is the survival function) evaluated at time
#'      \code{end_of_study} is calculated. If a single arm study, then this
#'      summarizes the treatment effect, else, if a two-armed study, the
#'      independent posteriors are used to estimate the posterior distribution
#'      of the difference. A posterior probability is calculated according to
#'      the specification of the test type (\code{alternative}) and the value of
#'      the null hypothesis (\code{h0}).
#'
#'  * Imputed final analysis (\code{imputed_final}).
#'      The overall final analysis conducted after accrual is suspended and
#'      follow-up is complete can be analyzed on imputed datasets (default) or
#'      on the non-imputed dataset. Since the imputations/predictions used
#'      during the interim analyses assume all subjects are imputed (since loss
#'      to follow-up is not yet known), it would seem most appropriate to
#'      conduct the trial in the same manner, especially if loss to follow-up
#'      rates are appreciable. Note, this only applies to subjects who are
#'      right-censored due to loss to follow-up, which we assume is a
#'      non-informative process. This can be used with any \code{method}.
#'
#' @return A data frame containing some input parameters (arguments) as well as
#'   statistics from the analysis, including:
#'
#'   \describe{
#'     \item{\code{N_treatment}}{
#'       integer. The number of patients enrolled in the treatment arm for
#'       each simulation.}
#'     \item{\code{N_control}}{
#'       integer. The number of patients enrolled in the control arm for
#'       each simulation.}
#'     \item{\code{est_interim}}{
#'       scalar. The treatment effect that was estimated at the time of the
#'       interim analysis. Note this is not actually used in the final
#'       analysis.}
#'     \item{\code{est_final}}{
#'       scalar. The treatment effect that was estimated at the final analysis.
#'       Final analysis occurs when either the maximum sample size is reached
#'       and follow-up complete, or the interim analysis triggered an early
#'       stopping of enrollment/accrual and follow-up for those subjects is
#'       complete.}
#'     \item{\code{post_prob_ha}}{
#'       scalar. The corresponding posterior probability from the final
#'       analysis. If \code{imputed_final} is true, this is calculated as the
#'       posterior probability of efficacy (or equivalent, depending on how
#'       \code{alternative} and \code{h0} were specified) for each imputed
#'       final analysis dataset, and then averaged over the \code{N_impute}
#'       imputations. If \code{method = "logrank"}, \code{post_prob_ha} is
#'       calculated in the same fashion, but value represents \eqn{1 - P},
#'       where \eqn{P} denotes the frequentist \eqn{P}-value.}
#'     \item{\code{stop_futility}}{
#'       integer. A logical indicator of whether the trial was stopped early for
#'       futility.}
#'     \item{\code{stop_expected_success}}{
#'       integer. A logical indicator of whether the trial was stopped early for
#'       expected success.}
#'  }
#'
#' @references
#' Broglio KR, Connor JT, Berry SM. Not too big, not too small: a Goldilocks
#' approach to sample size selection. \emph{Journal of Biopharmaceutical
#' Statistics}, 2014; 24(3): 685–705.
#'
#' @importFrom stats pexp runif sd coef
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
  imputed_final     = FALSE,
  debug             = FALSE
) {


  # Check: 'interim_look' bounded by maximum sample size
  if (!is.null(interim_look)) {
    stopifnot(all(N_total > interim_look))
  }

  # Check: 'alternative' is correctly specified
  if (alternative != "two.sided" & alternative != "greater" & alternative != "less") {
    stop("The input for alternative is wrong")
  }

  # Check: none of the 'cutpoints' are not more than 'end_of_study'
  if (!is.null(cutpoints)) {
    stopifnot(any(cutpoints < end_of_study))
  }

  # Check: Bayesian test only available as a one-sided test
  if (alternative == "two.sided" & method == "bayes") {
    stop("The Bayes test can only be used with alternative equal to 'greater' or 'less'")
  }

  # Check: log-rank test only available as a two.sided test
  if (alternative != "two.sided" & method %in% c("logrank", "cox")) {
    stop("The log-rank and Cox PH tests can only be applied with a two-sided test")
  }

  # Assign: indicator of whether single-arm study
  single_arm <- is.null(hazard_control)

  # Check: log-rank test only available for two-armed trials
  if (single_arm & method %in% c("logrank", "cox")) {
    stop("The log-rank and Cox PH tests can only be used for two-armed trials")
  }

  # Assign: interim look and final look
  analysis_at_enrollnumber <- c(interim_look, N_total)

  # Assign: futility assessment required?
  if (is.null(Fn) | all(Fn == 0)) {
    check_futility <- FALSE
  } else {
    check_futility <- TRUE
  }

  # Check: thresholds of consistent dimension
  if (length(Sn) != length(Fn)) {
    stop("Probability threshold vectors need to be the same length")
  }

  # Check: thresholds available for each interim look
  N_looks <- length(analysis_at_enrollnumber)
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

  # Assign: if not using Bayesian test, then set N_mcmc = 1
  if (method != "bayes") {
    N_mcmc <- 1
  }

  # Simulate enrollment times
  enrollment <- enrollment(lambda      = lambda,
                           N_total     = N_total,
                           lambda_time = lambda_time)
  enrollment <- enrollment + runif(length(enrollment))
  enrollment <- sort(enrollment)

  # Simulate treatment arm assignment
  if (!single_arm) {
    group <- randomization(N_total    = N_total,
                           block      = block,
                           allocation = rand_ratio)
  } else {
    group <- rep(1, N_total)
  }

  time  <- rep(NA, length = N_total)
  event <- rep(NA, length = N_total)

  # Simulate TTE outcome
  # - Note: time = time *from* enrollment
  if (!single_arm) {
    sim_control <- pwe_sim(hazard     = hazard_control,
                           n          = sum(!group),
                           maxtime    = end_of_study,
                           cutpoints  = cutpoints)
    time[group == 0]  <- sim_control$time
    event[group == 0] <- sim_control$event
  }

  sim_treatment <- pwe_sim(hazard    = hazard_treatment,
                           n         = sum(group),
                           maxtime   = end_of_study,
                           cutpoints = cutpoints)
  time[group == 1]  <- sim_treatment$time
  event[group == 1] <- sim_treatment$event

  # Simulate loss to follow-up
  loss_to_fu   <- rep(FALSE, N_total)
  if (prop_loss > 0) {
    n_loss_to_fu <- ceiling(prop_loss * N_total)
    loss_to_fu[sample(1:N_total, n_loss_to_fu)] <- TRUE
  }

  # Creating a new data.frame for all the variables
  data_total <- data.frame(
    time       = time,
    treatment  = group,
    event      = event,
    enrollment = enrollment,
    id         = 1:N_total,
    loss_to_fu = loss_to_fu)

  # Subjects lost are uniformly distributed
  if (prop_loss > 0) {
    data_total$time[data_total$loss_to_fu]  <- runif(
      n_loss_to_fu, 0, data_total$time[data_total$loss_to_fu])
    data_total$event[data_total$loss_to_fu] <- rep(0, n_loss_to_fu)
  }

  # KM plot for actual simulated data
  if (debug) {
    if (!single_arm) {
      plot(survfit(Surv(time, event) ~ treatment, data = data_total),
           col  = c(1, 2),
           main = "Simulated Data: Complete",
           xlab = "Time",
           ylab = "Freedom from event")
    } else if (single_arm) {
      plot(survfit(Surv(time, event) ~ 1, data = data_total),
           main = "Simulated Data: Complete",
           xlab = "Time",
           ylab = "Freedom from event")
    }
  }

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
      # - time_from_rand_at_look:  time from randomization to sample size look
      #                            e.g. if patient enrolled month 3, but look occurs month 7,
      #                            then patient could potentially be observed for 4 months

      if (debug) {
        print(paste("Look:", i))
      }

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
      data_interim <- within(data_interim, {
        time  = pmin(time, time_from_rand_at_look) + sd(time) / 1e4
        event = ifelse(subject_impute_success, 0, event)
      })

      # Carry out interim analysis on patients with complete data only
      # - Set-up new 'data' data frame
      data <- subset(data_interim,
                     subset = subject_enrolled,
                     select = c(time, event, treatment))

      # KM plot for masked data at interim analysis
      if (debug) {
        if (!single_arm) {
          plot(survfit(Surv(time, event) ~ treatment, data = data),
               col  = c(1, 2),
               main = "Simulated Data: Masked @ IA",
               xlab = "Time",
               ylab = "Freedom from event")
        } else if (single_arm) {
          plot(survfit(Surv(time, event) ~ 1, data = data),
               main = "Simulated Data: Masked @ IA",
               xlab = "Time",
               ylab = "Freedom from event")
        }
      }

      # Posterior distribution of lambdas: current data
      post_lambda <- posterior(data       = data,
                               cutpoints  = cutpoints,
                               prior      = prior,
                               N_mcmc     = N_impute,
                               single_arm = single_arm)

      # Imputation phase futility and expected success - initialize counters
      # for the current imputation phase
      futility_test         <- 0
      expected_success_test <- 0

      for (j in 1:N_impute) {

        if (debug) {
          print(paste("Imputation:", j))
        }

        ##########################################################################
        ### Expected success computations
        ##########################################################################

        # Single imputed data set
        data_success_impute <- impute_data(
          data_in      = data_interim,
          hazard       = post_lambda[j, , , drop = FALSE],
          end_of_study = end_of_study,
          cutpoints    = cutpoints,
          type         = "success",
          single_arm   = single_arm)

        # Create enrolled subject data frame for analysis
        data <- subset(data_success_impute,
                       subset = subject_enrolled,
                       select = c(time, event, treatment))

        # KM plot for imputed data at interim analysis (expected success)
        if (debug & (j == 1)) {
          if (!single_arm) {
            plot(survfit(Surv(time, event) ~ treatment, data = data),
                 col  = c(1, 2),
                 main = "Simulated Data: Imputed For Expected Success",
                 xlab = "Time",
                 ylab = "Freedom from event")
          } else if (single_arm) {
            plot(survfit(Surv(time, event) ~ 1, data = data),
                 main = "Simulated Data: Imputed For Expected Success",
                 xlab = "Time",
                 ylab = "Freedom from event")
          }
        }

        if (debug) {
          print("Testing: expected success")
        }

        # Apply primary analysis to imputed data
        success <- analyse_data(data         = data,
                                cutpoints    = cutpoints,
                                end_of_study = end_of_study,
                                prior        = prior,
                                N_mcmc       = N_mcmc,
                                single_arm   = single_arm,
                                method       = method,
                                alternative  = alternative,
                                h0           = h0)$success

        # Increase success counter by 1 if P(efficacy | data) > prob_ha
        if (success > prob_ha) {
          expected_success_test <- expected_success_test + 1
        }

        ##########################################################################
        ### Futility computations
        ##########################################################################

        if (check_futility) {
          # Take the already imputed data for expected success and append on
          # imputed event times for subjects not yet enrolled

          # Single imputed data set
          data_futility_impute <- impute_data(
            data_in      = data_success_impute,
            hazard       = post_lambda[j, , , drop = FALSE],
            end_of_study = end_of_study,
            cutpoints    = cutpoints,
            type         = "futility",
            single_arm   = single_arm)

          # Create data frame for analysis
          data <- subset(data_futility_impute,
                         select = c(time, event, treatment))

          # KM plot for imputed data at interim analysis (futility)
          if (debug & (j == 1)) {
            if (!single_arm) {
              plot(survfit(Surv(time, event) ~ treatment, data = data),
                   col  = c(1, 2),
                   main = "Simulated Data: Imputed For Futility",
                   xlab = "Time",
                   ylab = "Freedom from event")
            } else if (single_arm) {
              plot(survfit(Surv(time, event) ~ 1, data = data),
                   main = "Simulated Data: Imputed For Futility",
                   xlab = "Time",
                   ylab = "Freedom from event")
            }
          }

          if (debug) {
            print("Testing: futility")
          }

          # Apply primary analysis to imputed data
          success <- analyse_data(data         = data,
                                  cutpoints    = cutpoints,
                                  end_of_study = end_of_study,
                                  prior        = prior,
                                  N_mcmc       = N_mcmc,
                                  single_arm   = single_arm,
                                  method       = method,
                                  alternative  = alternative,
                                  h0           = h0)$success

          # Increase futility counter by 1 if P(efficacy | data) > prob_ha
          if (success > prob_ha) {
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

    ##############################################################################
    ### Effect at interim analysis (where trial stopped)
    ##############################################################################

    # TBC
    effect_int <- NA

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
  ### Effect at final analysis (after enrollment complete)
  ##############################################################################

  # All patients that have made it to the end of study
  # - complete follow-up (except any censoring)
  data_final <- subset(data_total, id <= stage_trial_stopped)
  data_final <- within(data_final, {
    time_from_rand_at_look = enrollment[analysis_at_enrollnumber[i]] - enrollment
    subject_impute_success = ((event == 0) & (time < end_of_study))
  })

  if (imputed_final) {
    # Posterior distribution of lambdas: final data
    post_lambda_final <- posterior(data       = data_final,
                                   cutpoints  = cutpoints,
                                   prior      = prior,
                                   N_mcmc     = N_impute,
                                   single_arm = single_arm)
    # Effect matrix + posterior probability
    effect_final_mat <- matrix(nrow = ifelse(method == "bayes", N_mcmc, 1),
                               ncol = N_impute)
    post_paa <- vector(length = N_impute)
    # Impute multiple data sets
    for (j in 1:N_impute) {
      # Single imputed data set
      data_success_impute <- impute_data(
        data_in      = data_final,
        hazard       = post_lambda_final[j, , , drop = FALSE],
        end_of_study = end_of_study,
        cutpoints    = cutpoints,
        type         = "success",
        single_arm   = single_arm)

      # Create enrolled subject data frame for analysis
      data <- subset(data_success_impute,
                     select = c(time, event, treatment))

      # Apply primary analysis to imputed data
      success <- analyse_data(data         = data,
                              cutpoints    = cutpoints,
                              end_of_study = end_of_study,
                              prior        = prior,
                              N_mcmc       = N_mcmc,
                              single_arm   = single_arm,
                              method       = method,
                              alternative  = alternative,
                              h0           = h0)

      post_paa[j] <- success$success
      if (method %in% c("cox", "bayes")) {
        effect_final_mat[, j] <- success$effect # See Gelman et al. (2004, p. 520)
      }
    }
    # Average over imputations
    post_paa  <- mean(post_paa)
    est_final <- mean(effect_final_mat)
  } else {
    # Apply primary analysis to final data (without imputation)
    success <- analyse_data(data         = data_final,
                            cutpoints    = cutpoints,
                            end_of_study = end_of_study,
                            prior        = prior,
                            N_mcmc       = N_mcmc,
                            single_arm   = single_arm,
                            method       = method,
                            alternative  = alternative,
                            h0           = h0)

    post_paa <- success$success
    est_final <- mean(success$effect)
  }

  N_treatment  <- sum(data_final$treatment == 1) # Total sample size analyzed: test group
  N_control    <- sum(data_final$treatment == 0) # Total sample size analyzed: control group

  if (length(analysis_at_enrollnumber) > 1) {
    est_interim <- mean(effect_int) # Posterior treatment effect at the interim analysis (if any)
  } else {
    est_interim <- NA
  }

  # Output
  results <- data.frame(
    prob_threshold        = prob_ha,
    margin                = h0,                       # Margin for error
    alternative           = alternative,              # Alternative hypothesis
    N_treatment           = N_treatment,
    N_control             = N_control,
    N_enrolled            = N_treatment + N_control,
    N_max                 = N_total, 				          # Total potential sample size
    post_prob_ha          = post_paa,                 # Posterior probability that alternative hypothesis is true
    est_final             = est_final,                # Posterior treatment effect at final analysis
    ppp_success           = ppp_success,              # Posterior predictive probability of eventual success when trial stopped
    est_interim           = est_interim,              # Posterior treatment effect at the interim analysis
    stop_futility         = stop_futility,            # Did the trial stop for futility
    stop_expected_success = stop_expected_success     # Did the trial stop for expected success
  )

  return(results)

}
