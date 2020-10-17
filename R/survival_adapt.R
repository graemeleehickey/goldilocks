#' @title Simulate a single adaptive clinical trial with a time-to-event
#'   endpoint
#'
#' @param hazard_treatment vector. Constant hazard rates under the treatment arm.
#' @param hazard_control vector. Constant hazard rates under the control arm.
#' @param cutpoint vector. Default is \code{cutpoint = 0}.
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
#' @param prop_loss_to_followup scalar. Overall proportion of subjects lost to
#'   follow-up.
#' @param alternative character. The string specifying the alternative
#'   hypothesis, must be one of \code{"greater"} (default), \code{"less"} or
#'   \code{"two-sided"}.
#' @param h0 scalar. Null hypothesis value of \eqn{p_\text{treatment} -
#'   p_\text{control}}. Default is \code{h0 = 0}.
#' @param futility_prob scalar \code{[0, 1]}. Probability threshold to stop
#'   early for futility.
#' @param expected_success_prob scalar \code{[0, 1]}. Probability threshold to
#'   stop early for expected success.
#' @param prob_ha scalar \code{[0, 1]}. Probability threshold of alternative
#'   hypothesis.
#' @param N_impute integer. Number of imputations for Monte Carlo simulation of
#'   missing data.
#' @param N_mcmc integer. Number of samples to from the posterior
#'   distribution.
#' @param imputed_final logical. Should the final analysis (after all subjects
#'   have been followed-up to the study end) be based on imputed outcomes for
#'   subjects who were LTFU (i.e. right-censored with time
#'   \code{<end_of_study})? Default is \code{TRUE}. Setting to \code{FALSE}
#'   means that the final analysis would incorporate right-censoring.
#' @param debug logical. If \code{TRUE} can be used to debug aspects of the
#'   code. Default is \code{debug = FALSE}.
#'
#' @return A list containing key input parameters (arguments) as well as
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
#'       analysis.}
#'     \item{\code{stop_futility}}{
#'     integer. A logical indicator of whether the
#'       trial was stopped early for futility.}
#'     \item{\code{stop_expected_success}}{
#'     integer. A logical indicator of whether the
#'       trial was stopped early for expected success.}
#'  }
#'
#' @importFrom stats pexp runif sd
#' @import dplyr
#' @export
#'
#' @examples
#' # RCT with exponential hazard (no piecewise breaks)
#' # Note: the number of imputations is small to enable this example to run
#' #       quickly on CRAN tests. In practice, much larger values are needed.
#' survival_adapt(
#'  hazard_treatment = -log(0.85) / 36,
#'  hazard_control = -log(0.7) / 36,
#'  cutpoint = 0,
#'  N_total = 600,
#'  lambda = 20,
#'  lambda_time = NULL,
#'  interim_look = 400,
#'  end_of_study = 36,
#'  prior = c(.1, .1),
#'  block = 2,
#'  rand_ratio = c(1, 1),
#'  prop_loss_to_followup = 0.30,
#'  alternative = "less",
#'  h0 = 0,
#'  futility_prob = 0.05,
#'  expected_success_prob = 0.9,
#'  prob_ha = 0.975,
#'  N_impute = 10,
#'  N_mcmc = 10)
survival_adapt <- function(
  hazard_treatment,
  hazard_control        = NULL,
  cutpoint              = 0,
  N_total,
  lambda                = 0.3,
  lambda_time           = NULL,
  interim_look          = NULL,
  end_of_study,
  prior                 = c(0.1, 0.1),
  block                 = 2,
  rand_ratio            = c(1, 1),
  prop_loss_to_followup = 0.10,
  alternative           = "greater",
  h0                    = 0,
  futility_prob         = 0.05,
  expected_success_prob = 0.9,
  prob_ha               = 0.95,
  N_impute              = 10,
  N_mcmc                = 100,
  imputed_final         = TRUE,
  debug = FALSE
) {

  # Check: 'interim_look' bounded by maximum sample size
  if (!is.null(interim_look)) {
    stopifnot(all(N_total > interim_look))
  }

  # Check: 'alternative' is correctly specified
  if (alternative != "two-sided" & alternative != "greater" & alternative != "less") {
    stop("The input for alternative is wrong!")
  }

  # Check: none of the 'cutpoints' are not more than 'end_of_study'
  if (!is.null(cutpoint)) {
    stopifnot(any(cutpoint < end_of_study))
  }

  # Assigning interim look and final look
  analysis_at_enrollnumber <- c(interim_look, N_total)

  # Assignment of enrollment based on the enrollment function
  enrollment <- enrollment(param   = lambda,
                           N_total = N_total,
                           time    = lambda_time)
  enrollment <- enrollment + runif(length(enrollment))
  enrollment <- sort(enrollment)

  # Indicator of single-arm study
  single_arm <- is.null(hazard_control)

  # Simulating group and treatment group assignment
  if (!single_arm) {
    group <- randomization(N_total = N_total, block = block,
                           allocation = rand_ratio)
  } else {
    group <- rep(1, N_total)
  }

  time  <- rep(NA, length = N_total)
  event <- rep(NA, length = N_total)

  # Simulate TTE outcome
  # - Note: time = time *from* enrollment
  if (!single_arm) {
    sim_control <- pwe_sim(hazard    = hazard_control,
                           n         = sum(!group),
                           maxtime   = end_of_study,
                           cutpoint  = cutpoint)
    time[which(!group)]  <- sim_control$time
    event[which(!group)] <- sim_control$event
  }

  sim_treatment <- pwe_sim(hazard   = hazard_treatment,
                           n        = sum(group),
                           maxtime  = end_of_study,
                           cutpoint = cutpoint)
  time[which(!!group)]  <- sim_treatment$time
  event[which(!!group)] <- sim_treatment$event

  # Simulate loss to follow-up
  n_loss_to_fu <- ceiling(prop_loss_to_followup * N_total)
  loss_to_fu   <- rep(FALSE, N_total)
  loss_to_fu[sample(1:N_total, n_loss_to_fu)] <- TRUE

  # Creating a new data.frame for all the variables
  data_total <- data.frame(
    time       = time,
    treatment  = group,
    event      = event,
    enrollment = enrollment,
    id         = 1:N_total,
    loss_to_fu = loss_to_fu)

  # Subjects lost are uniformly distributed
  data_total$time[data_total$loss_to_fu]  <- runif(n_loss_to_fu,
                                                   0, data_total$time[data_total$loss_to_fu])
  data_total$event[data_total$loss_to_fu] <- rep(0, n_loss_to_fu)

  # KM plot for actual simulated data
  if (debug & !single_arm) {
    plot(survfit(Surv(time, event) ~ treatment, data = data_total),
         col = c(1, 2),
         main = "Simulated Data: Complete",
         xlab = "Time",
         ylab = "Freedom from event")
  }

  # Assigning stop_futility and expected success
  stop_futility         <- 0
  stop_expected_success <- 0

  if (length(analysis_at_enrollnumber) > 1) {
    for (i in 1:(length(analysis_at_enrollnumber) - 1)) {

      # Analysis at the `analysis_at_enrollnumber` look
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
      data_interim <- data_total %>%
        mutate(
          subject_enrolled = id <= analysis_at_enrollnumber[i],
          subject_impute_futility = !subject_enrolled,
          time_from_rand_at_look = enrollment[analysis_at_enrollnumber[i]] - enrollment,
          subject_impute_success =
            # Had event, but has not occurred yet (based on interim look)
            ((event == 1) * (time_from_rand_at_look < time) & subject_enrolled) |
            # Event-free and not had opportunity to complete full follow
            ((event == 0) * (time_from_rand_at_look < end_of_study) & subject_enrolled) |
            (loss_to_fu & subject_enrolled))

      # Mask the data at time of look
      data_interim <- data_interim %>%
        mutate(
          time = pmin(time, time_from_rand_at_look) + sd(time) / 1e4,
          event = if_else(subject_impute_success, 0, event))

      # Carry out interim analysis on patients with complete data only
      # - Set-up new 'data' data frame
      data <- data_interim %>%
        filter(subject_enrolled) %>%
        select(time, event, treatment)

      # KM plot for masked data at interim analysis
      if (debug & !single_arm) {
        plot(survfit(Surv(time, event) ~ treatment, data = data),
             col = c(1, 2),
             main = "Simulated Data: Masked @ IA",
             xlab = "Time",
             ylab = "Freedom from event")
      }

      # Posterior distribution of lambdas: current data
      post_lambda <- posterior(data       = data,
                               cutpoint   = cutpoint,
                               prior      = prior,
                               N_mcmc     = N_mcmc,
                               single_arm = single_arm)

      # Imputation phase futility and expected success - initialize counters
      # for the current imputation phase
      futility_test         <- 0
      expected_success_test <- 0

      for (j in 1:N_impute) {

        ##########################################################################
        ### Expected success computations
        ##########################################################################

        data_success_impute <- impute_data(data_in = data_interim,
                                           hazard = post_lambda$post_treatment[j, ],
                                           end_of_study = end_of_study,
                                           cutpoint = cutpoint,
                                           type = "success",
                                           single_arm = single_arm)

        # Create enrolled subject data frame for analysis
        data <- data_success_impute %>%
          filter(subject_enrolled)  %>%
          ungroup() %>%
          select(time, event, treatment)

        # KM plot for imputed data at interim analysis (expected success)
        if (debug & !single_arm & (j == 1)) {
          plot(survfit(Surv(time, event) ~ treatment, data = data),
               col = c(1, 2),
               main = "Simulated Data: Imputed For Expected Success",
               xlab = "Time",
               ylab = "Freedom from event")
        }

        # Posterior distribution of lambdas: imputed data
        post_lambda_imp <- posterior(data       = data,
                                     cutpoint   = cutpoint,
                                     prior      = prior,
                                     N_mcmc     = N_mcmc,
                                     single_arm = single_arm)

        # Posterior distribution of event proportions: imputed data
        post_imp <- haz_to_prop(post         = post_lambda_imp,
                                cutpoint     = cutpoint,
                                end_of_study = end_of_study,
                                single_arm   = single_arm)

        # Apply statistical test to declare success (e.g. efficacy)
        success <- test(post_probs  = post_imp,
                        alternative = alternative,
                        h0          = h0,
                        single_arm  = single_arm)

        # Increase success counter by 1 if P(efficacy | data) > prob_ha
        if (success > prob_ha) {
          expected_success_test <- expected_success_test + 1
        }

        ##########################################################################
        ### Futility computations
        ##########################################################################

        # Take the already imputed data for expected success and append on
        # imputed event times for subjects not yet enrolled

        data_futility_impute <- impute_data(data_in = data_success_impute,
                                            hazard = post_lambda$post_treatment[j, ],
                                            end_of_study = end_of_study,
                                            cutpoint = cutpoint,
                                            type = "futility",
                                            single_arm = single_arm)

        # Create enrolled subject data frame for analysis
        data <- data_futility_impute %>%
          select(time, event, treatment)

        # KM plot for imputed data at interim analysis (futility)
        if (debug & !single_arm & (j == 1)) {
          plot(survfit(Surv(time, event) ~ treatment, data = data),
               col = c(1, 2),
               main = "Simulated Data: Imputed For Futility",
               xlab = "Time",
               ylab = "Freedom from event")
        }

        # Posterior distribution of lambdas: imputed data
        post_lambda_imp <- posterior(data       = data,
                                     cutpoint   = cutpoint,
                                     prior      = prior,
                                     N_mcmc     = N_mcmc,
                                     single_arm = single_arm)

        # Posterior distribution of event proportions: imputed data
        post_imp <- haz_to_prop(post         = post_lambda_imp,
                                cutpoint     = cutpoint,
                                end_of_study = end_of_study,
                                single_arm   = single_arm)

        # Apply statistical test to declare success (e.g. efficacy)
        success <- test(post_probs  = post_imp,
                        alternative = alternative,
                        h0          = h0,
                        single_arm  = single_arm)

        # Increase futility counter by 1 if P(efficacy | data) > prob_ha
        if (success > prob_ha) {
          futility_test <- futility_test + 1
        }

      }

      # Test if expected success criteria met
      if (expected_success_test / N_impute > expected_success_prob) {
        stop_expected_success <- 1
        stage_trial_stopped   <- analysis_at_enrollnumber[i]
        break # No further SS looks
      }

      # Test if futility success criteria is met
      if (futility_test / N_impute < futility_prob) {
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

    # Note 1: not imputed
    # Note 2: do not need to calculate posterior probability of success here,
    #         we are just interested in the effect size

    # Posterior distribution of event proportions: non-imputed data
    post_imp <- haz_to_prop(post         = post_lambda,
                            cutpoint     = cutpoint,
                            end_of_study = end_of_study,
                            single_arm   = single_arm)

    # Final interim analysis effect size
    effect_int <- mean(post_imp$effect)

    # Number of patients enrolled at trial stop
    N_enrolled <- nrow(data_interim[data_interim$id <= stage_trial_stopped, ])
  } else {
    # Assigning stage trial stopped given no interim look
    N_enrolled            <- N_total
    stage_trial_stopped   <- N_total
    stop_futility         <- 0
    stop_expected_success <- 0
  }

  ##############################################################################
  ### Effect at final analysis (after enrollment complete)
  ##############################################################################

  # All patients that have made it to the end of study
  # - complete follow-up (except any censoring)
  data_final <- data_total %>%
    filter(id <= stage_trial_stopped) %>%
    mutate(
      time_from_rand_at_look = enrollment[analysis_at_enrollnumber[i]] - enrollment,
      subject_impute_success =
        ((event == 0) & (time < end_of_study)))

  # Posterior distribution of lambdas: final data
  post_lambda_final <- posterior(data       = data_final,
                                 cutpoint   = cutpoint,
                                 prior      = prior,
                                 N_mcmc     = N_mcmc,
                                 single_arm = single_arm)

  if (imputed_final) {
    # Final analysis will multiply-impute subjects LTFU
    # Based on argument in Gelman et al. (2004, p. 520)
    effect_final_mat <- matrix(nrow = N_mcmc, ncol = N_impute)
    for (j in 1:N_impute) {
      # Step 1: Draw from post_lambda_final & impute a data set
      data_success_impute <- impute_data(
        data_in      = data_final,
        hazard       = post_lambda_final$post_treatment[j, ],
        end_of_study = end_of_study,
        cutpoint     = cutpoint,
        type         = "success",
        single_arm   = single_arm)
      data <- data_success_impute %>%
        select(time, event, treatment)
      # Step 2: Get posterior distribution of lambdas on imputed data
      post_lambda_final_imp <- posterior(data = data,
                                         cutpoint   = cutpoint,
                                         prior      = prior,
                                         N_mcmc     = N_mcmc,
                                         single_arm = single_arm)
      # Step 3: Transform to posterior sample of treatment effect
      post_final <- haz_to_prop(post         = post_lambda_final_imp,
                                cutpoint     = cutpoint,
                                end_of_study = end_of_study,
                                single_arm   = single_arm)
      # Step 4. Pool effects together
      effect_final_mat[, j] <- post_final$effect
    }
    # Step 5(i): Calculate average effect
    est_final <- mean(effect_final_mat)
    # Step 5(ii): Calculate average probability
    if (alternative == "two-sided") {
      post_paa <- max(c(mean(effect_final_mat > h0), mean(effect_final_mat < h0)))
    } else if (alternative == "greater") {
      post_paa <- mean(effect_final_mat > h0)
    } else {
      post_paa <- mean(effect_final_mat < h0)
    }
  } else {
    # Posterior distribution of event proportions: final data
    post_final <- haz_to_prop(post         = post_lambda_final,
                              cutpoint     = cutpoint,
                              end_of_study = end_of_study,
                              single_arm   = single_arm)

    # Apply statistical test to declare success (e.g. efficacy)
    post_paa <- test(post_final, alternative, h0, single_arm)

    # Final interim analysis effect size
    est_final <- mean(post_final$effect)
  }

  N_treatment  <- sum(data_final$treatment == 1) # Total sample size analyzed - test group
  N_control    <- sum(data_final$treatment == 0) # Total sample size analyzed - control group

  if (length(analysis_at_enrollnumber) > 1) {
    est_interim <- mean(effect_int) # Posterior treatment effect at the interim analysis (if any)
  } else {
    est_interim <- NA
  }

  # Output
  results_list <- list(
    hazard_treatment      = hazard_treatment,         # Probability of treatment in binomial
    hazard_control        = hazard_control,           # Probability of control in binomial
    cutpoint              = cutpoint,
    prob_threshold        = prob_ha,
    margin                = h0,                       # Margin for error
    alternative           = alternative,              # Alternative hypothesis
    interim_look          = interim_look,             # Print interim looks
    N_treatment           = N_treatment,
    N_control             = N_control,
    N_enrolled            = N_treatment + N_control,
    N_max                 = N_total, 				          # Total potential sample size
    post_prob_ha          = post_paa,                 # Posterior probability that alternative hypothesis is true
    est_final             = est_final,                # Posterior treatment effect at final analysis
    est_interim           = est_interim,              # Posterior treatment effect at the interim analysis
    stop_futility         = stop_futility,            # Did the trial stop for futility
    stop_expected_success = stop_expected_success     # Did the trial stop for expected success
  )

  # Return results
  return(results_list)

}
