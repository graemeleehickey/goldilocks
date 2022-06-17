#' @title Run test for whether we should stop enrolment at the interim analysis
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param check_futility Logical. Does the adaptive design include a test for
#'   assessment of futility?
#'
#' @return See analyse_data
#' @noRd
test_stop_success <- function(data, hazard, end_of_study, cutpoints, single_arm,
                              prior, N_mcmc, method, alternative, h0,
                              check_futility) {

  ##############################################################################
  ### Test for success at current sample size (-> stop for success)
  ##############################################################################

  # Single imputed data set
  data_success_impute <- impute_data(
    data_in      = data,
    hazard       = hazard,
    end_of_study = end_of_study,
    cutpoints    = cutpoints,
    type         = "success",
    single_arm   = single_arm)

  # Create enrolled subject data frame for analysis
  time             <- NULL
  event            <- NULL
  subject_enrolled <- NULL
  data <- subset(data_success_impute,
                 subset = subject_enrolled,
                 select = c(time, event, treatment))

  # Apply primary analysis to imputed data
  success_now <- analyse_data(data         = data,
                              cutpoints    = cutpoints,
                              end_of_study = end_of_study,
                              prior        = prior,
                              N_mcmc       = N_mcmc,
                              single_arm   = single_arm,
                              method       = method,
                              alternative  = alternative,
                              h0           = h0)

  ##############################################################################
  ### Test for success at maximum sample size (-> stop for futility)
  ##############################################################################

  if (check_futility) {

    # Take the already imputed data for expected success and append on
    # imputed event times for subjects not yet enrolled

    # Single imputed data set
    data_futility_impute <- impute_data(
      data_in      = data_success_impute,
      hazard       = hazard,
      end_of_study = end_of_study,
      cutpoints    = cutpoints,
      type         = "futility",
      single_arm   = single_arm)

    # Create data frame for analysis
    data <- subset(data_futility_impute,
                   select = c(time, event, treatment))

    # Apply primary analysis to imputed data
    success_max <- analyse_data(data         = data,
                                cutpoints    = cutpoints,
                                end_of_study = end_of_study,
                                prior        = prior,
                                N_mcmc       = N_mcmc,
                                single_arm   = single_arm,
                                method       = method,
                                alternative  = alternative,
                                h0           = h0)

  } else{
    success_max <- NA
  }

  return(list(
    success_now = success_now,
    success_max = success_max
  ))

}
