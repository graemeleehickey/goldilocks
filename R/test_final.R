#' @title Run final test once enrollment is complete
#'
#' @description Once enrollment is complete, the final analysis must be
#'   conducted. The final analysis can take two forms depending on whether the
#'   data is "complete" or not (see Details).
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param data_in data frame. The time-to-event data that requires imputation,
#'   with columns for treatment arm (\code{treatment}, which can be all 1s if
#'   single-arm), event time (\code{time}), event indicator (\code{event}),
#'   indicator of whether the subject requires imputation for expected success
#'   (\code{subject_impute_success}).
#'
#' @details If \code{prop_loss} is $>0$, then there is drop-out (loss to
#'   follow-up) for some subjects. During the Goldilocks algorithm, all subjects
#'   who are event-free and still eligible for follow-up for the primary
#'   endpoint are imputed. This means subjects who are lost to follow-up will be
#'   imputed. In the final analysis we need to decide whether we mimic that, or
#'   whether we will exclude or right-censor such subjects (depending on
#'   analysis \code{method} choice). This is controlled by \code{imputed_final},
#'   which is a Boolean: \code{TRUE} implies the final analysis will be based on
#'   fully imputed data, whereas \code{FALSE} implies the final analysis will be
#'   based on observed (non-imputed) data only.
#'
#'   If \code{imputed_final = FALSE} then intuitively, we might expect to see a
#'   marginal increase in the proportion of studies that stop for expected
#'   success, but which then go on to fail. We have not verified this aspect,
#'   but it should be noted.
#'
#' @return Vector with 1) posterior probability (or P-value equivalent) for
#'   alternative hypothesis, and 2) mean posterior treatment effect.
#' @noRd
test_final <- function(data_in, cutpoints, prior, N_mcmc, single_arm,
                       imputed_final, method, N_impute, alternative,
                       h0, end_of_study) {

  if (imputed_final) {
    # Posterior distribution of lambdas: final data
    post_lambda_final <- posterior(data       = data_in,
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
        data_in      = data_in,
        hazard       = post_lambda_final[j, , , drop = FALSE],
        end_of_study = end_of_study,
        cutpoints    = cutpoints,
        type         = "success",
        single_arm   = single_arm)

      # Create enrolled subject data frame for analysis
      time      <- NULL
      event     <- NULL
      treatment <- NULL
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
    success <- analyse_data(data         = data_in,
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

  return(c(post_paa, est_final))

}
