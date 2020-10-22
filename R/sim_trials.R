#' @title Simulate one or more clinical trials subject to known design
#'   parameters and treatment effect
#'
#' @description Simulate multiple clinical trials with fixed input parameters,
#'   and tidly extract the relevant data to generate operating characteristics.
#'
#' @inheritParams survival_adapt
#' @param N_trials integer. Number of trials to simulate.
#' @param ncores integer. Number of cores to use for parallel processing.
#'
#' @details This is basically a wrapper function for
#'   \code{\link{survival_adapt}}, whereby we repeatedly run the function for a
#'   independent number of trials (all with the same input design parameters and
#'   treatment effect).
#'
#'   To use will mutiple cores (where available), the argument \code{ncores}
#'   can be increased from the default of 1. Note: on Windows machines, it is not
#'   possible to use the \code{\link[parallel]{mclappy}} function with \code{ncores}
#'   \eqn{>1}.
#'
#' @return Data frame with 1 row per simulated trial and columns for key summary
#'   statistics.
#'
#' @importFrom parallel mclapply detectCores
#' @export
#'
#' @examples
#' hc <- haz_est(c(0.20, 0.30), 12, 36)
#' ht <- haz_est(c(0.05, 0.15), 12, 36)
#'
#' out <- sim_trials(
#'   hazard_treatment = ht,
#'   hazard_control = hc,
#'   cutpoint = c(0, 12),
#'   N_total = 600,
#'   lambda = 20,
#'   lambda_time = NULL,
#'   interim_look = c(400, 500),
#'   end_of_study = 36,
#'   prior = c(.1, .1),
#'   block = 2,
#'   rand_ratio = c(1, 1),
#'   prop_loss_to_followup = 0.30,
#'   alternative = "less",
#'   h0 = 0,
#'   futility_prob = 0.05,
#'   expected_success_prob = 0.9,
#'   prob_ha = 0.975,
#'   N_impute = 5,
#'   N_mcmc = 5,
#'   N_trials = 2,
#'   method = "logrank",
#'   ncores = 1)

sim_trials <- function(
  hazard_treatment,
  hazard_control        = NULL,
  cutpoint              = 0,
  N_total,
  lambda                = 0.3,
  lambda_time           = NULL,
  interim_look          = NULL,
  end_of_study,
  prior                 = c(.1, .1),
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
  N_trials              = 10,
  method                = "logrank",
  imputed_final         = TRUE,
  ncores                = 1L
) {

  # Check: missing 'ncores' defaults to max available (spare 1)
  if (missing(ncores)) {
    ncores <- min(1, parallel::detectCores() - 1)
  }

  # Check: cannot specify <1 core
  if (ncores < 1) {
    warning("Must use at least 1 core... setting ncores = 1")
  }

  # Check: if Windows and if ncores = 1
  if (.Platform$OS.type == "Windows" & ncores > 1L) {
    stop("On Windows machines it is required that ncores = 1L")
  }

  out <- mclapply(1:N_trials, survival_adapt, mc.cores = ncores,
                  hazard_treatment      = hazard_treatment,
                  hazard_control        = hazard_control,
                  cutpoint              = cutpoint,
                  N_total               = N_total,
                  lambda                = lambda,
                  lambda_time           = lambda_time,
                  interim_look          = interim_look,
                  end_of_study          = end_of_study,
                  prior                 = prior,
                  block                 = block,
                  rand_ratio            = rand_ratio,
                  prop_loss_to_followup = prop_loss_to_followup,
                  alternative           = alternative,
                  h0                    = h0,
                  futility_prob         = futility_prob,
                  expected_success_prob = expected_success_prob,
                  prob_ha               = prob_ha,
                  N_impute              = N_impute,
                  N_mcmc                = N_mcmc,
                  method                = method,
                  imputed_final         = imputed_final)

  out <- do.call("rbind", out)

  return(out)

}
