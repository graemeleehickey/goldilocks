#' @title Simulate one or more clinical trials subject to known design
#'   parameters and treatment effect
#'
#' @description Simulate multiple clinical trials with fixed input parameters,
#'   and tidily extract the relevant data to generate operating characteristics.
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
#'   To use will multiple cores (where available), the argument \code{ncores}
#'   can be increased from the default of 1. Note: on Windows machines, it is
#'   not possible to use the \code{\link[parallel]{mclapply}} function with
#'   \code{ncores} \eqn{>1}.
#'
#' @return Data frame with 1 row per simulated trial and columns for key summary
#'   statistics. See \code{\link{survival_adapt}} for details of what is
#'   returned in each row.
#'
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmclapply
#' @export
#'
#' @examples
#' hc <- prop_to_haz(c(0.20, 0.30), c(0, 12), 36)
#' ht <- prop_to_haz(c(0.05, 0.15), c(0, 12), 36)
#'
#' out <- sim_trials(
#'   hazard_treatment = ht,
#'   hazard_control = hc,
#'   cutpoints = c(0, 12),
#'   N_total = 600,
#'   lambda = 20,
#'   lambda_time = 0,
#'   interim_look = c(400, 500),
#'   end_of_study = 36,
#'   prior = c(0.1, 0.1),
#'   block = 2,
#'   rand_ratio = c(1, 1),
#'   prop_loss = 0.30,
#'   alternative = "two.sided",
#'   h0 = 0,
#'   Fn = 0.05,
#'   Sn = 0.9,
#'   prob_ha = 0.975,
#'   N_impute = 5,
#'   N_mcmc = 5,
#'   N_trials = 2,
#'   method = "logrank",
#'   ncores = 1)

sim_trials <- function(
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
  alternative       = "two.sided",
  h0                = 0,
  Fn                = 0.1,
  Sn                = 0.9,
  prob_ha           = 0.95,
  N_impute          = 10,
  N_mcmc            = 10,
  N_trials          = 10,
  method            = "logrank",
  imputed_final     = FALSE,
  ncores            = 1L
) {

  Call <- match.call()

  # Check: missing 'ncores' defaults to maximum available (spare 1)
  if (missing(ncores)) {
    ncores <- min(1, parallel::detectCores() - 1)
  }

  # Check: cannot specify <1 core
  if (ncores < 1) {
    warning("Must use at least 1 core... setting ncores = 1")
  }

  # Check: if Windows and if ncores = 1
  if (.Platform$OS.type == "Windows" & ncores > 1L) {
    message("On Windows machines it is required that ncores = 1L")
    ncores <- 1
  }

  survival_adapt_wrapper <- function(x) {
    survival_adapt(
      hazard_treatment = hazard_treatment,
      hazard_control   = hazard_control,
      cutpoints        = cutpoints,
      N_total          = N_total,
      lambda           = lambda,
      lambda_time      = lambda_time,
      interim_look     = interim_look,
      end_of_study     = end_of_study,
      prior            = prior,
      block            = block,
      rand_ratio       = rand_ratio,
      prop_loss        = prop_loss,
      alternative      = alternative,
      h0               = h0,
      Fn               = Fn,
      Sn               = Sn,
      prob_ha          = prob_ha,
      N_impute         = N_impute,
      N_mcmc           = N_mcmc,
      method           = method,
      imputed_final    = imputed_final,
      debug            = FALSE)
  }

  sims <- pbmclapply(1:N_trials, survival_adapt_wrapper, mc.cores = ncores)

  sims <- do.call("rbind", sims)
  out <- list(sims = sims, call = Call)

  return(out)

}
