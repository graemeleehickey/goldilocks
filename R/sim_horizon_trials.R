#' Simulate trials with a prespecified endpoint-horizon schedule
#'
#' @description Simulates only reached interim looks under a prespecified schedule, using separate success and futility prediction data cuts. Use [sim_horizon_paths()] when comparing many schedules or stopping thresholds.
#' @inheritParams sim_horizon_paths
#' @inheritParams apply_horizon_rules
#' @param endpoint_horizons Endpoint horizons used to allocate reproducible simulation streams. Defaults to the unique scheduled horizons plus `final_horizon`, in increasing order. To reproduce a larger set of saved paths exactly, supply the full set of horizons used to generate those paths.
#' @details At a 12-month look with an eventual 18-month final endpoint, success imputation is fitted only to available data through month 12; futility imputation is fitted to available data through month 18. Both respect the look's calendar date. The schedule changes only if the trial reaches the corresponding look. It is not selected using interim outcomes.
#'
#' Later interim calculations are skipped after a stop. The selected cohort's final analysis and the maximum-\eqn{N} final analysis are retained, including the counterfactual futility comparison. For the same seed, set of horizons, and other simulation settings, results agree with [sim_horizon_paths()] followed by [apply_horizon_rules()]. Event-time imputation, dropout interpretation, horizon-specific analysis settings, and the endpoint lock after a stop are described in [sim_horizon_paths()].
#' @return A simulation result as described in [apply_horizon_rules()]. The `calculation_failures` component additionally records failed simulation stages. The complete trial paths are not retained.
#' @export
#' @examples
#' schedule <- data.frame(N = c(20, 40), horizon = c(12, 18))
#' result <- sim_horizon_trials(
#'   hazard_treatment = 0.04, hazard_control = 0.07,
#'   N_total = 60, endpoint_schedule = schedule, final_horizon = 18,
#'   lambda = 3, method = "bayes-bin", bin_method = "quadrature",
#'   alternative = "less", N_impute = 5, N_trials = 2, seed = 123
#' )
#' summarise_sims(result)
sim_horizon_trials <- function(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = NULL,
  N_total,
  lambda = 0.3,
  lambda_time = NULL,
  endpoint_schedule,
  final_horizon,
  endpoint_horizons = sort(unique(c(endpoint_schedule$horizon, final_horizon))),
  followup_horizon = final_horizon,
  dropout_horizon = followup_horizon,
  horizon_settings = NULL,
  prior_surv = c(0.1, 0.1),
  prior_bin = c(1, 1),
  bin_method = "mc",
  block = 2,
  rand_ratio = c(control = 1, treatment = 1),
  prop_loss = 0,
  alternative = "greater",
  h0 = 0,
  N_impute = 500,
  N_mcmc = 1000,
  mc_conf_level = 0.95,
  N_trials = 10,
  method = "logrank",
  imputed_final = FALSE,
  empty_interval = c("prior", "propagate", "error"),
  ncores = 1L,
  backend = c("auto", "fork", "psock", "sequential"),
  seed = NULL,
  prior_surv_final = prior_surv,
  generation_cutpoints = cutpoints,
  trial_offset = 0L,
  Fn = 0.05,
  Sn = 0.9,
  Qn = 1,
  prob_ha = 0.95,
  return_trace = FALSE
) {
  Call <- match.call()
  if (is.function(Call[[1L]])) {
    Call[[1L]] <- quote(sim_horizon_trials)
  }
  args <- capture_arguments(sim_horizon_trials, environment())
  if (
    !is.data.frame(endpoint_schedule) ||
      !all(c("N", "horizon") %in% names(endpoint_schedule))
  ) {
    stop("'endpoint_schedule' must be a data frame with columns N and horizon")
  }
  args["interim_look"] <- list(
    if (nrow(endpoint_schedule)) endpoint_schedule$N else NULL
  )
  args[c(
    "endpoint_schedule",
    "Fn",
    "Sn",
    "Qn",
    "prob_ha",
    "return_trace"
  )] <- NULL
  args <- validate_horizon_arguments(args)
  schedule <- validate_endpoint_schedule(endpoint_schedule, args)
  rules <- normalize_horizon_rules(args, Fn, Sn, Qn, prob_ha)
  validate_logical_scalar(return_trace, "return_trace")
  paths <- run_horizon_simulations(args, Call, schedule, rules)
  out <- apply_horizon_rules(paths, schedule, Fn, Sn, Qn, prob_ha, return_trace)
  out$call <- Call
  out$calculation_failures <- paths$failures
  out
}
