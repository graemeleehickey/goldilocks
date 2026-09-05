#' Simulation results for the prob_ha calibration vignette
#'
#' Run from the package root after loading the development version of
#' goldilocks. The saved object contains summaries rather than trial-level
#' simulations so the installed vignette remains small.

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(goldilocks)
}

target_type1 <- 0.025
prob_ha_grid <- seq(0.965, 0.995, by = 0.005)
screening_seed <- 67201
validation_seed <- 67202

null_hazard <- prop_to_haz(0.30, endtime = 12)

calibration_design <- list(
  hazard_treatment = null_hazard,
  hazard_control = null_hazard,
  cutpoints = NULL,
  N_total = 300,
  lambda = 5,
  lambda_time = NULL,
  interim_look = seq(100, 275, 25),
  end_of_study = 12,
  prior_surv = c(0.1, 0.1),
  block = 2,
  rand_ratio = c(control = 1, treatment = 1),
  prop_loss = 0,
  alternative = "less",
  h0 = 0,
  Fn = rep(0.10, 8),
  Sn = c(1, rep(0.90, 7)),
  Qn = 1,
  N_impute = 100,
  method = "logrank",
  return_trace = FALSE,
  ncores = 8,
  backend = "auto"
)

screening_results <- lapply(prob_ha_grid, function(threshold) {
  do.call(
    sim_trials,
    c(
      calibration_design,
      list(
        prob_ha = threshold,
        N_trials = 2000,
        seed = screening_seed
      )
    )
  )
})
names(screening_results) <- sprintf("prob_ha = %.3f", prob_ha_grid)

calibration_screening <- summarise_sims(screening_results)
calibration_screening$prob_ha <- prob_ha_grid
calibration_screening$type1_status <- ifelse(
  calibration_screening$power_mc_upper <= target_type1,
  "controlled",
  ifelse(
    calibration_screening$power_mc_lower > target_type1,
    "not controlled",
    "inconclusive"
  )
)

controlled <- calibration_screening$prob_ha[
  calibration_screening$type1_status == "controlled"
]
if (length(controlled) == 0L) {
  stop("The screening grid did not identify a candidate for validation")
}
selected_prob_ha <- min(controlled)

validation_result <- do.call(
  sim_trials,
  c(
    calibration_design,
    list(
      prob_ha = selected_prob_ha,
      N_trials = 10000,
      seed = validation_seed
    )
  )
)
calibration_validation <- summarise_sims(list(
  "fresh-seed validation" = validation_result
))
calibration_validation$prob_ha <- selected_prob_ha
calibration_validation$type1_status <- ifelse(
  calibration_validation$power_mc_upper <= target_type1,
  "controlled",
  ifelse(
    calibration_validation$power_mc_lower > target_type1,
    "not controlled",
    "inconclusive"
  )
)

calibration_settings <- list(
  target_type1 = target_type1,
  prob_ha_grid = prob_ha_grid,
  screening_trials = 2000L,
  validation_trials = 10000L,
  N_impute = calibration_design$N_impute,
  screening_seed = screening_seed,
  validation_seed = validation_seed,
  selected_prob_ha = selected_prob_ha
)

save(
  calibration_screening,
  calibration_validation,
  calibration_settings,
  file = "vignettes/prob-ha-calibration.rda"
)
