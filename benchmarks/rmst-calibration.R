# Run from the package root:
# Rscript benchmarks/rmst-calibration.R 500 /tmp/rmst-calibration.csv
# Increase the trial/imputation counts for design-specific threshold selection.
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
n_trials <- if (length(args)) as.integer(args[1]) else 500L
output <- if (length(args) >= 2L) args[2] else tempfile(fileext = ".csv")
stopifnot(is.finite(n_trials), n_trials > 0)
tau <- 9
# Analytic two-interval RMST, independent of the fitted analysis.
rmst_truth <- function(h, split = 3, horizon = tau) {
  -expm1(-h[1] * split) /
    h[1] +
    exp(-h[1] * split) * -expm1(-h[2] * (horizon - split)) / h[2]
}
control <- c(0.1, 0.1)
crossing_late <- uniroot(
  function(h) {
    rmst_truth(c(0.04, h)) - rmst_truth(control)
  },
  c(0.1, 1),
  tol = 1e-12
)$root
scenarios <- list(
  equal_survival = list(h = control, h0 = 0),
  equal_rmst_crossing = list(h = c(0.04, crossing_late), h0 = 0),
  nonzero_null = list(
    h = c(0.12, 0.12),
    h0 = rmst_truth(c(0.12, 0.12)) - rmst_truth(control)
  ),
  delayed_benefit = list(h = c(0.1, 0.04), h0 = 0)
)
variants <- list(
  fixed = list(
    interim_look = NULL,
    prop_loss = 0,
    imputed_final = FALSE,
    cutpoints = 3
  ),
  adaptive = list(
    interim_look = 100,
    prop_loss = 0,
    imputed_final = FALSE,
    cutpoints = 3
  ),
  adaptive_dropout = list(
    interim_look = 100,
    prop_loss = c(control = 0.1, treatment = 0.2),
    imputed_final = FALSE,
    cutpoints = 3
  ),
  adaptive_imputed = list(
    interim_look = 100,
    prop_loss = c(control = 0.1, treatment = 0.2),
    imputed_final = TRUE,
    cutpoints = 3
  ),
  adaptive_misspecified = list(
    interim_look = 100,
    prop_loss = 0,
    imputed_final = FALSE,
    cutpoints = NULL
  )
)
rows <- list()
for (scenario in names(scenarios)) {
  for (variant in names(variants)) {
    design <- scenarios[[scenario]]
    seed <- 1430L + length(rows)
    result <- do.call(
      sim_trials,
      c(
        list(
          hazard_control = control,
          hazard_treatment = design$h,
          generation_cutpoints = 3,
          N_total = 200,
          lambda = 8,
          end_of_study = 12,
          rmst_tau = tau,
          method = "rmst",
          h0 = design$h0,
          alternative = "greater",
          prob_ha = 0.975,
          N_impute = 100,
          N_trials = n_trials,
          seed = seed,
          backend = "sequential"
        ),
        variants[[variant]]
      )
    )
    success <- result$sims$trial_success
    rate <- mean(success)
    ci <- binom.test(sum(success), length(success))$conf.int
    row <- data.frame(
      scenario = scenario,
      variant = variant,
      seed = seed,
      true_difference = rmst_truth(design$h) - rmst_truth(control),
      h0 = design$h0,
      requested = n_trials,
      completed = length(success),
      failed = nrow(result$failures),
      success_rate = rate,
      mcse = sqrt(rate * (1 - rate) / length(success)),
      lower = ci[1],
      upper = ci[2],
      mean_N = mean(result$sims$N_enrolled)
    )
    rows[[length(rows) + 1L]] <- row
    # Write after each scenario so a long calibration can be inspected.
    write.csv(do.call(rbind, rows), output, row.names = FALSE)
    print(row, row.names = FALSE)
  }
}
message("Calibration results: ", output)
