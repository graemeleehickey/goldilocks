# Graeme L. Hickey
# Run from the package root:
# Rscript benchmarks/trial-paths.R 100 20 /tmp/trial-path-benchmarks.csv
# These are implementation and performance comparisons, not design calibration.
pkgload::load_all(".", quiet = TRUE)
cli <- commandArgs(trailingOnly = TRUE)
n_trials <- if (length(cli)) as.integer(cli[1]) else 100L
n_impute <- if (length(cli) >= 2L) as.integer(cli[2]) else 20L
output <- if (length(cli) >= 3L) cli[3] else tempfile(fileext = ".csv")
stopifnot(n_trials > 1, n_impute > 1)
grid <- expand.grid(
  Fn = c(.05, .1),
  Sn = c(.85, .95),
  prob_ha = c(.975, .99),
  Qn = 1
)
rows <- list()
for (method in c("logrank", "bayes-bin")) {
  for (scenario in c("null", "alternative")) {
    fixed <- list(
      hazard_treatment = if (scenario == "null") .08 else .05,
      hazard_control = .08,
      N_total = 80,
      interim_look = c(40, 60),
      end_of_study = 12,
      lambda = 4,
      method = method,
      bin_method = "mc",
      alternative = "less",
      prop_loss = 0,
      N_trials = n_trials,
      N_impute = n_impute,
      N_mcmc = 100,
      backend = "sequential",
      seed = 7201
    )
    build_seconds <- system.time(bank <- do.call(sim_trial_paths, fixed))[[
      "elapsed"
    ]]
    for (i in seq_len(nrow(grid))) {
      rule <- as.list(grid[i, ])
      replay_seconds <- system.time(
        replayed <- do.call(apply_stopping_rules, c(list(paths = bank), rule))
      )[["elapsed"]]
      direct_args <- fixed
      direct_args$seed <- 7301
      direct_seconds <- system.time(
        direct <- do.call(sim_trials, c(direct_args, rule))
      )[["elapsed"]]
      bank_oc <- summarise_sims(replayed)
      direct_oc <- summarise_sims(direct)
      # Independent runs: display differences relative to their combined MCSE.
      difference <- bank_oc$power - direct_oc$power
      combined_mcse <- sqrt(bank_oc$power_mcse^2 + direct_oc$power_mcse^2)
      rows[[length(rows) + 1L]] <- data.frame(
        method,
        scenario,
        candidate = i,
        grid[i, ],
        n_trials,
        n_impute,
        build_seconds,
        replay_seconds,
        direct_seconds,
        bank_bytes = as.numeric(object.size(bank)),
        bank_power = bank_oc$power,
        direct_power = direct_oc$power,
        power_difference = difference,
        combined_mcse,
        bank_futility = bank_oc$stop_futility,
        direct_futility = direct_oc$stop_futility,
        bank_mean_N = bank_oc$mean_N,
        direct_mean_N = direct_oc$mean_N,
        bank_failures = bank_oc$n_failed,
        direct_failures = direct_oc$n_failed,
        failed_bank_calculations = nrow(bank$failures),
        package_version = as.character(utils::packageVersion("goldilocks")),
        R_version = as.character(getRversion())
      )
    }
    group <- do.call(rbind, tail(rows, nrow(grid)))
    cat(
      method,
      scenario,
      ": bank + eight replays =",
      round(build_seconds + sum(group$replay_seconds), 3),
      "s; eight direct runs =",
      round(sum(group$direct_seconds), 3),
      "s\n"
    )
  }
}
results <- do.call(rbind, rows)
write.csv(results, output, row.names = FALSE)
cat("Results:", output, "\n")
