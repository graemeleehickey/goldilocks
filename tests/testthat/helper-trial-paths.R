path_test_args <- function(...) {
  args <- list(
    hazard_treatment = 0.05,
    hazard_control = 0.08,
    N_total = 80,
    interim_look = c(40, 60),
    lambda = 4,
    end_of_study = 12,
    alternative = "less",
    method = "bayes-bin",
    bin_method = "quadrature",
    N_impute = 4,
    N_mcmc = 30,
    N_trials = 2,
    seed = 914
  )
  extra <- list(...)
  args[names(extra)] <- extra
  args
}

# Direct early-stopping reference: generate the data and evaluate only reached
# looks, using the bank's documented substream layout but no saved evidence.
direct_path_reference <- function(args, rules, trial = 1L) {
  restore_rng <- preserve_path_rng()
  on.exit(restore_rng(), add = TRUE)
  config <- args
  config$single_arm <- is.null(config$hazard_control)
  if (!config$method %in% c("bayes-bin", "bayes-surv")) {
    config$N_mcmc <- 1L
  }
  stream <- make_rng_streams(args$seed, trial)[[trial]]
  assign(".Random.seed", stream, envir = .GlobalEnv)
  data <- call_path_function(sim_comp_data, config)
  stopped <- config$N_total
  decision <- "continue"
  trace <- list()
  for (i in seq_along(config$interim_look)) {
    n <- config$interim_look[i]
    stream <- parallel::nextRNGSubStream(stream)
    assign(".Random.seed", stream, envir = .GlobalEnv)
    result <- call_path_function(
      evaluate_interim_decision,
      config,
      data_interim = prepare_simulated_interim(data, n, config$end_of_study),
      look = i,
      planned_N = n,
      calendar_time = data$enrollment[n],
      active_followup = active_followup_at(data, data$enrollment[n]),
      Fn = rules$Fn[i],
      Sn = rules$Sn[i],
      Qn = rules$Qn[i],
      prob_ha = rules$prob_ha,
      check_futility = TRUE
    )
    trace[[i]] <- result$trace
    decision <- result$decision
    if (decision != "continue") {
      stopped <- n
      break
    }
    # Reserve the cohort-final substream without running that unused analysis.
    stream <- parallel::nextRNGSubStream(stream)
  }
  cohort <- data[data$id <= stopped, , drop = FALSE]
  cohort$subject_impute_success <- cohort$event == 0 &
    cohort$time < config$end_of_study
  stream <- parallel::nextRNGSubStream(stream)
  assign(".Random.seed", stream, envir = .GlobalEnv)
  final <- if (decision == "stop_immediate_success") {
    c(NA_real_, NA_real_)
  } else {
    call_path_function(analyse_final, config, data_in = cohort)
  }
  list(
    N = stopped,
    score = final[1],
    effect = final[2],
    success = if (decision == "stop_immediate_success") {
      TRUE
    } else if (decision == "stop_futility") {
      FALSE
    } else {
      final[1] > rules$prob_ha
    },
    trace = new_trial_trace(trace),
    calendar = trial_calendar_metrics(cohort, config$end_of_study)
  )
}
