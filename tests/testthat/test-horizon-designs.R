horizon_test_args <- function(...) {
  args <- list(
    hazard_treatment = 0.05,
    hazard_control = 0.08,
    N_total = 80,
    endpoint_schedule = horizon_test_schedule(),
    lambda = 3,
    final_horizon = 18,
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

horizon_test_schedule <- function() {
  data.frame(N = c(40, 60), horizon = c(12, 18))
}

test_that("conjugate survival horizon paths ignore the final imputation flag", {
  args <- horizon_test_args(
    method = "bayes-surv",
    cutpoints = c(6, 15),
    generation_cutpoints = NULL,
    prior_surv = c(0.1, 0.1),
    prior_surv_final = c(1, 5),
    prop_loss = 0.35,
    N_trials = 1
  )
  args$imputed_final <- FALSE
  direct <- do.call(sim_horizon_paths, args)
  args$imputed_final <- TRUE
  flagged <- do.call(sim_horizon_paths, args)
  expect_identical(flagged$trials, direct$trials)
  expect_identical(flagged$failures, direct$failures)
  expect_equal(nrow(direct$failures), 0L)
})

test_that("endpoint recensoring preserves event and dropout ordering", {
  data <- data.frame(
    time = c(15, 15, 8, 12, 18),
    event = c(1, 0, 0, 1, 0),
    loss_to_fu = c(FALSE, TRUE, TRUE, FALSE, FALSE)
  )
  cut <- censor_horizon_data(data, 12)
  expect_equal(cut$time, c(12, 12, 8, 12, 12))
  expect_equal(cut$event, c(0, 0, 0, 1, 0))
  expect_identical(cut$loss_to_fu, c(FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(censor_horizon_data(data, 18), data)
})

test_that("late observed data inform futility but cannot inform earlier success", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  args <- attr(bank, "arguments")
  config12 <- horizon_analysis_config(args, 12)
  config18 <- horizon_analysis_config(args, 18)
  # A single shared hazard interval crosses month 12: mere re-analysis of an
  # 18-month posterior would fail this test even without an extra cutpoint.
  data <- data.frame(
    id = 1:80,
    treatment = rep(c(0, 1), 40),
    time = rep(18, 80),
    event = rep(0, 80),
    enrollment = c(rep(0, 39), 20:60),
    loss_to_fu = FALSE
  )
  late <- data
  late$time[1:20] <- 15
  late$event[1:20] <- 1
  evaluate <- function(data, config) {
    set.seed(45)
    horizon_interim_evidence(data, config, 1, 40)
  }
  expect_identical(evaluate(data, config12), evaluate(late, config12))
  original18 <- evaluate(data, config18)
  changed18 <- evaluate(late, config18)
  expect_false(identical(original18$scores$maximum, changed18$scores$maximum))
  expect_equal(sum(changed18$diagnostics$posterior$observed_events), 20)
  expect_equal(sum(original18$diagnostics$posterior$observed_events), 0)
  # Outcomes after the calendar cut, including outcomes of future recruits,
  # cannot affect either predictive calculation.
  future <- late
  future$time[40:80] <- 2
  future$event[40:80] <- 1
  expect_identical(evaluate(late, config12), evaluate(future, config12))
  expect_identical(evaluate(late, config18), evaluate(future, config18))
})

test_that("each horizon uses the leading master-prior intervals", {
  args <- horizon_test_args(
    cutpoints = c(6, 12, 15),
    generation_cutpoints = NULL,
    prior_surv = matrix(c(1, 10, 2, 20, 3, 30, 4, 40), nrow = 2),
    N_trials = 1
  )
  bank <- do.call(sim_horizon_paths, args)
  config <- horizon_analysis_config(attr(bank, "arguments"), 12)
  expect_equal(config$cutpoints, 6)
  expect_equal(
    unname(config$prior_surv[,, "control"]),
    matrix(c(1, 10, 2, 20), 2)
  )
  posterior <- bank$trials[[1]]$interims[[1]][[1]]$diagnostics$posterior
  expect_setequal(posterior$interval_end, c(6, 12))
  expect_setequal(posterior$prior_shape, c(1, 2))
  args$prior_surv[, 3:4] <- 100
  other <- do.call(sim_horizon_paths, args)
  expect_identical(
    bank$trials[[1]]$interims[[1]][[1]],
    other$trials[[1]]$interims[[1]][[1]]
  )
})

test_that("an endpoint inside an interval retains its prior but truncates observations", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(
    N_total = 8,
    endpoint_schedule = data.frame(N = 6, horizon = 12),
    cutpoints = c(6, 15),
    generation_cutpoints = NULL,
    prior_surv = matrix(c(1, 10, 2, 20, 3, 30), nrow = 2),
    N_trials = 1
  ))
  args <- attr(bank, "arguments")
  config12 <- horizon_analysis_config(args, 12)
  config18 <- horizon_analysis_config(args, 18)
  data <- data.frame(
    id = 1:8,
    treatment = rep(c(0, 1), 4),
    time = c(10, 14, 8, 18, 18, 18, 18, 18),
    event = c(1, 1, 0, 0, 0, 0, 0, 0),
    enrollment = c(0, 0, 0, 0, 0, 20, 21, 22),
    loss_to_fu = c(FALSE, FALSE, TRUE, rep(FALSE, 5))
  )
  evaluate <- function(data, config) {
    set.seed(45)
    horizon_interim_evidence(data, config, 1, 6)
  }
  success <- evaluate(data, config12)
  futility <- evaluate(data, config18)
  interval <- function(evidence) {
    posterior <- evidence$diagnostics$posterior
    posterior <- posterior[posterior$interval_start == 6, ]
    posterior[match(c("control", "treatment"), posterior$arm), ]
  }
  early <- interval(success)
  later <- interval(futility)
  # Through month 12, control exposure is (10-6) + (8-6) + (12-6),
  # and treatment exposure is (12-6) + (12-6). The month-14 event is hidden.
  expect_equal(early$interval_end, c(12, 12))
  expect_equal(early$observed_events, c(1, 0))
  expect_equal(early$observed_exposure, c(12, 12))
  expect_equal(early$prior_shape, c(2, 2))
  expect_equal(early$prior_rate, c(20, 20))
  expect_equal(early$posterior_shape, c(3, 2))
  expect_equal(early$posterior_rate, c(32, 32))
  # The final-horizon model uses the full (6,15] interval where observed.
  expect_equal(later$interval_end, c(15, 15))
  expect_equal(later$observed_events, c(1, 1))
  expect_equal(later$observed_exposure, c(15, 17))
  expect_equal(later$prior_shape, early$prior_shape)
  expect_equal(later$prior_rate, early$prior_rate)
  expect_equal(later$posterior_shape, c(3, 3))
  expect_equal(later$posterior_rate, c(35, 37))
  changed <- data
  changed$time[2] <- 16
  expect_identical(evaluate(changed, config12), success)
  expect_false(identical(interval(evaluate(changed, config18)), later))
})

test_that("direct horizon trials agree with full-bank replay for every method", {
  for (method in c(
    "logrank",
    "cox",
    "rmst",
    "bayes-surv",
    "bayes-bin",
    "riskdiff-wald",
    "riskdiff-fm"
  )) {
    args <- horizon_test_args(
      method = method,
      N_trials = 1,
      alternative = if (method == "rmst") "greater" else "less"
    )
    bank <- do.call(sim_horizon_paths, args)
    schedule <- horizon_test_schedule()
    for (rules in list(
      list(Fn = 0, Sn = 1),
      list(Fn = 1, Sn = 1, prob_ha = 1),
      list(Fn = 0, Sn = 0, Qn = 1, prob_ha = 0),
      list(Fn = 0, Sn = 0, Qn = 0, prob_ha = 0)
    )) {
      replay <- do.call(
        apply_horizon_rules,
        c(
          list(paths = bank, endpoint_schedule = schedule, return_trace = TRUE),
          rules
        )
      )
      direct_args <- args
      direct <- do.call(
        sim_horizon_trials,
        c(
          direct_args,
          list(return_trace = TRUE),
          rules
        )
      )
      expect_identical(direct$sims, replay$sims, info = method)
      expect_identical(direct$traces, replay$traces, info = method)
    }
  }
})

test_that("constant horizons reproduce fixed-horizon evidence and trial outcomes", {
  for (method in c(
    "logrank", "cox", "rmst", "bayes-surv", "bayes-bin",
    "riskdiff-wald", "riskdiff-fm"
  )) {
    args <- horizon_test_args(
      endpoint_schedule = data.frame(N = c(40, 60), horizon = 18),
      method = method,
      N_trials = 1,
      N_impute = 8,
      cutpoints = c(6, 12),
      generation_cutpoints = NULL,
      prior_surv = c(2, 40),
      prop_loss = 0.1,
      imputed_final = !method %in% c("logrank", "riskdiff-fm"),
      alternative = if (method == "rmst") "greater" else "less"
    )
    bank <- do.call(sim_horizon_paths, args)
    old_args <- args[setdiff(
      names(args),
      c("endpoint_schedule", "final_horizon")
    )]
    old_args$interim_look <- args$endpoint_schedule$N
    old_args$end_of_study <- 18
    fixed <- do.call(sim_trial_paths, old_args)
    for (i in 1:2) {
      expect_identical(
        bank$trials[[1]]$interims[[i]][[1]],
        fixed$trials[[1]]$interims[[i]],
        info = method
      )
    }
    expect_identical(
      bank$trials[[1]]$finals[names(fixed$trials[[1]]$finals)],
      fixed$trials[[1]]$finals,
      info = method
    )
    expect_equal(nrow(bank$failures), 0L, info = method)
    schedule <- args$endpoint_schedule
    rule_sets <- list(
      maximum_sample_size = list(Fn = 0, Sn = 1),
      futility = list(Fn = 1, Sn = 1, prob_ha = 1),
      expected_success = list(Fn = 0, Sn = 0, Qn = 1, prob_ha = 0),
      immediate_success = list(Fn = 0, Sn = 0, Qn = 0, prob_ha = 0)
    )
    for (reason in names(rule_sets)) {
      info <- paste(method, reason)
      rules <- rule_sets[[reason]]
      ordinary <- do.call(apply_stopping_rules, c(
        list(paths = fixed, return_trace = TRUE), rules
      ))
      constant <- do.call(apply_horizon_rules, c(
        list(paths = bank, endpoint_schedule = schedule, return_trace = TRUE),
        rules
      ))
      expect_identical(
        constant$sims[names(ordinary$sims)], ordinary$sims, info = info
      )
      expect_identical(
        constant$traces[names(ordinary$traces)], ordinary$traces, info = info
      )
      expect_equal(constant$sims$stopping_reason, reason, info = info)
      expect_equal(constant$sims$endpoint_horizon, 18, info = info)
    }
  }
})

test_that("single-arm horizon-specific hypotheses and imputed finals are retained", {
  for (method in c("bayes-bin", "bayes-surv")) {
    settings <- data.frame(horizon = c(12, 18), h0 = c(0.4, 0.6))
    args <- horizon_test_args(
      hazard_control = NULL,
      method = method,
      horizon_settings = settings,
      prop_loss = 0.2,
      dropout_horizon = 12,
      imputed_final = TRUE,
      N_trials = 1
    )
    bank <- do.call(sim_horizon_paths, args)
    replay <- apply_horizon_rules(bank, horizon_test_schedule(), Fn = 0, Sn = 1)
    expect_equal(replay$sims$margin, 0.6)
    expect_equal(replay$sims$N_control, 0)
    expect_true(is.finite(replay$sims$post_prob_ha))
    direct <- do.call(
      sim_horizon_trials,
      c(args, list(Fn = 0, Sn = 1))
    )
    expect_identical(direct$sims, replay$sims)
  }
})

test_that("a stop locks the endpoint and keeps counterfactual futility separate", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  bank$trials[[1]]$interims[[1]][[1]]$scores$current[] <- 0.96
  bank$trials[[1]]$interims[[1]][[2]]$scores$maximum[] <- 0.1
  bank$trials[[1]]$finals$post_prob_ha <- c(0.1, 0.99, 0.1, 0.99, 0.1, 0.99)
  schedule <- horizon_test_schedule()
  out <- apply_horizon_rules(
    bank,
    schedule,
    Sn = 0.9,
    Fn = 1,
    return_trace = TRUE
  )
  expect_equal(out$sims$N_enrolled, 40)
  expect_equal(out$sims$endpoint_horizon, 12)
  expect_false(out$sims$trial_success)
  expect_true(out$sims$success_at_max)
  expect_equal(out$sims$stop_expected_success, 1)
  expect_equal(nrow(out$traces), 1)
  expect_equal(out$traces$success_data_horizon, 12)
  expect_equal(out$traces$futility_data_horizon, 18)
  expect_equal(
    out$sims$planned_completion_time,
    out$sims$accrual_stop_time + 12
  )
  expect_identical(out$sims$decision_time, out$sims$analysis_ready_time)
  expect_false(out$sims$futility_and_success_at_max)
  # Equality does not cross the success boundary; futility then applies.
  out <- apply_horizon_rules(bank, schedule, Sn = 1, Fn = 0.1)
  expect_equal(out$sims$stop_futility, 1)
  expect_true(out$sims$futility_and_success_at_max)
  expect_equal(summarise_sims(out)$futility_and_success_at_max, 1)
})

test_that("a continued trial changes endpoint before evaluating the scheduled look", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  bank$trials[[1]]$interims[[1]][[1]]$scores$current[] <- 0.1
  bank$trials[[1]]$interims[[2]][[1]]$scores$current[] <- 0.1
  bank$trials[[1]]$interims[[2]][[2]]$scores$current[] <- 0.99
  out <- apply_horizon_rules(
    bank,
    horizon_test_schedule(),
    Fn = 0,
    return_trace = TRUE
  )
  expect_equal(out$sims$N_enrolled, 60)
  expect_equal(out$sims$endpoint_horizon, 18)
  expect_equal(out$traces$success_horizon, c(12, 18))
  never <- apply_horizon_rules(bank, horizon_test_schedule(), Fn = 0, Sn = 1)
  expect_equal(never$sims$N_enrolled, 80)
  expect_equal(
    never$sims$planned_completion_time,
    never$sims$accrual_stop_time + 18
  )
})

test_that("disabled horizon looks do not require saved interim results", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  bank$trials[[1]]$interims[] <- list(NULL)
  expect_silent(out <- apply_horizon_rules(
    bank, horizon_test_schedule(), Fn = 0, Sn = 1, Qn = 1,
    return_trace = TRUE
  ))
  expect_equal(nrow(out$failures), 0)
  expect_equal(out$sims$N_enrolled, 80)
  expect_equal(out$sims$endpoint_horizon, 18)
  expect_equal(out$sims$stopping_reason, "maximum_sample_size")
  expect_true(is.na(out$sims$ppp_success))
  expect_equal(nrow(out$traces), 0)
})

test_that("unused failures do not invalidate known horizon outcomes", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  bank$trials[[1]]$interims[[1]][[1]]$scores$current[] <- 0.99
  bank$trials[[1]]$interims[[1]][[2]] <- NULL
  # Preserve the horizon slot when recording a failed model.
  bank$trials[[1]]$interims[[1]][2] <- list(NULL)
  bank$trials[[1]]$finals$post_prob_ha[6] <- NA_real_
  out <- apply_horizon_rules(bank, horizon_test_schedule(), Qn = 0.9)
  expect_true(out$sims$trial_success)
  expect_true(is.na(out$sims$success_at_max))
  expect_false(out$sims$futility_and_success_at_max)
  expect_error(
    apply_horizon_rules(bank, horizon_test_schedule(), Sn = 1, Fn = 1),
    "F.*prediction|Maximum-cohort"
  )
  bank$trials[[1]]$interims[[1]][[1]]$scores$current[] <- 0.1
  bank$trials[[1]]$interims[[1]][2] <- list(list(
    scores = list(maximum = rep(0.1, 4)),
    diagnostics = list()
  ))
  out <- apply_horizon_rules(bank, horizon_test_schedule(), Sn = 1, Fn = 1)
  expect_true(is.na(out$sims$futility_and_success_at_max))
  expect_equal(summarise_sims(out)$futility_and_success_at_max_n_missing, 1)
})

test_that("direct simulation skips unreached interim fits", {
  original <- horizon_interim_evidence
  calls <- integer()
  local_mocked_bindings(horizon_interim_evidence = function(
    data,
    config,
    look,
    n
  ) {
    calls <<- c(calls, look)
    original(data, config, look, n)
  })
  args <- horizon_test_args(N_trials = 1)
  out <- do.call(
    sim_horizon_trials,
    c(
      args,
      list(
        Sn = 1,
        Fn = 1,
        prob_ha = 1
      )
    )
  )
  expect_equal(out$sims$stop_futility, 1)
  expect_identical(calls, c(1L, 1L))
})

test_that("horizon banks are reproducible and serializable", {
  set.seed(23)
  before <- .Random.seed
  bank <- do.call(sim_horizon_paths, horizon_test_args())
  expect_identical(.Random.seed, before)
  expect_named(attr(bank, "arguments"), names(formals(sim_horizon_paths)))
  regenerated <- do.call(sim_horizon_paths, attr(bank, "arguments"))
  expect_identical(regenerated$trials, bank$trials)
  first <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  second <- do.call(
    sim_horizon_paths,
    horizon_test_args(N_trials = 1, seed = 915)
  )
  expect_identical(bank$trials[1], first$trials)
  expect_identical(first$trials[[1]]$trial, 1L)
  expect_identical(second$trials[[1]]$trial, 1L)
  expect_false(identical(first$trials, second$trials))
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  saveRDS(bank, file)
  restored <- readRDS(file)
  original <- apply_horizon_rules(
    bank,
    horizon_test_schedule(),
    return_trace = TRUE
  )
  replay <- apply_horizon_rules(
    restored,
    horizon_test_schedule(),
    return_trace = TRUE
  )
  expect_identical(replay$sims, original$sims)
  expect_identical(replay$traces, original$traces)
  expect_identical(.Random.seed, before)
  expect_output(print(bank), "2 trials; 2 interim looks; horizons: 12, 18")
})

test_that("horizon banks agree across computation methods", {
  sequential <- do.call(
    sim_horizon_paths,
    horizon_test_args(backend = "sequential")
  )
  for (backend in c("psock", if (.Platform$OS.type != "windows") "fork")) {
    parallel <- do.call(
      sim_horizon_paths,
      horizon_test_args(backend = backend, ncores = 2)
    )
    expect_identical(parallel$trials, sequential$trials, info = backend)
    expect_identical(parallel$failures, sequential$failures, info = backend)
  }
})

test_that("the schedule derives all saved looks and includes the final horizon", {
  schedule <- data.frame(N = c(40, 60), horizon = c(12, 12))
  args <- horizon_test_args(endpoint_schedule = schedule, N_trials = 1)
  bank <- do.call(sim_horizon_paths, args)
  expect_identical(bank$schema_version, 2L)
  expect_identical(attr(bank, "arguments")$endpoint_schedule, schedule)
  for (look in bank$trials[[1]]$interims) {
    expect_length(look, 2)
    expect_equal(look[[1]]$diagnostics$posterior$interval_end, c(12, 12))
    expect_equal(look[[2]]$diagnostics$posterior$interval_end, c(18, 18))
  }
  expect_equal(bank$trials[[1]]$finals$horizon, rep(c(12, 18), 3))
  replay <- apply_horizon_rules(bank, schedule, Fn = 0, Sn = 1)
  direct <- do.call(sim_horizon_trials, c(args, list(Fn = 0, Sn = 1)))
  expect_identical(direct$sims, replay$sims)
  # A different switch time uses the same saved evidence, without new fits.
  switched <- apply_horizon_rules(bank, horizon_test_schedule(), Fn = 0, Sn = 1)
  expect_identical(switched$sims, replay$sims)
  expect_equal(replay$sims$followup_horizon, 18)
})

test_that("previously saved horizon paths retain trial numbers and omit collection outputs", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  removed <- c(
    "collection_ready_time", "planned_collection_time", "collection_person_time"
  )
  expect_false(any(removed %in% names(bank$trials[[1]]$finals)))
  bank$trials[[1]]$trial <- 101L
  expected <- apply_horizon_rules(bank, horizon_test_schedule(), return_trace = TRUE)
  expect_false(any(removed %in% names(expected$sims)))
  expect_identical(expected$sims$trial, 101L)
  previous <- bank
  # Earlier paths stored continued-collection summaries at each stopping cohort.
  finals <- previous$trials[[1]]$finals
  final_horizon_rows <- finals[finals$horizon == 18, ]
  row <- match(finals$planned_N, final_horizon_rows$planned_N)
  finals$collection_ready_time <- final_horizon_rows$analysis_ready_time[row]
  finals$planned_collection_time <- final_horizon_rows$planned_completion_time[row]
  finals$collection_person_time <- final_horizon_rows$followup_person_time[row]
  previous$trials[[1]]$finals <- finals
  attr(previous, "arguments")$trial_offset <- 100L
  actual <- apply_horizon_rules(previous, horizon_test_schedule(), return_trace = TRUE)
  expect_identical(actual$sims, expected$sims)
  expect_identical(actual$traces, expected$traces)

  # The original format also stored looks, horizons, and collection duration.
  previous$schema_version <- 1L
  args <- validate_horizon_arguments(attr(previous, "arguments"))
  old_names <- c(
    setdiff(names(formals(sim_horizon_paths)), "endpoint_schedule"),
    "interim_look", "endpoint_horizons", "followup_horizon", "trial_offset"
  )
  attr(previous, "arguments") <- args[old_names]
  actual <- apply_horizon_rules(previous, horizon_test_schedule(), return_trace = TRUE)
  expect_identical(actual$sims, expected$sims)
  expect_identical(actual$traces, expected$traces)
  expect_output(print(previous), "2 interim looks; horizons: 12, 18")
})

test_that("simulation functions no longer accept trial offsets", {
  for (simulate in list(sim_trial_paths, sim_horizon_paths, sim_horizon_trials)) {
    expect_error(simulate(trial_offset = 1L), "unused argument")
  }
})

test_that("dropout probability retains its fixed reference horizon", {
  a <- horizon_test_args(prop_loss = 0.2, dropout_horizon = 12, N_trials = 1)
  b <- horizon_test_args(
    prop_loss = -expm1(log1p(-0.2) * 18 / 12),
    dropout_horizon = 18,
    N_trials = 1
  )
  expect_identical(
    do.call(sim_horizon_paths, a)$trials,
    do.call(sim_horizon_paths, b)$trials
  )
})

test_that("horizon-specific risk-difference margins exclude boundary nulls", {
  for (method in c("riskdiff-wald", "riskdiff-fm")) {
    for (margin in c(-1, 1)) {
      args <- horizon_test_args(
        method = method,
        h0 = 0,
        horizon_settings = data.frame(horizon = c(12, 18), h0 = c(0, margin))
      )
      for (simulate in list(sim_horizon_paths, sim_horizon_trials)) {
        expect_error(
          do.call(simulate, args),
          "'h0' must lie in (-1, 1)",
          fixed = TRUE
        )
      }
    }
  }
})

test_that("invalid horizon settings fail before simulation", {
  for (override in list(
    list(endpoint_schedule = data.frame(N = c(40, 60), horizon = c(18, 12))),
    list(endpoint_schedule = data.frame(N = c(40, 60), horizon = c(0, 18))),
    list(endpoint_schedule = data.frame(N = c(40, 60), horizon = c(12, Inf))),
    list(endpoint_schedule = data.frame(N = c(40, 60), horizon = c(NA, 18))),
    list(endpoint_schedule = data.frame(N = c(40, 60), horizon = c("12", "18"))),
    list(endpoint_schedule = data.frame(N = c(60, 40), horizon = c(12, 18))),
    list(endpoint_schedule = data.frame(N = c(40, 40), horizon = c(12, 18))),
    list(endpoint_schedule = data.frame(N = c(40, 60), horizon = c(12, 18), x = 1)),
    list(endpoint_schedule = NULL),
    list(final_horizon = 12),
    list(dropout_horizon = 0),
    list(endpoint_schedule = data.frame(N = c(40, 80), horizon = c(12, 18))),
    list(cutpoints = 18),
    list(horizon_settings = data.frame(horizon = 12)),
    list(horizon_settings = data.frame(horizon = c(12, 18), h0 = c(0, NA))),
    list(
      horizon_settings = data.frame(horizon = c(12, 18), rmst_tau = c(18, 18))
    ),
    list(method = "riskdiff-fm", imputed_final = TRUE, prop_loss = 0.1)
  )) {
    expect_error(do.call(
      sim_horizon_paths,
      do.call(horizon_test_args, override)
    ))
  }
  bank <- do.call(sim_horizon_paths, horizon_test_args())
  for (schedule in list(
    data.frame(N = c(40, 60), horizon = c(18, 12)),
    data.frame(N = c(40, 60), horizon = c(6, 18)),
    data.frame(N = c(40, 80), horizon = c(12, 18))
  )) {
    expect_error(apply_horizon_rules(bank, schedule), "endpoint_schedule")
  }
  malformed <- bank
  malformed$trials[[1]]$finals$horizon[1] <- 18
  expect_error(
    apply_horizon_rules(malformed, horizon_test_schedule()),
    "Inconsistent"
  )
})

test_that("designs without interims use the final endpoint only", {
  schedule <- data.frame(N = numeric(), horizon = numeric())
  bank <- do.call(sim_horizon_paths, horizon_test_args(endpoint_schedule = schedule))
  out <- apply_horizon_rules(bank, schedule, return_trace = TRUE)
  expect_equal(out$sims$endpoint_horizon, c(18, 18))
  expect_equal(nrow(out$traces), 0)
  args <- horizon_test_args(endpoint_schedule = schedule)
  direct <- do.call(
    sim_horizon_trials,
    args
  )
  expect_identical(out$sims, direct$sims)
})
