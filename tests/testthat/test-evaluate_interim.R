interim_fixture <- function() {
  data.frame(
    id = 1:8,
    treatment = c(0, 1, 0, 1, 0, 1, 0, 1),
    enrollment = 0:7,
    time = c(8, 7, 6, 5, 4, 3, 2, 1),
    event = c(1, 0, 0, 1, 0, 0, 0, 0),
    status = c(
      "event",
      "pending",
      "censored",
      "event",
      "pending",
      "pending",
      "pending",
      "pending"
    ),
    stringsAsFactors = FALSE
  )
}

interim_args <- function() {
  list(
    data = interim_fixture(),
    data_cut = 8,
    look = 2,
    N_total = 12,
    end_of_study = 10,
    rand_ratio = c(control = 1, treatment = 1),
    alternative = "less",
    N_impute = 3,
    N_mcmc = 3,
    seed = 8401
  )
}

test_that("evaluate_interim returns an auditable one-look result", {
  input <- interim_fixture()
  before <- input
  args <- interim_args()
  args$data <- input

  out <- do.call(evaluate_interim, args)

  expect_s3_class(out, "goldilocks_interim")
  expect_identical(input, before)
  expect_named(
    out,
    c(
      "decision",
      "probabilities",
      "monte_carlo",
      "diagnostics",
      "trace",
      "metadata"
    )
  )
  expect_equal(nrow(out$decision), 1)
  expect_equal(nrow(out$trace), 1)
  expect_equal(out$trace$look, 2)
  expect_equal(out$trace$planned_N, 8)
  expect_equal(out$trace$calendar_time, 8)
  expect_equal(out$trace$active_followup, 5)
  expect_equal(out$trace$N_pending, 6)
  expect_equal(out$trace$N_not_enrolled, 4)
  expect_identical(
    out$diagnostics$participant_counts$status,
    c("event", "complete", "pending", "censored", "not_enrolled")
  )
  expect_identical(
    out$diagnostics$participant_counts$count,
    c(2L, 0L, 5L, 1L, 4L)
  )
  expect_identical(
    out$diagnostics$target_allocation,
    c(control = 6L, treatment = 6L)
  )
  expect_identical(
    out$diagnostics$current_allocation,
    c(control = 4L, treatment = 4L)
  )
  expect_identical(
    out$diagnostics$potential_accruals,
    c(control = 2L, treatment = 2L)
  )
  expect_identical(
    out$metadata$design$rand_ratio,
    c(control = 1, treatment = 1)
  )
  expect_identical(out$metadata$rng$seed, 8401)
  expect_identical(
    out$metadata$time_origin,
    "first participant randomization (time 0)"
  )
  expect_false("data" %in% names(out$metadata))
  grDevices::pdf(tempfile(fileext = ".pdf"))
  expect_identical(plot_trial_trace(out), out$trace)
  grDevices::dev.off()
  expect_equal(summarise_trial_trace(out)$last_look, 2)
})

test_that("evaluate_interim derives future accrual from rand_ratio", {
  data <- interim_fixture()[1:6, ]
  data$treatment <- c(0, 1, 1, 0, 1, 1)

  out <- evaluate_interim(
    data = data,
    data_cut = 8,
    look = 1,
    N_total = 12,
    end_of_study = 10,
    rand_ratio = c(treatment = 2, control = 1),
    alternative = "less",
    N_impute = 2,
    seed = 8402
  )

  expect_identical(
    out$metadata$design$rand_ratio,
    c(control = 1, treatment = 2)
  )
  expect_identical(
    out$diagnostics$target_allocation,
    c(control = 4L, treatment = 8L)
  )
  expect_identical(
    out$diagnostics$current_allocation,
    c(control = 2L, treatment = 4L)
  )
  expect_identical(
    out$diagnostics$potential_accruals,
    c(control = 2L, treatment = 4L)
  )
  expect_equal(out$trace$N_not_enrolled, 6)
})

test_that("evaluate_interim resolves arm-specific hazard priors", {
  args <- interim_args()
  args$cutpoints <- 4
  args$prior_surv <- list(
    treatment = rbind(shape = c(5, 6), rate = c(10, 12)),
    control = c(2, 4)
  )

  out <- do.call(evaluate_interim, args)

  expect_identical(
    out$diagnostics$prior$arm,
    rep(c("control", "treatment"), each = 2)
  )
  expect_equal(out$diagnostics$prior$shape, c(2, 2, 5, 6))
  expect_identical(out$metadata$prior_design, out$diagnostics$prior)
  expect_equal(
    out$diagnostics$posterior$posterior_shape,
    out$diagnostics$posterior$prior_shape +
      out$diagnostics$posterior$effective_events
  )
  expect_equal(
    out$diagnostics$posterior$posterior_rate,
    out$diagnostics$posterior$prior_rate +
      out$diagnostics$posterior$effective_exposure
  )
})

test_that("evaluate_interim supports single-arm Bayesian analyses", {
  data <- interim_fixture()[1:4, ]
  data$treatment <- 1

  out <- evaluate_interim(
    data = data,
    data_cut = 8,
    look = 1,
    N_total = 7,
    end_of_study = 10,
    single_arm = TRUE,
    method = "bayes-bin",
    binary_imputation = "bernoulli",
    N_impute = 2,
    N_mcmc = 2,
    seed = 8403
  )

  expect_identical(out$diagnostics$target_allocation, c(treatment = 7L))
  expect_identical(out$diagnostics$potential_accruals, c(treatment = 3L))
  expect_equal(out$trace$N_control, 0)
  expect_equal(out$trace$N_not_enrolled, 3)
})

test_that("evaluate_interim supports every two-arm analysis method", {
  methods <- c(
    "logrank",
    "cox",
    "riskdiff-wald",
    "riskdiff-fm",
    "bayes-surv",
    "bayes-bin"
  )
  results <- lapply(seq_along(methods), function(index) {
    evaluate_interim(
      data = interim_fixture(),
      data_cut = 8,
      look = 1,
      N_total = 12,
      end_of_study = 10,
      rand_ratio = c(control = 1, treatment = 1),
      alternative = "less",
      method = methods[[index]],
      N_impute = 2,
      N_mcmc = 3,
      seed = 8410 + index
    )
  })

  expect_true(all(vapply(results, inherits, logical(1), "goldilocks_interim")))
  expect_true(all(vapply(
    results,
    function(result) nrow(result$trace) == 1L,
    logical(1)
  )))
})

test_that("evaluate_interim maps deprecated riskdiff to riskdiff-wald", {
  expect_warning(
    result <- evaluate_interim(
      data = interim_fixture(),
      data_cut = 8,
      look = 1,
      N_total = 12,
      end_of_study = 10,
      rand_ratio = c(control = 1, treatment = 1),
      alternative = "less",
      method = "riskdiff",
      N_impute = 2,
      N_mcmc = 1,
      seed = 8463
    ),
    "deprecated; use `method = \"riskdiff-wald\"`"
  )

  expect_identical(result$metadata$design$method, "riskdiff-wald")
})

test_that("Fn zero disables the maximum-sample calculation", {
  args <- interim_args()
  args$Fn <- 0
  out <- do.call(evaluate_interim, args)

  expect_equal(nrow(out$probabilities), 1)
  expect_equal(nrow(out$monte_carlo), 1)
  expect_true(is.na(out$trace$ppp_success_at_max))
  expect_true(is.na(out$trace$futility_threshold))
})

test_that("evaluate_interim permits an already reached maximum sample size", {
  out <- evaluate_interim(
    data = interim_fixture(),
    data_cut = 8,
    look = 3,
    N_total = 8,
    end_of_study = 10,
    alternative = "less",
    N_impute = 2,
    seed = 8409
  )

  expect_identical(
    out$diagnostics$potential_accruals,
    c(control = 0L, treatment = 0L)
  )
  expect_equal(out$trace$N_not_enrolled, 0)
  expect_true(is.finite(out$trace$ppp_success_at_max))
})

test_that("a supplied seed is reproducible and preserves caller RNG state", {
  set.seed(8404)
  before_kind <- RNGkind()
  before_seed <- .Random.seed

  first <- do.call(evaluate_interim, interim_args())
  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_seed)

  second <- do.call(evaluate_interim, interim_args())
  expect_identical(first, second)
})

test_that("external and simulated interim paths use the same calculation", {
  set.seed(8405)
  N_total <- 20
  look_N <- 10
  end_of_study <- 12
  total <- sim_comp_data(
    hazard_treatment = -log(0.8) / end_of_study,
    hazard_control = -log(0.65) / end_of_study,
    N_total = N_total,
    lambda = 2,
    end_of_study = end_of_study,
    block = 2,
    rand_ratio = c(control = 1, treatment = 1),
    prop_loss = 0
  )
  look_time <- total$enrollment[[look_N]]
  simulated <- within(total, {
    subject_enrolled <- id <= look_N
    subject_impute_futility <- !subject_enrolled
    time_from_rand_at_look <- look_time - enrollment
    subject_impute_success <-
      ((event == 1) * (time_from_rand_at_look < time) & subject_enrolled) |
      ((event == 0) *
        (time_from_rand_at_look < end_of_study) &
        subject_enrolled) |
      (loss_to_fu & subject_enrolled)
  })
  simulated <- within(simulated, {
    time <- pmax(pmin(time, time_from_rand_at_look), .Machine$double.eps)
    event <- ifelse(subject_impute_success, 0, event)
  })

  observed <- simulated[simulated$subject_enrolled, ]
  external <- data.frame(
    id = observed$id,
    treatment = observed$treatment,
    enrollment = observed$enrollment,
    time = observed$time,
    event = observed$event,
    status = ifelse(
      observed$event == 1,
      "event",
      ifelse(observed$subject_impute_success, "pending", "complete")
    ),
    stringsAsFactors = FALSE
  )

  set.seed(8406)
  expected <- goldilocks:::evaluate_interim_decision(
    data_interim = simulated,
    look = 1,
    planned_N = look_N,
    calendar_time = look_time,
    active_followup = goldilocks:::active_followup_at(total, look_time),
    end_of_study = end_of_study,
    cutpoints = NULL,
    single_arm = FALSE,
    prior_surv = c(0.1, 0.1),
    prior_bin = c(1, 1),
    bin_method = "mc",
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.95,
    N_impute = 4,
    N_mcmc = 1,
    mc_conf_level = 0.95,
    empty_interval = "prior",
    method = "logrank",
    binary_imputation = "event-time",
    check_futility = TRUE
  )
  set.seed(8406)
  actual <- evaluate_interim(
    data = external,
    data_cut = look_time,
    look = 1,
    N_total = N_total,
    end_of_study = end_of_study,
    rand_ratio = c(control = 1, treatment = 1),
    alternative = "less",
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.95,
    N_impute = 4,
    N_mcmc = 1,
    method = "logrank"
  )

  expect_identical(actual$trace, expected$trace)
  expect_identical(actual$monte_carlo, expected$monte_carlo)
})

test_that("evaluate_interim validates allocation assumptions", {
  args <- interim_args()
  args$N_total <- 10
  args$rand_ratio <- c(control = 1, treatment = 2)
  expect_error(
    do.call(evaluate_interim, args),
    "must divide exactly according to 'rand_ratio'"
  )

  data <- interim_fixture()[1:6, ]
  data$treatment <- c(0, 0, 0, 0, 1, 1)
  expect_error(
    evaluate_interim(
      data,
      data_cut = 8,
      look = 1,
      N_total = 6,
      end_of_study = 10,
      N_impute = 1
    ),
    "Observed enrollment exceeds"
  )

  data$treatment <- 1
  expect_error(
    evaluate_interim(
      data,
      data_cut = 8,
      look = 1,
      N_total = 8,
      end_of_study = 10,
      N_impute = 1
    ),
    "must contain subjects in both arms"
  )
})

test_that("evaluate_interim validates the observed data schema", {
  args <- interim_args()

  invalid <- args
  invalid$data$status <- NULL
  expect_error(do.call(evaluate_interim, invalid), "containing columns")

  invalid <- args
  invalid$data$id[[2]] <- invalid$data$id[[1]]
  expect_error(do.call(evaluate_interim, invalid), "unique, non-missing")

  invalid <- args
  invalid$data$enrollment <- invalid$data$enrollment + 1
  expect_error(do.call(evaluate_interim, invalid), "time 0")

  invalid <- args
  invalid$data$enrollment[[8]] <- 9
  expect_error(do.call(evaluate_interim, invalid), "after 'data_cut'")

  invalid <- args
  invalid$data$time[[8]] <- 2
  expect_error(do.call(evaluate_interim, invalid), "follow-up available")

  invalid <- args
  invalid$data$event[[2]] <- 1
  expect_error(do.call(evaluate_interim, invalid), "Only rows with")

  invalid <- args
  invalid$data$status[[2]] <- "complete"
  expect_error(do.call(evaluate_interim, invalid), "time = end_of_study")

  invalid <- args
  invalid$seed <- -1
  expect_error(do.call(evaluate_interim, invalid), "non-negative integer")
})
