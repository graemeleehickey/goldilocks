binary_prediction_fixture <- function() {
  data.frame(
    time = c(4, 12, 3, 5, 0, 0, 0, 0),
    event = c(1L, 0L, 0L, 0L, 1L, 1L, 0L, 0L),
    treatment = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L),
    subject_enrolled = c(rep(TRUE, 4), rep(FALSE, 4)),
    subject_impute_success = c(
      FALSE,
      FALSE,
      TRUE,
      TRUE,
      FALSE,
      FALSE,
      FALSE,
      FALSE
    ),
    subject_impute_futility = c(rep(FALSE, 4), rep(TRUE, 4))
  )
}

binary_prediction_imputations <- function() {
  list(
    n_draws = 4L,
    current = list(
      rows = c(3L, 4L),
      time = matrix(12, nrow = 2L, ncol = 4L),
      event = matrix(
        c(
          0L,
          1L,
          1L,
          1L,
          0L,
          1L,
          1L,
          1L
        ),
        nrow = 2L
      )
    ),
    future = list(
      rows = 5:8,
      time = matrix(12, nrow = 4L, ncol = 4L),
      event = matrix(
        c(
          0L,
          1L,
          0L,
          1L,
          1L,
          0L,
          0L,
          1L,
          0L,
          1L,
          0L,
          1L,
          1L,
          0L,
          0L,
          1L
        ),
        nrow = 4L
      )
    )
  )
}

data_frame_binary_reference <- function(
  data,
  imputations,
  method,
  bin_method,
  N_mcmc,
  prob_ha = 0.6
) {
  current_successes <- 0L
  maximum_successes <- 0L
  uncertain_now <- 0L
  uncertain_max <- 0L
  for (draw in seq_len(imputations$n_draws)) {
    current <- goldilocks:::complete_predictive_data(
      data_in = data,
      imputations = imputations,
      draw = draw
    )
    current <- current[
      current$subject_enrolled,
      c("time", "event", "treatment"),
      drop = FALSE
    ]
    maximum <- goldilocks:::complete_predictive_data(
      data_in = data,
      imputations = imputations,
      draw = draw,
      include_future = TRUE
    )
    maximum <- maximum[, c("time", "event", "treatment"), drop = FALSE]
    success_now <- goldilocks:::analyse_data(
      data = current,
      end_of_study = 12,
      cutpoints = NULL,
      single_arm = FALSE,
      prior_surv = c(0.1, 0.1),
      N_mcmc = N_mcmc,
      method = method,
      alternative = "greater",
      h0 = 0,
      prior_bin = c(1, 1),
      bin_method = bin_method,
      empty_interval = "prior"
    )
    success_max <- goldilocks:::analyse_data(
      data = maximum,
      end_of_study = 12,
      cutpoints = NULL,
      single_arm = FALSE,
      prior_surv = c(0.1, 0.1),
      N_mcmc = N_mcmc,
      method = method,
      alternative = "greater",
      h0 = 0,
      prior_bin = c(1, 1),
      bin_method = bin_method,
      empty_interval = "prior"
    )
    classification_now <- goldilocks:::classify_completed_analysis(
      success_now,
      prob_ha = prob_ha,
      mc_conf_level = 0.95
    )
    classification_max <- goldilocks:::classify_completed_analysis(
      success_max,
      prob_ha = prob_ha,
      mc_conf_level = 0.95
    )
    current_successes <- current_successes +
      as.integer(classification_now$crossed)
    maximum_successes <- maximum_successes +
      as.integer(classification_max$crossed)
    uncertain_now <- uncertain_now +
      as.integer(classification_now$uncertain)
    uncertain_max <- uncertain_max +
      as.integer(classification_max$uncertain)
  }
  list(
    current_successes = current_successes,
    maximum_successes = maximum_successes,
    inner_mc_uncertain_now = uncertain_now,
    inner_mc_uncertain_max = uncertain_max
  )
}

test_that("predictive binary counts match completed patient data", {
  data <- binary_prediction_fixture()
  imputations <- binary_prediction_imputations()
  counts <- goldilocks:::predictive_binary_counts(
    data,
    imputations,
    single_arm = FALSE,
    check_futility = TRUE
  )

  for (draw in seq_len(imputations$n_draws)) {
    current <- goldilocks:::complete_predictive_data(
      data,
      imputations,
      draw
    )
    current <- current[current$subject_enrolled, ]
    maximum <- goldilocks:::complete_predictive_data(
      data,
      imputations,
      draw,
      include_future = TRUE
    )
    expected_current <- data.frame(
      events_control = sum(current$event[current$treatment == 0L]),
      n_control = sum(current$treatment == 0L),
      events_treatment = sum(current$event[current$treatment == 1L]),
      n_treatment = sum(current$treatment == 1L)
    )
    expected_maximum <- data.frame(
      events_control = sum(maximum$event[maximum$treatment == 0L]),
      n_control = sum(maximum$treatment == 0L),
      events_treatment = sum(maximum$event[maximum$treatment == 1L]),
      n_treatment = sum(maximum$treatment == 1L)
    )

    actual_current <- counts$current[draw, , drop = FALSE]
    actual_maximum <- counts$maximum[draw, , drop = FALSE]
    rownames(actual_current) <- NULL
    rownames(actual_maximum) <- NULL
    expect_identical(actual_current, expected_current)
    expect_identical(actual_maximum, expected_maximum)
  }
})

test_that("count analyses match general analyses at sparse and boundary counts", {
  cases <- data.frame(
    events_control = c(0L, 0L, 1L, 9L, 10L),
    events_treatment = c(0L, 10L, 9L, 1L, 10L),
    n = 10L
  )

  for (row in seq_len(nrow(cases))) {
    events <- c(
      rep.int(1L, cases$events_control[row]),
      rep.int(0L, cases$n[row] - cases$events_control[row]),
      rep.int(1L, cases$events_treatment[row]),
      rep.int(0L, cases$n[row] - cases$events_treatment[row])
    )
    data <- data.frame(
      event = events,
      treatment = rep(0:1, each = cases$n[row])
    )

    for (bin_method in c("normal", "quadrature", "mc")) {
      set.seed(9100 + row)
      checked <- goldilocks:::bayes_binomial_test(
        data = data,
        single_arm = FALSE,
        alternative = "greater",
        h0 = 0.05,
        prior_bin = c(0.5, 0.5),
        bin_method = bin_method,
        N_mcmc = 50L
      )
      checked_seed <- .Random.seed
      set.seed(9100 + row)
      counted <- goldilocks:::bayes_binomial_from_counts(
        events_control = cases$events_control[row],
        n_control = cases$n[row],
        events_treatment = cases$events_treatment[row],
        n_treatment = cases$n[row],
        single_arm = FALSE,
        alternative = "greater",
        h0 = 0.05,
        prior_bin = c(0.5, 0.5),
        bin_method = bin_method,
        N_mcmc = 50L
      )
      counted_seed <- .Random.seed

      expect_identical(counted, checked)
      expect_identical(counted_seed, checked_seed)
    }
  }

  estimable <- cases[c(3L, 4L), ]
  for (row in seq_len(nrow(estimable))) {
    data <- data.frame(
      time = rep(12, 20),
      event = c(
        rep.int(1L, estimable$events_control[row]),
        rep.int(0L, 10L - estimable$events_control[row]),
        rep.int(1L, estimable$events_treatment[row]),
        rep.int(0L, 10L - estimable$events_treatment[row])
      ),
      treatment = rep(0:1, each = 10L)
    )
    checked <- goldilocks:::risk_difference_wald_test_checked(
      data,
      end_of_study = 12,
      alternative = "two.sided",
      h0 = 0.1
    )
    counted <- goldilocks:::risk_difference_wald_from_counts(
      events_control = estimable$events_control[row],
      n_control = 10L,
      events_treatment = estimable$events_treatment[row],
      n_treatment = 10L,
      alternative = "two.sided",
      h0 = 0.1
    )
    expect_identical(counted, checked)
  }
})

test_that("binary count prediction matches the patient-data reference", {
  data <- binary_prediction_fixture()
  imputations <- binary_prediction_imputations()
  counts <- goldilocks:::predictive_binary_counts(
    data,
    imputations,
    single_arm = FALSE,
    check_futility = TRUE
  )
  configurations <- list(
    list(method = "riskdiff", bin_method = "mc", N_mcmc = 1L),
    list(method = "bayes-bin", bin_method = "normal", N_mcmc = 1L),
    list(method = "bayes-bin", bin_method = "quadrature", N_mcmc = 1L),
    list(method = "bayes-bin", bin_method = "mc", N_mcmc = 40L)
  )

  for (configuration in configurations) {
    set.seed(9201)
    expected <- data_frame_binary_reference(
      data = data,
      imputations = imputations,
      method = configuration$method,
      bin_method = configuration$bin_method,
      N_mcmc = configuration$N_mcmc
    )
    expected_seed <- .Random.seed

    set.seed(9201)
    actual <- goldilocks:::analyse_predictive_binary_counts(
      completed_counts = counts,
      single_arm = FALSE,
      N_mcmc = configuration$N_mcmc,
      method = configuration$method,
      alternative = "greater",
      h0 = 0,
      prior_bin = c(1, 1),
      bin_method = configuration$bin_method,
      check_futility = TRUE,
      prob_ha = 0.6,
      mc_conf_level = 0.95
    )
    actual_seed <- .Random.seed

    expect_identical(
      actual[names(expected)],
      expected,
      info = paste(configuration$method, configuration$bin_method)
    )
    expect_identical(
      actual_seed,
      expected_seed,
      info = paste(configuration$method, configuration$bin_method)
    )
    expect_gt(actual$reuse$repeated_count_summaries, 0)
    if (
      configuration$method == "bayes-bin" &&
        configuration$bin_method == "mc"
    ) {
      expect_false(actual$reuse$reuse_enabled)
      expect_identical(actual$reuse$reused_analyses, 0L)
    } else {
      expect_true(actual$reuse$reuse_enabled)
      expect_gt(actual$reuse$reused_analyses, 0)
    }
  }
})

test_that("binary analysis groups contain every analysis input", {
  base <- list(
    counts = data.frame(
      events_control = 2L,
      n_control = 10L,
      events_treatment = 4L,
      n_treatment = 10L
    ),
    method = "bayes-bin",
    single_arm = FALSE,
    alternative = "greater",
    h0 = 0,
    prior_bin = c(1, 1),
    bin_method = "quadrature",
    N_mcmc = 100L
  )
  variants <- list(
    within(base, counts$events_control <- 3L),
    within(base, counts$n_control <- 11L),
    within(base, counts$events_treatment <- 5L),
    within(base, counts$n_treatment <- 11L),
    within(base, method <- "riskdiff"),
    within(base, single_arm <- TRUE),
    within(base, alternative <- "less"),
    within(base, h0 <- 0.1),
    within(base, prior_bin <- c(0.5, 0.5)),
    within(base, bin_method <- "normal"),
    within(base, N_mcmc <- 101L)
  )
  keys <- vapply(
    c(list(base), variants),
    function(arguments) {
      do.call(goldilocks:::group_binary_analyses, arguments)
    },
    character(1L)
  )

  expect_length(unique(keys), length(keys))
})

test_that("binary interim diagnostics report bounded summary reuse", {
  data <- data.frame(
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
    )
  )
  out <- evaluate_interim(
    data = data,
    data_cut = 8,
    look = 1,
    N_total = 12,
    end_of_study = 10,
    method = "bayes-bin",
    bin_method = "quadrature",
    binary_imputation = "bernoulli",
    alternative = "less",
    N_impute = 20,
    N_mcmc = 1,
    seed = 9301
  )
  reuse <- out$diagnostics$binary_count_reuse

  expect_true(reuse$reuse_enabled)
  expect_identical(reuse$analysis_requests, 40L)
  expect_identical(
    reuse$analysis_requests,
    reuse$unique_count_summaries + reuse$reused_analyses
  )
  expect_equal(
    reuse$reused_analysis_rate,
    reuse$reused_analyses / reuse$analysis_requests
  )
})
