test_that("impute_data updates success rows by index for two-arm data", {
  data_in <- data.frame(
    time = c(2, 4, 6, 8, 10, 12),
    treatment = c(1, 0, 1, 0, 1, 0),
    event = c(0, 0, 1, 0, 0, 1),
    id = 1:6,
    subject_enrolled = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    subject_impute_success = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
    subject_impute_futility = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
  )
  hazard <- array(c(0.04, 0.06), dim = c(1, 1, 2))
  treatment_idx <- data_in$treatment == 1 & data_in$subject_impute_success
  control_idx <- data_in$treatment == 0 & data_in$subject_impute_success
  no_impute_idx <- !data_in$subject_impute_success

  set.seed(3121)
  expected_treatment <- pwe_impute(
    time = data_in$time[treatment_idx],
    hazard = hazard[1, , 1],
    maxtime = 36,
    cutpoints = NULL
  )
  expected_control <- pwe_impute(
    time = data_in$time[control_idx],
    hazard = hazard[1, , 2],
    maxtime = 36,
    cutpoints = NULL
  )

  set.seed(3121)
  out <- goldilocks:::impute_data(
    data_in = data_in,
    hazard = hazard,
    end_of_study = 36,
    cutpoints = NULL,
    type = "success",
    single_arm = FALSE
  )

  expect_equal(nrow(out), nrow(data_in))
  expect_named(out, names(data_in))
  expect_equal(out$id, data_in$id)
  expect_equal(out$time[no_impute_idx], data_in$time[no_impute_idx])
  expect_equal(out$event[no_impute_idx], data_in$event[no_impute_idx])
  expect_equal(out$time[treatment_idx], expected_treatment$time)
  expect_equal(out$event[treatment_idx], expected_treatment$event)
  expect_equal(out$time[control_idx], expected_control$time)
  expect_equal(out$event[control_idx], expected_control$event)
})

test_that("impute_data requires exactly one posterior hazard draw", {
  data_in <- data.frame(
    time = 1,
    treatment = 1,
    event = 0,
    subject_impute_success = TRUE,
    subject_impute_futility = FALSE
  )

  expect_error(
    goldilocks:::impute_data(
      data_in = data_in,
      hazard = array(0.04, dim = c(2, 1, 1)),
      end_of_study = 36,
      cutpoints = NULL,
      type = "success",
      single_arm = TRUE
    ),
    "exactly one posterior draw"
  )
})

test_that("impute_data updates futility rows by index for two-arm data", {
  data_in <- data.frame(
    time = c(2, 4, 6, 8, 10, 12),
    treatment = c(1, 0, 1, 0, 1, 0),
    event = c(0, 0, 1, 0, 0, 1),
    id = 1:6,
    subject_enrolled = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    subject_impute_success = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
    subject_impute_futility = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
  )
  hazard <- array(c(0.04, 0.06), dim = c(1, 1, 2))
  treatment_idx <- data_in$treatment == 1 & data_in$subject_impute_futility
  control_idx <- data_in$treatment == 0 & data_in$subject_impute_futility
  no_impute_idx <- !data_in$subject_impute_futility

  set.seed(5128)
  expected_treatment <- pwe_sim(
    n = sum(treatment_idx),
    hazard = hazard[1, , 1],
    maxtime = 36,
    cutpoints = NULL
  )
  expected_control <- pwe_sim(
    n = sum(control_idx),
    hazard = hazard[1, , 2],
    maxtime = 36,
    cutpoints = NULL
  )

  set.seed(5128)
  out <- goldilocks:::impute_data(
    data_in = data_in,
    hazard = hazard,
    end_of_study = 36,
    cutpoints = NULL,
    type = "futility",
    single_arm = FALSE
  )

  expect_equal(nrow(out), nrow(data_in))
  expect_named(out, names(data_in))
  expect_equal(out$id, data_in$id)
  expect_equal(out$time[no_impute_idx], data_in$time[no_impute_idx])
  expect_equal(out$event[no_impute_idx], data_in$event[no_impute_idx])
  expect_equal(out$time[treatment_idx], expected_treatment$time)
  expect_equal(out$event[treatment_idx], expected_treatment$event)
  expect_equal(out$time[control_idx], expected_control$time)
  expect_equal(out$event[control_idx], expected_control$event)
})

test_that("impute_data updates success rows by index for single-arm data", {
  data_in <- data.frame(
    time = c(2, 4, 6, 8),
    treatment = rep(1, 4),
    event = c(0, 1, 0, 0),
    id = 1:4,
    subject_enrolled = c(TRUE, TRUE, TRUE, FALSE),
    subject_impute_success = c(TRUE, FALSE, TRUE, FALSE),
    subject_impute_futility = c(FALSE, FALSE, FALSE, TRUE)
  )
  hazard <- array(c(0.04, 0.06), dim = c(1, 2, 1))
  treatment_idx <- data_in$treatment == 1 & data_in$subject_impute_success
  no_impute_idx <- !data_in$subject_impute_success

  set.seed(7824)
  expected_treatment <- pwe_impute(
    time = data_in$time[treatment_idx],
    hazard = hazard[1, , 1],
    maxtime = 36,
    cutpoints = 12
  )

  set.seed(7824)
  out <- goldilocks:::impute_data(
    data_in = data_in,
    hazard = hazard,
    end_of_study = 36,
    cutpoints = 12,
    type = "success",
    single_arm = TRUE
  )

  expect_equal(nrow(out), nrow(data_in))
  expect_named(out, names(data_in))
  expect_equal(out$id, data_in$id)
  expect_equal(out$time[no_impute_idx], data_in$time[no_impute_idx])
  expect_equal(out$event[no_impute_idx], data_in$event[no_impute_idx])
  expect_equal(out$time[treatment_idx], expected_treatment$time)
  expect_equal(out$event[treatment_idx], expected_treatment$event)
})

test_that("impute_data updates futility rows by index for single-arm data", {
  data_in <- data.frame(
    time = c(2, 4, 6, 8),
    treatment = rep(1, 4),
    event = c(0, 1, 0, 0),
    id = 1:4,
    subject_enrolled = c(TRUE, TRUE, FALSE, FALSE),
    subject_impute_success = c(TRUE, FALSE, FALSE, FALSE),
    subject_impute_futility = c(FALSE, FALSE, TRUE, TRUE)
  )
  hazard <- array(c(0.04, 0.06), dim = c(1, 2, 1))
  treatment_idx <- data_in$treatment == 1 & data_in$subject_impute_futility
  no_impute_idx <- !data_in$subject_impute_futility

  set.seed(9184)
  expected_treatment <- pwe_sim(
    n = sum(treatment_idx),
    hazard = hazard[1, , 1],
    maxtime = 36,
    cutpoints = 12
  )

  set.seed(9184)
  out <- goldilocks:::impute_data(
    data_in = data_in,
    hazard = hazard,
    end_of_study = 36,
    cutpoints = 12,
    type = "futility",
    single_arm = TRUE
  )

  expect_equal(nrow(out), nrow(data_in))
  expect_named(out, names(data_in))
  expect_equal(out$id, data_in$id)
  expect_equal(out$time[no_impute_idx], data_in$time[no_impute_idx])
  expect_equal(out$event[no_impute_idx], data_in$event[no_impute_idx])
  expect_equal(out$time[treatment_idx], expected_treatment$time)
  expect_equal(out$event[treatment_idx], expected_treatment$event)
})

test_that("Bernoulli and event-time imputation agree within tolerance", {
  n <- 20000
  data_in <- data.frame(
    time = rep(seq(0, 23, length.out = 200), length.out = n),
    treatment = rep(1, n),
    event = rep(0, n),
    subject_impute_success = rep(TRUE, n),
    subject_impute_futility = rep(FALSE, n)
  )
  hazard <- array(c(0.02, 0.08, 0.01), dim = c(1, 3, 1))

  set.seed(2101)
  event_time <- goldilocks:::impute_data(
    data_in = data_in,
    hazard = hazard,
    end_of_study = 24,
    cutpoints = c(5, 12),
    type = "success",
    single_arm = TRUE,
    binary_imputation = "event-time"
  )
  set.seed(2101)
  bernoulli <- goldilocks:::impute_data(
    data_in = data_in,
    hazard = hazard,
    end_of_study = 24,
    cutpoints = c(5, 12),
    type = "success",
    single_arm = TRUE,
    binary_imputation = "bernoulli"
  )

  event_rate_tolerance <- 0.005
  expect_lt(
    abs(mean(event_time$event) - mean(bernoulli$event)),
    event_rate_tolerance
  )
  expect_true(all(bernoulli$time == 24))
})

test_that("batched event-time imputations match the scalar reference", {
  data_in <- data.frame(
    time = c(2, 4, 6, 8, 0, 0),
    treatment = c(1, 0, 1, 0, 1, 0),
    event = c(0, 0, 1, 0, 0, 0),
    id = seq_len(6),
    subject_enrolled = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    subject_impute_success = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
    subject_impute_futility = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
  )
  hazards <- array(
    c(
      0.04,
      0.05,
      0.03,
      0.04,
      0.05,
      0.06,
      0.07,
      0.05,
      0.06,
      0.04,
      0.05,
      0.07
    ),
    dim = c(3, 2, 2)
  )

  set.seed(4201)
  expected <- lapply(seq_len(dim(hazards)[1]), function(draw) {
    current <- goldilocks:::impute_data(
      data_in = data_in,
      hazard = hazards[draw, , , drop = FALSE],
      end_of_study = 36,
      cutpoints = 12,
      type = "success",
      single_arm = FALSE
    )
    maximum <- goldilocks:::impute_data(
      data_in = current,
      hazard = hazards[draw, , , drop = FALSE],
      end_of_study = 36,
      cutpoints = 12,
      type = "futility",
      single_arm = FALSE
    )
    list(current = current, maximum = maximum)
  })

  set.seed(4201)
  batched <- goldilocks:::impute_predictive_draws(
    data_in = data_in,
    hazards = hazards,
    end_of_study = 36,
    cutpoints = 12,
    single_arm = FALSE,
    binary_imputation = "event-time",
    check_futility = TRUE
  )

  expect_identical(batched$n_draws, 3L)
  expect_identical(batched$current$rows, c(1L, 2L, 4L))
  expect_identical(batched$future$rows, c(5L, 6L))
  for (draw in seq_len(batched$n_draws)) {
    current <- goldilocks:::complete_predictive_data(
      data_in,
      batched,
      draw
    )
    maximum <- goldilocks:::complete_predictive_data(
      data_in,
      batched,
      draw,
      include_future = TRUE
    )
    expect_equal(current, expected[[draw]]$current, tolerance = 1e-12)
    expect_equal(maximum, expected[[draw]]$maximum, tolerance = 1e-12)
    expect_true(all(
      current$time[batched$current$rows] >= data_in$time[batched$current$rows]
    ))
  }
})

test_that("batched Bernoulli imputations match the scalar reference", {
  data_in <- data.frame(
    time = c(2, 4, 6, 8, 0, 0),
    treatment = c(1, 0, 1, 0, 1, 0),
    event = c(0, 0, 1, 0, 0, 0),
    subject_enrolled = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    subject_impute_success = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
    subject_impute_futility = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
  )
  hazards <- array(
    c(
      0.04,
      0.05,
      0.03,
      0.04,
      0.05,
      0.06,
      0.07,
      0.05,
      0.06,
      0.04,
      0.05,
      0.07
    ),
    dim = c(3, 2, 2)
  )

  set.seed(4202)
  expected <- lapply(seq_len(dim(hazards)[1]), function(draw) {
    current <- goldilocks:::impute_data(
      data_in = data_in,
      hazard = hazards[draw, , , drop = FALSE],
      end_of_study = 36,
      cutpoints = 12,
      type = "success",
      single_arm = FALSE,
      binary_imputation = "bernoulli"
    )
    maximum <- goldilocks:::impute_data(
      data_in = current,
      hazard = hazards[draw, , , drop = FALSE],
      end_of_study = 36,
      cutpoints = 12,
      type = "futility",
      single_arm = FALSE,
      binary_imputation = "bernoulli"
    )
    list(current = current, maximum = maximum)
  })

  set.seed(4202)
  batched <- goldilocks:::impute_predictive_draws(
    data_in = data_in,
    hazards = hazards,
    end_of_study = 36,
    cutpoints = 12,
    single_arm = FALSE,
    binary_imputation = "bernoulli",
    check_futility = TRUE
  )

  for (draw in seq_len(batched$n_draws)) {
    current <- goldilocks:::complete_predictive_data(
      data_in,
      batched,
      draw
    )
    maximum <- goldilocks:::complete_predictive_data(
      data_in,
      batched,
      draw,
      include_future = TRUE
    )
    expect_identical(current, expected[[draw]]$current)
    expect_identical(maximum, expected[[draw]]$maximum)
  }
})

test_that("batched imputation handles single-arm and empty blocks", {
  data_in <- data.frame(
    time = c(4, 12),
    treatment = c(1L, 1L),
    event = c(1, 0),
    subject_enrolled = c(TRUE, TRUE),
    subject_impute_success = c(FALSE, FALSE),
    subject_impute_futility = c(FALSE, FALSE)
  )
  hazards <- array(NA_real_, dim = c(4, 2, 2))
  hazards[,, 1] <- matrix(
    c(0.04, 0.05, 0.06, 0.07, 0.03, 0.04, 0.05, 0.06),
    nrow = 4
  )

  set.seed(4203)
  before <- .Random.seed
  out <- goldilocks:::impute_predictive_draws(
    data_in = data_in,
    hazards = hazards,
    end_of_study = 24,
    cutpoints = 12,
    single_arm = TRUE,
    check_futility = FALSE
  )

  expect_identical(dim(out$current$time), c(0L, 4L))
  expect_identical(dim(out$current$event), c(0L, 4L))
  expect_null(out$future)
  expect_identical(.Random.seed, before)
  expect_identical(
    goldilocks:::complete_predictive_data(data_in, out, 2),
    data_in
  )
  expect_error(
    goldilocks:::complete_predictive_data(
      data_in,
      out,
      2,
      include_future = TRUE
    ),
    "future draws unavailable"
  )
})

test_that("batched constant-hazard single-arm draws match scalar imputation", {
  data_in <- data.frame(
    time = c(3, 8, 0),
    treatment = rep(1L, 3),
    event = c(0, 1, 0),
    subject_enrolled = c(TRUE, TRUE, FALSE),
    subject_impute_success = c(TRUE, FALSE, FALSE),
    subject_impute_futility = c(FALSE, FALSE, TRUE)
  )
  hazards <- array(NA_real_, dim = c(4, 1, 2))
  hazards[, 1, 1] <- c(0.02, 0.04, 0.06, 0.08)

  set.seed(4204)
  expected <- lapply(seq_len(dim(hazards)[1]), function(draw) {
    current <- goldilocks:::impute_data(
      data_in,
      hazards[draw, , , drop = FALSE],
      12,
      NULL,
      "success",
      TRUE
    )
    maximum <- goldilocks:::impute_data(
      current,
      hazards[draw, , , drop = FALSE],
      12,
      NULL,
      "futility",
      TRUE
    )
    list(current = current, maximum = maximum)
  })

  set.seed(4204)
  batched <- goldilocks:::impute_predictive_draws(
    data_in,
    hazards,
    12,
    NULL,
    TRUE,
    check_futility = TRUE
  )

  for (draw in seq_len(batched$n_draws)) {
    expect_equal(
      goldilocks:::complete_predictive_data(data_in, batched, draw),
      expected[[draw]]$current,
      tolerance = 1e-12
    )
    expect_equal(
      goldilocks:::complete_predictive_data(
        data_in,
        batched,
        draw,
        include_future = TRUE
      ),
      expected[[draw]]$maximum,
      tolerance = 1e-12
    )
  }
})
