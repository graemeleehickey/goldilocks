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
