test_that("sim_comp_data returns data frame with correct columns (two-arm)", {
  set.seed(5927)
  out <- sim_comp_data(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = 0,
    N_total = 50,
    lambda = 5,
    lambda_time = 0,
    end_of_study = 36
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 50)
  expect_named(
    out,
    c("time", "treatment", "event", "enrollment", "id", "loss_to_fu")
  )
})

test_that("sim_comp_data returns correct columns for single-arm", {
  set.seed(8341)
  out <- sim_comp_data(
    hazard_treatment = 0.01,
    hazard_control = NULL,
    cutpoints = 0,
    N_total = 30,
    lambda = 5,
    lambda_time = 0,
    end_of_study = 36
  )
  expect_equal(nrow(out), 30)
  expect_true(all(out$treatment == 1))
})

test_that("sim_comp_data applies loss to follow-up", {
  set.seed(2764)
  out <- sim_comp_data(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = 0,
    N_total = 200,
    lambda = 20,
    lambda_time = 0,
    end_of_study = 36,
    prop_loss = 0.30
  )
  n_lost <- sum(out$loss_to_fu)
  expect_equal(n_lost, ceiling(0.30 * 200))
  # Subjects lost to follow-up should be censored
  expect_true(all(out$event[out$loss_to_fu] == 0))
})

test_that("sim_comp_data has no loss to follow-up by default", {
  set.seed(4192)
  out <- sim_comp_data(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = 0,
    N_total = 50,
    lambda = 5,
    lambda_time = 0,
    end_of_study = 36
  )
  expect_true(all(!out$loss_to_fu))
})

test_that("sim_comp_data works with piecewise hazard", {
  set.seed(6083)
  out <- sim_comp_data(
    hazard_treatment = c(0.005, 0.01),
    hazard_control = c(0.01, 0.02),
    cutpoints = c(0, 12),
    N_total = 80,
    lambda = 10,
    lambda_time = 0,
    end_of_study = 36
  )
  expect_equal(nrow(out), 80)
  expect_true(all(out$time > 0))
})

test_that("sim_comp_data validates counts and loss-to-follow-up probability", {
  expect_error(
    sim_comp_data(
      hazard_treatment = 0.01,
      hazard_control = 0.02,
      cutpoints = 0,
      N_total = 50.5,
      lambda = 5,
      lambda_time = 0,
      end_of_study = 36
    ),
    "N_total"
  )

  expect_error(
    sim_comp_data(
      hazard_treatment = 0.01,
      hazard_control = 0.02,
      cutpoints = 0,
      N_total = 50,
      lambda = 5,
      lambda_time = 0,
      end_of_study = 36,
      prop_loss = -0.1
    ),
    "prop_loss"
  )
})
