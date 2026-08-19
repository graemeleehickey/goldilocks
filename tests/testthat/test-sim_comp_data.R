test_that("sim_comp_data returns data frame with correct columns (two-arm)", {
  set.seed(5927)
  out <- sim_comp_data(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = NULL,
    N_total = 50,
    lambda = 5,
    lambda_time = NULL,
    end_of_study = 36
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 50)
  expect_named(
    out,
    c("time", "treatment", "event", "enrollment", "id", "loss_to_fu")
  )
  expect_equal(out$enrollment[1], 0)
  expect_true(all(diff(out$enrollment) > 0))
  expect_true(any(out$enrollment[-1] != floor(out$enrollment[-1])))
})

test_that("sim_comp_data accepts internal enrollment-rate knots", {
  set.seed(5419)
  out <- sim_comp_data(
    hazard_treatment = 0.01,
    N_total = 40,
    lambda = c(2, 5),
    lambda_time = 2.5,
    end_of_study = 36
  )

  expect_equal(out$enrollment[1], 0)
  expect_true(all(diff(out$enrollment) > 0))
})

test_that("sim_comp_data returns correct columns for single-arm", {
  set.seed(8341)
  out <- sim_comp_data(
    hazard_treatment = 0.01,
    hazard_control = NULL,
    cutpoints = NULL,
    N_total = 30,
    lambda = 5,
    lambda_time = NULL,
    end_of_study = 36,
    prop_loss = 0.20
  )
  expect_equal(nrow(out), 30)
  expect_true(all(out$treatment == 1))
  expect_equal(sum(out$loss_to_fu), ceiling(0.20 * 30))
})

test_that("sim_comp_data applies loss to follow-up", {
  set.seed(2764)
  out <- sim_comp_data(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = NULL,
    N_total = 200,
    lambda = 20,
    lambda_time = NULL,
    end_of_study = 36,
    prop_loss = 0.30
  )
  n_lost <- sum(out$loss_to_fu)
  expect_equal(n_lost, ceiling(0.30 * 200))
  # Subjects lost to follow-up should be censored
  expect_true(all(out$event[out$loss_to_fu] == 0))
})

test_that("scalar and equal arm-specific loss proportions are identical", {
  common_args <- list(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = NULL,
    N_total = 99,
    lambda = 20,
    lambda_time = NULL,
    end_of_study = 36,
    block = 3,
    rand_ratio = c(1, 2)
  )

  set.seed(2061)
  scalar <- do.call(
    sim_comp_data,
    c(common_args, list(prop_loss = 0.20))
  )
  set.seed(2061)
  arm_specific <- do.call(
    sim_comp_data,
    c(
      common_args,
      list(prop_loss = c(treatment = 0.20, control = 0.20))
    )
  )

  expect_identical(arm_specific, scalar)
})

test_that("sim_comp_data applies differential loss within randomized arms", {
  common_args <- list(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = NULL,
    N_total = 99,
    lambda = 20,
    lambda_time = NULL,
    end_of_study = 36,
    block = 3,
    rand_ratio = c(1, 2)
  )

  set.seed(6901)
  canonical <- do.call(
    sim_comp_data,
    c(
      common_args,
      list(prop_loss = c(control = 0.10, treatment = 0.25))
    )
  )
  set.seed(6901)
  reversed <- do.call(
    sim_comp_data,
    c(
      common_args,
      list(prop_loss = c(treatment = 0.25, control = 0.10))
    )
  )

  expect_identical(reversed, canonical)
  expect_equal(
    sum(canonical$loss_to_fu[canonical$treatment == 0L]),
    ceiling(0.10 * sum(canonical$treatment == 0L))
  )
  expect_equal(
    sum(canonical$loss_to_fu[canonical$treatment == 1L]),
    ceiling(0.25 * sum(canonical$treatment == 1L))
  )
  expect_true(all(canonical$event[canonical$loss_to_fu] == 0L))
})

test_that("sim_comp_data validates arm-specific loss proportions", {
  common_args <- list(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = NULL,
    N_total = 30,
    lambda = 5,
    lambda_time = NULL,
    end_of_study = 36
  )

  expect_error(
    do.call(
      sim_comp_data,
      c(common_args, list(prop_loss = c(0.10, 0.20)))
    ),
    "named 'control' and 'treatment'"
  )
  expect_error(
    do.call(
      sim_comp_data,
      c(
        common_args,
        list(prop_loss = c(control = 0.10, experimental = 0.20))
      )
    ),
    "named 'control' and 'treatment'"
  )
  expect_error(
    do.call(
      sim_comp_data,
      c(
        common_args,
        list(prop_loss = c(control = 0.10, treatment = 1.20))
      )
    ),
    "finite probabilities"
  )
  expect_error(
    do.call(
      sim_comp_data,
      modifyList(
        common_args,
        list(
          hazard_control = NULL,
          prop_loss = c(control = 0.10, treatment = 0.20)
        )
      )
    ),
    "single probability"
  )
})

test_that("sim_comp_data has no loss to follow-up by default", {
  set.seed(4192)
  out <- sim_comp_data(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = NULL,
    N_total = 50,
    lambda = 5,
    lambda_time = NULL,
    end_of_study = 36
  )
  expect_true(all(!out$loss_to_fu))
})

test_that("sim_comp_data works with piecewise hazard", {
  set.seed(6083)
  out <- sim_comp_data(
    hazard_treatment = c(0.005, 0.01),
    hazard_control = c(0.01, 0.02),
    cutpoints = 12,
    N_total = 80,
    lambda = 10,
    lambda_time = NULL,
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
      cutpoints = NULL,
      N_total = 50.5,
      lambda = 5,
      lambda_time = NULL,
      end_of_study = 36
    ),
    "N_total"
  )

  expect_error(
    sim_comp_data(
      hazard_treatment = 0.01,
      hazard_control = 0.02,
      cutpoints = NULL,
      N_total = 50,
      lambda = 5,
      lambda_time = NULL,
      end_of_study = 36,
      prop_loss = -0.1
    ),
    "prop_loss"
  )
})

test_that("sim_comp_data validates treatment-arm hazards", {
  common_args <- list(
    hazard_treatment = 0.01,
    hazard_control = 0.02,
    cutpoints = NULL,
    N_total = 50,
    lambda = 5,
    lambda_time = NULL,
    end_of_study = 36
  )

  expect_error(
    do.call(
      sim_comp_data,
      modifyList(common_args, list(hazard_treatment = Inf))
    ),
    "hazard_treatment"
  )
  expect_error(
    do.call(
      sim_comp_data,
      modifyList(common_args, list(hazard_control = NA_real_))
    ),
    "hazard_control"
  )
})

test_that("sim_comp_data validates the complete piecewise model", {
  common_args <- list(
    hazard_treatment = c(0.01, 0.02),
    hazard_control = c(0.02, 0.03),
    cutpoints = 12,
    N_total = 50,
    lambda = 5,
    lambda_time = NULL,
    end_of_study = 36
  )

  expect_error(
    do.call(sim_comp_data, modifyList(common_args, list(cutpoints = 100))),
    "end_of_study"
  )
  expect_error(
    do.call(
      sim_comp_data,
      modifyList(common_args, list(cutpoints = c(12, 12)))
    ),
    "strictly increasing"
  )
})
