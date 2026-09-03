test_that("complete-data simulation is invariant to the time unit", {
  scale <- 30
  month_arguments <- list(
    hazard_treatment = c(0.03, 0.01),
    hazard_control = c(0.04, 0.02),
    generation_cutpoints = 6,
    N_total = 40,
    lambda = c(2, 3),
    lambda_time = 4,
    end_of_study = 12,
    block = 2,
    rand_ratio = c(control = 1, treatment = 1),
    prop_loss = 0
  )
  day_arguments <- month_arguments
  day_arguments$hazard_treatment <-
    day_arguments$hazard_treatment / scale
  day_arguments$hazard_control <- day_arguments$hazard_control / scale
  day_arguments$generation_cutpoints <-
    day_arguments$generation_cutpoints * scale
  day_arguments$lambda <- day_arguments$lambda / scale
  day_arguments$lambda_time <- day_arguments$lambda_time * scale
  day_arguments$end_of_study <- day_arguments$end_of_study * scale

  set.seed(806)
  months <- do.call(sim_comp_data, month_arguments)
  set.seed(806)
  days <- do.call(sim_comp_data, day_arguments)

  expect_equal(days$time, months$time * scale, tolerance = 1e-12)
  expect_equal(
    days$enrollment,
    months$enrollment * scale,
    tolerance = 1e-12
  )
  expect_identical(
    days[c("treatment", "event", "id", "loss_to_fu")],
    months[c("treatment", "event", "id", "loss_to_fu")]
  )
})

test_that("predictive interim decisions are invariant to the time unit", {
  scale <- 30
  months_data <- data.frame(
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
  days_data <- months_data
  days_data$enrollment <- days_data$enrollment * scale
  days_data$time <- days_data$time * scale

  common_arguments <- list(
    look = 2,
    N_total = 12,
    rand_ratio = c(control = 1, treatment = 1),
    alternative = "less",
    N_impute = 200,
    N_mcmc = 1,
    method = "logrank",
    seed = 807
  )
  months <- do.call(
    evaluate_interim,
    c(
      list(
        data = months_data,
        data_cut = 8,
        end_of_study = 10,
        prior_surv = c(0.1, 0.1)
      ),
      common_arguments
    )
  )
  days <- do.call(
    evaluate_interim,
    c(
      list(
        data = days_data,
        data_cut = 8 * scale,
        end_of_study = 10 * scale,
        prior_surv = c(0.1, 0.1 * scale)
      ),
      common_arguments
    )
  )

  expect_identical(days$probabilities, months$probabilities)
  expect_identical(days$monte_carlo, months$monte_carlo)
  expect_identical(days$decision$decision, months$decision$decision)
})
