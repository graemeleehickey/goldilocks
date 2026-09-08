test_that("RMST Wald tests calibrate under zero and nonzero RMST nulls", {
  skip_on_cran()
  truth <- function(h) {
    integrate(
      function(t) exp(-h[1] * pmin(t, 3) - h[2] * pmax(t - 3, 0)),
      0,
      9,
      subdivisions = 200L,
      rel.tol = 1e-10
    )$value
  }
  control <- c(0.1, 0.1)
  late <- uniroot(
    function(h) truth(c(0.04, h)) - truth(control),
    c(0.1, 1),
    tol = 1e-10
  )$root
  scenarios <- list(
    equal = control,
    crossing = c(0.04, late),
    nonzero = c(0.12, 0.12)
  )
  for (i in seq_along(scenarios)) {
    h <- scenarios[[i]]
    margin <- truth(h) - truth(control)
    out <- sim_trials(
      hazard_treatment = h,
      hazard_control = control,
      generation_cutpoints = 3,
      N_total = 300,
      lambda = 100,
      end_of_study = 12,
      rmst_tau = 9,
      method = "rmst",
      h0 = margin,
      alternative = "greater",
      prob_ha = 0.975,
      prop_loss = c(control = 0.1, treatment = 0.2),
      N_trials = 500,
      N_impute = 1,
      backend = "sequential",
      seed = 1450 + i
    )
    expect_equal(nrow(out$failures), 0)
    expect_mc_rate(
      out$sims$trial_success,
      0.025,
      paste("RMST one-sided type-I error:", names(scenarios)[i])
    )
  }
})

test_that("censored RMST Wald confidence intervals cover the true difference", {
  skip_on_cran()
  tau <- 9
  true_difference <- -expm1(-0.06 * tau) / 0.06 - -expm1(-0.1 * tau) / 0.1
  set.seed(1454)
  covered <- replicate(500, {
    d <- sim_comp_data(
      hazard_control = 0.1,
      hazard_treatment = 0.06,
      N_total = 300,
      end_of_study = 12,
      prop_loss = c(control = 0.1, treatment = 0.2)
    )
    fit <- rmst_estimate(d$time, d$event, d$treatment, tau)
    abs(fit$estimate - true_difference) <= qnorm(0.975) * fit$std_error
  })
  expect_mc_rate(covered, 0.95, "RMST 95% Wald confidence-interval coverage")
})
