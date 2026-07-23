test_that("survival_adapt-bayes-surv", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "bayes-surv"
  )

  expect_s3_class(out, "data.frame")
})

test_that("survival_adapt-bayes-bin", {
  set.seed(1)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    bin_prior = c(1, 1),
    bin_method = "normal",
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "bayes-bin"
  )

  expect_s3_class(out, "data.frame")
  expect_true(!is.na(out$est_final))
})

test_that("survival_adapt supports Bernoulli binary imputation", {
  set.seed(2102)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    N_total = 100,
    lambda = 20,
    interim_look = 50,
    end_of_study = 36,
    bin_method = "normal",
    alternative = "less",
    N_impute = 2,
    N_mcmc = 2,
    method = "bayes-bin",
    binary_imputation = "bernoulli"
  )

  expect_s3_class(out, "data.frame")
  expect_true(out$ppp_success >= 0 && out$ppp_success <= 1)
})

test_that("survival_adapt-logrank", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "two.sided",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "logrank"
  )

  expect_s3_class(out, "data.frame")
})

test_that("survival_adapt-cox", {
  set.seed(2)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "two.sided",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "cox"
  )

  expect_s3_class(out, "data.frame")
  expect_true(!is.na(out$est_final))
})

test_that("survival_adapt-riskdiff", {
  set.seed(1)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "two.sided",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "riskdiff"
  )

  expect_s3_class(out, "data.frame")
  expect_true(!is.na(out$est_final))
})

test_that("survival_adapt-riskdiff excludes LTFU when imputed_final = FALSE", {
  set.seed(3927)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "two.sided",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "riskdiff",
    imputed_final = FALSE
  )

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
})

test_that("survival_adapt pools imputed Cox final analyses", {
  set.seed(8263)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 200,
    lambda = 20,
    lambda_time = NULL,
    interim_look = NULL,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0,
    prob_ha = 0.95,
    N_impute = 5,
    method = "cox",
    imputed_final = TRUE
  )

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
  expect_true(is.finite(out$est_final))
})

test_that("survival_adapt pools imputed risk-difference final analyses", {
  set.seed(2084)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 200,
    lambda = 20,
    lambda_time = NULL,
    interim_look = NULL,
    end_of_study = 36,
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0,
    prob_ha = 0.95,
    N_impute = 5,
    method = "riskdiff",
    imputed_final = TRUE
  )

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
  expect_true(is.finite(out$est_final))
})

test_that("survival_adapt requires multiple imputations for Rubin pooling", {
  for (method in c("cox", "riskdiff")) {
    expect_error(
      survival_adapt(
        hazard_treatment = -log(0.85) / 36,
        hazard_control = -log(0.7) / 36,
        cutpoints = NULL,
        N_total = 200,
        lambda = 20,
        lambda_time = NULL,
        interim_look = NULL,
        end_of_study = 36,
        prop_loss = 0.30,
        alternative = "less",
        method = method,
        N_impute = 1,
        imputed_final = TRUE
      ),
      "at least two imputations"
    )
  }
})

test_that("survival_adapt still rejects imputed log-rank final analyses", {
  expect_error(
    survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = NULL,
      N_total = 200,
      lambda = 20,
      lambda_time = NULL,
      interim_look = NULL,
      end_of_study = 36,
      alternative = "two.sided",
      method = "logrank",
      imputed_final = TRUE
    ),
    "no supported frequentist pooling rule"
  )
})

test_that("survival_adapt-complex", {
  skip_on_cran()

  hc <- prop_to_haz(c(0.20, 0.25, 0.30), c(3, 6), 24)
  ht <- prop_to_haz(c(0.075, 0.125, 0.175), c(3, 6), 24)

  set.seed(12345)
  out <- survival_adapt(
    hazard_treatment = ht,
    hazard_control = hc,
    cutpoints = c(3, 6),
    N_total = 500,
    lambda = c(6, 8),
    lambda_time = 6,
    interim_look = c(200, 350),
    end_of_study = 24,
    prior = c(0.1, 0.1),
    block = 3,
    rand_ratio = c(2, 1),
    prop_loss = 0,
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.95,
    N_impute = 20,
    N_mcmc = 20,
    method = "bayes-surv",
    imputed_final = TRUE
  )

  expect_s3_class(out, "data.frame")
})

test_that("error-interim-looks", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = NULL,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(200, 300, 400, 500),
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank"
    )
  )
})

test_that("error-interim-look-below-block-size", {
  # A two-arm interim look smaller than the block size can enroll one treatment
  # group only, leaving the interim posterior undefined for the missing group.
  # This should be caught as an input error before any simulation runs.
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = NULL,
      N_total = 100,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 3,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 4,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "less",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "bayes-surv"
    ),
    "must be at least the block size"
  )
})

test_that("error-alternative", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = NULL,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two-sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank"
    )
  )
})

test_that("error-cutpoint", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = 37,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank"
    )
  )
})

test_that("error-alternative-bayes-surv", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = NULL,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "bayes-surv"
    )
  )
})

test_that("survival_adapt-logrank-one-sided", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "logrank"
  )

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
})

test_that("survival_adapt-cox-one-sided-less", {
  set.seed(2)

  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "cox"
  )

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
  expect_true(!is.na(out$est_final))
})

test_that("survival_adapt-cox-one-sided-greater", {
  set.seed(2)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "greater",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "cox"
  )

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
})

test_that("survival_adapt-riskdiff supports one-sided margins", {
  set.seed(908)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0.05,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "riskdiff"
  )

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
  expect_true(out$est_final >= -1 && out$est_final <= 1)
})

test_that("survival_adapt no longer accepts method = 'chisq'", {
  expect_error(
    survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      N_total = 100,
      end_of_study = 36,
      method = "chisq"
    ),
    "or 'riskdiff'"
  )
})

test_that("survival_adapt validates risk-difference margins", {
  expect_error(
    survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      N_total = 100,
      end_of_study = 36,
      method = "riskdiff",
      h0 = 1.1
    ),
    "\\[-1, 1\\]"
  )
})

test_that("error-logrank-single_arm", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = NULL,
      cutpoints = NULL,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank"
    )
  )
})

test_that("error-prob-thresholds-length_v1", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = NULL,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(100, 200),
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = c(0.05, 0.05, 0.05),
      Sn = c(0.99, 0.9),
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank"
    )
  )
})

test_that("survival_adapt works with no interim looks and default thresholds", {
  set.seed(4817)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 200,
    lambda = 20,
    lambda_time = NULL,
    interim_look = NULL,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0,
    alternative = "two.sided",
    h0 = 0,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "logrank"
  )

  expect_s3_class(out, "data.frame")
  expect_equal(out$N_enrolled, 200)
  expect_equal(out$stop_futility, 0)
  expect_equal(out$stop_expected_success, 0)
  expect_true(is.na(out$ppp_success))
})

test_that("survival_adapt validates probability, count, and prior arguments", {
  common_args <- list(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 200,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 100,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "two.sided",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "logrank"
  )

  expect_error(
    do.call(survival_adapt, modifyList(common_args, list(prob_ha = -0.1))),
    "prob_ha"
  )
  expect_error(
    do.call(survival_adapt, modifyList(common_args, list(Sn = 1.2))),
    "Sn"
  )
  expect_error(
    do.call(
      survival_adapt,
      modifyList(common_args, list(prior = c(-0.1, 0.1)))
    ),
    "prior"
  )
  expect_error(
    do.call(survival_adapt, modifyList(common_args, list(N_impute = 1.5))),
    "N_impute"
  )
  expect_error(
    do.call(survival_adapt, modifyList(common_args, list(N_mcmc = 1.5))),
    "N_mcmc"
  )
  expect_error(
    do.call(survival_adapt, modifyList(common_args, list(N_total = 200.5))),
    "N_total"
  )
  expect_error(
    do.call(survival_adapt, modifyList(common_args, list(prop_loss = -0.1))),
    "prop_loss"
  )
})

test_that("survival_adapt validates interim ordering and h0", {
  common_args <- list(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 200,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 100,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0,
    alternative = "less",
    h0 = 0,
    Fn = 0,
    Sn = 1,
    prob_ha = 0.95,
    N_impute = 1,
    N_mcmc = 1,
    method = "bayes-surv"
  )

  expect_error(
    do.call(
      survival_adapt,
      modifyList(common_args, list(interim_look = c(150, 100)))
    ),
    "strictly increasing"
  )
  expect_error(
    do.call(
      survival_adapt,
      modifyList(common_args, list(interim_look = 200))
    ),
    "strictly less than 'N_total'"
  )
  expect_error(
    do.call(
      survival_adapt,
      modifyList(common_args, list(interim_look = c(100, 100)))
    ),
    "strictly increasing"
  )
  expect_error(
    do.call(survival_adapt, modifyList(common_args, list(h0 = Inf))),
    "single finite"
  )
  expect_error(
    do.call(survival_adapt, modifyList(common_args, list(h0 = 1.1))),
    "\\[-1, 1\\]"
  )
  expect_error(
    do.call(
      survival_adapt,
      modifyList(
        common_args,
        list(hazard_control = NULL, h0 = -0.1)
      )
    ),
    "\\[0, 1\\]"
  )
})

test_that("survival_adapt keeps futility disabled when Fn = 0", {
  set.seed(1)
  out <- survival_adapt(
    hazard_treatment = -log(0.3) / 36,
    hazard_control = -log(0.2) / 36,
    cutpoints = NULL,
    N_total = 200,
    lambda = 20,
    lambda_time = NULL,
    interim_look = c(100, 150),
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "two.sided",
    h0 = 0,
    Fn = c(0, 0),
    Sn = c(1, 1),
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "logrank"
  )

  expect_s3_class(out, "data.frame")
  expect_equal(out$stop_futility, 0)
})

test_that("error-prob-thresholds-length_v2", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = NULL,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(100, 200),
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = c(0.05, 0.05, 0.05),
      Sn = c(0.99, 0.99, 0.99),
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank"
    )
  )
})
