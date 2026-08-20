test_that("sim_trials-logrank", {
  hc <- prop_to_haz(c(0.20, 0.30), 12, 36)
  ht <- prop_to_haz(c(0.05, 0.15), 12, 36)

  out <- sim_trials(
    hazard_treatment = ht,
    hazard_control = hc,
    cutpoints = 12,
    N_total = 600,
    lambda = 20,
    lambda_time = NULL,
    interim_look = c(400, 500),
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "two.sided",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 5,
    N_mcmc = 2,
    N_trials = 2,
    method = "logrank",
    ncores = 1
  )

  expect_type(out, "list")
  expect_s3_class(out$sims, "data.frame")
  expect_named(
    out$failures,
    c("trial", "error_class", "message")
  )
  expect_equal(nrow(out$failures), 0L)

  summ_out <- summarise_sims(out$sims)
  expect_s3_class(summ_out, "data.frame")
})

test_that("sim_trials-bayes-surv", {
  hc <- prop_to_haz(c(0.20, 0.30), 12, 36)
  ht <- prop_to_haz(c(0.05, 0.15), 12, 36)

  out <- sim_trials(
    hazard_treatment = ht,
    hazard_control = hc,
    cutpoints = 12,
    N_total = 600,
    lambda = 20,
    lambda_time = NULL,
    interim_look = c(400, 500),
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 5,
    N_mcmc = 2,
    N_trials = 2,
    method = "bayes-surv",
    ncores = 1
  )

  expect_type(out, "list")
  expect_s3_class(out$sims, "data.frame")

  summ_out <- summarise_sims(out$sims)
  expect_s3_class(summ_out, "data.frame")
})

test_that("sim_trials-bayes-bin", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  out <- sim_trials(
    hazard_treatment = ht,
    hazard_control = hc,
    cutpoints = NULL,
    N_total = 200,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 100,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    prior_bin = c(1, 1),
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
    N_trials = 1,
    method = "bayes-bin",
    binary_imputation = "bernoulli",
    ncores = 1,
    seed = 9241
  )

  expect_type(out, "list")
  expect_s3_class(out$sims, "data.frame")
})

test_that("sim_trials rejects invalid ncores", {
  hc <- prop_to_haz(c(0.20, 0.30), 12, 36)
  ht <- prop_to_haz(c(0.05, 0.15), 12, 36)

  expect_error(
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = 12,
      N_total = 600,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(400, 500),
      end_of_study = 36,
      prior_surv = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "less",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 5,
      N_mcmc = 2,
      N_trials = 2,
      method = "bayes-surv",
      ncores = NULL
    )
  )
})

test_that("sim_trials defaults to serial execution", {
  expect_identical(formals(sim_trials)$ncores, 1L)
  expect_identical(formals(sim_trials)$return_trace, FALSE)
  expect_equal(
    as.character(formals(sim_trials)$backend)[-1],
    c("auto", "fork", "psock", "sequential")
  )
})

test_that("sim_trials retains complete evaluated arguments", {
  out <- sim_trials(
    hazard_treatment = 0.02,
    hazard_control = 0.03,
    N_total = 20,
    lambda = 10,
    end_of_study = 12,
    rand_ratio = c(treatment = 1, control = 1),
    prop_loss = c(treatment = 0.20, control = 0.10),
    N_trials = 2,
    method = "logrank",
    backend = "sequential",
    seed = 8601
  )
  arguments <- attr(out, "arguments", exact = TRUE)

  expect_named(arguments, names(formals(sim_trials)))
  expect_identical(arguments$hazard_treatment, 0.02)
  expect_identical(arguments$empty_interval, "prior")
  expect_identical(arguments$binary_imputation, "event-time")
  expect_identical(arguments$backend, "sequential")
  expect_identical(arguments$seed, 8601)
  expect_identical(
    arguments$rand_ratio,
    c(control = 1, treatment = 1)
  )
  expect_identical(
    arguments$prop_loss,
    c(control = 0.10, treatment = 0.20)
  )
  expect_identical(
    unserialize(serialize(arguments, NULL)),
    arguments
  )

  replayed <- do.call(sim_trials, arguments)
  expect_identical(replayed$sims, out$sims)
  expect_identical(replayed$failures, out$failures)
})

test_that("sim_trials optionally retains traces without changing summaries", {
  args <- list(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 80,
    lambda = 20,
    lambda_time = NULL,
    interim_look = c(40, 60),
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.05,
    alternative = "less",
    h0 = 0,
    Fn = c(0.05, 0.05),
    Sn = c(0.95, 0.9),
    prob_ha = 0.95,
    N_impute = 2,
    N_mcmc = 2,
    N_trials = 2,
    method = "bayes-surv",
    ncores = 1,
    seed = 4102
  )

  compact <- do.call(sim_trials, args)
  traced <- do.call(sim_trials, c(args, list(return_trace = TRUE)))

  expect_equal(traced$sims, compact$sims)
  expect_named(traced, c("sims", "traces", "failures", "call"))
  expect_s3_class(traced$traces, "data.frame")
  expect_true(all(
    c(
      "trial",
      "look",
      "ppp_stop_now",
      "ppp_success_at_max",
      "decision"
    ) %in%
      names(traced$traces)
  ))
  expect_setequal(unique(traced$traces$trial), 1:2)
})

test_that("sim_trials isolates failed trials and warns once", {
  run_fake_trials <- function(backend, fail_trial = 2L) {
    trial <- 0L
    fake_survival_adapt <- function(..., return_trace = FALSE) {
      trial <<- trial + 1L
      value <- stats::runif(1)
      if (identical(trial, fail_trial)) {
        rlang::abort(
          "simulated non-estimable analysis",
          class = "simulated_trial_error"
        )
      }
      summary <- data.frame(
        prob_threshold = 0.95,
        margin = 0,
        alternative = "less",
        N_treatment = 5L,
        N_control = 5L,
        N_enrolled = 10L,
        N_max = 10L,
        post_prob_ha = 0.99,
        est_final = value,
        ppp_success = 0.5,
        stop_futility = FALSE,
        stop_expected_success = FALSE
      )
      attr(summary, "decision_design") <- list(interim_look = NULL)
      if (!return_trace) {
        return(summary)
      }
      result <- list(
        summary = summary,
        trace = data.frame(look = 1L, N_enrolled = 10L),
        call = quote(fake_survival_adapt())
      )
      attr(result, "decision_design") <- attr(
        summary,
        "decision_design",
        exact = TRUE
      )
      result
    }

    fake_cluster <- structure(vector("list", 2L), class = "cluster")
    local_mocked_bindings(
      survival_adapt = fake_survival_adapt,
      make_psock_callable = function(...) fake_survival_adapt,
      make_sim_cluster = function(...) fake_cluster,
      initialize_sim_cluster = function(...) invisible(NULL),
      run_sim_cluster = function(cluster, trial_index, fun) {
        lapply(trial_index, fun)
      },
      stop_sim_cluster = function(...) invisible(NULL),
      .package = "goldilocks"
    )

    warnings <- character()
    result <- withCallingHandlers(
      sim_trials(
        hazard_treatment = 0.02,
        hazard_control = 0.03,
        N_total = 10,
        lambda = 10,
        end_of_study = 12,
        N_trials = 3,
        method = "logrank",
        backend = backend,
        ncores = 2,
        return_trace = TRUE,
        seed = 7201
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
    list(result = result, warnings = warnings)
  }

  sequential <- run_fake_trials("sequential")
  psock <- run_fake_trials("psock")
  complete <- run_fake_trials("sequential", fail_trial = NA_integer_)

  expect_length(sequential$warnings, 1L)
  expect_length(psock$warnings, 1L)
  expect_length(complete$warnings, 0L)
  expect_match(
    sequential$warnings,
    "1 of 3 simulated trials failed and were excluded"
  )
  expect_identical(sequential$result$failures, psock$result$failures)
  expect_identical(
    sequential$result$failures,
    data.frame(
      trial = 2L,
      error_class = "simulated_trial_error",
      message = "simulated non-estimable analysis"
    )
  )
  expect_equal(sequential$result$sims, psock$result$sims)
  expect_equal(
    sequential$result$sims$est_final,
    complete$result$sims$est_final[c(1L, 3L)]
  )
  expect_setequal(sequential$result$traces$trial, c(1L, 3L))

  summary <- summarise_sims(sequential$result)
  expect_identical(summary$n_requested, 3L)
  expect_identical(summary$n_analyzed, 2L)
  expect_identical(summary$n_failed, 1L)
  expect_identical(
    attr(summary, "simulation_metadata", exact = TRUE)[[1]]$failures,
    sequential$result$failures
  )
})

test_that("sim_trials stops with structured details when every trial fails", {
  fake_survival_adapt <- function(...) {
    rlang::abort(
      "simulated terminal failure",
      class = "simulated_trial_error"
    )
  }
  local_mocked_bindings(
    survival_adapt = fake_survival_adapt,
    .package = "goldilocks"
  )

  error <- tryCatch(
    sim_trials(
      hazard_treatment = 0.02,
      hazard_control = 0.03,
      N_total = 10,
      lambda = 10,
      end_of_study = 12,
      N_trials = 2,
      method = "logrank",
      backend = "sequential",
      seed = 7202
    ),
    error = identity
  )

  expect_s3_class(error, "goldilocks_all_trials_failed")
  expect_match(
    conditionMessage(error),
    "All 2 simulated trials failed.*simulated terminal failure"
  )
  expect_identical(error$failures$trial, 1:2)
  expect_identical(
    error$failures$error_class,
    rep("simulated_trial_error", 2)
  )
})

test_that("sim_trials validates return_trace", {
  expect_error(
    sim_trials(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      N_total = 40,
      interim_look = 20,
      end_of_study = 36,
      return_trace = NA
    ),
    "return_trace"
  )
})

test_that("sim_trials uses survival_adapt decision-rule defaults", {
  expect_identical(
    formals(sim_trials)$alternative,
    formals(survival_adapt)$alternative
  )
  expect_identical(formals(sim_trials)$Fn, formals(survival_adapt)$Fn)
})

test_that("sim_trials validates N_trials", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  expect_error(
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = NULL,
      N_total = 50,
      lambda = 20,
      lambda_time = NULL,
      interim_look = NULL,
      end_of_study = 36,
      N_trials = 1.5,
      method = "logrank",
      ncores = 1
    ),
    "N_trials"
  )
})

test_that("sim_trials is reproducible with seed and ncores = 1", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  run_once <- function() {
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = NULL,
      N_total = 200,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 100,
      end_of_study = 36,
      prior_surv = c(0.1, 0.1),
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
      N_trials = 2,
      method = "logrank",
      ncores = 1,
      seed = 4101
    )
  }

  expect_equal(run_once()$sims, run_once()$sims)
})

test_that("sim_trials uses reproducible per-trial streams in parallel", {
  skip_on_os("windows")

  hc <- rep(-log(0.7) / 36, 2)
  ht <- rep(-log(0.85) / 36, 2)

  run_with_cores <- function(ncores) {
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = 12,
      N_total = 200,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 100,
      end_of_study = 36,
      prior_surv = rbind(
        shape = c(0.1, 2),
        rate = c(0.1, 20)
      ),
      prior_surv_final = c(0.1, 0.1),
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
      N_trials = 4,
      method = "logrank",
      empty_interval = "prior",
      return_trace = TRUE,
      ncores = ncores,
      seed = 4101
    )
  }

  sequential <- run_with_cores(1)
  parallel <- run_with_cores(2)
  expect_equal(sequential$sims, parallel$sims)
  expect_equal(sequential$traces, parallel$traces)
})

test_that("sim_trials produces identical seeded PSOCK results", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  run_with_backend <- function(backend, ncores) {
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = NULL,
      N_total = 100,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 50,
      end_of_study = 36,
      prior_surv = c(0.1, 0.1),
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
      N_trials = 2,
      method = "logrank",
      return_trace = TRUE,
      ncores = ncores,
      backend = backend,
      seed = 4102
    )
  }

  sequential <- run_with_backend("sequential", 1)
  psock <- run_with_backend("psock", 2)
  psock_repeated <- run_with_backend("psock", 2)
  expect_equal(sequential$sims, psock$sims)
  expect_equal(sequential$traces, psock$traces)
  expect_equal(psock_repeated$sims, psock$sims)
  expect_identical(
    attr(psock_repeated, "parallel_metadata", exact = TRUE),
    attr(psock, "parallel_metadata", exact = TRUE)
  )
})

test_that("PSOCK workers load the package namespace and compiled routines", {
  namespace <- environment(sim_trials)
  package_path <- getNamespaceInfo(namespace, "path")
  skip_if_not(
    file.exists(file.path(package_path, "Meta", "package.rds")),
    "requires an installed package namespace"
  )

  cluster <- make_sim_cluster(1L)
  on.exit(stop_sim_cluster(cluster), add = TRUE)
  initialize_sim_cluster(cluster)

  worker_state <- parallel::clusterCall(cluster, function() {
    namespace <- asNamespace("goldilocks")
    c(
      namespace_loaded = "goldilocks" %in% loadedNamespaces(),
      dll_loaded = "goldilocks" %in% names(getLoadedDLLs()),
      native_symbol_loaded = exists(
        "_goldilocks_logrank_instance",
        envir = namespace,
        inherits = FALSE
      )
    )
  })[[1L]]

  expect_true(all(worker_state))
})

test_that("unseeded PSOCK simulations derive streams from caller RNG state", {
  run_psock <- function(ncores) {
    sim_trials(
      hazard_treatment = 0.02,
      hazard_control = 0.03,
      N_total = 30,
      lambda = 10,
      end_of_study = 12,
      N_trials = 3,
      method = "logrank",
      backend = "psock",
      ncores = ncores,
      seed = NULL
    )
  }

  set.seed(7101)
  first <- run_psock(2)
  set.seed(7101)
  repeated <- run_psock(2)
  set.seed(7101)
  one_worker <- run_psock(1)

  expect_equal(repeated$sims, first$sims)
  expect_equal(one_worker$sims, first$sims)

  metadata <- attr(first, "rng_metadata", exact = TRUE)
  expect_identical(metadata$seed_policy, "caller_derived_psock")
  expect_identical(metadata$backend, "psock")
  expect_identical(metadata$ncores, 2L)
  expect_identical(metadata$stream_kind, "L'Ecuyer-CMRG")
  expect_length(metadata$caller_kind, 3)
  expect_true(is.numeric(metadata$stream_seed))
})

test_that("unseeded PSOCK advances the parent RNG by one documented draw", {
  set.seed(7102)
  invisible(sample.int(.Machine$integer.max - 1L, size = 1L))
  expected_state <- .Random.seed

  set.seed(7102)
  sim_trials(
    hazard_treatment = 0.02,
    hazard_control = 0.03,
    N_total = 20,
    lambda = 10,
    end_of_study = 12,
    N_trials = 2,
    method = "logrank",
    backend = "psock",
    ncores = 2,
    seed = NULL
  )

  expect_identical(.Random.seed, expected_state)
})

test_that("consecutive unseeded sequential simulations advance caller streams", {
  set.seed(7103)
  initial_state <- .Random.seed
  first <- sim_trials(
    hazard_treatment = 0.02,
    hazard_control = 0.03,
    N_total = 30,
    lambda = 10,
    end_of_study = 12,
    N_trials = 2,
    method = "logrank",
    backend = "sequential",
    seed = NULL
  )
  first_state <- .Random.seed
  second <- sim_trials(
    hazard_treatment = 0.02,
    hazard_control = 0.03,
    N_total = 30,
    lambda = 10,
    end_of_study = 12,
    N_trials = 2,
    method = "logrank",
    backend = "sequential",
    seed = NULL
  )

  expect_false(identical(first_state, initial_state))
  expect_false(identical(.Random.seed, first_state))
  expect_false(identical(first$sims, second$sims))
  expect_identical(
    attr(first, "rng_metadata", exact = TRUE)$seed_policy,
    "caller_state"
  )
})

test_that("sim_trials auto backend selects the platform default", {
  expected <- if (.Platform$OS.type == "windows") "psock" else "fork"
  expect_identical(resolve_sim_backend("auto", 2L), expected)
  expect_identical(resolve_sim_backend("auto", 1L), "sequential")
})

test_that("automatic execution uses a documented workload crossover", {
  expected_parallel <- if (.Platform$OS.type == "windows") "psock" else "fork"

  small <- resolve_sim_execution("auto", 8L, N_trials = 3L)
  expect_identical(small$backend, "sequential")
  expect_identical(small$workers, 1L)
  expect_identical(small$reason, "auto_small_workload")

  medium <- resolve_sim_execution("auto", 8L, N_trials = 10L)
  expect_identical(medium$backend, expected_parallel)
  expect_identical(medium$workers, 5L)
  expect_identical(medium$reason, "auto_parallel")

  explicit <- resolve_sim_execution("psock", 8L, N_trials = 3L)
  expect_identical(explicit$backend, "psock")
  expect_identical(explicit$workers, 3L)
})

test_that("small automatic workloads do not start parallel workers", {
  local_mocked_bindings(
    make_sim_cluster = function(...) stop("PSOCK cluster was started"),
    pbmclapply = function(...) stop("fork workers were started"),
    .package = "goldilocks"
  )

  out <- sim_trials(
    hazard_treatment = 0.02,
    hazard_control = 0.03,
    N_total = 20,
    lambda = 10,
    end_of_study = 12,
    N_trials = 3,
    method = "logrank",
    backend = "auto",
    ncores = 8,
    seed = 9301
  )

  metadata <- attr(out, "parallel_metadata", exact = TRUE)
  expect_identical(metadata$backend, "sequential")
  expect_identical(metadata$workers, 1L)
  expect_identical(metadata$selection_reason, "auto_small_workload")
})

test_that("package-owned PSOCK clusters are capped and stopped on errors", {
  made_workers <- NULL
  initialized <- FALSE
  stopped <- FALSE
  fake_cluster <- structure(vector("list", 2L), class = "cluster")

  local_mocked_bindings(
    make_sim_cluster = function(workers) {
      made_workers <<- workers
      fake_cluster
    },
    initialize_sim_cluster = function(cluster) {
      initialized <<- identical(cluster, fake_cluster)
      invisible(cluster)
    },
    run_sim_cluster = function(...) stop("simulated worker failure"),
    stop_sim_cluster = function(cluster) {
      stopped <<- identical(cluster, fake_cluster)
    },
    .package = "goldilocks"
  )

  expect_error(
    sim_trials(
      hazard_treatment = 0.02,
      hazard_control = 0.03,
      N_total = 20,
      lambda = 10,
      end_of_study = 12,
      N_trials = 2,
      method = "logrank",
      backend = "psock",
      ncores = 8,
      seed = 9302
    ),
    "simulated worker failure"
  )
  expect_identical(made_workers, 2L)
  expect_true(initialized)
  expect_true(stopped)
})

test_that("sim_trials preserves the caller RNG state when seeded", {
  set.seed(9012)
  before_kind <- RNGkind()
  before <- .Random.seed

  sim_trials(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = NULL,
    N_total = 50,
    lambda = 20,
    lambda_time = NULL,
    interim_look = NULL,
    end_of_study = 36,
    N_trials = 1,
    method = "logrank",
    ncores = 1,
    seed = 4101
  )

  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before)
})

test_that("sim_trials validates seed", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  expect_error(
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = NULL,
      N_total = 200,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 100,
      end_of_study = 36,
      prior_surv = c(0.1, 0.1),
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
      N_trials = 2,
      method = "logrank",
      ncores = 1,
      seed = c(1, 2)
    ),
    "single integer"
  )
})
