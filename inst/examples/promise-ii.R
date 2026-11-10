# PROMISE II: the planned 120-patient design.
# Author: Graeme L. Hickey
# Sources: protocol Rev. 2, Tables 6-7; SAP Rev. 5, Sections 6.2-6.7,
# Tables 3 and 5, and Appendix 2. See the PROMISE II vignette for references.

promise_ii_design <- function() {
  list(
    N_total = 120L,
    interim_look = c(60L, 75L, 90L, 105L),
    end_of_study = 180,
    cutpoints = c(30, 60, 90),
    base_hazard = c(0.005417298, 0.00647187, 0.004045362, 0.0015350038),
    prior_surv = c(shape = 0.001, rate = 0.001),
    prior_bin = c(1, 1),
    h0 = 0.46,
    prob_ha = 0.977,
    Sn = c(0.98, 0.95, 0.90, 0.90),
    Fn = c(0.05, 0.10, 0.15, 0.15),
    prop_loss = 0.10,
    days_per_month = 30,
    accrual = list(
      expected = c(1, 2, 2, 3, 3, 4, 4, 5),
      slower = c(1, 1, 2, 2, 3),
      faster = c(2, 3, 4, 5, 6, 7)
    )
  )
}

promise_ii_hazard <- function(afs, design = promise_ii_design()) {
  stopifnot(length(afs) == 1L, is.finite(afs), afs > 0, afs < 1)
  widths <- diff(c(0, design$cutpoints, design$end_of_study))
  design$base_hazard * (-log(afs) / sum(design$base_hazard * widths))
}

promise_ii_enrollment <- function(accrual, enrollment_model, design) {
  rates <- design$accrual[[accrual]] / design$days_per_month
  knots <- seq_len(length(rates) - 1L) * design$days_per_month
  if (enrollment_model == "poisson") {
    return(goldilocks::enrollment(rates, design$N_total, knots))
  }
  # Sensitivity analysis: fixed spacing on the cumulative accrual-rate scale.
  # The first patient enters at zero under both enrollment assumptions.
  starts <- c(0, knots)
  cumulative <- c(0, cumsum(diff(starts) * head(rates, -1L)))
  target <- seq_len(design$N_total) - 1L
  interval <- findInterval(target, cumulative)
  starts[interval] + (target - cumulative[interval]) / rates[interval]
}

promise_ii_generate <- function(
  afs, accrual = "expected", enrollment_model = "poisson",
  design = promise_ii_design()
) {
  stopifnot(accrual %in% names(design$accrual),
            enrollment_model %in% c("poisson", "fixed"))
  enrollment <- promise_ii_enrollment(accrual, enrollment_model, design)
  lost <- stats::rbinom(design$N_total, 1, design$prop_loss) == 1L
  dropout <- rep(Inf, design$N_total)
  dropout[lost] <- stats::runif(sum(lost), 0, design$end_of_study)
  event_time <- goldilocks::pwe_sim(
    design$N_total, promise_ii_hazard(afs, design), design$cutpoints
  )$time
  data.frame(id = seq_len(design$N_total), treatment = 1L,
             enrollment = enrollment, event_time = event_time,
             dropout_time = dropout)
}

promise_ii_cut <- function(data, n, final = FALSE, design = promise_ii_design()) {
  data <- data[seq_len(n), , drop = FALSE]
  data_cut <- data$enrollment[n] + if (final) design$end_of_study else 0
  available <- pmin(data_cut - data$enrollment, design$end_of_study)
  time <- pmin(data$event_time, data$dropout_time, available)
  event <- as.integer(data$event_time <= time)
  status <- ifelse(event == 1L, "event",
                   ifelse(time >= design$end_of_study, "complete",
                          ifelse(data$dropout_time <= time, "censored", "pending")))
  out <- data.frame(id = data$id, treatment = data$treatment,
                    enrollment = data$enrollment, time = time,
                    event = event, status = status)
  attr(out, "data_cut") <- data_cut
  out
}

promise_ii_sufficient <- function(data, design = promise_ii_design()) {
  lower <- c(0, design$cutpoints)
  upper <- c(design$cutpoints, design$end_of_study)
  exposure <- vapply(seq_along(lower), function(j) {
    sum(pmax(0, pmin(data$time, upper[j]) - lower[j]))
  }, numeric(1))
  events <- vapply(seq_along(lower), function(j) {
    sum(data$event == 1L & data$time <= upper[j] &
          if (j == 1L) data$time >= 0 else data$time > lower[j])
  }, numeric(1))
  list(shape = design$prior_surv[1] + events,
       rate = design$prior_surv[2] + exposure)
}

# One column per hazard interval and one row per completed-data draw.
# A single hazard vector is shared by all subjects within each draw.
promise_ii_counts <- function(
  data, N_impute, n_future = 0L, design = promise_ii_design()
) {
  stopifnot(N_impute >= 1L, N_impute == as.integer(N_impute), n_future >= 0L)
  missing <- data$event == 0L & data$time < design$end_of_study
  counts <- rep(sum(data$event), N_impute)
  if (!any(missing) && n_future == 0L) {
    return(list(current = counts, maximum = counts))
  }
  post <- promise_ii_sufficient(data, design)
  hazards <- vapply(seq_along(post$shape), function(j) {
    stats::rgamma(N_impute, post$shape[j], rate = post$rate[j])
  }, numeric(N_impute))
  hazards <- matrix(hazards, nrow = N_impute)
  lower <- c(0, design$cutpoints)
  upper <- c(design$cutpoints, design$end_of_study)
  for (q in data$time[missing]) {
    remaining <- pmax(0, upper - pmax(q, lower))
    p_event <- -expm1(-drop(hazards %*% remaining))
    counts <- counts + stats::rbinom(N_impute, 1L, p_event)
  }
  maximum <- counts
  if (n_future > 0L) {
    # Future subjects are conditionally exchangeable given the hazard draw.
    p_event <- -expm1(-drop(hazards %*% (upper - lower)))
    maximum <- counts + stats::rbinom(N_impute, n_future, p_event)
  }
  list(current = counts, maximum = maximum)
}

promise_ii_score <- function(events, n, design = promise_ii_design()) {
  # Work on the event scale: Pr(p_event < 0.46 | completed binary data).
  lookup <- stats::pbeta(design$h0, design$prior_bin[1] + 0:n,
                        design$prior_bin[2] + n - 0:n)
  lookup[events + 1L]
}

promise_ii_predict <- function(data, N_impute, design = promise_ii_design()) {
  counts <- promise_ii_counts(data, N_impute, design$N_total - nrow(data), design)
  c(current = mean(promise_ii_score(counts$current, nrow(data), design) >
                     design$prob_ha),
    maximum = mean(promise_ii_score(counts$maximum, design$N_total, design) >
                     design$prob_ha))
}

promise_ii_package_interim <- function(
  data, look, N_impute, design = promise_ii_design(), seed = NULL
) {
  goldilocks::evaluate_interim(
    data = data, data_cut = attr(data, "data_cut"), look = look,
    N_total = design$N_total, end_of_study = design$end_of_study,
    cutpoints = design$cutpoints, prior_surv = design$prior_surv,
    prior_bin = design$prior_bin, bin_method = "quadrature",
    single_arm = TRUE, alternative = "less", h0 = design$h0,
    Fn = design$Fn[look], Sn = design$Sn[look], Qn = 1,
    prob_ha = design$prob_ha, N_impute = N_impute, method = "bayes-bin",
    binary_imputation = "bernoulli", empty_interval = "prior", seed = seed
  )
}

promise_ii_final <- function(data, N_impute, design = promise_ii_design()) {
  counts <- promise_ii_counts(data, N_impute, design = design)$current
  # Final inference averages posterior probabilities, not success indicators.
  mean(promise_ii_score(counts, nrow(data), design))
}

promise_ii_decision <- function(probabilities, look, design = promise_ii_design()) {
  if (probabilities[1] > design$Sn[look]) return("promising")
  if (probabilities[2] < design$Fn[look]) return("futile")
  "continue"
}

promise_ii_trial <- function(
  afs, N_impute, accrual, enrollment_model, stream,
  design = promise_ii_design(), keep_trace = FALSE
) {
  assign(".Random.seed", stream, envir = .GlobalEnv)
  data <- promise_ii_generate(afs, accrual, enrollment_model, design)
  reason <- "maximum"
  n <- design$N_total
  trace <- vector("list", length(design$interim_look))
  for (look in seq_along(design$interim_look)) {
    stream <- parallel::nextRNGSubStream(stream)
    assign(".Random.seed", stream, envir = .GlobalEnv)
    cut <- promise_ii_cut(data, design$interim_look[look], design = design)
    pp <- promise_ii_predict(cut, N_impute, design)
    decision <- promise_ii_decision(pp, look, design)
    trace[[look]] <- data.frame(n = nrow(cut), pp_current = pp[1],
                              pp_maximum = pp[2], decision = decision)
    if (decision != "continue") {
      reason <- decision
      n <- nrow(cut)
      break
    }
  }
  stream <- parallel::nextRNGSubStream(stream)
  assign(".Random.seed", stream, envir = .GlobalEnv)
  final <- promise_ii_cut(data, n, final = TRUE, design = design)
  probability <- promise_ii_final(final, N_impute, design)
  result <- data.frame(n = n, reason = reason, post_prob = probability,
                       success = probability > design$prob_ha,
                       analysis_day = attr(final, "data_cut"))
  if (keep_trace) {
    return(list(result = result, trace = do.call(rbind, trace), data = data))
  }
  result
}

promise_ii_simulate <- function(
  afs, N_trials, N_impute = 10000L, seed = 11012026L,
  accrual = "expected", enrollment_model = "poisson", ncores = 1L,
  design = promise_ii_design()
) {
  stopifnot(N_trials >= 2L, N_trials == as.integer(N_trials),
            N_impute >= 1L, N_impute == as.integer(N_impute), ncores >= 1L,
            accrual %in% names(design$accrual),
            enrollment_model %in% c("poisson", "fixed"))
  old_kind <- RNGkind()
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  streams <- vector("list", N_trials)
  streams[[1]] <- .Random.seed
  for (i in 2:N_trials) streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
  run <- function(i) {
    promise_ii_trial(afs, N_impute, accrual, enrollment_model, streams[[i]], design)
  }
  started <- proc.time()[[3]]
  trials <- if (ncores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(N_trials), run, mc.cores = ncores, mc.set.seed = FALSE)
  } else lapply(seq_len(N_trials), run)
  failed <- vapply(trials, inherits, logical(1), what = "try-error")
  if (any(failed)) stop("PROMISE II simulation failed; no failed trials are omitted")
  trials <- do.call(rbind, trials)
  rownames(trials) <- NULL
  trials$trial <- seq_len(N_trials)
  list(
    trials = trials,
    settings = list(afs = afs, N_trials = N_trials, N_impute = N_impute,
                    seed = seed, accrual = accrual, enrollment_model = enrollment_model,
                    ncores = ncores, rng = "L'Ecuyer-CMRG: trial streams, look substreams",
                    design = design),
    elapsed_seconds = proc.time()[[3]] - started
  )
}

promise_ii_summary <- function(simulation) {
  d <- simulation$trials
  p <- mean(d$success)
  ci <- stats::binom.test(sum(d$success), nrow(d))$conf.int
  data.frame(afs = simulation$settings$afs, N_trials = nrow(d), power = p,
             mcse = sqrt(p * (1 - p) / nrow(d)), lower = ci[1], upper = ci[2],
             mean_n = mean(d$n), mean_n_mcse = stats::sd(d$n) / sqrt(nrow(d)),
             futile_then_success = mean(d$reason == "futile" & d$success))
}

promise_ii_reference <- function() {
  # SAP Rev. 5, Table 3: pairs are futile first, promising second.
  data.frame(
    afs = c(.54, .65, .66, .67, .68, .69, .70, .75),
    power = c(.0244, .6336, .7130, .7810, .8423, .8879, .9253, .9945),
    mean_n = c(76.0, 95.2, 94.3, 92.3, 90.3, 87.9, 84.6, 71.7),
    futile_60 = c(.4647, .0529, .0355, .0260, .0208, .0112, .0093, .0005),
    promising_60 = c(.0019, .0911, .1181, .1486, .1726, .2190, .2654, .5219),
    futile_75 = c(.2317, .0594, .0422, .0314, .0228, .0158, .0099, .0004),
    promising_75 = c(.0033, .1159, .1336, .1578, .1865, .2049, .2287, .2699),
    futile_90 = c(.1379, .0594, .0478, .0374, .0235, .0171, .0133, .0011),
    promising_90 = c(.0049, .1287, .1475, .1654, .1786, .1783, .1846, .1316),
    futile_105 = c(.0688, .0466, .0357, .0238, .0180, .0118, .0080, .0001),
    promising_105 = c(.0058, .1294, .1429, .1482, .1550, .1529, .1424, .0518),
    maximum_120 = c(.0810, .3166, .2967, .2614, .2222, .1890, .1384, .0227)
  )
}
