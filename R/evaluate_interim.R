#' Evaluate an externally observed interim data cut
#'
#' @description Applies the same posterior-predictive decision calculation used
#'   by [survival_adapt()] to subject-level data observed at one interim look.
#'   The function does not simulate a trial or modify `data`. Participants who
#'   could still be enrolled before the maximum sample size are constructed
#'   internally from `N_total`, the observed arm counts, and `rand_ratio`.
#'
#' @param data Data frame with one row per enrolled subject and columns `id`,
#'   `treatment`, `enrollment`, `time`, `event`, and `status`. Treatment is
#'   coded `1` for treatment and `0` for control; single-arm data use `1`.
#'   `enrollment` is measured from first participant randomization, which must
#'   be zero, and `time` is follow-up from that subject's randomization.
#' @param data_cut Non-negative calendar time of the interim data cut, measured
#'   from the same origin and in the same units as `enrollment`, `time`,
#'   `end_of_study`, and `cutpoints`.
#' @param look Positive integer identifying the prespecified interim look.
#' @param N_total Positive integer maximum sample size.
#' @param end_of_study Positive subject-level follow-up horizon.
#' @param rand_ratio Length-two positive integer maximum-sample allocation
#'   ratio. Name the values `control` and `treatment`; either order is accepted.
#'   A legacy unnamed vector is interpreted as `c(control, treatment)`. The
#'   maximum sample size must divide exactly according to this ratio. Ignored
#'   for single-arm designs.
#' @param single_arm Logical indicating whether the design has one treatment
#'   arm and no control arm.
#' @param Fn Single probability threshold for stopping for futility at this
#'   look. Futility is declared when predictive success at the maximum sample
#'   size is strictly less than `Fn`. Set `Fn = 0` or `NULL` to disable the
#'   maximum-sample calculation.
#' @param Sn Single probability threshold for stopping for expected success at
#'   this look. Expected success is declared when predictive success among the
#'   currently enrolled participants is strictly greater than `Sn`.
#' @param seed `NULL`, or a non-negative integer used for the predictive Monte
#'   Carlo calculation. A supplied seed makes the result reproducible and
#'   preserves the caller's random-number state. With `seed = NULL`, the call
#'   uses and advances the current random-number state.
#' @inheritParams survival_adapt
#'
#' @details The arm-specific maximum enrollment is
#'   `N_total * rand_ratio / sum(rand_ratio)`. Potential future accrual in each
#'   arm is its maximum enrollment minus its observed enrollment. Randomization
#'   block information and concealed future assignment order are unnecessary:
#'   the maximum-sample predictive calculation requires only the number of
#'   potential future participants in each arm.
#'
#'   The `data$status` values distinguish observed follow-up: `"event"` for an
#'   observed endpoint event, `"complete"` for event-free completion of
#'   `end_of_study`, `"pending"` for a subject still under follow-up, and
#'   `"censored"` for permanent early censoring. Pending and censored outcomes
#'   are predictively imputed conditional on `time`.
#'
#'   `Sn` and `Fn` are scalar thresholds for this look. Expected success is
#'   declared when the estimated probability of completed-data success among
#'   the current participants is strictly greater than `Sn`. Futility is
#'   declared when the corresponding maximum-sample probability is strictly
#'   less than `Fn`. Set `Fn = 0` or `NULL` to disable futility. Exact
#'   one-sided Monte Carlo bounds are returned as diagnostics and do not drive
#'   either decision.
#'
#'   The function requires treatment assignments to perform the arm-specific
#'   posterior and completed-data analyses. In a blinded trial, an independent
#'   unblinded statistician or service should join the treatment assignments,
#'   run this function, and return the aggregate decision without distributing
#'   the subject-level input.
#'
#' @return An object of class `goldilocks_interim`, containing:
#'
#'   - `decision`: a one-row decision summary;
#'   - `probabilities`: current- and maximum-sample predictive probabilities;
#'   - `monte_carlo`: estimates, standard errors, and diagnostic bounds;
#'   - `diagnostics`: observed status counts, potential accruals, warnings,
#'     imputation diagnostics, and resolved Gamma prior and posterior parameters
#'     by arm and interval;
#'   - `trace`: a one-row decision trace compatible with
#'     [plot_trial_trace()] and [summarise_trial_trace()];
#'   - `metadata`: the evaluated design, resolved prior design, package version,
#'     time-origin, data-cut, and random-number policy. For `method =
#'     "bayes-bin"`, `metadata$design` retains the normalized imputation prior
#'     (`prior_surv`), completed-data analysis prior (`prior_bin`), and
#'     imputation horizon (`end_of_study`).
#'
#' @examples
#' interim_data <- data.frame(
#'   id = 1:6,
#'   treatment = c(0, 1, 0, 1, 0, 1),
#'   enrollment = c(0, 1, 2, 3, 4, 5),
#'   time = c(6, 5, 4, 3, 2, 1),
#'   event = c(1, 0, 0, 1, 0, 0),
#'   status = c("event", "pending", "pending", "event", "pending", "pending")
#' )
#'
#' evaluate_interim(
#'   data = interim_data,
#'   data_cut = 6,
#'   look = 1,
#'   N_total = 10,
#'   end_of_study = 12,
#'   rand_ratio = c(control = 1, treatment = 1),
#'   alternative = "less",
#'   N_impute = 5,
#'   seed = 2026
#' )
#'
#' @export
evaluate_interim <- function(
  data,
  data_cut,
  look,
  N_total,
  end_of_study,
  cutpoints = NULL,
  prior_surv = c(0.1, 0.1),
  prior_bin = c(1, 1),
  bin_method = "mc",
  rand_ratio = c(control = 1, treatment = 1),
  single_arm = FALSE,
  alternative = "greater",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.9,
  prob_ha = 0.95,
  N_impute = 500,
  N_mcmc = 1000,
  mc_conf_level = 0.95,
  empty_interval = c("prior", "propagate", "error"),
  method = "logrank",
  binary_imputation = c("event-time", "bernoulli"),
  seed = NULL
) {
  Call <- match.call()
  caller_rng_kind <- RNGkind()
  empty_interval <- match.arg(empty_interval)
  binary_imputation <- match.arg(binary_imputation)

  validate_positive_integer_scalar(look, "look")
  validate_positive_integer_scalar(N_total, "N_total")
  validate_logical_scalar(single_arm, "single_arm")
  validate_cutpoints(cutpoints)
  validate_endpoint_time(end_of_study, cutpoints, "end_of_study")
  validate_single_probability(prob_ha, "prob_ha")
  validate_positive_integer_scalar(N_impute, "N_impute")
  validate_positive_integer_scalar(N_mcmc, "N_mcmc")
  validate_single_probability(
    mc_conf_level,
    "mc_conf_level",
    upper_open = TRUE
  )
  if (mc_conf_level <= 0.5) {
    stop("'mc_conf_level' must be greater than 0.5 and less than 1")
  }
  validate_analysis_configuration(
    method,
    alternative,
    single_arm,
    imputed_final = FALSE
  )
  validate_h0(h0, method, single_arm)
  if (method == "bayes-bin") {
    validate_bayes_binomial_args(prior_bin, bin_method, N_mcmc)
  }

  Sn <- normalize_interim_threshold(Sn, 1L, "Sn")[[1L]]
  Fn <- normalize_interim_threshold(
    Fn,
    1L,
    "Fn",
    null_disables = TRUE
  )[[1L]]
  check_futility <- Fn != 0

  n_intervals <- length(cutpoints) + 1L
  prior_surv <- normalize_gamma_prior(
    prior_surv,
    n_intervals = n_intervals,
    single_arm = single_arm,
    name = "prior_surv"
  )
  requested_N_mcmc <- N_mcmc
  if (!method %in% c("bayes-surv", "bayes-bin")) {
    N_mcmc <- 1L
  }

  prepared <- prepare_external_interim_data(
    data = data,
    data_cut = data_cut,
    N_total = N_total,
    end_of_study = end_of_study,
    rand_ratio = rand_ratio,
    single_arm = single_arm
  )

  validate_interim_seed(seed)
  if (!is.null(seed)) {
    old_kind <- RNGkind()
    old_seed_exists <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
    if (old_seed_exists) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit(
      {
        do.call(RNGkind, as.list(old_kind))
        if (old_seed_exists) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (
          exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        ) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
    set.seed(seed)
  }

  core <- evaluate_interim_core(
    data_interim = prepared$data,
    look = look,
    planned_N = nrow(data),
    calendar_time = data_cut,
    active_followup = prepared$active_followup,
    end_of_study = end_of_study,
    cutpoints = cutpoints,
    single_arm = single_arm,
    prior_surv = prior_surv,
    prior_bin = prior_bin,
    bin_method = bin_method,
    alternative = alternative,
    h0 = h0,
    Fn = Fn,
    Sn = Sn,
    prob_ha = prob_ha,
    N_impute = N_impute,
    N_mcmc = N_mcmc,
    mc_conf_level = mc_conf_level,
    empty_interval = empty_interval,
    method = method,
    binary_imputation = binary_imputation,
    check_futility = check_futility
  )

  probabilities <- data.frame(
    estimand = c(
      "success_if_stop_now",
      if (check_futility) "success_at_maximum" else character()
    ),
    probability = core$monte_carlo$estimate,
    threshold = core$monte_carlo$threshold,
    direction = core$monte_carlo$direction,
    threshold_crossed = core$monte_carlo$point_crossed,
    stringsAsFactors = FALSE
  )
  decision <- core$trace[
    c(
      "look",
      "planned_N",
      "calendar_time",
      "decision",
      "decision_reason",
      "ppp_stop_now",
      "success_threshold",
      "ppp_success_at_max",
      "futility_threshold"
    )
  ]

  diagnostics <- c(
    list(
      status_counts = prepared$status_counts,
      status_counts_by_arm = prepared$status_counts_by_arm,
      participant_counts = prepared$participant_counts,
      target_allocation = prepared$target_allocation,
      current_allocation = prepared$current_allocation,
      potential_accruals = prepared$potential_accruals
    ),
    core$diagnostics
  )
  design <- list(
    look = as.integer(look),
    N_total = as.integer(N_total),
    end_of_study = end_of_study,
    cutpoints = cutpoints,
    prior_surv = prior_surv,
    prior_bin = prior_bin,
    bin_method = bin_method,
    rand_ratio = prepared$rand_ratio,
    single_arm = single_arm,
    alternative = alternative,
    h0 = h0,
    Fn = Fn,
    Sn = Sn,
    prob_ha = prob_ha,
    N_impute = as.integer(N_impute),
    N_mcmc = as.integer(requested_N_mcmc),
    mc_conf_level = mc_conf_level,
    empty_interval = empty_interval,
    method = method,
    binary_imputation = binary_imputation
  )

  out <- list(
    decision = decision,
    probabilities = probabilities,
    monte_carlo = core$monte_carlo,
    diagnostics = diagnostics,
    trace = core$trace,
    metadata = list(
      package_version = as.character(utils::packageVersion("goldilocks")),
      call = Call,
      design = design,
      prior_design = core$diagnostics$prior,
      data_cut = data_cut,
      time_origin = "first participant randomization (time 0)",
      rng = list(
        kind = caller_rng_kind,
        seed_policy = if (is.null(seed)) {
          "caller_state"
        } else {
          "explicit_preserve_caller"
        },
        seed = seed
      )
    )
  )
  class(out) <- "goldilocks_interim"
  out
}

#' Prepare externally observed data for predictive imputation
#'
#' @keywords internal
#' @noRd
prepare_external_interim_data <- function(
  data,
  data_cut,
  N_total,
  end_of_study,
  rand_ratio,
  single_arm
) {
  required_columns <- c(
    "id",
    "treatment",
    "enrollment",
    "time",
    "event",
    "status"
  )
  missing_columns <- setdiff(required_columns, names(data))
  if (!is.data.frame(data) || length(missing_columns) > 0L) {
    stop(
      "'data' must be a data frame containing columns: ",
      paste(required_columns, collapse = ", ")
    )
  }
  if (nrow(data) == 0L) {
    stop("'data' must contain at least one enrolled subject")
  }
  if (nrow(data) > N_total) {
    stop("The number of rows in 'data' must not exceed 'N_total'")
  }
  if (!is.atomic(data$id) || anyNA(data$id) || anyDuplicated(data$id)) {
    stop("'data$id' must contain unique, non-missing identifiers")
  }
  validate_binary_vector(data$treatment, "data$treatment")
  validate_binary_vector(data$event, "data$event")
  validate_nonnegative_numeric_vector(data$enrollment, "data$enrollment")
  validate_nonnegative_numeric_vector(data$time, "data$time")
  if (
    length(data_cut) != 1L ||
      !is.numeric(data_cut) ||
      is.na(data_cut) ||
      !is.finite(data_cut) ||
      data_cut < 0
  ) {
    stop("'data_cut' must be a single finite non-negative value")
  }

  tolerance <- sqrt(.Machine$double.eps) *
    max(
      1,
      abs(data_cut),
      abs(end_of_study),
      data$enrollment,
      data$time
    )
  if (abs(min(data$enrollment)) > tolerance) {
    stop(
      "'data$enrollment' must use first participant randomization as time 0"
    )
  }
  if (any(data$enrollment > data_cut + tolerance)) {
    stop("'data$enrollment' cannot occur after 'data_cut'")
  }
  available_followup <- pmin(data_cut - data$enrollment, end_of_study)
  if (any(data$time > available_followup + tolerance)) {
    stop(
      "'data$time' cannot exceed follow-up available at 'data_cut' or ",
      "'end_of_study'"
    )
  }

  status <- as.character(data$status)
  allowed_status <- c("event", "complete", "pending", "censored")
  if (
    length(status) != nrow(data) ||
      anyNA(status) ||
      any(!status %in% allowed_status)
  ) {
    stop(
      "'data$status' must contain only: ",
      paste(allowed_status, collapse = ", ")
    )
  }
  if (any((status == "event") != (data$event == 1))) {
    stop("Only rows with 'status = \"event\"' may have 'event = 1'")
  }
  complete <- status == "complete"
  if (any(complete & abs(data$time - end_of_study) > tolerance)) {
    stop(
      "Rows with 'status = \"complete\"' must have 'time = end_of_study'"
    )
  }
  incomplete <- status %in% c("pending", "censored")
  if (any(incomplete & data$time >= end_of_study - tolerance)) {
    stop(
      "Rows with pending or censored status must have 'time' below ",
      "'end_of_study'"
    )
  }

  treatment <- as.integer(data$treatment)
  event <- as.numeric(data$event)
  if (single_arm) {
    if (any(treatment != 1L)) {
      stop(
        "A single-arm interim data cut must have 'treatment = 1' in every row"
      )
    }
    rand_ratio <- c(treatment = 1L)
    target_allocation <- c(treatment = as.integer(N_total))
    current_allocation <- c(treatment = sum(treatment == 1L))
  } else {
    rand_ratio <- normalize_allocation(rand_ratio, "rand_ratio")
    if (!all(c(0L, 1L) %in% treatment)) {
      stop("A two-arm interim data cut must contain subjects in both arms")
    }
    target_raw <- N_total * rand_ratio / sum(rand_ratio)
    target_tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(target_raw))
    if (any(abs(target_raw - round(target_raw)) > target_tolerance)) {
      stop(
        "'N_total' must divide exactly according to 'rand_ratio'; observed ",
        "targets: control=",
        format(target_raw[["control"]]),
        ", treatment=",
        format(target_raw[["treatment"]])
      )
    }
    target_allocation <- as.integer(round(target_raw))
    names(target_allocation) <- names(rand_ratio)
    current_allocation <- c(
      control = sum(treatment == 0L),
      treatment = sum(treatment == 1L)
    )
  }
  potential_accruals <- target_allocation - current_allocation
  if (any(potential_accruals < 0L)) {
    offending <- names(potential_accruals)[potential_accruals < 0L]
    stop(
      "Observed enrollment exceeds the maximum implied by 'N_total' and ",
      "'rand_ratio' for arm(s): ",
      paste(offending, collapse = ", ")
    )
  }

  observed_internal <- data.frame(
    time = pmax(data$time, .Machine$double.eps),
    event = event,
    treatment = treatment,
    subject_enrolled = TRUE,
    subject_impute_success = incomplete,
    subject_impute_futility = FALSE
  )
  future_treatment <- if (single_arm) {
    rep.int(1L, potential_accruals[["treatment"]])
  } else {
    c(
      rep.int(0L, potential_accruals[["control"]]),
      rep.int(1L, potential_accruals[["treatment"]])
    )
  }
  future_internal <- data.frame(
    time = rep.int(0, length(future_treatment)),
    event = rep.int(0L, length(future_treatment)),
    treatment = future_treatment,
    subject_enrolled = rep.int(FALSE, length(future_treatment)),
    subject_impute_success = rep.int(FALSE, length(future_treatment)),
    subject_impute_futility = rep.int(TRUE, length(future_treatment))
  )
  internal <- rbind(observed_internal, future_internal)

  status_factor <- factor(status, levels = allowed_status)
  status_counts <- data.frame(
    status = allowed_status,
    count = as.integer(table(status_factor)),
    stringsAsFactors = FALSE
  )
  arm <- ifelse(treatment == 1L, "treatment", "control")
  status_counts_by_arm <- as.data.frame(
    table(
      arm = factor(
        arm,
        levels = if (single_arm) "treatment" else c("control", "treatment")
      ),
      status = status_factor
    ),
    stringsAsFactors = FALSE
  )
  names(status_counts_by_arm)[[3L]] <- "count"
  participant_counts <- rbind(
    status_counts,
    data.frame(
      status = "not_enrolled",
      count = as.integer(sum(potential_accruals)),
      stringsAsFactors = FALSE
    )
  )

  list(
    data = internal,
    active_followup = sum(status == "pending"),
    status_counts = status_counts,
    status_counts_by_arm = status_counts_by_arm,
    participant_counts = participant_counts,
    rand_ratio = rand_ratio,
    target_allocation = target_allocation,
    current_allocation = current_allocation,
    potential_accruals = potential_accruals
  )
}

#' Validate an interim-analysis seed
#'
#' @keywords internal
#' @noRd
validate_interim_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(TRUE))
  }
  if (
    length(seed) != 1L ||
      !is.numeric(seed) ||
      is.na(seed) ||
      !is.finite(seed) ||
      seed < 0 ||
      seed > .Machine$integer.max ||
      seed != floor(seed)
  ) {
    stop(
      "'seed' must be NULL or a single non-negative integer no greater than ",
      ".Machine$integer.max"
    )
  }
  invisible(TRUE)
}

#' Print an externally evaluated interim analysis
#'
#' @param x A `goldilocks_interim` object returned by [evaluate_interim()].
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.goldilocks_interim <- function(x, ...) {
  decision <- x$decision[1L, , drop = FALSE]
  cat("Goldilocks interim evaluation\n")
  cat("Look: ", decision$look, " (N = ", decision$planned_N, ")\n", sep = "")
  cat("Decision: ", decision$decision, "\n", sep = "")
  cat(
    "Predictive success if accrual stops now: ",
    format(decision$ppp_stop_now, digits = 4L),
    "\n",
    sep = ""
  )
  if (is.finite(decision$ppp_success_at_max)) {
    cat(
      "Predictive success at maximum sample size: ",
      format(decision$ppp_success_at_max, digits = 4L),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}
