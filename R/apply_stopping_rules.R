#' Apply stopping rules to a saved trial-path bank
#'
#' @description Reconstructs adaptive trial outcomes from [sim_trial_paths()]
#'   without predictive imputation, model fitting, or random-number generation.
#'
#' @param paths A complete `goldilocks_paths` object from [sim_trial_paths()].
#' @inheritParams survival_adapt
#' @inheritParams sim_trials
#'
#' @details At each look, predictive success is the fraction of saved scores
#'   strictly greater than `prob_ha`. Decisions retain the ordered immediate
#'   success (`Qn`), expected success (`Sn`), and binding futility (`Fn`) rules.
#'   The first stopping decision is terminal. Expected success requires a final
#'   analysis of the stopping cohort; immediate success and futility do not.
#'
#'   Only failures needed by the candidate invalidate a trial. Failures at
#'   unreached looks, unused maximum-cohort predictions, or optional final
#'   analyses remain in the bank but do not remove an evaluable trial. A missing
#'   maximum-N final analysis makes the counterfactual comparison unavailable,
#'   without invalidating an otherwise known adaptive outcome.
#'
#'   Candidate results from the same bank are paired by `trial`. Use paired
#'   uncertainty calculations for differences; marginal Monte Carlo intervals
#'   do not describe uncertainty in a difference. Screen candidate rules on one
#'   bank and validate selected rules on a fresh independent bank.
#'
#' @return A list compatible with [summarise_sims()], [summarise_calendar_time()],
#'   and simulation plots, containing `sims`, `failures`, and `call`, plus
#'   `traces` when requested. Each simulation row retains its original `trial`
#'   identifier and additionally reports `success_at_max` (the final outcome
#'   if enrollment continued to maximum N without interim stopping) and
#'   `futility_and_success_at_max` (binding futility together with that
#'   counterfactual success). The latter is `FALSE` for non-futility trials and
#'   `NA` for futility trials whose maximum-N final analysis failed. These
#'   quantities use the candidate's `prob_ha` and the bank's final analysis.
#'
#'   `summarise_sims()` reports the joint probability as
#'   `futility_and_success_at_max`, with its available and unavailable counts.
#'   Its ordinary `stop_futility` estimate is the futility probability under
#'   the simulated scenario; under an alternative this is the requested
#'   alternative-scenario futility rate. The joint probability is not the net
#'   power loss relative to a fixed maximum-N design.
#'
#'   Attributes include the evaluated design, bank provenance (`path_metadata`),
#'   and `stopping_rules`. The bank itself is not copied into the result.
#'   If no trial is evaluable, the function stops with a
#'   `goldilocks_all_trials_failed` error carrying the candidate failure table.
#'
#' @seealso [sim_trial_paths()], [summarise_sims()]
#' @export
#' @examples
#' paths <- sim_trial_paths(
#'   hazard_treatment = 0.05, hazard_control = 0.08,
#'   N_total = 60, interim_look = c(20, 40), end_of_study = 12,
#'   lambda = 5, method = "bayes-bin", bin_method = "quadrature",
#'   alternative = "less", N_impute = 5, N_trials = 2, seed = 123
#' )
#' candidates <- list(
#'   original = list(Fn = 0.05, Sn = 0.90, Qn = 1, prob_ha = 0.975),
#'   stricter = list(Fn = 0.10, Sn = 0.95, Qn = 1, prob_ha = 0.99)
#' )
#' results <- lapply(candidates, function(rule) {
#'   do.call(apply_stopping_rules, c(list(paths = paths), rule))
#' })
#' summarise_sims(results)
apply_stopping_rules <- function(
  paths,
  Fn = 0.05,
  Sn = 0.9,
  Qn = 1,
  prob_ha = 0.95,
  return_trace = FALSE
) {
  Call <- match.call()
  # do.call() may otherwise embed the whole bank in the stored call.
  Call$paths <- quote(paths)
  validate_trial_paths(paths)
  validate_single_probability(prob_ha, "prob_ha")
  validate_logical_scalar(return_trace, "return_trace")
  args <- attr(paths, "arguments", exact = TRUE)
  n_looks <- length(args$interim_look)
  Fn <- normalize_interim_threshold(Fn, n_looks, "Fn", null_disables = TRUE)
  Sn <- normalize_interim_threshold(Sn, n_looks, "Sn")
  Qn <- normalize_interim_threshold(Qn, n_looks, "Qn")
  validate_success_threshold_order(Sn, Qn)
  rules <- list(Fn = Fn, Sn = Sn, Qn = Qn, prob_ha = prob_ha)
  evaluated <- lapply(paths$trials, function(path) {
    tryCatch(
      replay_trial_path(path, args, rules, return_trace),
      error = function(error) {
        list(
          failure = data.frame(
            trial = path$trial,
            error_class = class(error)[1L],
            message = conditionMessage(error)
          )
        )
      }
    )
  })
  failed <- vapply(evaluated, function(x) !is.null(x$failure), logical(1))
  failures <- if (any(failed)) {
    dplyr::bind_rows(lapply(evaluated[failed], `[[`, "failure"))
  } else {
    data.frame(
      trial = integer(),
      error_class = character(),
      message = character()
    )
  }
  if (all(failed)) {
    error <- simpleError(paste0(
      "All ",
      length(evaluated),
      " trials failed under these stopping rules. First error: ",
      failures$message[1L]
    ))
    error$failures <- failures
    class(error) <- c("goldilocks_all_trials_failed", class(error))
    stop(error)
  }
  out <- list(
    sims = dplyr::bind_rows(lapply(evaluated[!failed], `[[`, "summary")),
    failures = failures,
    call = Call
  )
  if (return_trace) {
    out$traces <- dplyr::bind_rows(lapply(evaluated[!failed], `[[`, "trace"))
  }
  for (name in c(
    "enrollment_design",
    "prior_design",
    "rng_metadata",
    "parallel_metadata"
  )) {
    attr(out, name) <- attr(paths, name, exact = TRUE)
  }
  attr(out, "arguments") <- c(args, rules, list(return_trace = return_trace))
  attr(out, "decision_design") <- c(
    list(interim_look = args$interim_look),
    rules[c("Fn", "Sn", "Qn")],
    args[c("N_impute", "N_mcmc", "mc_conf_level")]
  )
  attr(out, "stopping_rules") <- rules
  attr(out, "path_metadata") <- list(
    schema_version = paths$schema_version,
    runtime = attr(paths, "runtime_metadata", exact = TRUE),
    arguments = args,
    trial_ids = vapply(paths$trials, `[[`, integer(1), "trial"),
    failed_calculations = nrow(paths$failures)
  )
  if (nrow(failures)) {
    warning(
      nrow(failures),
      " of ",
      length(evaluated),
      " trials could not be evaluated under these rules and were excluded. See 'result$failures'.",
      call. = FALSE
    )
  }
  out
}

validate_trial_paths <- function(paths) {
  if (
    !inherits(paths, "goldilocks_paths") || !identical(paths$schema_version, 1L)
  ) {
    stop("'paths' must be a supported goldilocks_paths bank (schema version 1)")
  }
  args <- attr(paths, "arguments", exact = TRUE)
  if (
    !is.list(args) ||
      !all(names(formals(sim_trial_paths)) %in% names(args)) ||
      !is.list(paths$trials) ||
      length(paths$trials) != args$N_trials ||
      !is.data.frame(paths$failures)
  ) {
    stop("Incomplete trial-path bank or metadata")
  }
  ids <- vapply(paths$trials, function(x) as.integer(x$trial), integer(1))
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop("Trial-path identifiers must be unique and non-missing")
  }
  for (path in paths$trials) {
    if (length(path$interims) != length(args$interim_look)) {
      stop("Incomplete interim look schedule in bank")
    }
    if (
      !is.null(path$finals) &&
        !identical(
          as.numeric(path$finals$planned_N),
          as.numeric(c(args$interim_look, args$N_total))
        )
    ) {
      stop("Inconsistent final cohort schedule in bank")
    }
    for (evidence in path$interims) {
      if (!is.null(evidence)) {
        for (scores in evidence$scores) {
          if (
            is.null(attr(scores, "failure", exact = TRUE)) &&
              length(scores) != args$N_impute
          ) {
            stop(
              "Incomplete predictive replicates in bank; no replicates may be dropped"
            )
          }
        }
      }
    }
  }
  invisible(TRUE)
}

require_path_scores <- function(scores, stage) {
  failure <- attr(scores, "failure", exact = TRUE)
  if (!is.null(failure)) {
    stop(stage, ": ", failure$message, call. = FALSE)
  }
  if (
    !length(scores) || any(!is.finite(scores)) || any(scores < 0 | scores > 1)
  ) {
    stop(stage, ": predictive scores are unavailable", call. = FALSE)
  }
  invisible(TRUE)
}

replay_trial_path <- function(path, args, rules, return_trace) {
  if (is.null(path$finals)) {
    stop(
      "Trial data generation failed; see the bank failure table",
      call. = FALSE
    )
  }
  stop_index <- nrow(path$finals)
  decision <- "continue"
  ppp_success <- NA_real_
  interim_time <- NA_real_
  trace_rows <- vector("list", length(path$interims))
  check_futility <- any(rules$Fn != 0)
  for (i in seq_along(path$interims)) {
    evidence <- path$interims[[i]]
    if (is.null(evidence)) {
      stop(
        "Interim calculation failed at look ",
        i,
        "; see the bank failure table",
        call. = FALSE
      )
    }
    require_path_scores(
      evidence$scores$current,
      paste("Current-cohort prediction at look", i)
    )
    ppp_success <- sum(evidence$scores$current > rules$prob_ha) / args$N_impute
    maximum <- evidence$scores$maximum
    has_maximum <- length(maximum) == args$N_impute &&
      is.null(attr(maximum, "failure", exact = TRUE)) &&
      all(is.finite(maximum)) &&
      all(maximum >= 0 & maximum <= 1)
    # Success takes precedence; a disabled or unused futility calculation is optional.
    if (ppp_success <= rules$Sn[i] && rules$Fn[i] > 0) {
      require_path_scores(
        maximum,
        paste("Maximum-cohort prediction at look", i)
      )
    }
    ppp_maximum <- if (has_maximum && check_futility) {
      sum(maximum > rules$prob_ha) / args$N_impute
    } else {
      NA_real_
    }
    selected <- select_interim_decision(
      ppp_success,
      ppp_maximum,
      rules$Sn[i],
      rules$Qn[i],
      rules$Fn[i],
      check_futility && has_maximum
    )
    decision <- selected$decision
    if (return_trace) {
      interim <- summarise_interim_evidence(
        evidence,
        rules$Fn[i],
        rules$Sn[i],
        rules$Qn[i],
        rules$prob_ha,
        args$mc_conf_level,
        check_futility && has_maximum
      )
      trace_rows[[i]] <- interim$trace
      if (check_futility) trace_rows[[i]]$futility_threshold <- rules$Fn[i]
    }
    if (decision != "continue") {
      stop_index <- i
      interim_time <- evidence$context$calendar_time
      break
    }
  }
  final <- path$finals[stop_index, , drop = FALSE]
  immediate <- decision == "stop_immediate_success"
  expected <- decision == "stop_expected_success"
  futility <- decision == "stop_futility"
  if (!immediate && !futility && !is.finite(final$post_prob_ha)) {
    stop(
      "Required final analysis failed for cohort N = ",
      final$planned_N,
      call. = FALSE
    )
  }
  trial_success <- if (immediate) {
    TRUE
  } else if (futility) {
    FALSE
  } else {
    final$post_prob_ha > rules$prob_ha
  }
  maximum_score <- path$finals$post_prob_ha[nrow(path$finals)]
  success_at_max <- if (is.finite(maximum_score)) {
    maximum_score > rules$prob_ha
  } else {
    NA
  }
  summary <- data.frame(
    trial = path$trial,
    prob_threshold = rules$prob_ha,
    margin = args$h0,
    alternative = args$alternative,
    N_treatment = final$N_treatment,
    N_control = final$N_control,
    N_enrolled = final$planned_N,
    N_max = args$N_total,
    post_prob_ha = if (immediate) NA_real_ else final$post_prob_ha,
    est_final = if (immediate) NA_real_ else final$est_final,
    ppp_success = ppp_success,
    stop_futility = as.numeric(futility),
    stop_immediate_success = as.numeric(immediate),
    stop_expected_success = as.numeric(expected),
    trial_success = trial_success,
    stopping_reason = calendar_stopping_reason(futility, expected, immediate),
    decision_time = if (immediate || futility) {
      interim_time
    } else {
      final$analysis_ready_time
    },
    final[c(
      "accrual_stop_time",
      "analysis_ready_time",
      "planned_completion_time",
      "followup_person_time",
      "peak_active_followup"
    )],
    success_at_max = success_at_max,
    futility_and_success_at_max = futility & success_at_max
  )
  trace <- if (return_trace) {
    trace <- new_trial_trace(trace_rows)
    trace$trial <- rep.int(path$trial, nrow(trace))
    trace[c("trial", setdiff(names(trace), "trial"))]
  } else {
    NULL
  }
  list(summary = summary, trace = trace)
}
