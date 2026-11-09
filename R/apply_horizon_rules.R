#' Apply a prespecified endpoint schedule to saved trial paths
#'
#' @description Applies endpoint horizons and stopping thresholds without further imputation, fitting, or random-number generation.
#' @param paths A complete `goldilocks_horizon_paths` object from [sim_horizon_paths()].
#' @param endpoint_schedule Data frame with exactly two columns, `N` and `horizon`, containing every planned interim enrollment size in increasing order and a nondecreasing endpoint horizon at each look. Horizons must be available in the saved paths. Use a zero-row data frame with numeric columns when there are no interims.
#' @inheritParams apply_stopping_rules
#' @param Fn Futility threshold: one probability repeated at every interim or one per scheduled look. `NULL` or zero disables futility; the default is `0.05`.
#' @param Sn Expected-success threshold: one probability repeated at every interim or one per scheduled look. The default is `0.9`; one disables expected-success stopping.
#' @param Qn Immediate-success threshold: one probability repeated at every interim or one per scheduled look. It must be at least `Sn` at every look. The default, one, disables immediate-success stopping.
#' @param return_trace Logical; whether to include a `traces` data frame containing the reached interim decisions in the returned list. Defaults to `FALSE`.
#' @details Success uses the scheduled endpoint and its own prediction data cut. Futility always uses the final endpoint specified when simulating the paths and its separate data cut, regardless of the scheduled success endpoint. Immediate success, expected success, and binding futility retain their existing priority and strict threshold comparisons. The first stop is terminal. An expected-success stop fixes the endpoint for that cohort's final analysis; failure there does not trigger extension.
#'
#' Only calculations needed by a candidate determine whether its trial is evaluable. A failed maximum-\eqn{N} counterfactual makes that comparison unavailable without discarding a known adaptive outcome. Both the ordinary futility probability and joint futility with counterfactual success are available through [summarise_sims()]. The latter always refers to maximum \eqn{N} at `final_horizon`.
#'
#' Under a generating scenario where every endpoint null hypothesis is true, the fraction of trials declaring success at any selected endpoint estimates the overall type I error. Calibrate the schedule and thresholds together. If some endpoint hypotheses are false, overall success includes legitimate findings; assessing false positives then requires identifying success at a truly null endpoint. A success at an earlier endpoint does not by itself reject the final-horizon null.
#'
#' @return A list compatible with [summarise_sims()] and [summarise_calendar_time()], containing `sims`, `failures`, and `call`, and optionally `traces`. Simulation rows additionally identify `endpoint_horizon`, `final_horizon`, and `followup_horizon`, and the null `margin` used for the selected endpoint. Treatment-effect estimates from different endpoint horizons should be reported separately. Overall success is success under the complete prespecified adaptive design, which may declare success at different endpoints.
#'
#' Calendar columns `analysis_ready_time`, `planned_completion_time`, and `followup_person_time` refer to the selected endpoint. Separate `collection_ready_time`, `planned_collection_time`, and `collection_person_time` refer to continued collection through `followup_horizon` among the enrolled cohort, even after a stop. They can differ after a 12-month stop even when `followup_horizon = final_horizon = 18`. `decision_time` equals the interim date for immediate success or futility and `analysis_ready_time` for an expected-success or maximum-N final analysis. An expected-success stop therefore ends recruitment at `accrual_stop_time` but determines final success or failure later. Event and dropout times can precede planned completion. Trace columns identify both prediction data horizons and futility event and pending-outcome counts; the unprefixed counts refer to success prediction. Saved posterior diagnostics remain separate by horizon.
#' @seealso [sim_horizon_paths()], [sim_horizon_trials()]
#' @export
apply_horizon_rules <- function(
  paths,
  endpoint_schedule,
  Fn = 0.05,
  Sn = 0.9,
  Qn = 1,
  prob_ha = 0.95,
  return_trace = FALSE
) {
  Call <- match.call()
  if (is.function(Call[[1L]])) {
    Call[[1L]] <- quote(apply_horizon_rules)
  }
  Call$paths <- quote(paths)
  validate_horizon_paths(paths)
  args <- attr(paths, "arguments", exact = TRUE)
  schedule <- validate_endpoint_schedule(endpoint_schedule, args)
  projected <- project_horizon_paths(paths, schedule)
  out <- apply_stopping_rules(projected, Fn, Sn, Qn, prob_ha, return_trace)
  out$call <- Call
  selected_horizons <- c(schedule$horizon, args$final_horizon)
  out$sims$endpoint_horizon <- selected_horizons[match(
    out$sims$N_enrolled,
    c(args$interim_look, args$N_total)
  )]
  out$sims$final_horizon <- args$final_horizon
  out$sims$followup_horizon <- args$followup_horizon
  out$sims$margin <- args$horizon_settings$h0[match(
    out$sims$endpoint_horizon,
    args$horizon_settings$horizon
  )]
  collection <- dplyr::bind_rows(lapply(projected$trials, function(path) {
    if (is.null(path$finals)) {
      return(NULL)
    }
    data.frame(
      trial = path$trial,
      path$finals[c(
        "planned_N",
        "collection_ready_time",
        "planned_collection_time",
        "collection_person_time"
      )]
    )
  }))
  out$sims <- dplyr::left_join(
    out$sims,
    collection,
    by = c("trial", "N_enrolled" = "planned_N")
  )
  attr(out, "endpoint_schedule") <- schedule
  attr(out, "arguments") <- c(
    args,
    attr(out, "stopping_rules", exact = TRUE),
    list(endpoint_schedule = schedule, return_trace = return_trace)
  )
  attr(out, "path_metadata")$arguments <- args
  attr(out, "path_metadata")$bank_type <- "goldilocks_horizon_paths"
  out
}

validate_horizon_paths <- function(paths) {
  if (
    !inherits(paths, "goldilocks_horizon_paths") ||
      !identical(paths$schema_version, 1L)
  ) {
    stop(
      "'paths' must be a supported goldilocks_horizon_paths object (schema version 1)"
    )
  }
  args <- attr(paths, "arguments", exact = TRUE)
  if (
    !is.list(args) ||
      !all(names(formals(sim_horizon_paths)) %in% names(args)) ||
      !is.list(paths$trials) ||
      length(paths$trials) != args$N_trials ||
      !is.data.frame(paths$failures)
  ) {
    stop("Incomplete horizon paths or metadata")
  }
  args <- validate_horizon_arguments(args)
  horizons <- args$endpoint_horizons
  expected_N <- rep(c(args$interim_look, args$N_total), each = length(horizons))
  expected_H <- rep(horizons, length(args$interim_look) + 1L)
  for (path in paths$trials) {
    if (length(path$interims) != length(args$interim_look)) {
      stop("Incomplete interim schedule in horizon paths")
    }
    if (
      !is.null(path$finals) &&
        (!identical(
          as.numeric(path$finals$planned_N),
          as.numeric(expected_N)
        ) ||
          !identical(as.numeric(path$finals$horizon), as.numeric(expected_H)))
    ) {
      stop("Inconsistent final cohort or horizon schedule in saved paths")
    }
    for (look in path$interims) {
      if (is.null(look)) {
        next
      }
      if (length(look) != length(horizons)) {
        stop("Incomplete endpoint horizons at an interim look")
      }
      for (evidence in look) {
        if (is.null(evidence)) {
          next
        }
        for (scores in evidence$scores) {
          if (
            !is.null(scores) &&
              is.null(attr(scores, "failure", exact = TRUE)) &&
              length(scores) != args$N_impute
          ) {
            stop(
              "Incomplete predictive replicates in horizon paths; no replicates may be dropped"
            )
          }
        }
      }
    }
  }
  invisible(TRUE)
}

# Project the chosen horizon at each look into the existing stopping-rule
# evaluator, keeping its decision priority and candidate-specific failures.
project_horizon_paths <- function(paths, schedule) {
  args <- validate_horizon_arguments(attr(paths, "arguments", exact = TRUE))
  horizons <- args$endpoint_horizons
  final_index <- match(args$final_horizon, horizons)
  selected <- match(c(schedule$horizon, args$final_horizon), horizons)
  out <- paths
  attr(out, "arguments") <- args
  class(out) <- "goldilocks_paths"
  out$trials <- lapply(paths$trials, function(path) {
    path$interims <- lapply(seq_along(args$interim_look), function(i) {
      evidence <- path$interims[[i]]
      if (is.null(evidence)) {
        return(NULL)
      }
      combine_horizon_evidence(
        evidence[[selected[i]]],
        evidence[[final_index]],
        schedule$horizon[i],
        args$final_horizon
      )
    })
    if (!is.null(path$finals)) {
      rows <- (seq_along(selected) - 1L) * length(horizons) + selected
      path$finals <- path$finals[rows, , drop = FALSE]
    }
    path
  })
  out
}
