#' Calculate a log-rank result from outcome vectors
#'
#' @param time A numeric vector of follow-up times.
#' @param event A binary vector of event indicators.
#' @param treatment A binary vector of treatment assignments.
#' @param alternative A character value specifying the direction of the
#'   alternative hypothesis.
#'
#' @return A list containing `1 - P` for the log-rank test and an unavailable
#'   treatment-effect estimate.
#'
#' @keywords internal
#' @noRd
analyse_logrank <- function(time, event, treatment, alternative) {
  control <- treatment == 0
  lr <- logrank_test(
    groupa = time[control],
    groupb = time[!control],
    groupacensored = event[control],
    groupbcensored = event[!control]
  )
  assert_logrank_estimable(lr)

  if (alternative == "two.sided") {
    success <- 1 - lr[3]
  } else {
    # Log-rank z > 0 when control has excess events (treatment beneficial).
    # This is opposite to the Cox convention.
    # "less" => treatment beneficial => large success when z >> 0
    # "greater" => treatment harmful => large success when z << 0
    z <- lr[2]
    if (alternative == "less") {
      success <- pnorm(z)
    } else if (alternative == "greater") {
      success <- 1 - pnorm(z)
    }
  }

  list(success = success, effect = NA)
}

#' @title Calculate a two-sample log-rank test
#'
#' @description Calculates a log-rank test from the follow-up times and event
#'   indicators in two independent treatment groups.
#'
#' @param groupa A numeric vector of non-negative follow-up times for group A.
#' @param groupb A numeric vector of non-negative follow-up times for group B.
#' @param groupacensored A zero-one integer vector indicating events in group A,
#'   where `1` denotes an event and `0` denotes censoring.
#' @param groupbcensored A zero-one integer vector indicating events in group B,
#'   where `1` denotes an event and `0` denotes censoring.
#' @param onlyz A single logical value indicating whether to return only the
#'   standardized log-rank statistic. The default is `FALSE`.
#'
#' @return A numeric vector containing the chi-squared statistic, standardized
#'   statistic, and two-sided *P*-value. When `onlyz = TRUE`, only the
#'   standardized statistic is returned.
#'
#' @examples
#' T1 <- c(6, 6, 6, 6, 7, 9, 10, 10, 11, 13, 16, 17, 19, 20, 22, 23, 25, 32, 32, 34, 35)
#' E1 <- c(1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0)
#' T2 <- c(1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22, 23)
#' E2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' logrank_test(T1, T2, E1, E2)
#' #1.679294e+01 -4.097919e+00, 4.168809e-05
#'
#' @noRd
#' @keywords internal
logrank_test <- function(
  groupa,
  groupb,
  groupacensored,
  groupbcensored,
  onlyz = FALSE
) {
  logrank_instance(groupa, groupb, groupacensored, groupbcensored, onlyz)
}

#' @title Assert that a log-rank result is estimable
#'
#' @description Stops with a diagnostic error when a log-rank statistic cannot
#'   be computed as a finite comparison between treatment groups.
#'
#' @noRd
assert_logrank_estimable <- function(lr) {
  if (
    length(lr) < 3 ||
      any(is.na(lr[1:3])) ||
      any(!is.finite(lr[1:3]))
  ) {
    stop(
      "Log-rank analysis is non-estimable: the test statistic is not finite. ",
      "This can occur when there are no events or insufficient information ",
      "to compare treatment groups."
    )
  }

  invisible(TRUE)
}
