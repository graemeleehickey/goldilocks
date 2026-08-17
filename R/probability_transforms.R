#' @title Convert cumulative hazards to event probabilities
#'
#' @description Uses `expm1()` to retain precision for cumulative hazards near
#'   zero. Zero maps to zero, positive infinity maps to one, and missing values
#'   propagate without warnings.
#'
#' @noRd
cumulative_hazard_to_probability <- function(cumulative_hazard) {
  -expm1(-cumulative_hazard)
}

#' @title Convert event probabilities to cumulative hazards
#'
#' @description Uses `log1p()` to retain precision for probabilities near
#'   zero. Zero maps to zero, one maps to positive infinity, and missing values
#'   propagate without warnings.
#'
#' @noRd
probability_to_cumulative_hazard <- function(probability) {
  -log1p(-probability)
}
