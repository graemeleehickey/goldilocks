#' @title Convert cumulative hazards to event probabilities
#'
#' @description Uses `expm1()` to retain precision for cumulative hazards near
#'   zero. Zero maps to zero, positive infinity maps to one, and missing values
#'   propagate without warnings.
#'
#' @param cumulative_hazard A numeric vector of non-negative cumulative hazards.
#'
#' @return A numeric vector of event probabilities in `[0, 1]`.
#'
#' @noRd
cum_haz_to_prob <- function(cumulative_hazard) {
  -expm1(-cumulative_hazard)
}

#' @title Convert event probabilities to cumulative hazards
#'
#' @description Uses `log1p()` to retain precision for probabilities near zero.
#'   Zero maps to zero, one maps to positive infinity, and missing values
#'   propagate without warnings.
#'
#' @param probability A numeric vector of event probabilities in `[0, 1]`.
#'
#' @return A numeric vector of non-negative cumulative hazards.
#'
#' @noRd
prob_to_cum_haz <- function(probability) {
  -log1p(-probability)
}
