#' Locate times using the survival counting-process convention
#'
#' @description Assigns positive times to intervals that are open on the left
#'   and closed on the right. A value equal to an internal boundary is therefore
#'   assigned to the interval ending at that boundary, as in [survival::Surv()]
#'   counting-process data and [survival::survSplit()]. Time zero returns index
#'   zero because it does not belong to a positive-length interval under this
#'   convention; callers with a fixed origin can handle it separately.
#'
#' @param x A numeric vector of finite, non-negative times to assign to
#'   intervals.
#' @param interval_starts A numeric vector of finite, strictly increasing
#'   interval start times beginning at zero.
#'
#' @return An integer vector of interval indices.
#'
#' @keywords internal
#' @noRd
right_closed_interval_index <- function(x, interval_starts) {
  findInterval(x, interval_starts, left.open = TRUE)
}
