#' @title Summarize simulations to get operating characteristics
#'
#' @param data list (of data frames) or a single data frame. If summarizing a
#'   single run of simulations, \code{data} will be a \code{data.frame} object
#'   returned from \code{\link{survival_adapt}}. If summarizing multiple
#'   simulation scenarios, \code{data} will be a \code{list} object, with each
#'   element being a \code{data.frame} object.
#'
#' @return Data frame reporting the operating characteristics, including the
#'   power (which will be equal to the type I error in the null case); the
#'   proportion of trials that stopped for early expected success, futility, or
#'   went to the maximum sample size. The average stopping sample size (and
#'   standard deviation) are also recorded. The proportion of trials that
#'   stopped early for expected success, yet went to ultimately fail are also
#'   reported.
#' @export
summarise_sims <- function(data) {

  if (inherits(data, "list")) {
    fnames <- names(data)
    if (is.null(fnames)) {
      fnames <- 1:length(data)
    }
    for (i in 1:length(data)) {
      data[[i]]$scenario <- fnames[i]
      data[[i]]
    }
    data <- do.call("rbind", data)
  } else {
    data$scenario <- 1
  }

  out <- data %>%
    group_by(.data$scenario) %>%
    summarise(
      "power"         = mean(!.data$stop_futility & .data$post_prob_ha > .data$prob_threshold),
      "stop_success"  = mean(.data$stop_expected_success),
      "stop_futility" = mean(.data$stop_futility),
      "stop_max_N"    = 1 - mean(.data$stop_success) - mean(.data$stop_futility),
      "mean_N"        = mean(.data$N_enrolled),
      "sd_N"          = sd(.data$N_enrolled),
      "stop_and_fail" = mean(.data$stop_expected_success & .data$post_prob_ha <= .data$prob_threshold)
    )

  return(out)

}
