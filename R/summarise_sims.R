#' @title Summarize simulations to get operating characteristics
#'
#' @param data list (of data frames) or a single data frame. If summarizing a
#'   single run of simulations, \code{data} will be a \code{data.frame} object
#'   returned from \code{\link{survival_adapt}}. If summarizing multiple
#'   simulation scenarios, \code{data} will be a \code{list} object, with each
#'   element being a \code{data.frame} object.
#'
#' @return Data frame reporting the operating characteristics.
#' @export
summarise_sims <- function(data) {

  if (is.list(data)) {
    multi <- TRUE
    if (!all(sapply(data, is.data.frame))) {
      stop("Each element of list must be a data.frame")
    }
  } else if (!is.data.frame(results)) {
    multi <- FALSE
    stop("Input data should be a data.frame")
  }

  out <- data %>%
    summarise(
      "type_2_error"  = mean(!stop_futility & post_prob_ha > prob_threshold),
      "stop_success"  = mean(stop_expected_success),
      "stop_futility" = mean(stop_futility),
      "stop_max_N"    = mean(!(stop_expected_success | stop_futility)),
      "mean_N"        = mean(N_enrolled),
      "sd_N"          = sd(N_enrolled),
      "stop_and_fail" = mean(stop_expected_success & post_prob_ha <= prob_threshold)
    )

  return(out)

}
