#' @title Posterior distribution of piecewise exponential constant hazard rates
#'
#' @description Using the Beta-Gamma conjugacy property, the posterior
#'   distribution of the piecewise hazard rates (\eqn{\lambda_j}, for \code{j=1,
#'   \dots, J}) is calculated and sampled from.
#'
#' @inheritParams survival_adapt
#' @inheritParams haz_to_prop
#' @param data data frame. Minimum requirements are 3 columns: event time
#'   (\code{time}), indicator of the event (\code{event}), and indicator for
#'   treatment arm (\code{treatment}). Other columns can be included in the data
#'   frame and will be handled in the split.
#'
#' @return An array of dimension 3. The first dimension is of length
#'   \code{N_mcmc}, the second dimension is of length \eqn{J} (one column for
#'   each hazard piece), and the third dimension is of length 2, with the first
#'   slice including posterior samples from \code{post_treatment}, and the
#'   second slice including posterior samples from \code{post_control}.
#'
#' @importFrom stats rgamma
#' @importFrom dplyr %>% summarise group_by ungroup
#' @importFrom rlang .data
#' @import survival
#'
#' @noRd
posterior <- function(data, cutpoints, prior, N_mcmc, single_arm) {

  n_intervals <- length(cutpoints)

  ik <- cutpoints[-1] # Note: survSplit() doesn't like cuts at 0, so ik = inner knots
  if (length(ik) == 0) {
    ik <- max(data$time)
  }

  data_survsplit <- survSplit(
    Surv(time, event) ~ .,
    data = data,
    cut = ik,
    episode = "interval")

  # In case interim data doesn't span all intervals possible
  data_survsplit$interval <- factor(data_survsplit$interval,
                                    levels = 1:n_intervals)

  data_summ <- data_survsplit %>%
    group_by(.data$treatment, .data$interval, .drop = FALSE) %>%
    summarise(n = length(.data$time),
              tot_time = sum(.data$time - .data$tstart),
              tot_events = sum(.data$event)) %>%
    ungroup()

  # If a time-interval has zero subjects at a given interim analysis, it will
  # mean that there is zero exposure time and events for that stratum. To avoid
  # unrealistic hazard parameter estimates, we propagate the exposure time and
  # event counts from the last non-zero stratum. This is equivalent to
  # independent draws from Gamma posterior from the last non-zero stratum.
  if (any(data_summ$n == 0)) {
    get_i <- which(data_summ$n == 0)
    for (i in get_i) {
      if (i == 1) {
        stop("No subjects in first strata")
      }
      data_summ[i, c("tot_time", "tot_events")] <- data_summ[(i - 1), c("tot_time", "tot_events")]
    }
  }

  post <- array(dim = c(N_mcmc, n_intervals, 2))
  post_treatment <- matrix(nrow = N_mcmc, ncol = n_intervals)
  post_control <- matrix(nrow = N_mcmc, ncol = n_intervals)

  for (j in 1:n_intervals) {
    post[, j, 1] <- with(
      subset(data_summ, treatment == 1),
      rgamma(N_mcmc, prior[1] + tot_events[j], prior[2] + tot_time[j])
    )
  }

  if (!single_arm) { # If control patients present
    for (j in 1:n_intervals) {
      post[, j, 2] <- with(
        subset(data_summ, treatment == 0),
        rgamma(N_mcmc, prior[1] + tot_events[j], prior[2] + tot_time[j]))
    }
  }

  return(post)

}
