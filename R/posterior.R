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
#' @importFrom data.table data.table
#' @import survival
#' @export
posterior <- function(data, cutpoint, prior, N_mcmc, single_arm) {

  cutpoint <- cutpoint[-1] # Note: survSplit() doesn't like cuts at 0
  if (length(cutpoint) == 0) {
    cutpoint <- max(data$time)
  }

  data_survsplit <- survSplit(
    Surv(time, event) ~ .,
    data = data,
    cut = cutpoint,
    episode = "interval")

  time <- tstart <- event <- treatment <- interval <- NULL
  data_survsplit <- data.table(data_survsplit)
  data_summ <- data_survsplit[, .(n = length(time),
                                  tot_time = sum(time - tstart),
                                  tot_events = sum(event)),
                              by = list(treatment, interval)]

  # # dplyr implementation
  # data_summ <- data_survsplit %>%
  #   group_by(.data$treatment, .data$interval) %>%
  #   summarise(n = length(.data$time),
  #             tot_time = sum(.data$time - .data$tstart),
  #             tot_events = sum(.data$event)) %>%
  #   ungroup()

  nbreaks <- max(data_survsplit$interval) - 1

  post <- array(dim = c(N_mcmc, nbreaks + 1, 2))
  post_treatment <- matrix(nrow = N_mcmc, ncol = nbreaks + 1)
  post_control <- matrix(nrow = N_mcmc, ncol = nbreaks + 1)

  for (j in 1:(nbreaks + 1)) {
    post[, j, 1] <- with(
      subset(data_summ, treatment == 1),
      rgamma(N_mcmc, prior[1] + tot_events[j], prior[2] + tot_time[j])
    )
  }

  if (!single_arm) { # If control patients present
    for (j in 1:(nbreaks + 1)) {
      post[, j, 2] <- with(
        subset(data_summ, treatment == 0),
        rgamma(N_mcmc, prior[1] + tot_events[j], prior[2] + tot_time[j]))
    }
  }

  return(post)

}
