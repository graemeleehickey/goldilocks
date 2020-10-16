#' @title Posterior distribution of piecewise exponential constant hazard rates
#'
#' @description Using the Beta-Gamma conjugacy property, the posterior
#'   distribution of the piecewise hazard rates (\eqn{\lambda_j}, for \code{j=1,
#'   \dots, J}) is calculated and sampled from.
#'
#' @inheritParams survival_adapt
#' @param data data frame. Minimum requirements are 3 columns: event time
#'   (\code{time}), indicator of the event (\code{event}), and indicator for
#'   treatment arm (\code{treatment}). Other columns can be included in the data
#'   frame and will be handled in the split.
#'
#' @return A list of two data frame elements: \code{post_treatment} and
#'   \code{post_control}. Each data frame is of length \code{N_mcmc} and has
#'   \eqn{J} columns -- one column per each constant hazard piece.
#'
#' @importFrom stats rgamma
#' @import survival
#' @export
posterior <- function(data, cutpoint, prior, N_mcmc) {

  cutpoint <- cutpoint[-1] # survSplit doesn't like cuts at 0
  if (length(cutpoint) == 0) {
    cutpoint <- max(data$time)
  }

  data_survsplit <- survSplit(
    Surv(time, event) ~ .,
    data = data,
    cut = cutpoint,
    episode = "interval")

  data_summ <- data_survsplit %>%
    group_by(treatment, interval) %>%
    summarise(n = n(),
              tot_time = sum(time - tstart),
              tot_events = sum(event))

  nbreaks <- max(data_survsplit$interval) - 1
  post_treatment <- matrix(nrow = N_mcmc, ncol = nbreaks + 1)
  post_control <- matrix(nrow = N_mcmc, ncol = nbreaks + 1)

  for (j in 1:(nbreaks + 1)) {
    post_treatment[, j] <- with(
      dplyr::filter(data_summ, treatment == 1),
      rgamma(N_mcmc, prior[1] + tot_events[j], prior[2] + tot_time[j])
    )
  }

  if (any(data$treatment == 0)) { # If control patients present
    for (j in 1:(nbreaks + 1)) {
      post_control[, j] <- with(
        dplyr::filter(data_summ, treatment == 0),
        rgamma(N_mcmc, prior[1] + tot_events[j], prior[2] + tot_time[j]))
    }
  }

  return(list("post_treatment" = post_treatment,
              "post_control" = post_control))

}
