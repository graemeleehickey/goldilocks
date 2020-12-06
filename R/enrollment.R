#' @title Simulate enrollment times
#'
#' @description Simulate enrollment time using a piecewise Poisson distribution.
#'
#' @param param vector. Rate parameter(s) for Poisson distribution.
#' @param time vector. Knots (of \code{length(param)} - 1) indicating end of
#'   time when a specific lambda is used. Note: final element of \code{param} is
#'   assumed to be constant as \code{time} tends to infinity.
#' @param N_total integer. Value of total sample size.
#'
#' @details Subject recruitment is assumed to follow a (piecewise stationary)
#'   Poisson process. We assume trial recruitment to be an independent process,
#'   thus the 'memoryless' property modeling of subject recruitment is used.
#'   Since the subject recruitment rate can vary over time, we can account for
#'   differential rates over time. Note that the first trial enrollment is
#'   assumed to occur at time zero.
#'
#'   To illustrate, suppose we use a piecewise function to specify the change in
#'   enrollment rate over time:
#'
#'   \deqn{
#'     \lambda = \left\{
#'       \begin{array}{ll}
#'         0.3 & \textrm{time} \in [0, 5) \\
#'         0.7 & \textrm{time} \in [5, 10) \\
#'         0.9 & \textrm{time} \in [10, 15) \\
#'         1.2 & \textrm{time} \in [15, \infty) \\
#'       \end{array}
#'     \right.
#' }
#'
#'   Then, to simulate individual patient enrollment dates with a sample size
#'   (\code{N_total}) of 50, we use
#'
#'   \code{enrollment(param = c(0.3, 0.7, 0.9, 1.2), N_total = 50, time = c(5, 10, 15))}
#'
#'   This function was ported from the
#'   [\code{bayesCT}](https://cran.r-project.org/web/packages/bayesCT/index.html)
#'   R package.
#'
#' @return A vector of enrollment times (from time of first patient enrollment)
#'   in days.
#'
#' @importFrom stats rpois
#' @export enrollment
#'
#' @examples
#' enrollment(param = c(0.003, 0.7), 100, time = 10)
#' enrollment(param = c(0.3, 0.5, 0.9, 1.2, 2.1), 200, c(20, 30, 40, 60))
enrollment <- function(param, N_total, time = NULL) {

  if (any(param <= 0)) {
    stop("The lambda(s) for poisson enrollment rate should be non-negative")
  }

  if (N_total <= 0) {
    stop("The sample size for enrollment needs to be greater than 0")
  }

  if (length(param) > 1  & (is.null(time) | (length(param) != (length(time) + 1)))) {
    stop("The cutoff time for lambda is not correct")
  }

  if (length(param) == 1 & !is.null(time)) {
    warning("The time input is not being used")
  }

  output <- NULL
  count <- 0
  # For constant lambda in Poisson distribution
  if (length(param) == 1) {
    while (length(output) < N_total) {
      count <- count + 1
      output <- c(output, rep(count, rpois(1, param)))
    }
  } else {
    # For different lambda values in Poisson distribution as a function of time
    while (length(output) < N_total) {
      count <- count + 1
      index <- min(c(which(time >= (count - 1)), length(param)))
      output <- c(output, rep(count, rpois(1, param[index])))
    }
  }

  # Adjusting trial and outputs
  output <- output[1:N_total]
  output <- output - output[1]
  return(output)

}
