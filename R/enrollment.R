#' @title Simulate enrollment times
#'
#' @description Simulate enrollment time using a piecewise Poisson distribution.
#'
#' @param lambda vector. Rate parameter(s) for Poisson distribution.
#' @param lambda_time vector. Knots (of \code{length(lambda)}) indicating
#'   regions where a specific hazard rate (\code{lambda}) applies. The first
#'   element is always \code{lambda_time = 0}, denoting the trial start time.
#'   Note: final element of \code{lambda} is assumed to be constant as
#'   \code{lambda_time} tends to infinity.
#' @param N_total integer. Value of total sample size.
#'
#' @details Subject recruitment is assumed to follow a (piecewise stationary)
#'   Poisson process. We assume trial recruitment to be an independent process,
#'   thus the 'memoryless' property modelling of subject recruitment is used.
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
#'   \code{enrollment(lambda = c(0.3, 0.7, 0.9, 1.2), N_total = 50,
#'         lambda_time = c(0, 5, 10, 15))}
#'
#' @return A vector of enrollment times (from time of first patient enrollment)
#'   in unit time (e.g. days).
#'
#' @seealso This function is based on the \code{enrollment} function from the
#'   [\code{bayesCT}](https://cran.r-project.org/package=bayesCT)
#'    R package.
#'
#' @importFrom stats rpois
#' @export
#'
#' @examples
#' enrollment(lambda = c(0.003, 0.7), N_total = 100, lambda_time = c(0, 10))
#' enrollment(lambda = c(0.3, 0.5, 0.9, 1.2, 2.1), N_total = 200,
#'            lambda_time = c(0, 20, 30, 40, 60))
enrollment <- function(lambda = 1, N_total, lambda_time = 0) {

  if (any(lambda <= 0)) {
    stop("The lambda(s) for Poisson enrollment rate should be non-negative")
  }

  if (N_total <= 0) {
    stop("The sample size for enrollment needs to be greater than 0")
  }

  if (length(lambda) != length(lambda_time)) {
    stop("The length of rates should match the length of knots")
  }

  if (is.null(lambda_time) | any(is.na(lambda_time))) {
    stop("The lambda_time argument is required")
  }

  if (lambda_time[1] != 0) {
    stop("The first cutpoint should always 0")
  }

  output <- NULL
  count <- 0
  # For constant lambda in Poisson distribution
  if (length(lambda) == 1) {
    while (length(output) < N_total) {
      count <- count + 1
      output <- c(output, rep(count, rpois(1, lambda)))
    }
  } else {
    # For different lambda values in Poisson distribution as a function of
    # lambda_time
    while (length(output) < N_total) {
      count <- count + 1
      index <- min(c(which(lambda_time[-1] >= (count - 1)), length(lambda)))
      output <- c(output, rep(count, rpois(1, lambda[index])))
    }
  }

  # Adjusting trial and outputs
  output <- output[1:N_total]
  output <- output - output[1]
  return(output)

}
