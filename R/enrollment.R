#' @title Simulate enrollment times
#'
#' @description Simulate enrollment time using a piecewise Poisson distribution.
#'
#' @param lambda finite positive rate parameter(s) for the Poisson distribution.
#' @param lambda_time finite, strictly increasing knots (of `length(lambda)`)
#'   indicating regions where a specific rate (`lambda`) applies. The first
#'   element must be `lambda_time = 0`, denoting the trial start time. The final
#'   element of `lambda` is assumed to be constant as `lambda_time` tends to
#'   infinity.
#' @param N_total positive integer total sample size.
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
#'   (`N_total`) of 50, we use
#'
#'   ```
#'   enrollment(
#'     lambda = c(0.3, 0.7, 0.9, 1.2),
#'     N_total = 50,
#'     lambda_time = c(0, 5, 10, 15)
#'   )
#'   ```
#'
#' @return A vector of enrollment times (from time of first patient enrollment)
#'   in unit time (e.g. days).
#'
#' @seealso This function is based on the `enrollment` function from the
#'   [`bayesCT`](https://cran.r-project.org/package=bayesCT) R package.
#'
#' @importFrom stats rpois
#' @export
#'
#' @examples
#' enrollment(lambda = c(0.003, 0.7), N_total = 100, lambda_time = c(0, 10))
#' enrollment(lambda = c(0.3, 0.5, 0.9, 1.2, 2.1), N_total = 200,
#'            lambda_time = c(0, 20, 30, 40, 60))
enrollment <- function(lambda = 1, N_total, lambda_time = 0) {
  validate_enrollment_schedule(lambda, lambda_time, N_total)

  chunks <- list()
  n_enrolled <- 0
  n_chunks <- 0
  count <- 0
  # For constant lambda in Poisson distribution
  if (length(lambda) == 1) {
    while (n_enrolled < N_total) {
      count <- count + 1
      n_new <- rpois(1, lambda)
      if (n_new > 0) {
        n_chunks <- n_chunks + 1
        chunks[[n_chunks]] <- rep(count, n_new)
        n_enrolled <- n_enrolled + n_new
      }
    }
  } else {
    # For different lambda values in Poisson distribution as a function of
    # lambda_time
    while (n_enrolled < N_total) {
      count <- count + 1
      index <- min(c(which(lambda_time[-1] > (count - 1)), length(lambda)))
      n_new <- rpois(1, lambda[index])
      if (n_new > 0) {
        n_chunks <- n_chunks + 1
        chunks[[n_chunks]] <- rep(count, n_new)
        n_enrolled <- n_enrolled + n_new
      }
    }
  }

  # Adjusting trial and outputs
  output <- unlist(chunks, use.names = FALSE)
  output <- output[1:N_total]
  output <- output - output[1]
  return(output)
}
