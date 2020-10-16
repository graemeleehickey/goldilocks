#' @title Calculate posterior probability of the alternative hypothesis
#'
#' @description For a given null hypothesis value (\code{h0}), the posterior
#'   probability of the alternative hypothesis is calculated. This can be
#'   compared to threshold values for declaring expected success or futility.
#'
#' @inheritParams survival_adapt
#' @param post_probs matrix (3 columns). Sampled draws from the posterior
#'   distribution of the event rate for the \code{treatment} arm and the
#'   \code{control} arm. If the trial is single-armed, the \code{control} arm
#'   probabilities will all be \code{NA}s. The third column is the derived
#'   treatment \code{effect}.
#'
#' @return A posterior probability (scalar \eqn{[0, 1]}).
#'
#' @export
test <- function(post_probs, alternative, h0) {

  p_treatment <- post_probs$p_treatment
  p_control <- post_probs$p_control
  control <- all(!is.na(p_control)) # Is there a control arm?

  if (control) {
    # Two-arm
    effect <- p_treatment - p_control
  } else {
    # One-arm
    effect <- p_treatment
  }

  # Apply test with direction
  if (alternative == "two-sided") {
    success <- max(c(mean(effect > h0), mean(effect < h0)))
  } else if (alternative == "greater") {
    success <- mean(effect > h0)
  } else {
    success <- mean(effect < h0)
  }

  return(success)

}
