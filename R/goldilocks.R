#' @keywords internal
"_PACKAGE"
#' @name goldilocks
#' @title goldilocks
#'
#' @description The goal of `goldilocks` is to implement the Goldilocks Bayesian
#'   adaptive design proposed by Broglio et al. (2014) for time-to-event
#'   endpoint trials, both one- and two-arm, with an underlying piecewise
#'   exponential hazard model. The method can be used for a confirmatory trial
#'   to select a trial's sample size based on accumulating data. During accrual,
#'   frequent sample size selection analyses are made and predictive
#'   probabilities are used to determine whether the current sample size is
#'   sufficient or whether continuing accrual would be futile. The algorithm
#'   explicitly accounts for complete follow-up of all patients before the
#'   primary analysis is conducted. Broglio et al. (2014) refer to this as a
#'   Goldilocks trial design, as it is constantly asking the question, **“Is the
#'   sample size too big, too small, or just right?”**
#'
#' @references
#' Broglio KR, Connor JT, Berry SM. Not too big, not too small: a Goldilocks
#' approach to sample size selection. *Journal of Biopharmaceutical Statistics*,
#' 2014; **24(3)**: 685–705.
#'
#' @importFrom utils globalVariables
NULL


# Quiets concerns of R CMD check re: no visible binding
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("id", "subject_impute_futility",
                           "subject_impute_success", "treatment"))
}
