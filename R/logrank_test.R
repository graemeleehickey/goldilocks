#' @title Calculate the log-rank test very quickly
#'
#' @param groupa vector of group a's survival times
#' @param groupb vector of group b's survival times
#' @param groupacensored vector of censored information of group a's survival
#'   times
#' @param groupbcensored vector of censored information of group b's survival
#'   times
#' @param onlyz (optional) calculate only z-statistic
#'
#' @return chi-squared statistic, z-statistic, p-value
#'
#' @examples
#' T1 <- c(6, 6, 6, 6, 7, 9, 10, 10, 11, 13, 16, 17, 19, 20, 22, 23, 25, 32, 32, 34, 35)
#' E1 <- c(1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0)
#' T2 <- c(1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22, 23)
#' E2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' logrank_test(T1, T2, E1, E2)
#' #1.679294e+01 -4.097919e+00, 4.168809e-05
#'
#' @import Rcpp
#' @noRd
#' @keywords internal
#' @useDynLib goldilocks, .registration = TRUE
logrank_test <- function(
    groupa,
    groupb,
    groupacensored,
    groupbcensored,
    onlyz = FALSE) {

    logrank_instance(groupa, groupb, groupacensored, groupbcensored, onlyz)

}
