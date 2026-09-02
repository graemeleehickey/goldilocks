#' Capture the evaluated arguments of a function call
#'
#' @description Collects every formal argument from a function's evaluation
#'   environment. The resulting named list can be saved and passed back to the
#'   function with [do.call()].
#'
#' @param fn A function whose evaluated arguments should be collected.
#' @param env An environment containing the evaluated arguments of the current
#'   function call.
#'
#' @return A named list with one element for every formal argument.
#'
#' @keywords internal
#' @noRd
capture_arguments <- function(fn, env) {
  argument_names <- names(formals(fn))
  mget(argument_names, envir = env, inherits = FALSE)
}
