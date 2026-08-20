#' @title Randomization allocation
#'
#' @description Generates a randomized treatment assignment sequence for
#'   control and treatment arms with different randomization ratios and block
#'   sizes.
#'
#' @param N_total integer. Total sample size for randomization allocation.
#' @param block vector. One or more positive integer block sizes. Every block
#'   size must be a multiple of `sum(allocation)`.
#' @param allocation length-two positive integer vector. Name the values
#'   `control` and `treatment`; either supplied order is accepted and normalized
#'   internally. A legacy unnamed vector remains accepted in
#'   `c(control, treatment)` order. Unequal unnamed values produce a warning
#'   because names may be required in a future major release.
#'
#' @details Complete randomization may not always be ideal due to the chance of
#'   drawing a large block assigned to one treatment group, potentially
#'   impacting the time to enrollment completion. Therefore, a block
#'   randomization allocation may be preferable. The block randomization
#'   allocation specification allows for different two-arm randomization
#'   ratios, but they must be given in integer form. For every value `b` in
#'   `block`, the required relationship is `b %% sum(allocation) == 0`; see the
#'   equal- and unequal-allocation examples below.
#'
#' @return An integer treatment assignment vector, coded `0` for control and
#'   `1` for treatment.
#'
#' @export
#'
#' @examples
#' # Implementing treatment allocation for control to treatment with 1:1.5
#' # randomization ratio
#' randomization(
#'   N_total = 100,
#'   block = 5,
#'   allocation = c(control = 2, treatment = 3)
#' )
#'
#' # Treatment allocation with 2:1 for control to treatment
#' randomization(
#'   N_total = 70,
#'   block = 9,
#'   allocation = c(treatment = 1, control = 2)
#' )
#'
#' # Treatment allocation for control to treatment with 1:2 for control
#' # to treatment with multiple block sizes c(3, 9, 6)
#' randomization(
#'   N_total = 100,
#'   block = c(3, 9, 6),
#'   allocation = c(control = 1, treatment = 2)
#' )
#'
#' # For complete randomization set the N_total to block size
#' randomization(
#'   N_total = 100,
#'   block = 100,
#'   allocation = c(control = 1, treatment = 1)
#' )
#'
#' # randomization() is a two-arm helper; a multi-arm allocation is rejected.
#' try(randomization(
#'   N_total = 60,
#'   block = 6,
#'   allocation = c(1, 1, 1)
#' ), silent = TRUE)
randomization <- function(
  N_total,
  block = 2,
  allocation = c(control = 1, treatment = 1)
) {
  allocation <- validate_randomization_args(N_total, block, allocation)

  next_block <- NULL

  # Creating different block sizes for multiple blocks
  blocking <- rep(block, N_total %/% sum(block))
  n_blocking <- length(blocking)
  extra_blocks <- vector(mode = typeof(block), length = length(block))
  n_extra_blocks <- 0
  blocking_total <- sum(blocking)

  for (k in 1:length(block)) {
    if ((blocking_total + block[k]) < N_total) {
      n_extra_blocks <- n_extra_blocks + 1
      extra_blocks[n_extra_blocks] <- block[k]
      blocking_total <- blocking_total + block[k]
      next_block <- block[(k %% length(block)) + 1]
    } else {
      break
    }
  }

  if (n_extra_blocks > 0) {
    blocking <- c(blocking, extra_blocks[seq_len(n_extra_blocks)])
    n_blocking <- length(blocking)
  }

  # Making sure the next block is assigned
  if (is.null(next_block)) {
    next_block <- block[1]
  }

  # Within each block, randomize with the correct allocation
  sampling <- integer(N_total)
  start <- 1
  for (m in seq_len(n_blocking)) {
    item <- rep(
      rep(0:1, times = allocation),
      each = blocking[m] / sum(allocation)
    )
    end <- start + blocking[m] - 1
    sampling[start:end] <- sample(item)
    start <- end + 1
  }

  # Fill up the remainder of the allocation using next block
  if (N_total >= start) {
    item <- rep(
      rep(0:1, times = allocation),
      each = next_block / sum(allocation)
    )
    sampling[start:N_total] <- sample(item, size = N_total - start + 1)
  }

  return(sampling)
}
