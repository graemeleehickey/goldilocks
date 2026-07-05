#' @title Randomization allocation
#'
#' @description Generates a randomized treatment assignment sequence for
#'   control and treatment arms with different randomization ratios and block
#'   sizes.
#'
#' @param N_total integer. Total sample size for randomization allocation.
#' @param block vector. Block size for randomization. Note that it needs to be a
#'   multiple of the sum of `allocation`.
#' @param allocation vector. The randomization allocation in the order
#'   `c(control, treatment)`.
#'
#' @details Complete randomization may not always be ideal due to the chance of
#'   drawing a large block assigned to one treatment group, potentially
#'   impacting the time to enrollment completion. Therefore, a block
#'   randomization allocation may be preferable. The block randomization
#'   allocation specification allows for different randomization ratios, but
#'   they must be given in integer form. Additionally, the block size should be
#'   an integer that is divisible by the sum of the randomization allocation;
#'   see the examples.
#'
#' @return An integer treatment assignment vector, coded `0` for control and
#'   `1` for treatment.
#'
#' @export
#'
#' @examples
#' # Implementing treatment allocation for control to treatment with 1:1.5
#' # randomization ratio
#' randomization(N_total = 100, block = 5, allocation = c(2, 3))
#'
#' # Treatment allocation with 2:1 for control to treatment
#' randomization(N_total = 70, block = 9, allocation = c(2, 1))
#'
#' # Treatment allocation for control to treatment with 1:2 for control
#' # to treatment with multiple block sizes c(3, 9, 6)
#' randomization(N_total = 100, block = c(3, 9, 6), allocation = c(1, 2))
#'
#' # For complete randomization set the N_total to block size
#' randomization(N_total = 100, block = 100, allocation = c(1, 1))
randomization <- function(N_total, block = 2, allocation = c(1, 1)) {
  if (any(block %% 1 != 0) | any(block <= 0)) {
    stop("'block' must be a non-negative integer")
  }

  if (any(block %% sum(allocation) != 0)) {
    stop("Sum of 'allocation' must be a multiple of 'block'")
  }

  if (N_total < sum(block)) {
    stop("Number of subjects must be at least the size of 'block'")
  }

  if (any(allocation %% 1 != 0)) {
    stop("All values of 'allocation' must be integer values")
  }

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
