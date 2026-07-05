test_that("randomization returns correct length", {
  out <- randomization(N_total = 100, block = 2, allocation = c(1, 1))
  expect_length(out, 100)
})

test_that("randomization produces only 0s and 1s", {
  set.seed(3847)
  out <- randomization(N_total = 100, block = 4, allocation = c(1, 1))
  expect_true(all(out %in% c(0, 1)))
})

test_that("1:1 randomization gives balanced allocation within blocks", {
  set.seed(5612)
  out <- randomization(N_total = 100, block = 2, allocation = c(1, 1))
  # With block size 2 and 1:1, 50 complete blocks of 2 = exactly balanced
  expect_equal(sum(out), 50)
})

test_that("randomization respects unequal allocation ratio", {
  set.seed(9214)
  out <- randomization(N_total = 90, block = 3, allocation = c(1, 2))
  # 30 complete blocks: 30 control, 60 treatment
  expect_equal(sum(out == 0), 30)
  expect_equal(sum(out == 1), 60)
})

test_that("randomization works with multiple block sizes", {
  set.seed(1053)
  out <- randomization(N_total = 100, block = c(3, 9, 6), allocation = c(1, 2))
  expect_length(out, 100)
  expect_true(all(out %in% c(0, 1)))
})

test_that("randomization errors on non-integer block", {
  expect_error(
    randomization(N_total = 100, block = 2.5, allocation = c(1, 1)),
    "non-negative integer"
  )
})

test_that("randomization errors when block not divisible by allocation sum", {
  expect_error(
    randomization(N_total = 100, block = 3, allocation = c(1, 1)),
    "multiple"
  )
})

test_that("randomization errors when N_total < block", {
  expect_error(
    randomization(N_total = 1, block = 4, allocation = c(1, 1)),
    "at least the size"
  )
})

test_that("randomization works with multiple block sizes needing remainder fill", {
  set.seed(7238)
  # sum(block) = 18, N_total = 25 forces the loop to exhaust all block elements
  # and then fill a remainder using next_block
  out <- randomization(N_total = 25, block = c(6, 6, 6), allocation = c(1, 2))
  expect_length(out, 25)
  expect_true(all(out %in% c(0, 1)))
})

test_that("randomization preserves allocation in final partial block setup", {
  set.seed(2894)
  out <- randomization(N_total = 25, block = c(6, 9, 3), allocation = c(1, 2))
  complete_blocks <- split(out[1:24], rep(seq_len(4), c(6, 9, 3, 6)))

  expect_equal(
    unname(vapply(complete_blocks, length, integer(1))),
    c(6, 9, 3, 6)
  )
  expect_equal(unname(vapply(complete_blocks, sum, integer(1))), c(4, 6, 2, 4))
  expect_true(out[25] %in% c(0, 1))
})

test_that("randomization errors on non-integer allocation", {
  expect_error(
    randomization(N_total = 100, block = 2, allocation = c(0.5, 0.5)),
    "integer"
  )
})

test_that("randomization validates allocation length and positivity", {
  expect_error(
    randomization(N_total = 100, block = 2, allocation = c(1, 1, 1)),
    "two positive"
  )

  expect_error(
    randomization(N_total = 100, block = 2, allocation = c(0, 1)),
    "two positive"
  )

  expect_error(
    randomization(N_total = 100, block = 2, allocation = c(-1, 2)),
    "two positive"
  )

  expect_error(
    randomization(N_total = 100, block = 2, allocation = c(1, Inf)),
    "integer"
  )
})
