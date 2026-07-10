test_that("haz_to_prop preserves dimensions for one piecewise posterior draw", {
  post <- array(
    c(0.02, 0.03, 0.03, 0.04),
    dim = c(1, 2, 2)
  )

  out <- haz_to_prop(
    post = post,
    cutpoints = c(0, 12),
    end_of_study = 24,
    single_arm = FALSE
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1)
  expect_true(all(is.finite(unlist(out))))
})
