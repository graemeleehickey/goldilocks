sim_oc_data <- function() {
  data.frame(
    scenario = c("target", "null", "small"),
    effect = c(0.7, 1, 0.85),
    power = c(0.9, 0.025, 0.55),
    stop_success = c(0.75, 0.01, 0.35),
    stop_futility = c(0.05, 0.7, 0.25),
    stop_max_N = c(0.2, 0.29, 0.4),
    mean_N = c(210, 180, 250),
    sd_N = c(55, 45, 60),
    stop_and_fail = c(0.02, 0.001, 0.01)
  )
}

test_that("operating-characteristic plot orders scenarios by effect", {
  operating_characteristics <- sim_oc_data()
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    matplot = function(x, y, ...) {
      captured$probability_x <- x
      captured$probabilities <- y
    },
    legend = function(...) invisible(NULL),
    .package = "graphics"
  )
  local_mocked_bindings(
    plot = function(x, y, ...) {
      captured$sample_size_x <- x
      captured$mean_N <- y
    },
    .package = "base"
  )

  returned <- plot_sim_ocs(
    operating_characteristics,
    effect = "effect",
    xlab = "True hazard ratio"
  )

  expect_identical(returned, operating_characteristics)
  expect_equal(captured$probability_x, c(0.7, 0.85, 1))
  expect_equal(
    unname(captured$probabilities[, "power"]),
    c(0.9, 0.55, 0.025)
  )
  expect_equal(captured$sample_size_x, c(0.7, 0.85, 1))
  expect_equal(captured$mean_N, c(210, 250, 180))
})

test_that("operating-characteristic plot accepts an effect vector", {
  operating_characteristics <- sim_oc_data()
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(
    plot_sim_ocs(
      operating_characteristics,
      effect = operating_characteristics$effect
    )
  )
})

test_that("operating-characteristic plot validates its inputs", {
  operating_characteristics <- sim_oc_data()

  expect_error(
    plot_sim_ocs(operating_characteristics[, -4], "effect"),
    "simulation summary"
  )
  expect_error(
    plot_sim_ocs(operating_characteristics, "unknown"),
    "not found"
  )
  expect_error(
    plot_sim_ocs(operating_characteristics, c(0.7, 0.85)),
    "one finite numeric value"
  )
  operating_characteristics$power[1] <- 1.1
  expect_error(
    plot_sim_ocs(operating_characteristics, "effect"),
    "between 0 and 1"
  )
})
