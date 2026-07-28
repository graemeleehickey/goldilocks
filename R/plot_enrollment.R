#' @title Plot an enrollment projection
#'
#' @description Draws the expected cumulative enrollment curve for a
#'   Goldilocks trial design, together with optional random enrollment
#'   trajectories and projected interim and maximum-sample-size milestones.
#'
#' @param x `NULL`, or a result returned by [survival_adapt()] or
#'   [sim_trials()]. Results created by current versions of `goldilocks` retain
#'   the evaluated enrollment design needed by this function.
#' @param lambda finite positive enrollment rates per unit time. Required when
#'   `x = NULL`; otherwise it can override the rate stored in `x`.
#' @param N_total positive integer maximum sample size. Required when
#'   `x = NULL`; otherwise it can override the value stored in `x`.
#' @param lambda_time `NULL`, or the finite, positive, strictly increasing
#'   internal times at which the enrollment rate changes. See [enrollment()].
#' @param interim_look `NULL`, or the enrollment counts at interim looks.
#' @param end_of_study optional positive follow-up duration. When available and
#'   `annotate = TRUE`, it is reported beneath the plot.
#' @param n_sim non-negative integer number of random enrollment trajectories
#'   to draw.
#' @param seed `NULL`, or a single non-negative integer seed for the random
#'   trajectories. A supplied seed does not alter the caller's random-number
#'   state.
#' @param time_unit `NULL`, or a character label for the design's unit of time,
#'   such as `"months"` or `"days"`.
#' @param xlab,ylab,main axis and main-title labels. With `xlab = NULL`, a label
#'   is constructed from `time_unit`.
#' @param annotate logical. Should the follow-up and simulation notes be drawn
#'   beneath the plot?
#' @param projection_col,simulation_col,milestone_col colours for the expected
#'   projection, random trajectories, and milestone guides.
#'
#' @details The blue projection is
#' \eqn{1 + \Lambda(t)}, where \eqn{\Lambda(t)} is the cumulative intensity of
#' the piecewise-constant Poisson enrollment process. The first patient is
#' fixed at time zero, consistently with [enrollment()]. A milestone's
#' projected time solves \eqn{1 + \Lambda(t) = N}. With a constant enrollment
#' rate this is also the mean arrival time, \eqn{(N - 1) / \lambda}. With a
#' piecewise rate it is an expected-count projection rather than the mean of
#' the corresponding arrival-time distribution.
#'
#' If `x` supplies a stored design, explicitly supplied design arguments
#' override the corresponding stored values. This makes it possible, for
#' example, to compare a fitted design with a different enrollment rate.
#'
#' @return Invisibly, a list containing the evaluated `design`, the
#'   `projection` data frame, the `milestones` data frame, and the simulated
#'   enrollment-time vectors in `simulations`.
#'
#' @examples
#' plot_enrollment(
#'   lambda = 20,
#'   N_total = 600,
#'   interim_look = 400,
#'   end_of_study = 12,
#'   n_sim = 20,
#'   seed = 20260727,
#'   time_unit = "months"
#' )
#'
#' # Piecewise enrollment rates are supported.
#' plot_enrollment(
#'   lambda = c(8, 20),
#'   lambda_time = 6,
#'   N_total = 200,
#'   interim_look = c(100, 150),
#'   n_sim = 5,
#'   seed = 1,
#'   time_unit = "months"
#' )
#'
#' @export
plot_enrollment <- function(
  x = NULL,
  lambda = NULL,
  N_total = NULL,
  lambda_time = NULL,
  interim_look = NULL,
  end_of_study = NULL,
  n_sim = 20L,
  seed = NULL,
  time_unit = NULL,
  xlab = NULL,
  ylab = "Cumulative number of enrolled patients",
  main = NULL,
  annotate = TRUE,
  projection_col = "#276E9B",
  simulation_col = "#777777",
  milestone_col = "#C8682A"
) {
  design <- enrollment_design_from_result(x)

  if (!missing(lambda)) {
    design$lambda <- lambda
  }
  if (!missing(N_total)) {
    design$N_total <- N_total
  }
  if (!missing(lambda_time)) {
    design$lambda_time <- lambda_time
  }
  if (!missing(interim_look)) {
    design$interim_look <- interim_look
  }
  if (!missing(end_of_study)) {
    design$end_of_study <- end_of_study
  }

  if (is.null(design$lambda) || is.null(design$N_total)) {
    stop(
      "Supply an object with a stored enrollment design, or supply both ",
      "'lambda' and 'N_total'"
    )
  }
  if (!"lambda_time" %in% names(design)) {
    design$lambda_time <- NULL
  }
  if (!"interim_look" %in% names(design)) {
    design$interim_look <- NULL
  }
  if (!"end_of_study" %in% names(design)) {
    design$end_of_study <- NULL
  }

  validate_enrollment_schedule(
    design$lambda,
    design$lambda_time,
    design$N_total
  )
  if (!is.null(design$interim_look)) {
    validate_positive_integer_vector(design$interim_look, "interim_look")
    if (
      any(diff(design$interim_look) <= 0) ||
        any(design$interim_look >= design$N_total)
    ) {
      stop(
        "'interim_look' must be strictly increasing and less than 'N_total'"
      )
    }
  }
  if (!is.null(design$end_of_study)) {
    if (
      length(design$end_of_study) != 1L ||
        !is.numeric(design$end_of_study) ||
        is.na(design$end_of_study) ||
        !is.finite(design$end_of_study) ||
        design$end_of_study <= 0
    ) {
      stop("'end_of_study' must be NULL or a single finite positive value")
    }
  }
  if (
    length(n_sim) != 1L ||
      !is.numeric(n_sim) ||
      is.na(n_sim) ||
      !is.finite(n_sim) ||
      n_sim < 0 ||
      n_sim != floor(n_sim)
  ) {
    stop("'n_sim' must be a single non-negative integer")
  }
  if (
    !is.null(seed) &&
      (length(seed) != 1L ||
        !is.numeric(seed) ||
        is.na(seed) ||
        !is.finite(seed) ||
        seed < 0 ||
        seed > .Machine$integer.max ||
        seed != floor(seed))
  ) {
    stop(
      "'seed' must be NULL or a single non-negative integer no greater than ",
      ".Machine$integer.max"
    )
  }
  if (
    !is.null(time_unit) &&
      (!is.character(time_unit) ||
        length(time_unit) != 1L ||
        is.na(time_unit) ||
        !nzchar(time_unit))
  ) {
    stop("'time_unit' must be NULL or a non-empty character value")
  }
  validate_plot_label(xlab, "xlab", allow_null = TRUE)
  validate_plot_label(ylab, "ylab")
  validate_plot_label(main, "main", allow_null = TRUE)
  validate_logical_scalar(annotate, "annotate")
  validate_plot_colour(projection_col, "projection_col")
  validate_plot_colour(simulation_col, "simulation_col")
  validate_plot_colour(milestone_col, "milestone_col")

  design <- new_enrollment_design(
    lambda = design$lambda,
    N_total = design$N_total,
    lambda_time = design$lambda_time,
    interim_look = design$interim_look,
    end_of_study = design$end_of_study
  )

  if (!is.null(seed)) {
    old_kind <- RNGkind()
    old_seed_exists <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
    if (old_seed_exists) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit(
      {
        do.call(RNGkind, as.list(old_kind))
        if (old_seed_exists) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (
          exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        ) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
    set.seed(seed)
  }

  simulations <- replicate(
    n_sim,
    enrollment(
      lambda = design$lambda,
      N_total = design$N_total,
      lambda_time = design$lambda_time
    ),
    simplify = FALSE
  )

  milestone_n <- c(design$interim_look, design$N_total)
  milestone_time <- inverse_enrollment_intensity(
    milestone_n - 1,
    lambda = design$lambda,
    lambda_time = design$lambda_time
  )
  milestones <- data.frame(
    N = milestone_n,
    projected_time = milestone_time
  )

  projection_time <- if (max(milestone_time) > 0) {
    seq(0, max(milestone_time), length.out = 500L)
  } else {
    0
  }
  projection <- data.frame(
    time = projection_time,
    enrolled = 1 +
      cumulative_enrollment_intensity(
        projection_time,
        lambda = design$lambda,
        lambda_time = design$lambda_time
      )
  )

  simulation_max <- if (length(simulations) > 0L) {
    max(vapply(simulations, max, numeric(1)))
  } else {
    0
  }
  plot_max_time <- max(milestone_time, simulation_max)
  xlim <- if (plot_max_time > 0) {
    c(0, 1.04 * plot_max_time)
  } else {
    c(0, 1)
  }
  y_padding <- max(1, 0.13 * design$N_total)
  ylim <- c(0, design$N_total + y_padding)
  if (is.null(xlab)) {
    xlab <- if (is.null(time_unit)) {
      "Time since first patient in"
    } else {
      paste0("Enrollment time (", time_unit, ")")
    }
  }

  has_annotation <- annotate &&
    (!is.null(design$end_of_study) || n_sim > 0L)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mar = c(if (has_annotation) 6.1 else 4.4, 4.8, 3.0, 0.8),
    mgp = c(2.6, 0.8, 0)
  )

  graphics::plot(
    projection$time,
    projection$enrolled,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    axes = FALSE,
    bty = "l",
    cex.lab = 1.15
  )
  graphics::axis(1)
  y_ticks <- pretty(c(0, design$N_total))
  y_ticks <- y_ticks[y_ticks >= 0 & y_ticks <= design$N_total]
  graphics::axis(2, at = y_ticks, las = 1)
  grid_ticks <- y_ticks[y_ticks > 0]
  if (length(grid_ticks) > 0L) {
    graphics::abline(
      h = grid_ticks,
      col = "#E6E6E6",
      lwd = 0.8
    )
  }

  for (simulated_enrollment in simulations) {
    graphics::lines(
      simulated_enrollment,
      seq_len(design$N_total),
      type = "s",
      col = grDevices::adjustcolor(simulation_col, alpha.f = 0.22),
      lwd = 0.9
    )
  }
  graphics::segments(
    x0 = milestones$projected_time,
    y0 = 0,
    x1 = milestones$projected_time,
    y1 = design$N_total,
    col = milestone_col,
    lty = 2,
    lwd = 1.5
  )
  graphics::lines(
    projection$time,
    projection$enrolled,
    col = projection_col,
    lwd = 3
  )
  graphics::points(
    milestones$projected_time,
    milestones$N,
    pch = 19,
    col = projection_col,
    cex = 0.9
  )
  time_label <- if (length(design$lambda) == 1L) "Avg." else "Projected"
  graphics::text(
    x = milestones$projected_time,
    y = design$N_total + 0.45 * y_padding,
    labels = sprintf(
      "N = %d\n%s %.1f%s",
      milestones$N,
      time_label,
      milestones$projected_time,
      if (is.null(time_unit)) "" else paste0(" ", time_unit)
    ),
    col = grDevices::adjustcolor(milestone_col, alpha.f = 0.9),
    cex = 0.78,
    xpd = NA
  )
  graphics::box(bty = "l")

  if (has_annotation) {
    if (!is.null(design$end_of_study)) {
      unit_label <- if (is.null(time_unit)) "time units" else time_unit
      graphics::mtext(
        sprintf(
          "Final analysis after the last patient completes %s %s of follow-up",
          format(design$end_of_study, trim = TRUE),
          unit_label
        ),
        side = 1,
        line = 4.7,
        adj = 0,
        cex = 0.72,
        col = "#555555"
      )
    }
    if (n_sim > 0L) {
      graphics::mtext(
        sprintf(
          "%d simulated enrollment %s",
          n_sim,
          if (n_sim == 1L) "trajectory" else "trajectories"
        ),
        side = 1,
        line = 4.7,
        adj = 1,
        cex = 0.72,
        col = "#555555"
      )
    }
  }

  invisible(list(
    design = design,
    projection = projection,
    milestones = milestones,
    simulations = simulations
  ))
}

#' @keywords internal
#' @noRd
new_enrollment_design <- function(
  lambda,
  N_total,
  lambda_time,
  interim_look,
  end_of_study
) {
  list(
    lambda = lambda,
    N_total = N_total,
    lambda_time = lambda_time,
    interim_look = interim_look,
    end_of_study = end_of_study
  )
}

#' @keywords internal
#' @noRd
enrollment_design_from_result <- function(x) {
  if (is.null(x)) {
    return(list())
  }

  design <- attr(x, "enrollment_design", exact = TRUE)
  if (is.null(design)) {
    stop(
      "'x' does not contain a stored enrollment design; re-fit it with the ",
      "current version of goldilocks or supply the design arguments directly"
    )
  }
  if (!is.list(design)) {
    stop("the enrollment design stored in 'x' is invalid")
  }

  design
}

#' @keywords internal
#' @noRd
cumulative_enrollment_intensity <- function(time, lambda, lambda_time) {
  interval_starts <- c(0, lambda_time)
  cumulative_at_start <- c(
    0,
    cumsum(diff(interval_starts) * lambda[seq_along(lambda_time)])
  )
  interval <- findInterval(time, interval_starts)

  cumulative_at_start[interval] +
    (time - interval_starts[interval]) * lambda[interval]
}

#' @keywords internal
#' @noRd
inverse_enrollment_intensity <- function(intensity, lambda, lambda_time) {
  interval_starts <- c(0, lambda_time)
  cumulative_at_start <- c(
    0,
    cumsum(diff(interval_starts) * lambda[seq_along(lambda_time)])
  )
  interval <- findInterval(intensity, cumulative_at_start)

  interval_starts[interval] +
    (intensity - cumulative_at_start[interval]) / lambda[interval]
}

#' @keywords internal
#' @noRd
validate_plot_label <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(invisible(TRUE))
  }
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be a single character value")
  }

  invisible(TRUE)
}

#' @keywords internal
#' @noRd
validate_plot_colour <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("'", name, "' must be a single colour value")
  }
  tryCatch(
    grDevices::col2rgb(x),
    error = function(e) {
      stop("'", name, "' must be a valid colour", call. = FALSE)
    }
  )

  invisible(TRUE)
}
