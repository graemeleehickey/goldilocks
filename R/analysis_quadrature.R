# Numerical checks for two-arm Beta quadrature. These diagnose probabilities;
# they do not replace the integration result or change the stopping rules.
check_beta_quadrature <- function(
  probability,
  treatment,
  control,
  h0,
  alternative
) {
  tolerance <- 1e-6
  symmetry <- h0 == 0 &&
    treatment$alpha == control$alpha &&
    treatment$beta == control$beta &&
    abs(probability - 0.5) > tolerance

  # Cantelli's inequality for D = p_treatment - p_control - h0.
  # With negative mean, Pr(D < 0) >= mean^2 / (variance + mean^2);
  # with positive mean, Pr(D < 0) <= variance / (variance + mean^2).
  mean_difference <- treatment$mean - control$mean - h0
  variance <- treatment$variance + control$variance
  denominator <- variance + mean_difference^2
  moments <- FALSE
  if (
    is.finite(mean_difference) &&
      is.finite(variance) &&
      variance >= 0 &&
      is.finite(denominator) &&
      denominator > 0
  ) {
    probability_less <- if (alternative == "less") {
      probability
    } else {
      1 - probability
    }
    if (mean_difference < 0) {
      moments <- probability_less < mean_difference^2 / denominator - tolerance
    } else if (mean_difference > 0) {
      moments <- probability_less > variance / denominator + tolerance
    }
  }
  if (symmetry || moments) {
    warn_quadrature_diagnostics(c(
      flagged_calculations = 1,
      symmetry = as.numeric(symmetry),
      moments = as.numeric(moments)
    ))
  }
  invisible(NULL)
}

new_quadrature_diagnostics <- function() {
  diagnostics <- new.env(parent = emptyenv())
  diagnostics$counts <- c(flagged_calculations = 0, symmetry = 0, moments = 0)
  diagnostics
}

warn_quadrature_diagnostics <- function(counts) {
  if (counts[["flagged_calculations"]] == 0) {
    return(invisible(NULL))
  }
  message <- paste0(
    "Numerical accuracy checks flagged ",
    counts[["flagged_calculations"]],
    " two-arm quadrature probability calculation(s) (symmetry: ",
    counts[["symmetry"]],
    "; posterior mean-variance bounds: ",
    counts[["moments"]],
    "). Success and stopping decisions may be unreliable. ",
    "Check these probabilities using sufficiently precise Monte Carlo Beta ",
    "draws (bin_method = \"mc\"). The quadrature results have not been changed."
  )
  warning(structure(
    list(message = message, call = NULL, counts = counts),
    class = c("goldilocks_quadrature_warning", "warning", "condition")
  ))
}

# Catch only this diagnostic warning; leave other warning policies unchanged.
with_quadrature_warnings <- function(expr, diagnostics) {
  withCallingHandlers(expr, goldilocks_quadrature_warning = function(w) {
    diagnostics$counts <- diagnostics$counts + w$counts
    invokeRestart("muffleWarning")
  })
}

finish_quadrature_warnings <- function(result, diagnostics, warn = TRUE) {
  if (diagnostics$counts[["flagged_calculations"]] > 0) {
    attr(result, "quadrature_diagnostics") <- c(
      diagnostics$counts,
      absolute_tolerance = 1e-6
    )
    if (warn) warn_quadrature_diagnostics(diagnostics$counts)
  }
  result
}

# Workers return counts to the calling process, which issues the run warning.
capture_quadrature_warnings <- function(expr) {
  diagnostics <- new_quadrature_diagnostics()
  value <- with_quadrature_warnings(expr, diagnostics)
  finish_quadrature_warnings(value, diagnostics, warn = FALSE)
}

combine_quadrature_warnings <- function(results) {
  diagnostics <- new_quadrature_diagnostics()
  for (result in results) {
    counts <- attr(result, "quadrature_diagnostics", exact = TRUE)
    if (!is.null(counts)) {
      diagnostics$counts <- diagnostics$counts +
        counts[names(diagnostics$counts)]
    }
  }
  diagnostics
}
