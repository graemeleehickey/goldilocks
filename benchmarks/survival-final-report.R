# Summarize saved Bayesian survival final-analysis benchmarks without rerunning.
# Author: Graeme L. Hickey
# Rscript benchmarks/survival-final-report.R benchmarks/survival-final

cli <- commandArgs(trailingOnly = TRUE)
prefix <- if (length(cli)) cli[1] else "benchmarks/survival-final"
saved <- readRDS(paste0(prefix, "-results.rds"))
stopifnot(length(saved$results) == length(saved$metadata$scenarios))
summary <- saved$summary
comparison <- do.call(
  rbind,
  lapply(names(saved$results), function(name) {
    rows <- summary[summary$scenario == name, ]
    legacy <- rows[rows$method == "legacy", ]
    direct <- rows[rows$method == "direct_1000", ]
    higher <- rows[rows$method == "direct_50000", ]
    matched <- rows[rows$method == "direct_matched", ]
    data.frame(
      scenario = name,
      probability = legacy$reference_probability,
      observed_dropout = legacy$observed_dropouts / legacy$n,
      legacy_ms = legacy$median_ms,
      direct_1000_ms = direct$median_ms,
      direct_50000_ms = higher$median_ms,
      default_speedup = legacy$median_ms / direct$median_ms,
      higher_speedup = legacy$median_ms / higher$median_ms,
      legacy_mcse = legacy$empirical_probability_sd,
      direct_1000_mcse = direct$theoretical_probability_sd,
      direct_50000_mcse = higher$theoretical_probability_sd,
      default_mcse_ratio = direct$theoretical_probability_sd /
        legacy$empirical_probability_sd,
      equivalent_draws = matched$n_mcmc,
      matched_speedup = legacy$median_ms / matched$median_ms
    )
  })
)
print(comparison, digits = 5, row.names = FALSE)

decisions <- do.call(
  rbind,
  lapply(names(saved$results), function(name) {
    result <- saved$results[[name]]
    probability <- result$reference$probability
    do.call(
      rbind,
      lapply(saved$metadata$thresholds, function(threshold) {
        true_success <- probability > threshold
        data.frame(
          scenario = name,
          threshold = threshold,
          probability = probability,
          legacy_error = mean(
            (result$legacy_estimates$probability > threshold) != true_success
          ),
          direct_1000_error = mean(
            (result$direct_estimates$probability > threshold) != true_success
          ),
          direct_1000_error_exact = if (true_success) {
            pbinom(floor(threshold * 1000), 1000, probability)
          } else {
            pbinom(
              floor(threshold * 1000),
              1000,
              probability,
              lower.tail = FALSE
            )
          },
          direct_50000_error_exact = if (true_success) {
            pbinom(floor(threshold * 50000), 50000, probability)
          } else {
            pbinom(
              floor(threshold * 50000),
              50000,
              probability,
              lower.tail = FALSE
            )
          }
        )
      })
    )
  })
)
cat("\nNumerical classification error for the fixed datasets:\n")
print(decisions, digits = 5, row.names = FALSE)

observed <- summary[summary$empirical_repetitions > 0, ]
observed$probability_bias_z <- observed$probability_bias /
  sqrt(
    observed$empirical_probability_sd^2 /
      observed$empirical_repetitions +
      observed$reference_mcse^2
  )
observed$effect_bias_z <- observed$effect_bias /
  (observed$empirical_effect_sd / sqrt(observed$empirical_repetitions))
cat("\nConsistency with the reference (bias / its standard error):\n")
print(
  observed[c("scenario", "method", "probability_bias_z", "effect_bias_z")],
  digits = 4,
  row.names = FALSE
)
cat("\nDirect empirical / theoretical probability SD:\n")
direct <- observed[observed$method == "direct_1000", ]
print(
  setNames(
    direct$empirical_probability_sd / direct$theoretical_probability_sd,
    direct$scenario
  ),
  digits = 4
)

# Conditional decision errors follow an exact binomial distribution. A
# probability on the threshold has no nonzero separation to resolve, so the
# plotted curves omit the point at zero and join each side separately.
grid <- expand.grid(
  threshold = saved$metadata$thresholds,
  draws = c(1000L, 10000L, 50000L, 500000L),
  gap = c(
    seq(-0.015, -0.00001, length.out = 500L),
    seq(0.00001, 0.015, length.out = 500L)
  )
)
grid$probability <- grid$threshold + grid$gap
grid$success_probability <- pbinom(
  floor(grid$threshold * grid$draws),
  grid$draws,
  grid$probability,
  lower.tail = FALSE
)
grid$error <- ifelse(
  grid$gap > 0,
  1 - grid$success_probability,
  grid$success_probability
)
grid$draws <- factor(
  grid$draws,
  levels = c(1000, 10000, 50000, 500000),
  labels = c("1,000 (default)", "10,000", "50,000", "500,000")
)
grid$panel <- paste("Success threshold", grid$threshold)
grid$side <- ifelse(grid$gap > 0, "above", "below")
library(ggplot2)
figure <- ggplot(
  grid,
  aes(
    x = 100 * gap,
    y = 100 * error,
    colour = draws,
    group = interaction(draws, side)
  )
) +
  geom_hline(yintercept = 1, colour = "grey55", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey75") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~panel) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73", "#5B4B8A")) +
  scale_x_continuous(breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) +
  scale_y_continuous(limits = c(0, 55), breaks = seq(0, 50, 10)) +
  labs(
    title = "Monte Carlo error can change a near-threshold final decision",
    subtitle = "Exact binomial probabilities for independent posterior draws",
    x = "True posterior probability minus success threshold (percentage points)",
    y = "Probability of the wrong numerical decision (%)",
    colour = "Posterior draws",
    caption = paste(
      "Dashed line: 1% numerical decision error. This is not trial type I error.",
      "Author: Graeme L. Hickey",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )
ggsave(
  paste0(prefix, "-precision.png"),
  figure,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 180,
  bg = "white"
)
