# Plot stopping outcomes from trial simulations

Draws a stacked bar chart of stopping outcomes by enrolled sample size,
with colours distinguishing expected-success, futility, and
maximum-sample-size outcomes. The `type` argument controls whether the
function draws marginal, conditional, or cumulative bars, or a flowchart
through successive interim looks. Bar-chart subtitles state the
denominator used by the selected view. The input can be the `sims`
element returned by
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
or the complete
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
result.

## Usage

``` r
plot_sim_stopping(
  x,
  type = c("marginal", "conditional", "cumulative", "flowchart")
)
```

## Arguments

- x:

  A simulation result data frame or the list returned by sim_trials.

- type:

  Character string specifying the percentages to plot. `"marginal"`
  shows the percentage of all simulated trials ending at each sample
  size; its bars sum to 100 percent across sample sizes. `"conditional"`
  shows the percentage stopping at each look among trials still active
  at the start of that look. `"cumulative"` shows the status of all
  simulated trials after each look; every bar sums to 100 percent and
  includes trials continuing to the next look. `"flowchart"` starts with
  all simulated trials and branches at each look into futility,
  continued enrollment, and expected-success nodes labelled with trial
  counts.

## Value

For bar-chart types, the simulation result data frame, invisibly. For
`type = "flowchart"`, a `DiagrammeR` `grViz` htmlwidget.

## Details

The flowchart requires the `N_max` column and uses interim sample sizes
observed in `N_enrolled`. When the complete result from
`sim_trials(return_trace = TRUE)` is supplied, sample sizes recorded in
`traces` are also included, so looks at which no trial stopped still
appear. The flowchart is rendered with
[`DiagrammeR::grViz()`](https://rich-iannone.github.io/DiagrammeR/reference/grViz.html).
