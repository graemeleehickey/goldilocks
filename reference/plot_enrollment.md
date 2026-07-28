# Plot an enrollment projection

Draws the expected cumulative enrollment curve for a Goldilocks trial
design, together with optional random enrollment trajectories and
projected interim and maximum-sample-size milestones.

## Usage

``` r
plot_enrollment(
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
)
```

## Arguments

- x:

  `NULL`, or a result returned by
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  or
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md).
  Results created by current versions of `goldilocks` retain the
  evaluated enrollment design needed by this function.

- lambda:

  finite positive enrollment rates per unit time. Required when
  `x = NULL`; otherwise it can override the rate stored in `x`.

- N_total:

  positive integer maximum sample size. Required when `x = NULL`;
  otherwise it can override the value stored in `x`.

- lambda_time:

  `NULL`, or the finite, positive, strictly increasing internal times at
  which the enrollment rate changes. See
  [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md).

- interim_look:

  `NULL`, or the enrollment counts at interim looks.

- end_of_study:

  optional positive follow-up duration. When available and
  `annotate = TRUE`, it is reported beneath the plot.

- n_sim:

  non-negative integer number of random enrollment trajectories to draw.

- seed:

  `NULL`, or a single non-negative integer seed for the random
  trajectories. A supplied seed does not alter the caller's
  random-number state.

- time_unit:

  `NULL`, or a character label for the design's unit of time, such as
  `"months"` or `"days"`.

- xlab, ylab, main:

  axis and main-title labels. With `xlab = NULL`, a label is constructed
  from `time_unit`.

- annotate:

  logical. Should the follow-up and simulation notes be drawn beneath
  the plot?

- projection_col, simulation_col, milestone_col:

  colours for the expected projection, random trajectories, and
  milestone guides.

## Value

Invisibly, a list containing the evaluated `design`, the `projection`
data frame, the `milestones` data frame, and the simulated
enrollment-time vectors in `simulations`.

## Details

The blue projection is \\1 + \Lambda(t)\\, where \\\Lambda(t)\\ is the
cumulative intensity of the piecewise-constant Poisson enrollment
process. The first patient is fixed at time zero, consistently with
[`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md).
A milestone's projected time solves \\1 + \Lambda(t) = N\\. With a
constant enrollment rate this is also the mean arrival time, \\(N - 1) /
\lambda\\. With a piecewise rate it is an expected-count projection
rather than the mean of the corresponding arrival-time distribution.

If `x` supplies a stored design, explicitly supplied design arguments
override the corresponding stored values. This makes it possible, for
example, to compare a fitted design with a different enrollment rate.

## Examples

``` r
plot_enrollment(
  lambda = 20,
  N_total = 600,
  interim_look = 400,
  end_of_study = 12,
  n_sim = 20,
  seed = 20260727,
  time_unit = "months"
)


# Piecewise enrollment rates are supported.
plot_enrollment(
  lambda = c(8, 20),
  lambda_time = 6,
  N_total = 200,
  interim_look = c(100, 150),
  n_sim = 5,
  seed = 1,
  time_unit = "months"
)

```
