# Simulate a complete clinical trial with event data drawn from a piecewise exponential distribution

Simulate a complete clinical trial with event data drawn from a
piecewise exponential distribution

## Usage

``` r
sim_comp_data(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = NULL,
  N_total,
  lambda = 0.3,
  lambda_time = NULL,
  end_of_study,
  block = 2,
  rand_ratio = c(1, 1),
  prop_loss = 0
)
```

## Arguments

- hazard_treatment:

  vector. Finite non-negative constant hazard rates under the treatment
  arm.

- hazard_control:

  vector. Finite non-negative constant hazard rates under the control
  arm.

- cutpoints:

  finite, positive, strictly increasing interior times at which the
  baseline hazard changes. The number of hazards for each arm must be
  one greater than the number of cutpoints. Default is `NULL`, which
  corresponds to a simple (non-piecewise) exponential model.

- N_total:

  integer. Maximum sample size allowable

- lambda:

  finite positive enrollment rates per unit time. Supply one rate for
  each interval defined by `lambda_time`. See
  [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  for the precise continuous-time process and time-origin convention.

- lambda_time:

  `NULL`, or finite, positive, strictly increasing internal times at
  which the enrollment rate changes. The initial boundary at zero is
  implicit, so `length(lambda)` must equal `length(lambda_time) + 1`.

- end_of_study:

  finite study endpoint, strictly greater than the last cutpoint.

- block:

  scalar. Block size for generating the randomization schedule.

- rand_ratio:

  vector. Randomization allocation for the ratio of control to
  treatment. Integer values mapping the size of the block. See
  [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  for more details.

- prop_loss:

  scalar. Overall proportion of subjects lost to follow-up. Subjects are
  selected at random for LTFU regardless of treatment assignment or
  event status. Each LTFU subject's observed time is drawn from a
  `Uniform(0, t)` distribution, where `t` is their potential event or
  censoring time. Since the LTFU time is always less than `t`, the event
  has not yet occurred at dropout and the subject is right-censored.
  Defaults to zero.

## Value

A data frame with 1 row per subject and columns:

- `time`: Time of event or censoring time.

- `treatment`: Treatment assignment, coded `1L` for the treatment arm
  and `0L` for the control arm. Single-arm designs have `treatment = 1L`
  for every subject.

- `event`: Indicator of whether event occurred (`1L` if occurred and
  `0L` if right-censored).

- `enrollment`: Time of patient enrollment relative to the time the
  trial enrolled the first patient. The package treats enrollment and
  randomization as occurring at the same time.

- `id`: Identification number for each patient.

- `loss_to_fu`: Indicator of whether the patient was lost to follow-up
  during observation.

## Details

Enrollment is simulated directly in continuous time by
[`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md).
The first patient is placed at time zero and all subsequent enrollment
times are measured from first patient in. No uniform jitter is added in
`sim_comp_data()`.

`lambda_time` and `cutpoints` both contain internal change times, but
they describe different clocks. `lambda_time` describes changes in the
trial's calendar-time enrollment rate measured from first patient in.
`cutpoints` describes changes in an individual subject's event hazard
measured from that subject's enrollment. They need not have the same
values or length. All time quantities supplied to a simulation should
nevertheless use one common unit, such as days or months.
