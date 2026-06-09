# Simulate a complete clinical trial with event data drawn from a piecewise exponential distribution

Simulate a complete clinical trial with event data drawn from a
piecewise exponential distribution

## Usage

``` r
sim_comp_data(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = 0,
  N_total,
  lambda = 0.3,
  lambda_time = 0,
  end_of_study,
  block = 2,
  rand_ratio = c(1, 1),
  prop_loss = 0
)
```

## Arguments

- hazard_treatment:

  vector. Constant hazard rates under the treatment arm.

- hazard_control:

  vector. Constant hazard rates under the control arm.

- cutpoints:

  vector. Times at which the baseline hazard changes. Default is
  `cutpoints = 0`, which corresponds to a simple (non-piecewise)
  exponential model.

- N_total:

  integer. Maximum sample size allowable

- lambda:

  vector. Enrollment rates across simulated enrollment times. See
  [`enrollment`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  for more details.

- lambda_time:

  vector. Enrollment time(s) at which the enrollment rates change. Must
  be same length as lambda. See
  [`enrollment`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  for more details.

- end_of_study:

  scalar. Length of the study; i.e. time at which endpoint will be
  evaluated.

- block:

  scalar. Block size for generating the randomization schedule.

- rand_ratio:

  vector. Randomization allocation for the ratio of control to
  treatment. Integer values mapping the size of the block. See
  [`randomization`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  for more details.

- prop_loss:

  scalar. Overall proportion of subjects lost to follow-up. Subjects are
  selected at random for LTFU regardless of treatment arm or event
  status. Each LTFU subject's observed time is drawn from a
  `Uniform(0, t)` distribution, where `t` is their potential event or
  censoring time. Since the LTFU time is always less than `t`, the event
  has not yet occurred at dropout and the subject is right-censored.
  Defaults to zero.

## Value

A data frame with 1 row per subject and columns:

- `time:`:

  numeric. Time of event or censoring time.

- `treatment:`:

  integer. Treatment arm with values `1L` for experimental arm, and `0L`
  for control arm (only if `hazard_control` is given).

- `event:`:

  integer. Indicator of whether event occurred (`=1L` if occurred and
  `=0L` if right-censored).

- `enrollment:`:

  numeric. Time of patient enrollment relative to time trial enrolled
  first patient.

- `id:`:

  integer. Identification number for each patient.

- `loss_to_fu:`:

  logical. Indicator of whether the patient was lost to follow-up during
  the course of observation.
