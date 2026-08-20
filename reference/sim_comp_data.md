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
  rand_ratio = c(control = 1, treatment = 1),
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
  corresponds to a simple (non-piecewise) exponential model. Realized
  event times are assigned to analysis intervals using the survival
  counting-process convention, open on the left and closed on the right;
  an event exactly at a cutpoint therefore belongs to the interval
  ending at that cutpoint.

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

  length-two positive integer randomization allocation. Name the values
  `control` and `treatment`; either supplied order is accepted and
  normalized internally. A legacy unnamed vector remains accepted in
  `c(control, treatment)` order. Unequal unnamed values produce a
  warning because names may be required in a future major release. See
  [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  for more details.

- prop_loss:

  one or two probabilities. A scalar applies the same LTFU proportion to
  every arm. For a two-arm design, differential attrition can be
  specified with a length-two vector named `control` and `treatment`;
  the supplied order does not matter. Within each arm,
  `ceiling(prop_loss * arm size)` subjects are selected at random
  regardless of event status. Each selected subject's observed time is
  drawn from a `Uniform(0, t)` distribution, where `t` is their
  potential event or censoring time. Since the LTFU time is always less
  than `t`, the event has not yet occurred at dropout and the subject is
  right-censored. Single-arm designs require one probability. Defaults
  to zero.

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

PWEALL represents the continuous generating hazard with pieces closed on
the left and open on the right. This differs from the package's
open-left, closed-right convention for assigning realized times only at
the cutpoints themselves, which have probability zero under the
continuous model. The cumulative hazard, event-time distribution, and
generated simulations are therefore unchanged.
