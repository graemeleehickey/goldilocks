# Summarize operating characteristics on the calendar-time scale

Summarizes trial duration, follow-up burden, and (when retained)
interim-look timing without adding any new simulation arguments.
Calendar time is measured from the first patient's enrollment at time
zero. `analysis_ready_time` is the time at which the last enrolled
subject's observed event or censoring becomes available; it does not
include external data-cleaning or database-lock delays.

Pass a complete result from
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
to retain requested, analyzed, and failed simulation counts. Interim
timing requires simulations run with `return_trace = TRUE`. A simulation
data frame can also be supplied, but only the trial-duration table can
then be calculated.

## Usage

``` r
summarise_calendar_time(data)
```

## Arguments

- data:

  A complete result returned by
  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md),
  a simulation `data.frame`, or a list of either form. Named list
  elements identify scenarios. Existing `scenario` columns and grouping
  variables are preserved.

## Value

An object of class `goldilocks_calendar_summary`, containing two wide
data frames:

- `trial_duration`: one row per scenario and stopping reason, plus an
  overall row. It reports simulation denominators, stopping counts,
  sample size, accrual-stop time, analysis-ready time, planned
  completion time, total person-time under follow-up, and peak
  concurrent follow-up.

- `interim_timing`: one row per scenario and interim look. It reports
  how often the look was reached, its calendar time, and the number of
  subjects actively under follow-up. This table has zero rows when
  traces are not available.

Continuous quantities are reported as means, Monte Carlo standard
errors, and the 10th, 50th, and 90th percentiles. Percentages use the
requested number of simulations as their denominator, so excluded
failures remain visible.
