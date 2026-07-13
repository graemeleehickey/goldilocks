# Plot stopping outcomes from trial simulations

Draws a stacked bar chart of final enrolled sample sizes, with colours
distinguishing expected-success, futility, and maximum-sample-size
outcomes. Each bar is labelled with its marginal percentage of simulated
trials. The input can be the sims element returned by sim_trials or the
complete sim_trials result.

## Usage

``` r
plot_sim_stopping(x)
```

## Arguments

- x:

  A simulation result data frame or the list returned by sim_trials.

## Value

The simulation result data frame, invisibly.
