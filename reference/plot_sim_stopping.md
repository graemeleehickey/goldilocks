# Plot stopping outcomes from trial simulations

Draws a bar chart of expected-success, futility, and maximum-sample-size
outcomes together with a histogram of enrolled sample sizes. The input
can be the sims element returned by sim_trials or the complete
sim_trials result.

## Usage

``` r
plot_sim_stopping(x)
```

## Arguments

- x:

  A simulation result data frame or the list returned by sim_trials.

## Value

The simulation result data frame, invisibly.
