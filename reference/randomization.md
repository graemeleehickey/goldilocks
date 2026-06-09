# Randomization allocation

Implements a randomization allocation for control and treatment arms
with different randomization ratios and block sizes.

## Usage

``` r
randomization(N_total, block = 2, allocation = c(1, 1))
```

## Arguments

- N_total:

  integer. Total sample size for randomization allocation.

- block:

  vector. Block size for randomization. Note that it needs to be a
  multiple of the sum of `allocation`.

- allocation:

  vector. The randomization allocation in the order
  `c(control, treatment)`.

## Value

The randomization allocation with 0, 1 for control and treatment,
respectively.

## Details

Complete randomization may not always be ideal due to the chance of
drawing a large block of a single treatment arm, potentially impacting
the time to enrollment completion. Therefore, a block randomization
allocation may be preferable. The block randomization allocation
specification allows for different randomization ratios, but they must
be given in integer form. Additionally, the block size should be an
integer that is divisible by the sum of the randomization allocation;
see the examples.

## Examples

``` r
# Implementing treatment allocation for control to treatment with 1:1.5
# randomization ratio
randomization(N_total = 100, block = 5, allocation = c(2, 3))
#>   [1] 1 0 1 1 0 0 1 1 1 0 1 1 0 1 0 1 0 0 1 1 1 0 1 0 1 1 1 1 0 0 0 1 0 1 1 1 0
#>  [38] 1 1 0 0 1 0 1 1 0 1 1 0 1 0 1 0 1 1 1 0 1 0 1 0 1 0 1 1 1 1 0 0 1 1 1 0 0
#>  [75] 1 1 1 1 0 0 1 1 0 0 1 0 1 0 1 1 1 1 1 0 0 1 1 0 1 0

# Treatment allocation with 2:1 for control to treatment
randomization(N_total = 70, block = 9, allocation = c(2, 1))
#>  [1] 0 0 1 1 0 0 0 1 0 0 1 1 0 0 0 0 0 1 1 1 0 0 0 0 0 0 1 1 0 1 0 0 0 0 1 0 0 1
#> [39] 1 0 0 0 0 1 0 1 0 0 0 0 1 1 0 0 0 1 1 0 1 0 0 0 0 1 0 0 1 1 0 0

# Treatment allocation for control to treatment with 1:2 for control
# to treatment with multiple block sizes c(3, 9, 6)
randomization(N_total = 100, block = c(3, 9, 6), allocation = c(1, 2))
#>   [1] 1 0 1 0 1 0 0 1 1 1 1 1 0 1 1 1 1 0 1 1 0 1 0 1 1 0 1 1 0 1 1 0 1 1 1 0 1
#>  [38] 1 0 0 1 1 0 0 1 1 1 1 1 1 0 1 1 0 1 1 0 1 1 0 1 1 1 0 0 1 1 0 1 1 1 0 0 1
#>  [75] 1 1 0 1 0 1 0 1 1 1 1 0 1 1 0 1 0 1 1 0 1 1 1 1 0 1

# For complete randomization set the N_total to block size
randomization(N_total = 100, block = 100, allocation = c(1, 1))
#>   [1] 0 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 1 1 0 1 0 1 1 1 0 1 0 1 1 0 0 0 0 0 0 1 1
#>  [38] 0 0 0 0 0 0 0 0 1 1 1 1 0 1 1 1 0 1 0 0 1 1 0 0 0 0 1 0 1 1 0 1 1 1 0 0 0
#>  [75] 0 1 0 0 1 1 1 0 1 1 1 0 1 1 1 0 1 0 1 1 1 0 1 0 0 0
```
