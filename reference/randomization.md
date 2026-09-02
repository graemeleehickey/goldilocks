# Generate a block-randomized treatment sequence

Generates a randomized treatment assignment sequence for control and
treatment arms using a fixed allocation ratio and one or more permitted
block sizes.

## Usage

``` r
randomization(N_total, block = 2, allocation = c(control = 1, treatment = 1))
```

## Arguments

- N_total:

  A required positive integer giving the total number of treatment
  assignments.

- block:

  A positive integer vector containing one or more permitted block
  sizes. Every block size must be a multiple of `sum(allocation)`. The
  default is `2`.

- allocation:

  A length-two positive integer vector giving the control to treatment
  allocation ratio. The default is `c(control = 1, treatment = 1)`. Name
  the values `control` and `treatment`; either supplied order is
  accepted and matched by name. A legacy unnamed vector remains accepted
  in `c(control, treatment)` order. Unequal unnamed values produce a
  warning because names may be required in a future major release.

## Value

An integer treatment assignment vector, coded `0` for control and `1`
for treatment.

## Details

Complete randomization may not always be ideal due to the chance of
drawing a large block assigned to one treatment group, potentially
impacting the time to enrollment completion. Therefore, a block
randomization allocation may be preferable. The block randomization
allocation specification allows for different two-arm randomization
ratios, but they must be given in integer form. For every value `b` in
`block`, the required relationship is `b %% sum(allocation) == 0`; see
the equal- and unequal-allocation examples below.

## Examples

``` r
# Implementing treatment allocation for control to treatment with 1:1.5
# randomization ratio
randomization(
  N_total = 100,
  block = 5,
  allocation = c(control = 2, treatment = 3)
)
#>   [1] 0 1 1 1 0 1 1 1 0 0 0 1 1 0 1 1 0 0 1 1 1 0 0 1 1 0 1 1 0 1 0 1 1 1 0 0 1
#>  [38] 1 1 0 1 1 0 1 0 0 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0 1 1 1 0 1 1 1 0 0 1 1 0 0
#>  [75] 1 0 1 0 1 1 0 0 1 1 1 1 1 1 0 0 1 0 1 1 0 0 1 1 1 0

# Treatment allocation with 2:1 for control to treatment
randomization(
  N_total = 70,
  block = 9,
  allocation = c(treatment = 1, control = 2)
)
#>  [1] 0 1 1 0 1 0 0 0 0 0 0 0 1 0 0 1 0 1 0 1 0 1 0 0 0 1 0 1 1 0 0 1 0 0 0 0 0 0
#> [39] 0 1 1 0 0 0 1 1 0 0 1 0 0 0 1 0 0 0 1 0 1 0 0 1 0 1 0 0 1 0 0 0

# Treatment allocation for control to treatment with 1:2 for control
# to treatment with multiple block sizes c(3, 9, 6)
randomization(
  N_total = 100,
  block = c(3, 9, 6),
  allocation = c(control = 1, treatment = 2)
)
#>   [1] 1 0 1 1 1 0 1 1 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1 0 1 0 1 1 1 1 1 0 1 1 0 1 0
#>  [38] 1 1 0 1 1 1 0 1 1 1 0 0 0 1 1 1 1 1 0 1 0 1 1 1 1 0 0 1 1 0 0 1 1 1 1 1 1
#>  [75] 0 0 1 1 1 1 1 0 1 0 0 1 1 1 0 1 0 1 1 0 1 0 0 1 1 1

# For complete randomization set the N_total to block size
randomization(
  N_total = 100,
  block = 100,
  allocation = c(control = 1, treatment = 1)
)
#>   [1] 0 1 0 0 0 1 1 1 0 0 1 1 0 0 1 1 0 1 1 0 0 0 0 0 0 1 0 0 0 1 1 1 1 1 1 0 1
#>  [38] 1 1 0 1 1 0 0 1 1 0 0 1 0 0 1 0 1 0 0 0 0 0 1 1 0 1 1 1 0 0 1 0 1 1 0 1 0
#>  [75] 1 0 0 0 1 0 1 1 0 1 1 0 1 0 0 1 1 1 1 0 1 1 1 0 0 0

# randomization() is a two-arm helper; a multi-arm allocation is rejected.
try(randomization(
  N_total = 60,
  block = 6,
  allocation = c(1, 1, 1)
), silent = TRUE)
```
