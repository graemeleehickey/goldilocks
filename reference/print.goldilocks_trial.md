# Print a Goldilocks adaptive trial result

Prints the final trial summary and reports how many interim looks were
completed for a `goldilocks_trial` object returned by
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md).

## Usage

``` r
# S3 method for class 'goldilocks_trial'
print(x, ...)
```

## Arguments

- x:

  A `goldilocks_trial` result returned by
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  when `return_trace = TRUE`.

- ...:

  Additional arguments passed to
  [`base::print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html)
  when printing the final trial summary.

## Value

The input object, invisibly.
