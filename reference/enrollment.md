# Simulate enrollment times

Simulate enrollment time using a piecewise Poisson distribution.

## Usage

``` r
enrollment(lambda = 1, N_total, lambda_time = 0)
```

## Arguments

- lambda:

  vector. Rate parameter(s) for Poisson distribution.

- N_total:

  integer. Value of total sample size.

- lambda_time:

  vector. Knots (of `length(lambda)`) indicating regions where a
  specific hazard rate (`lambda`) applies. The first element is always
  `lambda_time = 0`, denoting the trial start time. Note: final element
  of `lambda` is assumed to be constant as `lambda_time` tends to
  infinity.

## Value

A vector of enrollment times (from time of first patient enrollment) in
unit time (e.g. days).

## Details

Subject recruitment is assumed to follow a (piecewise stationary)
Poisson process. We assume trial recruitment to be an independent
process, thus the 'memoryless' property modelling of subject recruitment
is used. Since the subject recruitment rate can vary over time, we can
account for differential rates over time. Note that the first trial
enrollment is assumed to occur at time zero.

To illustrate, suppose we use a piecewise function to specify the change
in enrollment rate over time:

\$\$ \lambda = \left\\ \begin{array}{ll} 0.3 & \textrm{time} \in \[0, 5)
\\ 0.7 & \textrm{time} \in \[5, 10) \\ 0.9 & \textrm{time} \in \[10, 15)
\\ 1.2 & \textrm{time} \in \[15, \infty) \\ \end{array} \right. \$\$

Then, to simulate individual patient enrollment dates with a sample size
(`N_total`) of 50, we use

`enrollment(lambda = c(0.3, 0.7, 0.9, 1.2), N_total = 50, lambda_time = c(0, 5, 10, 15))`

## See also

This function is based on the `enrollment` function from the
[`bayesCT`](https://cran.r-project.org/package=bayesCT) R package.

## Examples

``` r
enrollment(lambda = c(0.003, 0.7), N_total = 100, lambda_time = c(0, 10))
#>   [1]   0   0   0   2   3   5   5   5   6   8   9  10  14  17  18  18  21  22
#>  [19]  23  24  26  27  28  28  28  29  29  29  34  36  38  41  43  45  47  47
#>  [37]  48  49  58  59  59  61  62  66  66  69  71  71  71  71  73  74  74  75
#>  [55]  77  78  81  82  86  87  88  89  90  91  92  92  94  94  95  95  95  95
#>  [73]  96  98  99 100 101 101 102 104 104 105 105 108 112 113 113 114 114 115
#>  [91] 116 117 120 120 122 122 122 122 123 125
enrollment(lambda = c(0.3, 0.5, 0.9, 1.2, 2.1), N_total = 200,
           lambda_time = c(0, 20, 30, 40, 60))
#>   [1]   0   1  11  12  12  13  14  16  17  17  18  21  22  25  25  27  28  29
#>  [19]  31  32  34  34  36  36  38  40  40  42  42  42  43  44  45  45  45  46
#>  [37]  46  48  48  48  49  49  50  51  51  51  51  51  51  52  53  53  53  57
#>  [55]  57  58  59  61  61  61  61  61  62  63  63  63  64  65  65  65  65  66
#>  [73]  67  67  67  67  67  67  68  68  69  69  69  70  70  70  71  71  71  72
#>  [91]  72  73  75  75  75  76  76  77  77  77  78  79  79  79  79  79  80  80
#> [109]  80  81  81  81  81  82  82  82  82  82  82  83  84  84  84  85  85  87
#> [127]  87  88  88  89  89  90  90  90  90  91  91  92  92  92  92  92  94  94
#> [145]  94  95  96  96  96  96  97  97  97  98  98  98  99  99  99 100 100 101
#> [163] 102 102 103 103 104 104 104 105 105 106 106 107 107 107 107 107 108 110
#> [181] 110 111 112 113 113 113 113 113 114 114 114 115 115 115 116 116 117 117
#> [199] 119 119
```
