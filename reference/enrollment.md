# Simulate enrollment times

Simulate enrollment time using a piecewise Poisson distribution.

## Usage

``` r
enrollment(lambda = 1, N_total, lambda_time = 0)
```

## Arguments

- lambda:

  finite positive rate parameter(s) for the Poisson distribution.

- N_total:

  positive integer total sample size.

- lambda_time:

  finite, strictly increasing knots (of `length(lambda)`) indicating
  regions where a specific rate (`lambda`) applies. The first element
  must be `lambda_time = 0`, denoting the trial start time. The final
  element of `lambda` is assumed to be constant as `lambda_time` tends
  to infinity.

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

    enrollment(
      lambda = c(0.3, 0.7, 0.9, 1.2),
      N_total = 50,
      lambda_time = c(0, 5, 10, 15)
    )

## See also

This function is based on the `enrollment` function from the
[`bayesCT`](https://cran.r-project.org/package=bayesCT) R package.

## Examples

``` r
enrollment(lambda = c(0.003, 0.7), N_total = 100, lambda_time = c(0, 10))
#>   [1]   0   0   9   9   9  11  12  14  14  14  15  17  18  19  23  26  27  27
#>  [19]  30  31  32  33  35  36  37  37  37  38  38  38  43  45  47  50  52  54
#>  [37]  56  56  57  58  67  68  68  70  71  75  75  78  80  80  80  80  82  83
#>  [55]  83  84  86  87  90  91  95  96  97  98  99 100 101 101 103 103 104 104
#>  [73] 104 104 105 107 108 109 110 110 111 113 113 114 114 117 121 122 122 123
#>  [91] 123 124 125 126 129 129 131 131 131 131
enrollment(lambda = c(0.3, 0.5, 0.9, 1.2, 2.1), N_total = 200,
           lambda_time = c(0, 20, 30, 40, 60))
#>   [1]   0   4   5  15  16  16  17  18  20  21  21  21  22  25  26  29  29  29
#>  [19]  30  31  32  33  35  36  38  38  40  40  42  44  44  46  46  46  47  48
#>  [37]  49  49  49  50  50  52  52  52  53  53  54  55  55  55  55  55  55  56
#>  [55]  57  57  57  58  60  61  61  61  61  62  63  65  65  65  65  65  66  67
#>  [73]  67  67  68  69  69  69  69  70  71  71  71  71  71  71  72  72  73  73
#>  [91]  73  74  74  74  75  75  75  76  76  77  79  79  79  80  80  81  81  81
#> [109]  82  83  83  83  83  83  84  84  84  85  85  85  85  86  86  86  86  86
#> [127]  86  87  88  88  88  89  89  91  91  92  92  93  93  94  94  94  94  95
#> [145]  95  96  96  96  96  96  98  98  98  99 100 100 100 100 101 101 101 102
#> [163] 102 102 103 103 103 104 104 105 106 106 107 107 108 108 108 109 109 110
#> [181] 110 111 111 111 111 111 112 114 114 115 116 117 117 117 117 117 118 118
#> [199] 118 119
```
