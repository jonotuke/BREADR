# test_degree

Test if a degree of relatedness is consistent with an observed PMR

## Usage

``` r
test_degree(in_tibble, row, degree, verbose = TRUE)
```

## Arguments

- in_tibble:

  a tibble that is the output of the callRelatedness() function.

- row:

  either the row number or pair name for which the posterior
  distribution is to be plotted.

- degree:

  the degree of relatedness to be tested.

- verbose:

  a logical (boolean) for whether all test output should be printed to
  screen.

## Value

the associated p-value for the test

## Examples

``` r
test_degree(relatedness_example, 1, 1)
#> Testing H0       : "Ind1 - Ind2" are 1st-degree relatives.
#> Expected PMR     : 0.1635
#> Observed PMR     : 0.2042
#> Estimated degree : 2.9826
#> p-value          : 2.666909e-05
#> Decision         : Reject H0
#> [1] 2.666909e-05
```
