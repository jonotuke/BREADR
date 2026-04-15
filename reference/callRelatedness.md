# callRelatedness

A function that takes PMR observations, and (given a prior distribution
for degrees of relatedness) returns the posterior probabilities of all
pairs of individuals being (a) the same individual/twins, (b)
first-degree related, (c) second-degree related or (d) "unrelated"
(third-degree or higher). The highest posterior probability degree of
relatedness is also returned as a hard classification. Options include
setting the background relatedness (or using the sample median), a
minimum number of overlapping SNPs if one uses the sample median for
background relatedness, and a minimum number of overlapping SNPs for
including pairs in the analysis.

## Usage

``` r
callRelatedness(
  pmr_tibble,
  class_prior = rep(0.25, 4),
  average_relatedness = NULL,
  median_co = 500,
  filter_n = 1
)
```

## Arguments

- pmr_tibble:

  a tibble that is the output of the processEigenstrat function.

- class_prior:

  the prior probabilities for same/twin, 1st-degree, 2nd-degree,
  unrelated, respectively.

- average_relatedness:

  a single numeric value, or a vector of numeric values, to use as the
  average background relatedness. If NULL, the sample median is used.

- median_co:

  if average_relatedness is left NULL, then the minimum cutoff for the
  number of overlapping snps to be included in the median calculation is
  500.

- filter_n:

  the minimum number of overlapping SNPs for which pairs are removed
  from the entire analysis. If NULL, default is 1.

## Value

results_tibble: A tibble containing 13 columns:

- row: The row number

- pair: the pair of individuals that are compared.

- relationship: the highest posterior probability estimate of the degree
  of relatedness.

- pmr: the pairwise mismatch rate (mismatch/nsnps).

- sd: the estimated standard deviation of the pmr.

- mismatch: the number of sites which did not match for each pair.

- nsnps: the number of overlapping snps that were compared for each
  pair.

- ave_re;: the value for the background relatedness used for
  normalisation.

- Same_Twins: the posterior probability associated with a same
  individual/twins classification.

- First_Degree: the posterior probability associated with a first-degree
  classification.

- Second_Degree: the posterior probability associated with a
  second-degree classification.

- Unrelated: the posterior probability associated with an unrelated
  classification.

- BF: A strength of confidence in the Bayes Factor associated with the
  highest posterior probability classification compared to the 2nd
  highest. (No longer included)

## Examples

``` r
callRelatedness(counts_example,
  class_prior=rep(0.25,4),
  average_relatedness=NULL,
  median_co=5e2,filter_n=1
)
#> # A tibble: 15 × 12
#>      row pair       relationship   pmr      sd mismatch nsnps ave_rel Same_Twins
#>    <int> <chr>      <fct>        <dbl>   <dbl>    <dbl> <dbl>   <dbl>      <dbl>
#>  1     1 Ind1 - In… Unrelated    0.204 0.0103       310  1518   0.218  6.71e- 26
#>  2     2 Ind1 - In… Unrelated    0.222 0.00428     2093  9435   0.218  1.22e-214
#>  3     3 Ind1 - In… Unrelated    0.224 0.00458     1854  8283   0.218  2.00e-194
#>  4     4 Ind1 - In… Unrelated    0.235 0.0116       314  1336   0.218  2.68e- 37
#>  5     5 Ind1 - In… Unrelated    0.215 0.00867      481  2242   0.218  9.82e- 46
#>  6     6 Ind2 - In… Unrelated    0.218 0.00432     1988  9119   0.218  5.06e-195
#>  7     7 Ind2 - In… Unrelated    0.213 0.00458     1699  7984   0.218  4.95e-156
#>  8     8 Ind2 - In… Unrelated    0.229 0.0122       270  1179   0.218  1.80e- 30
#>  9     9 Ind2 - In… Unrelated    0.215 0.00927      423  1965   0.218  1.10e- 40
#> 10    10 Ind3 - In… Same_Twins   0.108 0.00214     2253 20952   0.218  1   e+  0
#> 11    11 Ind3 - In… Unrelated    0.213 0.00458     1703  7994   0.218  6.83e-157
#> 12    12 Ind3 - In… Unrelated    0.218 0.00394     2398 10994   0.218  2.05e-235
#> 13    13 Ind4 - In… Unrelated    0.210 0.00489     1451  6924   0.218  1.92e-127
#> 14    14 Ind4 - In… Unrelated    0.220 0.00419     2141  9745   0.218  2.95e-214
#> 15    15 Ind5 - In… Unrelated    0.220 0.00994      383  1739   0.218  3.64e- 39
#> # ℹ 3 more variables: First_Degree <dbl>, Second_Degree <dbl>, Unrelated <dbl>
```
