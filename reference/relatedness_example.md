# relatedness_example

this is an example of the tibble made by callRelatedness()

## Usage

``` r
relatedness_example
```

## Format

### `relatedness_example`

A data frame with 15 rows and 13 columns:

- row:

  The row number

- pair:

  the pair of individuals that are compared.

- relationship:

  the highest posterior probability estimate of the degree of
  relatedness.

- pmr:

  the pairwise mismatch rate (mismatch/nsnps).

- sd:

  the estimated standard deviation of the pmr.

- mismatch:

  the number of sites which did not match for each pair.

- nsnps:

  the number of overlapping snps that were compared for each pair.

- ave_re:

  the value for the background relatedness used for normalisation.

- Same_Twins:

  the posterior probability associated with a same individual/twins
  classification.

- First_Degree:

  the posterior probability associated with a first-degree
  classification.

- Second_Degree:

  the posterior probability associated with a second-degree
  classification.

- Unrelated:

  the posterior probability associated with an unrelated classification.

- BF:

  A strength of confidence in the Bayes Factor associated with the
  highest posterior probability classification compared to the 2nd
  highest.
