# plotDOUGH

plotDOUGH

## Usage

``` r
plotDOUGH(
  in_tibble,
  indfile,
  nsnps_cutoff = NULL,
  fntsize = 7,
  verbose = TRUE,
  removeBetween = FALSE
)
```

## Arguments

- in_tibble:

  A tibble containing the output from processEigenstrat().

- indfile:

  Path to an ind file with columns for the individuals IDs, genetic sex
  (not used) and population.

- nsnps_cutoff:

  Minimum number of overlapping sites for a pair to be included.

- fntsize:

  Font size for x axis tick labels.

- verbose:

  Controls the printing of progress to console.

- removeBetween:

  TRUE/FALSE for whether to include the “between population” PMR values
  as their own group (FALSE), or to exclude them (TRUE).

## Value

plots
