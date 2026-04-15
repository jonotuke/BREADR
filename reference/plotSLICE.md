# plotSLICE

A function for plotting the diagnostic information when classifying a
specific pair (defined by the row number or pair name) of individuals.
Output includes the PDFs for each degree of relatedness (given the
number of overlapping SNPs) in panel A, and the normalised posterior
probabilities for each possible degree of relatedness.

## Usage

``` r
plotSLICE(
  in_tibble,
  row,
  title = NULL,
  class_prior = rep(1/4, 4),
  showPlot = TRUE,
  which_plot = 0,
  labels = NULL
)
```

## Arguments

- in_tibble:

  a tibble that is the output of the callRelatedness() function.

- row:

  either the row number or pair name for which the posterior
  distribution is to be plotted.

- title:

  an optional title for the plot. If NULL, the pair from the
  user-defined row is used.

- class_prior:

  the prior probabilities for same/twin, 1st-degree, 2nd-degree,
  unrelated, respectively.

- showPlot:

  If TRUE, display plot. If FALSE, just pass plot as a variable.

- which_plot:

  if 1, returns just the plot of the posterior distributions, if 2
  returns just the normalised posterior values. Anything else returns
  both plots.

- labels:

  a length two character vector of labels for plots. Default is no
  labels.

## Value

a two-panel diagnostic ggplot object

## Examples

``` r
plotSLICE(relatedness_example, row=1)
#> Warning: Removed 4428 rows containing missing values or values outside the scale range
#> (`geom_ribbon()`).
```
