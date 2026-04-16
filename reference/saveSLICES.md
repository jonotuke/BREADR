# saveSLICES

Plots all pairwise diagnostic plots (in a tibble as output by
callRelatedness), as produced by plotSLICE, to a folder. Options include
the width and height of the output files, and the units in which these
dimensions are measured.

## Usage

``` r
saveSLICES(
  in_tibble,
  outFolder = NULL,
  width = 297,
  height = 210,
  units = "mm",
  verbose = TRUE
)
```

## Arguments

- in_tibble:

  a tibble that is the output of the callRelatedness() function.

- outFolder:

  the folder into which all diagnostic plots will be saved

- width:

  the width of the output PDFs.

- height:

  the height of the output PDFs.

- units:

  the units for the height and width of the output PDFs.

- verbose:

  Controls the printing of progress to console.

## Value

nothing

## Examples

``` r
# \donttest{
saveSLICES(relatedness_example[1:3, ], outFolder=tempdir())
#> Starting pairwise plot creation.
#>   |                                                                              |                                                                      |   0%
#> Warning: Removed 4428 rows containing missing values or values outside the scale range
#> (`geom_ribbon()`).
#> Warning: Removed 15212 rows containing missing values or values outside the scale range
#> (`geom_ribbon()`).
#>   |                                                                              |===================================                                   |  50%
#> Warning: Removed 15184 rows containing missing values or values outside the scale range
#> (`geom_ribbon()`).
#>   |                                                                              |======================================================================| 100%
#> Completed
#> All plots in: /tmp/RtmpjORzCt/
#> [1] TRUE
# }
```
