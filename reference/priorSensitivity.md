# priorSensitivity

priorSensitivity

## Usage

``` r
priorSensitivity(
  in_BREADR,
  row,
  degrees = NULL,
  grid_space = 0.05,
  maxPrior = 0.5
)
```

## Arguments

- in_BREADR:

  A tibble containing the output of callRelatedness().

- row:

  Either the row number or pair name for which the sensitivity analysis
  should be run.

- degrees:

  A vector of integers identifying which degrees of relatedness degrees
  to plot. Same/Twins (1), First-Degree (2), Second-Degree (3) and
  Unrelated (4).

- grid_space:

  The space between prior probability values that are tested (smaller
  values mean a finer grid size).

- maxPrior:

  Maximum value a prior probability can take (the right-hand limit of
  the x-axis).

## Value

plot
