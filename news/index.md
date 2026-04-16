# Changelog

## BREADR 1.1.0

- Added
  [`plotDOUGH()`](https://jonotuke.github.io/BREADR/reference/plotDOUGH.md)
  function
- Added
  [`priorSensitivity()`](https://jonotuke.github.io/BREADR/reference/priorSensitivity.md)
  function
- Changed R requirement for package to make it more accessible
- Removed commented code from callRelatedness.
- Made prior sum to one a warning and less restrictive for
  [`callRelatedness()`](https://jonotuke.github.io/BREADR/reference/callRelatedness.md).
- Made sure column names were in input data for functions, not exact.
- Added “by population” option to
  [`processEigenstrat()`](https://jonotuke.github.io/BREADR/reference/processEigenstrat.md)
  to make it faster.
- Added a new plot called `plotROLL()` to compare PMR values between
  groups.

## BREADR 1.0.3

CRAN release: 2025-04-14

## BREADR 1.0.4

## BREADR 1.0.3

CRAN release: 2025-04-14

- Reviewer for JOSS noticed that needed a later version than in
  description due to a call to the Matrix package. So added later R
  version requirement.

## BREADR 1.0.2

CRAN release: 2024-09-09

- Problem with processEigenstrat as it used a bash cut function native
  to linux and mac, but not available on Windows. The old version is
  still available as processEigenstrat_old.
- Change callRelatedness to use any function
- Problem with plotLoaf not showing legend is level is missing - fixed.

## BREADR 1.0.1

CRAN release: 2023-04-12

Patch as error in read files causing miscalculation of PMR.

## BREADR 1.0.0

CRAN release: 2023-04-07

This is the version that we released to CRAN.

## BREADR 0.0.0.9000

### Major changes

- This is the first development version. It has the primary functions.
- Also set up a pkgdown website to support.
