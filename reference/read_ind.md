# read_ind

read_ind

## Usage

``` r
read_ind(filename)
```

## Arguments

- filename:

  a IND text file.

## Value

tibble with column headings: ind (CHR), sex (CHR), pop (CHR)

## Examples

``` r
ind_snpfile <- system.file("extdata", "example.ind.txt", package = "BREADR")
read_ind(ind_snpfile)
#> # A tibble: 6 × 3
#>   ind   sex   pop  
#>   <chr> <chr> <chr>
#> 1 Ind1  U     Pop1 
#> 2 Ind2  U     Pop1 
#> 3 Ind3  U     Pop1 
#> 4 Ind4  U     Pop1 
#> 5 Ind5  U     Pop2 
#> 6 Ind6  U     Pop2 
```
