# split line

takes a line for a SNP file and splits into parts.

## Usage

``` r
split_line(x)
```

## Arguments

- x:

  line from SNP file

## Value

tibble with 6 columns.

## Examples

``` r
split_line("1_14.570829090394763     1        0.000000              14 A X")
#> # A tibble: 1 × 6
#>   snp                    chr   pos  site anc   der  
#>   <chr>                <dbl> <dbl> <dbl> <chr> <chr>
#> 1 1_14.570829090394763     1     0    14 A     X    
split_line("rs3094315  1  0.0  752566  G  A")
#> # A tibble: 1 × 6
#>   snp         chr   pos   site anc   der  
#>   <chr>     <dbl> <dbl>  <dbl> <chr> <chr>
#> 1 rs3094315     1     0 752566 G     A    
```
