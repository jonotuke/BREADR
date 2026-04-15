# read_snp

read_snp

## Usage

``` r
read_snp(filename)
```

## Arguments

- filename:

  a SNP text file.

## Value

tibble with column headings: snp (CHR), chr (DBL), pos (DBL), site
(DBL), anc (CHR), and der (CHR).

## Examples

``` r
std_snpfile <- system.file("extdata", "example.snp.txt", package = "BREADR")
broken_snpfile <- system.file("extdata", "broken.snp.txt", package = "BREADR")
read_snp(std_snpfile)
#> # A tibble: 1,000 × 6
#>    snp         chr     pos    site anc   der  
#>    <chr>       <chr> <dbl>   <int> <chr> <chr>
#>  1 rs3094315   1         0  752566 G     A    
#>  2 rs7419119   1         0  842013 T     G    
#>  3 rs13302957  1         0  891021 G     A    
#>  4 rs6696609   1         0  903426 C     T    
#>  5 rs8997      1         0  949654 A     G    
#>  6 rs9442372   1         0 1018704 A     G    
#>  7 rs147606383 1         0 1045331 G     A    
#>  8 rs4970405   1         0 1048955 A     G    
#>  9 rs11807848  1         0 1061166 T     C    
#> 10 rs4970421   1         0 1108637 G     A    
#> # ℹ 990 more rows
read_snp(broken_snpfile)
#> # A tibble: 100 × 6
#>    snp                  chr        pos  site anc   der  
#>    <chr>                <chr>    <dbl> <int> <chr> <chr>
#>  1 1_14.570829090394763 1     0           14 A     X    
#>  2 1_31.762980536196316 1     0           31 X     X    
#>  3 1_58.40133360272569  1     0.000001    58 A     X    
#>  4 1_146.13704767420893 1     0.000001   146 X     X    
#>  5 1_350.25559508884965 1     0.000003   350 A     X    
#>  6 1_373.54266589096653 1     0.000004   373 X     X    
#>  7 1_376.1963053157634  1     0.000004   376 X     X    
#>  8 1_402.00713564772997 1     0.000004   402 A     X    
#>  9 1_473.15441796651334 1     0.000005   473 A     X    
#> 10 1_608.1074747673024  1     0.000006   608 X     X    
#> # ℹ 90 more rows
```
