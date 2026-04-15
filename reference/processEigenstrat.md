# process Eigenstrat data - alternative version

A function that takes paths to an eigenstrat trio (ind, snp and geno
file) and returns the pairwise mismatch rate for all pairs on a thinned
set of SNPs. Options include choosing thinning parameter, subsetting by
population names, and filtering out SNPs for which deamination is
possible.

## Usage

``` r
processEigenstrat(
  indfile,
  genofile,
  snpfile,
  filter_length = NULL,
  pop_pattern = NULL,
  filter_deam = FALSE,
  outfile = NULL,
  chromosomes = NULL,
  verbose = TRUE,
  byPop = FALSE,
  minMerge = 2
)
```

## Arguments

- indfile:

  path to eigenstrat ind file

- genofile:

  path to eigenstrat geno file.

- snpfile:

  path to eigenstrat snp file.

- filter_length:

  the minimum distance between sites to be compared (to reduce the
  effect of LD).

- pop_pattern:

  a character vector of population names to filter the ind file if only
  some populations are to compared.

- filter_deam:

  a TRUE/FALSE for if C-\>T and G-\>A sites should be ignored.

- outfile:

  (OPTIONAL) a path and filename to which we can save the output of the
  function as a TSV, if NULL, no back up saved. If no outfile, then a
  tibble is returned.

- chromosomes:

  the chromosome to filter the data on.

- verbose:

  controls printing of messages to console

- byPop:

  boolean

- minMerge:

  TBA

## Value

out_tibble: A tibble containing four columns:

## Examples

``` r
# Use internal files to the package as an example
indfile <- system.file("extdata", "example.ind.txt", package="BREADR")
genofile <- system.file("extdata", "example.geno.txt", package="BREADR")
snpfile <- system.file("extdata", "example.snp.txt", package="BREADR")
processEigenstrat(
indfile, genofile, snpfile,
filter_length=1e5,
pop_pattern=NULL,
filter_deam=FALSE
)
#> Reading in SNP data.
#> Analysing chromosomes:
#> 1
#> Starting to read in genotype data.
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%  Complete.
#> 
#> Starting to compare genotypes and calculate PMR.
#>   |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
#> Complete.
#> 
#> # A tibble: 15 × 4
#>    pair        nsnps mismatch    pmr
#>    <chr>       <dbl>    <dbl>  <dbl>
#>  1 Ind1 - Ind2     2        1  0.5  
#>  2 Ind1 - Ind3    14        2  0.143
#>  3 Ind1 - Ind4    14        3  0.214
#>  4 Ind1 - Ind5     1        0  0    
#>  5 Ind1 - Ind6     4        0  0    
#>  6 Ind2 - Ind3    13        3  0.231
#>  7 Ind2 - Ind4    11        2  0.182
#>  8 Ind2 - Ind5     0        0 NA    
#>  9 Ind2 - Ind6     4        2  0.5  
#> 10 Ind3 - Ind4    35        5  0.143
#> 11 Ind3 - Ind5    10        0  0    
#> 12 Ind3 - Ind6    19        2  0.105
#> 13 Ind4 - Ind5    10        1  0.1  
#> 14 Ind4 - Ind6    20        7  0.35 
#> 15 Ind5 - Ind6     2        1  0.5  
```
