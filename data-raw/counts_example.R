# load files ----
indfile <- "data-raw/example.ind.txt"
genofile <- "data-raw/example.geno.txt"
snpfile <- "data-raw/example.snp.txt"

# load package ----
pacman::p_load(tidyverse, targets, BREAD)

# use function processEigenstrat to get counts_example ----
counts_example <- processEigenstrat(
  indfile, genofile, snpfile,
  filter_length=1e5,
  pop_pattern=NULL,
  filter_deam=FALSE
)

# save data ----
usethis::use_data(counts_example, overwrite = TRUE)
