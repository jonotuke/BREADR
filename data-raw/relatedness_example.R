## code to prepare `kin_example` dataset goes here
relatedness_example <- callRelatedness(
  counts_example,
  class_prior=rep(0.25,4),
  average_relatedness=NULL,
  median_co=5e2,
  filter_n=1
)
usethis::use_data(relatedness_example, overwrite = TRUE)
