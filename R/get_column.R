get_column_old <- function(genofile, col = 1){
  system(paste0('cut -c ',col,' ',genofile),intern=T)
}
get_column_new <- function(genofile, col = 1){
  read_lines(genofile) |>
    map(\(x) stringr::str_sub(x, col, col)) |>
    unlist()
}
pacman::p_load(conflicted, tidyverse, targets, data.table)
genofile <- system.file("extdata", "example.geno.txt", package = "BREADR")
identical(get_column_old(genofile, 1), get_column_new(genofile, 1))

microbenchmark::microbenchmark(
  get_column_old(genofile, 1),
  get_column_new(genofile, 1)
)
