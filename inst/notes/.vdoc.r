#
#
#
#
#
#
#
#
#
#
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
#
#
#
library(BREADR)
pacman::p_load(conflicted, tidyverse, targets, data.table)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
sim_geno(n_ind = 5, n_snp = 10, here::here("inst", "notes", "geno.txt"))
#
#
#
#
#
#
#| eval: false
geno_list[[i]] <- (system(paste0('cut -c ',ind$row_number[i],' ',genofile),intern=T) %>% 
  dplyr::na_if("9"))[snps$relative_pos]
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
x <- 12345
n <- 1
for(i in 1:5){
  print(x %/% 10^(trunc(log10(x)) - i + 1) %% 10)
}
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
x <- 00345
n <- 1
for(i in 1:5){
  print(x %/% 10^(trunc(log10(x)) - i + 1) %% 10)
}
#
#
#
#
#
x <- 00345
n <- 1
for(i in 1:5){
  print(x %/% 10^(5 - i) %% 10)
}
#
#
#
#
#
#
n_snp <- 500
n_ind <- 10
file <- here::here("inst", "notes", "small-geno.txt")
sim_geno(n_ind = n_ind, n_snp = n_snp, filename = file)
#
#
#
#
library(tictoc)
tic()
geno_list_1 <- vector(mode='list',length=n_ind)
for(i in 1:n_ind){
  geno_list_1[[i]] <- (system(paste0('cut -c ',i,' ',file),intern=T)) 
}
toc()
#
#
#
#
#
#
tic()
geno_list_2 <- list(length=n_ind)
dt <- fread(file)
for(i in 1:n_ind){
  dt[, digit := V1 %/% 10^(n_ind - i) %% 10,]
  geno_list_2[[i]] <- as.character(dt[[2]])
}
toc()
#
#
#
#
identical(
  geno_list_1[[1]],
  geno_list_2[[1]]
  )
#
#
#
#
#
#
n_snp <- 10000
n_ind <- 200
file <- here::here("inst", "notes", "large-geno.txt")
sim_geno(n_ind = n_ind, n_snp = n_snp, filename = file)
#
#
#
#
library(tictoc)
tic()
geno_list_1 <- vector(mode='list',length=n_ind)
for(i in 1:n_ind){
  geno_list_1[[i]] <- (system(paste0('cut -c ',i,' ',file),intern=T)) 
}
toc()
#
#
#
#
#
#
tic()
geno_list_2 <- list(length=n_ind)
dt <- fread(file)
for(i in 1:n_ind){
  dt[, digit := V1 %/% 10^(n_ind - i) %% 10,]
  geno_list_2[[i]] <- as.character(dt[[2]])
}
toc()
#
#
#
#
identical(
  geno_list_1[[1]],
  geno_list_2[[1]]
  )
#
#
#
