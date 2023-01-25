test_that("processEigenstart errors", {
  # load files
  indfile <- system.file("extdata", "example.ind.txt", package = "BREAD")
  genofile <- system.file("extdata", "example.geno.txt", package = "BREAD")
  snpfile <- system.file("extdata", "example.snp.txt", package = "BREAD")
  # check input
  expect_error(
    processEigenstrat(indfile = "no-file.txt")
  )
  expect_error(
    processEigenstrat(indfile = indfile, snpfile = "no-file.txt")
  )
  expect_error(
    processEigenstrat(indfile = indfile,
                      snpfile = snpfile,
                      genofile = "no-file.txt")
  )
  expect_error(
    processEigenstrat(indfile = indfile,
                      snpfile = snpfile,
                      genofile = genofile,
                      filter_length = -1)
  )
  expect_error(
    processEigenstrat(indfile = indfile,
                      snpfile = snpfile,
                      genofile = genofile,
                      filter_deam = "NO")
  )
})
