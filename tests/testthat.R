# This file is used by R CMD check and devtools::test()
library(testthat)

# Should work from any directory now!
library(testthat)
test_file(here::here("tests", "testthat", "test-physmap.R"))

