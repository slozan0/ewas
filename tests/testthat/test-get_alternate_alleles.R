# Tests for get_alternate_alleles function
source("../../scripts/functions/easy_chi2_fun.r")

test_that("get_alternate_alleles marks reference as non-alternate", {
  # Reference is A, have A and C
  observed <- c(100, 50, 0, 0, 0, 0)
  result <- get_alternate_alleles("A", observed)

  # A should be FALSE (it's the reference)
  expect_false(result[1])
  # C should be TRUE (it's alternate and present)
  expect_true(result[2])
  # Rest should be FALSE (not present)
  expect_false(result[3])
  expect_false(result[4])
  expect_false(result[5])
  expect_false(result[6])
})

test_that("get_alternate_alleles handles all nucleotides", {
  # Test each reference nucleotide
  observed <- c(10, 10, 10, 10, 10, 10)

  for (ref in c("A", "C", "G", "T", "I", "D")) {
    result <- get_alternate_alleles(ref, observed)

    # Should have exactly 5 TRUE values (all except reference)
    expect_equal(sum(result), 5)
  }
})

test_that("get_alternate_alleles handles zero observations", {
  # Reference is A, only A observed
  observed <- c(100, 0, 0, 0, 0, 0)
  result <- get_alternate_alleles("A", observed)

  # All should be FALSE (no alternates present)
  expect_equal(sum(result), 0)
})

test_that("get_alternate_alleles rejects invalid reference", {
  observed <- c(10, 10, 0, 0, 0, 0)

  expect_error(
    get_alternate_alleles("X", observed),
    "Unknown character"
  )
})

test_that("get_alternate_alleles ignores zero-count alleles", {
  # Reference is A, have C but not G/T/I/D
  observed <- c(100, 50, 0, 0, 0, 0)
  result <- get_alternate_alleles("A", observed)

  # Only C should be TRUE
  expect_equal(sum(result), 1)
  expect_true(result[2]) # C position
})
