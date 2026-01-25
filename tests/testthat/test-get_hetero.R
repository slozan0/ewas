# Tests for get_hetero function
source("../../scripts/functions/easy_chi2_fun.r")

test_that("get_hetero returns value between 0 and 1", {
  # Simple test case: all same nucleotide (homozygous)
  counts_homo <- matrix(c(100, 0, 0, 0, 0, 0), nrow = 1)
  result <- get_hetero(counts_homo, n_groups = 1)

  expect_gte(result, 0)
  expect_lte(result, 1)
  expect_equal(result, 0) # Homozygous should be 0
})

test_that("get_hetero handles perfect heterozygosity", {
  # Equal distribution across all 6 nucleotides
  counts_hetero <- matrix(c(10, 10, 10, 10, 10, 10), nrow = 1)
  result <- get_hetero(counts_hetero, n_groups = 1)

  expect_gte(result, 0)
  expect_lte(result, 1)
  expect_gt(result, 0.8) # Should be high heterozygosity
})

test_that("get_hetero handles two alleles", {
  # Two alleles with 50/50 split
  counts_two <- matrix(c(50, 50, 0, 0, 0, 0), nrow = 1)
  result <- get_hetero(counts_two, n_groups = 1)

  expected <- 1 - (0.5^2 + 0.5^2) # 1 - (0.25 + 0.25) = 0.5
  expect_equal(result, expected)
})

test_that("get_hetero handles multiple groups", {
  # Multiple groups (rows)
  counts_multi <- matrix(
    c(50, 0, 0, 0, 0, 0,
      0, 50, 0, 0, 0, 0),
    nrow = 2, byrow = TRUE
  )
  result <- get_hetero(counts_multi, n_groups = 2)

  # Combined: 50 A + 50 C = 0.5 each
  expected <- 1 - (0.5^2 + 0.5^2)
  expect_equal(result, expected)
})

test_that("get_hetero handles edge case: all zeros", {
  counts_zero <- matrix(c(0, 0, 0, 0, 0, 0), nrow = 1)

  # Division by zero should produce NaN
  expect_warning(result <- get_hetero(counts_zero, n_groups = 1))
  expect_true(is.nan(result) || is.na(result))
})
