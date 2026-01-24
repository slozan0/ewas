# Tests for get_chi function
source("../../scripts/functions/easy_chi2_fun.r")

test_that("get_chi returns valid chi-square and degrees of freedom", {
  # Simple case: two alleles, observed vs expected
  observed <- c(50, 50, 0, 0, 0, 0, 0)
  w_obs1 <- c(30, 20, 0, 0, 0, 0, 0)
  w_obs2 <- c(20, 30, 0, 0, 0, 0, 0)

  result <- get_chi(
    nuc_position = 12345,
    observed = observed,
    w_obs1 = w_obs1,
    w_obs2 = w_obs2
  )

  # Should return a list with chiSqr and degFreedom
  expect_type(result, "list")
  expect_named(result, c("chiSqr", "degFreedom"))

  # Chi-square should be non-negative
  expect_gte(result$chiSqr, 0)

  # Degrees of freedom should be (n_alleles - 1)
  expect_equal(result$degFreedom, 1) # 2 alleles - 1 = 1
})

test_that("get_chi handles perfect fit", {
  # When observed matches expected, chi-square should be ~0
  observed <- c(50, 50, 0, 0, 0, 0, 0)
  w_obs1 <- c(25, 25, 0, 0, 0, 0, 0)
  w_obs2 <- c(25, 25, 0, 0, 0, 0, 0)

  result <- get_chi(
    nuc_position = 12345,
    observed = observed,
    w_obs1 = w_obs1,
    w_obs2 = w_obs2
  )

  # Chi-square should be very small (near 0)
  expect_lt(result$chiSqr, 0.01)
})

test_that("get_chi handles multiple alleles", {
  # Three alleles
  observed <- c(30, 30, 30, 0, 0, 0, 0)
  w_obs1 <- c(15, 15, 15, 0, 0, 0, 0)
  w_obs2 <- c(15, 15, 15, 0, 0, 0, 0)

  result <- get_chi(
    nuc_position = 12345,
    observed = observed,
    w_obs1 = w_obs1,
    w_obs2 = w_obs2
  )

  # Degrees of freedom should be 2 (3 alleles - 1)
  expect_equal(result$degFreedom, 2)
})
