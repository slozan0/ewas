# Regression tests - ensure refactoring doesn't change results
# This tests the entire pipeline with real data

source("../../scripts/functions/easy_chi2_fun.r")

test_that("get_easy_chi_estimates produces consistent results", {
  # Use sample data (committed to git)
  sample_file <- "data/sample/test_chr1.rds"

  skip_if_not(file.exists(sample_file),
    message = paste(
      "Sample data not found. Run:",
      "source('scripts/create_sample_data.r')"
    )
  )

  # Load sample data
  poly_sites <- readRDS(sample_file)

  # Take just the first row
  poly_site <- poly_sites[1, ]
  poly_site$ref <- poly_site$ref1

  # Run the function
  result <- get_easy_chi_estimates(poly_site)

  # Check structure
  expect_type(result, "double")
  expect_named(result, c(
    "nucPosition", "refNuc", "group1AltAllFreq", "group2AltAllFreq",
    "lod", "group1Heteroz", "group2Heteroz", "totalHeteroz",
    "As", "Cs", "Gs", "Ts", "Is", "Ds",
    "group1ChiSqr", "group2ChiSqr", "totalChiSqr",
    "group1DegFreedom", "group2DegFreedom", "totalDegFreedom"
  ))

  # Check value ranges
  expect_gte(result["lod"], 0)
  expect_gte(result["group1Heteroz"], 0)
  expect_lte(result["group1Heteroz"], 1)
  expect_gte(result["group2Heteroz"], 0)
  expect_lte(result["group2Heteroz"], 1)

  # Save the result for future comparison (first time only)
  golden_file <- "tests/testthat/golden-output.rds"
  if (!file.exists(golden_file)) {
    saveRDS(result, golden_file)
    skip("Creating golden output file for future comparison")
  }

  # Compare to golden output
  golden <- readRDS(golden_file)
  expect_equal(result, golden, tolerance = 1e-6)
})
