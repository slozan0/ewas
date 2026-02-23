# Regression tests - ensure refactoring doesn't change results
# This tests the entire pipeline with real data

source(here("scripts", "functions", "easy_chi2_fun.r"))

sample_file <- here("data", "sample", "5feb_tem_chr3_avd_sample.rds")

test_that("get_easy_chi_estimates produces consistent results", {
  # Use sample data (committed to git)

  skip_if_not(file.exists(sample_file),
    message = paste(
      "Sample data not found. Run:",
      "source('scripts/create_sample_data.r')"
    )
  )

  # Load sample data
  poly_sites <- readRDS(sample_file)

  # Take just the first row and apply the same preprocessing as easy_chi2.r
  poly_site <- poly_sites[1, ]
  cols_to_remove <- c(
    "refnuc1", "refnuc2", "refnuc3", "refnuc4",
    "chrom1", "chrom2", "chrom3", "chrom4",
    "sumDepth1", "sumDepth2", "sumDepth3", "sumDepth4"
  )
  poly_site <- poly_site[, !(names(poly_site) %in% cols_to_remove), drop = FALSE]

  # Run the function
  result <- get_easy_chi_estimates(poly_site)

  # Check structure
  expect_type(result, "double")
  expect_named(result, c(
    "nuc_position", "ref_nuc", "group1_alt_all_freq", "group2_alt_all_freq",
    "lod", "group1_heteroz", "group2_heteroz", "total_heteroz",
    "a_s", "c_s", "g_s", "t_s", "i_s", "d_s",
    "group1_chi_sqr", "group2_chi_sqr", "total_chi_sqr",
    "group1_deg_freedom", "group2_deg_freedom", "total_deg_freedom"
  ))

  # Check value ranges
  expect_gte(result["lod"], 0)
  expect_gte(result["group1_heteroz"], 0)
  expect_lte(result["group1_heteroz"], 1)
  expect_gte(result["group2_heteroz"], 0)
  expect_lte(result["group2_heteroz"], 1)

  # Save the result for future comparison (first time only)
  golden_file <- here("tests", "testthat", "fixtures",
                      "ezchi_golden_output.rds")

  if (!file.exists(golden_file)) {
    dir.create(here("tests", "testthat", "fixtures"),
      showWarnings = FALSE, recursive = TRUE
    )
    saveRDS(result, golden_file)
    skip("Creating golden output file for future comparison")
  }

  # Compare to golden output
  golden <- readRDS(golden_file)
  expect_equal(result, golden, tolerance = 1e-6)
})
