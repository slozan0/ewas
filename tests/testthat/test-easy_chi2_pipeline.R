# Test for easy_chi2.r - Full Pipeline Test
# Tests the complete chi-square analysis pipeline with sample data

library(here)

test_that("easy_chi2 pipeline processes sample data correctly", {
  # Skip if sample data doesn't exist
  sample_file <- here("data", "sample", "5feb_tem_chr3_avd_sample.rds")
  skip_if_not(file.exists(sample_file),
    message = "Sample data not found. Run:
    source('scripts/create_sample_data.r')"
  )

  # Load required libraries
  library(data.table)
  source(here("scripts", "functions", "easy_chi2_fun.r"))

  # Load sample data
  cat("Loading sample data...\n")
  poly_sites <- readRDS(sample_file)

  # Test on first 100 rows (faster, but still comprehensive)
  n_test_rows <- min(100, nrow(poly_sites))
  test_data <- poly_sites[1:n_test_rows, ]

  cat(sprintf("Processing %d test rows...\n", n_test_rows))

  # Process sequentially (no parallel for testing)
  results_list <- list()
  for (i in 1:n_test_rows) {
    results_list[[i]] <- get_easy_chi_estimates(poly_site = test_data[i, ])
  }

  # Combine results
  raw_results <- do.call(rbind, results_list)
  ezchi_results <- as.data.table(raw_results, key = "nuc_position")

  # === TESTS START HERE ===

  # Test 1: Output structure
  expect_s3_class(ezchi_results, "data.table")

  expected_cols <- c(
    "nuc_position", "ref_nuc", "group1_alt_all_freq", "group2_alt_all_freq",
    "lod", "group1_heteroz", "group2_heteroz", "total_heteroz",
    "a_s", "c_s", "g_s", "t_s", "i_s", "d_s",
    "group1_chi_sqr", "group2_chi_sqr", "total_chi_sqr",
    "group1_deg_freedom", "group2_deg_freedom", "total_deg_freedom"
  )

  expect_true(
    all(expected_cols %in% names(ezchi_results)),
    info = "All expected columns should be present"
  )

  # Test 2: Data types
  expect_type(ezchi_results$nuc_position, "double")
  expect_type(ezchi_results$ref_nuc, "double")
  expect_type(ezchi_results$lod, "double")
  expect_type(ezchi_results$group1_heteroz, "double")

  # Test 3: Value ranges - LOD scores
  expect_true(
    all(ezchi_results$lod >= 0),
    info = "All LOD scores should be non-negative"
  )

  # Test 4: Value ranges - Heterozygosity (0 to 1)
  expect_true(
    all(ezchi_results$group1_heteroz >= 0 &
      ezchi_results$group1_heteroz <= 1),
    info = "Group 1 heterozygosity should be between 0 and 1"
  )
  expect_true(
    all(ezchi_results$group2_heteroz >= 0 &
      ezchi_results$group2_heteroz <= 1),
    info = "Group 2 heterozygosity should be between 0 and 1"
  )
  expect_true(
    all(ezchi_results$total_heteroz >= 0 &
      ezchi_results$total_heteroz <= 1),
    info = "Total heterozygosity should be between 0 and 1"
  )

  # Test 5: Value ranges - Frequencies (0 to 1)
  expect_true(
    all(ezchi_results$group1_alt_all_freq >= 0 &
      ezchi_results$group1_alt_all_freq <= 1),
    info = "Group 1 allele frequencies should be between 0 and 1"
  )
  expect_true(
    all(ezchi_results$group2_alt_all_freq >= 0 &
      ezchi_results$group2_alt_all_freq <= 1),
    info = "Group 2 allele frequencies should be between 0 and 1"
  )

  # Test 6: Chi-square values are non-negative
  expect_true(
    all(ezchi_results$group1_chi_sqr >= 0),
    info = "Group 1 chi-square should be non-negative"
  )
  expect_true(
    all(ezchi_results$group2_chi_sqr >= 0),
    info = "Group 2 chi-square should be non-negative"
  )
  expect_true(
    all(ezchi_results$total_chi_sqr >= 0),
    info = "Total chi-square should be non-negative"
  )

  # Test 7: Degrees of freedom are reasonable
  expect_true(
    all(ezchi_results$group1_deg_freedom >= 0 &
      ezchi_results$group1_deg_freedom <= 5),
    info = "Degrees of freedom should be 0-5 (max 6 alleles - 1)"
  )
  expect_true(
    all(ezchi_results$group2_deg_freedom >= 0 &
      ezchi_results$group2_deg_freedom <= 5),
    info = "Degrees of freedom should be 0-5"
  )

  # Test 8: Nucleotide counts are non-negative integers
  expect_true(
    all(ezchi_results$a_s >= 0),
    info = "A counts should be non-negative"
  )
  expect_true(
    all(ezchi_results$c_s >= 0),
    info = "C counts should be non-negative"
  )
  expect_true(
    all(ezchi_results$g_s >= 0),
    info = "G counts should be non-negative"
  )
  expect_true(
    all(ezchi_results$t_s >= 0),
    info = "T counts should be non-negative"
  )

  # Test 9: Total nucleotide count makes sense
  total_nuc <- ezchi_results$a_s + ezchi_results$c_s +
    ezchi_results$g_s + ezchi_results$t_s +
    ezchi_results$i_s + ezchi_results$d_s

  expect_true(
    all(total_nuc > 0),
    info = "Each position should have at least some nucleotides"
  )

  # Test 10: Reference nucleotide is valid
  valid_refs <- c(1, 2, 3, 4) # A=1, C=2, G=3, T=4
  expect_true(
    all(ezchi_results$ref_nuc %in% valid_refs),
    info = "Reference nucleotide should be 1-4 (A, C, G, T)"
  )

  # Test 11: Positions are in expected range (1M-2M for sample data)
  expect_true(
    all(ezchi_results$nuc_position >= 1e6 &
      ezchi_results$nuc_position <= 2e6),
    info = "Sample data positions should be in 1M-2M range"
  )

  # Test 12: No NAs in critical columns
  expect_false(any(is.na(ezchi_results$nuc_position)),
    info = "No NAs in nuc_position"
  )
  expect_false(any(is.na(ezchi_results$lod)),
    info = "No NAs in LOD scores"
  )
  expect_false(any(is.na(ezchi_results$ref_nuc)),
    info = "No NAs in reference nucleotide"
  )
})

test_that("easy_chi2 probability calculations are correct", {
  sample_file <- here("data", "sample", "5feb_tem_chr3_avd_sample.rds")
  skip_if_not(file.exists(sample_file), message = "Sample data not found")

  library(data.table)
  source(here("scripts", "functions", "easy_chi2_fun.r"))

  poly_sites <- readRDS(sample_file)

  # Test first row
  result <- get_easy_chi_estimates(poly_site = poly_sites[1, ])
  result_dt <- as.data.table(t(result))

  # Calculate probability from chi-square and df
  total_prob <- 1 - pchisq(
    q = result_dt$total_chi_sqr,
    df = result_dt$total_deg_freedom
  )

  group1_prob <- 1 - pchisq(
    q = result_dt$group1_chi_sqr,
    df = result_dt$group1_deg_freedom
  )

  group2_prob <- 1 - pchisq(
    q = result_dt$group2_chi_sqr,
    df = result_dt$group2_deg_freedom
  )

  # Test probability calculations
  expect_true(total_prob >= 0 && total_prob <= 1,
    info = "Total probability should be between 0 and 1"
  )
  expect_true(group1_prob >= 0 && group1_prob <= 1,
    info = "Group 1 probability should be between 0 and 1"
  )
  expect_true(group2_prob >= 0 && group2_prob <= 1,
    info = "Group 2 probability should be between 0 and 1"
  )

  # LOD should be -log10(p-value) for chi < 30
  if (result_dt$total_chi_sqr < 30) {
    expected_lod <- -log10(total_prob)
    expect_equal(result_dt$lod, expected_lod,
      tolerance = 1e-6,
      info = "LOD calculation should match expected value"
    )
  }
})

test_that("Benjamini-Hochberg correction works correctly", {
  sample_file <- here("data", "sample", "5feb_tem_chr3_avd_sample.rds")
  skip_if_not(file.exists(sample_file), message = "Sample data not found")

  library(data.table)
  source(here("scripts", "functions", "easy_chi2_fun.r"))

  # Test with known p-values
  test_pvalues <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.9)

  bh_thresh <- get_benjamini_hochber_thresh(test_pvalues)

  # Threshold should be a probability between 0 and 1
  expect_true(bh_thresh >= 0 && bh_thresh <= 1,
    info = "BH threshold should be a valid probability"
  )

  # Should be one of the input p-values
  expect_true(bh_thresh %in% test_pvalues,
    info = "BH threshold should be one of the input p-values"
  )

  # Test with uniform p-values (edge case - no p-values pass BH at alpha=0.01)
  # seq(0.01, 0.99) values all exceed their j_alpha = i * 0.01/n thresholds,
  # so the function correctly returns numeric(0) (nothing significant).
  uniform_pvals <- seq(0.01, 0.99, length.out = 100)
  bh_thresh_uniform <- get_benjamini_hochber_thresh(uniform_pvals)

  expect_equal(
    length(bh_thresh_uniform), 0,
    info = "BH threshold should be empty when no p-values pass the threshold"
  )
})

test_that("mark_inconsistency function works correctly", {
  source(here("scripts", "functions", "easy_chi2_fun.r"))

  # Test case 1: No inconsistency (high p-values)
  mark1 <- mark_inconsistency(
    chi1 = 0.5, deg_freedom1 = 1, inconsistency_mark1 = "1*",
    chi2 = 0.3, deg_freedom2 = 1, inconsistency_mark2 = "2*",
    mark_threshold = 0.05
  )
  expect_equal(mark1, "",
    info = "No mark for non-significant differences"
  )

  # Test case 2: Group 1 inconsistent (low p-value)
  mark2 <- mark_inconsistency(
    chi1 = 10, deg_freedom1 = 1, inconsistency_mark1 = "1*",
    chi2 = 0.3, deg_freedom2 = 1, inconsistency_mark2 = "2*",
    mark_threshold = 0.05
  )
  expect_equal(mark2, "1*",
    info = "Should mark group 1 inconsistency"
  )

  # Test case 3: Group 2 inconsistent
  mark3 <- mark_inconsistency(
    chi1 = 0.3, deg_freedom1 = 1, inconsistency_mark1 = "1*",
    chi2 = 10, deg_freedom2 = 1, inconsistency_mark2 = "2*",
    mark_threshold = 0.05
  )
  expect_equal(mark3, "2*",
    info = "Should mark group 2 inconsistency"
  )
})

test_that("get_alleles_label creates correct labels", {
  source(here("scripts", "functions", "easy_chi2_fun.r"))

  # Test case 1: A reference with C and G alternates
  label1 <- get_alleles_label(
    nuc_position = 1000000,
    ref_nucleotide = 1, # A
    a_s = 100, c_s = 50, g_s = 30, t_s = 0, i_s = 0, d_s = 0
  )
  expect_true(grepl("^A", label1),
    info = "Label should start with reference allele A"
  )
  expect_true(grepl("C", label1),
    info = "Label should contain alternate allele C"
  )
  expect_true(grepl("G", label1),
    info = "Label should contain alternate allele G"
  )
  expect_false(grepl("T", label1),
    info = "Label should not contain T (count = 0)"
  )

  # Test case 2: T reference with A alternate
  label2 <- get_alleles_label(
    nuc_position = 1000001,
    ref_nucleotide = 4, # T
    a_s = 20, c_s = 0, g_s = 0, t_s = 80, i_s = 0, d_s = 0
  )
  expect_true(grepl("^T", label2),
    info = "Label should start with reference allele T"
  )
  expect_true(grepl("A", label2),
    info = "Label should contain alternate allele A"
  )
})

test_that("Golden file comparison - full pipeline", {
  sample_file <- here("data", "sample", "5feb_tem_chr3_avd_sample.rds")
  skip_if_not(file.exists(sample_file), message = "Sample data not found")

  library(data.table)
  source(here("scripts", "functions", "easy_chi2_fun.r"))

  # Process first 5 rows for golden file
  poly_sites <- readRDS(sample_file)

  results_list <- list()
  for (i in 1:5) {
    results_list[[i]] <- get_easy_chi_estimates(poly_site = poly_sites[i, ])
  }

  raw_results <- do.call(rbind, results_list)
  ezchi_results <- as.data.table(raw_results, key = "nuc_position")

  # Add probability calculations
  ezchi_results[, total_prob := 1 -
    pchisq(q = total_chi_sqr, df = total_deg_freedom),
  by = nuc_position
  ]

  golden_file <- here(
    "tests", "testthat", "fixtures",
    "ezchi_pipeline_golden.rds"
  )

  if (!file.exists(golden_file)) {
    dir.create(here("tests", "testthat", "fixtures"),
      showWarnings = FALSE, recursive = TRUE
    )
    saveRDS(ezchi_results, golden_file)
    skip("Created golden output file for future comparison")
  }

  # Compare to golden file
  golden_output <- readRDS(golden_file)

  expect_equal(nrow(ezchi_results), nrow(golden_output),
    info = "Number of rows should match golden output"
  )

  expect_equal(ezchi_results, golden_output,
    info = "Full pipeline output should match golden output",
    tolerance = 1e-8
  )
})
