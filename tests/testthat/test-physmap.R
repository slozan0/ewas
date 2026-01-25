# Test for physmap.r script
# Tests the complete physmap pipeline with sample data
library(here)

test_that("physmap processes sample data correctly", {
  # Skip if C++ code not compiled
  # Install "BH" library if you haven't install.packages("BH")

  cpp_file <- here("scripts", "functions", "physmap.cpp")

  skip_if_not(file.exists(cpp_file),
    message = "physmap.cpp not found"
  )

  # Skip if sample data doesn't exist
  sample_file <- here("data", "sample", "5feb_tem_a1_c3_sample.rds")
  skip_if_not(file.exists(sample_file),
    message =
      "Sample data not found. Run: source('scripts/create_sample_data.r')"
  )

  # Load required libraries
  library(Rcpp)
  library(data.table)

  # Compile C++ code (if not already compiled)
  # This is safe to run multiple times
  suppressMessages(sourceCpp(cpp_file))

  # Load sample data
  raw_data <- readRDS(sample_file)

  # Run physmap processing (same logic as physmap.r)
  mono_sites <- raw_data[raw_data$snp1 == "", ]
  mono_sites <- as.matrix(mono_sites)

  poly_sites <- raw_data[raw_data$snp1 != "", ]
  poly_sites <- as.matrix(poly_sites)

  # Process sites
  my_mono_lines <- MapMonoSites(mono_sites)
  my_poly_lines <- MapPolySites(poly_sites)

  # Convert to data.table
  dt_mono_sites <- data.table(
    chrom = as.numeric(my_mono_lines[, 1]),
    pos = as.numeric(my_mono_lines[, 2]),
    ref = as.character(my_mono_lines[, 3]),
    a = as.numeric(my_mono_lines[, 4]),
    c = as.numeric(my_mono_lines[, 5]),
    g = as.numeric(my_mono_lines[, 6]),
    t = as.numeric(my_mono_lines[, 7]),
    i = as.numeric(my_mono_lines[, 8]),
    d = as.numeric(my_mono_lines[, 9]),
    key = c("chrom", "pos")
  )
  dt_mono_sites[is.na(dt_mono_sites)] <- 0

  dt_poly_sites <- data.table(
    chrom = as.numeric(my_poly_lines[, 1]),
    pos = as.numeric(my_poly_lines[, 2]),
    ref = as.character(my_poly_lines[, 3]),
    a = as.numeric(my_poly_lines[, 4]),
    c = as.numeric(my_poly_lines[, 5]),
    g = as.numeric(my_poly_lines[, 6]),
    t = as.numeric(my_poly_lines[, 7]),
    i = as.numeric(my_poly_lines[, 8]),
    d = as.numeric(my_poly_lines[, 9]),
    key = c("chrom", "pos")
  )
  dt_poly_sites[is.na(dt_poly_sites)] <- 0

  # Combine
  dt_chrom <- rbind(dt_poly_sites, dt_mono_sites)

  # Add sumDepth
  dt_chrom <- dt_chrom[, sumDepth := sum(a, c, g, t, i, d),
    by = seq_len(NROW(dt_chrom))
  ]

  # Filter by depth
  dt_chrom <- dt_chrom[sumDepth >= 25]
  dt_chrom <- dt_chrom[sumDepth < 1000]

  # === TESTS START HERE ===

  # Test 1: Output structure
  expect_s3_class(dt_chrom, "data.table")
  expect_named(dt_chrom, c("chrom", "pos", "ref",
                           "a", "c", "g", "t", "i", "d", "sumDepth"))

  # Test 2: Data types
  expect_type(dt_chrom$chrom, "double")
  expect_type(dt_chrom$pos, "double")
  expect_type(dt_chrom$ref, "character")
  expect_type(dt_chrom$sumDepth, "double")

  # Test 3: Depth filtering worked
  expect_true(all(dt_chrom$sumDepth >= 25),
    info = "All rows should have sumDepth >= 25"
  )
  expect_true(all(dt_chrom$sumDepth < 1000),
    info = "All rows should have sumDepth < 1000"
  )

  # Test 4: sumDepth calculation is correct
  calculated_sum <- dt_chrom$a + dt_chrom$c + dt_chrom$g +
    dt_chrom$t + dt_chrom$i + dt_chrom$d
  expect_equal(dt_chrom$sumDepth, calculated_sum,
    info = "sumDepth should equal sum of nucleotide counts"
  )

  # Test 5: No NAs in critical columns
  expect_false(any(is.na(dt_chrom$pos)),
    info = "Position column should have no NAs"
  )
  expect_false(any(is.na(dt_chrom$ref)),
    info = "Reference column should have no NAs"
  )

  # Test 6: Reference nucleotides are valid
  valid_refs <- c("A", "C", "G", "T")
  expect_true(all(dt_chrom$ref %in% valid_refs),
    info = "All reference nucleotides should be A, C, G, or T"
  )

  # Test 7: Positions are in expected range (1M-2M for sample data)
  expect_true(all(dt_chrom$pos >= 1e6 & dt_chrom$pos <= 2e6),
    info = "Sample data positions should be in 1M-2M range"
  )

  # Test 8: Golden file comparison (regression test)
  golden_file <- here("tests", "testthat", "fixtures",
                      "physmap_golden_output.rds")

  if (!file.exists(golden_file)) {
    # First run - create golden file
    dir.create(here("tests", "testthat", "fixtures"),
               showWarnings = FALSE, recursive = TRUE)
    saveRDS(dt_chrom, golden_file)
    skip("Created golden output file for future comparison")
  } else {
    # Compare to golden file
    golden_output <- readRDS(golden_file)

    expect_equal(nrow(dt_chrom), nrow(golden_output),
      info = "Number of rows should match golden output"
    )

    expect_equal(dt_chrom, golden_output,
      info = "Output should exactly match golden output",
      tolerance = 1e-10
    )
  }
})

test_that("physmap handles edge cases", {
  cpp_file <- here("scripts", "functions", "physmap.cpp")

  skip_if_not(file.exists(cpp_file),
              message = "physmap.cpp not found")

  library(Rcpp)
  library(data.table)
  suppressMessages(sourceCpp(cpp_file))

  # Test with minimal data structure
  minimal_data <- data.frame(
    chrom = "3",
    position = 1000000,
    refnuc = "A",
    depth = 50,
    q30_depth = 48,
    refQA = "A:48:1:36:1:48:0:0",
    snp1 = "",
    snp2 = "",
    snp3 = "",
    snp4 = "",
    snp5 = "",
    snp6 = "",
    snp7 = "",
    snp8 = "",
    snp9 = "",
    snp10 = ""
  )

  # Should not crash
  expect_error(
    {
      mono_sites <- as.matrix(minimal_data)
      result <- MapMonoSites(mono_sites)
    },
    NA
  )
})
