# Helper Script to Create/Update Golden Files
# Run this when you need to create or update expected test outputs

rm(list = ls())

library(testthat)
library(here)

cat("\n")
cat("========================================\n")
cat("Golden File Manager\n")
cat("========================================\n")
cat("\n")

# Create fixtures directory if it doesn't exist
fixtures_dir <- here("tests", "testthat", "fixtures")
if (!dir.exists(fixtures_dir)) {
  dir.create(fixtures_dir, recursive = TRUE)
  cat("Created fixtures directory\n\n")
}

# === Function to create golden files ===

create_physmap_golden <- function() {
  cat("Creating physmap golden file...\n")

  # Check if sample data exists
  sample_file <- here("data", "sample", "5feb_tem_a1_c3_sample.rds")
  if (!file.exists(sample_file)) {
    cat("✗ Sample data not found:", sample_file, "\n")
    cat("  Run: source('scripts/create_sample_data.r')\n\n")
    return(FALSE)
  }

  # Check if C++ code exists
  cpp_file <- here("scripts", "functions", "physmap.cpp")
  if (!file.exists(cpp_file)) {
    cat("✗ physmap.cpp not found\n\n")
    return(FALSE)
  }

  # Load libraries
  library(Rcpp)
  library(data.table)

  # Compile C++ code
  cat("  Compiling C++ code...\n")
  suppressMessages(sourceCpp(cpp_file))

  # Run physmap logic
  cat("  Processing sample data...\n")
  raw_data <- readRDS(sample_file)

  mono_sites <- raw_data[raw_data$snp1 == "", ]
  mono_sites <- as.matrix(mono_sites)

  poly_sites <- raw_data[raw_data$snp1 != "", ]
  poly_sites <- as.matrix(poly_sites)

  # nolint start: object_usage_linter (C++ functions)
  my_mono_lines <- MapMonoSites(mono_sites)
  my_poly_lines <- MapPolySites(poly_sites)
  # nolint end

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

  dt_chrom <- rbind(dt_poly_sites, dt_mono_sites)

  # nolint start: object_usage_linter
  dt_chrom <- dt_chrom[, sumDepth := sum(a, c, g, t, i, d),
    by = seq_len(NROW(dt_chrom))
  ]

  dt_chrom <- dt_chrom[sumDepth >= 25]
  dt_chrom <- dt_chrom[sumDepth < 1000]
  # nolint end

  # Save golden file
  golden_file <- here("tests", "testthat", "fixtures",
                      "physmap_golden_output.rds")

  # Check if file already exists
  if (file.exists(golden_file)) {
    cat("  ⚠️  Golden file already exists!\n")
    cat("     Compare old vs new before overwriting\n\n")

    old_golden <- readRDS(golden_file)

    cat("  Old golden file:\n")
    cat(sprintf("    Rows: %d\n", nrow(old_golden)))
    cat(sprintf("    Columns: %s\n", paste(names(old_golden), collapse = ", ")))

    cat("\n  New output:\n")
    cat(sprintf("    Rows: %d\n", nrow(dt_chrom)))
    cat(sprintf("    Columns: %s\n", paste(names(dt_chrom), collapse = ", ")))

    if (identical(old_golden, dt_chrom)) {
      cat("\n  ✓ New output is identical to existing golden file\n")
      cat("    No update needed\n\n")
      return(TRUE)
    } else {
      cat("\n  ✗ New output differs from existing golden file\n")
      response <- readline("  Overwrite? (yes/no): ")

      if (tolower(response) != "yes") {
        cat("  Cancelled. Golden file not updated.\n\n")
        return(FALSE)
      }
    }
  }

  saveRDS(dt_chrom, golden_file, compress = "xz")

  cat("\n✓ Golden file created:\n")
  cat(" ", golden_file, "\n")
  cat(sprintf("  Rows: %d\n", nrow(dt_chrom)))
  cat(sprintf("  Size: %.2f KB\n",
    file.info(golden_file)$size / 1024
  ))
  cat("\n")

  TRUE
}

# === Main Menu ===

cat("What would you like to do?\n\n")
cat("1. Create physmap golden file\n")
cat("2. List existing golden files\n")
cat("3. Remove all golden files (start fresh)\n")
cat("4. Exit\n\n")

choice <- readline("Enter choice (1-4): ")

if (choice == "1") {
  create_physmap_golden()
} else if (choice == "2") {
  cat("\nExisting golden files:\n")
  golden_files <- list.files(
    fixtures_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )

  if (length(golden_files) == 0) {
    cat("  No golden files found\n")
  } else {
    for (f in golden_files) {
      size_kb <- file.info(f)$size / 1024
      cat(sprintf("  %s (%.2f KB)\n", basename(f), size_kb))
    }
  }
  cat("\n")
} else if (choice == "3") {
  golden_files <- list.files(fixtures_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )

  if (length(golden_files) == 0) {
    cat("\nNo golden files to remove\n\n")
  } else {
    cat("\nThis will delete:\n")
    for (f in golden_files) {
      cat("  -", basename(f), "\n")
    }

    response <- readline("\nAre you sure? (yes/no): ")

    if (tolower(response) == "yes") {
      unlink(golden_files)
      cat("✓ Golden files removed\n\n")
    } else {
      cat("Cancelled\n\n")
    }
  }
} else if (choice == "4") {
  cat("\nExiting\n\n")
} else {
  cat("\nInvalid choice\n\n")
}

cat("Done!\n")
cat("========================================\n\n")
