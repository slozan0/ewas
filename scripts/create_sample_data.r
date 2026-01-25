# Create Sample Data for Testing
# This script creates small sample datasets from real data for testing/sharing

rm(list = ls())

library(data.table)

# Configuration ----
# Position range to sample (there isn't much data before 1M)
pos_start <- 1000000
pos_end <- 2000000
position_col <- "position"  # Column name for nucleotide position

input_dir <- "data/output"  # These are the split chromosome files
sample_dir <- "data/sample"

# Create sample directory if it doesn't exist
if (!dir.exists(sample_dir)) {
  dir.create(sample_dir, recursive = TRUE)
}

# Function to create sample from RDS file by position range ----
create_sample_by_position <- function(input_file, output_file,
                                      pos_start, pos_end,
                                      pos_col = "position") {
  if (!file.exists(input_file)) {
    warning(sprintf("Input file not found: %s", input_file))
    return(FALSE)
  }

  cat(sprintf("Reading: %s\n", input_file))
  full_data <- readRDS(input_file)

  # Convert to data.table if needed
  if (!is.data.table(full_data)) {
    setDT(full_data)
  }

  # Check if position column exists
  if (!pos_col %in% names(full_data)) {
    warning(sprintf(
      "Position column '%s' not found in %s. Available columns: %s",
      pos_col,
      basename(input_file),
      paste(names(full_data), collapse = ", ")
    ))
    return(FALSE)
  }

  cat(sprintf(
    "  Total rows: %s\n",
    format(nrow(full_data), big.mark = ",")
  ))
  cat(sprintf(
    "  Position range in file: %s - %s\n",
    format(min(full_data[[pos_col]], na.rm = TRUE), big.mark = ","),
    format(max(full_data[[pos_col]], na.rm = TRUE), big.mark = ",")
  ))

  # Filter by position range
  sample_data <- full_data[
    get(pos_col) >= pos_start & get(pos_col) <= pos_end
  ]

  cat(sprintf(
    "  Filtered to positions %s - %s: %s rows\n",
    format(pos_start, big.mark = ","),
    format(pos_end, big.mark = ","),
    format(nrow(sample_data), big.mark = ",")
  ))

  if (nrow(sample_data) == 0) {
    warning("No data in specified position range!")
    return(FALSE)
  }

  # Save sample
  cat(sprintf("Saving to: %s\n", output_file))
  saveRDS(sample_data, output_file, compress = "xz")

  # Report size
  size_mb <- file.info(output_file)$size / 1024^2
  cat(sprintf("Sample file size: %.2f MB\n\n", size_mb))

  if (size_mb > 10) {
    warning(sprintf(
      "Sample file is large (%.2f MB). Consider narrowing position range.",
      size_mb
    ))
  }

  return(TRUE)
}

# Create sample data from chromosome 3 files ----
# These are the 4 files you specified
sample_files <- list(
  list(
    input = file.path(input_dir, "5feb_tem_a1_c3_rc.rds"),
    output = file.path(sample_dir, "5feb_tem_a1_c3_sample.rds")
  ),
  list(
    input = file.path(input_dir, "5feb_tem_a2_c3_rc.rds"),
    output = file.path(sample_dir, "5feb_tem_a2_c3_sample.rds")
  ),
  list(
    input = file.path(input_dir, "5feb_tem_d1_c3_rc.rds"),
    output = file.path(sample_dir, "5feb_tem_d1_c3_sample.rds")
  ),
  list(
    input = file.path(input_dir, "5feb_tem_d2_c3_rc.rds"),
    output = file.path(sample_dir, "5feb_tem_d2_c3_sample.rds")
  )
)

cat("Creating sample data from chromosome 3 files\n")
cat(sprintf("Position range: %s to %s\n\n",
  format(pos_start, big.mark = ","),
  format(pos_end, big.mark = ",")
))

success_count <- 0
for (file_pair in sample_files) {
  if (create_sample_by_position(
    file_pair$input,
    file_pair$output,
    pos_start = pos_start,
    pos_end = pos_end,
    pos_col = position_col
  )) {
    success_count <- success_count + 1
  }
}

# Summary ----
cat("\n")
cat("========================================\n")
cat("Sample Data Creation Summary\n")
cat("========================================\n")
cat(sprintf(
  "Successfully created: %d / %d files\n",
  success_count,
  length(sample_files)
))
cat(sprintf(
  "Position range: %s - %s\n",
  format(pos_start, big.mark = ","),
  format(pos_end, big.mark = ",")
))

# List created files and sizes
cat("\nCreated files:\n")
sample_files_created <- list.files(
  sample_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)
if (length(sample_files_created) > 0) {
  for (f in sample_files_created) {
    size_mb <- file.info(f)$size / 1024^2
    cat(sprintf("  %s (%.2f MB)\n", basename(f), size_mb))
  }
  total_size <- sum(sapply(
    sample_files_created,
    function(f) file.info(f)$size
  )) / 1024^2
  cat(sprintf("\nTotal size: %.2f MB\n", total_size))
} else {
  cat("  No sample files created.\n")
}

cat("\nNext steps:\n")
cat("1. Review sample files in data/sample/\n")
cat("2. Verify file sizes are reasonable (<10 MB each)\n")
cat("3. Test with: source('run_tests.R')\n")
cat("4. Commit sample data to git:\n")
cat("   git add data/sample/\n")
cat("   git commit -m 'Add sample test data (1M-2M)'\n")
cat("========================================\n")
