# Quick Check - Load and Compare 4 Files
# Fast overview of file structure and position ranges

rm(list = ls())
library(data.table)

# Load all 4 files ----
cat("Loading files...\n\n")

files <- list(
  a1 = "data/output/5feb_tem_a1_c3_rc.rds",
  a2 = "data/output/5feb_tem_a2_c3_rc.rds",
  d1 = "data/output/5feb_tem_d1_c3_rc.rds",
  d2 = "data/output/5feb_tem_d2_c3_rc.rds"
)

data_list <- list()

for (name in names(files)) {
  cat(sprintf("Loading %s... ", name))
  if (file.exists(files[[name]])) {
    data_list[[name]] <- readRDS(files[[name]])
    cat("✓\n")
  } else {
    cat("✗ NOT FOUND\n")
  }
}

cat("\n")

# Quick comparison ----
if (length(data_list) > 0) {
  cat("Quick Comparison:\n")
  cat(sprintf("%-10s %12s %10s %15s %15s\n",
    "File", "Rows", "Cols", "Min Pos", "Max Pos"
  ))
  cat(paste(rep("-", 65), collapse = ""), "\n")

  for (name in names(data_list)) {
    d <- data_list[[name]]
    cat(sprintf(
      "%-10s %12s %10d %15s %15s\n",
      name,
      format(nrow(d), big.mark = ","),
      ncol(d),
      format(min(d$position, na.rm = TRUE), big.mark = ","),
      format(max(d$position, na.rm = TRUE), big.mark = ",")
    ))
  }

  cat("\n")

  # Check position in 1M-2M range
  cat("Positions in 1M-2M range:\n")
  for (name in names(data_list)) {
    d <- data_list[[name]]
    count <- sum(d$position >= 1e6 & d$position <= 2e6, na.rm = TRUE)
    cat(sprintf("  %s: %s rows\n", name, format(count, big.mark = ",")))
  }

  # Column comparison
  cat("\nColumn names:\n")
  all_cols <- lapply(data_list, names)

  # Check if all have same columns
  first_cols <- all_cols[[1]]
  same_cols <- all(sapply(all_cols, function(x) identical(x, first_cols)))

  if (same_cols) {
    cat("  ✓ All files have identical column names\n")
    cat(sprintf("  Columns (%d): %s\n",
      length(first_cols),
      paste(head(first_cols, 10), collapse = ", ")
    ))
    if (length(first_cols) > 10) {
      cat(sprintf("    ... and %d more\n", length(first_cols) - 10))
    }
  } else {
    cat("  ✗ Files have different columns:\n")
    for (name in names(all_cols)) {
      cat(sprintf("    %s: %s\n",
        name,
        paste(head(all_cols[[name]], 5), collapse = ", ")
      ))
    }
  }
}

cat("\n")
cat("For detailed analysis, run:\n")
cat("  source('scripts/diagnose_data_files.r')\n")
