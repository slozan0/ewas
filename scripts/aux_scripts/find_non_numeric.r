# Find Non-Numeric Values in Position Column
# This helps diagnose why a column is character instead of integer

rm(list = ls())

# Load the problematic file
cat("Loading a2...\n")
data_a2 <- readRDS("data/output/5feb_tem_a2_c3_rc.rds")

cat(sprintf("Total rows: %s\n\n", format(nrow(data_a2), big.mark = ",")))

# Check data type
cat("Position column type:", class(data_a2$position), "\n\n")

# Find non-numeric values
cat("Finding non-numeric positions...\n")

# Try to convert to numeric - NAs will be created for non-numeric values
numeric_positions <- suppressWarnings(as.numeric(data_a2$position))

# Find which ones became NA (these were non-numeric)
non_numeric_idx <- which(is.na(numeric_positions) & !is.na(data_a2$position))

cat(sprintf("Non-numeric positions found: %d\n",
  length(non_numeric_idx)
))

if (length(non_numeric_idx) > 0) {
  cat("\nFirst 20 non-numeric position values:\n")
  cat("----------------------------------------\n")

  # Show examples
  examples <- head(data_a2$position[non_numeric_idx], 20)
  for (i in seq_along(examples)) {
    cat(sprintf(
      "Row %s: '%s'\n",
      format(non_numeric_idx[i], big.mark = ","),
      examples[i]
    ))
  }

  # Check if they're empty strings
  cat("\n--- Analysis ---\n")
  empty_strings <- sum(data_a2$position[non_numeric_idx] == "")
  if (empty_strings > 0) {
    cat(sprintf("Empty strings: %d\n", empty_strings))
  }

  # Check for whitespace
  whitespace_only <- sum(grepl("^\\s+$", data_a2$position[non_numeric_idx]))
  if (whitespace_only > 0) {
    cat(sprintf("Whitespace only: %d\n", whitespace_only))
  }

  # Check for special characters
  special_chars <- data_a2$position[non_numeric_idx]
  special_chars <- special_chars[!special_chars %in% c("", " ")]
  if (length(special_chars) > 0) {
    cat("Special/unusual values:\n")
    cat(paste("  ", head(unique(special_chars), 10)), sep = "\n")
  }

  # Summary by pattern
  cat("\n--- Pattern Summary ---\n")
  cat("Unique non-numeric values:\n")
  unique_bad <- unique(data_a2$position[non_numeric_idx])
  cat(sprintf("  Total unique: %d\n", length(unique_bad)))

  if (length(unique_bad) <= 20) {
    cat("  All unique values:\n")
    for (val in unique_bad) {
      count <- sum(data_a2$position[non_numeric_idx] == val)
      cat(sprintf("    '%s' (occurs %d times)\n", val, count))
    }
  } else {
    cat("  Top 20 by frequency:\n")
    freq_table <- sort(table(data_a2$position[non_numeric_idx]),
                       decreasing = TRUE)
    print(head(freq_table, 20))
  }

  # Show full rows with problematic positions
  cat("\n--- Sample Problematic Rows ---\n")
  problem_rows <- head(data_a2[non_numeric_idx, ], 5)
  print(problem_rows)

} else {
  cat("✓ All positions can be converted to numeric!\n")
  cat("The issue might be elsewhere.\n")
}

# Also check for decimal numbers (would be numeric but not integer)
cat("\n--- Checking for Decimal Numbers ---\n")
numeric_positions_clean <- numeric_positions[!is.na(numeric_positions)]
has_decimals <- any(numeric_positions_clean != floor(numeric_positions_clean))

if (has_decimals) {
  decimal_idx <- which(numeric_positions != floor(numeric_positions))
  cat(sprintf("Decimal positions found: %d\n", length(decimal_idx)))
  cat("Examples:\n")
  print(head(data_a2$position[decimal_idx], 10))
} else {
  cat("No decimal numbers found\n")
}

cat("\n========================================\n")
cat("Recommendation:\n")
cat("========================================\n")

if (length(non_numeric_idx) > 0) {
  pct <- (length(non_numeric_idx) / nrow(data_a2)) * 100

  cat(sprintf("%.4f%% of positions are non-numeric\n\n", pct))

  if (pct < 0.01) {
    cat("This is a small fraction. You can:\n")
    cat("1. Remove these rows:\n")
    cat("   data_clean <- data_a2[!is.na(as.numeric(data_a2$position)), ]\n\n")
    cat("2. Or fix them if they follow a pattern\n")
  } else {
    cat("This is significant. Investigate why.\n")
    cat("Possible causes:\n")
    cat("- Data corruption\n")
    cat("- Mixed data sources\n")
    cat("- Processing error upstream\n")
  }
} else {
  cat("All positions are numeric.\n")
  cat("The character type might be from:\n")
  cat("- read.table() with stringsAsFactors=FALSE\n")
  cat("- Explicit conversion somewhere\n")
  cat("\nYou can safely convert:\n")
  cat("  data_a2$position <- as.integer(data_a2$position)\n")
}
