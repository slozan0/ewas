# Diagnostic Script for Data Files
# Checks for missing positions, duplicates, and compares file structures

rm(list = ls())

library(data.table)

# Configuration ----
input_dir <- "data/output"

# Files to diagnose
files_to_check <- c(
  "5feb_tem_a1_c3_rc.rds",
  "5feb_tem_a2_c3_rc.rds",
  "5feb_tem_d1_c3_rc.rds",
  "5feb_tem_d2_c3_rc.rds"
)

position_col <- "position"  # Adjust if different

# Function to analyze a single file ----
analyze_file <- function(file_path, pos_col = "position") {
  cat("\n")
  cat("========================================\n")
  cat(sprintf("Analyzing: %s\n", basename(file_path)))
  cat("========================================\n")

  if (!file.exists(file_path)) {
    cat("ERROR: File not found!\n")
    return(NULL)
  }

  # Load data
  cat("Loading data...\n")
  data <- readRDS(file_path)

  if (!is.data.table(data)) {
    setDT(data)
  }

  # Basic info
  cat(sprintf("  Rows: %s\n", format(nrow(data), big.mark = ",")))
  cat(sprintf("  Columns: %d\n", ncol(data)))
  cat(sprintf("  Column names: %s\n",
    paste(head(names(data), 10), collapse = ", ")
  ))
  if (ncol(data) > 10) {
    cat(sprintf("    ... and %d more\n", ncol(data) - 10))
  }

  # Check if position column exists
  if (!pos_col %in% names(data)) {
    cat(sprintf("\nERROR: Position column '%s' not found!\n", pos_col))
    cat("Available columns:\n")
    cat(paste("  -", names(data), collapse = "\n"))
    cat("\n")
    return(NULL)
  }

  # Position statistics
  cat("\n--- Position Statistics ---\n")
  positions <- data[[pos_col]]

  cat(sprintf("  Min position: %s\n",
    format(min(positions, na.rm = TRUE), big.mark = ",")
  ))
  cat(sprintf("  Max position: %s\n",
    format(max(positions, na.rm = TRUE), big.mark = ",")
  ))
  cat(sprintf("  Range: %s\n",
    format(
      max(positions, na.rm = TRUE) - min(positions, na.rm = TRUE),
      big.mark = ","
    )
  ))

  # Check for NAs
  n_na <- sum(is.na(positions))
  if (n_na > 0) {
    cat(sprintf("  WARNING: %s NA positions found!\n",
      format(n_na, big.mark = ",")
    ))
  } else {
    cat("  No NA positions\n")
  }

  # Check for duplicates
  cat("\n--- Duplicate Check ---\n")
  n_duplicates <- sum(duplicated(positions))
  if (n_duplicates > 0) {
    cat(sprintf("  WARNING: %s duplicate positions!\n",
      format(n_duplicates, big.mark = ",")
    ))

    # Show some examples
    dup_positions <- positions[duplicated(positions)]
    cat("  Example duplicates (first 10):\n")
    cat(paste("   ",
      format(head(unique(dup_positions), 10), big.mark = ","),
      collapse = "\n"
    ))
    cat("\n")
  } else {
    cat("  No duplicate positions\n")
  }

  # Check for gaps
  cat("\n--- Gap Analysis ---\n")
  sorted_pos <- sort(unique(positions))
  gaps <- diff(sorted_pos)
  cat(sprintf("  Unique positions: %s\n",
    format(length(sorted_pos), big.mark = ",")
  ))
  cat(sprintf("  Min gap: %s\n", format(min(gaps), big.mark = ",")))
  cat(sprintf("  Max gap: %s\n", format(max(gaps), big.mark = ",")))
  cat(sprintf("  Median gap: %s\n", format(median(gaps), big.mark = ",")))

  # Find large gaps (> 10000)
  large_gaps <- which(gaps > 10000)
  if (length(large_gaps) > 0) {
    cat(sprintf("  Large gaps (>10,000): %d\n", length(large_gaps)))
    cat("  First 5 large gaps:\n")
    for (i in head(large_gaps, 5)) {
      cat(sprintf(
        "    Position %s to %s (gap: %s)\n",
        format(sorted_pos[i], big.mark = ","),
        format(sorted_pos[i + 1], big.mark = ","),
        format(gaps[i], big.mark = ",")
      ))
    }
  }

  # Position density in 1M bins
  cat("\n--- Position Distribution (1M bins) ---\n")
  min_pos <- min(positions, na.rm = TRUE)
  max_pos <- max(positions, na.rm = TRUE)

  bins <- seq(
    floor(min_pos / 1e6) * 1e6,
    ceiling(max_pos / 1e6) * 1e6,
    by = 1e6
  )

  if (length(bins) <= 20) {
    for (i in 1:(length(bins) - 1)) {
      count <- sum(positions >= bins[i] & positions < bins[i + 1])
      if (count > 0) {
        cat(sprintf(
          "  %sM - %sM: %s positions\n",
          format(bins[i] / 1e6, big.mark = ","),
          format(bins[i + 1] / 1e6, big.mark = ","),
          format(count, big.mark = ",")
        ))
      }
    }
  } else {
    cat("  Too many bins to display (>20)\n")
    cat(sprintf(
      "  Try: range %s to %s\n",
      format(min_pos, big.mark = ","),
      format(max_pos, big.mark = ",")
    ))
  }

  # Return summary for comparison
  summary <- list(
    file = basename(file_path),
    n_rows = nrow(data),
    n_cols = ncol(data),
    col_names = names(data),
    min_pos = min(positions, na.rm = TRUE),
    max_pos = max(positions, na.rm = TRUE),
    n_unique_pos = length(unique(positions)),
    n_duplicates = n_duplicates,
    n_na = n_na
  )

  return(summary)
}

# Analyze all files ----
cat("\n")
cat("████████████████████████████████████████\n")
cat("  DATA FILE DIAGNOSTIC REPORT\n")
cat("████████████████████████████████████████\n")

summaries <- list()

for (file_name in files_to_check) {
  file_path <- file.path(input_dir, file_name)
  summary <- analyze_file(file_path, pos_col = position_col)
  if (!is.null(summary)) {
    summaries[[file_name]] <- summary
  }
}

# Comparison across files ----
if (length(summaries) > 1) {
  cat("\n")
  cat("========================================\n")
  cat("CROSS-FILE COMPARISON\n")
  cat("========================================\n")

  # Compare column names
  cat("\n--- Column Name Comparison ---\n")
  all_cols <- unique(unlist(lapply(summaries, function(x) x$col_names)))

  cat("Common columns across all files:\n")
  common_cols <- Reduce(
    intersect,
    lapply(summaries, function(x) x$col_names)
  )
  cat(paste("  -", common_cols, collapse = "\n"))
  cat("\n")

  if (length(common_cols) < length(all_cols)) {
    cat("\nColumns not in all files:\n")
    for (file_name in names(summaries)) {
      unique_cols <- setdiff(summaries[[file_name]]$col_names, common_cols)
      if (length(unique_cols) > 0) {
        cat(sprintf("  %s: %s\n",
          file_name,
          paste(unique_cols, collapse = ", ")
        ))
      }
    }
  }

  # Compare position ranges
  cat("\n--- Position Range Comparison ---\n")
  cat(sprintf(
    "%-35s %15s %15s %10s\n",
    "File", "Min Pos", "Max Pos", "Unique"
  ))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (file_name in names(summaries)) {
    s <- summaries[[file_name]]
    cat(sprintf(
      "%-35s %15s %15s %10s\n",
      substr(file_name, 1, 35),
      format(s$min_pos, big.mark = ","),
      format(s$max_pos, big.mark = ","),
      format(s$n_unique_pos, big.mark = ",")
    ))
  }

  # Check for position overlaps
  cat("\n--- Position Overlap Check ---\n")
  if (length(summaries) >= 2) {
    file_names <- names(summaries)
    for (i in 1:(length(file_names) - 1)) {
      for (j in (i + 1):length(file_names)) {
        f1 <- file_names[i]
        f2 <- file_names[j]
        s1 <- summaries[[f1]]
        s2 <- summaries[[f2]]

        # Check if ranges overlap
        overlap <- !(s1$max_pos < s2$min_pos || s2$max_pos < s1$min_pos)

        if (overlap) {
          overlap_start <- max(s1$min_pos, s2$min_pos)
          overlap_end <- min(s1$max_pos, s2$max_pos)
          cat(sprintf(
            "  %s & %s: OVERLAP %s to %s\n",
            substr(f1, 1, 20),
            substr(f2, 1, 20),
            format(overlap_start, big.mark = ","),
            format(overlap_end, big.mark = ",")
          ))
        } else {
          cat(sprintf(
            "  %s & %s: No overlap\n",
            substr(f1, 1, 20),
            substr(f2, 1, 20)
          ))
        }
      }
    }
  }

  # Summary table
  cat("\n--- Summary Statistics ---\n")
  cat(sprintf(
    "%-35s %10s %10s %10s\n",
    "File", "Rows", "Dups", "NAs"
  ))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (file_name in names(summaries)) {
    s <- summaries[[file_name]]
    cat(sprintf(
      "%-35s %10s %10s %10s\n",
      substr(file_name, 1, 35),
      format(s$n_rows, big.mark = ","),
      format(s$n_duplicates, big.mark = ","),
      format(s$n_na, big.mark = ",")
    ))
  }
}

# Save report
cat("\n")
cat("========================================\n")
cat("Saving detailed report...\n")
report_file <- "data/diagnostic_report.txt"

# Redirect output to file
sink(report_file)
cat("DATA FILE DIAGNOSTIC REPORT\n")
cat("Generated:", format(Sys.time()), "\n\n")

for (file_name in files_to_check) {
  file_path <- file.path(input_dir, file_name)
  analyze_file(file_path, pos_col = position_col)
}

sink()

cat(sprintf("Report saved to: %s\n", report_file))
cat("========================================\n\n")

cat("✓ Diagnostic complete!\n")
