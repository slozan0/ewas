#!/usr/bin/env Rscript
# Run all tests for the EWAS PhysMap project

library(testthat)

# Run tests
cat("Running tests...\n\n")
test_results <- test_dir(
  "tests/testthat",
  reporter = "progress"  # Can also use "summary", "minimal", or "check"
)

# Print summary
cat("\n")
cat("========================================\n")
cat("Test Summary\n")
cat("========================================\n")
cat(sprintf("Passed: %d\n", sum(test_results$passed)))
cat(sprintf("Failed: %d\n", sum(test_results$failed)))
cat(sprintf("Warnings: %d\n", sum(test_results$warning)))
cat(sprintf("Skipped: %d\n", sum(test_results$skipped)))
cat("========================================\n")

# Exit with error code if any tests failed
if (sum(test_results$failed) > 0) {
  quit(status = 1)
}
