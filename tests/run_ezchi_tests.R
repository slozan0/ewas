# Run Easy Chi2 Tests
# Quick script to run only the easy_chi2 related tests

library(testthat)
library(here)

cat("\n")
cat("========================================\n")
cat("Running Easy Chi2 Tests\n")
cat("========================================\n\n")

# Test files to run
test_files <- c(
  "test-get_hetero.R",
  "test-get_chi.R",
  "test-get_alternate_alleles.R",
  "test-regression.R",
  "test-easy_chi2_pipeline.R"
)

total_passed <- 0
total_failed <- 0
total_skipped <- 0

for (test_file in test_files) {
  test_path <- here("tests", "testthat", test_file)

  if (file.exists(test_path)) {
    cat(sprintf("\n--- Running %s ---\n", test_file))

    result <- test_file(test_path, reporter = "progress")

    total_passed <- total_passed + sum(result$passed)
    total_failed <- total_failed + sum(result$failed)
    total_skipped <- total_skipped + sum(result$skipped)
  } else {
    cat(sprintf("⚠ Skipping %s (not found)\n", test_file))
  }
}

cat("\n")
cat("========================================\n")
cat("Easy Chi2 Test Summary\n")
cat("========================================\n")
cat(sprintf("✓ Passed:  %d\n", total_passed))
cat(sprintf("✗ Failed:  %d\n", total_failed))
cat(sprintf("⊘ Skipped: %d\n", total_skipped))
cat("========================================\n\n")

if (total_failed > 0) {
  cat("❌ Some tests failed! Review the output above.\n\n")
  quit(status = 1)
} else {
  cat("✅ All easy_chi2 tests passed!\n\n")
}
