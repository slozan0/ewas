# Testing Quick Start

New to automated testing? Follow these steps!

## Setup (One-Time)

### 1. Install testthat
```r
renv::install("testthat")
renv::snapshot()
```

### 2. Create sample data (if not already done)
```r
source("scripts/create_sample_data.r")
```

### 3. Create golden files (expected outputs)
```r
source("tests/create_golden_files.R")
# Choose option 1: Create physmap golden file
```

---

## Daily Workflow

### Run all tests
```r
source("run_tests.R")
```

**Expected output:**
```
✔ | 8 | test-physmap
✔ | 5 | test-get_hetero
✔ | 3 | test-get_chi

========================================
Test Summary
========================================
Passed: 16
Failed: 0
========================================
```

✅ **Green = Good!** Your code works.
❌ **Red = Problem!** You broke something or need to update golden files.

---

## Files You'll Work With

```
tests/
├── README.md                    ← You are here
├── TESTING_GUIDE.md             ← Full guide (read this!)
├── create_golden_files.R        ← Helper to create expected outputs
├── run_tests.R                  ← Run all tests
├── testthat.R                   ← Auto-runs tests
└── testthat/
    ├── fixtures/                ← Expected outputs ("golden files")
    │   └── physmap_golden_output.rds
    ├── test-physmap.R           ← Tests for physmap.r
    ├── test-get_hetero.R        ← Tests for get_hetero()
    └── test-get_chi.R           ← Tests for get_chi()
```

---

## Common Tasks

### ✅ Test a single file
```r
library(testthat)
test_file("tests/testthat/test-physmap.R")
```

### ✅ Create/update golden files
```r
source("tests/create_golden_files.R")
```

### ✅ List golden files
```r
list.files("tests/testthat/fixtures")
```

### ✅ Inspect a golden file
```r
golden <- readRDS("tests/testthat/fixtures/physmap_golden_output.rds")
str(golden)
head(golden)
```

### ✅ Delete and recreate golden files
```r
# Delete
unlink("tests/testthat/fixtures/physmap_golden_output.rds")

# Recreate
source("tests/create_golden_files.R")
```

---

## When Tests Fail

### Option 1: You broke something (most common)
```r
# Fix your code
# Run tests again
source("run_tests.R")
```

### Option 2: Output intentionally changed
```r
# 1. Verify new output is correct
result <- my_function(input)
head(result)  # Looks good!

# 2. Update golden file
source("tests/create_golden_files.R")

# 3. Run tests again
source("run_tests.R")  # Should pass now
```

---

## Learn More

**New to testing?** Read [TESTING_GUIDE.md](TESTING_GUIDE.md) - it explains everything!

**Key concepts:**
- **Test files** (`test-*.R`) - Code that checks if your functions work
- **Golden files** (in `fixtures/`) - Expected correct outputs
- **Assertions** (`expect_equal()`) - Checks that verify results

---

## Quick Reference

| Task | Command |
|------|---------|
| Run all tests | `source("run_tests.R")` |
| Run one test file | `test_file("tests/testthat/test-physmap.R")` |
| Create golden files | `source("tests/create_golden_files.R")` |
| List golden files | `list.files("tests/testthat/fixtures")` |
| View golden file | `readRDS("tests/testthat/fixtures/physmap_golden_output.rds")` |
| Delete golden file | `unlink("tests/testthat/fixtures/physmap_golden_output.rds")` |

---

## Workflow Example

```r
# 1. Make changes to physmap.r
# ... edit code ...

# 2. Run tests
source("run_tests.R")

# 3a. If tests pass - you're done! ✓
# 3b. If tests fail:

# Check what failed
test_file("tests/testthat/test-physmap.R")

# If you broke something: fix and retest
# If output intentionally changed: update golden file
source("tests/create_golden_files.R")

# 4. Commit your changes
# git add .
# git commit -m "Update physmap with new filtering logic"
```

---

## Getting Help

1. **Read the guide**: [TESTING_GUIDE.md](TESTING_GUIDE.md)
2. **Look at test examples**: Open `test-physmap.R` to see how tests work
3. **Check test output**: It tells you exactly what failed
4. **Run tests interactively**: Copy test code into R console to debug

---

## Pro Tips

✅ **Run tests before committing** - Catch bugs early
✅ **Keep golden files small** - Use sample data, not full datasets
✅ **Update golden files carefully** - Make sure new output is actually correct
✅ **Write descriptive test names** - Future you will thank you
✅ **One test per concept** - Easier to understand failures

---

Happy testing! 🧪

Questions? Check [TESTING_GUIDE.md](TESTING_GUIDE.md) for detailed explanations.
