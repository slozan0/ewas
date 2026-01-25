# Automated Testing Guide for Beginners

If you've never done automated testing before, don't worry! This guide explains everything.

## What is Automated Testing?

**The Problem**: You make changes to your code and want to make sure nothing broke.

**Manual Testing** (old way):
1. Run the script
2. Look at the output
3. Check if it "looks right"
4. Hope you didn't miss anything

**Automated Testing** (better way):
1. Write tests once
2. Run `source('run_tests.R')` anytime
3. Tests tell you if something broke
4. Takes seconds, catches bugs you'd miss

---

## Key Concepts

### 1. **Test Files**

Tests live in `tests/testthat/test-*.R`:
- `test-physmap.R` - Tests for physmap.r
- `test-get_hetero.R` - Tests for get_hetero function
- Each file contains multiple test cases

### 2. **Golden Files** (Expected Outputs)

Golden files are **the correct answer** that your code should produce.

**Location**: `tests/testthat/fixtures/`

**How they work**:
```
First Run:
  Your code → Output → Saved as golden file ✓

Future Runs:
  Your code → New output → Compare to golden file

  If identical: ✓ Test passes
  If different: ✗ Test fails (you broke something!)
```

**Example**:
```r
# Run your code
result <- my_function(input_data)

# First time: create golden file
if (!file.exists("golden.rds")) {
  saveRDS(result, "golden.rds")
} else {
  # Compare to golden
  golden <- readRDS("golden.rds")
  expect_equal(result, golden)  # Test passes if identical
}
```

### 3. **Test Assertions**

Assertions check if something is true:

```r
# Check if equal
expect_equal(actual, expected)

# Check if TRUE
expect_true(x > 0)

# Check column names
expect_named(data, c("a", "b", "c"))

# Check no errors
expect_error(my_function(), NA)
```

---

## Testing physmap.r - Step by Step

### **Step 1: Run the Test for the First Time**

This creates the golden file:

```r
# Install testthat if needed
renv::install("testthat")

# Run physmap test
library(testthat)
test_file("tests/testthat/test-physmap.R")
```

**Expected output**:
```
✔ | 7 | test-physmap
Skip (test-physmap.R:115:5): Created golden output file
```

**What happened**:
- Test ran successfully
- No golden file existed
- Test created `tests/testthat/fixtures/physmap_golden_output.rds`
- This is now **your reference** for future runs

---

### **Step 2: Verify the Golden File is Correct**

**IMPORTANT**: The golden file should be from working, correct code!

```r
# Load and inspect the golden file
golden <- readRDS("tests/testthat/fixtures/physmap_golden_output.rds")

# Check structure
str(golden)

# Check sample rows
head(golden)

# Check if it looks correct
summary(golden)
```

**If it looks good**: You're done! This is your baseline.

**If it looks wrong**: Delete it and fix your code first!
```r
# Delete bad golden file
unlink("tests/testthat/fixtures/physmap_golden_output.rds")

# Fix your code
# Run test again to create new golden file
```

---

### **Step 3: Run Tests After Making Changes**

Anytime you modify code:

```r
# Run all tests
source("run_tests.R")

# Or just physmap tests
library(testthat)
test_file("tests/testthat/test-physmap.R")
```

**If tests pass** (green ✓):
```
✔ | 8 | test-physmap
```
Your changes didn't break anything!

**If tests fail** (red ✗):
```
✖ | 1 8 | test-physmap
Failure (test-physmap.R:120:5): Output should match golden
nrow(dt_chrom) not equal to nrow(golden_output)
1234 != 1189
```
Something changed! Either:
- You broke something (fix it)
- The change is intentional (update golden file)

---

## Where Golden Files Go

```
tests/
└── testthat/
    ├── fixtures/                    ← Golden files here!
    │   ├── physmap_golden_output.rds
    │   ├── ezchi_golden_output.rds
    │   └── tiny_sample.rds          ← Small test data
    ├── test-physmap.R               ← Test code
    ├── test-get_hetero.R
    └── ...
```

**Why `fixtures/`?**
- Standard naming convention
- Clearly separates test code from test data
- Easy to find and manage

---

## Common Testing Workflows

### **Workflow 1: Adding a New Feature**

```r
# 1. Write the feature code
# 2. Run tests to make sure you didn't break existing features
source("run_tests.R")

# 3. Write tests for the new feature
# (add to test-*.R file)

# 4. Run tests again
source("run_tests.R")
```

### **Workflow 2: Fixing a Bug**

```r
# 1. Write a test that FAILS (reproduces the bug)
test_that("bug is fixed", {
  result <- buggy_function(input)
  expect_equal(result, correct_answer)
})

# 2. Run test - should fail
test_file("tests/testthat/test-mybug.R")

# 3. Fix the bug

# 4. Run test - should pass now!
test_file("tests/testthat/test-mybug.R")
```

### **Workflow 3: Intentional Change to Output**

If you intentionally change how code works:

```r
# 1. Make your code changes

# 2. Run tests - they will fail
source("run_tests.R")
# ✗ Output doesn't match golden

# 3. Verify new output is correct
result <- my_function(input)
head(result)  # Looks good!

# 4. Update golden file
saveRDS(result, "tests/testthat/fixtures/my_golden.rds")

# 5. Run tests again - should pass
source("run_tests.R")
# ✓ All tests pass
```

---

## Example: Complete Test for physmap

Here's what the test does:

```r
test_that("physmap processes sample data correctly", {

  # === SETUP ===
  # Load sample data
  raw_data <- readRDS("data/sample/5feb_tem_a1_c3_sample.rds")

  # === RUN THE CODE ===
  # (Same logic as physmap.r)
  # ... processing steps ...
  result <- dt_chrom

  # === CHECK RESULTS ===

  # Test 1: Structure is correct
  expect_s3_class(result, "data.table")
  expect_named(result, c("chrom", "pos", "ref", ...))

  # Test 2: Values are in valid range
  expect_true(all(result$sumDepth >= 25))
  expect_true(all(result$sumDepth < 1000))

  # Test 3: No invalid data
  expect_false(any(is.na(result$pos)))
  expect_true(all(result$ref %in% c("A", "C", "G", "T")))

  # Test 4: Compare to golden file
  golden <- readRDS("tests/testthat/fixtures/physmap_golden.rds")
  expect_equal(result, golden)
})
```

---

## Interpreting Test Results

### ✅ **All Tests Pass**
```
✔ | 8 | test-physmap
✔ | 5 | test-get_hetero
✔ | 3 | test-get_chi

========================================
Test Summary
========================================
Passed: 16
Failed: 0
Warnings: 0
Skipped: 0
========================================
```
**Meaning**: Everything works! Safe to commit.

---

### ⚠️ **Test Skipped**
```
Skip (test-physmap.R:10:3): Sample data not found
```
**Meaning**: Test couldn't run (missing file, etc.). Not a failure, but test didn't check anything.

**Fix**: Create the missing file, then run again.

---

### ❌ **Test Failed**
```
✖ | 1 7 | test-physmap

Failure (test-physmap.R:95:3): sumDepth filtering worked
all(dt_chrom$sumDepth >= 25) is not TRUE

Actual value: FALSE
```

**Meaning**: Code produced wrong output.

**How to debug**:
```r
# Run the test code manually
source("tests/testthat/test-physmap.R")

# Or run interactively
library(testthat)
test_file("tests/testthat/test-physmap.R")

# Check what went wrong
dt_chrom[sumDepth < 25]  # These shouldn't exist!
```

---

## Managing Golden Files

### **When to Update Golden Files**

✅ **Update when**:
- You intentionally changed output format
- You fixed a bug (old output was wrong)
- You improved an algorithm (new output is better)

❌ **Don't update when**:
- Test fails and you don't know why
- You're not sure if new output is correct
- You're just trying to make tests pass

### **How to Update**

```r
# Option 1: Delete and regenerate
unlink("tests/testthat/fixtures/physmap_golden_output.rds")
test_file("tests/testthat/test-physmap.R")  # Creates new one

# Option 2: Manually save new output
new_result <- my_corrected_function(input)
saveRDS(new_result, "tests/testthat/fixtures/my_golden.rds")
```

### **Golden File Best Practices**

1. **Keep them small** - Use sample data, not full datasets
2. **Use compression** - `saveRDS(data, file, compress = "xz")`
3. **Document them** - Add comments explaining what they are
4. **Version control** - Commit golden files to git
5. **Review changes** - When updating, check the diff

---

## Tips for Beginners

### **Start Small**
Don't try to test everything at once. Start with:
1. One simple test
2. Get it working
3. Add more tests gradually

### **Use Descriptive Test Names**
```r
# Good
test_that("physmap filters out low-coverage sites", { ... })

# Bad
test_that("test1", { ... })
```

### **One Test, One Concept**
```r
# Good - separate tests
test_that("sumDepth is calculated correctly", { ... })
test_that("positions are filtered by depth", { ... })

# Bad - too much in one test
test_that("everything works", { ... })
```

### **Use Helpful Error Messages**
```r
# Good - tells you what went wrong
expect_true(all(x > 0), info = "All values should be positive")

# Less helpful
expect_true(all(x > 0))
```

---

## Next Steps

1. **Create golden files**:
   ```r
   test_file("tests/testthat/test-physmap.R")
   ```

2. **Inspect golden files**:
   ```r
   golden <- readRDS("tests/testthat/fixtures/physmap_golden_output.rds")
   str(golden)
   summary(golden)
   ```

3. **Run all tests**:
   ```r
   source("run_tests.R")
   ```

4. **Make changes, run tests again**:
   ```r
   # Edit code
   source("run_tests.R")  # Did you break anything?
   ```

5. **Commit tests and golden files**:
   ```bash
   git add tests/
   git commit -m "Add automated tests for physmap"
   ```

---

## Questions?

**Q: What if my golden file is huge?**
A: Use smaller sample data. The golden file should be based on `data/sample/*`, not full datasets.

**Q: Tests pass on my machine but fail on colleague's machine?**
A: Probably different package versions or random number generation. Use `set.seed()` and ensure same package versions via `renv`.

**Q: How do I test functions that use randomness?**
A: Use `set.seed()` before calling the function to make results reproducible.

**Q: Should I commit golden files to git?**
A: Yes! They're part of your test suite.

**Q: How many tests do I need?**
A: Enough to feel confident. Start with main functionality, add edge cases later.

---

Happy testing! 🧪
