# Creating Sample Data

## Quick Start

```r
# Run from project root
source("scripts/create_sample_data.r")
```

This will create sample files from positions **1,000,000 to 2,000,000** on chromosome 3.

---

## What Gets Created

The script processes these 4 files:
- `data/output/5feb_tem_a1_c3_rc.rds` → `data/sample/5feb_tem_a1_c3_sample.rds`
- `data/output/5feb_tem_a2_c3_rc.rds` → `data/sample/5feb_tem_a2_c3_sample.rds`
- `data/output/5feb_tem_d1_c3_rc.rds` → `data/sample/5feb_tem_d1_c3_sample.rds`
- `data/output/5feb_tem_d2_c3_rc.rds` → `data/sample/5feb_tem_d2_c3_sample.rds`

---

## Customizing the Sample

### Change Position Range

Edit `scripts/create_sample_data.r`:

```r
# Configuration ----
pos_start <- 1000000  # Change this
pos_end <- 2000000    # Change this
```

**Why 1M-2M?** There isn't much data before position 1,000,000 in your dataset.

### Change Position Column Name

If your data uses a different column name:

```r
position_col <- "position"  # Change to "pos", "Position", etc.
```

### Add Different Files

Add more files to the `sample_files` list:

```r
sample_files <- list(
  # ... existing files ...
  list(
    input = file.path(input_dir, "your_file.rds"),
    output = file.path(sample_dir, "your_sample.rds")
  )
)
```

---

## Verifying Sample Data

### Check What Was Created

```r
# List sample files
list.files("data/sample", pattern = "\\.rds$")

# Check sizes
sapply(
  list.files("data/sample", pattern = "\\.rds$", full.names = TRUE),
  function(f) file.info(f)$size / 1024^2  # Size in MB
)
```

### Inspect Sample Content

```r
# Load a sample
sample <- readRDS("data/sample/5feb_tem_a1_c3_sample.rds")

# Check structure
str(sample)

# Check position range
range(sample$position)  # Should be ~1M to 2M

# Number of rows
nrow(sample)
```

---

## File Size Guidelines

**Target**: Keep each file under 10 MB

**If files are too large:**

1. **Narrow the position range:**
   ```r
   pos_start <- 1000000
   pos_end <- 1500000  # Reduced from 2M to 1.5M
   ```

2. **Sample fewer rows:**
   ```r
   # After filtering by position, sample randomly
   sample_data <- sample_data[sample(.N, min(.N, 10000))]
   ```

3. **Use higher compression:**
   ```r
   saveRDS(sample_data, output_file, compress = "xz")
   ```

---

## Troubleshooting

### "Position column not found"

Check your column names:
```r
data <- readRDS("data/output/5feb_tem_a1_c3_rc.rds")
names(data)  # List all column names
```

Then update `position_col` in the script.

### "No data in specified position range"

Your data might be in a different range:
```r
data <- readRDS("data/output/5feb_tem_a1_c3_rc.rds")
range(data$position, na.rm = TRUE)  # Check actual range
```

Adjust `pos_start` and `pos_end` accordingly.

### "Input file not found"

Check the file path:
```r
# Verify file exists
file.exists("data/output/5feb_tem_a1_c3_rc.rds")

# List available files
list.files("data/output", pattern = "c3.*\\.rds$")
```

---

## After Creating Samples

1. **Review files:**
   ```bash
   ls -lh data/sample/*.rds
   ```

2. **Test with sample data:**
   ```r
   source("run_tests.R")
   ```

3. **Commit to git:**
   ```bash
   git add data/sample/
   git commit -m "Add sample data (chr3, positions 1M-2M)"
   ```

---

## Why This Sampling Strategy?

✅ **Position-based sampling** is better than just taking first N rows because:
- Captures real biological variation
- More representative of actual data
- Better tests edge cases in specific genomic regions
- Ensures data quality in your region of interest

✅ **Position 1M-2M** was chosen because:
- Your data has sparse coverage before 1M
- This range has good data density
- Still small enough for fast testing
- Large enough to be meaningful

---

## Advanced: Creating Test Fixtures

For unit tests, you might want even smaller samples:

```r
# After running create_sample_data.r
sample <- readRDS("data/sample/5feb_tem_a1_c3_sample.rds")

# Create tiny test fixture (100 rows)
test_fixture <- sample[1:100, ]
saveRDS(test_fixture, "tests/testthat/fixtures/tiny_sample.rds")
```
