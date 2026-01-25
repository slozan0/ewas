# Data Diagnostics Scripts

Two scripts to inspect and compare your chromosome 3 data files.

## Quick Check (Fast)

For a fast overview:

```r
source("scripts/quick_check_files.r")
```

**What it shows:**
- File sizes (rows/columns)
- Position ranges (min/max)
- How many positions in 1M-2M range
- Whether columns match across files

**Output example:**
```
File         Rows       Cols         Min Pos         Max Pos
-----------------------------------------------------------------
a1         41,971      50         547,912      41,963,435
a2         42,105      50         548,001      42,104,789
...

Positions in 1M-2M range:
  a1: 1,234 rows
  a2: 1,189 rows
...
```

---

## Full Diagnostic (Detailed)

For comprehensive analysis:

```r
source("scripts/diagnose_data_files.r")
```

**What it checks:**

### Per File:
- ✅ Basic stats (rows, columns, column names)
- ✅ Position range (min, max, median gap)
- ✅ **Duplicate positions** ⚠️
- ✅ **Missing values (NAs)** ⚠️
- ✅ **Large gaps** in positions
- ✅ Position distribution by 1M bins

### Cross-File:
- ✅ Column name comparison
- ✅ Position range overlaps
- ✅ Summary statistics table

**Output example:**
```
Analyzing: 5feb_tem_a1_c3_rc.rds
========================================
  Rows: 41,971,435
  Columns: 50
  ...

--- Duplicate Check ---
  WARNING: 1,234 duplicate positions!
  Example duplicates:
    1,234,567
    2,345,678
    ...

--- Gap Analysis ---
  Min gap: 1
  Max gap: 45,123
  Median gap: 15
  Large gaps (>10,000): 5
    Position 1,000,000 to 1,050,000 (gap: 50,000)
    ...

--- Position Distribution (1M bins) ---
  0M - 1M: 123,456 positions
  1M - 2M: 234,567 positions
  2M - 3M: 345,678 positions
  ...
```

**Also saves**: Full report to `data/diagnostic_report.txt`

---

## Common Issues to Look For

### 1. Duplicate Positions ⚠️

If you see duplicates, investigate why:
```r
# Load one file
data <- readRDS("data/output/5feb_tem_a1_c3_rc.rds")

# Find duplicates
library(data.table)
setDT(data)
duplicates <- data[duplicated(position) | duplicated(position, fromLast = TRUE)]
duplicates[order(position)]
```

**Possible causes:**
- Multiple reads at same position (expected)
- Data processing error (unexpected)
- Different samples merged (expected)

### 2. Large Gaps

Gaps > 10,000 bp might indicate:
- Natural low-coverage regions
- Filtering removed positions
- Assembly gaps

### 3. Different Column Names

If files have different columns:
```r
# Compare
col_a1 <- names(readRDS("data/output/5feb_tem_a1_c3_rc.rds"))
col_a2 <- names(readRDS("data/output/5feb_tem_a2_c3_rc.rds"))

setdiff(col_a1, col_a2)  # In a1 but not a2
setdiff(col_a2, col_a1)  # In a2 but not a1
```

### 4. Unexpected Position Ranges

If position ranges don't make sense:
- Check chromosome (should all be chr3)
- Verify files weren't mixed up
- Check if positions are 0-based or 1-based

---

## Customizing the Scripts

### Check different files:

Edit `files_to_check` in either script:
```r
files_to_check <- c(
  "your_file_1.rds",
  "your_file_2.rds"
)
```

### Different position column:

```r
position_col <- "pos"  # Or whatever your column is named
```

### Different position range:

In `quick_check_files.r`:
```r
# Check 2M-3M instead of 1M-2M
count <- sum(d$position >= 2e6 & d$position <= 3e6, na.rm = TRUE)
```

---

## Workflow

**Step 1: Quick check**
```r
source("scripts/quick_check_files.r")
```

**Step 2: If something looks wrong, run full diagnostic**
```r
source("scripts/diagnose_data_files.r")
```

**Step 3: Review the saved report**
```r
file.show("data/diagnostic_report.txt")
```

**Step 4: Fix issues** (if any)

**Step 5: Create sample data**
```r
source("scripts/create_sample_data.r")
```

---

## Interpreting Results

### Good signs ✅
- No duplicates (or very few)
- No NAs
- Position ranges make sense
- All files have same columns
- Reasonable gaps between positions

### Warning signs ⚠️
- Many duplicates (>1% of data)
- Lots of NAs
- Huge gaps (>100,000 bp)
- Different columns across files
- Position ranges don't overlap (if they should)

### Red flags 🚨
- All positions are duplicated
- Position range is 0-0
- File won't load
- Wrong number of columns
- Positions outside expected chromosome range

---

## Example: Finding the Best Position Range for Sampling

After running the diagnostic:

```r
# Load the diagnostic report
report <- readLines("data/diagnostic_report.txt")

# Or manually check each file
data <- readRDS("data/output/5feb_tem_a1_c3_rc.rds")

# Find densest region
library(data.table)
setDT(data)

# Count positions by 100k bins
data[, bin := floor(position / 1e5) * 1e5]
density <- data[, .N, by = bin][order(-N)]

# Top 10 densest bins
head(density, 10)

# Use densest region for sampling
best_start <- density[1, bin]
best_end <- best_start + 1e6  # 1M range
```

---

## Troubleshooting

**Script runs very slow:**
- Large files take time
- Try quick_check first
- Comment out parts you don't need

**"Column not found" error:**
- Check column names with `names(data)`
- Update `position_col` variable

**Out of memory:**
- Process one file at a time
- Use `data.table::fread()` if CSV
- Increase R memory limit

---

## Next Steps

After diagnostics are clean:
1. ✅ Run `create_sample_data.r`
2. ✅ Run tests: `source('run_tests.R')`
3. ✅ Commit sample data to git
