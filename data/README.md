# Data Directory Structure

This directory contains datasets for the EWAS PhysMap project. Data is organized to separate shareable sample data from large real datasets.

## Directory Layout

```
data/
├── sample/          # Small test datasets (committed to git)
├── input/           # Real/complete data (ignored by git)
├── output/          # Generated results (ignored by git)
└── cache/           # Intermediate results (ignored by git)
```

---

## 📁 `sample/` - Sample Data (Committed)

**Purpose**: Small datasets for testing, CI/CD, and sharing with collaborators

**Contents**:
- Small subset of real data (~100-1000 rows)
- Unit test fixtures
- Expected outputs for regression tests
- Example data showing structure/schema

**Size limit**: Keep individual files < 10 MB

**Usage**:
```r
# Load sample data for testing
test_data <- readRDS("data/sample/test_chr1.rds")

# Run tests
source("run_tests.R")
```

---

## 📁 `input/` - Real Data (NOT Committed)

**Purpose**: Full datasets for actual analysis

**Contents**:
- Raw sequencing data (`.readcounts`, `.fastq`, etc.)
- Pre-processed RDS files
- Annotation databases (DuckDB, SQLite)
- Large reference files

**Location on your machine**:
- Local: `data/input/`
- Server/HPC: (document your data location here)

**How to populate**:
```bash
# Option 1: Copy from server/external drive
cp /path/to/real/data/*.rds data/input/

# Option 2: Download from cloud storage
# (add your download script here)

# Option 3: Symlink to avoid duplication
ln -s /path/to/central/data data/input
```

**⚠️ Important**:
- Never commit these files (they're in .gitignore)
- Keep a backup in a separate location
- Document data sources in this README

---

## 📁 `output/` - Generated Results (NOT Committed)

**Purpose**: Store analysis results and figures

**Contents**:
- Chi-square test results (`.chi` files)
- Annotated datasets (`.rds` files)
- Figures and plots
- Summary statistics

**Cleaning**:
```bash
# Remove all output files (BE CAREFUL!)
rm -rf data/output/*

# Or use R
unlink("data/output", recursive = TRUE)
dir.create("data/output")
```

---

## 📁 `cache/` - Cached Results (NOT Committed)

**Purpose**: Store intermediate results to speed up reruns

**Contents**:
- Partially processed data
- Expensive computation results
- Temporary files

**Cleaning**:
```bash
# Safe to delete - will be regenerated
rm -rf data/cache/*
```

---

## Creating Sample Data

To create sample data from your real datasets:

```r
# Example: Create a small sample from chromosome 1
full_data <- readRDS("data/input/5feb_tem_chr1_avd.rds")

# Take first 500 rows
sample_data <- full_data[1:500, ]

# Save to sample directory
saveRDS(sample_data, "data/sample/test_chr1.rds")

# Verify size
file.info("data/sample/test_chr1.rds")$size / 1024^2  # MB
```

### Anonymizing Sensitive Data

If your data contains sensitive information:

```r
# Remove or scramble sensitive columns
sample_data$patient_id <- NULL
sample_data$location <- sample(sample_data$location)  # Shuffle

# Or generate synthetic data with same structure
synthetic <- data.frame(
  pos = sample(1:1000000, 500),
  ref = sample(c("A", "C", "G", "T"), 500, replace = TRUE),
  # ... match your real data structure
)
```

---

## Data Requirements

### Input Data Format

**For `poly_sites`**:
- Columns: `pos`, `ref`, `ref1-ref4`, `chrom1-chrom4`, nucleotide counts (A, C, G, T, I, D) for each group
- Format: RDS (data.table or data.frame)

**For annotations**:
- DuckDB database with `chrom_data` table
- Columns: `identifier` (INTEGER PRIMARY KEY), `info` (TEXT)

### Output Data Format

**Chi-square results (`.chi`)**:
- CSV format with header
- Columns: SNPID, MUTATION, FREQ(ALIVE), FREQ(DEAD), LOD, etc.

**RDS results**:
- data.table with chi-square statistics
- Columns: nucPosition, refNuc, group1AltAllFreq, group2AltAllFreq, lod, etc.

---

## Troubleshooting

### "File not found" errors

```r
# Check if data directory exists
dir.exists("data/input")

# List available files
list.files("data/input", recursive = TRUE)
```

### Large file issues

If sample files are too large:
```r
# Reduce sample size
sample_data <- full_data[1:100, ]  # Even smaller sample

# Or compress
saveRDS(sample_data, "data/sample/test.rds", compress = "xz")
```

### Missing real data

If you're setting up on a new machine:
1. Check .gitignore to see what's excluded
2. Copy real data from backup/server
3. Verify file permissions
4. Update paths in scripts if needed

---

## Data Provenance

### Data Sources

| File | Source | Date | Description |
|------|--------|------|-------------|
| `5feb_tem_a1.readcounts` | [Lab/Source] | 2023-XX-XX | Raw read counts from... |
| `chrom1_3.duckdb` | [Database/Source] | 2023-XX-XX | Annotation data for... |

### Processing History

Document major processing steps:
1. Raw `.readcounts` → Split by chromosome → `_c1.rds`, `_c2.rds`, `_c3.rds`
2. Apply variant filters → `_avd.rds`
3. Chi-square analysis → `_ezchi.rds`
4. Annotation → `annotated_*.rds`

---

## Best Practices

✅ **DO**:
- Keep sample data small (<10 MB per file)
- Document data sources and processing steps
- Test with sample data before running on full dataset
- Use relative paths (`data/input/...`) not absolute paths

❌ **DON'T**:
- Commit large files to git
- Hard-code absolute file paths
- Commit sensitive/private data
- Delete `data/input/` without backup!

---

## Questions?

If you're collaborating on this project and need access to the full datasets, contact:
- [Your name/email]
- Data location: [Server path or cloud storage link]
