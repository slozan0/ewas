library(data.table)

mat_to_dt <- function(m) {
  nuc_cols <- c("a", "c", "g", "t", "i", "d")
  dt <- data.table(
    chrom    = as.integer(m[, 1]),
    position = as.integer(m[, 2]),
    refnuc   = as.character(m[, 3]),
    a        = as.integer(m[, 4]),
    c        = as.integer(m[, 5]),
    g        = as.integer(m[, 6]),
    t        = as.integer(m[, 7]),
    i        = as.integer(m[, 8]),
    d        = as.integer(m[, 9])
  )
  dt[, (nuc_cols) := lapply(.SD, function(x) ifelse(is.na(x), 0L, x)),
     .SDcols = nuc_cols]
  dt
}

# process_sample: run C++ physmap on one splitter output ----
# Requires MapMonoSites and MapPolySites from physmap.cpp via sourceCpp.
process_sample <- function(input_file, min_depth, max_depth) {
  raw_data <- readRDS(input_file)

  mono_sites <- raw_data[is.na(raw_data$snp1) | raw_data$snp1 == "", ]
  poly_sites  <- raw_data[!is.na(raw_data$snp1) & raw_data$snp1 != "", ]
  rm(raw_data)

  mapped_mono <-
    MapMonoSites(as.matrix(mono_sites)) # nolint: object_usage_linter.
  mapped_poly <-
    MapPolySites(as.matrix(poly_sites)) # nolint: object_usage_linter.
  rm(mono_sites, poly_sites)

  dt <- rbind(mat_to_dt(mapped_mono), mat_to_dt(mapped_poly))
  rm(mapped_mono, mapped_poly)

  dt[, sumDepth := a + c + g + t + i + d]
  dt <- dt[sumDepth >= min_depth & sumDepth < max_depth]
  setkey(dt, position)

  dt
}
