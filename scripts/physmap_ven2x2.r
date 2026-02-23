# physmap_ven2x2.r
# Replaces running physmap.r x4 then ven2x2.r separately.
# Processes all 4 samples through physmap (C++), saves per-sample outputs,
# then merges and filters to polymorphic sites.

library(Rcpp)
library(data.table)
library(here)

sourceCpp(here("scripts", "functions", "physmap.cpp"))

# set i/o ----
samples <- list(
  list(
    input  = here("data", "sample", "5feb_tem_a1_c3_rc_sample.rds"),
    output = here("data", "output", "5feb_tem_a1_c3_sample.rds"),
    label  = "5feb_tem_a1"
  ),
  list(
    input  = here("data", "sample", "5feb_tem_a2_c3_rc_sample.rds"),
    output = here("data", "output", "5feb_tem_a2_c3_sample.rds"),
    label  = "5feb_tem_a2"
  ),
  list(
    input  = here("data", "sample", "5feb_tem_d1_c3_rc_sample.rds"),
    output = here("data", "output", "5feb_tem_d1_c3_sample.rds"),
    label  = "5feb_tem_d1"
  ),
  list(
    input  = here("data", "sample", "5feb_tem_d2_c3_rc_sample.rds"),
    output = here("data", "output", "5feb_tem_d2_c3_sample.rds"),
    label  = "5feb_tem_d2"
  )
)

rds_out  <- here("data", "output", "5feb_tem_chr3_avd_sample.rds")
text_out <- here("data", "output", "5feb_tem_chr3_avd_sample.tbl")

min_depth <- 25
max_depth <- 1000

# process_sample: run C++ physmap on one splitter output ----
process_sample <- function(input_file) {
  raw_data <- readRDS(input_file)

  mono_sites <- raw_data[is.na(raw_data$snp1) | raw_data$snp1 == "", ]
  poly_sites  <- raw_data[!is.na(raw_data$snp1) & raw_data$snp1 != "", ]
  rm(raw_data)

  mapped_mono <- MapMonoSites(as.matrix(mono_sites))
  mapped_poly <- MapPolySites(as.matrix(poly_sites))
  rm(mono_sites, poly_sites)

  nuc_cols <- c("a", "c", "g", "t", "i", "d")

  dt_mono <- data.table(
    chrom    = as.numeric(mapped_mono[, 1]),
    position = as.numeric(mapped_mono[, 2]),
    refnuc   = as.character(mapped_mono[, 3]),
    a        = as.numeric(mapped_mono[, 4]),
    c        = as.numeric(mapped_mono[, 5]),
    g        = as.numeric(mapped_mono[, 6]),
    t        = as.numeric(mapped_mono[, 7]),
    i        = as.numeric(mapped_mono[, 8]),
    d        = as.numeric(mapped_mono[, 9])
  )
  dt_mono[, (nuc_cols) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)),
          .SDcols = nuc_cols]
  rm(mapped_mono)

  dt_poly <- data.table(
    chrom    = as.numeric(mapped_poly[, 1]),
    position = as.numeric(mapped_poly[, 2]),
    refnuc   = as.character(mapped_poly[, 3]),
    a        = as.numeric(mapped_poly[, 4]),
    c        = as.numeric(mapped_poly[, 5]),
    g        = as.numeric(mapped_poly[, 6]),
    t        = as.numeric(mapped_poly[, 7]),
    i        = as.numeric(mapped_poly[, 8]),
    d        = as.numeric(mapped_poly[, 9])
  )
  dt_poly[, (nuc_cols) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)),
          .SDcols = nuc_cols]
  rm(mapped_poly)

  dt <- rbind(dt_mono, dt_poly)
  rm(dt_mono, dt_poly)

  dt[, sumDepth := a + c + g + t + i + d]
  dt <- dt[sumDepth >= min_depth & sumDepth < max_depth]
  setkey(dt, position)

  dt
}

# process all samples ----
start <- Sys.time()

processed <- lapply(samples, function(s) {
  message("Processing: ", s$label)
  dt <- process_sample(s$input)
  saveRDS(dt, s$output)
  dt
})

# merge by position (inner join — only positions present in all 4 samples) ----
p1_merged <- merge(processed[[1]], processed[[2]],
                   by = "position", suffixes = c("1", "2"))
p2_merged <- merge(processed[[3]], processed[[4]],
                   by = "position", suffixes = c("3", "4"))
pol_all <- merge(p1_merged, p2_merged, by = "position")
rm(p1_merged, p2_merged, processed)

# filter polymorphic sites ----
# a monomorphic position has all reads in one nucleotide per sample:
# 5 zeros per sample x 4 samples = 20 zeros across nucleotide columns.
# keeping rowSums <= 19 retains only positions with at least one alt allele.
nuc_cols_all <- paste0(
  rep(c("a", "c", "g", "t", "i", "d"), times = 4),
  rep(1:4, each = 6)
)
pol <- pol_all[rowSums(pol_all[, ..nuc_cols_all] == 0) <= 19]
rm(pol_all)

# add refnuc column for downstream scripts ----
pol[, refnuc := refnuc1]

end <- Sys.time()
message("Done. Elapsed: ", round(end - start, 2), " ", units(end - start))
message("Polymorphic sites: ", nrow(pol))

# save outputs ----
saveRDS(pol, rds_out)

n_rows <- nrow(pol)
file_conn <- file(text_out, open = "wt")

for (i in seq_len(n_rows)) {
  line1 <- sprintf("%10i", pol$position[i])
  line2 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
    samples[[1]]$label,
    pol$a1[i], pol$c1[i], pol$g1[i], pol$t1[i], pol$i1[i], pol$d1[i])
  line3 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
    samples[[2]]$label,
    pol$a2[i], pol$c2[i], pol$g2[i], pol$t2[i], pol$i2[i], pol$d2[i])
  line4 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
    samples[[3]]$label,
    pol$a3[i], pol$c3[i], pol$g3[i], pol$t3[i], pol$i3[i], pol$d3[i])
  line5 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
    samples[[4]]$label,
    pol$a4[i], pol$c4[i], pol$g4[i], pol$t4[i], pol$i4[i], pol$d4[i])

  writeLines(c(line1, line2, line3, line4, line5), file_conn)
}

close(file_conn)
