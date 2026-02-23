# physmap.r
# Replaces running physmap.r x4 then ven2x2.r separately.
# Processes all 4 samples through physmap (C++), saves per-sample outputs,
# then merges and filters to polymorphic sites.

library(Rcpp)
library(data.table)
library(here)

sourceCpp(here("scripts", "functions", "physmap.cpp"))
source(here("scripts", "functions", "physmap_fun.r"))

# set i/o ----
samples <- list(
  list(
    input  = here("data", "input", "5feb_tem_a1_c3_rc.rds"),
    output = here("data", "output", "5feb_tem_a1_c3_sample.rds"),
    label  = "5feb_tem_a1"
  ),
  list(
    input  = here("data", "input", "5feb_tem_a2_c3_rc.rds"),
    output = here("data", "output", "5feb_tem_a2_c3_sample.rds"),
    label  = "5feb_tem_a2"
  ),
  list(
    input  = here("data", "input", "5feb_tem_d1_c3_rc.rds"),
    output = here("data", "output", "5feb_tem_d1_c3_sample.rds"),
    label  = "5feb_tem_d1"
  ),
  list(
    input  = here("data", "input", "5feb_tem_d2_c3_rc.rds"),
    output = here("data", "output", "5feb_tem_d2_c3_sample.rds"),
    label  = "5feb_tem_d2"
  )
)

rds_out  <- here("data", "output", "5feb_tem_chr3_avd.rds")
text_out <- here("data", "output", "5feb_tem_chr3_avd.tbl")

min_depth <- 25
max_depth <- 1000

dir.create(here("data", "output"), recursive = TRUE, showWarnings = FALSE)

# process all samples ----
start <- Sys.time()

processed <- lapply(samples, function(s) {
  message("Processing: ", s$label)
  dt <- process_sample(s$input, min_depth, max_depth)
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

# add refnuc column and drop per-sample auxiliary columns ----
pol[, refnuc := refnuc1]
cols_to_drop <- c(
  "refnuc1", "refnuc2", "refnuc3", "refnuc4",
  "chrom1", "chrom2", "chrom3", "chrom4",
  "sumDepth1", "sumDepth2", "sumDepth3", "sumDepth4"
)
pol[, (cols_to_drop) := NULL]

# save outputs ----
saveRDS(pol, rds_out)

line1 <- sprintf("%10i", pol$position)
line2 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
  samples[[1]]$label,
  pol$a1, pol$c1, pol$g1, pol$t1, pol$i1, pol$d1
)
line3 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
  samples[[2]]$label,
  pol$a2, pol$c2, pol$g2, pol$t2, pol$i2, pol$d2
)
line4 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
  samples[[3]]$label,
  pol$a3, pol$c3, pol$g3, pol$t3, pol$i3, pol$d3
)
line5 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
  samples[[4]]$label,
  pol$a4, pol$c4, pol$g4, pol$t4, pol$i4, pol$d4
)

file_conn <- file(text_out, open = "wt")
writeLines(as.vector(rbind(line1, line2, line3, line4, line5)), file_conn)
close(file_conn)

end <- Sys.time()
message("Done. Elapsed: ", round(end - start, 2), " ", units(end - start))
message("Polymorphic sites: ", nrow(pol))
