# set up env ----
rm(list = ls())

library(Rcpp)
library(data.table)
library(here)

# Install "BH" library if you haven't install.packages("BH")
sourceCpp(here("scripts", "functions", "physmap.cpp"))

# set i/o ----
input_file_name <- "data/sample/5feb_tem_d2_c3_sample.rds"
out_file_name <- "data/output/5feb_tem_d2_c3.rds"

# program starts here ----
raw_data <- readRDS(file = input_file_name)

mono_sites <- raw_data[raw_data$snp1 == "", ]
mono_sites <- as.matrix(mono_sites)

poly_sites <- raw_data[raw_data$snp1 != "", ]
poly_sites <- as.matrix(poly_sites)

start <- Sys.time()

my_mono_lines <- MapMonoSites(mono_sites)

my_poly_lines <- MapPolySites(poly_sites)
rm(raw_data)

# convert to dt object and change NAs to zeros
dt_mono_sites <- data.table(
  chrom = as.numeric(my_mono_lines[, 1]),
  pos = as.numeric(my_mono_lines[, 2]),
  ref = as.character(my_mono_lines[, 3]),
  a = as.numeric(my_mono_lines[, 4]),
  c = as.numeric(my_mono_lines[, 5]),
  g = as.numeric(my_mono_lines[, 6]),
  t = as.numeric(my_mono_lines[, 7]),
  i = as.numeric(my_mono_lines[, 8]),
  d = as.numeric(my_mono_lines[, 9]),
  key = c("chrom", "pos")
)

dt_mono_sites[is.na(dt_mono_sites)] <- 0
rm(mono_sites, my_mono_lines)

dt_poly_sites <- data.table(
  chrom = as.numeric(my_poly_lines[, 1]),
  pos = as.numeric(my_poly_lines[, 2]),
  ref = as.character(my_poly_lines[, 3]),
  a = as.numeric(my_poly_lines[, 4]),
  c = as.numeric(my_poly_lines[, 5]),
  g = as.numeric(my_poly_lines[, 6]),
  t = as.numeric(my_poly_lines[, 7]),
  i = as.numeric(my_poly_lines[, 8]),
  d = as.numeric(my_poly_lines[, 9]),
  key = c("chrom", "pos")
)
dt_poly_sites[is.na(dt_poly_sites)] <- 0
rm(poly_sites, my_poly_lines)

end <- Sys.time()
elapse <- end - start
elapse

dt_chrom <- rbind(dt_poly_sites, dt_mono_sites)

# sum the number of reads per nucleotide and create a new field
dt_chrom <- dt_chrom[, sumDepth := sum(a, c, g, t, i, d),
  by = seq_len(NROW(dt_chrom))
]

# remove positions with less than 25 reads or more than 1000

dt_chrom <- dt_chrom[sumDepth >= 25]
dt_chrom <- dt_chrom[sumDepth < 1000]

saveRDS(dt_chrom, out_file_name)
