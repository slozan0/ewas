# set up environment ----
rm(list = ls())
library(data.table)

# set i/o ----
## input
file_path <- "data/input/5feb_tem_a2.readcounts"
## output
chr1_file <- "data/output/5feb_tem_a2_c1_rc.rds"
chr2_file <- "data/output/5feb_tem_a2_c2_rc.rds"
chr3_file <- "data/output/5feb_tem_a2_c3_rc.rds"

# program starts here ----
##
column_separator_character <- " "
## find the maximum number of columns
n_col <- max(count.fields(file_path, sep = column_separator_character))
n_snps <- n_col - 7
snp_col <- vector()

for (i in 1:n_snps) {
  snp_col[i] <- paste("snp", i, sep = "")
}

start <- Sys.time()

raw_data <- read.table(
  file = file_path,
  header = TRUE,
  sep = column_separator_character,
  fill = TRUE,
  col.names = c("chrom", "position", "refnuc",
                "depth", "q30_depth", "refQA", snp_col),
  stringsAsFactors = FALSE,
  colClasses = c(
    "character",  # chrom
    "integer",    # position (force integer!)
    "character",  # refnuc
    "integer",    # depth
    "integer",    # q30_depth
    "character",  # refQA
    rep("character", length(snp_col))  # snp columns
  )
)


end <- Sys.time()
elapse <- end - start
elapse

raw_data |>
  dplyr::filter(chrom == "1") |>
  saveRDS(file = chr1_file)

raw_data |>
  dplyr::filter(chrom == "2") |>
  saveRDS(file = chr2_file)

raw_data |>
  dplyr::filter(chrom == "3") |>
  saveRDS(file = chr3_file)
