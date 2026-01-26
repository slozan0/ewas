# setup environment
rm(list = ls())

library(data.table)
library(doParallel)
library(foreach)
library(here)

function_file <- here("scripts", "functions", "easy_chi2_fun.r")
poly_sites_file <- "data/sample/5feb_tem_chr3_avd_sample.rds"

if (!file.exists(function_file)) {
  stop("Cannot find easy_chi2_fun.r")
}
source(function_file)

# set number of processor for parallel processing
number_of_procesors <- 8
# set within group marking threshold
mark_threshold <- 0.05

# set i/o ----
## input files
poly_sites <- readRDS(poly_sites_file)

## output file names
### text output
txt_output_file <- "data/output/5feb_tem_ezchi_c3_sample.chi"
### r binary output
raw_ezchi_results_file <- "data/output/5feb_tem_raw_ezchi_c3_sample.rds"
rds_ezchi_results_file <- "data/output/5feb_tem_ezchi_c3_sample.rds"

# program starts here ----
poly_sites$refnuc <- poly_sites$refnuc1

# remove extra columns ----
cols_to_remove <- c(
  "refnuc1", "refnuc2", "refnuc3", "refnuc4",
  "chrom1", "chrom2", "chrom3", "chrom4",
  "sumDepth1", "sumDepth2", "sumDepth3", "sumDepth4"
)
poly_sites[, (cols_to_remove) := NULL]

# nolint start:
# rowPosition <- which(polySites == 20423119, arr.ind=TRUE)[,"row"]
# polySites[6, ]
# polySite <- (polySites[1, ])
# print(dbug:on)
# nolint end:

# parallel ----
## parameters
n_lines <- nrow(poly_sites)

cl <- parallel::makeCluster(number_of_procesors)
doParallel::registerDoParallel(cl)

raw_ezchi_results <- foreach(i = 1:n_lines, .combine = "rbind") %dopar% {
  get_easy_chi_estimates(poly_site = poly_sites[i, ])
}

stopCluster(cl)

saveRDS(raw_ezchi_results, file = raw_ezchi_results_file, compress = FALSE)

# let's make the matrix a data.table and
# make the nucleotide position a local key
ezchi_results <- as.data.table(raw_ezchi_results, key = "nucPosition")
rm(raw_ezchi_results)
rm(poly_sites)

# let's crate the x2 probs for the contingency tables
ezchi_results[, totalProb := 1 - pchisq(q = totalChiSqr, df = totalDegFreedom),
  by = nucPosition
]

ezchi_results[, group1Prob := 1 -
    pchisq(q = group1ChiSqr, df = group1DegFreedom),
  by = nucPosition
]

ezchi_results[, group2Prob := 1 -
    pchisq(q = group2ChiSqr, df = group2DegFreedom),
  by = nucPosition
]

# estimate Benjamini
bh_threshold <- -log10(get_benjamini_hochber_thresh(ezchi_results$totalProb))

# let's remove the sites with a global probability less than reject threshold
ezchi_results <- ezchi_results[ezchi_results$lod > bh_threshold]

# mark what is inconsistent within the group
ezchi_results[, inconsistency := mark_inconsistency(
  chi1 = group1ChiSqr,
  deg_freedom1 = group1DegFreedom,
  inconsistency_mark1 = "1*",
  chi2 = group2ChiSqr,
  deg_freedom2 = group2DegFreedom,
  inconsistency_mark2 = "2*",
  mark_threshold = mark_threshold
),
by = nucPosition
]

# let's add the allele list the first character is the reference allele
ezchi_results[, alleles := get_alleles_label(
  nuc_position = nucPosition,
  ref_nucleotide = refNuc,
  a_s = As,
  c_s = Cs,
  g_s = Gs,
  t_s = Ts,
  i_s = Is,
  d_s = Ds
),
by = nucPosition
]

saveRDS(ezchi_results, file = rds_ezchi_results_file, compress = FALSE)

n_remaining_sites <- nrow(ezchi_results)

ezchi_results$alleles <- sprintf("%-6s", ezchi_results$alleles)
ezchi_results$alleles <- gsub(" ", "_", ezchi_results$alleles)

# Open a connection to the file
fileconn <- file(txt_output_file, open = "wt")

# Header
header <- paste0(
  "SNPID,MUTATION,FREQ(ALIVE),FREQ(DEAD),LOD,",
  "HET(ALL),HET(ALIVE),HET(DEAD),FREQ(A),FREQ(C),",
  "FREQ(G),FREQ(T),FREQ(I),FREQ(D),CHISQ(ALL),",
  "CHISQ(ALIVE),CHISQ(DEAD),DF(ALL),DF(ALIVE),DF(DEAD)"
)

writeLines(header, fileconn)

print_format_temp <- data.frame(
  nucPosition = "%10.0f",
  alleles = "%6s",
  group1AltAllFreq = "%8.5f",
  group2AltAllFreq = "%8.5f",
  lod = "%7.2f",
  group1Heteroz = "%8.5f",
  group2Heteroz = "%8.5f",
  totalHeteroz = "%8.5f",
  As = "%8.0f",
  Cs = "%7.0f",
  Gs = "%8.0f",
  Ts = "%8.0f",
  Is = "%8.0f",
  Ds = "%8.0f",
  group1ChiSqr = "%11.5f",
  group2ChiSqr = "%10.5f",
  totalChiSqr = "%10.5f",
  group1DegFree = "%5.0f",
  group2DegFree = "%5.0f",
  totalDegFree = "%5.0f",
  inconsistency = "%3s"
)


# Prepare format string
print_format <- paste(print_format_temp[, ], collapse = ",")

# Create all lines at once (vectorized - much faster!)
all_lines <- sprintf(
  print_format,
  ezchi_results$nucPosition,
  ezchi_results$alleles,
  ezchi_results$group1AltAllFreq,
  ezchi_results$group2AltAllFreq,
  ezchi_results$lod,
  ezchi_results$totalHeteroz,
  ezchi_results$group1Heteroz,
  ezchi_results$group2Heteroz,
  ezchi_results$As,
  ezchi_results$Cs,
  ezchi_results$Gs,
  ezchi_results$Ts,
  ezchi_results$Is,
  ezchi_results$Ds,
  ezchi_results$group1ChiSqr,
  ezchi_results$group2ChiSqr,
  ezchi_results$totalChiSqr,
  ezchi_results$group1DegFree,
  ezchi_results$group2DegFree,
  ezchi_results$totalDegFree,
  ezchi_results$inconsistency
)

# Write all lines at once
writeLines(all_lines, fileconn)

# Close the file connection
close(fileconn)
