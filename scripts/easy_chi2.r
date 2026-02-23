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
ezchi_results <- as.data.table(raw_ezchi_results, key = "nuc_position")
rm(raw_ezchi_results)
rm(poly_sites)

# let's crate the x2 probs for the contingency tables
ezchi_results[,
  total_prob := 1 - pchisq(q = total_chi_sqr, df = total_deg_freedom),
  by = nuc_position
]

ezchi_results[,
  group1_prob := 1 - pchisq(q = group1_chi_sqr, df = group1_deg_freedom),
  by = nuc_position
]

ezchi_results[,
  group2_prob := 1 - pchisq(q = group2_chi_sqr, df = group2_deg_freedom),
  by = nuc_position
]

# estimate Benjamini
bh_threshold <- -log10(get_benjamini_hochber_thresh(ezchi_results$total_prob))

# let's remove the sites with a global probability less than reject threshold
ezchi_results <- ezchi_results[ezchi_results$lod > bh_threshold]

# mark what is inconsistent within the group
ezchi_results[, inconsistency := mark_inconsistency(
  chi1 = group1_chi_sqr,
  deg_freedom1 = group1_deg_freedom,
  inconsistency_mark1 = "1*",
  chi2 = group2_chi_sqr,
  deg_freedom2 = group2_deg_freedom,
  inconsistency_mark2 = "2*",
  mark_threshold = mark_threshold
),
by = nuc_position
]

# let's add the allele list the first character is the reference allele
ezchi_results[, alleles := get_alleles_label(
  nuc_position = nuc_position,
  ref_nucleotide = ref_nuc,
  a_s = a_s,
  c_s = c_s,
  g_s = g_s,
  t_s = t_s,
  i_s = i_s,
  d_s = d_s
),
by = nuc_position
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
  nuc_position = "%10.0f",
  alleles = "%6s",
  group1_alt_all_freq = "%8.5f",
  group2_alt_all_freq = "%8.5f",
  lod = "%7.2f",
  group1_heteroz = "%8.5f",
  group2_heteroz = "%8.5f",
  total_heteroz = "%8.5f",
  a_s = "%8.0f",
  c_s = "%7.0f",
  g_s = "%8.0f",
  t_s = "%8.0f",
  i_s = "%8.0f",
  d_s = "%8.0f",
  group1_chi_sqr = "%11.5f",
  group2_chi_sqr = "%10.5f",
  total_chi_sqr = "%10.5f",
  group1_deg_freedom = "%5.0f",
  group2_deg_freedom = "%5.0f",
  total_deg_freedom = "%5.0f",
  inconsistency = "%3s"
)


# Prepare format string
print_format <- paste(print_format_temp[, ], collapse = ",")

# Create all lines at once (vectorized - much faster!)
all_lines <- sprintf(
  print_format,
  ezchi_results$nuc_position,
  ezchi_results$alleles,
  ezchi_results$group1_alt_all_freq,
  ezchi_results$group2_alt_all_freq,
  ezchi_results$lod,
  ezchi_results$total_heteroz,
  ezchi_results$group1_heteroz,
  ezchi_results$group2_heteroz,
  ezchi_results$a_s,
  ezchi_results$c_s,
  ezchi_results$g_s,
  ezchi_results$t_s,
  ezchi_results$i_s,
  ezchi_results$d_s,
  ezchi_results$group1_chi_sqr,
  ezchi_results$group2_chi_sqr,
  ezchi_results$total_chi_sqr,
  ezchi_results$group1_deg_freedom,
  ezchi_results$group2_deg_freedom,
  ezchi_results$total_deg_freedom,
  ezchi_results$inconsistency
)

# Write all lines at once
writeLines(all_lines, fileconn)

# Close the file connection
close(fileconn)
