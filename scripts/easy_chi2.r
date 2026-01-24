# setup environment
rm(list = ls())

library(data.table)
library(doParallel)
library(foreach)
source("easy_chi2_fun.r")

# set number of processor for parallel processing
number_of_procesors <- 8
# set within group marking threshold
mark_threshold <- 0.05

# set i/o ----
## input files
poly_sites <- readRDS("data/5feb_tem_chr1_avd.rds")
## input file names
### text output
txt_output_file <- "data/output/5feb_tem_ezchi_c1.chi"
### r binary output
raw_ezchi_results <- "data/output/5feb_tem_raw_ezchi_c1.rds"
rds_ezchi_results <- "data/output/5feb_tem_ezchi_c1.rds"

# program starts here ----
poly_sites$ref <- poly_sites$ref1

# remove extra columns ----
poly_sites$ref1 <- NULL
poly_sites$ref2 <- NULL
poly_sites$ref3 <- NULL
poly_sites$ref4 <- NULL

poly_sites$chrom1 <- NULL
poly_sites$chrom2 <- NULL
poly_sites$chrom3 <- NULL
poly_sites$chrom4 <- NULL

poly_sites$sumDepth1 <- NULL
poly_sites$sumDepth2 <- NULL
poly_sites$sumDepth3 <- NULL
poly_sites$sumDepth4 <- NULL

# nolint start:
# rowPosition <- which(polySites == 20423119, arr.ind=TRUE)[,"row"]
# polySites[6, ]
# polySite <- (polySites[1, ])
# print(dbug:on)
# nolint end:

# parallel ----
# parameters ----
n_lines <- nrow(poly_sites)

cl <- parallel::makeCluster(number_of_procesors)
doParallel::registerDoParallel(cl)

raw_ezchi_results <- foreach(i = 1:n_lines, .combine = "rbind") %dopar% {
  get_easy_chi_estimates(poly_site = poly_sites[i, ])
}

stopCluster(cl)

saveRDS(raw_ezchi_results, file = raw_ezchi_results, compress = FALSE)

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
  inconsistency_bark1 = "1*",
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

saveRDS(ezchi_results, file = rds_ezchi_results, compress = FALSE)

n_remaining_sites <- nrow(ezchi_results)

ezchi_results$alleles <- sprintf("%-6s", ezchi_results$alleles)
ezchi_results$alleles <- gsub(" ", "_", ezchi_results$alleles)

# Open a connection to the file
fileconn <- file(txt_output_file, open = "wt")

# Header
header <-
  "SNPID,MUTATION,FREQ(ALIVE),FREQ(DEAD),LOD,HET(ALL),HET(ALIVE),HET(DEAD),FREQ(A),FREQ(C),FREQ(G),FREQ(T),FREQ(I),FREQ(D),CHISQ(ALL),CHISQ(ALIVE),CHISQ(DEAD),DF(ALL),DF(ALIVE),DF(DEAD)" # nolint: line_length_linter.
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

# Write each row to the file
for (i in 1:n_remaining_sites) {
  line <- sprintf(
    print_format,
    ezchi_results$nucPosition[i], # 1
    ezchi_results$alleles[i], # 2
    ezchi_results$group1AltAllFreq[i], # 3
    ezchi_results$group2AltAllFreq[i], # 4
    ezchi_results$lod[i], # 5
    ezchi_results$totalHeteroz[i], # 6
    ezchi_results$group1Heteroz[i], # 7
    ezchi_results$group2Heteroz[i], # 8
    ezchi_results$As[i], # 9
    ezchi_results$Cs[i], # 10
    ezchi_results$Gs[i], # 11
    ezchi_results$Ts[i], # 12
    ezchi_results$Is[i], # 13
    ezchi_results$Ds[i], # 14
    ezchi_results$group1ChiSqr[i], # 15
    ezchi_results$group2ChiSqr[i], # 16
    ezchi_results$totalChiSqr[i], # 17
    ezchi_results$group1DegFree[i], # 18
    ezchi_results$group2DegFree[i], # 19
    ezchi_results$totalDegFree[i], # 20
    ezchi_results$inconsistency[i] # 21
  )
  writeLines(line, fileconn)
}

# Close the file connection
close(fileconn)
