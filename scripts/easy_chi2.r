# setup environment
rm(list = ls())

library(data.table)
library(doParallel)
library(foreach)
source("easy_chi2_fun.r")

## set number of processor for parallel processing
numberOfProcesors <- 8 
## set within group marking threshold 
markThreshold <- 0.05

#set i/o ----
## input files
polySites <- readRDS("data/5feb_tem_chr1_avd.rds")
## input file names
### text output
txtOutputFile   <- "data/output/5feb_tem_ezchi_c1.chi"
### r binary output
rawEzchiResults <- "data/output/5feb_tem_raw_ezchi_c1.rds"
rdsEzChiResults <- "data/output/5feb_tem_ezchi_c1.rds"

# program starts here ----
polySites$ref <- polySites$ref1

# remove extra columns ----
polySites$ref1 <- NULL
polySites$ref2 <- NULL
polySites$ref3 <- NULL
polySites$ref4 <- NULL

polySites$chrom1 <- NULL
polySites$chrom2 <- NULL
polySites$chrom3 <- NULL
polySites$chrom4 <- NULL

polySites$sumDepth1 <- NULL
polySites$sumDepth2 <- NULL
polySites$sumDepth3 <- NULL
polySites$sumDepth4 <- NULL

#dbug ----
#rowPosition <- which(polySites == 20423119, arr.ind=TRUE)[,"row"]
#polySites[6, ]
#polySite <- (polySites[1, ])
#print(dbug:on)
#dbug ----

# parallel ----
# parameters ----
nLines <- nrow(polySites)

cl <- parallel::makeCluster(numberOfProcesors)
doParallel::registerDoParallel(cl)

rawEzChiResults <- foreach( i = 1:nLines, .combine = 'rbind') %dopar% {
  GetEasyChiEstimates(polySite = polySites[i, ])
}

stopCluster(cl)



saveRDS(rawEzChiResults, file = rawEzchiResults, compress = FALSE)

#let's make the matrix a data.table and make the nucleotide position a local key
ezChiResults <- as.data.table(rawEzChiResults, key = "nucPosition")
rm(rawEzChiResults)
rm(polySites)

#let's crate the x2 probs for the contingency tables
ezChiResults[ , totalProb := 1 - pchisq( q = totalChiSqr, df = totalDegFreedom ) ,
              by = nucPosition ]

ezChiResults[ , group1Prob := 1 - pchisq( q = group1ChiSqr, df = group1DegFreedom ) ,
              by = nucPosition ]

ezChiResults[ , group2Prob := 1 - pchisq( q = group2ChiSqr, df = group2DegFreedom ) ,
              by = nucPosition ]


#estimate Benjamini
bhThreshold <- -log10( GetBenjaminiHochberThreshold(ezChiResults$totalProb))

#let's remove the sites with a global probability less than reject threshold
ezChiResults <- ezChiResults[ ezChiResults$lod > bhThreshold ]

# mark what is inconsistent within the group 
ezChiResults[ , inconsistency := MarkInconsistency(chi1 = group1ChiSqr,
                                                   degFreedom1 = group1DegFreedom,
                                                   inconsistencyMark1 = "1*",
                                                   chi2 = group2ChiSqr ,
                                                   degFreedom2 = group2DegFreedom,
                                                   inconsistencyMark2 = "2*",
                                                   markThreshold = markThreshold),
              by = nucPosition ]


#let's add the allele list the first character is the reference allele
ezChiResults[ , alleles := GetAllelesLabel(
  nucPosition = nucPosition,
  refNucleotide = refNuc, 
  As = As,
  Cs = Cs,
  Gs = Gs,
  Ts = Ts,
  Is = Is,
  Ds = Ds),
  by = nucPosition ]

saveRDS(ezChiResults, file = rdsEzChiResults, compress = FALSE)

nRemainingSites <- nrow(ezChiResults)

ezChiResults$alleles <- sprintf("%-6s", ezChiResults$alleles)
ezChiResults$alleles <- gsub(" ", "_", ezChiResults$alleles)

# Open a connection to the file
fileConn <- file(txtOutputFile, open = "wt")

# Header
header <- "SNPID,MUTATION,FREQ(ALIVE),FREQ(DEAD),LOD,HET(ALL),HET(ALIVE),HET(DEAD),FREQ(A),FREQ(C),FREQ(G),FREQ(T),FREQ(I),FREQ(D),CHISQ(ALL),CHISQ(ALIVE),CHISQ(DEAD),DF(ALL),DF(ALIVE),DF(DEAD)"
writeLines(header, fileConn)

printFormatTemp <- data.frame(nucPosition  = "%10.0f",
                              alleles          = "%6s",
                              group1AltAllFreq = "%8.5f",
                              group2AltAllFreq = "%8.5f",
                              lod              = "%7.2f",
                              group1Heteroz    = "%8.5f",
                              group2Heteroz    = "%8.5f",
                              totalHeteroz     = "%8.5f",
                              As               = "%8.0f",
                              Cs               = "%7.0f",
                              Gs               = "%8.0f",
                              Ts               = "%8.0f",
                              Is               = "%8.0f",
                              Ds               = "%8.0f",
                              group1ChiSqr     = "%11.5f",
                              group2ChiSqr     = "%10.5f",
                              totalChiSqr      = "%10.5f",
                              group1DegFree    = "%5.0f",
                              group2DegFree    = "%5.0f",
                              totalDegFree     = "%5.0f",
                              inconsistency    = "%3s"
)


# Prepare format string
printFormat <- paste(printFormatTemp[, ], collapse = ",")

# Write each row to the file
for (i in 1:nRemainingSites) {
  line <- sprintf(printFormat,
                  ezChiResults$nucPosition[i],     #1
                  ezChiResults$alleles[i],         #2
                  ezChiResults$group1AltAllFreq[i],#3
                  ezChiResults$group2AltAllFreq[i],#4
                  ezChiResults$lod[i],             #5
                  ezChiResults$totalHeteroz[i],    #6
                  ezChiResults$group1Heteroz[i],   #7
                  ezChiResults$group2Heteroz[i],   #8
                  ezChiResults$As[i],              #9
                  ezChiResults$Cs[i],              #10
                  ezChiResults$Gs[i],              #11
                  ezChiResults$Ts[i],              #12
                  ezChiResults$Is[i],              #13
                  ezChiResults$Ds[i],              #14
                  ezChiResults$group1ChiSqr[i],    #15
                  ezChiResults$group2ChiSqr[i],    #16
                  ezChiResults$totalChiSqr[i],     #17
                  ezChiResults$group1DegFree[i],   #18
                  ezChiResults$group2DegFree[i],   #19
                  ezChiResults$totalDegFree[i],    #20
                  ezChiResults$inconsistency[i]    #21
  )
  writeLines(line, fileConn)
}

# Close the file connection
close(fileConn)


