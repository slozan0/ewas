#set up environment ----
rm(list = ls())
library(data.table)

#set i/o ----
## input
filePath <- "data/input/5feb_tem_a1.readcounts"
## output
chr1File <- "data/output/5feb_tem_a1_c1_rc.rds"
chr2File <- "data/output/5feb_tem_a1_c2_rc.rds"
chr3File <- "data/output/5feb_tem_a1_c3_rc.rds"

# program starts here ----
## 
columnSeparatorCharacter <- " "
##find the maximum number of columns
nCol   <- max( count.fields(filePath, sep = columnSeparatorCharacter) )
nSnps  <- nCol - 7
snpCol <- vector()

for (i in 1:nSnps) {
  snpCol[i] <- paste("snp", i, sep = "")  
}

start <- Sys.time()
rawData <- read.table(file = filePath, 
                      sep = columnSeparatorCharacter, 
                      fill = TRUE,
                      col.names = c("chrom", "position", "refnuc", "depth", "q30_depth", "refQA", snpCol ),
                      stringsAsFactors = FALSE)

end <- Sys.time()
elapse <- end - start
elapse

rawData |> dplyr::filter(chrom == "1") |>
  saveRDS( file = chr1File)

rawData |> dplyr::filter(chrom == "2") |>
  saveRDS( file = chr2File)

rawData |> dplyr::filter(chrom == "3") |>
  saveRDS( file = chr3File)
