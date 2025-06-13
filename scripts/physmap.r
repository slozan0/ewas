#set up env ----
rm(list = ls())

library(Rcpp)
library(data.table)
sourceCpp("physmap.cpp")

#set i/o ----
inputFileName <- "data/input/5feb_tem_a1_c2_rc.rds"
outFileName   <- "data/output/5feb_tem_a1_c2.rds"

#program starts here ----
rawData <- readRDS(file = inputFileName)

monoSites <- rawData[ rawData$snp1 ==  "", ]
monoSites <- as.matrix( monoSites)

polySites <- rawData[ rawData$snp1 != "", ]
polySites <- as.matrix( polySites)

start <- Sys.time()

myMonoLines <- MapMonoSites(monoSites)

myPolyLines <- MapPolySites(polySites)
rm(rawData)

#convert to dt object and change NAs to zeros
dtMonoSites <- data.table(chrom = as.numeric(myMonoLines[ , 1]),
                          pos   = as.numeric(myMonoLines[ , 2]),
                          ref   = as.character(myMonoLines[ , 3]),
                          a     = as.numeric(myMonoLines[ , 4]),
                          c     = as.numeric(myMonoLines[ , 5]),
                          g     = as.numeric(myMonoLines[ , 6]),
                          t     = as.numeric(myMonoLines[ , 7]),
                          i     = as.numeric(myMonoLines[ , 8]),
                          d     = as.numeric(myMonoLines[ , 9]),
                          key = c("chrom","pos")
                                         )
dtMonoSites[is.na(dtMonoSites)] <- 0
rm(monoSites, myMonoLines)

dtPolySites <- data.table(chrom = as.numeric(myPolyLines[,1]),
                          pos   = as.numeric(myPolyLines[,2]),
                          ref   = as.character(myPolyLines[,3]),
                          a     = as.numeric(myPolyLines[,4]),
                          c     = as.numeric(myPolyLines[,5]),
                          g     = as.numeric(myPolyLines[,6]),
                          t     = as.numeric(myPolyLines[,7]),
                          i     = as.numeric(myPolyLines[,8]),
                          d     = as.numeric(myPolyLines[,9]),
                          key = c("chrom","pos")
)
dtPolySites[is.na(dtPolySites)] <- 0
rm(polySites, myPolyLines)

end <- Sys.time()
elapse <- end - start
elapse

dtChrom <- rbind(dtPolySites, dtMonoSites)

#sum the number of reads per nucleotide and create a new field
dtChrom <- dtChrom[ , sumDepth := sum(a,c,g,t,i,d) , by = 1:NROW(dtChrom) ]

# remove positions with less than 25 reads or more than 1000

dtChrom <- dtChrom[ sumDepth >= 25 ] 
dtChrom <- dtChrom[ sumDepth < 1000 ]

saveRDS(dtChrom, outFileName)
