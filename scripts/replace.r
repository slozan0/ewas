rm(list = ls())

library(dplyr)
source("replacer_functions.r")

chiResults <- readRDS("data/output/annotated_5feb_chrom2.rds")

replaceInput <- chiResults |>
  dplyr::select(nucPosition, alleles, 
                group1AltAllFreq,
                group2AltAllFreq,
                lod,
                totalHeteroz,
                group1Heteroz, 
                group2Heteroz,
                info)

replaceInput$info <- gsub("\\s+", " ", replaceInput$info)

replaceInput <- replaceInput %>%
  tidyr::separate(info, 
           into = c("LOC", "POL", "REFN", "GTYPE", "CODPOS", "RESID", "CODON1", "AA1", "CODON2", "AA2"), 
           sep = " ", 
           extra = "merge")

replaceInput$info <- NULL
replaceInput$REPLACEMENT <- NA #make space in ram for new data
replaceInput$NREPS <- NA #make space in ram for new data

# replaceInput1000 <- head(replaceInput, n = 15000)
# nRows <- nrow(replaceInput1000)

# i <- 766 # one indel
# i <- 2735 # triple variant everyting syn
# i <- 204 # triple variant double rep

for(i in 1:nrow(replaceInput)) {
  
  if(replaceInput$GTYPE[i] == "CDS" & !is.na(replaceInput$GTYPE[i])) {
    nucPos <- replaceInput$nucPosition[i]
    
    # print( paste("CDS mutation in nuc pos:", nucPos) )

    refCodon <- replaceInput$CODON1[i]
    refAmino <- replaceInput$AA1[i]
    refNuc <- replaceInput$REFN[i]
    
    varPosition <- as.numeric(replaceInput$CODPOS[i])
    variants <- replaceInput$alleles[i]
    polarity <- replaceInput$POL[i]
      
    variantsTable <- GetSubstitutions(refCodon = refCodon,
                             refAmino = refAmino,
                             refNuc = refNuc,
                             varPosition = varPosition, 
                             variants = variants,
                             polarity = polarity)
    
    if ( nrow(variantsTable ) > 2 ) { #if more than 2 alleles
    
      temp <- ReportTripleVariant( variantsTable )
      replaceInput$CODON2[i] <- temp$codon[2]
      replaceInput$AA2[i] <- temp$amino[2]
      replaceInput$REPLACEMENT[i] <- temp$replacement[2]
      replaceInput$NREPS[i] <- temp$nRep[2]
    
    } else{
      replaceInput$CODON2[i] <- variantsTable$codon[2]
      replaceInput$AA2[i] <- variantsTable$amino[2]
      replaceInput$REPLACEMENT[i] <- variantsTable$replacement[2]
    }
    
    } else {
    # variation not in CDS; do nothing. 
  }
  
}

#save merge results to file
