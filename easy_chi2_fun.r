GetHetero <- function(counts, nGroups) {
  #dbug----
  #groups <- allGroup
  #nGroups <- 4
  #dbug----

  a <- sum( counts[ , 1] )
  c <- sum( counts[ , 2] )
  g <- sum( counts[ , 3] )
  t <- sum( counts[ , 4] )
  i <- sum( counts[ , 5] )
  d <- sum( counts[ , 6] )
  
  SUMNUC = a + c + g + t + i + d;
  
  freqA = a / SUMNUC;
  freqC = c / SUMNUC;
  freqG = g / SUMNUC;
  freqT = t / SUMNUC;
  freqI = i / SUMNUC;
  freqD = d / SUMNUC;
  
  homozygocity = (freqA * freqA) + (freqC * freqC) + (freqG * freqG) + (freqT * freqT) + 
                 (freqI * freqI) + (freqD * freqD);
  
  heterozygosity = 1 - homozygocity;
  
  return(heterozygosity);
}

GetFrequencies <- function(observed, total, altAlleles) {
  #debug----
  #observed <- ROBS
  #total <- T1 + T2
  #altAlleles <- altAlleles
  #end debug----
  
  nAllelesG1 <- 0
  nAllelesG2 <- 0

  freqAll <- vector( mode = "numeric", length = 6 )

  freqGroupOneTemp <- observed[3, ] / sum( observed[3, ] )
  freqGroupTwoTemp <- observed[6, ] / sum( observed[6, ] )

    #number of alleles in nuc position
  for (i in 1:6) {
    if( freqGroupOneTemp[i] > 0.0 ){ 
      nAllelesG1 <- nAllelesG1 + 1 
      }
    if( freqGroupTwoTemp[i] > 0.0 ){ 
      nAllelesG2 <- nAllelesG2 + 1 
      }
  }
  
  if(nAllelesG1 != nAllelesG2) warning("Freq. Estimation: The number of alleles is not equal between the groups.")
  
  freqGroupOne <- sum(freqGroupOneTemp[altAlleles])
  freqGroupTwo <- sum(freqGroupTwoTemp[altAlleles])
  
  results <- list(
    RF1 = freqGroupOne,
    RF2 = freqGroupTwo
  )
  
  return(results)
}

GetChi <- function(nucPosition, observed, wObs1, wObs2) {
  
  #dbug
   #rObs <- ROBS[3, ]
   #wObs1 = WOBS[1, ]
   #wObs2 = WOBS[2, ]
  #dbug end
  
  removeChar <- ""
  nAlleles <- 0
  
  expectedValuesGroupOne <- vector(mode = "numeric", length = 7)
  expectedValuesGroupTwo <- vector(mode = "numeric", length = 7)
  
  chiGroupOne <- vector(mode = "numeric", length = 7)
  chiGroupTwo <- vector(mode = "numeric", length = 7)
  
  columnProportion <- vector(mode = "numeric", length = 7)
  
  for (i in 1:6) {
    if(observed[i] > 0) nAlleles = nAlleles + 1
  }

  total <- sum(observed)
  
  columnProportion <- observed / total

  freqGroupOne <- sum(wObs1) / total 
  freqGroupTwo <- sum(wObs2) / total 
  
  expectedValuesGroupOne <- freqGroupOne * columnProportion * total
  expectedValuesGroupTwo <- freqGroupTwo * columnProportion * total

  chiGroupOne <- (( wObs1 - expectedValuesGroupOne ) ^ 2 ) / expectedValuesGroupOne
  chiGroupTwo <- (( wObs2 - expectedValuesGroupTwo ) ^ 2 ) / expectedValuesGroupTwo

  chiSqr <- sum(chiGroupOne, chiGroupTwo, na.rm = TRUE)
  degFreedom <- nAlleles - 1
  
  # results ----  
  results <- list(chiSqr = chiSqr,
                  degFreedom = degFreedom)
  
  return(results)
}

GetChiAll <- function(observed, total) {
  
  #dbug
  #observed <- ROBS
  #total <- T1 + T2
  #dbug end
  
  freqTotal <- vector(mode = "numeric", length = 7)
  expectedValueGroup1 <- vector(mode = "numeric", length = 7)
  expectedValueGroup2 <- vector(mode = "numeric", length = 7)
  chiGroup1 <- vector(mode = "numeric", length = 7)
  chiGroup2 <- vector(mode = "numeric", length = 7)
  
  nAlleles = 0
  
  if( observed[7, 1] > 0 ) nAlleles = nAlleles + 1
  if( observed[7, 2] > 0 ) nAlleles = nAlleles + 1
  if( observed[7, 3] > 0 ) nAlleles = nAlleles + 1
  if( observed[7, 4] > 0 ) nAlleles = nAlleles + 1
  if( observed[7, 5] > 0 ) nAlleles = nAlleles + 1
  if( observed[7, 6] > 0 ) nAlleles = nAlleles + 1
  
  freqTotal <- observed[7, ] / total

  freqGroupOne = sum( observed[3, ] ) / total
  freqGroupTwo = sum( observed[6, ] ) / total
  
  expectedValueGroup1 <- freqGroupOne * freqTotal * total
  expectedValueGroup2 <- freqGroupTwo * freqTotal * total
  
  chiGroup1 <- ( (observed[3, ] - expectedValueGroup1 )^2 ) / expectedValueGroup1
  chiGroup2 <- ( (observed[6, ] - expectedValueGroup2 )^2 ) / expectedValueGroup2
  
  chiAll = sum(chiGroup1 + chiGroup2, na.rm = TRUE)
  
  degFreedom = nAlleles - 1
  
  results <- list(
    TCHISQ3 = chiAll,
    IDF3 = degFreedom
  )

  return(results)
}

GetAlternateAlleles <- function(refNuc, observed) {
  
  #debug
  #refNuc <- refNuc
  #observed <- ROBS[7,]
  #end dbug
  
  altNucleotides <- matrix(data = FALSE, nrow = 1, ncol = 6)
  colnames(altNucleotides) <- c("A", "C", "G", "T", "I", "D")
  #mark the reference nuc/allele as not being the alternate
  if(refNuc == "A"){
    altNucleotides <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  } else if(refNuc == "C"){
    altNucleotides <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  } else if(refNuc == "G"){
    altNucleotides <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
  } else if(refNuc == "T"){
    altNucleotides <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  } else if(refNuc == "I"){
    altNucleotides <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  } else if(refNuc == "D"){
    altNucleotides <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  } else {
    stop("Freq. Estimation: Unknown character used as Ref. nucleotide")
  }
  #mark the non-existing nuc/allele as not being alternate
  for (i in 1:6) {
    if( altNucleotides[i] == TRUE & observed[i] > 0.0 ){
      altNucleotides[i] = TRUE
    } else{
      altNucleotides[i] = FALSE
    }
  }
  
  return(altNucleotides)
}

MarkInconsistency <- function(chi1, degFreedom1, inconsistencyMark1,
                              chi2, degFreedom2, inconsistencyMark2, 
                              markThreshold) {
  
  mark <- ""
  
  pval1 <- 1 - pchisq( q = chi1, df = degFreedom1 ) 
  pval2 <- 1 - pchisq( q = chi2, df = degFreedom2 ) 
  
  if( pval1 < markThreshold) {
    mark <- inconsistencyMark1
  } else{
    #do nothing because the difference is not statistically significant
  }
  
  if( pval2 < markThreshold) {
    mark <- inconsistencyMark2
  } else{
    #do nothing because the difference is not statistically significant
  }
  
  return(mark)
}

GetAllelesLabel <- function(nucPosition, refNucleotide, As, Cs, Gs, Ts, Is, Ds) {
  
  #dbug
  # i <- 14
  # refNucleotide <- ezChiResults$As[i]
  # As = ezChiResults$As[i]
  # Cs = ezChiResults$Cs[i]
  # Gs = ezChiResults$Gs[i]
  # Ts = ezChiResults$Ts[i]
  # Is = ezChiResults$Is[i]
  # Ds = ezChiResults$Ds[i]
  # #observed <-  c(125, 0, 0, 104, 0, 0)
  # #dbug
  
  observed <- c(As, Cs, Gs, Ts, Is, Ds)
  
  refNucChar <- " "

  allelesCharacter <- c("A", "C", "G", "T", "I", "D")
  
  altNucleotides <- matrix(data = FALSE, nrow = 1, ncol = 6)
  colnames(altNucleotides) <- allelesCharacter
  
  if(refNucleotide == 1){
    refNucChar <- "A"
    altNucleotides <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  } else if(refNucleotide == 2){
    refNucChar <- "C"
    altNucleotides <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  } else if(refNucleotide == 3){
    refNucChar <- "G"
    altNucleotides <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
  } else if(refNucleotide == 4){
    refNucChar <- "T"
    altNucleotides <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  } else if(refNucleotide == 5){
    refNucChar <- "I"
    altNucleotides <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  } else if(refNucleotide == 6){
    refNucChar <- "D"
    altNucleotides <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  } else {
    stop("Freq. Estimation: Unknown character used as Ref. nucleotide")
  }
  
  for (i in 1:6) {
    
    if( altNucleotides[i] == TRUE & observed[i] > 0.0 ){
      altNucleotides[i] = TRUE
    } else{
      altNucleotides[i] = FALSE
    }
  }
  
  temp1 <- paste(allelesCharacter[altNucleotides], collapse = "")

  allelesLabel <- paste(refNucChar, temp1, collapse = "", sep = "")
  
  return(allelesLabel)

}

GetEasyChiEstimates <- function(polySite) {
  
  refNuc <- polySite$ref
  nucPosition <- polySite$pos
  #nAllelesPerRepeat(group1_repeat_1 (g1_1), g1_2, NA , g2_1, g2_2, NA, NA)
  nAllelesPerRepeat <- c( 25, 25, 0, 25, 25, 0, 0)
  
  #observed[g1_1, g1_2, NA , g2_1, g2_2, NA, NA), (a, c, g, t, i, d)]
  observed <- matrix( 0, nrow = 7, ncol = 6)
  
  observed[1, ] <- as.vector(polySite[ , 2:7], mode = "numeric")
  observed[2, ] <- as.vector(polySite[ , 8:13], mode = "numeric")
  observed[4, ] <- as.vector(polySite[ , 14:19], mode = "numeric")
  observed[5, ] <- as.vector(polySite[ , 20:25], mode = "numeric")
  
  #count nucleotides 
  As <- sum( observed[ , 1] )
  Cs <- sum( observed[ , 2] )
  Gs <- sum( observed[ , 3] )
  Ts <- sum( observed[ , 4] )
  Is <- sum( observed[ , 5] )
  Ds <- sum( observed[ , 6] )
  
  #do hetero calculations
  totalHeterozygosity  <- GetHetero( counts = rbind( observed[1, ],observed[2, ],observed[4, ], observed[5, ]), nGroups = 4 )
  group1Heterozygosity <- GetHetero( counts = rbind( observed[1, ],observed[2, ]), nGroups = 2 )
  group2Heterozygosity <- GetHetero( counts = rbind( observed[4, ],observed[5, ]), nGroups = 2 )

  #++++++++++++++++++++++++++++++++++++++++++++
  #Convert observed values based upon coverage to observed values based 
  #upon actual numbers of individuals analyzed in each of the four groups
  
  rowSum <- vector(mode = "numeric", length = 7)
  
  WOBS <- matrix(0, nrow = 7, ncol = 6)
  
  for (ir in 1:7) {
    for (ic in 1:6) {
      rowSum[ir] <- observed[ir, ic] + rowSum[ir]
    }
  }
  
  for (ir in 1:7) {
    for (ic in 1:6) {
      if( rowSum[ir] == 0) break
      TBL <- observed[ir, ic] / rowSum[ir]
      #normalize the number of alleles based on having 25 mosquitoes
      WOBS[ir, ic] <- TBL * nAllelesPerRepeat[ir] * 2
    }
  }
  
  observed[3, 1] = WOBS[1, 1] + WOBS[2, 1]
  observed[3, 2] = WOBS[1, 2] + WOBS[2, 2]
  observed[3, 3] = WOBS[1, 3] + WOBS[2, 3]
  observed[3, 4] = WOBS[1, 4] + WOBS[2, 4]
  observed[3, 5] = WOBS[1, 5] + WOBS[2, 5]
  observed[3, 6] = WOBS[1, 6] + WOBS[2, 6]
  
  observed[6, 1] = WOBS[4, 1] + WOBS[5, 1]
  observed[6, 2] = WOBS[4, 2] + WOBS[5, 2]
  observed[6, 3] = WOBS[4, 3] + WOBS[5, 3]
  observed[6, 4] = WOBS[4, 4] + WOBS[5, 4]
  observed[6, 5] = WOBS[4, 5] + WOBS[5, 5]
  observed[6, 6] = WOBS[4, 6] + WOBS[5, 6]
  
  observed[7, 1] = observed[3, 1] + observed[6, 1]
  observed[7, 2] = observed[3, 2] + observed[6, 2]
  observed[7, 3] = observed[3, 3] + observed[6, 3]
  observed[7, 4] = observed[3, 4] + observed[6, 4]
  observed[7, 5] = observed[3, 5] + observed[6, 5]
  observed[7, 6] = observed[3, 6] + observed[6, 6]
  
  #find alternate alleles and create label
  altAlleles <- GetAlternateAlleles(refNuc = refNuc, observed = observed[7, ])
  allelesCharacter <- c("A", "C", "G", "T", "I", "D")
  
  if(refNuc == "A") { 
    refNuc <- 1
  } else if (refNuc == "C") {
    refNuc <- 2 
  } else if (refNuc == "G") {
    refNuc <- 3 
  } else if (refNuc == "T") {
    refNuc <- 4
  } else {
    #do nothing
  }

#  LABEL <- paste(refNuc, allelesCharacter[altAlleles], sep = "")
  #end of alternate alleles
  
  remove = 0
  
  if(group1Heterozygosity != 0) { # if het > 0 
    #do chi for 
    myChiResults <- GetChi(
      nucPosition = nucPosition,
      observed = observed[3, ],
      wObs1 = WOBS[1, ],
      wObs2 = WOBS[2, ]) 
    
     group1ChiSqr <- myChiResults$chiSqr
     group1DegFreedom <- myChiResults$degFreedom

  } else {
    group1ChiSqr <- 0
    group1DegFreedom <- 0
  }
  
  if(group2Heterozygosity != 0) { # if het > 0 
    #do chi for 
    myChiResults <- GetChi(
      nucPosition = nucPosition,
      observed = observed[6, ],
      wObs1 = WOBS[4, ], 
      wObs2 = WOBS[5, ])
    
    group2ChiSqr <- myChiResults$chiSqr
    group2DegFreedom <- myChiResults$degFreedom

  } else {
    group2ChiSqr <- 0
    group2DegFreedom <- 0
  }  

  totalChiSqrResults <- GetChiAll( observed = observed, 
                                   total = sum( c(observed[3, ], observed[6, ]) ) ) 
  
  totalChiSqr <- totalChiSqrResults$TCHISQ3
  totalDegFreedom <- totalChiSqrResults$IDF3
  
  myFreq <- GetFrequencies(observed = observed, 
                           total =  sum( c(observed[3, ], sum( observed[6, ] )) ), 
                           altAlleles = altAlleles)
  
  group1AltAllFreq <- myFreq$RF1
  group2AltAllFreq <- myFreq$RF2
  
  PROB <- 0
  RLOD <- 0
  
  if(totalChiSqr < 30){
    PROB <- pchisq( q = totalChiSqr, df = totalDegFreedom )  
    RLOD <- -log10( 1 - PROB )
  } else if (totalChiSqr >= 30){
    if(totalDegFreedom == 1) RLOD = (0.219269476 * totalChiSqr) + 0.864404467
    if(totalDegFreedom == 2) RLOD = (0.217147241 * totalChiSqr) + 7.10543E-15
    if(totalDegFreedom == 3) RLOD = (0.215026228 * totalChiSqr) - 0.668488951
    if(totalDegFreedom == 4) RLOD = (0.212906488 * totalChiSqr) - 1.232280047
    if(totalDegFreedom == 5) RLOD = (0.210788076 * totalChiSqr) - 1.725110612
  }

  estimates <- c(nucPosition = nucPosition,
                 refNuc = refNuc,
                 group1AltAllFreq = group1AltAllFreq, #frequency of alternate Allele in group 1
                 group2AltAllFreq = group2AltAllFreq,
                 lod = RLOD,
                 group1Heteroz = group1Heterozygosity,
                 group2Heteroz = group2Heterozygosity,
                 totalHeteroz = totalHeterozygosity,
                 As = As,
                 Cs = Cs,
                 Gs = Gs,
                 Ts = Ts,
                 Is = Is,
                 Ds = Ds,
                 group1ChiSqr = group1ChiSqr,
                 group2ChiSqr = group2ChiSqr,
                 totalChiSqr = totalChiSqr,
                 group1DegFreedom = group1DegFreedom,
                 group2DegFreedom = group2DegFreedom,
                 totalDegFreedom = totalDegFreedom)
  
  return(estimates)
  
}

###Benjamini Hochberg to get the cut off for significant SNPs
###11/14/19 Original algorithm William C. Black IV
## modified as a function by lozano.saul@gmail.com 2020-06-22
GetBenjaminiHochberThreshold <- function(probabilities) {
  
  probabilities <- ezChiResults$totalProb
  
  n <- length(probabilities)
  sortedPvalues <- sort(probabilities)
  
  #create vector from 1 to n of probabilities vector increasing by 0.01/length
  jAlpha <- (1:n) * ( 0.01 / n )
  
  differences <- sortedPvalues - jAlpha
  
  negativeDifferences <- differences[ differences < 0 ]
  positiveDifferences <- negativeDifferences[ length(negativeDifferences) ]
  index <- differences == positiveDifferences
  
  benjaminiHochberThreshold <- sortedPvalues[index]
  
  return(benjaminiHochberThreshold)

}

