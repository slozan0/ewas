library(dplyr)

#180470 pol "+"
#477822 pol "-"

GetSubstitutions <- function(refCodon, refAmino, refNuc,
                             varPosition, variants, 
                             polarity) {
  require(dplyr)
  #get ref aminoacid
  nAlleles <- nchar(variants)
  
  variantsTable <- data.frame(
    nuc = character(nAlleles),
    codon = character(nAlleles), 
    amino = character(nAlleles),
    reference = integer(nAlleles),
    replacement = character(nAlleles)
  )
  
  variantsTable$nuc <- unlist(strsplit(variants, ""))
  
  #check if there are indels, if there are then do something

  if(polarity == "+"){
    variantsTable <- variantsTable %>%
      mutate(nuc = case_when(
        nuc == "I" ~ "X",
        nuc == "D" ~ "X",
        TRUE ~ nuc  # Keep the value unchanged if it doesn't match any of the conditions
      ))
    
  } else if(polarity == "-"){
    # Use mutate to change values in place
    variantsTable <- variantsTable %>%
      mutate(nuc = case_when(
        nuc == "A" ~ "T",
        nuc == "C" ~ "G",
        nuc == "G" ~ "C",
        nuc == "T" ~ "A",
        nuc == "I" ~ "X",
        nuc == "D" ~ "X",
        TRUE ~ nuc  # Keep the value unchanged if it doesn't match any of the conditions
      ))
    
  } else {
    print("something is wrong, you should be reading this!")
  }
  
  #get variant aminoacids
  for (i in 1:nAlleles) {
    variantsTable$codon[i] <- MakeVariantCodon(refCodon, 
                                               variantsTable$nuc[i], 
                                               varPosition)
    
    if( grepl("X", variantsTable$codon[i]) ){ 
    variantsTable$amino[i] <- "UNK"
    } else { 
      variantsTable$amino[i] <- GetAmino(variantsTable$codon[i])
    }
  
  }
  
  # mark reference  
  variantsTable$reference <- ifelse(variantsTable$codon == refCodon, 1, 0)
  
  for (i in 1:nAlleles) {
    
    if(variantsTable$amino[i] == refAmino) {
      variantsTable$replacement[i] <- "syn"
    }else if(variantsTable$amino[i] == "UNK" ){
      variantsTable$replacement[i] <- "unk"
    } else {
      variantsTable$replacement[i] <- "rep"
    }
  }
  
  return(variantsTable)
  
}

MakeVariantCodon <- function(refCod, newNuc, positionsToChange) {
  # Check if positionsToChange is a vector
  if (!is.vector(positionsToChange)) {
    stop("positionsToChange must be a vector of integers.")
  }
  
  # Check if positions are within the reference codon length
  if (any(positionsToChange < 1 | positionsToChange > nchar(refCod))) {
    stop("positionsToChange must be within the reference codon length.")
  }
  
  # Sort positions in increasing order
  positionsToChange <- sort(positionsToChange)
  
  # Split the reference codon into substrings
  codonParts <- strsplit(refCod, split = "")[[1]]
  
  # Replace nucleotides at specified positions
  for (i in positionsToChange) {
    codonParts[[i]] <- newNuc
  }
  
  # Combine the parts back into a single string
  newCodon <- paste(codonParts, collapse = "")
  
  # Return the updated codon
  return(newCodon)
}

GetAmino <- function(codon) {
  codon <- toupper(codon)  # Convert to uppercase for case-insensitive matching
  lookup_table <- c("AAA" = "K", "AAG" = "K", "AAT" = "N", "AAC" = "N",
                    "AGA" = "R", "AGG" = "R", "AGT" = "S", "AGC" = "S",
                    "ACA" = "T", "ACG" = "T", "ACT" = "T", "ACC" = "T",
                    "ATA" = "I", "ATG" = "M", "ATT" = "I", "ATC" = "I",
                    "GAA" = "E", "GAG" = "E", "GAT" = "D", "GAC" = "D",
                    "GGA" = "G", "GGG" = "G", "GGT" = "G", "GGC" = "G",
                    "GCA" = "A", "GCG" = "A", "GCT" = "A", "GCC" = "A",
                    "GTA" = "V", "GTG" = "V", "GTT" = "V", "GTC" = "V",
                    "CAA" = "Q", "CAG" = "Q", "CAT" = "H", "CAC" = "H",
                    "CGA" = "R", "CGG" = "R", "CGT" = "R", "CGC" = "R",
                    "CCA" = "P", "CCG" = "P", "CCT" = "P", "CCC" = "P",
                    "CTA" = "L", "CTG" = "L", "CTT" = "L", "CTC" = "L",
                    "TAA" = "*", "TAG" = "*", "TAT" = "Y", "TAC" = "Y",
                    "TGA" = "*", "TGG" = "W", "TGT" = "C", "TGC" = "C",
                    "TCA" = "S", "TCG" = "S", "TCT" = "S", "TCC" = "S",
                    "TTA" = "L", "TTG" = "L", "TTT" = "F", "TTC" = "F")
  
  result <- lookup_table[codon]
  return(result)
}

ReportTripleVariant <- function(variantsTable) {
  require(dplyr)
  #if any of the variants are replacement
  # pick that if not, pick the second
  
  #count the number of replacements
  #count the number of unks 
  #count the number of syn 
  
  counts <- table(variantsTable$replacement)
  # str(counts["rep"])
  
  variantsTable$nRep <- NA
  
  if("rep" %in% variantsTable$replacement) {
    #count how many replacements
    #make the second row the replacement
    replacement_row_index <- which(variantsTable$replacement == "rep")[1]
    
    dfReordered <- variantsTable %>%
      # Remove the row to be moved
      slice(-replacement_row_index) %>%
      # Add it as the second row
      add_row(slice(variantsTable, replacement_row_index), .before = 2)
    
    dfReordered$nRep[2] <- counts["rep"]
    
    return(dfReordered)
  } else {
    print(variantsTable)
    return(variantsTable)
  }

  
  #make notation on the replacement type 
  
  # + single syn
  # ++ single syn
  
  #* single replacement
  #** double replacement
  
  
}


