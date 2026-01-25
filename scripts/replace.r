rm(list = ls())

library(dplyr)
source("replacer_functions.r")

chi_results <- readRDS("data/output/annotated_5feb_chrom2.rds")

replace_input <- chi_results |>
  dplyr::select(
    nucPosition, alleles,
    group1AltAllFreq,
    group2AltAllFreq,
    lod,
    totalHeteroz,
    group1Heteroz,
    group2Heteroz,
    info
  )

replace_input$info <- gsub("\\s+", " ", replace_input$info)

replace_input <- replace_input |>
  tidyr::separate(info,
    into = c("LOC", "POL", "REFN",
             "GTYPE", "CODPOS", "RESID",
             "CODON1", "AA1", "CODON2", "AA2"),
    sep = " ",
    extra = "merge"
  )

replace_input$info <- NULL
replace_input$REPLACEMENT <- NA # make space in ram for new data
replace_input$NREPS <- NA # make space in ram for new data

# replaceInput1000 <- head\(replaceInput, n = 15000\)
# nRows <- nrow\(replaceInput1000\)

# \i <- 766 # one indel
# \i <- 2735 # triple variant everyting syn
# \i <- 204 # triple variant double rep

for (i in seq_len(nrow(replace_input))) {
  if (replace_input$GTYPE[i] == "CDS" && !is.na(replace_input$GTYPE[i])) {
    nuc_pos <- replace_input$nucPosition[i]

    # print\( paste("CDS mutation in nuc pos:", nucPos) \)

    ref_codon <- replace_input$CODON1[i]
    ref_amino <- replace_input$AA1[i]
    ref_nuc <- replace_input$REFN[i]

    var_position <- as.numeric(replace_input$CODPOS[i])
    variants <- replace_input$alleles[i]
    polarity <- replace_input$POL[i]

    variants_table <- get_substitutions(
      ref_codon = ref_codon,
      ref_amino = ref_amino,
      ref_nuc = ref_nuc,
      var_position = var_position,
      variants = variants,
      polarity = polarity
    )

    if (nrow(variants_table) > 2) { # if more than 2 alleles

      temp <- report_triple_variant(variants_table)
      replace_input$CODON2[i] <- temp$codon[2]
      replace_input$AA2[i] <- temp$amino[2]
      replace_input$REPLACEMENT[i] <- temp$replacement[2]
      replace_input$NREPS[i] <- temp$nRep[2]
    } else {
      replace_input$CODON2[i] <- variants_table$codon[2]
      replace_input$AA2[i] <- variants_table$amino[2]
      replace_input$REPLACEMENT[i] <- variants_table$replacement[2]
    }
  } else {
    # variation not in CDS; do nothing.
  }
}
