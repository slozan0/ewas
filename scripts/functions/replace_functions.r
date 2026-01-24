library(dplyr)

# 180470 pol "+"
# 477822 pol "-"

get_substitutions <- function(ref_codon, ref_amino, ref_nuc,
                              var_position, variants, polarity) {
  require(dplyr)
  # get ref aminoacid
  n_alleles <- nchar(variants)

  variants_table <- data.frame(
    nuc = character(n_alleles),
    codon = character(n_alleles),
    amino = character(n_alleles),
    reference = integer(n_alleles),
    replacement = character(n_alleles)
  )

  variants_table$nuc <- unlist(strsplit(variants, ""))

  # check if there are indels, if there are then do something

  if (polarity == "+") {
    variants_table <- variants_table |>
      mutate(nuc = case_when(
        nuc == "I" ~ "X",
        nuc == "D" ~ "X",
        # Keep the value unchanged if it doesn't match any of the conditions
        TRUE ~ nuc
      ))
  } else if (polarity == "-") {
    # Use mutate to change values in place
    variants_table <- variants_table |>
      mutate(nuc = case_when(
        nuc == "A" ~ "T",
        nuc == "C" ~ "G",
        nuc == "G" ~ "C",
        nuc == "T" ~ "A",
        nuc == "I" ~ "X",
        nuc == "D" ~ "X",
        # Keep the value unchanged if it doesn't match any of the conditions
        TRUE ~ nuc
      ))
  } else {
    print("something is wrong, you should be reading this!")
  }

  # get variant aminoacids
  for (i in 1:n_alleles) {
    variants_table$codon[i] <- make_variant_codon(
      ref_codon,
      variants_table$nuc[i],
      var_position
    )

    if (grepl("X", variants_table$codon[i])) {
      variants_table$amino[i] <- "UNK"
    } else {
      variants_table$amino[i] <- get_amino(variants_table$codon[i])
    }
  }

  # mark reference
  variants_table$reference <- ifelse(variants_table$codon == ref_codon, 1, 0)

  for (i in 1:n_alleles) {
    if (variants_table$amino[i] == ref_amino) {
      variants_table$replacement[i] <- "syn"
    } else if (variants_table$amino[i] == "UNK") {
      variants_table$replacement[i] <- "unk"
    } else {
      variants_table$replacement[i] <- "rep"
    }
  }

  variants_table
}

make_variant_codon <- function(ref_cod, new_nuc, positions_to_change) {
  # Check if positionsToChange is a vector
  if (!is.vector(positions_to_change)) {
    stop("positionsToChange must be a vector of integers.")
  }

  # Check if positions are within the reference codon length
  if (any(positions_to_change < 1 | positions_to_change > nchar(ref_cod))) {
    stop("positionsToChange must be within the reference codon length.")
  }

  # Sort positions in increasing order
  positions_to_change <- sort(positions_to_change)

  # Split the reference codon into substrings
  codon_parts <- strsplit(ref_cod, split = "")[[1]]

  # Replace nucleotides at specified positions
  for (i in positions_to_change) {
    codon_parts[[i]] <- new_nuc
  }

  # Combine the parts back into a single string
  new_codon <- paste(codon_parts, collapse = "")

  # Return the updated codon
  new_codon
}

get_amino <- function(codon) {
  codon <- toupper(codon) # Convert to uppercase for case-insensitive matching
  lookup_table <- c(
    "AAA" = "K", "AAG" = "K", "AAT" = "N", "AAC" = "N",
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
    "TTA" = "L", "TTG" = "L", "TTT" = "F", "TTC" = "F"
  )

  result <- lookup_table[codon]
  result
}

report_triple_variant <- function(variants_table) {
  require(dplyr)
  # if any of the variants are replacement
  # pick that, if not, pick the second

  # count the number of replacements
  # count the number of unks
  # count the number of syn

  counts <- table(variants_table$replacement)
  # str\(counts["rep"]\)

  variants_table$nRep <- NA

  if ("rep" %in% variants_table$replacement) {
    # count how many replacements
    # make the second row the replacement
    replacement_row_index <- which(variants_table$replacement == "rep")[1]

    df_reordered <- variants_table |>
      # Remove the row to be moved
      slice(-replacement_row_index) |>
      # Add it as the second row
      add_row(slice(variants_table, replacement_row_index), .before = 2)

    df_reordered$nRep[2] <- counts["rep"]

    df_reordered
  } else {
    print(variants_table)
    variants_table
  }

  # make notation on the replacement type
  # + single syn
  # ++ single syn
  #* single replacement
  #** double replacement
}
