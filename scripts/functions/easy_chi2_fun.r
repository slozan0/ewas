get_hetero <- function(counts, n_groups) {
  a <- sum(counts[, 1])
  c <- sum(counts[, 2])
  g <- sum(counts[, 3])
  t <- sum(counts[, 4])
  i <- sum(counts[, 5])
  d <- sum(counts[, 6])

  sum_nuc <- a + c + g + t + i + d

  freq_a <- a / sum_nuc
  freq_c <- c / sum_nuc
  freq_g <- g / sum_nuc
  freq_t <- t / sum_nuc
  freq_i <- i / sum_nuc
  freq_d <- d / sum_nuc

  homozygocity <-
    (freq_a * freq_a) + (freq_c * freq_c) +
    (freq_g * freq_g) + (freq_t * freq_t) +
    (freq_i * freq_i) + (freq_d * freq_d)

  heterozygosity <- 1 - homozygocity

  heterozygosity
}

get_frequencies <- function(observed, total, alt_alleles) {
  # nolint start:
  # observed <- ROBS
  # total <- T1 + T2
  # altAlleles <- altAlleles
  # nolint end:

  n_alleles_g1 <- 0
  n_alleles_g2 <- 0

  freq_group_one_temp <- observed[3, ] / sum(observed[3, ])
  freq_group_two_temp <- observed[6, ] / sum(observed[6, ])

  # number of alleles in nuc position
  for (i in 1:6) {
    if (freq_group_one_temp[i] > 0.0) {
      n_alleles_g1 <- n_alleles_g1 + 1
    }
    if (freq_group_two_temp[i] > 0.0) {
      n_alleles_g2 <- n_alleles_g2 + 1
    }
  }

  warning_msg <-
    "Freq. Estimation: The number of alleles is not equal between the groups."
  if (n_alleles_g1 != n_alleles_g2) warning(warning_msg)

  freq_group_one <- sum(freq_group_one_temp[alt_alleles])
  freq_group_two <- sum(freq_group_two_temp[alt_alleles])

  results <- list(
    RF1 = freq_group_one,
    RF2 = freq_group_two
  )

  results
}

get_chi <- function(nuc_position, observed, w_obs1, w_obs2) {
  # nolint start:
  # rObs <- ROBS[3, ]
  # wObs1 = WOBS[1, ]
  # wObs2 = WOBS[2, ]
  # nolint end:

  n_alleles <- 0

  expected_values_group_one <- vector(mode = "numeric", length = 7)
  expected_values_group_two <- vector(mode = "numeric", length = 7)

  chi_group_one <- vector(mode = "numeric", length = 7)
  chi_group_two <- vector(mode = "numeric", length = 7)

  column_proportion <- vector(mode = "numeric", length = 7)

  for (i in 1:6) {
    if (observed[i] > 0) n_alleles <- n_alleles + 1
  }

  total <- sum(observed)

  column_proportion <- observed / total

  freq_group_one <- sum(w_obs1) / total
  freq_group_two <- sum(w_obs2) / total

  expected_values_group_one <- freq_group_one * column_proportion * total
  expected_values_group_two <- freq_group_two * column_proportion * total

  chi_group_one <- ((w_obs1 - expected_values_group_one)^2) /
    expected_values_group_one
  chi_group_two <- ((w_obs2 - expected_values_group_two)^2) /
    expected_values_group_two

  chi_sqr <- sum(chi_group_one, chi_group_two, na.rm = TRUE)
  deg_freedom <- n_alleles - 1

  results <- list(
    chiSqr = chi_sqr,
    degFreedom = deg_freedom
  )

  results
}

get_chi_all <- function(observed, total) {
  # nolint start:
  # observed <- ROBS
  # total <- T1 + T2
  # nolint end:

  freq_total <- vector(mode = "numeric", length = 7)
  expected_value_group1 <- vector(mode = "numeric", length = 7)
  expected_value_group2 <- vector(mode = "numeric", length = 7)
  chi_group1 <- vector(mode = "numeric", length = 7)
  chi_group2 <- vector(mode = "numeric", length = 7)

  n_alleles <- 0

  if (observed[7, 1] > 0) n_alleles <- n_alleles + 1
  if (observed[7, 2] > 0) n_alleles <- n_alleles + 1
  if (observed[7, 3] > 0) n_alleles <- n_alleles + 1
  if (observed[7, 4] > 0) n_alleles <- n_alleles + 1
  if (observed[7, 5] > 0) n_alleles <- n_alleles + 1
  if (observed[7, 6] > 0) n_alleles <- n_alleles + 1

  freq_total <- observed[7, ] / total

  freq_group_one <- sum(observed[3, ]) / total
  freq_group_two <- sum(observed[6, ]) / total

  expected_value_group1 <- freq_group_one * freq_total * total
  expected_value_group2 <- freq_group_two * freq_total * total

  chi_group1 <- ((observed[3, ] - expected_value_group1)^2) /
    expected_value_group1
  chi_group2 <- ((observed[6, ] - expected_value_group2)^2) /
    expected_value_group2

  chi_all <- sum(chi_group1 + chi_group2, na.rm = TRUE)

  deg_freedom <- n_alleles - 1

  results <- list(
    TCHISQ3 = chi_all,
    IDF3 = deg_freedom
  )

  results
}

get_alternate_alleles <- function(ref_nuc, observed) {
  # nolint start:
  # refNuc <- refNuc
  # observed <- ROBS[7,]
  # nolint end:

  alt_nucleotides <- matrix(data = FALSE, nrow = 1, ncol = 6)
  colnames(alt_nucleotides) <- c("A", "C", "G", "T", "I", "D")

  # mark the reference nuc/allele as not being the alternate
  if (ref_nuc == "A") {
    alt_nucleotides <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  } else if (ref_nuc == "C") {
    alt_nucleotides <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  } else if (ref_nuc == "G") {
    alt_nucleotides <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
  } else if (ref_nuc == "T") {
    alt_nucleotides <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  } else if (ref_nuc == "I") {
    alt_nucleotides <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  } else if (ref_nuc == "D") {
    alt_nucleotides <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  } else {
    stop("Freq. Estimation: Unknown character used as Ref. nucleotide")
  }
  # mark the non-existing nuc/allele as not being alternate
  for (i in 1:6) {
    if (alt_nucleotides[i] == TRUE && observed[i] > 0.0) {
      alt_nucleotides[i] <- TRUE
    } else {
      alt_nucleotides[i] <- FALSE
    }
  }

  alt_nucleotides
}

mark_inconsistency <- function(chi1, deg_freedom1, inconsistency_mark1,
                               chi2, deg_freedom2, inconsistency_mark2,
                               mark_threshold) {
  mark <- ""

  pval1 <- 1 - pchisq(q = chi1, df = deg_freedom1)
  pval2 <- 1 - pchisq(q = chi2, df = deg_freedom2)

  if (pval1 < mark_threshold) {
    mark <- inconsistency_mark1
  } else {
    # do nothing because the difference is not statistically significant
  }

  if (pval2 < mark_threshold) {
    mark <- inconsistency_mark2
  } else {
    # do nothing because the difference is not statistically significant
  }

  mark
}

get_alleles_label <- function(nuc_position, ref_nucleotide,
                              a_s, c_s, g_s, t_s, i_s, d_s) {
  # nolint start:
  # i <- 14
  # refNucleotide <- ezChiResults$As[i]
  # As = ezChiResults$As[i]
  # Cs = ezChiResults$Cs[i]
  # Gs = ezChiResults$Gs[i]
  # Ts = ezChiResults$Ts[i]
  # Is = ezChiResults$Is[i]
  # Ds = ezChiResults$Ds[i]
  # #observed <-  c(125, 0, 0, 104, 0, 0)
  # nolint end:

  observed <- c(a_s, c_s, g_s, t_s, i_s, d_s)

  ref_nuc_char <- " "

  alleles_character <- c("A", "C", "G", "T", "I", "D")

  alt_nucleotides <- matrix(data = FALSE, nrow = 1, ncol = 6)
  colnames(alt_nucleotides) <- alleles_character

  if (ref_nucleotide == 1) {
    ref_nuc_char <- "A"
    alt_nucleotides <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  } else if (ref_nucleotide == 2) {
    ref_nuc_char <- "C"
    alt_nucleotides <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  } else if (ref_nucleotide == 3) {
    ref_nuc_char <- "G"
    alt_nucleotides <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
  } else if (ref_nucleotide == 4) {
    ref_nuc_char <- "T"
    alt_nucleotides <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  } else if (ref_nucleotide == 5) {
    ref_nuc_char <- "I"
    alt_nucleotides <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  } else if (ref_nucleotide == 6) {
    ref_nuc_char <- "D"
    alt_nucleotides <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  } else {
    stop("Freq. Estimation: Unknown character used as Ref. nucleotide")
  }

  for (i in 1:6) {
    if (alt_nucleotides[i] == TRUE && observed[i] > 0.0) {
      alt_nucleotides[i] <- TRUE
    } else {
      alt_nucleotides[i] <- FALSE
    }
  }

  temp1 <- paste(alleles_character[alt_nucleotides], collapse = "")

  alleles_label <- paste(ref_nuc_char, temp1, collapse = "", sep = "")

  alleles_label
}

get_easy_chi_estimates <- function(poly_site) {
  ref_nuc <- poly_site$refnuc
  nuc_position <- poly_site$position
  # nAllelesPerRepeat \( group1_repeat_1 (g1_1), g1_2, NA , g2_1, g2_2, NA, NA\)
  # Number of mosquitoes per replicate (25 per group)
  n_alleles_per_repeat <- c(25, 25, 0, 25, 25, 0, 0)

  # observed\(g1_1, g1_2, NA , g2_1, g2_2, NA, NA), (a, c, g, t, i, d)]
  observed <- matrix(0, nrow = 7, ncol = 6)

  observed[1, ] <- as.vector(poly_site[, 2:7], mode = "numeric")
  observed[2, ] <- as.vector(poly_site[, 8:13], mode = "numeric")
  observed[4, ] <- as.vector(poly_site[, 14:19], mode = "numeric")
  observed[5, ] <- as.vector(poly_site[, 20:25], mode = "numeric")

  # count nucleotides
  a_s <- sum(observed[, 1])
  c_s <- sum(observed[, 2])
  g_s <- sum(observed[, 3])
  t_s <- sum(observed[, 4])
  i_s <- sum(observed[, 5])
  d_s <- sum(observed[, 6])

  # do hetero calculations
  total_heterozygosity <-
    get_hetero(counts = rbind(observed[1, ],
                              observed[2, ],
                              observed[4, ],
                              observed[5, ]),
               n_groups = 4)

  group1_heterozygosity <-
    get_hetero(counts = rbind(observed[1, ], observed[2, ]), n_groups = 2)
  group2_heterozygosity <-
    get_hetero(counts = rbind(observed[4, ], observed[5, ]), n_groups = 2)

  #++++++++++++++++++++++++++++++++++++++++++++
  # Convert observed values based upon coverage to observed values based
  # upon actual numbers of individuals analyzed in each of the four groups

  row_sum <- vector(mode = "numeric", length = 7)

  w_obs <- matrix(0, nrow = 7, ncol = 6)

  for (ir in 1:7) {
    for (ic in 1:6) {
      row_sum[ir] <- observed[ir, ic] + row_sum[ir]
    }
  }

  for (ir in 1:7) {
    for (ic in 1:6) {
      if (row_sum[ir] == 0) break
      tbl <- observed[ir, ic] / row_sum[ir]
      # normalize the number of alleles based on having 25 mosquitoes
      w_obs[ir, ic] <- tbl * n_alleles_per_repeat[ir] * 2
    }
  }

  observed[3, 1] <- w_obs[1, 1] + w_obs[2, 1]
  observed[3, 2] <- w_obs[1, 2] + w_obs[2, 2]
  observed[3, 3] <- w_obs[1, 3] + w_obs[2, 3]
  observed[3, 4] <- w_obs[1, 4] + w_obs[2, 4]
  observed[3, 5] <- w_obs[1, 5] + w_obs[2, 5]
  observed[3, 6] <- w_obs[1, 6] + w_obs[2, 6]

  observed[6, 1] <- w_obs[4, 1] + w_obs[5, 1]
  observed[6, 2] <- w_obs[4, 2] + w_obs[5, 2]
  observed[6, 3] <- w_obs[4, 3] + w_obs[5, 3]
  observed[6, 4] <- w_obs[4, 4] + w_obs[5, 4]
  observed[6, 5] <- w_obs[4, 5] + w_obs[5, 5]
  observed[6, 6] <- w_obs[4, 6] + w_obs[5, 6]

  observed[7, 1] <- observed[3, 1] + observed[6, 1]
  observed[7, 2] <- observed[3, 2] + observed[6, 2]
  observed[7, 3] <- observed[3, 3] + observed[6, 3]
  observed[7, 4] <- observed[3, 4] + observed[6, 4]
  observed[7, 5] <- observed[3, 5] + observed[6, 5]
  observed[7, 6] <- observed[3, 6] + observed[6, 6]

  # find alternate alleles and create label
  alt_alleles <-
    get_alternate_alleles(ref_nuc = ref_nuc, observed = observed[7, ])

  if (ref_nuc == "A") {
    ref_nuc <- 1
  } else if (ref_nuc == "C") {
    ref_nuc <- 2
  } else if (ref_nuc == "G") {
    ref_nuc <- 3
  } else if (ref_nuc == "T") {
    ref_nuc <- 4
  } else {
    # do nothing
  }

  #  LABEL equals paste(refNuc, allelesCharacter[altAlleles], sep = "")
  # end of alternate alleles

  if (group1_heterozygosity != 0) { # if het > 0
    # do chi for
    my_chi_results <- get_chi(
      nuc_position = nuc_position,
      observed = observed[3, ],
      w_obs1 = w_obs[1, ],
      w_obs2 = w_obs[2, ]
    )

    group1_chi_sqr <- my_chi_results$chiSqr
    group1_deg_freedom <- my_chi_results$degFreedom
  } else {
    group1_chi_sqr <- 0
    group1_deg_freedom <- 0
  }

  if (group2_heterozygosity != 0) { # if het > 0
    # do chi for
    my_chi_results <- get_chi(
      nuc_position = nuc_position,
      observed = observed[6, ],
      w_obs1 = w_obs[4, ],
      w_obs2 = w_obs[5, ]
    )

    group2_chi_sqr <- my_chi_results$chiSqr
    group2_deg_freedom <- my_chi_results$degFreedom
  } else {
    group2_chi_sqr <- 0
    group2_deg_freedom <- 0
  }

  total_chi_sqr_results <- get_chi_all(
    observed = observed,
    total = sum(c(observed[3, ], observed[6, ]))
  )

  total_chi_sqr <- total_chi_sqr_results$TCHISQ3
  total_deg_freedom <- total_chi_sqr_results$IDF3

  my_freq <- get_frequencies(
    observed = observed,
    total = sum(c(observed[3, ], sum(observed[6, ]))),
    alt_alleles = alt_alleles
  )

  group1_alt_all_freq <- my_freq$RF1
  group2_alt_all_freq <- my_freq$RF2

  prob <- 0
  rlod <- 0

  if (total_chi_sqr < 30) {
    prob <- pchisq(q = total_chi_sqr, df = total_deg_freedom)
    rlod <- -log10(1 - prob)
  } else if (total_chi_sqr >= 30) {
    if (total_deg_freedom == 1) rlod <-
      (0.219269476 * total_chi_sqr) + 0.864404467
    if (total_deg_freedom == 2) rlod <-
      (0.217147241 * total_chi_sqr) + 7.10543E-15
    if (total_deg_freedom == 3) rlod <-
      (0.215026228 * total_chi_sqr) - 0.668488951
    if (total_deg_freedom == 4) rlod <-
      (0.212906488 * total_chi_sqr) - 1.232280047
    if (total_deg_freedom == 5) rlod <-
      (0.210788076 * total_chi_sqr) - 1.725110612
  }

  estimates <- c(
    nucPosition = nuc_position,
    refNuc = ref_nuc,
    # frequency of alternate Allele in group 1
    group1AltAllFreq = group1_alt_all_freq,
    group2AltAllFreq = group2_alt_all_freq,
    lod = rlod,
    group1Heteroz = group1_heterozygosity,
    group2Heteroz = group2_heterozygosity,
    totalHeteroz = total_heterozygosity,
    As = a_s,
    Cs = c_s,
    Gs = g_s,
    Ts = t_s,
    Is = i_s,
    Ds = d_s,
    group1ChiSqr = group1_chi_sqr,
    group2ChiSqr = group2_chi_sqr,
    totalChiSqr = total_chi_sqr,
    group1DegFreedom = group1_deg_freedom,
    group2DegFreedom = group2_deg_freedom,
    totalDegFreedom = total_deg_freedom
  )

  estimates
}

### Benjamini Hochberg to get the cut off for significant SNPs
### 11/14/19 Original algorithm William C. Black IV
## modified as a function by lozano.saul@gmail.com 2020-06-22
get_benjamini_hochber_thresh <- function(probabilities) {
  n <- length(probabilities)
  sorted_p_values <- sort(probabilities)

  # create vector from 1 to n of probabilities vector increasing by 0.01/length
  j_alpha <- (1:n) * (0.01 / n)

  differences <- sorted_p_values - j_alpha

  negative_differences <- differences[differences < 0]
  positive_differences <- negative_differences[length(negative_differences)]
  index <- differences == positive_differences

  ben_hoc_threshold <- sorted_p_values[index]

  ben_hoc_threshold
}
