# setup environment
rm(list = ls())
library(data.table)

# set i/o ----
## input files
phen1_1 <- readRDS("data/output/5feb_tem_a1_c3_sample.rds")
phen1_2 <- readRDS("data/output/5feb_tem_a2_c3_sample.rds")
phen2_1 <- readRDS("data/output/5feb_tem_d1_c3_sample.rds")
phen2_2 <- readRDS("data/output/5feb_tem_d2_c3_sample.rds")

## output file names
rds_file <- "data/output/5feb_tem_chr3_avd_sample.rds"
text_file <- "data/output/5feb_tem_chr3_avd_sample.tbl"

## input file labels
phen1_1_name <- "5feb_tem_a1"
phen1_2_name <- "5feb_tem_a2"
phen2_1_name <- "5feb_tem_d1"
phen2_2_name <- "5feb_tem_d2"

# program starts here ----
p1_1xp1_2 <- merge(x = phen1_1, y = phen1_2,
                   by = "position", suffixes = c("1", "2"))
p2_1xp2_2 <- merge(x = phen2_1, y = phen2_2,
                   by = "position", suffixes = c("3", "4"))

# venn2x2X2
p1_v_p2 <- merge(x = p1_1xp1_2, y = p2_1xp2_2, by = "position")

rm(p1_1xp1_2)
rm(p2_1xp2_2)

pol <- p1_v_p2[rowSums(p1_v_p2 == 0) <= 19, ] # find the polymorphic sites
# monomorphic position
# pos   chrom   A  C  G  T  I  D
# 5995      1  64  0  0  0  0  0
# 5995      1  64  0  0  0  0  0
# 5995      1  64  0  0  0  0  0
# 5995      1  64  0  0  0  0  0

# polymorphic position
# pos   chrom   A  C  G  T  I  D
# 5995      1  64  0  0  0  0  0
# 5995      1  63  1  0  0  0  0
# 5995      1  64  0  0  0  0  0
# 5995      1  64  0  0  0  0  0

# find in the row all the zeros and give them true or false (p1Vp2 == 0)
# if the number of TRUEs is 20 then site is the row monomorphic

saveRDS(pol, file = rds_file)

n_rows <- nrow(pol)

file_conn <- file(text_file, open = "wt")

# Write each row to the file
for (i in 1:n_rows) {
  line1 <- sprintf("%10i", pol$pos[i])
  line2 <-
    sprintf("%s%20i%10i%10i%10i%10i%10i",
            phen1_1_name,
            pol$a1[i], pol$c1[i], pol$g1[i], pol$t1[i], pol$i1[i], pol$d1[i])
  line3 <-
    sprintf("%s%20i%10i%10i%10i%10i%10i",
            phen1_2_name,
            pol$a2[i], pol$c2[i], pol$g2[i], pol$t2[i], pol$i2[i], pol$d2[i])
  line4 <-
    sprintf("%s%20i%10i%10i%10i%10i%10i",
            phen2_1_name,
            pol$a3[i], pol$c3[i], pol$g3[i], pol$t3[i], pol$i3[i], pol$d3[i])
  line5 <-
    sprintf("%s%20i%10i%10i%10i%10i%10i",
            phen2_2_name,
            pol$a4[i], pol$c4[i], pol$g4[i], pol$t4[i], pol$i4[i], pol$d4[i])

  writeLines(c(line1, line2, line3, line4, line5), file_conn)
}

# Close the file connection
close(file_conn)
