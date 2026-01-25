# Fix d1
data_d1 <- readRDS("data/output/5feb_tem_d1_c3_rc.rds")
data_d1$position <- as.integer(data_d1$position)
data_d1 <- data_d1[!is.na(data_d1$position), ]
saveRDS(data_d1, "data/output/5feb_tem_d1_c3_rc.rds")

# Fix d2
data_d2 <- readRDS("data/output/5feb_tem_d2_c3_rc.rds")
data_d2$position <- as.integer(data_d2$position)
data_d2 <- data_d2[!is.na(data_d2$position), ]
saveRDS(data_d2, "data/output/5feb_tem_d2_c3_rc.rds")

unique(data_d2$snp11)
