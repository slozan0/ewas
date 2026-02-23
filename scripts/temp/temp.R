# Fix d1
data_d1 <- readRDS("data/sample/5feb_tem_d1_c3_rc_sample.rds")
data_d1$position <- as.integer(data_d1$position)
data_d1 <- data_d1[!is.na(data_d1$position), ]

unique(data_d2$snp11)
