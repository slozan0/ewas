rm(list = ls())

library(data.table)
library(DBI)
library(duckdb)

# read in your results (should be small enough to fit in memory)
chiResults <- readRDS("data/input/5feb_tem_ezchi_c1.rds")

# path to your DuckDB file
dbPath <- "data/input/annotations/chrom1_3.duckdb"

# path to your results file
resultsPath <- "data/output/annotated_5feb_chrom1_3.rds"


# connect (read-only)
con <- dbConnect(duckdb::duckdb(), dbdir = dbPath, read_only = TRUE)

# 1) push just the lookup keys (nucPosition) into a temporary DuckDB table
#    (no matter how big chrom_data is, this only uploads your chiResults$nucPosition vector)
dbWriteTable(con,
             name = "tmp_ids",
             value = data.frame(nucPosition = chiResults$nucPosition),
             temporary = TRUE,
             overwrite = TRUE)

# 2) use a single SQL join to grab the matching info values
#    DuckDB will use the index on chrom_data.identifier (if one exists) 
#    and stream through the 400M rows efficiently
lookup <- dbGetQuery(con, "
  SELECT 
    t.nucPosition,
    c.info
  FROM tmp_ids AS t
  LEFT JOIN chrom_data AS c
    ON t.nucPosition = c.identifier
")

# 3) clean up
dbDisconnect(con, shutdown = TRUE)

# 4) merge the info back into the original data.table
setDT(lookup)
setkey(lookup, nucPosition)
chiResults[, info := lookup[.SD, info, on = "nucPosition"]]

# 5) save
saveRDS(chiResults, resultsPath)
