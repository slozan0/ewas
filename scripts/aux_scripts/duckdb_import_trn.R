rm(list = ls())

library(duckdb)
library(dbplyr)

con <- dbConnect(duckdb("chrom3_4.duckdb")) # new DB on disk

# One-shot bulk import ──────────────────────────────────────────────────────
dbExecute(con, "
  PRAGMA memory_limit='20GB';              -- optional safety guard
  CREATE OR REPLACE TABLE chrom_data AS
  SELECT
      substr(line, 1, 10)  AS identifier, -- characters 1-10
      substr(line, 12)     AS info        -- from 12 to end
  FROM read_csv(
        'data/input/annotations/Chrom3_4.trn',
        delim      = '\n',                -- newline = record separator
        header     = false,
        quote      = '',                  -- ignore quotes
        escape     = '',                  -- ignore escapes
        columns    = {'line':'VARCHAR'},  -- exactly one column called line
        auto_detect = false,              -- skip sampling
        parallel = true,
        ignore_errors = true,       -- <-- skip bad lines
        rejects_table  = 'bad_rows' -- keep the rejected rows for inspection
  );
")

# check if contents of bad formatted lines
dbGetQuery(con, "SELECT * FROM bad_rows LIMIT 5")
# Quick check
dbGetQuery(con, "SELECT * FROM chrom_data LIMIT 10")

# Quick duplicate check before primary-key
# If duplicates exist the primary-key build will fail anyway,
# so it’s smart to look once.
dbGetQuery(con, "
  SELECT identifier, COUNT(*) AS n
  FROM chrom_data
  GROUP BY identifier HAVING n > 1
  LIMIT 10")

# Check the data types
dbGetQuery(con, "PRAGMA table_info('chrom_data');")

## change identifier VARCHAR → BIGINT in place ───────────────────────────
dbExecute(con, "
  ALTER TABLE chrom_data
  ALTER COLUMN identifier
  SET DATA TYPE BIGINT                 -- or INTEGER if you’re sure < 2.1 b
  USING CAST(identifier AS BIGINT);    -- single streaming pass
")

# dbGetQuery(con, "SELECT *
#                  FROM chrom_data
#                  WHERE identifier = 409333263 ")

## remove dupes from db
# dbGetQuery(con, "
# DELETE FROM chrom_data
# USING (
#   SELECT rowid,
#   ROW_NUMBER() OVER (PARTITION BY identifier ORDER BY rowid) AS rn
#   FROM chrom_data
# ) dup
# WHERE chrom_data.rowid = dup.rowid
# AND dup.rn > 1;")

# Add primary key
dbExecute(con, "
  ALTER TABLE chrom_data
  ADD CONSTRAINT chrom_pk PRIMARY KEY(identifier);
")

dbDisconnect(con, shutdown = TRUE)
