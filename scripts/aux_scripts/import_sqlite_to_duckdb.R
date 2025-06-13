rm(list = ls())

# Load packages ----
library(duckdb)
library(dbplyr)

# open / create the DuckDB file ----
con <- dbConnect(duckdb("chrom3_1.duckdb"))  # on-disk file

dbExecute(con, "SET memory_limit='25GB';      -- or e.g. '8GB'")

# Load duckdb sqlite extentions
dbExecute(con, "LOAD sqlite;")        

# Attach sqlite database
## Use normalized path otherwise it will not attach
sqlite_path <- normalizePath("data/input/annotations/chrom3_1.db")  

dbExecute(
  con,
  sprintf("CALL sqlite_attach('%s');", sqlite_path)
)

# Activate the schema
dbExecute(con, "USE chrom3_1")

# Test attachment to sqlite db
dbGetQuery(con, "
  SELECT *
  FROM chrom_data
  LIMIT 20;
")

# Create a deduplicated DuckDB table, already sorted by key
dbExecute(con, "
  CREATE OR REPLACE TABLE chrom_data_duck AS
  SELECT
      identifier,
      ANY_VALUE(info) AS info          -- all duplicate infos are identical
  FROM chrom_data                     -- <-- the SQLite view in 'main'
  GROUP BY identifier
  ORDER BY identifier;                -- keeps the table clustered
");


# If duplicates exist the primary-key build will fail anyway, so it’s smart to look once:
dbGetQuery(con, "
  SELECT identifier, COUNT(*) AS n
  FROM chrom_data_duck
  GROUP BY identifier HAVING n > 1
  LIMIT 10")

# dbExecute(con, "DELETE FROM chrom_data_duck WHERE identifier IS NULL;")

# Add the primary-key constraint (now memory-cheap)
dbExecute(con, "
  ALTER TABLE chrom_data_duck
  ADD CONSTRAINT chrom_pk PRIMARY KEY(identifier);
");

# Drop the original view so you don’t confuse the two
dbExecute(con, "DROP VIEW chrom_data;");

# Verify counts
dbGetQuery(con, "
  SELECT COUNT(*) AS n_rows
  FROM chrom_data_duck;
")

# Check db structure
dbGetQuery(con, "
  SELECT table_schema, table_name
  FROM information_schema.tables
  WHERE table_name LIKE 'chrom_data%';
")

# Remove the old SQLite view if it’s still around
dbExecute(con, "DROP VIEW IF EXISTS chrom_data;")   # no error if it’s gone

# Rename the native table 
dbExecute(con, "
  ALTER TABLE chrom_data_duck
  RENAME TO chrom_data;
")

# Check again
dbGetQuery(con, "
  SELECT table_schema, table_name
  FROM information_schema.tables
  WHERE table_name = 'chrom_data';
")

dbGetQuery(con, "
  SELECT COUNT(*) AS n_rows
  FROM chrom_data;
")

dbDisconnect(con, shutdown = TRUE) 
