library(RSQLite)

# Connect to the SQLite database
db <- dbConnect(SQLite(), dbname = "chrom3_1_c.db", timeout = 10000)

# Create table if it doesn't exist
sql <- paste(
  "CREATE TABLE IF NOT EXISTS chrom_data",
  "(identifier INTEGER PRIMARY KEY, info TEXT)"
)
dbExecute(db, sql)


# Function to handle bulk inserts
bulk_insert <- function(data) {
  # Convert list to data frame
  params <- do.call(rbind, data)

  # Prepare the SQL insert statement
  sql <- "INSERT OR IGNORE INTO chrom_data (identifier, info) VALUES (?, ?)"
  stmt <- dbSendQuery(db, sql)

  # Bind parameters and execute row by row
  for (i in seq_len(nrow(params))) {
    dbBind(stmt, params[i, , drop = FALSE])
  }

  # Clear the statement
  dbClearResult(stmt)
}

# Open the file connection
con <- file("data/input/annotations/chrom3_1_c.trn", "r")

# Initialize variables for bulk insert
data_buffer <- list()
buffer_size <- 10000 # Adjust based on performance and system memory
batch_counter <- 0
max_chars <- 61

# Begin a transaction
dbBegin(db)

# Read the file line by line
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  id <- as.integer(substr(line, 1, 10)) # Convert identifier to integer
  info <- trimws(substr(line, 12, max_chars)) # Trim leading/trailing spaces

  data_buffer[[length(data_buffer) + 1]] <- list(id, info)

  # When buffer reaches the specified size, do a bulk insert
  if (length(data_buffer) >= buffer_size) {
    bulk_insert(data_buffer)
    data_buffer <- list() # Reset buffer after insertion
    batch_counter <- batch_counter + 1
    if (batch_counter >= 10) {
      dbCommit(db) # Commit after every 10 batches
      dbBegin(db) # Begin a new transaction
      batch_counter <- 0
    }
  }
}

# Insert any remaining data in the buffer
if (length(data_buffer) > 0) {
  bulk_insert(data_buffer)
}

# Commit any remaining transaction
dbCommit(db)

# Close file and finalize database transaction
close(con)
dbDisconnect(db)