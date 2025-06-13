# Function to get info from the database for a given identifier
get_info <- function(identifier, dbPath) {
  require(RSQLite)
  conn <- dbConnect(RSQLite::SQLite(), dbname = dbPath)
  query <- sprintf("SELECT info FROM chrom_data WHERE identifier = %d", identifier)
  result <- dbGetQuery(conn, query)
  dbDisconnect(conn)  # Close the connection
  if (nrow(result) == 0) {
    return(NA)  # Return NA if no match is found
  }
  return(result$info[1])
}