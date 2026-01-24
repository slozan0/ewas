# Function to get info from the database for a given identifier
get_info <- function(identifier, db_path) {
  require(RSQLite)

  conn <- dbConnect(RSQLite::SQLite(), dbname = db_path)
  query <- "SELECT info FROM chrom_data WHERE identifier = ?"
  result <- dbGetQuery(conn, query, params = list(identifier))
  dbDisconnect(conn)  # Close the connection

  if (nrow(result) == 0) {
    return(NA)  # Return NA if no match is found
  }

  result$info[1]
}