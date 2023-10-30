#==============================================================================================================
# Establish the connection with the database
# Returns the connection object
#==============================================================================================================
connectDB <- function(dbname, host, port, user, password) {
  db <- dbConnect(RPostgres::Postgres(),
                  dbname = dbname, 
                  host = host, 
                  port = port, 
                  user = user,
                  password = password)
  return(db)
}
#==============================================================================================================


#==============================================================================================================
# Query the keyword table for keywords associated with a given attribute
# Returns a dataframe containing the options for each attribute
#==============================================================================================================
get_keywords <- function(attribute, db) {
  attribute <- paste0("'", attribute, "'")
  keywords_data <- dbGetQuery(db,
                              stringr::str_interp(
                                paste("SELECT keyword, old_keyword",
                                      "FROM keyword",
                                      "WHERE keyword_type LIKE ${attribute}")))
  keywords_data <- keywords_data %>%
    mutate(old_keyword = ifelse(sapply(old_keyword, function(x) is.null(x)), keyword, old_keyword))
  keywords <- keywords_data$old_keyword
  #names(keywords) <- keywords_data$keyword
  
  return(keywords)
}
#==============================================================================================================


#==============================================================================================================
# Function used to query the dataset table for all datasets
#==============================================================================================================
get_datasets <- function(db) {

  # Pull the datasets
  datasets <- dbGetQuery(db,
                         stringr::str_interp(
                           paste("SELECT dataset_acc, title, disease, organism, experiment_type, source, cell_type, treatment, timepoint, dose",
                                 "FROM dataset")))
  
  # Sort the data description by dataset ID, organism, and experiment type
  datasets <- arrange(datasets, disease, dataset_acc)
  
  # Wrap the titles/summaries (and add ellipsis for long titles/summaries)
  datasets <- datasets %>%
    mutate(title = str_trunc(title, width = 40, ellipsis = "..."),
           disease = str_trunc(disease, width = 25, ellipsis = "..."),
           source = str_trunc(source, width = 25, ellipsis = "..."),
           treatment = str_trunc(treatment, width = 25, ellipsis = "..."))
  
  return(datasets)
}
#==============================================================================================================