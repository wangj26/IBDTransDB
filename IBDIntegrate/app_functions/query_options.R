#==============================================================================================================
# Function used to dynamically change the query options presented to the user
#==============================================================================================================
change_query_options <- function(session, query_option, input, disease_input, organism_input, tissue_input, treatment_input, experiment_type_input, cell_type_input, timepoint_input,db) {
  
  # Define the query options
  options <- c("disease", "organism", "tissue", "treatment", "experiment_type", "cell_type", "timepoint")
  
  # Pull the data descriptions from datasets associated with the specified options
  data_desc <- get_datasets(disease_input, organism_input, tissue_input, treatment_input, experiment_type_input, cell_type_input, timepoint_input, db)
  
  # Update the query options
  for (current_option in options) {
    if (current_option == query_option) {
      next
    }
    if (is.null(input[[paste0(current_option, "_input")]])) {
      if (current_option == "tissue") {
        keywords <- strsplit(data_desc[["source"]], ";") 
        keywords <- unique(unlist(keywords))
      } else {
        keywords <- strsplit(data_desc[[current_option]], ";") 
        keywords <- unique(unlist(keywords))
      }
      updateSelectInput(session, paste0(current_option, "_input"), choices = keywords)
    }
  }
}
#==============================================================================================================