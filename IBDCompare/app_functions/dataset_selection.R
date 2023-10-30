#==============================================================================================================
# Generate a table allowing users to select specific datasets of interest (all datasets selected by default)
# Returns a DT object with dataset descriptions
#==============================================================================================================
generate_dataset_table <- function(dataset_table, select_all) {
  
  if (select_all) {
    # Create the DT object
    dataset_table <- datatable(dataset_table,
                               extensions = "Scroller",
                               selection = list(mode = "multiple", selected = 1:nrow(dataset_table)),
                               options = list(order = list(list(3, "asc"), list(1, "desc")),
                                              autoWidth = TRUE,
                                              scrollX = TRUE,
                                              scrollY = 825,
                                              scroller = TRUE,
                                              dom = '<"top">ift'
                               ),
                               rownames = FALSE)
  } else {
    # Create the DT object
    dataset_table <- datatable(dataset_table,
                               extensions = "Scroller",
                               selection = list(mode = "multiple"),
                               options = list(order = list(list(3, "asc"), list(1, "desc")),
                                              autoWidth = TRUE,
                                              scrollX = TRUE,
                                              scrollY = 825,
                                              scroller = TRUE,
                                              dom = '<"top">ift'
                               ),
                               rownames = FALSE)
  }
  dataset_table
  # Format the DT object columns
  dataset_table <- dataset_table %>%
    formatStyle(columns = c(1, 3, 4, 5, 6, 7, 8, 9), width = "150px") %>%
    formatStyle(columns = c(2), width = "300px") %>%
    formatStyle(columns = 1:10, fontSize = "75%")
  
  return(dataset_table)
}
#==============================================================================================================
