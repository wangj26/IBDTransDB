#==============================================================================================================
# Function used to create a DT table from a dataframe of datasets
# Returns a DT table
#==============================================================================================================
style_table <- function(dataset_table) {
  dataset_table <- datatable(dataset_table,
                             extensions = "Scroller",
                             selection = "single",
                             options = list(order = list(list(3, "asc"), list(1, "desc")),
                                            autoWidth = TRUE,
                                            scrollX = TRUE,
                                            scrollY = 670,
                                            scroller = TRUE,
                                            dom = '<"top">ift'
                             ),
                             rownames = FALSE) %>%
    formatStyle(columns = c(1, 3, 4, 5, 6, 7, 8, 9), width = "150px") %>%
    formatStyle(columns = c(2), width = "300px") %>%
    formatStyle(columns = 1:9, fontSize = "75%")
  
  return(dataset_table)
}
#==============================================================================================================


#==============================================================================================================
# Function used to filter the table with datasets when the user modifies an option
# Returns a DT table
#==============================================================================================================
update_table <- function(dataset_table, input) {
  
  # Create a named list of the options for attributes specified by the user
  attributes <- list("disease" = input$disease_input,
                     "source" = input$tissue_input, 
                     "cell_type" = input$cell_type_input, 
                     "treatment" = input$treatment_input, 
                     "timepoint" = input$timepoint_input)
  
  print(input$disease_input)
  print(input$tissue_input)
  print(input$cell_type_input)
  print(input$treatment_input)
  print(input$timepoint_input)
  # Determine the rows to keep in the dataset table based on the user specified options
  rows_to_keep <- lapply(names(attributes), function(x) {
    if (is.null(attributes[[x]])) {
      return(rep(TRUE, nrow(dataset_table)))
    } else {
      return(sapply(dataset_table[[x]], function(y) any(str_detect(y, attributes[[x]]))))
    }
  })
  rows_to_keep <- Reduce("&", rows_to_keep)
  
  # Filter the dataset table for datasets that match the values for the modified attributes
  dataset_table <- filter(dataset_table, rows_to_keep)
  print(dim(dataset_table))
  
  return(dataset_table)
}
#==============================================================================================================
