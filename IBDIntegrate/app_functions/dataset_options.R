#==============================================================================================================
# Function used to create a DT table displaying dataset descriptions and selecting datasets
# Returns a DT table
#==============================================================================================================
generate_dataset_table <- function(datasets) {
  
  # Wrap the titles, summaries, and treatments (and add ellipsis)
  datasets <- datasets %>%
    mutate(title = str_trunc(title, width = 35, ellipsis = "..."),
           source = str_trunc(source, width = 25, ellipsis = "..."),
           cell_type = str_trunc(cell_type, width = 25, ellipsis = "..."),
           treatment = str_trunc(treatment, width = 25, ellipsis = "..."))
  
  # Generate the DT table
  datasets <- datasets %>%
    rename("experiment type" = "experiment_type", "cell type" = "cell_type") %>%
    datatable(extensions = "Scroller",
              rownames = FALSE,
              selection = list(mode = "multiple", selected = 1:nrow(datasets)),
              options = list(dom = "ti",
                             autoWidth = FALSE,
                             scrollX = TRUE,
                             scrollY = 300,
                             scroller = TRUE)) %>%
    formatStyle(columns = c(1, 3, 4, 5, 6, 7, 8), width = "125px") %>%
    formatStyle(columns = c(2), width = "300px") %>%
    formatStyle(columns = seq_len(8), fontSize = "75%")
  
  return(datasets)
}
#==============================================================================================================