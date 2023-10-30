#==============================================================================================================
# Function used to create a DT table displaying dataset descriptions and selecting datasets
# Returns a DT table
#==============================================================================================================
generate_dataset_table <- function(datasets, select_all) {
  
  # Wrap the titles, summaries, and treatments (and add ellipsis) and rename columns
  datasets <- datasets %>%
    mutate(title = str_trunc(title, width = 35, ellipsis = "..."),
           source = str_trunc(source, width = 25, ellipsis = "..."),
           cell_type = str_trunc(cell_type, width = 25, ellipsis = "..."),
           treatment = str_trunc(treatment, width = 25, ellipsis = "...")) %>%
    rename("dataset ID" = "dataset_acc", "experiment type" = "experiment_type", "cell type" = "cell_type")
  
  # Generate the DT table
  if (select_all) {
    datasets <- datasets %>%
      datatable(class = "display nowrap",
                extensions = "Scroller",
                rownames = FALSE,
                selection = list(mode = "multiple", selected = 1:nrow(datasets)),
                options = list(dom = "ti",
                               autoWidth = FALSE,
                               scrollX = TRUE,
                               scrollY = 500,
                               scroller = TRUE)) %>%
      formatStyle(columns = seq_len(8), fontSize = "75%")
  } else {
    datasets <- datasets %>%
      datatable(class = "display nowrap",
                extensions = "Scroller",
                rownames = FALSE,
                selection = list(mode = "multiple"),
                options = list(dom = "ti",
                               autoWidth = FALSE,
                               scrollX = TRUE,
                               scrollY = 500,
                               scroller = TRUE)) %>%
      formatStyle(columns = seq_len(8), fontSize = "75%")
  }
  
  return(datasets)
}
#==============================================================================================================


#==============================================================================================================
# Function used to create a DT table displaying comparisons and selecting comparisons
# Returns a DT table
#==============================================================================================================
generate_comparison_table <- function(dataset_acc_list, comparisons) {
  
  # If there are no dataset IDs, return NULL
  if (length(dataset_acc_list) == 0) {
    return(NULL)
  }
  # If there are no comparisons, return NULL
  if (nrow(comparisons) == 0) {
    return(NULL)
  }
  
  # Filter the comparisons for the selected dataset IDs
  comparisons <- filter(comparisons, dataset_acc %in% dataset_acc_list)
  
  # Format the comparisons
  comparisons <- format_comparison(comparisons)
  
  # # Select columns from the comparison table and determine the number of datasets that contain each comparison
  # comparisons <- comparisons %>%
  #   select(c("group", "comparison", "case_ann", "control_ann")) %>%
  #   group_by(group, comparison, case_ann, control_ann) %>%
  #   summarise("number of datasets" = n()) %>%
  #   ungroup() %>%
  #   as.data.frame() %>%
  #   arrange(comparison, group)
  
  # Select columns from the comparison table and determine the number of datasets that contain each comparison
  comparisons <- comparisons %>%
    select(dataset_acc, group, comparison, full_comparison) %>%
    group_by(group, comparison, full_comparison) %>%
    mutate("number of comparisons" = n(),
           "number of datasets" = n_distinct(dataset_acc)) %>%
    ungroup() %>%
    as.data.frame() %>%
    select(-dataset_acc) %>%
    unique() %>%
    arrange(comparison)
  
  return(comparisons)
}
#==============================================================================================================


#==============================================================================================================
# Function used to create a DT table displaying comparisons and selecting comparisons
# Returns a DT table
#==============================================================================================================
style_comparison_table <- function(comparison_table, select_all) {
  
  # Select columns from the comparison table
  comparison_table <- select(comparison_table, c("comparison", "group", "number of datasets", "number of comparisons")) 
  
  # Generate the DT table
  if (select_all) {
    comparison_table <- comparison_table %>%
      datatable(class = "display nowrap",
                extensions = "Scroller",
                rownames = FALSE,
                selection = list(mode = "multiple", selected = 1:nrow(comparison_table)),
                options = list(dom = "ti",
                               autoWidth = FALSE,
                               scrollX = TRUE,
                               scrollY = 500,
                               scroller = TRUE)) %>%
      formatStyle(columns = seq_len(3), fontSize = "75%")
  } else {
    comparison_table <- comparison_table %>%
      datatable(class = "display nowrap",
                extensions = "Scroller",
                rownames = FALSE,
                selection = list(mode = "multiple"),
                options = list(dom = "ti",
                               autoWidth = FALSE,
                               scrollX = TRUE,
                               scrollY = 500,
                               scroller = TRUE)) %>%
      formatStyle(columns = seq_len(4), fontSize = "75%")
  }
  
  return(comparison_table)
}
#==============================================================================================================


#==============================================================================================================
# Function used to format the comparison columns (case vs control) in the output tables
#==============================================================================================================
format_comparison <- function(comparisons) {
  # Convert the case samples and control samples columns to lists of lists
  case_samples <- lapply(comparisons$case_ann, function(x) str_split(x, ";")[[1]])
  control_samples <- lapply(comparisons$control_ann, function(x) str_split(x, ";")[[1]])
  
  # Determine the attributes that are shared between the case and control sample annotations
  shared_attributes <- mapply(function(x, y) {
    paste(intersect(x, y), collapse = ", ")
  }, case_samples, control_samples)
  shared_attributes[shared_attributes == ""] <- "---"
  
  # Determine the attributes that vary between the case and control sample annotations
  case_unique_attributes <- mapply(function(x,y) {
    paste(x[!(x %in% y)], collapse = ", ")
  }, case_samples, control_samples)
  control_unique_attributes <- mapply(function(x,y) {
    paste(y[!(y %in% x)], collapse = ", ")
  }, case_samples, control_samples)
  
  # Define the comparison column
  unique_comparisons <- paste(case_unique_attributes, "vs", control_unique_attributes)
  
  # Add the shared attributes and comparison columns
  comparisons <- comparisons %>%
    mutate("group" = shared_attributes, 
           "comparison" = unique_comparisons,
           "full_comparison" = paste(case_ann, control_ann, sep = " vs "))
  
  return(comparisons)
}
#==============================================================================================================
