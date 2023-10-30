#==============================================================================================================
# Function used to format the comparison columns (case vs control) in the output tables
#==============================================================================================================
format_comparison <- function(comparisons) {
  # Convert the case samples and control samples columns to lists of lists
  case_samples <- lapply(unlist(comparisons['case_ann']), function(x) str_split(x, ";")[[1]])
  control_samples <- lapply(unlist(comparisons['control_ann']), function(x) str_split(x, ";")[[1]])
  
  # Determine the attributes that are shared between the case and control sample annotations
  shared_attributes <- sapply(1:nrow(comparisons), function(x) {
    paste(intersect(case_samples[[x]], control_samples[[x]]), collapse = ", ")
  })
  shared_attributes[shared_attributes == ""] <- "---"
  
  # Determine the attributes that vary between the case and control sample annotations
  case_unique_attributes <- sapply(1:nrow(comparisons), function(x) {
    paste(case_samples[[x]][!(case_samples[[x]] %in% control_samples[[x]])], collapse = ", ")
  })
  control_unique_attributes <- sapply(1:nrow(comparisons), function(x) {
    paste(control_samples[[x]][!(control_samples[[x]] %in% case_samples[[x]])], collapse = ", ")
  })
  
  # Define the comparison column
  unique_comparisons <- paste(case_unique_attributes, "vs", control_unique_attributes)
  
  # Add the shared attributes and comparison columns
  comparisons <- mutate(comparisons, "group" = shared_attributes, "comparison" = unique_comparisons)
  
  return(comparisons)
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
    select(dataset_acc, comparison) %>%
    group_by(comparison) %>%
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
  
  # # Select columns from the comparison table
  # comparison_table <- select(comparison_table, c("number of datasets", "group", "comparison")) 
  # 
  # # Generate the DT table
  # comparison_table <- comparison_table %>%
  #   datatable(extensions = "Scroller",
  #             rownames = FALSE,
  #             selection = list(mode = "multiple", selected = 1:nrow(comparison_table)),
  #             options = list(dom = "ti",
  #                            autoWidth = FALSE,
  #                            scrollX = TRUE,
  #                            scrollY = 300,
  #                            scroller = TRUE)) %>%
  #   formatStyle(columns = c(1), width = "125px") %>%
  #   formatStyle(columns = c(2), width = "200px") %>%
  #   formatStyle(columns = c(3), width = "300px") %>%
  #   formatStyle(columns = seq_len(3), fontSize = "75%")
  
  # Select columns from the comparison table
  comparison_table <- select(comparison_table, c("number of datasets", "number of comparisons", "comparison")) 
  
  if (select_all) {
    # Generate the DT table
    comparison_table <- comparison_table %>%
      datatable(extensions = "Scroller",
                rownames = FALSE,
                selection = list(mode = "multiple", selected = 1:nrow(comparison_table)),
                options = list(dom = "ti",
                               autoWidth = FALSE,
                               scrollX = TRUE,
                               scrollY = 300,
                               scroller = TRUE)) %>%
      formatStyle(columns = c(1, 2), width = "125px") %>%
      formatStyle(columns = c(3), width = "300px") %>%
      formatStyle(columns = seq_len(3), fontSize = "75%")
  } else {
    # Generate the DT table
    comparison_table <- comparison_table %>%
      datatable(extensions = "Scroller",
                rownames = FALSE,
                selection = list(mode = "multiple"),
                options = list(dom = "ti",
                               autoWidth = FALSE,
                               scrollX = TRUE,
                               scrollY = 300,
                               scroller = TRUE)) %>%
      formatStyle(columns = c(1, 2), width = "125px") %>%
      formatStyle(columns = c(3), width = "300px") %>%
      formatStyle(columns = seq_len(3), fontSize = "75%")
  }
  
  return(comparison_table)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate the comparisons the user can select from before ranking targets
# Returns a dataframe with the comparisons, associated comparison IDs, and associated diseases
#==============================================================================================================
generate_comparison_options <- function(datasets, comparisons) {

  # If the disease was not found in the DB, return NULL
  if (nrow(datasets) == 0) {
    return(NULL)
  }
  # Pull the dataset ID's from the data description (these will be used to order the dataframe after the join)
  unique_dataset_acc <- datasets$dataset_acc
  # If no comparisons match the user specified queries (genes, sources, treatments, cell_types), return NULL
  if (is.null(comparisons)) {
    return(NULL)
  }
  
  # Format the case and control samples columns (these are replaced by the group and comparison columns)
  comparisons <- format_comparison(comparisons)
  
  # Combine the data description and the comparison description
  comparisons <- comparisons %>%
    left_join(datasets, by = c('dataset_acc' = 'dataset_acc')) %>% 
    arrange(factor(dataset_acc, levels = unique_dataset_acc))
  
  # Remove the ID, title, and summary columns
  comparisons <- select(comparisons, -c("title", "summary",
                                        "organism", "experiment_type", "source",
                                        "case_sample", "control_sample"))
  
  return(comparisons)
}
#==============================================================================================================
