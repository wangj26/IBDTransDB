#==============================================================================================================
# Function used to prepare files associated with a given dataset for downloading
# Returns a named list of dataframes
#==============================================================================================================
prepare_files <- function(dataset_acc, selected_comparisons, dataset_description, comparisons, sample_ann) {
  
  # Pull the expression data and convert to wide format
  exp_data <- get_full_exp_data(dataset_acc, db)
  ### REMOVE DUPLICATE SAMPLE+GENE PAIRS ###
  ### THIS WILL BE REMOVED ###
  rows_to_keep <- !duplicated(paste0(exp_data$sample_acc, "_", exp_data$gene))
  exp_data <- filter(exp_data, rows_to_keep)
  ##########
  exp_data <- exp_data %>%
    pivot_wider(names_from = gene, values_from = expr) %>% 
    as.data.frame()
  rownames(exp_data) <- exp_data$sample_acc
  exp_data <- select(exp_data, -c("sample_acc"))
  exp_data <- t(exp_data)
  
  # Order the sample annotation dataframe by the sample order in the expression data
  sample_ann <- sample_ann[match(rownames(exp_data), sample_ann$sample_id),]
  
  # Convert the dataset description to a two column dataframe
  dataset_description <- data.frame("Attribute" = colnames(dataset_description),
                                    "Value" = unlist(unname(dataset_description[1,])))
  
  # Filter the comparisons and comparison data for the selected comparisons
  comparisons <- filter(comparisons, comparison %in% selected_comparisons) 
  selected_comparison_ids <- pull(comparisons, id)
  comparison_data <- get_comparison_data(selected_comparison_ids, db)

  # Create a named list used to store all files to be downloaded
  files_to_download <- list("data" = exp_data,
                            "data_description" = dataset_description,
                            "comparison_description" = select(comparisons, -c(id, comparison)),
                            "sample_annotations" = sample_ann)
  
  # Add selected comparison data to the named list
  for (current_comparison in selected_comparisons) {
    # Pull the comparison ID associated with the current comparison
    current_comparison_id <- comparisons %>%
      filter(comparison == current_comparison) %>%
      pull(id) %>%
      unique()
    
    # Filter the comparison data based on the current comparison ID
    current_comparison_data <- comparison_data %>%
      filter(comparison_id == current_comparison_id) %>%
      select(-comparison_id)
    
    # Add the comparison data associated with the current comparison to the named list
    files_to_download[[gsub(" ", "_", current_comparison)]] <- current_comparison_data
  }
  
  return(files_to_download)
}
#==============================================================================================================


#==============================================================================================================
# Function used to prepare files associated with a given dataset for downloading (single cell)
# Returns a named list of dataframes
#==============================================================================================================
prepare_files_sc <- function(dataset_acc, selected_comparisons, dataset_description, comparisons, sample_ann) {
  
  # Order the sample annotation dataframe by the sample order in the expression data
  sample_ann <- sample_ann[match(rownames(exp_data), sample_ann$sample_id),]
  
  # Convert the dataset description to a two column dataframe
  dataset_description <- data.frame("Attribute" = colnames(dataset_description),
                                    "Value" = unlist(unname(dataset_description[1,])))
  
  # Filter the comparisons and comparison data for the selected comparisons
  comparisons <- filter(comparisons, comparison %in% selected_comparisons) 
  selected_comparison_ids <- pull(comparisons, id)
  comparison_data <- get_comparison_data(selected_comparison_ids, db)
  
  # Create a named list used to store all files to be downloaded
  files_to_download <- list("data_description" = dataset_description,
                            "comparison_description" = select(comparisons, -c(id, comparison)),
                            "sample_annotations" = sample_ann)
  
  # Add selected comparison data to the named list
  for (current_comparison in selected_comparisons) {
    # Pull the comparison ID associated with the current comparison
    current_comparison_id <- comparisons %>%
      filter(comparison == current_comparison) %>%
      pull(id) %>%
      unique()
    
    # Filter the comparison data based on the current comparison ID
    current_comparison_data <- comparison_data %>%
      filter(comparison_id == current_comparison_id) %>%
      select(-comparison_id)
    
    # Add the comparison data associated with the current comparison to the named list
    files_to_download[[gsub(" ", "_", current_comparison)]] <- current_comparison_data
  }
  
  return(files_to_download)
}
#==============================================================================================================
