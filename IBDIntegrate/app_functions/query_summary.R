#==============================================================================================================
# Function used to generate a summary of the matching datasets and genesets for each specified disease
# Returns a kable table with the number of datasets (of each experiment type) and genesets for each disease
#==============================================================================================================
generate_summary_table <- function(datasets, comparisons, diseases, experiment_types, db) {

  # If no datasets match the user query, return NULL
  if (nrow(datasets) == 0) {
    return(NULL)
  }
  # If no comparisons match the user query, return NULL
  if (is.null(comparisons)) {
    return(NULL)
  }
  
  # Combine the data description and the comparison description
  comparisons <- left_join(comparisons, datasets, by = c("dataset_acc" = "dataset_acc"))
  
  #==========
  # DATAFRAME FOR THE DATASET SUMMARY
  
  # Define the summary table
  summary_table <- data.frame("disease" = diseases)
  
  # Determine the number of DATASETS/COMPARISONS for each experiment type
  if (is.null(experiment_types)) {
    experiment_types <- get_keywords("ExperimentType", db)
  }
  for (exp_type in experiment_types) {
    num_results <- sapply(diseases, function(current_disease) {
      num_datasets <- comparisons %>%
        filter(str_detect(disease, current_disease), experiment_type == exp_type) %>%
        pull(dataset_acc) %>%
        unique() %>%
        length()
      num_comparisons <- comparisons %>%
        filter(str_detect(disease, current_disease), experiment_type == exp_type) %>%
        nrow()
      return(paste0(num_datasets, " (", num_comparisons, ")"))
    })
    
    summary_table[exp_type] <- num_results
  }
  #==========
  
  
  #==========
  # STYLED DATASET/GENESET SUMMARY
  summary_table <- summary_table %>% 
    rename(" " = "disease") %>%
    kable("html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(font_size = 12, full_width = FALSE) %>%
    row_spec(0, font_size = 14, color = "white", background = "#2C3E4C") %>%
    row_spec(1:nrow(summary_table), background = "white") %>%
    row_spec(nrow(summary_table), extra_css = "border-bottom: 1px solid black;") %>%
    column_spec(1, width = "100px", extra_css = "border-left: 1px solid black;") %>%
    column_spec(2:ncol(summary_table), width = "125px") %>%
    column_spec(ncol(summary_table), extra_css = "border-right: 1px solid black;") %>%
    add_header_above(c(" " = 1, "# datasets (# comparisons)" = length(experiment_types)),
                     extra_css = "border-left: 1px solid black; border-right: 1px solid black")
  #==========
  
  return(summary_table)
}
#==============================================================================================================
