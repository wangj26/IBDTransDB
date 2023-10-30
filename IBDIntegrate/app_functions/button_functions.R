#==============================================================================================================
# Function used to generate options for the comparison filter dropdown menu
#==============================================================================================================
comparison_filter_options <- function(results_table) {
  filter_options <- paste0(results_table["data ID"], ": ", 
                           results_table["case samples"], " vs ", 
                           results_table["control samples"])
  
  return(filter_options)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate options for the comparisons filter dropdown menu
#==============================================================================================================
generate_comparison_table <- function(genes, diseases, organisms, sources, cell_types, treatments, experiment_types, platforms, pval_threshold, logfc_threshold, normalization_methods, db) {
  
  # Create a list used to store the dataframes per disease
  comparison_df_list <- list()
  
  # Iterate over every disease
  for (disease in diseases) {
    
    # Pull the list of data ID's
    dataDesc <- get_dataDesc(disease, organisms, sources, treatments, experiment_types, platforms, cell_types, normalization_methods, db)
    
    # If the disease was not found in the DB, return NULL
    if (nrow(dataDesc) == 0) {
      return(NULL)
    }
    
    # Pull the comparison description and comparison data
    comp_list <- get_comparisonData(unlist(dataDesc['dataset_acc']), genes, sources, treatments, cell_types, db)
    compDesc <- comp_list[[1]]
    compData <- comp_list[[2]]
    
    # If no comparisons match the user specified queries (genes, sources, treatments, cell_types), return NULL
    if (nrow(compData) == 0) {
      return(NULL)
    }
    
    # Add gene columns to the comparison description table 
    # Initialize values to NA - some of these will be changed later
    for (gene in genes) {
      compDesc[gene] <- 'NA'
    }
    
    # Combine the data description and the comparison description
    compDesc <- left_join(dataDesc, compDesc, by = c('dataset_acc' = 'dataset_acc'))
    
    # Pull the data ID's before the data_acc column is converted to cell_spec elements (NEEDED for plotting)
    dataID_list <- unique(compDesc$dataset_acc)
    
    # Iterate over every gene
    for (current_gene in genes) {
      # Define the list used to store the labels shown in the table
      gene_labels <- c()
      # Define the list used to store the comparison information displayed when hovering over labels
      gene_data <- c()
      
      # Iterate over every comparison
      for (id in compDesc$id) {
        included_genes <- compData[compData$comparison_id == id,]
        included_genes <- included_genes$gene
        
        # Determine if the given comparison contains data on the given gene
        if (current_gene %in% included_genes) {
          comparison_info <- compData[compData$comparison_id == id & compData$gene == current_gene,]
          
          # Determine if the given gene passes the user specified thresholds (i.e. p-value, q-value, logFC)
          significance_status <- comparison_info$p_value <= pval_threshold & abs(comparison_info$log_fc) >= logfc_threshold
          if (significance_status) {
            ### ADD INFO TO COPY ###
            compDesc[compDesc$id == id, current_gene] <- "info"
          } else {
            ### ADD INFO TO COPY ###
            compDesc[compDesc$id == id, current_gene] <- "NS"
          }
          
        } else {
          ### ADD INFO TO COPY ###
          compDesc[compDesc$id == id, current_gene] <- "NA"
        }
        
      }
      
    }
    
    # Remove the ID, title, and summary columns
    compDesc <- select(compDesc, -c('id', 'title', 'summary',
                                    'organism', 'experiment_type', 'source',
                                    'case_sample', 'control_sample'))
    # Change column names
    new_col <- c('data ID', 'case samples', 'control samples')
    colnames(compDesc)[1:3] <- new_col
    
    # Replace ; with , in the case samples and control samples columns
    compDesc$`case samples` <- gsub(';', ', ', compDesc$`case samples`)
    compDesc$`control samples` <- gsub(';', ', ', compDesc$`control samples`)
    
    # Append the comparison description to the list
    comparison_df_list <- append(comparison_df_list, compDesc)
  }

  # Merge the dataframes
  comparison_df <- do.call("rbind", comparison_df_list)
  
  return(comparison_df)
}
#==============================================================================================================
  