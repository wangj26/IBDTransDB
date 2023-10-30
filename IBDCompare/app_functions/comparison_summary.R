#==============================================================================================================
# Function used to generate a simple table for each disease used for the summary sections
#==============================================================================================================
generate_dataset_summary_table <- function(genes, selected_disease, datasets, comparisons, comparison_data) {
  # Filter the datasets for the given disease
  datasets <- datasets %>%
    select(-id) %>%
    filter(str_detect(disease, selected_disease))
  # If the disease was not found in the DB, return NULL
  if (nrow(datasets) == 0) {
    return(NULL)
  }
  
  # Filter the comparisons and comparison data
  comparisons <- filter(comparisons, dataset_acc %in% datasets$dataset_acc)
  selected_comparison_ids <- comparisons %>%
    filter(dataset_acc %in% datasets$dataset_acc) %>%
    pull(id)
  comparison_data <- filter(comparison_data, comparison_id %in% selected_comparison_ids)
  
  # Combine the data description and the comparison description
  comparisons <- left_join(comparisons, datasets, by = c("dataset_acc" = "dataset_acc"))
  
  # Add gene columns to the comparison description table 
  # Initialize values to NA - some of these will be changed later
  for (gene in genes) {
    comparisons[gene] <- "NA"
  }
  
  # Iterate over every gene
  for (current_gene in genes) {
    
    # Iterate over every comparison
    for (current_id in comparisons$id) {
      included_genes <- comparison_data %>%
        filter(comparison_id == current_id) %>%
        pull(gene)

      # Determine if the given comparison contains data on the given gene
      if (current_gene %in% included_genes) {
        comparisons[comparisons$id == current_id, current_gene] <- "info"
      }
    }
  }
  
  # Remove the unneeded columns, rename some columns, and format the group annotations
  comparisons <- comparisons %>%
    select(-c('id', 'title', 'summary', 'disease', 'organism', 'experiment_type', 'source', 'case_sample', 'control_sample')) %>%
    rename("case samples" = "case_ann", "control samples" = "control_ann") %>%
    mutate(`case samples` = gsub(';', ', ', `case samples`),
           `control samples` = gsub(';', ', ', `control samples`))
  
  return(comparisons)
}
#==============================================================================================================


#==============================================================================================================
# Function used to compute summary stats for each disease (e.g. number of datasets returned)
#==============================================================================================================
generate_dataset_summary_stats <- function(genes, results_table) {
  summary_df <- data.frame("datasets" = rep(0, length(genes)),
                           "comparisons" = rep(0, length(genes)))
  rownames(summary_df) <- genes
  
  # Calculate the total number of comparisons and datasets found for each gene
  if (is.null(results_table)) {
    total_num_comparisons <- 0
    total_num_datasets <- 0
  } else {
    # Calculate the number of comparisons found for each gene 
    num_found <- sapply(genes, function(gene) sum(results_table[gene] != "NA"))
    summary_df["comparisons"] <- num_found
    # Calculate the total number of comparisons that matched the user specified queries
    total_num_comparisons <- nrow(results_table)
    
    # Calculate the number of datasets found for each gene
    num_found <- sapply(genes, function(gene) {
      dataset_list <- sapply(unique(results_table[["dataset_acc"]]), function(dataID) {
        data <- results_table %>%
          filter(results_table$dataset_acc == dataID) %>%
          pull(gene)
        return(all(data != "NA"))
      })
      return(sum(dataset_list))
    })
    summary_df["datasets"] <- num_found
    # Calculate the total number of datasets that matched the user specified queries
    total_num_datasets <- length(unique(results_table[["dataset_acc"]]))
  }
  
  return(list("summary_df" = summary_df, "total_num_datasets" = total_num_datasets, "total_num_comparisons" = total_num_comparisons))
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate a simple table for each disease used for the summary sections
#==============================================================================================================
pathway_generate_dataset_summary_table <- function(pathways, selected_disease, datasets, comparisons, pathway_data) {
  # Filter the datasets for the given disease
  datasets <- datasets %>%
    select(-id) %>%
    filter(str_detect(disease, selected_disease))
  # If the disease was not found in the DB, return NULL
  if (nrow(datasets) == 0) {
    return(NULL)
  }
  
  # Filter the comparisons and comparison data
  comparisons <- filter(comparisons, dataset_acc %in% datasets$dataset_acc)
  selected_comparison_ids <- comparisons %>%
    filter(dataset_acc %in% datasets$dataset_acc) %>%
    pull(id)
  pathway_data <- filter(pathway_data, comparison_id %in% selected_comparison_ids)
  # If no comparisons match the user specified queries (genes, sources, treatments, cell_types), return NULL
  if (is.null(pathway_data)) {
    return(NULL)
  }
  # Combine the data description and the comparison description
  comparisons <- left_join(comparisons, datasets, by = c("dataset_acc" = "dataset_acc"))
  
  # Add gene columns to the comparison description table 
  # Initialize values to NA - some of these will be changed later
  for (pathway in pathways) {
    comparisons[pathway] <- "NA"
  }
  
  # Iterate over every gene
  for (current_pathway in pathways) {
    
    # Iterate over every comparison
    for (current_id in comparisons$id) {
      included_pathways <- pathway_data %>%
        filter(comparison_id == current_id) %>%
        pull(description)
      
      # Determine if the given comparison contains data on the given gene
      if (current_pathway %in% included_pathways) {
        comparisons[comparisons$id == current_id, current_pathway] <- "info"
      }
    }
  }
  
  # Remove the unneeded columns, rename some columns, and format the group annotations
  comparisons <- comparisons %>%
    select(-c('id', 'title', 'summary', 'disease', 'organism', 'experiment_type', 'source', 'case_sample', 'control_sample')) %>%
    rename("case samples" = "case_ann", "control samples" = "control_ann") %>%
    mutate(`case samples` = gsub(';', ', ', `case samples`),
           `control samples` = gsub(';', ', ', `control samples`))
  
  return(comparisons)
}
#==============================================================================================================


#==============================================================================================================
# Function used to compute summary stats for each disease (e.g. number of datasets returned)
#==============================================================================================================
pathway_generate_dataset_summary_stats <- function(pathways, results_table) {
  summary_df <- data.frame("datasets" = rep(0, length(pathways)), 
                           "comparisons" = rep(0, length(pathways)))
  rownames(summary_df) <- pathways
  
  # Calculate the total number of comparisons and datasets found for each gene
  if (is.null(results_table)) {
    total_num_comparisons <- 0
    total_num_datasets <- 0
  } else {
    # Calculate the number of comparisons found for each gene 
    num_found <- sapply(pathways, function(pathway) sum(results_table[pathway] != "NA"))
    summary_df["comparisons"] <- num_found
    # Calculate the total number of comparisons that matched the user specified queries
    total_num_comparisons <- nrow(results_table)
    
    # Calculate the number of datasets found for each gene
    num_found <- sapply(pathways, function(pathway) {
      dataset_list <- sapply(unique(results_table[["dataset_acc"]]), function(dataID) {
        data <- results_table %>%
          filter(dataset_acc == dataID) %>%
          pull(pathway)
        return(all(data != "NA"))
      })
      return(sum(dataset_list))
    })
    summary_df["datasets"] <- num_found
    # Calculate the total number of datasets that matched the user specified queries
    total_num_datasets <- length(unique(results_table[["dataset_acc"]]))
  }
  
  return(list("summary_df" = summary_df, "total_num_datasets" = total_num_datasets, "total_num_comparisons" = total_num_comparisons))
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate summary tables for the summary section
# Table 1 displays the proportion of datasets that contain at least one significant comparison per gene per disease
# Table 2 displays the proportion of comparisons that are significant per gene per disease
#==============================================================================================================
generate_results_summaries <- function(genes, diseases, experiment_types, datasets, comparisons, comparison_data, pval_option, pval_threshold, logfc_threshold, db) {
  
  datasets <- select(datasets, -id)
  
  metrics_df_list_1 <- list()
  metrics_df_list_2 <- list()
  
  # Determine if the user specified p-value or adjusted p-value under the SIGNIFICANCE THRESHOLDS option
  pval_col <- pval_option
  
  # If experiment_type is NULL, set experiment_type to all experiment types in the DB
  if (is.null(experiment_types)) {
    experiment_types <- get_keywords("ExperimentType", db)
  }
  
  #================
  # Iterate over every disease
  for (current_disease in diseases) {
    # Filter the datasets for the given disease
    current_datasets <- filter(datasets, str_detect(disease, current_disease))

    # If no datasets match the user specified queries, specify NA in the metrics dataframes and continue to the next iteration
    if (nrow(current_datasets) == 0) {
      metrics_df_1 <- data.frame("Disease" = rep(current_disease, length(genes)), "gene" = genes)
      for (type in experiment_types) {
        metrics_df_1[type] <- "NA"
      }
      metrics_df_list_1 <- append(metrics_df_list_1, list(metrics_df_1))
      
      metrics_df_2 <- data.frame("Disease" = rep(current_disease, length(genes)), "gene" = genes)
      for (type in experiment_types) {
        metrics_df_2[type] <- "NA"
      }
      metrics_df_list_2 <- append(metrics_df_list_2, list(metrics_df_2))
      
      next
    }
    
    # Filter the comparisons and comparison data
    current_comparisons <- filter(comparisons, dataset_acc %in% current_datasets$dataset_acc)
    current_comparison_ids <- comparisons %>%
      filter(dataset_acc %in% current_datasets$dataset_acc) %>%
      pull(id)
    current_comparison_data <- filter(comparison_data, comparison_id %in% current_comparison_ids)
    
    # If no comparisons match the user specified queries, specify NA in the metrics dataframes and continue to the next iteration
    if (is.null(current_comparison_data)) {
      metrics_df_1 <- data.frame("Disease" = rep(current_disease, length(genes)), "gene" = genes)
      for (type in experiment_types) {
        metrics_df_1[type] <- "NA"
      }
      metrics_df_list_1 <- append(metrics_df_list_1, list(metrics_df_1))
      
      metrics_df_2 <- data.frame("Disease" = rep(current_disease, length(genes)), "gene" = genes)
      for (type in experiment_types) {
        metrics_df_2[type] <- "NA"
      }
      metrics_df_list_2 <- append(metrics_df_list_2, list(metrics_df_2))
      
      next
    }
    
    # Combine the data description and the comparison description
    current_comparisons <- left_join(current_comparisons, current_datasets, by = c("dataset_acc" = "dataset_acc"))
    
    # Add gene columns to the comparison description table 
    # Initialize values to NA - some of these will be changed later
    for (gene in genes) {
      current_comparisons[gene] <- "NA"
    }
    
    # Pull the data ID's before the data_acc column is converted to cell_spec elements (NEEDED for plotting)
    dataset_acc_list <- unique(current_comparisons$dataset_acc)
    
    # Iterate over every gene
    for (current_gene in genes) {
      
      # Determine the significance status across all selected comparisons for the given gene 
      current_comparisons[current_gene] <- sapply(current_comparisons$id, function(current_id) {
        included_genes <- current_comparison_data %>%
          filter(comparison_id == current_id) %>%
          pull(gene)
        included_genes <- unique(included_genes) ##some comparisons may be uploaded twice
        # Determine if the given comparison contains data on the given gene
        if (current_gene %in% included_genes) {
          comparison_info <- unique(filter(current_comparison_data, comparison_id == current_id & gene == current_gene))
          # Determine if the given gene passes the user specified thresholds (i.e. p-value, q-value, logFC)
          significance_status <- comparison_info[[pval_col]] <= pval_threshold & abs(comparison_info$log_fc) >= logfc_threshold
          return(ifelse(significance_status, "info", "NS"))
        } else {
          return("NA")
        }
      })
      
    }
    
    
    ### Generate a table that displays the proportion of datasets that contain ###
    ### at least one significant comparison for each experiment type ###
    metrics_df_1 <- data.frame("Disease" = rep(current_disease, length(genes)), "gene" = genes)
    # Iterate over every experiment type
    for (exp_type in experiment_types) {
      
      # Pull the data associated with the given experiment type
      gene_data <- current_comparisons %>% 
        filter(experiment_type == exp_type) %>%
        select(-c("experiment_type"))
      
      metrics_df_1[exp_type] <- ""
      for (gene in genes) {
        # Handle the case when there are no datasets associated with the given disease and experiment type
        if (nrow(gene_data) == 0) {
          metrics_df_1[metrics_df_1$gene == gene, exp_type] <- "NA"
          next
        }
        
        num_significant <- sapply(unique(gene_data$dataset_acc), function(x) {
          dataset_data <- gene_data[gene_data$dataset_acc == x,]
          return(any(dataset_data[[gene]] == 'info'))
        })
        num_significant = sum(num_significant)
        
        num_total <- sapply(unique(gene_data$dataset_acc), function(x) {
          dataset_data <- gene_data[gene_data$dataset_acc == x,]
          return(any(dataset_data[[gene]] == 'info' | dataset_data[[gene]] == 'NS'))
        })
        num_total <- sum(num_total)
        
        percent_signifcant <- round(num_significant/num_total * 100, 2)
        
        # If the given gene was not found in any datasets, the associated entry in the table is NA
        if (num_total == 0) {
          metrics_df_1[metrics_df_1$gene == gene, exp_type] <- "NA"
        } else {
          metrics_df_1[metrics_df_1$gene == gene, exp_type] <- paste0(num_significant, "/", num_total, " (", percent_signifcant,"%)")
        }
        
      }
      
    }
    
    # Append the metrics dataframe to the list
    metrics_df_list_1 <- append(metrics_df_list_1, list(metrics_df_1))
    ###
    
    
    # Remove the unneeded columns from the comparisons
    current_comparisons <- select(current_comparisons, -c('id', 'dataset_acc', 'title', 'summary', 'disease',
                                                          'organism', 'source', 'case_ann', 'case_sample',
                                                          'control_ann', 'control_sample','platform','cell_type'))
    
    
    ### Generate a table that displays the proportion of significant comparisons for each experiment type ###
    metrics_df_2 <- data.frame("Disease" = rep(current_disease, length(genes)), "gene" = genes)
    # Iterate over every experiment type
    for (exp_type in experiment_types) {
      
      # Pull the data associated with the given experiment type
      gene_data <- current_comparisons %>%
        filter(experiment_type == exp_type) %>%
        select(-c("experiment_type"))
      
      metrics_df_2[exp_type] <- ""
      for (gene in genes) {
        # Handle the case when there are no datasets associated with the given disease and experiment type
        if (nrow(gene_data) == 0) {
          metrics_df_2[metrics_df_2$gene == gene, exp_type] <- "NA"
          next
        }
        
        num_significant <- sum(gene_data[[gene]] == "info")
        num_total <- sum(gene_data[[gene]] == "info" | gene_data[[gene]] == "NS")
        percent_signifcant <- round(num_significant/num_total * 100, 2)
        
        # If the given gene was not found in any comparisons, the associated entry in the table is NA
        if (num_total == 0) {
          metrics_df_2[metrics_df_2$gene == gene, exp_type] <- "NA"
        } else {
          metrics_df_2[metrics_df_2$gene == gene, exp_type] <- paste0(num_significant, "/", num_total, " (", percent_signifcant,"%)")
        }
        
      }
      
    }
    
    # Append the metrics dataframe to the list
    metrics_df_list_2 <- append(metrics_df_list_2, list(metrics_df_2))
    ###
    
  }
  #================
  
  
  # Merge the metrics dataframes
  merged_metrics_df_1 <- do.call("rbind", metrics_df_list_1)
  
  # Create a disease factor used to group the table
  row_headers <- factor(merged_metrics_df_1$Disease, unique(merged_metrics_df_1$Disease))
  merged_metrics_df_1 <- select(merged_metrics_df_1, -c("Disease"))
  # Create the table
  merged_metrics_df_1 <- merged_metrics_df_1 %>% 
    kable("html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(font_size = 12) %>%
    row_spec(0, font_size = 14, color = "white", background = "#2C3E4C") %>%
    row_spec(1:nrow(merged_metrics_df_1), background = "white") %>%
    column_spec(1, extra_css = "border-left: 1px solid black;") %>%
    column_spec(ncol(merged_metrics_df_1), extra_css = "border-right: 1px solid black;") %>%
    pack_rows(index = table(row_headers), background = '#e1e2e3', 
              label_row_css = "border-top: 2px solid black; border-left: 1px solid black; border-right: 1px solid black;") %>%
    row_spec(nrow(merged_metrics_df_1), extra_css = "border-bottom: 1px solid black;") %>%
    footnote(general = paste("The entries under the experiment type columns indicate the number of datasets that",
                             "contain at least one significant comparison out of the total number of datasets for",
                             "which the given gene was found for the given disease and experiment type.",
                             "\nFormat: <# significant>/<# total> (% significant)"))
  
  
  # Merge the metrics dataframes
  merged_metrics_df_2 <- do.call("rbind", metrics_df_list_2)
  
  # Create a disease factor used to group the table
  row_headers <- factor(merged_metrics_df_2$Disease, unique(merged_metrics_df_2$Disease))
  merged_metrics_df_2 <- select(merged_metrics_df_2, -c("Disease"))
  # Create the table
  merged_metrics_df_2 <- merged_metrics_df_2 %>% 
    kable("html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(font_size = 12) %>%
    row_spec(0, font_size = 14, color = "white", background = "#2C3E4C") %>%
    row_spec(1:nrow(merged_metrics_df_2), background = "white") %>%
    column_spec(1, extra_css = "border-left: 1px solid black;") %>%
    column_spec(ncol(merged_metrics_df_2), extra_css = "border-right: 1px solid black;") %>%
    pack_rows(index = table(row_headers), background = '#e1e2e3', 
              label_row_css = "border-top: 2px solid black; border-left: 1px solid black; border-right: 1px solid black;") %>%
    row_spec(nrow(merged_metrics_df_2), extra_css = "border-bottom: 1px solid black;") %>%
    footnote(general = paste("The entries under the experiment type columns indicate the number of significant",
                             "comparisons out of the total number of comparisons for which the given gene was",
                             "found for the given disease and experiment type.",
                             "\nFormat: <# significant>/<# total> (% significant)"))
  
  return(list("dataset_results_summary" = merged_metrics_df_1, "comparison_results_summary" = merged_metrics_df_2))
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate summary tables for the summary section
# Table 1 displays the proportion of datasets that contain at least one significant comparison per gene per disease
# Table 2 displays the proportion of comparisons that are significant per gene per disease
#==============================================================================================================
signature_generate_results_summaries <- function(diseases, experiment_types, datasets, comparisons, signature_data, pval_threshold, db) {
  
  metrics_df_list_1 <- list()
  metrics_df_list_2 <- list()
  
  # If experiment_type is NULL, set experiment_type to all experiment types in the DB
  if (is.null(experiment_types)) {
    experiment_types <- get_keywords("ExperimentType", db)
  }
  
  #================
  # Iterate over every disease
  for (current_disease in diseases) {
    # Filter the datasets for the given disease
    current_datasets <- filter(datasets, str_detect(disease, current_disease))
    
    # If no datasets match the user specified queries, specify NA in the metrics dataframes and continue to the next iteration
    if (nrow(current_datasets) == 0) {
      metrics_df_1 <- data.frame("Disease" = current_disease, "signature" = "signature")
      for (type in experiment_types) {
        metrics_df_1[type] <- "NA"
      }
      metrics_df_list_1 <- append(metrics_df_list_1, list(metrics_df_1))
      
      metrics_df_2 <- data.frame("Disease" = current_disease, "signature" = "signature")
      for (type in experiment_types) {
        metrics_df_2[type] <- "NA"
      }
      metrics_df_list_2 <- append(metrics_df_list_2, list(metrics_df_2))
      
      next
    }
    
    # Filter the comparisons and comparison data
    current_comparisons <- filter(comparisons, dataset_acc %in% current_datasets$dataset_acc)
    current_comparison_ids <- comparisons %>%
      filter(dataset_acc %in% current_datasets$dataset_acc) %>%
      pull(id)
    current_signature_data <- filter(signature_data, comparison_id %in% current_comparison_ids)
    
    # If no comparisons match the user specified queries, specify NA in the metrics dataframes and continue to the next iteration
    if (is.null(current_signature_data)) {
      metrics_df_1 <- data.frame("Disease" = current_disease, "signature" = "signature")
      for (type in experiment_types) {
        metrics_df_1[type] <- "NA"
      }
      metrics_df_list_1 <- append(metrics_df_list_1, list(metrics_df_1))
      
      metrics_df_2 <- data.frame("Disease" = current_disease, "signature" = "signature")
      for (type in experiment_types) {
        metrics_df_2[type] <- "NA"
      }
      metrics_df_list_2 <- append(metrics_df_list_2, list(metrics_df_2))
      
      next
    }
    
    # Combine the data description, comparison description, and signature data
    current_comparisons <- current_comparisons %>%
      left_join(current_datasets, by = c("dataset_acc" = "dataset_acc")) %>% 
      left_join(current_signature_data, by = c("id" = "comparison_id", "dataset_acc" = "dataset_acc"))

    # Determine the significance status across all selected comparisons for the given signature 
    current_comparisons$signature <- sapply(current_comparisons$id, function(current_id) {
      comparison_info <- filter(current_comparisons, id == current_id)
      if (is.na(comparison_info$num_genes) | comparison_info$num_genes != 0) {
        # Determine if the given gene passes the user specified thresholds (i.e. p-value)
        significance_status <- comparison_info$p_value <= pval_threshold
        return(ifelse(significance_status, "info", "NS"))
      } else {
        return("NA")
      }
    })
    
    
    ### Generate a table that displays the proportion of datasets that contain ###
    ### at least one significant comparison for each experiment type ###
    metrics_df_1 <- data.frame("Disease" = current_disease, "signature" = "signature")
    # Iterate over every experiment type
    for (exp_type in experiment_types) {
      
      # Pull the data associated with the given experiment type
      gene_data <- current_comparisons %>% 
        filter(experiment_type == exp_type) %>%
        select(-c("experiment_type"))
      
      metrics_df_1[[exp_type]] <- ""

      # Handle the case when there are no datasets associated with the given disease and experiment type
      if (nrow(gene_data) == 0) {
        metrics_df_1[[exp_type]] <- "NA"
        next
      }
      
      num_significant <- sapply(unique(gene_data$dataset_acc), function(x) {
        dataset_data <- gene_data[gene_data$dataset_acc == x,]
        return(any(dataset_data$signature == 'info'))
      })
      num_significant = sum(num_significant)
      
      num_total <- sapply(unique(gene_data$dataset_acc), function(x) {
        dataset_data <- gene_data[gene_data$dataset_acc == x,]
        return(any(dataset_data$signature == 'info' | dataset_data$signature == 'NS'))
      })
      num_total <- sum(num_total)
      
      percent_signifcant <- round(num_significant/num_total * 100, 2)
      
      # If the given gene was not found in any datasets, the associated entry in the table is NA
      if (num_total == 0) {
        metrics_df_1[[exp_type]] <- "NA"
      } else {
        metrics_df_1[[exp_type]] <- paste0(num_significant, "/", num_total, " (", percent_signifcant,"%)")
      }

    }
    
    # Append the metrics dataframe to the list
    metrics_df_list_1 <- append(metrics_df_list_1, list(metrics_df_1))
    ###
    
    
    # Remove the unneeded columns from the comparisons
    current_comparisons <- select(current_comparisons, -c('id', 'dataset_acc', 'title', 'summary', 'disease',
                                                          'organism', 'source', 'case_ann', 'case_sample',
                                                          'control_ann', 'control_sample'))
    
    
    ### Generate a table that displays the proportion of significant comparisons for each experiment type ###
    metrics_df_2 <- data.frame("Disease" = current_disease, "signature" = "signature")
    # Iterate over every experiment type
    for (exp_type in experiment_types) {
      
      # Pull the data associated with the given experiment type
      gene_data <- current_comparisons %>%
        filter(experiment_type == exp_type) %>%
        select(-c("experiment_type"))
      
      metrics_df_2[[exp_type]] <- ""
      
      # Handle the case when there are no datasets associated with the given disease and experiment type
      if (nrow(gene_data) == 0) {
        metrics_df_2[[exp_type]] <- "NA"
        next
      }
      
      num_significant <- sum(gene_data$signature == "info")
      num_total <- sum(gene_data$signature == "info" | gene_data$signature == "NS")
      percent_signifcant <- round(num_significant/num_total * 100, 2)
      
      # If the given gene was not found in any comparisons, the associated entry in the table is NA
      if (num_total == 0) {
        metrics_df_2[[exp_type]] <- "NA"
      } else {
        metrics_df_2[[exp_type]] <- paste0(num_significant, "/", num_total, " (", percent_signifcant,"%)")
      }
      
    }
    
    # Append the metrics dataframe to the list
    metrics_df_list_2 <- append(metrics_df_list_2, list(metrics_df_2))
    ###
    
  }
  #================
  
  
  # Merge the metrics dataframes
  merged_metrics_df_1 <- do.call("rbind", metrics_df_list_1)
  
  # Create a disease factor used to group the table
  row_headers <- factor(merged_metrics_df_1$Disease, unique(merged_metrics_df_1$Disease))
  merged_metrics_df_1 <- select(merged_metrics_df_1, -c("Disease"))
  # Create the table
  merged_metrics_df_1 <- merged_metrics_df_1 %>% 
    kable("html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(font_size = 12) %>%
    row_spec(0, font_size = 14, color = "white", background = "#2C3E4C") %>%
    row_spec(1:nrow(merged_metrics_df_1), background = "white") %>%
    column_spec(1, extra_css = "border-left: 1px solid black;") %>%
    column_spec(ncol(merged_metrics_df_1), extra_css = "border-right: 1px solid black;") %>%
    pack_rows(index = table(row_headers), background = '#e1e2e3', 
              label_row_css = "border-top: 2px solid black; border-left: 1px solid black; border-right: 1px solid black;") %>%
    row_spec(nrow(merged_metrics_df_1), extra_css = "border-bottom: 1px solid black;") %>%
    footnote(general = paste("The entries under the experiment type columns indicate the number of datasets that",
                             "contain at least one significant comparison out of the total number of datasets for",
                             "which the given gene was found for the given disease and experiment type.",
                             "\nFormat: <# significant>/<# total> (% significant)"))
  
  
  # Merge the metrics dataframes
  merged_metrics_df_2 <- do.call("rbind", metrics_df_list_2)
  
  # Create a disease factor used to group the table
  row_headers <- factor(merged_metrics_df_2$Disease, unique(merged_metrics_df_2$Disease))
  merged_metrics_df_2 <- select(merged_metrics_df_2, -c("Disease"))
  # Create the table
  merged_metrics_df_2 <- merged_metrics_df_2 %>% 
    kable("html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(font_size = 12) %>%
    row_spec(0, font_size = 14, color = "white", background = "#2C3E4C") %>%
    row_spec(1:nrow(merged_metrics_df_2), background = "white") %>%
    column_spec(1, extra_css = "border-left: 1px solid black;") %>%
    column_spec(ncol(merged_metrics_df_2), extra_css = "border-right: 1px solid black;") %>%
    pack_rows(index = table(row_headers), background = '#e1e2e3', 
              label_row_css = "border-top: 2px solid black; border-left: 1px solid black; border-right: 1px solid black;") %>%
    row_spec(nrow(merged_metrics_df_2), extra_css = "border-bottom: 1px solid black;") %>%
    footnote(general = paste("The entries under the experiment type columns indicate the number of significant",
                             "comparisons out of the total number of comparisons for which the given gene was",
                             "found for the given disease and experiment type.",
                             "\nFormat: <# significant>/<# total> (% significant)"))
  
  return(list("dataset_results_summary" = merged_metrics_df_1, "comparison_results_summary" = merged_metrics_df_2))
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate summary tables for the summary section
# Table 1 displays the proportion of datasets that contain at least one significant comparison per gene per disease
# Table 2 displays the proportion of comparisons that are significant per gene per disease
#==============================================================================================================
pathway_generate_results_summaries <- function(pathways, diseases, experiment_types, datasets, comparisons, pathway_data, pval_option, pval_threshold, es_option, es_threshold, db) {
  
  datasets <- select(datasets, -id)
  
  metrics_df_list_1 <- list()
  metrics_df_list_2 <- list()
  
  # Determine if the user specified p-value or adjusted p-value under the SIGNIFICANCE THRESHOLDS option
  pval_col <- pval_option
  # Determine if the user specified ES or NES under the SIGNIFICANCE THRESHOLDS option
  enrichment_score_col <- es_option
  
  # If experiment_type is NULL, set experiment_type to all experiment types in the DB
  if (is.null(experiment_types)) {
    experiment_types <- get_keywords("ExperimentType", db)
  }
  
  #================
  # Iterate over every disease
  for (current_disease in diseases) {
    # Filter the datasets for the given disease
    current_datasets <- filter(datasets, str_detect(disease, current_disease))
    
    # If no datasets match the user specified queries, specify NA in the metrics dataframes and continue to the next iteration
    if (nrow(current_datasets) == 0) {
      metrics_df_1 <- data.frame("Disease" = rep(current_disease, length(pathways)), "pathway" = pathways)
      for (type in experiment_types) {
        metrics_df_1[type] <- "NA"
      }
      metrics_df_list_1 <- append(metrics_df_list_1, list(metrics_df_1))
      
      metrics_df_2 <- data.frame("Disease" = rep(current_disease, length(pathways)), "pathway" = pathways)
      for (type in experiment_types) {
        metrics_df_2[type] <- "NA"
      }
      metrics_df_list_2 <- append(metrics_df_list_2, list(metrics_df_2))
      
      next
    }
    # print(current_datasets$dataset_acc)
    # print(comparisons)
    # print(pathway_data)
    # Filter the comparisons and comparison data
    current_comparisons <- filter(comparisons, dataset_acc %in% current_datasets$dataset_acc)
    current_comparison_ids <- comparisons %>%
      filter(dataset_acc %in% current_datasets$dataset_acc) %>%
      pull(id)
    current_pathway_data <- filter(pathway_data, comparison_id %in% current_comparison_ids)
    
    # If no comparisons match the user specified queries, specify NA in the metrics dataframes and continue to the next iteration
    if (is.null(current_pathway_data)) {
      metrics_df_1 <- data.frame("Disease" = rep(current_disease, length(pathways)), "gene" = pathways)
      for (type in experiment_types) {
        metrics_df_1[type] <- "NA"
      }
      metrics_df_list_1 <- append(metrics_df_list_1, list(metrics_df_1))
      
      metrics_df_2 <- data.frame("Disease" = rep(current_disease, length(pathways)), "gene" = pathways)
      for (type in experiment_types) {
        metrics_df_2[type] <- "NA"
      }
      metrics_df_list_2 <- append(metrics_df_list_2, list(metrics_df_2))
      
      next
    }
    
    # Combine the data description and the comparison description
    current_comparisons <- left_join(current_comparisons, current_datasets, by = c("dataset_acc" = "dataset_acc"))
    
    # Add gene columns to the comparison description table 
    # Initialize values to NA - some of these will be changed later
    for (pathway in pathways) {
      current_comparisons[pathway] <- "NA"
    }
    
    # Iterate over every pathway
    for (current_pathway in pathways) {
      
      # Determine the significance status across all selected comparisons for the given pathway 
      current_comparisons[current_pathway] <- sapply(current_comparisons$id, function(current_id) {
        included_pathways <- current_pathway_data %>%
          filter(comparison_id == current_id) %>%
          pull(description)
        
        # Determine if the given comparison contains data on the given gene
        if (current_pathway %in% included_pathways) {
          comparison_info <- filter(current_pathway_data, comparison_id == current_id & description == current_pathway)[1,]
          # Determine if the given gene passes the user specified thresholds (i.e. p-value, q-value, logFC)
          significance_status <- comparison_info[[pval_col]] <= pval_threshold & abs(comparison_info[[enrichment_score_col]]) >= es_threshold
          return(ifelse(significance_status, "info", "NS"))
        } else {
          return("NA")
        }
      })
      
    }
    
    
    ### Generate a table that displays the proportion of datasets that contain ###
    ### at least one significant comparison for each experiment type ###
    metrics_df_1 <- data.frame("Disease" = rep(current_disease, length(pathways)), "pathway" = pathways)
    # Iterate over every experiment type
    for (exp_type in experiment_types) {
      
      # Pull the data associated with the given experiment type
      pathway_info <- current_comparisons %>% 
        filter(experiment_type == exp_type) %>%
        select(-c("experiment_type"))
      
      metrics_df_1[exp_type] <- ""
      for (pathway in pathways) {
        # Handle the case when there are no datasets associated with the given disease and experiment type
        if (nrow(pathway_info) == 0) {
          metrics_df_1[metrics_df_1$pathway == pathway, exp_type] <- "NA"
          next
        }
        
        num_significant <- sapply(unique(pathway_info$dataset_acc), function(x) {
          dataset_data <- pathway_info[pathway_info$dataset_acc == x,]
          return(any(dataset_data[[pathway]] == 'info'))
        })
        num_significant = sum(num_significant)
        
        num_total <- sapply(unique(pathway_info$dataset_acc), function(x) {
          dataset_data <- pathway_info[pathway_info$dataset_acc == x,]
          return(any(dataset_data[[pathway]] == 'info' | dataset_data[[pathway]] == 'NS'))
        })
        num_total <- sum(num_total)
        
        percent_signifcant <- round(num_significant/num_total * 100, 2)
        
        # If the given gene was not found in any datasets, the associated entry in the table is NA
        if (num_total == 0) {
          metrics_df_1[metrics_df_1$pathway == pathway, exp_type] <- "NA"
        } else {
          metrics_df_1[metrics_df_1$pathway == pathway, exp_type] <- paste0(num_significant, "/", num_total, " (", percent_signifcant,"%)")
        }
        
      }
      
    }
    
    # Append the metrics dataframe to the list
    metrics_df_list_1 <- append(metrics_df_list_1, list(metrics_df_1))
    ###
    
    
    # Remove the unneeded columns from the comparisons
    current_comparisons <- select(current_comparisons, -c('id', 'dataset_acc', 'title', 'summary', 'disease',
                                                          'organism', 'source', 'case_ann', 'case_sample',
                                                          'control_ann', 'control_sample','platform','cell_type'))
    
    
    ### Generate a table that displays the proportion of significant comparisons for each experiment type ###
    metrics_df_2 <- data.frame("Disease" = rep(current_disease, length(pathways)), "pathway" = pathways)
    # Iterate over every experiment type
    for (exp_type in experiment_types) {
      
      # Pull the data associated with the given experiment type
      pathway_info <- current_comparisons %>%
        filter(experiment_type == exp_type) %>%
        select(-c("experiment_type"))
      
      metrics_df_2[exp_type] <- ""
      for (pathway in pathways) {
        # Handle the case when there are no datasets associated with the given disease and experiment type
        if (nrow(pathway_info) == 0) {
          metrics_df_2[metrics_df_2$pathway == pathway, exp_type] <- "NA"
          next
        }
        
        num_significant <- sum(pathway_info[[pathway]] == "info")
        num_total <- sum(pathway_info[[pathway]] == "info" | pathway_info[[pathway]] == "NS")
        percent_signifcant <- round(num_significant/num_total * 100, 2)
        
        # If the given gene was not found in any comparisons, the associated entry in the table is NA
        if (num_total == 0) {
          metrics_df_2[metrics_df_2$pathway == pathway, exp_type] <- "NA"
        } else {
          metrics_df_2[metrics_df_2$pathway == pathway, exp_type] <- paste0(num_significant, "/", num_total, " (", percent_signifcant,"%)")
        }
        
      }
      
    }
    
    # Append the metrics dataframe to the list
    metrics_df_list_2 <- append(metrics_df_list_2, list(metrics_df_2))
    ###
    
  }
  #================
  
  
  # Merge the metrics dataframes
  merged_metrics_df_1 <- do.call("rbind", metrics_df_list_1)
  
  # Create a disease factor used to group the table
  row_headers <- factor(merged_metrics_df_1$Disease, unique(merged_metrics_df_1$Disease))
  merged_metrics_df_1 <- select(merged_metrics_df_1, -c("Disease"))
  # Create the table
  merged_metrics_df_1 <- merged_metrics_df_1 %>% 
    kable("html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(font_size = 12) %>%
    row_spec(0, font_size = 14, color = "white", background = "#2C3E4C") %>%
    row_spec(1:nrow(merged_metrics_df_1), background = 'white') %>%
    column_spec(1, extra_css = "border-left: 1px solid black;") %>%
    column_spec(ncol(merged_metrics_df_1), extra_css = "border-right: 1px solid black;") %>%
    pack_rows(index = table(row_headers), background = '#e1e2e3', 
              label_row_css = "border-top: 2px solid black; border-left: 1px solid black; border-right: 1px solid black;") %>%
    row_spec(nrow(merged_metrics_df_1), extra_css = "border-bottom: 1px solid black;") %>%
    footnote(general = paste("The entries under the experiment type columns indicate the number of datasets that",
                             "contain at least one significant comparison out of the total number of datasets for",
                             "which the given gene was found for the given disease and experiment type.",
                             "\nFormat: <# significant>/<# total> (% significant)"))
  
  
  # Merge the metrics dataframes
  merged_metrics_df_2 <- do.call("rbind", metrics_df_list_2)
  
  # Create a disease factor used to group the table
  row_headers <- factor(merged_metrics_df_2$Disease, unique(merged_metrics_df_2$Disease))
  merged_metrics_df_2 <- select(merged_metrics_df_2, -c("Disease"))
  # Create the table
  merged_metrics_df_2 <- merged_metrics_df_2 %>% 
    kable("html", escape = FALSE, row.names = FALSE) %>%
    kable_styling(font_size = 12) %>%
    row_spec(0, font_size = 14, color = "white", background = "#2C3E4C") %>%
    row_spec(1:nrow(merged_metrics_df_2), background = 'white') %>%
    column_spec(1, extra_css = "border-left: 1px solid black;") %>%
    column_spec(ncol(merged_metrics_df_2), extra_css = "border-right: 1px solid black;") %>%
    pack_rows(index = table(row_headers), background = '#e1e2e3', 
              label_row_css = "border-top: 2px solid black; border-left: 1px solid black; border-right: 1px solid black;") %>%
    row_spec(nrow(merged_metrics_df_2), extra_css = "border-bottom: 1px solid black;") %>%
    footnote(general = paste("The entries under the experiment type columns indicate the number of significant",
                             "comparisons out of the total number of comparisons for which the given gene was",
                             "found for the given disease and experiment type.",
                             "\nFormat: <# significant>/<# total> (% significant)"))
  
  return(list("dataset_results_summary" = merged_metrics_df_1, "comparison_results_summary" = merged_metrics_df_2))
}
#==============================================================================================================
