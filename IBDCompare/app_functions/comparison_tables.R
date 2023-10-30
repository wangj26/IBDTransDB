#==============================================================================================================
# Function used to define the user query
# Returns a named list with the user specified options
#==============================================================================================================
define_user_query <- function(genes, diseases, organisms, sources, cell_types, treatments, experiment_types, platforms, pval_option, pval_threshold, logfc_threshold, normalization_methods, row_specification) {
  user_query <- list("genes" = genes,
                     "diseases" = diseases, 
                     "organisms" = organisms,
                     "sources" = sources,
                     "cell_types" = cell_types,
                     "treatments" = treatments,
                     "experiment_types" = experiment_types,
                     "platforms" = platforms,
                     "pval_option" = pval_option,
                     "pval_threshold" = pval_threshold,
                     "logfc_threshold" = logfc_threshold,
                     "normalization_methods" = normalization_methods,
                     "row_specification" = row_specification)
  
  return(user_query)
}
#==============================================================================================================


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
  
  # Remove the case and control sample columns and add the shared attributes and comparison columns
  comparisons <- comparisons %>%
    select(-c("case_ann", "control_ann")) %>%
    mutate("group" = shared_attributes, "comparison" = unique_comparisons)
  
  return(comparisons)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate a comparison table for each disease
# UNSTYLED
# Returns a list of dataframes
#==============================================================================================================
generate_comparison_tables <- function(genes, diseases, datasets, comparisons, comparison_data, pval_option, pval_threshold, logfc_threshold, remove_ns) {
  
  datasets <- select(datasets, -id)
  
  # Produce a set of table(s) for every user specified disease
  table_list <- list()
  
  # Format the case and control samples columns (these are replaced by the group and comparison columns)
  comparisons <- format_comparison(comparisons)
  
  # Pull the dataset ID's from the data description (these will be used to order the dataframe after the join)
  unique_dataset_acc <- unique(comparisons$dataset_acc)
  # Combine the data description and the comparison description
  comparisons <- comparisons %>%
    left_join(datasets, by = c("dataset_acc" = "dataset_acc")) %>% 
    arrange(factor(dataset_acc, levels = unique_dataset_acc))
  
  # Determine if the user specified p-value or adjusted p-value under the SIGNIFICANCE THRESHOLDS option
  pval_col <- pval_option
  
  ##########
  for (current_disease in diseases) {
    
    # Filter the comparisons for the given disease and add gene columns (initialized to NA)
    current_comparisons <- filter(comparisons, disease == current_disease)
    for (gene in genes) {
      current_comparisons[gene] <- "NA"
    }
    
    # Convert the data ID's to hyperlinks
    current_comparisons$dataset_acc <- unlist(lapply(1:nrow(current_comparisons), function(x) {
      dataset_acc <- current_comparisons$dataset_acc[x]
      experiment_type <- current_comparisons$experiment_type[x]
      if (experiment_type != "scRNASeq") {
        return(cell_spec(dataset_acc, link = paste0("https://abbviegrc.shinyapps.io/ibdexplore_dataset/?dataset_acc=", tolower(dataset_acc)), new_tab = TRUE))
      } else {
        return(cell_spec(dataset_acc, link = paste0("http://rconnect.abbvienet.com:8080/immunoexplore_singlecell_dataset/?dataset_acc=", tolower(dataset_acc)), new_tab = TRUE))
      }
    }))
    
    ### Iterate over every gene ###
    for (current_gene in genes) {
      
      # Iterate over every comparison
      for (current_id in current_comparisons$id) {
        included_genes <- comparison_data %>%
          filter(comparison_id == current_id) %>%
          pull(gene)
        
        # Determine if the given comparison contains data on the given gene
        if (current_gene %in% included_genes) {
          comparison_info <- filter(comparison_data, comparison_id == current_id & gene == current_gene)
          current_gene_data <- HTML(paste0("name: ", current_gene,
                                           "<br/> logFC: ", round(comparison_info$log_fc, 4),
                                           "<br/> p-value: ", round(comparison_info$p_value, 6),
                                           "<br/> adj p-value: ", round(comparison_info$p_value_adj, 6)))
          
          # Determine if the given gene passes the user specified thresholds (i.e. p-value, q-value, logFC)
          significance_status <- comparison_info[[pval_col]] <= pval_threshold & abs(comparison_info$log_fc) >= logfc_threshold
          if (significance_status) {
            # Add the popup for the given cell (blue background: positive logFC, red background: negative logFC)
            if (comparison_info$log_fc > 0) {
              current_comparisons[current_comparisons$id == current_id, current_gene] <- cell_spec("up",
                                                                                                   popover = current_gene_data,
                                                                                                   color = "#1f37f0",
                                                                                                   background = "#f2bebb")
            } else {
              current_comparisons[current_comparisons$id == current_id, current_gene] <- cell_spec("down",
                                                                                                   popover = current_gene_data,
                                                                                                   color = "#1f37f0",
                                                                                                   background = "#bad8f7")
            }
          } else {
            # Add the popup for the given cell
            current_comparisons[current_comparisons$id == current_id, current_gene] <- cell_spec("NS",
                                                                                                 popover = current_gene_data,
                                                                                                 color = "black")
          }
        } else {
          current_gene_data <- HTML(paste0("name: ", current_gene,
                                           "<br/> No associated data"))
          # Add the popup for the given cell
          current_comparisons[current_comparisons$id == current_id, current_gene] <- cell_spec("NA",
                                                                                               popover = current_gene_data,
                                                                                               color = "#db160f")
        }
      }
      
    }
    ##########
    
    # Determine the rows that contain all NA's, all NS's, or all of both
    # This is used to determine which rows to remove based on the user specified option
    gene_df <- select(current_comparisons, all_of(genes))
    row_labels <- apply(gene_df, 1, function(x) {
      row_text <- unname(sapply(x, function(y) ifelse(y == "NA", "NA", str_match(y, ">\\s*([a-zA-Z]+?)\\s*</span>$")[2])))
      if (all(row_text == 'NA')) {
        return("all NA")
      } else if (all(row_text == 'NA' | row_text == 'NS')) {
        return("all NA or NS")
      } else {
        return("keep")
      }
    })
    
    
    ### Remove rows based on the user specified preferences ###
    if (remove_ns) {
      rows_to_keep <- row_labels == "keep"
    } else {
      rows_to_keep <- !(row_labels == "all NA")
    }
    current_comparisons <- filter(current_comparisons, rows_to_keep)
    ### IF ALL ROWS ARE REMOVED, assign NULL for the given disease in the table_list ###
    if (nrow(current_comparisons) == 0) {
      table_list[[current_disease]] <- NULL
      next
    }
    
    # Pull the text from the dataset ID column (these are now hyperlinks so we need the text)
    id_text <- sapply(current_comparisons$dataset_acc, function(x) str_match(x, ">\\s*(.+?)\\s*</a>$")[2])
    # Add a column with checkboxes (this is used to create selections for filtering rows and plotting)
    checkbox_col <- paste0('<input type="checkbox" id="', id_text, "_", current_comparisons$id, '">')
    current_comparisons["selection"] <- checkbox_col
    # Remove the checkboxes from scRNA-seq datasets
    current_comparisons[current_comparisons$experiment_type == "scRNASeq", "selection"] <- ""
    
    # Assign the data to the table list
    table_list[[current_disease]] <- current_comparisons
  }
  ##########
  
  return(table_list)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate a comparison table for each disease
# UNSTYLED
# Returns a list of dataframes
#==============================================================================================================
signature_generate_comparison_tables <- function(diseases, datasets, comparisons, signature_data, pval_threshold, remove_ns) {
  
  # Produce a set of table(s) for every user specified disease
  table_list <- list()
  
  # Format the case and control samples columns (these are replaced by the group and comparison columns)
  comparisons <- format_comparison(comparisons)
    
  # Pull the dataset ID's from the data description (these will be used to order the dataframe after the join)
  unique_dataset_acc <- unique(comparisons$dataset_acc)
  # Combine the data description, comparison description, and signature data
  comparisons <- comparisons %>%
    left_join(datasets, by = c("dataset_acc" = "dataset_acc")) %>% 
    left_join(signature_data, by = c("id" = "comparison_id", "dataset_acc" = "dataset_acc")) %>%
    arrange(factor(dataset_acc, levels = unique_dataset_acc))

  
  ##########
  for (current_disease in diseases) {
    
    # Filter the comparisons for the given disease and add a column for the signature (initialized to NA)
    current_comparisons <- comparisons %>%
      filter(disease == current_disease) %>%
      mutate(signature = "NA")
    
    # Convert the data ID's to hyperlinks
    current_comparisons$dataset_acc <- lapply(1:nrow(current_comparisons), function(x) {
      dataset_acc <- current_comparisons$dataset_acc[x]
      experiment_type <- current_comparisons$experiment_type[x]
      if (experiment_type != "scRNASeq") {
        return(cell_spec(dataset_acc, link = paste0("http://rconnect.abbvienet.com:8080/immunoexplore_dataset/?dataset_acc=", tolower(dataset_acc)), new_tab = TRUE))
      } else {
        return(cell_spec(dataset_acc, link = paste0("http://rconnect.abbvienet.com:8080/immunoexplore_singlecell_dataset/?dataset_acc=", tolower(dataset_acc)), new_tab = TRUE))
      }
    })
    
      
    ### Iterate over every comparison ###
    for (current_id in current_comparisons$id) {

      comparison_info <- filter(current_comparisons, id == current_id)
      
      # Determine if the given comparison contains data on the given gene
      if (comparison_info$num_genes != 0) {
        current_comparison_data <- HTML(paste0("logFC: ", signif(comparison_info$log_fc, 4),
                                               "<br/> p-value: ", signif(comparison_info$p_value, 6),
                                               "<br/> number of genes: ", comparison_info$num_genes))
        
        # Determine if the given gene passes the user specified thresholds (i.e. p-value)
        significance_status <- comparison_info$p_value <= pval_threshold
        if (significance_status) {
          # Add the popup for the given cell (red background: positive logFC, blue background: negative logFC)
          if (comparison_info$log_fc > 0) {
            current_comparisons[current_comparisons$id == current_id, "signature"] <- cell_spec("up",
                                                                                                popover = current_comparison_data,
                                                                                                color = "#1f37f0",
                                                                                                background = "#f2bebb")
          } else {
            current_comparisons[current_comparisons$id == current_id, "signature"] <- cell_spec("down",
                                                                                                popover = current_comparison_data,
                                                                                                color = "#1f37f0",
                                                                                                background = "#bad8f7")
          }
        } else {
          # Add the popup for the given cell
          current_comparisons[current_comparisons$id == current_id, "signature"] <- cell_spec("NS",
                                                                                              popover = current_comparison_data,
                                                                                              color = "black")
        }
      } else {
        current_comparison_data <- HTML("No associated data")
        # Add the popup for the given cell
        current_comparisons[current_comparisons$id == current_id, "signature"] <- cell_spec("NA",
                                                                                            popover = current_comparison_data,
                                                                                            color = "#db160f")
      }
    }
    current_comparisons <- select(current_comparisons, -c(log_fc, p_value, num_genes))
    ##########
    
    # Determine the rows that contain all NA's, all NS's, or all of both
    # This is used to determine which rows to remove based on the user specified option
    gene_df <- select(current_comparisons, signature)
    row_labels <- apply(gene_df, 1, function(x) {
      row_text <- unname(sapply(x, function(y) ifelse(y == "NA", "NA", str_match(y, ">\\s*([a-zA-Z]+?)\\s*</span>$")[2])))
      if (all(row_text == 'NA')) {
        return("all NA")
      } else if (all(row_text == 'NA' | row_text == 'NS')) {
        return("all NA or NS")
      } else {
        return("keep")
      }
    })
    
    
    ### Remove rows based on the user specified preferences ###
    if (remove_ns) {
      rows_to_keep <- row_labels == "keep"
    } else {
      rows_to_keep <- !(row_labels == "all NA")
    }
    current_comparisons <- filter(current_comparisons, rows_to_keep)
    ### IF ALL ROWS ARE REMOVED, assign NULL for the given disease in the table_list ###
    if (nrow(current_comparisons) == 0) {
      table_list[[current_disease]] <- NULL
      next
    }
    
    # Pull the text from the dataset ID column (these are now hyperlinks so we need the text)
    id_text <- sapply(current_comparisons$dataset_acc, function(x) str_match(x, ">\\s*(.+?)\\s*</a>$")[2])
    # Add a column with checkboxes (this is used to create selections for filtering rows and plotting)
    checkbox_col <- paste0('<input type="checkbox" id="', id_text,
                           "_", current_comparisons$id, '" value="', 
                           1:nrow(current_comparisons), '">',"")
    # current_comparisons["selection"] <- checkbox_col
    # # Remove the checkboxes from scRNA-seq datasets
    # current_comparisons[current_comparisons$experiment_type == "scRNASeq", "selection"] <- ""
    
    # Assign the data to the table list
    table_list[[current_disease]] <- current_comparisons
  }
  ##########
  
  return(table_list)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate a comparison table for each disease
# UNSTYLED
# Returns a list of dataframes
#==============================================================================================================
pathway_generate_comparison_tables <- function(pathways, diseases, datasets, comparisons, pathway_data, pval_option, pval_threshold, es_option, es_threshold, remove_ns) {
  
  datasets <- select(datasets, -id)
  
  # Produce a set of table(s) for every user specified disease
  table_list <- list()
  
  # Format the case and control samples columns (these are replaced by the group and comparison columns)
  comparisons <- format_comparison(comparisons)
  
  # Pull the dataset ID's from the data description (these will be used to order the dataframe after the join)
  unique_dataset_acc <- unique(comparisons$dataset_acc)
  # Combine the data description and the comparison description
  comparisons <- comparisons %>%
    left_join(datasets, by = c("dataset_acc" = "dataset_acc")) %>% 
    arrange(factor(dataset_acc, levels = unique_dataset_acc))
  
  # Determine if the user specified p-value or adjusted p-value under the SIGNIFICANCE THRESHOLDS option
  pval_col <- pval_option
  # Determine if the user specified ES or NES under the SIGNIFICANCE THRESHOLDS option
  if (es_option == "ES") {
    enrichment_score_col <- "es"
  } else {
    enrichment_score_col <- "nes"
  }
  
  ##########
  for (current_disease in diseases) {
    
    # Filter the comparisons for the given disease and add pathway columns (initialized to NA)
    current_comparisons <- filter(comparisons, disease == current_disease)
    for (pathway in pathways) {
      current_comparisons[pathway] <- "NA"
    }
    
    # Convert the data ID's to hyperlinks
    current_comparisons$dataset_acc <- unlist(lapply(1:nrow(current_comparisons), function(x) {
      dataset_acc <- current_comparisons$dataset_acc[x]
      experiment_type <- current_comparisons$experiment_type[x]
      if (experiment_type != "scRNASeq") {
        return(cell_spec(dataset_acc, link = paste0("https://abbviegrc.shinyapps.io/ibdexplore_dataset/?dataset_acc=", tolower(dataset_acc)), new_tab = TRUE))
      } else {
        return(cell_spec(dataset_acc, link = paste0("http://rconnect.abbvienet.com:8080/immunoexplore_singlecell_dataset/?dataset_acc=", tolower(dataset_acc)), new_tab = TRUE))
      }
    }))
    
    
    ### Iterate over every pathway ###
    for (current_pathway in pathways) {
      
      # Iterate over every comparison
      for (current_id in current_comparisons$id) {
        included_pathways <- pathway_data %>%
          filter(comparison_id == current_id) %>%
          pull(description)
        
        # Determine if the given comparison contains data on the given gene
        if (current_pathway %in% included_pathways) {
          comparison_info <- filter(pathway_data, comparison_id == current_id & description == current_pathway)
          current_pathway_data <- HTML(paste("name:", current_pathway,
                                             "<br/> ES:", round(comparison_info$es, 4),
                                             "<br/> NES:", round(comparison_info$nes, 4),
                                             "<br/> p-value:", round(comparison_info$p_value, 6),
                                             "<br/> adj p-value:", round(comparison_info$p_value_adj, 6)))
          
          # Determine if the given gene passes the user specified thresholds (i.e. p-value, q-value, ES, NES)
          significance_status <- comparison_info[[pval_col]] <= pval_threshold & abs(comparison_info[[enrichment_score_col]]) >= es_threshold
          if (significance_status) {
            # # Add the popup for the given cell (blue background: positive logFC, red background: negative logFC)
            if (comparison_info$es > 0) {
              current_comparisons[current_comparisons$id == current_id, current_pathway] <- cell_spec("up",
                                                                                                      popover = current_pathway_data,
                                                                                                      color = "#1f37f0",
                                                                                                      background = "#f2bebb")
            } else {
              current_comparisons[current_comparisons$id == current_id, current_pathway] <- cell_spec("down",
                                                                                                      popover = current_pathway_data,
                                                                                                      color = "#1f37f0",
                                                                                                      background = "#bad8f7")
            }
          } else {
            # Add the popup for the given cell
            current_comparisons[current_comparisons$id == current_id, current_pathway] <- cell_spec("NS",
                                                                                                    popover = current_pathway_data,
                                                                                                    color = "black")
          }
        } else {
          current_pathway_data <- HTML(paste0("name: ", current_pathway,
                                              "<br/> No associated data"))
          # Add the popup for the given cell
          current_comparisons[current_comparisons$id == current_id, current_pathway] <- cell_spec("NA",
                                                                                                  popover = current_pathway_data,
                                                                                                  color = "#db160f")
        }
      }
      
    }
    ##########
    
    # Determine the rows that contain all NA's, all NS's, or all of both
    # This is used to determine which rows to remove based on the user specified option
    pathway_df <- select(current_comparisons, all_of(pathways))
    row_labels <- apply(pathway_df, 1, function(x) {
      row_text <- unname(sapply(x, function(y) ifelse(y == "NA", "NA", str_match(y, ">\\s*([a-zA-Z]+?)\\s*</span>$")[2])))
      if (all(row_text == 'NA')) {
        return("all NA")
      } else if (all(row_text == 'NA' | row_text == 'NS')) {
        return("all NA or NS")
      } else {
        return("keep")
      }
    })
    
    
    ### Remove rows based on the user specified preferences ###
    if (remove_ns) {
      rows_to_keep <- row_labels == "keep"
    } else {
      rows_to_keep <- !(row_labels == "all NA")
    }
    current_comparisons <- filter(current_comparisons, rows_to_keep)
    ### IF ALL ROWS ARE REMOVED, assign NULL for the given disease in the table_list ###
    if (nrow(current_comparisons) == 0) {
      table_list[[current_disease]] <- NULL
      next
    }
    
    # Pull the text from the dataset ID column (these are now hyperlinks so we need the text)
    id_text <- sapply(current_comparisons$dataset_acc, function(x) str_match(x, ">\\s*(.+?)\\s*</a>$")[2])
    # Add a column with checkboxes (this is used to create selections for filtering rows and plotting)
    checkbox_col <- paste0('<input type="checkbox" id="', id_text,
                           "_", current_comparisons$id, '" value="', 
                           1:nrow(current_comparisons), '">',"")
    current_comparisons["selection"] <- checkbox_col
    # Remove the checkboxes from scRNA-seq datasets
    current_comparisons[current_comparisons$experiment_type == "scRNASeq", "selection"] <- ""
    
    # Assign the data to the table list
    table_list[[current_disease]] <- current_comparisons
  }
  ##########
  
  return(table_list)
}
#==============================================================================================================


#==============================================================================================================
# Function used to style the comparison tables produced by the generate_comparison_tables function
# STYLED
# Returns a list of kable objects
#==============================================================================================================
style_comparison_tables <- function(diseases, table_list) {
  
  # Iterate over every disease in order to style the associated dataframe
  for (disease in diseases) {
    
    # Pull the associated dataframe - If null, continue to the next iteration
    comparison_table <- table_list[[disease]]
    if (is.null(comparison_table)) {
      next
    }
    
    # Create headers including the title, summary, organism, type, and source for every data ID
    data_header_df <- select(comparison_table, c("dataset_acc", "title", "organism", "experiment_type", "source", "treatment","cell_type"))
    header_list <- c()
    for (i in 1:nrow(data_header_df)) {
      # Specify the data description for the given data ID
      study_data <- paste0("Title: ", gsub('"', '', data_header_df$title[i]),
                           "\nSource: ", data_header_df$source[i],
                           "\nCellType: ", data_header_df$cell_type[i],
                           "\nTreatment: ", data_header_df$treatment[i],
                           "\nTimepoints: ", ifelse(is.null(data_header_df$timepoint[i]),NA,data_header_df$timepoint[i])
                           )
      # Append the data to the header list
      header_list[i] <- study_data
    }
    
    # Remove the ID, title, and summary columns and format the column names
    comparison_table <- comparison_table %>%
      select(-c("id", "disease", "title", "summary", "organism", "experiment_type", "source", "treatment", "timepoint", "case_sample", "control_sample","platform","cell_type")) %>%
      rename("data ID" = "dataset_acc")
    
    # Generate the styled table
    data_info <- factor(header_list, unique(header_list))
    comparison_table <- comparison_table %>% kable('html', escape = FALSE, row.names = FALSE) %>%
      kable_styling(font_size = 12) %>%
      row_spec(0, font_size = 14, color = "white", background = "#2C3E4C") %>%
      row_spec(1:nrow(comparison_table), background = 'white') %>%
      column_spec(1, width = "15em", extra_css = "text-align: left;") %>%
      column_spec(2, width = "20em") %>%
      column_spec(3, width = "30em") %>%
      collapse_rows(1, valign = 'top') %>%
      pack_rows(index = table(data_info), background = "#E1E2E3", 
                label_row_css = "border-top: 2px solid black;")
    
    table_list[[disease]] <- comparison_table
  }
  
  return(table_list)
}
#==============================================================================================================


#==============================================================================================================
# Function used to filter the comparison tables produced by the generate_comp_tables_unstyled function
# Returns a list of dataframes
#==============================================================================================================
filter_comparison_tables <- function(diseases, checked_comparisons, filter_option, table_list) {
  
  # Determine all unique comparisons selected across all tables
  unique_comparisons <- lapply(diseases, function(disease) {
    # Pull the associated dataframe - If null, continue to the next iteration
    comparison_table <- table_list[[disease]]
    if (is.null(comparison_table)) {
      return(NULL)
    }
    
    # Determine the unique comparisons based on the user selections if the user wants to filter for similar comparisons
    selected_comparisons <- comparison_table %>%
      mutate(comparison_ids = str_match(selection,'id=\"(.*?)\" value')[,2]) %>%
      filter(comparison_ids %in% checked_comparisons) %>%
      pull("comparison") %>%
      unique()
    
    return(selected_comparisons)
  })
  unique_comparisons <- unlist(unique_comparisons)
  
  # Iterate over every disease and organism in order to filter the associated tables
  for (disease in diseases) {

    # Pull the associated dataframe - If null, continue to the next iteration
    comparison_table <- table_list[[disease]]
    if (is.null(comparison_table)) {
      next
    }
    
    if (filter_option) {
      # Filter the comparison table for the user selected comparisons
      comparison_table <- filter(comparison_table, comparison %in% unique_comparisons)
    } else {
      # Filter the comparison table for the user selected comparisons
      comparison_table <- comparison_table %>%
        mutate(comparison_ids = str_match(selection,'id=\"(.*?)\" value')[,2]) %>%
        filter(comparison_ids %in% checked_comparisons) %>%
        select(-c("comparison_ids"))
    }
    
    # Necessary for the style_comparison_tables function to operate properly 
    # when no comparisons where selected for a particular table
    if (nrow(comparison_table) == 0) {
      comparison_table <- NULL
    }
    
    table_list[[disease]] <- comparison_table
  }
  
  return(table_list)
}
#==============================================================================================================
