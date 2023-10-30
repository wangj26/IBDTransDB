#==============================================================================================================
# Establish the connection with the database
# Returns the connection object
#==============================================================================================================
connectDB <- function(dbname, host, port, user, password) {
  db <- dbConnect(RPostgres::Postgres(),
                   dbname = dbname, 
                   host = host, 
                   port = port, 
                   user = user,
                   password = password)
  return(db)
}
#==============================================================================================================


#==============================================================================================================
# Query the keyword table for keywords associated with a given attribute
# Returns a dataframe containing the options for each attribute
#==============================================================================================================
get_keywords <- function(attribute, db) {
  attribute <- paste0("'", attribute, "'")
  keywords_data <- dbGetQuery(db,
                              stringr::str_interp(
                                paste("SELECT keyword, old_keyword",
                                      "FROM keyword",
                                      "WHERE keyword_type LIKE ${attribute}")))
  keywords_data <- keywords_data %>%
    mutate(old_keyword = ifelse(sapply(old_keyword, function(x) is.null(x)), keyword, old_keyword))
  keywords <- keywords_data$old_keyword
  #names(keywords) <- keywords_data$keyword
  
  return(keywords)
}
#==============================================================================================================


#==============================================================================================================
# Format user specified queries 
#==============================================================================================================
format_query_ExactMatch <- function(query) {
  if (is.null(query)) {
    return("'%'")
  } else if (length(query) == 1) {
    return(paste0("'", query, "'"))
  } else {
    return(paste0("'(", paste(query, collapse = "|"), ")'"))
  }
}

format_query_SoftMatch <- function(query) {
  if (is.null(query)) {
    return("'%'")
  } else if (length(query) == 1) {
    return(paste0("'%", query, "%'"))
  } else {
    return(paste0("'%(", paste(query, collapse = "|"), ")%'"))
  }
}
#==============================================================================================================


#==============================================================================================================
# Query the dataset table for diseases that match the user query
# Returns a list of diseases
#==============================================================================================================
get_diseases <- function(diseases, organisms, sources, treatments, experiment_types, cell_types, db) {

  # Import the filtered data description for each disease and append the dataframes
  datasets_full <- list()
  for (disease in diseases) {
    
    # Import the filtered data description
    datasets <- dbGetQuery(db,
                           stringr::str_interp(
                             paste("SELECT dataset_acc, disease",
                                   "FROM dataset")))
    # Filter the datasets
    if (!is.null(organisms)) {
      datasets <- filter(datasets, organism %in% organisms)
    }
    if (!is.null(sources)) {
      datasets <- filter(datasets, sapply(source, function(x) any(unlist(strsplit(x, ";")) %in% sources)))
    }
    if (!is.null(treatments)) {
      datasets <- filter(datasets, sapply(treatment, function(x) any(unlist(strsplit(x, ";")) %in% treatments)))
    }
    if (!is.null(experiment_types)) {
      datasets <- filter(datasets, experiment_type %in% experiment_types)
    }
    if (!is.null(cell_types)) {
      datasets <- filter(datasets, sapply(cell_type, function(x) any(unlist(strsplit(x, ";")) %in% cell_types)))
    }
    # Filter the datasets based on the selected diseases
    datasets <- datasets %>%
      filter(sapply(disease, function(x) any(unlist(strsplit(x, ";")) %in% diseases)))
    
    datasets_full <- append(datasets_full, list(datasets))
  }
  datasets_full <- bind_rows(datasets_full)
  
  # Sort the data description by disease and return the diseases as a vector
  diseases <- datasets_full %>% 
    arrange(disease) %>% 
    pull(disease) %>%
    unique
  
  return(diseases)
}
#==============================================================================================================


#==============================================================================================================
# Query the dataset table for the following matches
# Used for generating a summary of the user query and displaying datasets to choose from
# Returns a dataframe with dataset descriptions for datasets matching the user query 
#==============================================================================================================
get_datasets <- function(diseases, organisms, sources, treatments, experiment_types, platforms, cell_types, timepoints, db) {
  
  # Import the filtered data description
  datasets <- dbGetQuery(db,
                         stringr::str_interp(
                           paste("SELECT id, dataset_acc, disease, organism, experiment_type, source, treatment, title, summary, platform, cell_type, timepoint",
                                 "FROM dataset")))
  attributes <- list("disease" = diseases,
                     "source" = sources, 
                     "cell_type" = cell_types, 
                     "treatment" = treatments, 
                     "timepoint" = timepoints)
  print(timepoints)
  rows_to_keep <- lapply(names(attributes), function(x) {
    if (is.null(attributes[[x]])) {
      return(rep(TRUE, nrow(datasets)))
    } else {
      return(sapply(datasets[[x]], function(y) any(str_detect(y, attributes[[x]]))))
    }
  })
  rows_to_keep <- Reduce("&", rows_to_keep)
  
  # Filter the dataset table for datasets that match the values for the modified attributes
  datasets <- filter(datasets, rows_to_keep)
  print(dim(datasets))
  # Filter the datasets
  # if (!is.null(diseases)) {
  #   datasets <- filter(datasets, sapply(disease, function(x) any(unlist(strsplit(x, ";")) %in% diseases)))
  # }
  # 
  # if (!is.null(sources)) {
  #   datasets <- filter(datasets, sapply(source, function(x) any(unlist(strsplit(x, ";")) %in% sources)))
  # }
  # 
  # if (!is.null(treatments)) {
  #   datasets <- filter(datasets, sapply(treatment, function(x) any(unlist(strsplit(x, ";")) %in% treatments)))
  # }
  # 
  # if (!is.null(cell_types)) {
  #   datasets <- filter(datasets, sapply(cell_type, function(x) any(unlist(strsplit(x, ";")) %in% cell_types)))
  # }
  # 
  # if (!is.null(timepoints)) {
  #   datasets <- filter(datasets, sapply(timepoint, function(x) any(unlist(strsplit(x, ";")) %in% timepoints)))
  # }
  
  
  # Filter the datasets based on the selected diseases and sort the datasets
  # datasets <- datasets %>%
  #   dplyr::filter(sapply(disease, function(x) any(unlist(strsplit(x, ";")) %in% diseases))) %>%
  #   dplyr::arrange(organism, experiment_type, dataset_acc)
  
  return(datasets)
}

#==============================================================================================================


#==============================================================================================================
# Query the comparison table for the following matches
# Returns a dataframe with comparison descriptions
#==============================================================================================================
get_comparisons <- function(datasets, db) {
  # Format the data ID list and other features so that they can be used for the query
  formatted_data_acc <- paste0("('", paste(datasets$dataset_acc, collapse = "', '"), "')")
  
  # Import the comparisons
  comparisons <- dbGetQuery(db,
                            stringr::str_interp( 
                              paste("SELECT id, dataset_acc, case_ann, case_sample, control_ann, control_sample",
                                    "FROM comparison",
                                    "WHERE dataset_acc IN ${formatted_data_acc}")))
  # # Filter the comparisons
  # if (!is.null(timepoints)) {
  #   comparisons <- filter(comparisons, sapply(case_ann, function(x) any(unlist(strsplit(x, ";")) %in% timepoints)) | sapply(control_ann, function(x) any(unlist(strsplit(x, ";")) %in% timepoints)))
  # }
  # if (!is.null(cell_types)) {
  #   comparisons <- filter(comparisons, sapply(case_ann, function(x) any(unlist(strsplit(x, ";")) %in% cell_types)) | sapply(control_ann, function(x) any(unlist(strsplit(x, ";")) %in% cell_types)))
  # }
  
  # If no comparisons match the user specified queries, return None for both the comparisons and comparison data
  if (nrow(comparisons) == 0) {
    return(NULL)
  }
  #print(unique(comparisons$dataset_acc))
  # Filter out comparisons not associated with any of the specified diseases, sources, or treatments
  comparisons <- mutate(comparisons, annotation = paste0(case_ann, ";", control_ann))
  # rows_to_keep <- sapply(1:nrow(comparisons), function(x) {
  #   current_dataset_acc <- comparisons$dataset_acc[x]
  #   current_annotations <- unlist(strsplit(comparisons$annotation[x], ";"))
  #   
  #   # Determine if the comparison is relevant based on the specified diseases
  #   current_diseases <- datasets %>%
  #     filter(dataset_acc == current_dataset_acc) %>%
  #     pull(disease)
  #   
  #   current_diseases <- unlist(strsplit(current_diseases, ";"))
  # 
  #   disease_relevance <- ifelse(length(current_diseases) > 1, any(diseases %in% current_annotations), TRUE)
  #   # Determine if the comparison is relevant based on the specified sources
  #   current_sources <- datasets %>%
  #     filter(dataset_acc == current_dataset_acc) %>%
  #     pull(source)
  #   current_sources <- unlist(strsplit(current_sources, ";"))
  #   source_relevance <- ifelse(!is.null(sources) & length(current_sources) > 1, any(sources %in% current_annotations), TRUE)
  #   # Determine if the comparison is relevant based on the specified treatments
  #   current_treatments <- datasets %>%
  #     filter(dataset_acc == current_dataset_acc) %>%
  #     pull(treatment)
  #   current_treatments <- unlist(strsplit(current_treatments, ";"))
  #   treatment_relevance <- ifelse(!is.null(treatments) & length(current_treatments) > 1, any(treatments %in% current_annotations), TRUE)
  #   return(disease_relevance | source_relevance | treatment_relevance)
  # })
  
  # comparisons <- comparisons %>%
  #   filter(rows_to_keep) %>%
  #   dplyr::select(-annotation)
  # print(unique(comparisons$dataset_acc))
  return(comparisons)
}
#==============================================================================================================


#==============================================================================================================
# Query the comparison_data table for the following matches
# Returns a dataframe with comparison metrics
#==============================================================================================================
get_comparison_data <- function(comparisons, genes, gene_map, db) {
  # If no comparisons match the user specified queries, return None for both the comparisons and comparison data
  if (nrow(comparisons) == 0) {
    return(NULL)
  }
  
  # Format the primary key ID's and gene list so that they can be used for the query
  formatted_comparison_ids <- paste0("('", paste(comparisons$id, collapse = "', '"), "')")
  formatted_genes <- gene_map$id[match(genes, gene_map$gene)]
  formatted_genes <- paste0("(", paste(formatted_genes, collapse = ", "), ")")
  
  # Import the filtered comparison data
  comparison_data <- dbGetQuery(db,
                                stringr::str_interp(
                                  paste("SELECT comparison_id, gene, log_fc, p_value, p_value_adj",
                                        "FROM comparison_data",
                                        "WHERE (comparison_id IN ${formatted_comparison_ids}",
                                        "AND gene IN ${formatted_genes})")))
  # Convert the gene identifiers in the comparison data to gene names
  comparison_data <- comparison_data %>%
    merge(gene_map, by.x = "gene", by.y = "id", all.x = TRUE) %>%
    select(-gene) %>%
    rename("gene" = "gene.y") %>%
    relocate(comparison_id, gene)
  
  # If none of the specified genes are found for the selected comparisons, return None
  if (nrow(comparison_data) == 0) {
    return(NULL)
  }
  
  return(comparison_data)
}
#==============================================================================================================


#==============================================================================================================
# Query the dataset_data table for the given signature across matching datasets
# Returns a dataframe with expression data
#==============================================================================================================
get_signature_data <- function(genes, dataset_acc_list, db) {
  
  # Format the genes and dataset IDs so that they can be used for the query
  formatted_genes <- paste0("('", paste(genes, collapse = "', '"), "')")
  formatted_dataset_acc_list <- paste0("('", paste(dataset_acc_list, collapse = "', '"), "')")
  
  # Import the filtered expression data
  exp_data <- dbGetQuery(db,
                         stringr::str_interp(
                           paste("SELECT dataset_acc, gene, sample_acc, expr",
                                 "FROM dataset_data",
                                 "WHERE (dataset_acc IN ${formatted_dataset_acc_list}",
                                 "AND gene in ${formatted_genes})")))
  
  return(exp_data)
}
#==============================================================================================================


#==============================================================================================================
# Query the enrichment_results table for the selected comparisons and pathways
# Returns a dataframe with comparison metrics
#==============================================================================================================
get_pathway_data <- function(comparisons, pathway_database, db) {
  # If no comparisons match the user specified queries, return None for both the comparisons and comparison data
  if (nrow(comparisons) == 0) {
    return(NULL)
  }
  
  # Format the pathway database and comparison IDs so that they can be used for the query
  formatted_pathway_database <- paste0("'", pathway_database, "'")
  formatted_comparison_ids <- paste0("('", paste(comparisons$id, collapse = "', '"), "')")
  
  # Import the enrichment run ids corresponding to the selected comparisons and pathway database
  enrichment_runs <- dbGetQuery(db,
                                stringr::str_interp(
                                  paste("SELECT id, comparison_id",
                                        "FROM enrichment_run",
                                        "WHERE (comparison_id IN ${formatted_comparison_ids}",
                                        "AND database = ${formatted_pathway_database})")))
  enrichment_run_ids <- pull(enrichment_runs, id)
  
  # Format the enrichment run IDs so that they can be used for the query
  formatted_enrichment_ids <- paste0("('", paste(enrichment_run_ids, collapse = "', '"), "')")

  # Import the enrichment results corresponding to the selected comparisons
  pathway_data <- dbGetQuery(db,
                             stringr::str_interp(
                               paste("SELECT enrichment_run_id, description, size, leading_edge_number, es, nes, p_value, p_value_adj",
                                     "FROM enrichment_results",
                                     "WHERE enrichment_run_id IN ${formatted_enrichment_ids}")))
  
  # If none of the specified genes are found for the selected comparisons, return None for both the compDesc and compData dataframes
  if (nrow(pathway_data) == 0) {
    return(NULL)
  }
  
  # Append the comparison IDs to the enrichment results
  pathway_data <- pathway_data %>%
    left_join(enrichment_runs, by = c("enrichment_run_id" = "id")) %>%
    select(-c("enrichment_run_id"))
  
  return(pathway_data)
}
#==============================================================================================================


#==============================================================================================================
# Get the comparison data related to the selected comparisons (Used for expression plots in the visuals)
# Returns a dataframe with the comparison metrics
#==============================================================================================================
get_selected_comparison_results <- function(gene, gene_map, comparisons, db) {
  # Pull the dataset ID's and comparison ID's
  dataset_acc <- unique(comparisons$dataset_acc)
  comparison_ids <- comparisons$comparison_id
  
  # Format the dataset ID and comparison ID's so that they can be used for the query
  formatted_dataset_acc <- paste0("('", paste(dataset_acc, collapse = "', '"), "')")
  formatted_comparison_ids <- paste0("('", paste(comparison_ids, collapse = "', '"), "')")
  
  # Import the associated comparison descriptions
  selected_comparisons <- dbGetQuery(db,
                                     stringr::str_interp(
                                       paste("SELECT id, case_ann, control_ann, dataset_acc",
                                             "FROM comparison",
                                             "WHERE (dataset_acc IN ${formatted_dataset_acc}",
                                             "AND id IN ${formatted_comparison_ids})"))) 
  
  # Format the gene so that it can be used for the query
  formatted_gene <- gene_map$id[match(gene, gene_map$gene)]
  
  # Import the associated comparison results
  comparison_results <- dbGetQuery(db,
                                   stringr::str_interp(
                                     paste("SELECT comparison_id, p_value, p_value_adj",
                                           "FROM comparison_data",
                                           "WHERE (comparison_id IN ${formatted_comparison_ids}",
                                           "AND gene = ${formatted_gene})")))  
  
  # Join the selected comaprisons and comparison results dataframes
  combined_results <- left_join(selected_comparisons, comparison_results, 
                                by = c("id" = "comparison_id"))
  
  return(combined_results)
}
#==============================================================================================================


#==============================================================================================================
# Generate options for target selection option 2 (i.e. cytokines, TNF super family)
#==============================================================================================================
get_predefined_groups <- function(db) {
  # Import the gene groups
  gene_groups <- dbGetQuery(db,
                            stringr::str_interp(
                              paste("SELECT type",
                                    "FROM gene_type")))
  gene_groups <- unique(gene_groups$type)
  
  return(gene_groups)
}
#==============================================================================================================


#==============================================================================================================
# Get genes associated with selected groups if target selection option 2 is used
#==============================================================================================================
get_predefined_group_genes <- function(selected_groups, db) {
  # Format the selected group list so that it can be used for the query
  selected_groups <- paste0("('", paste(selected_groups, collapse = "', '"), "')")
  
  # Import the genes from selected groups
  genes <- dbGetQuery(db,
                      stringr::str_interp(
                        paste("SELECT gene",
                              "FROM gene_type",
                              "WHERE type IN ${selected_groups}")))
  genes <- unique(genes$gene)
  
  return(genes)
}
#==============================================================================================================


#==============================================================================================================
# Function used to pull the mapping of gene identifiers to gene names
#==============================================================================================================
get_gene_map <- function(db) {
  # Pull the gene map
  gene_map <- dbGetQuery(db,
                         stringr::str_interp(
                           paste("SELECT *",
                                 "FROM gene_map")))
  
  return(gene_map)
}
#==============================================================================================================
