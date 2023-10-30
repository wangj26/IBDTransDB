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
# Query the dataset data table for all searchable genes
# Returns a list of diseases
#==============================================================================================================
get_genes <- function(db) {
  # Import the gene names
  genes <- dbGetQuery(db,
                      stringr::str_interp(
                        paste("SELECT gene",
                              "FROM dataset_data"))) 
  genes <- unique(genes$gene)
  
  return(genes)
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
# Produce a list of exact diseases from the datasets matching the user query
# Returns a list of diseases
#==============================================================================================================
get_exact_diseases <- function(datasets) {
  # Remove keywords Healthy and nonIBD and pull the unique diseases
  diseases <- datasets %>%
    arrange(disease) %>%
    pull(disease) %>%
    unique()
  
  return(diseases)
}
#==============================================================================================================


#==============================================================================================================
# Query the dataset table for datasets that match the user query
# Returns a dataframe with dataset descriptions
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
  print(class(dataset_table))
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
  
  return(dataset_table)
}

#==============================================================================================================
## GET DATASET FUNCTION BASED OFF IMMUNOCOMPARE
#==============================================================================================================
# get_datasets <- function(diseases, organisms, sources, treatments, experiment_types, platforms, cell_types, db, exact_disease) {
#   # Format the queries
#   organisms <- format_query_ExactMatch(organisms)
#   sources <- format_query_SoftMatch(sources)
#   treatments <- format_query_SoftMatch(treatments)
#   experiment_types <- format_query_ExactMatch(experiment_types)
#   platforms <- ifelse(is.null(platforms), 
#                       "(platform SIMILAR TO '%' OR platform IS NULL)", 
#                       paste("platform SIMILAR TO", format_query_ExactMatch(platforms)))
#   cell_types <- ifelse(is.null(cell_types), 
#                        "(cell_type SIMILAR TO '%' OR cell_type IS NULL)", 
#                        paste("cell_type SIMILAR TO", format_query_SoftMatch(cell_types)))
#   
#   # Determine whether to use an exact match or soft match for the disease
#   if (exact_disease) {
#     diseases <- format_query_ExactMatch(diseases)
#   } else {
#     diseases <- format_query_SoftMatch(diseases)
#   }
#   
#   # Import the filtered data description
#   datasets <- dbGetQuery(db,
#                          stringr::str_interp(
#                            paste("SELECT dataset_acc, disease, organism, experiment_type, source, treatment, title, summary",
#                                  "FROM dataset",
#                                  "WHERE (disease SIMILAR TO ${diseases}",
#                                  "AND organism SIMILAR TO ${organisms}",
#                                  "AND source SIMILAR TO ${sources}",
#                                  "AND treatment SIMILAR TO ${treatments}",
#                                  "AND experiment_type SIMILAR TO ${experiment_types}",
#                                  "AND ${platforms}",
#                                  "AND ${cell_types})")))
#   
#   # Sort the data description
#   datasets <- datasets %>%
#     arrange(organism, experiment_type, dataset_acc)
#   
#   return(datasets)
# }
#==============================================================================================================


#==============================================================================================================
# Pull the dataset descriptions displayed in the dataset selection table from the DB
# Returns a dataframe
#==============================================================================================================
get_dataset_table_data <- function(datasets, db) {
  # Format the dataset IDs for the query
  formatted_dataset_acc_list <- paste0("('", paste(datasets$dataset_acc, collapse = "', '"), "')")
  # Pull the datasets
  datasets <- dbGetQuery(db,
                         stringr::str_interp(
                           paste("SELECT dataset_acc, title, disease, organism, experiment_type, source, treatment, timepoint, dose, sample_number",
                                 "FROM dataset",
                                  " WHERE dataset_acc IN ${formatted_dataset_acc_list}"
                                 )))
  
  # Sort the data description by several features and format long features (i.e. add ellipsis)
  datasets <- datasets %>%
    arrange(disease, organism, experiment_type, dataset_acc) %>%
    mutate(title = str_trunc(title, width = 40, ellipsis = "..."),
           disease = str_trunc(disease, width = 25, ellipsis = "..."),
           source = str_trunc(source, width = 25, ellipsis = "..."),
           treatment = str_trunc(treatment, width = 25, ellipsis = "..."))
  
  return(datasets)
}
#==============================================================================================================


#==============================================================================================================
# Query the comparison table for the following matches
# Returns a dataframes with comparison descriptions
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
# Function used to query the sample_ann table for the sample annotations of the given dataset
#==============================================================================================================
get_sample_ann <- function(datasets, db) {
  
  # Format the dataset IDs so that they can be used for the query
  formatted_dataset_acc <- paste0("('", paste(datasets$dataset_acc, collapse = "', '"), "')")
  
  # Pull the sample annotations
  sample_ann <- dbGetQuery(db,
                           stringr::str_interp(
                             paste("SELECT id, sample_id, sample_ann_type, sample_ann_value",
                                   "FROM sample_ann",
                                   "WHERE dataset_acc IN ${formatted_dataset_acc}")))
  
  return(sample_ann)
}
#==============================================================================================================


#==============================================================================================================
# Query the comparison_data table for the following matches
# Returns a dataframe with comparison metrics
#==============================================================================================================
get_comparison_data <- function(comparisons, genes, gene_map, db) {
  # If no comparisons match the user specified queries, return None for both the comparisons and comparison data
  if (nrow(comparisons) == 0 | is.null(comparisons)) {
    return(NULL)
  }
  
  # Format the genes so that they can be used for the query
  formatted_genes <- gene_map$id[match(genes, gene_map$gene)]
  formatted_genes <- formatted_genes[!is.na(formatted_genes)]
  if(length(formatted_genes)==0){
    return(NULL)
  }
  formatted_genes <- paste0("(", paste(formatted_genes, collapse = ", "), ")")
  
  comparison_data <- lapply(comparisons$id, function(comparison_id) {
    # Format the primary key ID so that it can be used for the query
    formatted_comparison_id <- paste0("'", comparison_id, "'")
    
    # Import the filtered comparison data
    results <- dbGetQuery(db,
                          stringr::str_interp(
                            paste("SELECT comparison_id, gene, log_fc, p_value, p_value_adj",
                                  "FROM comparison_data",
                                  "WHERE (comparison_id = ${formatted_comparison_id}",
                                  "AND gene IN ${formatted_genes})")))
    
    return(results)
  })
  comparison_data <- do.call(rbind, comparison_data)
  
  # Convert the gene identifiers in the comparison data to gene names
  comparison_data <- comparison_data %>%
    merge(gene_map, by.x = "gene", by.y = "id", all.x = TRUE) %>%
    dplyr::select(-gene) %>%
    dplyr::rename("gene" = "gene.y") %>%
    relocate(comparison_id, gene)
  return(comparison_data)
}
#==============================================================================================================


#==============================================================================================================
# Query the enrichment_results table for the following matches
# Returns a dataframe with comparison metrics
#==============================================================================================================
get_pathway_data <- function(comparisons, pathway_databases, pathways, db) {
  # If no comparisons match the user specified queries, return None for both the comparisons and comparison data
  if (nrow(comparisons) == 0) {
    return(NULL)
  }
  
  # Format the pathway databases and comparison IDs so that they can be used for the query
  formatted_pathway_databases <- paste0("('", paste(pathway_databases, collapse = "', '"), "')")
  formatted_comparison_ids <- paste0("('", paste(comparisons$id, collapse = "', '"), "')")

  # Import the enrichment run ids corresponding to the selected comparisons
  enrichment_runs <- dbGetQuery(db,
                                stringr::str_interp(
                                  paste("SELECT id, comparison_id",
                                        "FROM enrichment_run",
                                        "WHERE (comparison_id IN ${formatted_comparison_ids}",
                                        "AND database IN ${formatted_pathway_databases})")))
  enrichment_run_ids <- pull(enrichment_runs, id)
  
  # Format the enrichment run IDs and pathway list so that they can be used for the query
  formatted_enrichment_ids <- paste0("('", paste(enrichment_run_ids, collapse = "', '"), "')")
  formatted_pathways <- paste0("('", paste(pathways, collapse = "', '"), "')")
  
  # Import the enrichment results corresponding to the selected comparisons
  pathway_data <- dbGetQuery(db,
                             stringr::str_interp(
                               paste("SELECT enrichment_run_id, description, size, leading_edge_number, es, nes, p_value, p_value_adj",
                                     "FROM enrichment_results",
                                     "WHERE (enrichment_run_id IN ${formatted_enrichment_ids}",
                                     "AND description IN ${formatted_pathways})")))
  
  # Append the comparison IDs to the enrichment results
  pathway_data <- pathway_data %>%
    left_join(enrichment_runs, by = c("enrichment_run_id" = "id")) %>%
    dplyr::select(-c("enrichment_run_id"))
  
  return(pathway_data)
}
#==============================================================================================================


#==============================================================================================================
# Query the dataset_data table for the given signature across matching datasets
# Returns a dataframe with expression data
#==============================================================================================================
# get_signature_data <- function(genes, gene_map, dataset_acc_list, sample_ann, db) {
#   
#   # Format the genes and dataset IDs so that they can be used for the query
#   formatted_genes <- gene_map$id[match(genes, gene_map$gene)]
#   formatted_genes <- paste0("(", paste(formatted_genes, collapse = ", "), ")")
#   formatted_dataset_acc_list <- paste0("('", paste(dataset_acc_list, collapse = "', '"), "')")
#   
#   # Import the filtered expression data
#   exp_data <- dbGetQuery(db,
#                          stringr::str_interp(
#                            paste("SELECT dataset_acc, gene, sample_acc, expr",
#                                  "FROM dataset_data",
#                                  "WHERE (dataset_acc IN ${formatted_dataset_acc_list}",
#                                  "AND gene in ${formatted_genes})")))
#   
#   # Convert the dataset, gene, and sample identifiers to their respective names
#   sample_ann <- select(id, sample_id)
#   exp_data <- exp_data %>%
#     mutate(dataset_acc = current_dataset_acc) %>%
#     merge(gene_map, by.x = "gene", by.y = "id", all.x = TRUE) %>%
#     select(-gene) %>%
#     rename("gene" = "gene.y") %>%
#     merge(sample_ann, by.x = "sample_acc", by.y = "id", all.x = TRUE) %>%
#     select(-sample_acc) %>%
#     rename("sample_acc" = "sample_id")
#   
#   return(exp_data)
# }
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
  comparison_results <- unique(dbGetQuery(db,
                                   stringr::str_interp(
                                     paste("SELECT comparison_id, p_value, p_value_adj",
                                           "FROM comparison_data",
                                           "WHERE (comparison_id IN ${formatted_comparison_ids}",
                                           "AND gene = ${formatted_gene})"))))
  
  # Join the selected comparisons and comparison results dataframes
  combined_results <- left_join(selected_comparisons, comparison_results, 
                                by = c("id" = "comparison_id"))
  
  return(combined_results)
}
#==============================================================================================================


#==============================================================================================================
# Get the expression data related to the selected comparisons (Used for expression plots in the visuals tab)
# Returns a dataframe with the expression data
#==============================================================================================================
get_selected_comparison_expression <- function(gene, gene_map, datasets, comparisons, sample_ann, db) {
  # Pull the dataset ID's and comparison ID's
  dataset_acc <- unique(sapply(comparisons, function(x) unlist(str_split(x, "_"))[1]))
  comparison_ids <- unique(sapply(comparisons, function(x) unlist(str_split(x, "_"))[2]))
  
  # Format the dataset ID and comparison ID's so that they can be used for the query
  formatted_dataset_acc <- paste0("('", paste(dataset_acc, collapse = "', '"), "')")
  formatted_comparison_ids <- paste0("('", paste(comparison_ids, collapse = "', '"), "')")
  
  # Import the associated comparison descriptions
  selected_comparisons <- dbGetQuery(db,
                                     stringr::str_interp(
                                       paste("SELECT case_ann, case_sample, control_ann, control_sample, dataset_acc",
                                             "FROM comparison",
                                             "WHERE (dataset_acc IN ${formatted_dataset_acc}",
                                             "AND id IN ${formatted_comparison_ids})")))
  # Format the associated comparisons so that the dataframe contains only two columns ("annotation" and "sample")
  case_df <- selected_comparisons %>% 
    dplyr::select(-c("control_ann", "control_sample")) %>%
    dplyr::rename(annotation = case_ann, sample = case_sample)
  control_df <- selected_comparisons %>%
    dplyr::select(-c('case_ann', 'case_sample')) %>%
    dplyr::rename(annotation = control_ann, sample = control_sample)
  selected_comparisons <- rbind(case_df, control_df)
  print("format")
  # Explode the dataframe so that each row corresponds to one sample
  selected_comparisons <- separate_rows(selected_comparisons, sample, sep = ";")
  # Keep the unique rows in the dataframe
  selected_comparisons <- selected_comparisons[!duplicated(selected_comparisons),]
  print(head(selected_comparisons))
  # Format the genes dataset IDs, and sample IDs so that they can be used for the query
  formatted_gene <- gene_map$id[match(gene, gene_map$gene)]
  formatted_dataset_acc <- datasets$id[match(dataset_acc, datasets$dataset_acc)]
  formatted_dataset_acc <- paste0("(", paste(formatted_dataset_acc, collapse = ", "), ")")
  formatted_samples <- sample_ann$id[match(selected_comparisons$sample, sample_ann$sample_id)]
  formatted_samples <- paste0("(", paste(formatted_samples, collapse = ", "), ")")
  print("formatted_samples")
  print(formatted_dataset_acc)
  print(formatted_samples)
  print(formatted_gene)
  # Import the associated dataset data based on the selected samples
  exp_data <- dbGetQuery(db,
                         stringr::str_interp(
                           paste("SELECT gene, sample_acc, expr",
                                 "FROM dataset_data",
                                 "WHERE (dataset_acc IN ${formatted_dataset_acc}",
                                 "AND gene = ${formatted_gene}",
                                 "AND sample_acc IN ${formatted_samples})")))
  
  print("exp_data_Retr")
  print(head(exp_data))
  # Convert the dataset, gene, and sample identifiers to their respective names
  
  sample_ann <- select(sample_ann, c(id, sample_id))
  exp_data <- exp_data %>%
    merge(gene_map, by.x = "gene", by.y = "id", all.x = TRUE) %>%
    dplyr::select(-gene) %>%
    dplyr::rename("gene" = "gene.y") %>%
    merge(sample_ann, by.x = "sample_acc", by.y = "id", all.x = TRUE) %>%
    dplyr::select(-sample_acc) %>%
    dplyr::rename("sample" = "sample_id")
  print(head(exp_data))
  # Join the selected comparisons and expression data dataframes
  combined_data <- inner_join(selected_comparisons, exp_data, by = c("sample" = "sample"))
  
  return(combined_data)
}
#==============================================================================================================


#==============================================================================================================
# List all pathways that can be queried
# Returns a list of pathways
#==============================================================================================================
get_pathways <- function() {
  # Create a dataframe with the available databases
  enrichment_databases <- listGeneSet()
  
  # Remove unneeded databases
  enrichment_databases <- enrichment_databases %>%
    filter(!grepl("cancer|CPTAC|TCGA|community-contributed", name)) %>%
    filter(!xor(grepl("^geneontology", name), grepl("noRedundant$", name)))
  
  # Create a column with the database type and label for the select input
  enrichment_databases <- enrichment_databases %>%
    mutate(type = sapply(name, function(x) unlist(str_split(x, "_"))[1]),
           label = sapply(name, function(x) paste(unlist(str_split(x, "_"))[-1], collapse = "_"))) %>%
    dplyr::select(c(name, type, label))
  
  # Create a named list that will be used to for the grouped select input
  database_types <- unique(enrichment_databases$type)
  enrichment_database_list <- lapply(database_types, function(x) {
    databases <- filter(enrichment_databases, type == x)
    type_list <- split(databases$name, databases$label)
    
    return(type_list)
  })
  names(enrichment_database_list) <- database_types

  return(enrichment_database_list)
}
#==============================================================================================================


#==============================================================================================================
# Query the database table for all database info
# Returns a dataframe with database names and descriptions
#==============================================================================================================
get_db_info <- function(db) {
  # Import the filtered geneset description
  database_links <- dbGetQuery(db,
                               stringr::str_interp(
                                 paste("SELECT database_name, database_full_name, type, summary, link",
                                       "FROM database")))

  return(database_links)
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


