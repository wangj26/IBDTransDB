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
# Function used to query the dataset table for the description of the given dataset
#==============================================================================================================
get_dataset_description <- function(dataset_acc, db) {
  
  # Format the dataset ID so that it can be used for the query
  dataset_acc <- paste0("'", dataset_acc, "'")
  
  # Pull the datasets
  dataset <- dbGetQuery(db,
                        stringr::str_interp(
                          paste("SELECT id, dataset_acc, title,",
                                "disease, organism, experiment_type,",
                                "source, cell_type, sample_number,",
                                "treatment, timepoint, dose,",
                                "platform, normalization_method,",
                                "summary, design,", 
                                "contact_name, contact_email",
                                "FROM dataset",
                                "WHERE dataset_acc = ${dataset_acc}")))
  
  return(dataset)
}
#==============================================================================================================


#==============================================================================================================
# Function used to query the dataset_data table for the expression data of the given dataset
#==============================================================================================================
# get_full_exp_data <- function(current_dataset_acc, dataset_acc_identifier, gene_map, sample_ann, db) {
#   
#   # Pull the expression data
#   exp_data <- dbGetQuery(db,
#                          stringr::str_interp(
#                            paste("SELECT gene, sample_acc, expr",
#                                  "FROM dataset_data",
#                                  "WHERE dataset_acc = ${dataset_acc_identifier}")))
#   # Convert the dataset, gene, and sample identifiers to their respective names
#   sample_ann <- select(sample_ann, c(id, sample_id))
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
# Function used to query the dataset_data table for the expression data of the given dataset and selected genes
#==============================================================================================================
get_exp_data <- function(genes, gene_map, dataset_acc_identifier, sample_ann, db) {
  
  # Format the genes so that they can be used for the query
  genes <- gene_map$id[match(genes, gene_map$gene)]
  genes <- paste0("(", paste(genes, collapse = ", "), ")")
  
  # Pull the expression data
  exp_data <- dbGetQuery(db,
                         stringr::str_interp(
                           paste("SELECT gene, sample_acc, expr",
                                 "FROM dataset_data",
                                 "WHERE (dataset_acc = ${dataset_acc_identifier}",
                                 "AND gene IN ${genes})")))
  # Convert the dataset, gene, and sample identifiers to their respective names
  sample_ann <- select(sample_ann, c(id, sample_id))
  exp_data <- exp_data %>%
    merge(gene_map, by.x = "gene", by.y = "id", all.x = TRUE) %>%
    select(-gene) %>%
    rename("gene" = "gene.y") %>%
    merge(sample_ann, by.x = "sample_acc", by.y = "id", all.x = TRUE) %>%
    select(-sample_acc) %>%
    rename("sample_acc" = "sample_id")
  
  return(exp_data)
}
#==============================================================================================================


#==============================================================================================================
# Function used to query the sample_ann table for the sample annotations of the given dataset
#==============================================================================================================
get_sample_ann <- function(dataset_acc, db) {
  
  # Format the dataset ID so that it can be used for the query
  dataset_acc <- paste0("'", dataset_acc, "'")
  
  # Pull the sample annotations
  sample_ann <- dbGetQuery(db,
                           stringr::str_interp(
                             paste("SELECT id, sample_id, sample_ann_type, sample_ann_value",
                                   "FROM sample_ann",
                                   "WHERE dataset_acc = ${dataset_acc}")))
  
  return(sample_ann)
}
#==============================================================================================================


#==============================================================================================================
# Function used to query the comparison table for the comparisons for the given dataset
#==============================================================================================================
get_comparisons <- function(dataset_acc, db) {
  
  # Format the dataset ID so that it can be used for the query
  dataset_acc <- paste0("'", dataset_acc, "'")
  
  # Pull the comparisons
  comparisons <- dbGetQuery(db,
                            stringr::str_interp(
                              paste("SELECT id, case_ann, case_sample, control_ann, control_sample",
                                    "FROM comparison",
                                    "WHERE dataset_acc = ${dataset_acc}")))
  
  # Add a column that indicates the comparisons
  comparisons <- mutate(comparisons, comparison = paste(case_ann, "vs", control_ann))
  
  return(comparisons)
}
#==============================================================================================================


#==============================================================================================================
# Function used to query the comparison table for the comparisons for the given dataset (single cell)
#==============================================================================================================
get_sc_comparisons <- function(dataset_acc, db) {
  
  # Format the dataset ID so that it can be used for the query
  dataset_acc <- paste0("'", dataset_acc, "'")
  
  # Pull the comparisons
  comparisons <- dbGetQuery(db,
                            stringr::str_interp(
                              paste("SELECT id, case_ann, case_sample, control_ann, control_sample",
                                    "FROM comparison",
                                    "WHERE dataset_acc = ${dataset_acc}")))
  
  # Add a column that indicates the comparisons
  comparisons <- mutate(comparisons, comparison = paste(id, case_ann, "vs", control_ann))
  
  return(comparisons)
}
#==============================================================================================================


#==============================================================================================================
# Function used to query the comparison_data table for comparison metrics for the given comparison
#==============================================================================================================
get_comparison_data <- function(comparison_ids, gene_map, db) {
  
  # Format the comparison IDs so that it can be used for the query
  if (length(comparison_ids) == 1) {
    comparison_ids <- paste0("('", comparison_ids, "')")
  } else {
    comparison_ids <- paste0("('", paste(comparison_ids, collapse = "', '"), "')")
  }

  # Pull the comparison data
  comparison_data <- dbGetQuery(db,
                                stringr::str_interp(
                                  paste("SELECT comparison_id, gene, log_fc, p_value, p_value_adj",
                                        "FROM comparison_data",
                                        "WHERE comparison_id IN ${comparison_ids}")))
  
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
