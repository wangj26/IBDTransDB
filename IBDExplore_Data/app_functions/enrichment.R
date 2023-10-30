#==============================================================================================================
# Function used to run FGSEA
# Returns a DT table
#==============================================================================================================
# run_gsea <- function(selected_comparison, selected_collection, comparisons, comparison_data, active_enrichment_data) {
#   
#   # Use the selected comparison ID to filter the comparison data
#   selected_comparison_id <- comparisons %>%
#     filter(comparison == selected_comparison) %>%
#     pull(id)
#   comparison_data <- comparison_data %>%
#     filter(comparison_id == selected_comparison_id) %>%
#     select(c(gene, log_fc))
#   
#   # Create a named list of genes where the ranking metric is the logFC
#   ranked_genes <- comparison_data$log_fc
#   names(ranked_genes) <- comparison_data$gene
#   
#   # Pull the MSigDB gene set collection that is specified (e.g. hallmark geneset)
#   gene_set_collection <- msigdbr(species = "human", category = selected_collection)
#   
#   # Format the gene sets so that they are compatible with fgsea
#   gene_set_collection <- split(x = gene_set_collection$gene_symbol,
#                                f = gene_set_collection$gs_name)
#   
#   # Perform GSEA
#   gsea_results <- fgsea(pathways = gene_set_collection,
#                         stats = ranked_genes,
#                         minSize = 15,
#                         maxSize = 500)
#   
#   ### Update REACTIVE VALUES ###
#   active_enrichment_data$gene_set_colelction <- gene_set_collection
#   active_enrichment_data$ranked_genes <- ranked_genes
#   active_enrichment_data$gsea_results <- gsea_results
# 
#   # Keep the top 20 most significant gene sets
#   gsea_results <- gsea_results %>%
#     top_n(20, wt = -padj) %>%
#     arrange(desc(abs(NES)))
#   
#   ### Update REACTIVE VALUES ###
#   active_enrichment_data$top_gsea_results <- gsea_results
#   
#   # Generate the DT table
#   gsea_results <- gsea_results %>%
#     select(-c(nMoreExtreme, size, leadingEdge)) %>%
#     rename("p-value" = "pval", "adj p-value" = "padj") %>%
#     datatable(gsea_results,
#               extensions = "Scroller",
#               selection = "single",
#               options = list(autoWidth = TRUE,
#                              scrollX = TRUE,
#                              scrollY = 350,
#                              scroller = TRUE)) %>%
#     formatStyle(columns = 1, width = "175px") %>%
#     formatStyle(columns = c(2, 3, 4, 5), width = "100px") %>%
#     formatStyle(columns = c(1, 2, 3, 4, 5), fontSize = "90%")
#   
#   return(gsea_results)
# }
#==============================================================================================================


#==============================================================================================================
# Function used to create an enrichment score plot
#==============================================================================================================
# generate_enrichment_score_plot <- function(selected_row, active_enrichment_data) {
#   
#   # Determine the associated gene set given the selected row
#   selected_gene_set <- pull(active_enrichment_data[["top_gsea_results"]], pathway)
#   selected_gene_set <- selected_gene_set[selected_row]
#   
#   # Generate the enrichment score plot for the given gene set
#   enrichment_score_plot <- plotEnrichment(active_enrichment_data[["gene_set_collection"]][[selected_gene_set]],
#                                           active_enrichment_data[["ranked_genes"]])
#   
#   return(enrichment_score_plot)
# }
#==============================================================================================================


#==============================================================================================================
# Function used to list databases users can select for the enrichment analysis
# Returns a dataframe with the type and name of the available databases
#==============================================================================================================
list_enrichment_databases <- function() {
  # Create a dataframe with the available databases
  enrichment_databases <- listGeneSet()
  
  # Remove unneeded databases
  enrichment_databases <- enrichment_databases %>%
    filter(!grepl("Kinase_target|cancer|CPTAC|TCGA|community-contributed", name)) %>%
    filter(!(idType == "phosphositeSeq")) %>%
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
# Function used to generate a map of organism names (names in DB and names used by Webgestalt)
# Returns a named list
#==============================================================================================================
generate_organism_map <- function() {
  # Create the organism map
  organism_map <- list("Human" = "hsapiens",
                       "Mouse" = "mmusculus",
                       "Rat" = "rnorvegicus")
  
  return(organism_map)
}
#==============================================================================================================


#==============================================================================================================
# Function used to perform ORA using fora
# Returns a dataframe with pathways
#==============================================================================================================
perform_ora <- function(selected_comparison, organism, database, dataset_acc, pval_option, p_value_threshold, logfc_option, log_fc_threshold, comparisons, comparison_data) {
  
  # Filter the comparison data for the selected comparison and pull genes that meet the specified filters
  selected_comparison_id <- comparisons %>%
    filter(comparison == selected_comparison) %>%
    pull(id)
  if (logfc_option == "up") {
    genes <- comparison_data %>% 
      filter(comparison_id == selected_comparison_id,
             log_fc >= log_fc_threshold,
             comparison_data[[pval_option]] <= p_value_threshold) %>%
      pull(gene)
  } else {
    genes <- comparison_data %>% 
      filter(comparison_id == selected_comparison_id,
             log_fc <= -log_fc_threshold,
             comparison_data[[pval_option]] <= p_value_threshold) %>%
      pull(gene)
  }
  universe <- comparison_data %>% 
    filter(comparison_id == selected_comparison_id) %>%
    pull(gene)
  
  if (length(genes) == 0) {
    return(NULL)
  }
  
  # Define the output directory for the results
  project_name <- paste(dataset_acc, selected_comparison_id, database, sep = "_")
  project_name <- gsub("-", "_", project_name)
  
  # Perform ORA
  enrichment_results <- WebGestaltR(enrichMethod = "ORA", 
                                    organism = organism, 
                                    enrichDatabase = database,
                                    interestGene = genes,
                                    interestGeneType = "genesymbol",
                                    referenceGene = universe,
                                    referenceGeneType = "genesymbol",
                                    collapseMethod = "mean",
                                    minNum = 5,
                                    maxNum = 2000,
                                    fdrMethod = "BH",
                                    sigMethod = "fdr",
                                    fdrThr = 0.1,
                                    isOutput = TRUE,
                                    outputDirectory = "enrichment_results",
                                    projectName = project_name,hostName="https://www.webgestalt.org/")
  if(is.null(enrichment_results)){
    return(NULL)
  }else{
    return(project_name)
  }
}
#==============================================================================================================


#==============================================================================================================
# Function used to perform GSEA using Webgestalt
# Returns the project name and downloads an html report
#==============================================================================================================
perform_gsea <- function(selected_comparison, organism, database, dataset_acc, comparisons, comparison_data) {
  
  # Pull the comparison ID and pull the ranked gene list associated with the given comparison
  selected_comparison_id <- comparisons %>%
    filter(comparison == selected_comparison) %>%
    pull(id)
  ranked_genes <- comparison_data %>%
    filter(comparison_id == selected_comparison_id) %>%
    select(c(gene, log_fc))
  
  # Define the output directory for the results
  project_name <- paste(dataset_acc, selected_comparison_id, database, sep = "_")
  project_name <- gsub("-", "_", project_name)
  
  # Run Webgestalt
  enrichment_results <- WebGestaltR(enrichMethod = "GSEA", 
                                    organism = organism, 
                                    enrichDatabase = database,
                                    interestGene = ranked_genes,
                                    interestGeneType = "genesymbol",
                                    collapseMethod = "mean",
                                    minNum = 5,
                                    maxNum = 2000,
                                    isOutput = TRUE,
                                    outputDirectory = "enrichment_results",
                                    projectName = project_name,hostName="https://www.webgestalt.org/")
  
  return(project_name)
}
#==============================================================================================================


#==============================================================================================================
# Function used to modify the enrichment report produced by Webgestalt
# Writes the new report to the appropriate location
#==============================================================================================================
modify_enrichment_report <- function(project_name) {
  # Import the enrichment report and cast the rownames as numeric
  enrichment_report <- read.table(paste0("./enrichment_results/Project_", project_name, "/Report_", project_name, ".html"), 
                                  sep = "\n", quote = "")
  rownames(enrichment_report) <- as.numeric(rownames(enrichment_report))
  
  # Modify rows of the enrichment report
  row_index <- which(enrichment_report$V1 == "<hr><main>")
  enrichment_report$V1[row_index] <- gsub("<hr>", "", enrichment_report$V1[row_index])
  row_index <- which(grepl('<a href="$', enrichment_report$V1)) 
  enrichment_report$V1[row_index] <- '<a href="#" class="card-header-icon">'
  
  ### Remove the header and footer of the report ###
  title_indices <- which(startsWith(enrichment_report$V1, "<title>"))
  title_indices <- c(title_indices, title_indices + 1)
  #==========
  header_start <- which(startsWith(enrichment_report$V1, "<header>"))
  header_end <- which(startsWith(enrichment_report$V1, "</header>"))
  header_indices <- header_start:header_end
  #==========
  footer_start <- which(startsWith(enrichment_report$V1, "<footer"))
  footer_end <- which(startsWith(enrichment_report$V1, "</footer>"))
  footer_indices <- footer_start:(footer_end - 1)
  #==========
  enrichment_report$V1[footer_end] <- gsub("</footer>", "", enrichment_report$V1[footer_end])
  #==========
  rows_to_remove <- c(title_indices, header_indices, footer_indices)
  enrichment_report <- filter(enrichment_report, !(rownames(enrichment_report) %in% rows_to_remove))
  ##########
  
  # Escape double quotes and single quotes
  enrichment_report$V1 <- gsub('"', '\"', enrichment_report$V1)
  enrichment_report$V1 <- gsub("'", "\'", enrichment_report$V1)
  
  write.table(enrichment_report, file = paste0("./enrichment_results/Project_", project_name, "/Report_", project_name, ".html"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
#==============================================================================================================


#==============================================================================================================
# Function used to upload the enrichment results to the DB
# Writes results to the enrichment_run, enrichment_results, and enrichment_leading_edge tables
#==============================================================================================================
upload_enrichment_results <- function(selected_comparison, organism, selected_database, dataset_acc, comparisons, comparison_data, db) {
  
  # Pull the comparison ID and pull the ranked gene list associated with the given comparison
  selected_comparison_id <- comparisons %>%
    filter(comparison == selected_comparison) %>%
    pull(id)
  ranked_genes <- comparison_data %>%
    filter(comparison_id == selected_comparison_id) %>%
    mutate(rank = -1*log10(p_value_adj) * sign(log_fc)) %>%
    select(c(gene, rank))
  
  # Skip the upload if the enrichment analysis is already run for the given dataset, comparison, and enrichment database
  enrichment_results <- dbGetQuery(db,
                                   stringr::str_interp(
                                     paste("SELECT *",
                                           "FROM enrichment_run",
                                           "WHERE (dataset_acc = ${dataset_acc}",
                                           "AND comparison_id = ${selected_comparison_id}",
                                           "AND database = ${selected_database})")))
  
  if (nrow(enrichment_results) > 0) {
    # Run Webgestalt
    enrichment_results <- WebGestaltR(enrichMethod = "GSEA", 
                                      organism = organism, 
                                      enrichDatabase = selected_database,
                                      interestGene = ranked_genes,
                                      interestGeneType = "genesymbol",
                                      collapseMethod = "mean",
                                      minNum = 5,
                                      maxNum = 2000,
                                      isOutput = FALSE,
                                      fdrThr = 1)
    
    ### ENRICHMENT RUN UPLOAD ###
    enrichment_run <- data.frame("dataset_acc" = dataset_acc,
                                 "comparison_id" = selected_comparison_id,
                                 "database" = selected_database)
    dbAppendTable(db, "enrichment_run", enrichment_run)
    ##########
    
    ### ENRICHMENT RESULTS UPLOAD ###
    # Store the leading edge info
    leading_edge_data <- select(enrichment_results, c(geneSet, userId))
    # Select columns and modify column names
    enrichment_results <- enrichment_results %>%
      select(-c(link, plotPath, leadingEdgeId, userId)) %>%
      rename(geneset = geneSet,
             leading_edge_number = leadingEdgeNum,
             es = enrichmentScore,
             nes = normalizedEnrichmentScore,
             p_value = pValue,
             p_value_adj = FDR)
    # Pull the primary key ID from the enrichment run table
    dataset_id <- paste0("'", dataset_acc, "'")
    comparison_id <- paste0("'", selected_comparison_id, "'")
    database <- paste0("'", selected_database, "'")
    run_id <- dbGetQuery(db,
                         stringr::str_interp(
                           paste("SELECT id",
                                 "FROM enrichment_run",
                                 "WHERE (dataset_acc = ${dataset_id}",
                                 "AND comparison_id = ${comparison_id}",
                                 "AND database = ${database})")))
    run_id <- pull(run_id, id)
    enrichment_results$enrichment_run_id <- run_id
    # Upload the enrichment results
    dbAppendTable(db, "enrichment_results", enrichment_results)
    ##########
    
    ### ENRICHMENT LEADING EDGE UPLOAD ###
    # Pull the IDs associated with each gene set in the enrichment results table
    run_id <- paste0("'", run_id, "'")
    results_ids <- dbGetQuery(db,
                              stringr::str_interp(
                                paste("SELECT id, geneset",
                                      "FROM enrichment_results",
                                      "WHERE enrichment_run_id = ${run_id}")))
    # Pull the associated comparison data
    comparison_data <- filter(comparison_data, comparison_id == selected_comparison_id)
    
    # Create a list of dataframes containing the leading edge genes for each significant gene set
    leading_edge_list <- lapply(leading_edge_data$geneSet, function(x) {
      current_id <- results_ids %>%
        filter(geneset == x) %>%
        pull(id)
      
      genes <- leading_edge_data %>%
        filter(geneSet == x) %>% 
        pull(userId)
      genes <- unlist(strsplit(genes, split = ";"))
      
      leading_edge <- comparison_data %>%
        select(gene, log_fc) %>%
        filter(gene %in% genes) %>%
        arrange(desc(log_fc)) %>%
        mutate(enrichment_results_id = current_id) %>%
        rename(score = log_fc)
    })
    # Combine the dataframes 
    leading_edge_results <- do.call("rbind", leading_edge_list)
    
    # Upload the leading edge results
    dbAppendTable(db, "enrichment_leading_edge", leading_edge_results)
    ##########
  }
}
#==============================================================================================================
