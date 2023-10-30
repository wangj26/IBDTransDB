#==============================================================================================================
# Functions used to calculate one-tailed p-values
#==============================================================================================================
one_tailed_upreg <- function(p_value, log_fc) {
  if (log_fc > 0) {
    p_value <- p_value/2
    return(p_value)
  } else {
    p_value <- 1 - p_value/2
    return(p_value)
  }
}

one_tailed_downreg <- function(p_value, log_fc) {
  if (log_fc > 0) {
    p_value <- 1 - p_value/2
    return(p_value)
  } else {
    p_value <- p_value/2
    return(p_value)
  }
}
#==============================================================================================================


#==============================================================================================================
# Function used to rank targets based on comparison evidence
# Returns a dataframe of ranked targets
#==============================================================================================================
rank_targets <- function(genes, comparisons, comparison_data) {
  
  # Compute metrics for each gene - use parallel backend (8 cores)
  gene_scores <- mclapply(genes, function(current_gene) {
    
    ### Compute the meta p-value for the given gene using Fisher's method ###
    # Get the associated comparison data for the given gene
    gene_data <- comparison_data %>%
      filter(gene == current_gene) %>%
      select(comparison_id, p_value, log_fc)
    # Calculate the meta p-value
    if (nrow(gene_data) == 0) {
      meta_p <- "NA"
      regulation_dir <- "NA"
      num_comparisons <- "NA"
      num_datasets <- "NA"
      mean_fc <- "NA"
    } else {
      # Must convert two-tailed p-values to one-tailed p-values
      
      # Calculate the meta p-value using one-tailed p-values from alternative hypothesis: UPREGULATED
      p_values <- mapply(one_tailed_upreg, gene_data$p_value, gene_data$log_fc)
      meta_p_upreg <- fisher(p_values)$p
      # Calculate the meta p-value using one-tailed p-values from alternative hypothesis: DOWNREGULATED
      p_values <- mapply(one_tailed_downreg, gene_data$p_value, gene_data$log_fc)
      meta_p_downreg <- fisher(p_values)$p
      
      # Select the meta p-value that is smaller
      if (meta_p_upreg < meta_p_downreg) {
        meta_p <- meta_p_upreg
        regulation_dir <- "up"
      } else {
        meta_p <- meta_p_downreg
        regulation_dir <- "down"
      }
      
      # Determine the number of comparisons and datasets associated with the given gene
      num_comparisons <- comparisons %>%
        filter(id %in% gene_data$comparison_id) %>%
        nrow()
      num_datasets <- comparisons %>%
        filter(id %in% gene_data$comparison_id) %>%
        pull(dataset_acc) %>%
        unique() %>%
        length()
      
      mean_fc <- mean(gene_data$log_fc)
    }
    ##########
    
    return(c(current_gene, meta_p, mean_fc, regulation_dir, num_comparisons, num_datasets))
  }, mc.cores = 8)
  integrate_df <- as.data.frame(do.call(rbind, gene_scores))
  colnames(integrate_df) <- c("target", "meta_p_value", "mean_logFC","direction", "num_comparisons", "num_datasets")
  
  # Calculate adjusted meta p-values and format columns of the dataframe 
  # Filter out genes that don't have any associated comparisons - CHECK THIS
  integrate_df <- integrate_df %>%
    filter(meta_p_value != "NA") %>%
    mutate(meta_p_value = as.numeric(meta_p_value), mean_logFC=as.numeric(mean_logFC)) %>%
    mutate(adjusted_meta_p_value = p.adjust(meta_p_value, method = "BH")) %>%
    arrange(meta_p_value) %>%
    mutate(meta_p_value = signif(meta_p_value, 5),
           mean_logFC = signif(mean_logFC, 5),
           adjusted_meta_p_value = signif(adjusted_meta_p_value, 5)) %>%
    relocate(target, meta_p_value, adjusted_meta_p_value, mean_logFC,direction, num_comparisons, num_datasets) %>%
    rename("meta p-value" = "meta_p_value", 
           "adjusted meta p-value" = "adjusted_meta_p_value",
           "mean logFC" = "mean_logFC",
           "# of comparisons" = "num_comparisons",
           "# of datasets" = "num_datasets")
  
  return(integrate_df)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate extra info for targets
# Returns a list of DT objects
#==============================================================================================================
generate_target_extra_info <- function(genes, comparisons, comparison_data) {
  # Compute metrics for each gene - use parallel backend (8 cores)
  extra_info <- lapply(genes, function(current_gene) {
    # Get the associated comparison data for the given gene and merge with the comparison description
    gene_data <- comparison_data %>%
      filter(gene == current_gene) %>%
      select(comparison_id, p_value, p_value_adj, log_fc) %>%
      left_join(comparisons, by = c("comparison_id" = "id")) %>%
      select(dataset_acc, case_ann, control_ann, log_fc, p_value, p_value_adj) %>%
      mutate(p_value = signif(p_value, 5),
             log_fc = signif(log_fc, 5),
             p_value_adj = signif(p_value_adj, 5)) %>%
      rename("dataset ID" = "dataset_acc", "case group" = "case_ann", "control group" = "control_ann",
             "logFC" = "log_fc", "p-value" = "p_value", "adj. p-value" = "p_value_adj") %>%
      datatable(options = list(dom = "p", pageLength = 10))
    
    return(gene_data)
  })
  
  return(extra_info)
}
#==============================================================================================================


#==============================================================================================================
# Function used to rank a signature based on comparison evidence
# Returns a dataframe of metrics for the given signature
#==============================================================================================================
rank_signature <- function(genes, selected_comparison_ids, comparisons, db) {
  
  comparisons <- filter(comparisons, id %in% selected_comparison_ids)
  
  # Pull the expression data associated with the given datasets and specified signature
  exp_data <- get_signature_data(genes, comparisons$dataset_acc, db)

  # Perform a Wilcoxon rank sum test for each comparison and store the p-values
  p_values <- lapply(selected_comparison_ids, function(current_id) {
    # Pull the case and control samples for the given comparison
    case_samples <- comparisons %>%
      filter(id == current_id) %>%
      pull(case_sample) %>%
      strsplit(";") %>%
      unlist()
    control_samples <- comparisons %>%
      filter(id == current_id) %>%
      pull(control_sample) %>%
      strsplit(";") %>%
      unlist()
    
    # Determine the number of genes included in the comparison
    num_genes <- exp_data %>%
      filter(sample_acc %in% c(case_samples, control_samples)) %>%
      pull(gene) %>%
      unique() %>%
      length()
    
    if (num_genes == 0) {
      return(NULL)
    }
    
    # Run a Wilcoxon rank sum test between the two groups
    case_values <- exp_data %>%
      filter(sample_acc %in% case_samples) %>%
      select(-gene) %>%
      group_by(dataset_acc, sample_acc) %>%
      summarize(expr = mean(expr)) %>%
      pull(expr)
    control_values <- exp_data %>%
      filter(sample_acc %in% control_samples) %>%
      select(-gene) %>%
      group_by(dataset_acc, sample_acc) %>%
      summarize(expr = mean(expr)) %>%
      pull(expr)
    # Alternative is greater
    test_results <- wilcox.test(case_values, control_values, alternative = "greater")
    p_val_greater <- signif(test_results$p.value, 6)
    # Alternative is less
    test_results <- wilcox.test(case_values, control_values, alternative = "less")
    p_val_less <- signif(test_results$p.value, 6)
    
    return(c(p_val_greater, p_val_less))
  })
  p_val_df <- as.data.frame(do.call(rbind, p_values))
  colnames(p_val_df) <- c("p_val_upreg", "p_val_downreg")
  
  # Perform Fisher's method on the one-sided p-values
  meta_p_upreg <- fisher(p_val_df$p_val_upreg)$p
  meta_p_downreg <- fisher(p_val_df$p_val_downreg)$p
  
  # Select the meta p-value that is smaller
  if (meta_p_upreg < meta_p_downreg) {
    meta_p <- meta_p_upreg
    regulation_dir <- "upregulated"
  } else {
    meta_p <- meta_p_downreg
    regulation_dir <- "downregulated"
  }
  
  # Create a dataframe used to store the target ranks and associated metrics
  integrate_df <- data.frame("meta_p_value" = meta_p, "direction" = regulation_dir)
  
  # Modify columns of the table
  integrate_df <- integrate_df %>%
    rename("meta p-value" = "meta_p_value")
  
  return(integrate_df)
}
#==============================================================================================================


#==============================================================================================================
# Function used to rank pathways based on comparison evidence
# Returns a dataframe of ranked pathways
#==============================================================================================================
rank_pathways <- function(comparisons,pathway_data) {
  
  pathways <- unique(pathway_data$description)
  pathway_metrics <- mclapply(pathways, function(x) {
    current_pathway_data <- filter(pathway_data, description == x)
    
    # Calculate the meta p-value using one-tailed p-values from alternative hypothesis: UPREGULATED
    p_values <- mapply(one_tailed_upreg, current_pathway_data$p_value, current_pathway_data$es)
    meta_p_upreg <- fisher(p_values)$p
    # Calculate the meta p-value using one-tailed p-values from alternative hypothesis: UPREGULATED
    p_values <- mapply(one_tailed_downreg, current_pathway_data$p_value, current_pathway_data$es)
    meta_p_downreg <- fisher(p_values)$p
    
    # Select the meta p-value that is smaller
    if (meta_p_upreg < meta_p_downreg) {
      meta_p <- meta_p_upreg
      regulation_dir <- "upregulated"
    } else {
      meta_p <- meta_p_downreg
      regulation_dir <- "downregulated"
    }
    
    # Calculate the mean ES and mean NES
    mean_es <- mean(current_pathway_data$es)
    mean_nes <- mean(current_pathway_data$nes)
    
    # Determine the number of comparisons and datasets associated with the given gene
    num_comparisons <- comparisons %>%
      filter(id %in% pathway_data$comparison_id) %>%
      nrow()
    num_datasets <- comparisons %>%
      filter(id %in% pathway_data$comparison_id) %>%
      pull(dataset_acc) %>%
      unique() %>%
      length()
    
    return(c(x, meta_p, mean_es, mean_nes, regulation_dir, num_comparisons, num_datasets))
  }, mc.cores = 8)
  
  # Create a dataframe used to store the pathway integration metrics
  integrate_df <- as.data.frame(do.call(rbind,pathway_metrics))
  colnames(integrate_df) <- c("pathway", "meta_p_value", "mean_es", "mean_nes", "direction", "num_comparisons", "num_datasets")


  
  # Modify columns of the table
  integrate_df <- integrate_df %>%
    filter(meta_p_value != "NA") %>%
    mutate(meta_p_value = as.numeric(meta_p_value), mean_es = as.numeric(mean_es), mean_nes = as.numeric(mean_nes)) %>%
    arrange(meta_p_value) %>%
    mutate(adjusted_meta_p_value = p.adjust(meta_p_value, method = "BH")) %>%
    mutate(meta_p_value = signif(meta_p_value, 5),
           adjusted_meta_p_value = signif(adjusted_meta_p_value, 5),
           mean_es = signif(mean_es,5),
           mean_nes = signif(mean_nes,5)) %>%
    relocate(pathway, meta_p_value, adjusted_meta_p_value, mean_es, mean_nes, direction, num_comparisons, num_datasets) %>%
    rename("meta p-value" = "meta_p_value", 
           "adjusted meta p-value" = "adjusted_meta_p_value",
           "# of comparisons" = "num_comparisons",
           "# of datasets" = "num_datasets",
           "mean es" = "mean_es",
           "mean nes" = "mean_nes")
  
  return(integrate_df)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate extra info for pathways
# Returns a list of DT objects
#==============================================================================================================
generate_pathway_extra_info <- function(rank_df,comparisons, pathway_data) {
  # Compute metrics for each gene - use parallel backend (8 cores)
  pathways <- unique(rank_df$pathway)
  extra_info <- lapply(pathways, function(current_pathway) {
    # Get the associated comparison data for the given gene and merge with the comparison description
    pathway_data <- pathway_data %>%
      filter(description == current_pathway) %>%
      left_join(comparisons, by = c("comparison_id" = "id")) %>%
      select(dataset_acc, case_ann, control_ann, es, nes, p_value, p_value_adj) %>%
      mutate(es = signif(es, 5),
             nes = signif(nes, 5),
             p_value = signif(p_value,5),
             p_value_adj = signif(p_value_adj,5)) %>%
      rename("dataset ID" = "dataset_acc", "case group" = "case_ann", "control group" = "control_ann",
             "ES" = "es", "NES" = "nes", "p-value" = "p_value", "adj. p-value" = "p_value_adj") %>%
      datatable(options = list(dom = "p", pageLength = 10))
    
    return(pathway_data)
  })
  
  return(extra_info)
}
#==============================================================================================================


#==============================================================================================================
# Function used to create a DT table from a ranked targets dataframe
# Returns a reactable table
#==============================================================================================================
style_rank_df <- function(integrate_df, extra_info) {
  # integrate_df <- datatable(integrate_df,
  #                           options = list(pageLength = 15, 
  #                                          lengthChange = FALSE,
  #                                          dom = '<"top">ft<"bottom">ip'
  #                           ),
  #                           rownames = FALSE)
  
  # Generate the reactable table
  integrate_df <- reactable(
    integrate_df,
    defaultColDef = colDef(
      minWidth = 75,
      headerStyle = list(background = "#f7f7f8")
    ),
    bordered = TRUE,
    striped = TRUE,
    highlight = TRUE,
    height = 500,
    pagination = FALSE,
    sortable = TRUE,
    details = function(index) {
      htmltools::div(style = "padding: 1rem",
                     htmltools::p("Data included in integration analysis: ", style = "font-weight: bold;"),
                     htmltools::br(),
                     extra_info[[index]])
    }
  )
  
  return(integrate_df)
}
#==============================================================================================================