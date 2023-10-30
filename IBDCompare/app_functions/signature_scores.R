#==============================================================================================================
# Function used to perform signature-based comparisons
# Returns a dataframe containing p-values derived from Wilcoxon rank sum tests
#==============================================================================================================
perform_signature_comparisons <- function(genes, comparisons, db) {

  # Pull the expression data associated with the selected datasets and genes
  dataset_acc_list <- unique(comparisons$dataset_acc)
  exp_data <- get_signature_data(genes, dataset_acc_list, db)
  
  # Perform signature-based comparisons
  signature_data <- lapply(1:nrow(comparisons), function(i) {
    current_comparison_id <- comparisons$id[i]
    
    # Calculate signature scores for the case and control groups
    case_data <- exp_data %>%
      filter(sample_acc %in% unlist(strsplit(comparisons$case_sample[i], ";"))) %>%
      select(sample_acc, expr) %>%
      group_by(sample_acc) %>%
      summarize(expr = mean(expr)) %>%
      pull(expr)
    control_data <- exp_data %>%
      filter(sample_acc %in% unlist(strsplit(comparisons$control_sample[i], ";"))) %>%
      select(sample_acc, expr) %>%
      group_by(sample_acc) %>%
      summarize(expr = mean(expr)) %>%
      pull(expr)
    
    if (length(case_data) == 0 | length(control_data) == 0) {
      return(NULL)
    }
    
    # Determine the log2 fold change between the case and control groups
    log_fc <- log2(mean(case_data)/mean(control_data))
    # Perform a Wilcoxon rank sum test between the case and control signature scores
    p_val <- wilcox.test(case_data, control_data, alternative = "two.sided")$p.value
    # Determine the number of genes pertaining to the signature that are found in the given dataset
    num_genes <- exp_data %>%
      filter(dataset_acc == comparisons$dataset_acc[i]) %>%
      pull(gene) %>%
      unique() %>%
      length()
    
    return(c(current_comparison_id, log_fc, p_val, num_genes))
  })
  signature_data <- data.frame(do.call(rbind, signature_data))
  colnames(signature_data) <- c("id", "log_fc", "p_value", "num_genes")
  signature_data <- signature_data %>%
    left_join(comparisons, by = c("id" = "id")) %>%
    select(id, dataset_acc, log_fc, p_value, num_genes) %>%
    rename("comparison_id" = "id")
  
  return(signature_data)
}
#==============================================================================================================
