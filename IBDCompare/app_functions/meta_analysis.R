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
# Function used to perform meta analysis for selected comparisons
# Returns a reactable object
#==============================================================================================================
perform_meta_analysis <- function(datasets, comparisons, comparison_data, checked_comparisons) {
  # Filter for selected comparisons and compute one-tailed p-values (alternative hypotheses UPREGULATED and DOWNREGULATED)
  checked_comparisons <- lapply(checked_comparisons, function(x) unlist(str_split(x, "_"))[2])
  comparison_data <- comparison_data %>%
    filter(comparison_id %in% checked_comparisons) %>%
    mutate(up_p_value = mapply(one_tailed_upreg, p_value, log_fc),  
           down_p_value = mapply(one_tailed_downreg, p_value, log_fc))
  
  # Perform the meta analysis
  genes <- unique(comparison_data$gene)
  meta_analysis_results <- lapply(genes, function(current_gene) {
    current_comparison_data <- filter(comparison_data, gene == current_gene)
    # Perform meta analysis using Fisher's method
    meta_p_upreg <- fisher(current_comparison_data$up_p_value)$p
    meta_p_downreg <- fisher(current_comparison_data$down_p_value)$p
    # Select the meta p-value that is smaller
    if (meta_p_upreg < meta_p_downreg) {
      meta_p <- meta_p_upreg
      regulation_dir <- "up"
    } else {
      meta_p <- meta_p_downreg
      regulation_dir <- "down"
    }
    return(c(current_gene, meta_p, regulation_dir))
  })
  meta_analysis_results <- data.frame(do.call(rbind, meta_analysis_results))
  colnames(meta_analysis_results) <- c("gene", "meta_p", "regulation_dir")
  
  # Compute corrected meta p-values and format the columns
  meta_analysis_results <- meta_analysis_results %>%
    mutate(meta_p = as.numeric(meta_p)) %>%
    mutate(adj_meta_p = p.adjust(meta_p, method = "BH")) %>%
    mutate(meta_p = signif(meta_p, 5),
           adj_meta_p = signif(adj_meta_p, 5)) %>%
    relocate(gene, meta_p, adj_meta_p, regulation_dir) %>%
    rename("meta p" = "meta_p", "adj. meta p" = "adj_meta_p", "direction" = "regulation_dir")
  
  # Generate the extra info that will be displayed in expanded sections of the table
  extra_info <- lapply(genes, function(current_gene) {
    data <- comparison_data %>%
      filter(gene == current_gene) %>%
      left_join(comparisons, by = c("comparison_id" = "id")) %>%
      left_join(datasets, by = c("dataset_acc" = "dataset_acc")) %>%
      select(dataset_acc, case_ann, control_ann, log_fc, p_value, p_value_adj) %>%
      mutate(log_fc = signif(log_fc, 5),
             p_value = signif(p_value, 5),
             p_value_adj = signif(p_value_adj, 5)) %>%
      rename("dataset ID" = "dataset_acc", "case group" = "case_ann", "control group" = "control_ann",
             "logFC" = "log_fc", "p-value" = "p_value", "adj. p-value" = "p_value_adj") %>%
      datatable(options = list(dom = "p", pageLength = 10))
  })
  
  # Generate the reactable table
  meta_analysis_results <- reactable(
    meta_analysis_results,
    defaultColDef = colDef(
      minWidth = 75,
      headerStyle = list(background = "#f7f7f8")
    ),
    bordered = TRUE,
    striped = TRUE,
    highlight = TRUE,
    height = 310,
    pagination = FALSE,
    sortable = TRUE,
    details = function(index) {
      htmltools::div(style = "padding: 1rem",
                     htmltools::p("Data included in meta analysis: ", style = "font-weight: bold;"),
                     extra_info[[index]])
    }
  )
  
  return(meta_analysis_results)
}
#==============================================================================================================
