#==============================================================================================================
# Function used to generate volcano plots
# Returns a ggplot object
#==============================================================================================================
generate_volcano_plot <- function(comparisons, selected_comparison, comparison_data, pval_option, log_fc_threshold, p_value_threshold) {
  
  # Use the selected comparison ID to filter the comparison data
  selected_comparison_id <- comparisons %>%
    filter(comparison == selected_comparison) %>%
    pull(id)
  comparison_data <- comparison_data %>%
    filter(comparison_id == selected_comparison_id) %>%
    dplyr::select(-comparison_id)
  
  # Determine significant genes based on the fold change and p-value (these will be colored red)
  comparison_data$significant <- comparison_data[[pval_option]] <= p_value_threshold & abs(comparison_data$log_fc) >= log_fc_threshold
  # Significant data points are solid and all others have some transparency
  alpha <- ifelse(comparison_data$significant, 1, 0.5)
  
  # Create the volcano plot
  pval_label <- ifelse(pval_option == "p_value", "p-value:", "adj p-value:")
  volcano_plot <- ggplot(data = comparison_data, aes(x = log_fc, 
                                                     y = -log10(.data[[pval_option]]), 
                                                     color = significant,
                                                     text = paste("name:", gene,
                                                                  "\nlogFC:", formatC(log_fc, format = "f", digits = 3), 
                                                                  "\n", pval_label, formatC(p_value, format = "e", digits = 3)))) +
    geom_point(size = 1.5, alpha = alpha) +
    geom_hline(yintercept = -log10(p_value_threshold), color = "blue", size = 0.3) +
    geom_vline(xintercept = -log_fc_threshold, color = "blue", size = 0.3) +
    geom_vline(xintercept = log_fc_threshold, color = "blue", size = 0.3) +
    ggtitle(selected_comparison) +
    xlab("log2(fold change)") +
    ylab(ifelse(pval_option == "p_value", "-log10(p-value)", "-log10(adj p-value)")) +
    theme_bw() +
    scale_color_manual(values = c("black", "red")) +
    theme(plot.title = element_text(size = 14, hjust = 0.5), axis.title = element_text(size = 12),
          axis.text = element_text(size = 12, color = "black"), axis.text.x = element_text(vjust = .5), 
          legend.position = "none", aspect.ratio = 0.8, panel.grid = element_blank())
  # Convert to plotly object
  # volcano_plot <- ggplotly(volcano_plot, tooltip = "text") %>% 
  #   style(hoverinfo = "none") %>% 
  #   config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "zoom2d", 
  #                                     "pan2d", "select2d", "lasso2d", 
  #                                     "drawclosedpath", "drawopenpath",
  #                                     "drawline", "drawrect", "drawcircle"), 
  #          displaylogo = FALSE)

  return(volcano_plot)
}
#==============================================================================================================


#==============================================================================================================
# Function used to label the current volcano plot with selected genes
# Returns a ggplot object
#==============================================================================================================
update_volcano_plot <- function(selected_genes, comparisons, selected_comparison, comparison_data, pval_option, log_fc_threshold, p_value_threshold) {
  
  # Use the selected comparison ID to filter the comparison data
  selected_comparison_id <- comparisons %>%
    filter(comparison == selected_comparison) %>%
    pull(id)
  comparison_data <- comparison_data %>%
    filter(comparison_id == selected_comparison_id) %>%
    select(-comparison_id)
  
  # Determine significant genes based on the fold change and p-value (these will be colored red)
  comparison_data$significant <- comparison_data[[pval_option]] <= p_value_threshold & abs(comparison_data$log_fc) >= log_fc_threshold
  # Significant data points and labeled data points are solid and all others have some transparency
  alpha <- ifelse(comparison_data$significant | comparison_data$gene %in% selected_genes, 1, 0.5)
  # Change the values in the significant column
  comparison_data$significant <- ifelse(comparison_data$significant, "significant", "nonsignificant")
  
  # Determine the genes to label in the volcano plot
  comparison_data$label <- comparison_data$gene %in% selected_genes
  # Change the significant column to indicate genes to be labeled (used for coloring the associated data points)
  comparison_data <- comparison_data %>%
    mutate(significant = replace(significant, label, "to_label"))
  
  # Determine the max -log10(p-value) and the min and max log2(fold change) (used for restraining the labels to a specified area)
  min_x <- min(comparison_data$log_fc)
  max_x <- max(comparison_data$log_fc)
  max_y <- max(-log10(comparison_data$p_value_adj))
  
  # Create the volcano plot
  pval_label <- ifelse(pval_option == "p_value", "p-value:", "adj p-value:")
  volcano_plot <- ggplot(data = comparison_data, aes(x = log_fc, 
                                                     y = -log10(.data[[pval_option]]), 
                                                     color = significant, 
                                                     text = paste("name:", gene,
                                                                  "\nlogFC:", formatC(log_fc, format = "f", digits = 3), 
                                                                  "\n", pval_label, formatC(p_value, format = "e", digits = 3)))) +
    geom_point(size = 1.5, alpha = alpha) +
    geom_hline(yintercept = -log10(p_value_threshold), color = "blue", size = 0.3) +
    geom_vline(xintercept = -log_fc_threshold, color = "blue", size = 0.3) +
    geom_vline(xintercept = log_fc_threshold, color = "blue", size = 0.3) +
    ggtitle(selected_comparison) +
    xlab("log2(fold change)") +
    ylab(ifelse(pval_option == "p_value", "-log10(p-value)", "-log10(adj p-value)")) +
    theme_bw() +
    scale_color_manual(values = c("black", "red", "#8C7E10")) +
    geom_label(data = filter(comparison_data, label), 
               aes(label = gene),
               subset = comparison_data$label,
               nudge_x = 0.5, 
               nudge_y = 0.5) +
    # geom_text_repel(aes(label = ifelse(comparison_data$label, as.character(comparison_data$gene), ""), hjust = 0, vjust = 0), 
    #                 size = 5,
    #                 force = 1,
    #                 max.overlaps = Inf,
    #                 xlim = c(min_x/2, max_x/2),
    #                 ylim = c(max_y/2, Inf)) +
    theme(plot.title = element_text(size = 14, hjust = 0.5), axis.title = element_text(size = 12),
          axis.text = element_text(size = 12, color = "black"), axis.text.x = element_text(vjust = .5), 
          legend.position = "none", aspect.ratio = 0.8, panel.grid = element_blank())
  # Convert to plotly object
  # volcano_plot <- ggplotly(volcano_plot, tooltip = "text") %>% 
  #   style(hoverinfo = "none") %>% 
  #   config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "zoom2d", 
  #                                     "pan2d", "select2d", "lasso2d", 
  #                                     "drawclosedpath", "drawopenpath",
  #                                     "drawline", "drawrect", "drawcircle"), 
  #          displaylogo = FALSE)
  
  return(volcano_plot)
}
#==============================================================================================================


#==============================================================================================================
# Function used to label the current volcano plot with selected geneset 
# Color all unselected points the same color and do not label selected genes
# Returns a ggplot object
#==============================================================================================================
update_volcano_plot_2 <- function(selected_genes, comparisons, selected_comparison, comparison_data, pval_option, log_fc_threshold, p_value_threshold) {
  
  # Use the selected comparison ID to filter the comparison data
  selected_comparison_id <- comparisons %>%
    filter(comparison == selected_comparison) %>%
    pull(id)
  comparison_data <- comparison_data %>%
    filter(comparison_id == selected_comparison_id) %>%
    select(-comparison_id)
  
  # Determine the genes to be highlighted in the plot
  comparison_data$selected <- comparison_data$gene %in% selected_genes
  # Unselected data points have some transparency
  alpha <- ifelse(comparison_data$selected, 1, 0.1)
  
  # Create the volcano plot
  pval_label <- ifelse(pval_option == "p_value", "p-value:", "adj p-value:")
  volcano_plot <- ggplot(data = comparison_data, aes(x = log_fc, 
                                                     y = -log10(.data[[pval_option]]), 
                                                     color = selected, 
                                                     text = paste("name:", gene,
                                                                  "\nlogFC:", formatC(log_fc, format = "f", digits = 3), 
                                                                  "\n", pval_label, formatC(p_value, format = "e", digits = 3)))) +
    geom_point(size = 1.5, alpha = alpha) +
    geom_hline(yintercept = -log10(p_value_threshold), color = "blue", size = 0.3) +
    geom_vline(xintercept = -log_fc_threshold, color = "blue", size = 0.3) +
    geom_vline(xintercept = log_fc_threshold, color = "blue", size = 0.3) +
    ggtitle(selected_comparison) +
    xlab("log2(fold change)") +
    ylab(ifelse(pval_option == "p_value", "-log10(p-value)", "-log10(adj p-value)")) +
    theme_bw() +
    scale_color_manual(values = c("#D3D3D3", "#FFA500")) +
    theme(plot.title = element_text(size = 14, hjust = 0.5), axis.title = element_text(size = 12),
          axis.text = element_text(size = 12, color = "black"), axis.text.x = element_text(vjust = .5), 
          legend.position = "none", aspect.ratio = 0.8, panel.grid = element_blank())
  # Convert to plotly object
  # volcano_plot <- ggplotly(volcano_plot, tooltip = "text") %>% 
  #   style(hoverinfo = "none") %>% 
  #   config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "zoom2d", 
  #                                     "pan2d", "select2d", "lasso2d", 
  #                                     "drawclosedpath", "drawopenpath",
  #                                     "drawline", "drawrect", "drawcircle"), 
  #          displaylogo = FALSE)
  
  return(volcano_plot)
}
#==============================================================================================================


#==============================================================================================================
# Function used to determine the number of DE genes based on the p-value and logFC thresholds
# Returns a string specifying the number of DE genes
#==============================================================================================================
determine_num_signif_genes <- function(comparisons, selected_comparison, comparison_data, pval_option, log_fc_threshold, p_value_threshold) {
  
  # Use the selected comparison ID to filter the comparison data
  selected_comparison_id <- comparisons %>%
    filter(comparison == selected_comparison) %>%
    pull(id)
  comparison_data <- comparison_data %>%
    filter(comparison_id == selected_comparison_id) %>%
    dplyr::select(-comparison_id)
  
  # Filter the comparison data for genes that satisfy the fold change and p-value thresholds
  rows_to_keep <- comparison_data[[pval_option]] <= p_value_threshold & abs(comparison_data$log_fc) >= log_fc_threshold
  comparison_data <- comparison_data %>%
    filter(rows_to_keep) 
  # Determine the number of DE genes
  num_total <- nrow(comparison_data)
  num_up <- nrow(filter(comparison_data, log_fc > 0))
  num_down <- nrow(filter(comparison_data, log_fc < 0))
  
  # The statement specifying the number of DE genes
  if (pval_option == "p_value") {
    signif_genes_statement <- div(
      p("The number of differentially expressed genes based on a p-value threshold of", 
        p_value_threshold, 
        "and a |logFC| threshold of", 
        log_fc_threshold,
        ":"),
      p("Total:", strong(num_total)),
      p("Up:", strong(num_up)),
      p("Down:", strong(num_down))
    )
  } else {
    signif_genes_statement <- div(
      p("The number of differentially expressed genes based on an adj. p-value threshold of", 
        p_value_threshold, 
        "and a |logFC| threshold of", 
        log_fc_threshold,
        ":"),
      p("Total:", strong(num_total)),
      p("Up:", strong(num_up)),
      p("Down:", strong(num_down))
    )
  }

  return(signif_genes_statement)
}
#==============================================================================================================


#==============================================================================================================
# Function used to create a DT table from a dataframe of comparison metrics for a given comparison
# Returns a DT table
#==============================================================================================================
generate_comparison_table <- function(comparisons, selected_comparison, comparison_data, active_comparison_data) {
  
  # Use the selected comparison ID to filter the comparison data
  selected_comparison_id <- comparisons %>%
    filter(comparison == selected_comparison) %>%
    pull(id)
  comparison_data <- comparison_data %>%
    filter(comparison_id == selected_comparison_id) %>%
    select(-comparison_id) %>%
    arrange(p_value_adj) %>%
    mutate(log_fc = signif(log_fc, 4),
           p_value = signif(p_value, 4),
           p_value_adj = signif(p_value_adj, 4)) %>%
    head(50)
  
  ### Update reactive value for the active comparison data ###
  active_comparison_data(comparison_data)
  
  # Generate the DT table
  comparison_data <- comparison_data %>%
    rename("logFC" = "log_fc", "p-value" = "p_value", "adj p-value" = "p_value_adj") %>%
    datatable(comparison_data,
              extensions = "Scroller",
              selection = "multiple",
              options = list(dom = "ti",
                             order = list(list(4, "asc")),
                             autoWidth = TRUE,
                             scrollX = TRUE,
                             scrollY = 450,
                             scroller = TRUE)) %>%
    formatStyle(columns = c(1, 2, 3, 4), width = "150px") %>%
    formatStyle(columns = c(1, 2, 3, 4), fontSize = "90%")
  
  return(comparison_data)
}
#==============================================================================================================


#==============================================================================================================
# Function used to update the DT table with comparison metrics with selected genes
# Returns a DT table
#==============================================================================================================
update_comparison_table <- function(selected_genes, comparisons, selected_comparison, comparison_data, active_comparison_data) {
  
  # Use the selected comparison ID to filter the comparison data
  selected_comparison_id <- comparisons %>%
    filter(comparison == selected_comparison) %>%
    pull(id)
  comparison_data <- comparison_data %>%
    filter(comparison_id == selected_comparison_id) %>%
    select(-comparison_id) %>%
    filter(gene %in% selected_genes) %>%
    mutate(log_fc = signif(log_fc, 4),
           p_value = signif(p_value, 4),
           p_value_adj = signif(p_value_adj, 4))
  
  ### Update reactive value for the active comparison data ###
  active_comparison_data(comparison_data)
  
  # Generate the DT table
  comparison_data <- comparison_data %>%
    rename("logFC" = "log_fc", "p-value" = "p_value", "adj p-value" = "p_value_adj") %>%
    datatable(comparison_data,
              extensions = "Scroller",
              selection = "multiple",
              options = list(dom = "ti",
                             order = list(list(4, "asc")),
                             autoWidth = TRUE,
                             scrollX = TRUE,
                             scrollY = 450,
                             scroller = TRUE)) %>%
    formatStyle(columns = c(1, 2, 3, 4), width = "150px") %>%
    formatStyle(columns = c(1, 2, 3, 4), fontSize = "90%")
  
  return(comparison_data)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate expression boxplots based on the selected rows in the comparison metrics table
# Returns a list of ggplot objects
#==============================================================================================================
generate_gene_plots <- function(selected_rows, comparisons, selected_comparison, comparison_data, exp_data, sample_group_map, pval_option, jitter) {
  
  # Iterate over each selected row and produce an expression plot
  plot_list <- lapply(selected_rows, function(selected_row) {
    # Pull the gene associated with the selected row
    selected_gene <- comparison_data$gene[selected_row]

    # Use the selected comparison ID to filter the comparisons and comparison data
    comparisons <- filter(comparisons, comparison == selected_comparison)
    comparison_data <- comparison_data %>%
      filter(gene == selected_gene)
    
    # Join the selected comparisons and comparison data dataframes
    comparison_results <- cbind(comparisons, comparison_data)

    # Filter the sample to group map for the selected groups and explode the dataframe (each row corresponds to a sample)
    selected_groups <- c(comparison_results$case_ann, comparison_results$control_ann)
    sample_group_map <- sample_group_map %>%
      filter(group %in% selected_groups) %>%
      separate_rows(samples, sep = ";")

    # Filter the expression data for the given gene and groups and add a column specifying the group
    exp_data <- exp_data %>%
      filter(gene == selected_gene) %>%
      filter(sample_acc %in% sample_group_map$samples) %>%
      left_join(sample_group_map, by = c("sample_acc" = "samples")) %>%
      select(-gene) %>%
      mutate(group = factor(group, levels = c(comparison_results$case_ann, comparison_results$control_ann)),
             group = gsub(";", "; ", group),
             group = stringr::str_wrap(group, 10))

    # Generate the y coordinates for the p-value and p-value line
    comparison_exp <- pull(exp_data, expr)
    p_val_line_position <- max(comparison_exp) + (max(comparison_exp) - min(comparison_exp))*0.1
    p_val_position <- max(comparison_exp) + (max(comparison_exp) - min(comparison_exp))*0.15

    # Create a dataframe used to display the p-values in the expression boxplot
    ### Must str_wrap group1 and group2 annotations in the dataset comparisons dataframe ###
    ### Otherwise p-values will not display properly ###
    comparison_results <- comparison_results %>%
      select(case_ann, control_ann, p_value, p_value_adj) %>%
      mutate(p_value = formatC(p_value, format = "e", digits = 2), 
             p_value_adj = formatC(p_value_adj, format = "e", digits = 2))

    # Define the color palette
    plot_colors <- brewer.pal(8, "Dark2")[1:2]

    # Boxplot of expression for the given comparison and gene
    exp_plot <- ggplot(exp_data, aes(x = group, y = expr, fill = group)) +
      geom_boxplot() +
      annotate("text", x = 1.5, y = p_val_position, label = comparison_results[[pval_option]], size = 4) +
      geom_segment(aes(x = 1, y = p_val_line_position, xend = 2, yend = p_val_line_position)) +
      xlab(NULL) +
      ylab("Expression") +
      ggtitle(selected_gene) +
      scale_fill_manual(values = plot_colors) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      theme(panel.background = element_rect(fill = 'white'),
            plot.title = element_text(size = 14),
            axis.title.y = element_text(size = 10),
            axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
            axis.text.y = element_text(size = 12),
            legend.position = "None",
            panel.border = element_rect(colour = "black", fill = NA),
            plot.margin = margin(1, 1, 1, 40))
    if (jitter) {
      exp_plot <- exp_plot + geom_jitter()
    }
    
    return(exp_plot)
  })

  # Define the names for the plot list (necessary to create plotOutput objects)
  names(plot_list) <- paste0(comparison_data$gene[selected_rows], "_plot")

  return(plot_list)
}
#==============================================================================================================


#==============================================================================================================
# Function used to determine the plot layout for the gene expression plots
#==============================================================================================================
determine_plot_layout <- function(num_plots) {
  
  # Maximum number of plots to display is 9
  if (num_plots > 9) {
    num_plots <- 9
  }
  
  # Determine the plot layout based on the number of selected rows in the comparison metrics table
  num_row <- num_plots%/%3
  num_remaining <- num_plots%%3
  if (num_remaining == 0) {
    col_1_indices <- 1:num_row
    col_2_indices <- (num_row + 1):(2*num_row)
    col_3_indices <- (2*num_row + 1):num_plots
  } else if (num_remaining == 1) {
    col_1_indices <- 1:(num_row + 1)
    col_2_indices <- (num_row + 2):(2*num_row + 1)
    col_3_indices <- (2*num_row + 2):num_plots
  } else {
    col_1_indices <- 1:(num_row + 1)
    col_2_indices <- (num_row + 2):(2*num_row + 2)
    col_3_indices <- (2*num_row + 3):num_plots
  }
  
  return(list(col_1_indices, col_2_indices, col_3_indices))
}
#==============================================================================================================
