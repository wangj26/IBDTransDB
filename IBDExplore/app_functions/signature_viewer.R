#==============================================================================================================
# Function used to map samples to groups (used for gene expression group options)
#==============================================================================================================
map_samples_to_groups <- function(comparisons) {
  # Create a two column dataframe (first contains group annotation, second contains sample IDs)
  case_df <- comparisons %>%
    dplyr::select(c(case_ann, case_sample)) %>%
    dplyr::rename(group = case_ann, samples = case_sample)
  control_df <- comparisons %>%
    dplyr::select(c(control_ann, control_sample)) %>%
    dplyr::rename(group = control_ann, samples = control_sample)
  sample_group_map <- rbind(case_df, control_df)
  sample_group_map <- unique(sample_group_map)
  
  return(sample_group_map)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate an expression boxplot for each selected comparison (p-value is displayed)
# Returns a list of ggplot objects
#==============================================================================================================
# generate_comparison_exp_plots <- function(selected_genes, selected_comparisons, exp_data, sample_group_map, comparisons, comparison_data) {
# 
#   # Iterate over each selected comparison and produce an expression plot
#   plot_list <- lapply(selected_comparisons, function(selected_comparison) {
# 
#     # Use the selected comparison ID to filter the comparisons and comparison data
#     selected_comparison_id <- comparisons %>%
#       filter(comparison == selected_comparison) %>%
#       pull(id)
#     comparisons <- filter(comparisons, comparison == selected_comparison)
#     comparison_data <- comparison_data %>%
#       filter(comparison_id == selected_comparison_id) %>%
#       filter(gene %in% selected_genes)
# 
#     # Join the selected comparisons and comparison results dataframes
#     comparison_results <- left_join(comparisons, comparison_data, by = c("id" = "comparison_id")) %>%
#       filter(row_number() == 1)
# 
#     # Filter the sample to group map for the selected groups and explode the dataframe (each row corresponds to a sample)
#     selected_groups <- c(comparison_results$case_ann, comparison_results$control_ann)
#     sample_group_map <- sample_group_map %>%
#       filter(group %in% selected_groups) %>%
#       separate_rows(samples, sep = ";")
# 
#     # Filter the expression data for the given gene and groups and add a column specifying the group
#     # Average the expression values across genes for each sample (if multiple genes are specified)
#     # Include a space after every ; in the annotation column - necessary for the tick labels to wrap in the plot
#     exp_data <- exp_data %>%
#       filter(gene %in% selected_genes) %>%
#       filter(sample_acc %in% sample_group_map$samples) %>%
#       left_join(sample_group_map, by = c("sample_acc" = "samples")) %>%
#       select(-gene) %>%
#       group_by(sample_acc, group) %>%
#       summarize(expr = mean(expr)) %>%
#       as.data.frame() %>%
#       mutate(group = gsub(";", "; ", group)) %>%
#       mutate(group = stringr::str_wrap(group, 10))
# 
#     # Generate the y coordinates for each p-value
#     group_comparisons <- list(c(comparison_results$case_ann[1], comparison_results$control_ann[1]))
#     y_positions <- get_y_position(exp_data,
#                                   expr ~ group,
#                                   comparisons = group_comparisons,
#                                   step.increase = 0.1)
#     y_positions <- select(y_positions, -c("groups"))
#     # Increase the y position of the p-value
#     y_positions$y.position <- y_positions$y.position + 0.25
# 
#     if (length(selected_genes) == 1) {
#       # Create a dataframe used to display the p-values in the expression boxplot
#       ### Must str_wrap group1 and group2 annotations in the dataset comparisons dataframe ###
#       ### Otherwise p-values will not display properly ###
#       comparison_results <- comparison_results %>%
#         select(case_ann, control_ann, p_value, p_value_adj) %>%
#         rename(group1 = case_ann, group2 = control_ann) %>%
#         mutate(p_value = formatC(p_value, format = "e", digits = 2),
#                p_value_adj = formatC(p_value_adj, format = "e", digits = 2)) %>%
#         left_join(y_positions, by = c("group1" = "group1", "group2" = "group2")) %>%
#         mutate(group1 = gsub(";", "; ", group1), group2 = gsub(";", "; ", group2)) %>%
#         mutate(group1 = stringr::str_wrap(group1, 10), group2 = stringr::str_wrap(group2, 10))
#     } else {
#       # Perform a Wilcoxon rank sum test between the two groups for the given comparison (given a signature is specified)
#       group1_name <- unlist(strsplit(selected_comparison, split = " vs "))[1]
#       group1_name <- stringr::str_wrap(gsub(";", "; ", group1_name), 10)
#       group1_values <- exp_data %>%
#         filter(group == group1_name) %>%
#         pull(expr)
#       group2_name <- unlist(strsplit(selected_comparison, split = " vs "))[2]
#       group2_name <- stringr::str_wrap(gsub(";", "; ", group2_name), 10)
#       group2_values <- exp_data %>%
#         filter(group == group2_name) %>%
#         pull(expr)
#       test_results <- wilcox.test(group1_values, group2_values, alternative = "two.sided")
#       p_val <- signif(test_results$p.value, 4)
# 
#       # Create the comparison results dataframe
#       comparison_results <- comparison_results %>%
#         filter(row_number() == 1) %>%
#         select(case_ann, control_ann, p_value_adj) %>%
#         rename(group1 = case_ann, group2 = control_ann) %>%
#         mutate(p_value_adj = p_val) %>%
#         mutate(p_value_adj = formatC(p_value_adj, format = "e", digits = 2)) %>%
#         left_join(y_positions, by = c("group1" = "group1", "group2" = "group2")) %>%
#         mutate(group1 = gsub(";", "; ", group1), group2 = gsub(";", "; ", group2)) %>%
#         mutate(group1 = stringr::str_wrap(group1, 10), group2 = stringr::str_wrap(group2, 10))
#     }
# 
# 
#     # Define the color palette
#     plot_colors <- brewer.pal(8, "Dark2")[1:2]
# 
#     # Boxplot of expression for the given comparison and gene
#     exp_plot <- ggplot(exp_data, aes(x = group, y = expr, fill = group)) +
#       geom_boxplot() +
#       xlab(NULL) +
#       ylab("Expression") +
#       ggtitle(NULL) +
#       scale_fill_manual(values = plot_colors) +
#       stat_pvalue_manual(comparison_results, label = "p_value_adj", size = 4, inherit.aes = FALSE) +
#       scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
#       theme(panel.background = element_rect(fill = "white"),
#             plot.title = element_text(size = 14),
#             axis.title.y = element_text(size = 10),
#             axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
#             axis.text.y = element_text(size = 12),
#             legend.position = "None",
#             panel.border = element_rect(colour = "black", fill=NA),
#             plot.margin = margin(1, 1, 1, 40))
# 
#     return(exp_plot)
#   })
# 
#   # Define the names for the plot list (necessary to create plotOutput objects)
#   names(plot_list) <- paste0(gsub(" ", "_", selected_comparisons), "_plot")
# 
#   return(plot_list)
# }
#==============================================================================================================


#==============================================================================================================
# Function used to generate an expression boxplot for each selected comparison (p-value is displayed)
# Returns a list of ggplot objects
#==============================================================================================================
generate_comparison_exp_plots <- function(selected_genes, selected_comparisons, exp_data, sample_group_map, comparisons, comparison_data, pval_option, jitter) {

  # Iterate over each selected comparison and produce an expression plot
  plot_list <- lapply(selected_comparisons, function(selected_comparison) {

    # Use the selected comparison ID to filter the comparisons and comparison data
    selected_comparison_id <- comparisons %>%
      filter(comparison == selected_comparison) %>%
      pull(id)
    comparisons <- filter(comparisons, comparison == selected_comparison)
    comparison_data <- comparison_data %>%
      filter(comparison_id == selected_comparison_id) %>%
      filter(gene %in% selected_genes)

    # Join the selected comparisons and comparison results dataframes
    comparison_results <- left_join(comparisons, comparison_data, by = c("id" = "comparison_id")) %>%
      filter(row_number() == 1)

    # Filter the sample to group map for the selected groups and explode the dataframe (each row corresponds to a sample)
    selected_groups <- c(comparison_results$case_ann, comparison_results$control_ann)
    sample_group_map <- sample_group_map %>%
      filter(group %in% selected_groups) %>%
      separate_rows(samples, sep = ";")

    # Filter the expression data for the given gene and groups and add a column specifying the group
    # Average the expression values across genes for each sample (if multiple genes are specified)
    # Include a space after every ; in the annotation column - necessary for the tick labels to wrap in the plot
    exp_data <- exp_data %>%
      filter(gene %in% selected_genes) %>%
      filter(sample_acc %in% sample_group_map$samples) %>%
      left_join(sample_group_map, by = c("sample_acc" = "samples")) %>%
      select(-gene) %>%
      group_by(sample_acc, group) %>%
      summarize(expr = mean(expr)) %>%
      as.data.frame() %>%
      mutate(group = factor(group, levels = c(comparison_results$case_ann, comparison_results$control_ann)),
             group = gsub(";", "; ", group),
             group = stringr::str_wrap(group, 10))

    # Generate the y coordinates for the p-value and p-value line
    comparison_exp <- pull(exp_data, expr)
    p_val_line_position <- max(comparison_exp) + (max(comparison_exp) - min(comparison_exp))*0.1
    p_val_position <- max(comparison_exp) + (max(comparison_exp) - min(comparison_exp))*0.15

    if (length(selected_genes) == 1) {
      # Create a dataframe used to display the p-values in the expression boxplot
      ### Must str_wrap group1 and group2 annotations in the dataset comparisons dataframe ###
      ### Otherwise p-values will not display properly ###
      comparison_results <- comparison_results %>%
        select(case_ann, control_ann, p_value, p_value_adj) %>%
        mutate(p_value = formatC(p_value, format = "e", digits = 2),
               p_value_adj = formatC(p_value_adj, format = "e", digits = 2))
      pval <- ifelse(pval_option == "p_value", comparison_results$p_value, comparison_results$p_value_adj)
    } else {
      # Perform a Wilcoxon rank sum test between the two groups for the given comparison (given a signature is specified)
      group1_name <- unlist(strsplit(selected_comparison, split = " vs "))[1]
      group1_name <- stringr::str_wrap(gsub(";", "; ", group1_name), 10)
      group1_values <- exp_data %>%
        filter(group == group1_name) %>%
        pull(expr)
      group2_name <- unlist(strsplit(selected_comparison, split = " vs "))[2]
      group2_name <- stringr::str_wrap(gsub(";", "; ", group2_name), 10)
      group2_values <- exp_data %>%
        filter(group == group2_name) %>%
        pull(expr)
      test_results <- wilcox.test(group1_values, group2_values, alternative = "two.sided")
      pval <- signif(test_results$p.value, 4)

      # Create the comparison results dataframe
      comparison_results <- comparison_results %>%
        filter(row_number() == 1) %>%
        select(case_ann, control_ann, p_value_adj) %>%
        mutate(p_value = pval) %>%
        mutate(p_value = formatC(p_value, format = "e", digits = 2))
      pval <- comparison_results$p_value
    }


    # Define the color palette
    plot_colors <- brewer.pal(8, "Dark2")[1:2]

    # Boxplot of expression for the given comparison and gene
    exp_plot <- ggplot(exp_data, aes(x = group, y = expr, fill = group)) +
      geom_boxplot() +
      annotate("text", x = 1.5, y = p_val_position, label = pval, size = 4) +
      geom_segment(aes(x = 1, y = p_val_line_position, xend = 2, yend = p_val_line_position)) +
      xlab(NULL) +
      ylab("Expression") +
      ggtitle(NULL) +
      scale_fill_manual(values = plot_colors) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      theme(panel.background = element_rect(fill = "white"),
            plot.title = element_text(size = 14),
            axis.title.y = element_text(size = 10),
            axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
            axis.text.y = element_text(size = 12),
            legend.position = "None",
            panel.border = element_rect(colour = "black", fill=NA),
            plot.margin = margin(1, 1, 1, 40))
    if (jitter) {
      exp_plot <- exp_plot + geom_jitter()
    }
    
    return(exp_plot)
  })

  # Define the names for the plot list (necessary to create plotOutput objects)
  names(plot_list) <- paste0(gsub(" ", "_", selected_comparisons), "_plot")

  return(plot_list)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate an expression boxplot for selected groups
# Returns one ggplot object
#==============================================================================================================
generate_grouped_exp_plot <- function(selected_genes, selected_groups, exp_data, sample_group_map, jitter) {
  
  # Filter the sample to group map for the selected groups and explode the dataframe (each row corresponds to a sample)
  sample_group_map <- sample_group_map %>%
    filter(group %in% selected_groups) %>%
    separate_rows(samples, sep = ";")
  
  # Filter the expression data for the given gene and groups and add a column specifying the group
  # Average the expression values across genes for each sample (if multiple genes are specified)
  exp_data <- exp_data %>%
    filter(gene %in% selected_genes) %>%
    filter(sample_acc %in% sample_group_map$samples) %>%
    left_join(sample_group_map, by = c("sample_acc" = "samples")) %>%
    select(-gene) %>%
    group_by(sample_acc, group) %>%
    summarize(expr = mean(expr)) %>%
    as.data.frame()

  # Include a space after every ; in the annotation column - necessary for the tick labels to wrap in the plot
  exp_data$group <- gsub(";", "; ", exp_data$group)
  
  # Extend the color palette if necessary
  num_colors <- length(unique(exp_data$group))
  plot_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
  
  # Boxplot of expression for the given dataset and gene
  exp_plot <- ggplot(exp_data, aes(x = stringr::str_wrap(group, 10), y = expr, fill = group)) + 
    geom_boxplot() +
    xlab(NULL) +
    ylab("Expression") +
    ggtitle(NULL) +
    scale_fill_manual(values = plot_colors) +
    theme(panel.background = element_rect(fill = "white"), 
          plot.title = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14, angle = 30, hjust = 1), 
          axis.text.y = element_text(size = 14),
          legend.position = "None",
          panel.border = element_rect(colour = "black", fill=NA),
          plot.margin = margin(1, 1, 1, 40))
  if (jitter) {
    exp_plot <- exp_plot + geom_jitter()
  }
  
  return(exp_plot)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate a table of comparisons for the generated expression plot
# Returns one DT table object
#==============================================================================================================
generate_signature_plot_table <- function(selected_genes, selected_groups, comparisons, exp_data, sample_group_map) {
  
  # Filter the comparisons for the selected groups
  comparisons <- comparisons %>%
    filter(case_ann %in% selected_groups, control_ann %in% selected_groups) %>%
    select(case_ann, control_ann)
  
  if (nrow(comparisons) == 0) {
    return(NULL)
  }
  
  # Filter the sample to group map for the selected groups and explode the dataframe (each row corresponds to a sample)
  sample_group_map <- sample_group_map %>%
    filter(group %in% selected_groups) %>%
    separate_rows(samples, sep = ";")
  
  # Filter the expression data for the given gene and groups and add a column specifying the group
  # Average the expression values across genes for each sample (if multiple genes are specified)
  exp_data <- exp_data %>%
    filter(gene %in% selected_genes) %>%
    filter(sample_acc %in% sample_group_map$samples) %>%
    left_join(sample_group_map, by = c("sample_acc" = "samples")) %>%
    select(-gene) %>%
    group_by(sample_acc, group) %>%
    summarize(expr = mean(expr)) %>%
    as.data.frame()
  
  # Perform a Wilcoxon rank sum test between the two groups for each viable comparison (given a signature is specified)
  test_stats <- lapply(1:nrow(comparisons), function(x) {
    
    # Compute the Wilcoxon rank sum p-value
    group1_name <- comparisons$case_ann[x]
    group1_values <- exp_data %>%
      filter(group == group1_name) %>%
      pull(expr)
    group2_name <- comparisons$control_ann[x]
    group2_values <- exp_data %>%
      filter(group == group2_name) %>%
      pull(expr)
    test_results <- wilcox.test(group1_values, group2_values, alternative = "two.sided")
    p_val <- signif(test_results$p.value, 4)
    
    # Compute the group 1 mean and median
    group1_mean <- signif(mean(group1_values), 4)
    group1_median <- signif(median(group1_values), 4)
    
    # Compute the group 2 mean and median
    group2_mean <- signif(mean(group2_values), 4)
    group2_median <- signif(median(group2_values), 4)
    
    return(c(group1_mean, group1_median, group2_mean, group2_median, p_val))
  })
  test_stats <- data.frame(do.call(rbind, test_stats))
  colnames(test_stats) <- c("case mean", "case median", "control mean", "control median", "p-value")
  
  comparisons <- cbind(comparisons, test_stats)
  
  # Generate the DT table
  signature_plot_table <- comparisons %>%
    relocate("control_ann", .after = "case median") %>%
    mutate(case_ann = gsub(";", ", ", case_ann), control_ann = gsub(";", ", ", control_ann)) %>%
    rename("case group" = "case_ann", "control group" = "control_ann") %>%
    datatable(options = list(dom = "ti",
                             order = list(list(1, 4, "desc")),
                             autoWidth = TRUE),
              rownames = FALSE) %>%
    formatStyle(columns = c(1, 4), width = "175px") %>%
    formatStyle(columns = 1:ncol(comparisons), fontSize = "90%")
  
  return(signature_plot_table)
}
#==============================================================================================================
