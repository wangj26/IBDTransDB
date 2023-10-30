#==============================================================================================================
# Function used to generate a boxplot of gene expression for every selected comparison
# Groups from comparisons from the same dataset are included in the same plot
#==============================================================================================================
# plot_expression <- function(gene, comparisons, db) {
#   
#   # Split the comparisons into dataset ID's and the comparison ID's
#   comparisons <- lapply(comparisons, function(x) str_split(x, "_")[[1]])
#   # Convert the list of dataset ID's and comparison ID's to a dataframe
#   comparison_df <- as.data.frame(do.call(rbind, comparisons))
#   colnames(comparison_df) <- c("dataset_acc", "comparison_id")
#   
#   # Pull the expression data and comparison results associated with the selected comparisons
#   expression_data <- get_selected_comparison_expression(gene, comparison_df, db)
#   comparison_results <- get_selected_comparison_results(gene, comparison_df, db)
#   
#   # Include a space after every ; in the annotation column - necessary for the tick labels to wrap in the plot
#   expression_data$annotation <- gsub(";", "; ", expression_data$annotation)
#   comparison_results$case_ann <- gsub(";", "; ", comparison_results$case_ann)
#   comparison_results$control_ann <- gsub(";", "; ", comparison_results$control_ann)
#   
#   # Extend the color palette if necessary
#   num_colors <- max(aggregate(annotation ~ dataset_acc, expression_data, function(x) length(unique(x)))$annotation)
#   plot_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
#   
#   # Create a plot for every dataset 
#   plot_list <- lapply(unique(expression_data$dataset_acc), function(x) {
#     
#     # Filter the expression data and comparison results for the given gene 
#     dataset_data <- expression_data[expression_data$dataset_acc == x,]
#     dataset_comparisons <- comparison_results[comparison_results$dataset_acc == x,]
#     
#     # Generate the y coordinates for each p-value
#     group_comparisons <- lapply(1:nrow(dataset_comparisons), function(x) {
#       c(dataset_comparisons$case_ann[x], dataset_comparisons$control_ann[x])
#     })
#     y_positions <- get_y_position(dataset_data, 
#                                   expr ~ annotation, 
#                                   comparisons = group_comparisons,
#                                   step.increase = 0.1)
#     y_positions <- select(y_positions, -c("groups"))
#     
#     # Create a dataframe used to display the p-values in the expression boxplot
#     ### Must str_wrap group1 and group2 annotations in the dataset comparisons dataframe ###
#     ### Otherwise p-values will not display properly ###
#     dataset_comparisons <- dataset_comparisons %>% 
#       select(case_ann, control_ann, p_value, p_value_adj) %>%
#       rename(group1 = case_ann, group2 = control_ann) %>%
#       mutate(p_value = formatC(p_value, format = "e", digits = 2), p_value_adj = formatC(p_value_adj, format = "e", digits = 2)) %>%
#       left_join(y_positions, by = c("group1" = "group1", "group2" = "group2")) %>%
#       mutate(group1 = stringr::str_wrap(group1, 10), group2 = stringr::str_wrap(group2, 10))
#     
#     # Boxplot of expression for the given dataset and gene
#     exp_plot <- ggplot(dataset_data, aes(x = stringr::str_wrap(annotation, 10), y = expr, fill = annotation)) + 
#       geom_boxplot() +
#       ggtitle(x) +
#       xlab(NULL) +
#       ylab(NULL) +
#       scale_fill_manual(values = plot_colors) +
#       stat_pvalue_manual(dataset_comparisons, label = "p_value", inherit.aes = FALSE) +
#       theme(panel.background = element_rect(fill = 'white'), 
#             title = element_text(size = 16), 
#             axis.text.x = element_text(size = 10, angle = 30, hjust = 1), 
#             axis.text.y = element_text(size = 12),
#             legend.position = 'None',
#             panel.border = element_rect(colour = "black", fill=NA),
#             plot.margin = margin(1, 1, 1, 40))
#     
#     return(exp_plot)
#   })
#   
#   # Generate the grid of plots
#   grouped_plot <- grid.arrange(grobs = plot_list, ncol = 2)
#   
#   return(grouped_plot)
# }
# #==============================================================================================================  



#==============================================================================================================
# Function used to determine the plot layout for a given dataset in the Expression Visuals tab
# Each comparison is represented by a boxplot
#==============================================================================================================
# plot_layout <- function(current_dataset_acc, comparisons) {
#   # Split the comparisons into dataset ID's and the comparison ID's
#   comparisons <- lapply(comparisons, function(x) str_split(x, "_")[[1]])
#   # Convert the list of dataset ID's and comparison ID's to a dataframe and filter for the given dataset ID
#   comparison_df <- as.data.frame(do.call(rbind, comparisons))
#   comparison_df <- comparison_df %>% 
#     setNames(c("dataset_acc", "comparison_id")) %>%
#     filter(dataset_acc == current_dataset_acc)
#   
#   # Determine the plot layout based on the number of datasets chosen using splitLayout
#   num_row <- nrow(comparison_df)%/%3
#   num_remaining <- nrow(comparison_df)%%3
#   if (num_remaining == 0) {
#     col_1_indices <- 1:num_row
#     col_2_indices <- (num_row + 1):(2*num_row)
#     col_3_indices <- (2*num_row + 1):nrow(comparison_df)
#   } else if (num_remaining == 1) {
#     col_1_indices <- 1:(num_row + 1)
#     col_2_indices <- (num_row + 2):(2*num_row + 1)
#     col_3_indices <- (2*num_row + 2):nrow(comparison_df)
#   } else {
#     col_1_indices <- 1:(num_row + 1)
#     col_2_indices <- (num_row + 2):(2*num_row + 2)
#     col_3_indices <- (2*num_row + 3):nrow(comparison_df)
#   }
# 
#   return(list("col_1_indices" = col_1_indices, "col_2_indices" = col_2_indices, "col_3_indices" = col_3_indices))
# }
#==============================================================================================================


#==============================================================================================================
# Alternate function used to determine the plot layout under each gene tab in the Expression Visuals tab
# Plots are organized by columns for each dataset
# Each comparison is represented by a boxplot
#==============================================================================================================
# plot_layout_2 <- function(chosen_dataset_acc, comparisons) {
#   
#   # Split the comparisons into dataset ID's and the comparison ID's
#   comparisons <- lapply(comparisons, function(x) str_split(x, "_")[[1]])
#   # Convert the list of dataset ID's and comparison ID's to a dataframe
#   comparison_df <- as.data.frame(do.call(rbind, comparisons))
#   
#   # Iterate over every column in the plot layout (6 columns)
#   plot_index_list <- list()
#   for (i in 1:6) {
#     if (i > length(chosen_dataset_acc)) {
#       plot_index_list[[i]] <- list("", "")
#     } else {
#       current_dataset_acc <- chosen_dataset_acc[i]
#       
#       # Filter for the given dataset ID
#       dataset_comparisons <- comparison_df %>% 
#         setNames(c("dataset_acc", "comparison_id")) %>%
#         filter(dataset_acc == current_dataset_acc)
#       
#       # Generate the plotOutput names based on the current dataset and comparisons chosen (including the title output)
#       plot_index_list[[i]] <- list(paste0(current_dataset_acc, "_plot_title"), 
#                                    paste0(current_dataset_acc, "_expression_plot_", c(1:nrow(dataset_comparisons))))
#     }
#   }
#   
#   return(plot_index_list)
# }
#==============================================================================================================


#==============================================================================================================
# Function used to determine the plot layout in the Expression Visuals tab
# Each comparison is represented by a boxplot
#==============================================================================================================
plot_layout <- function(genes, chosen_dataset_acc, comparisons) {
  
  
  # Convert the list of dataset ID's and comparison ID's to a dataframe
  comparison_df <- as.data.frame(do.call(rbind, lapply(comparisons, function(x) str_split(x, "_")[[1]])))
  comparison_df <- comparison_df %>% setNames(c("dataset_acc", "comparison_id"))
  # Iterate over every column in the plot layout (6 columns)
  plot_index_list <- list()
  for (i in 1:6) {
    if (i > nrow(comparison_df)) {
      plot_index_list[[i]] <- list("", "")
    } else {
      #current_dataset_acc <- chosen_dataset_acc[i]
      
      # Filter for the given dataset ID
      dataset_comparisons <- comparison_df[i,]
      print(dataset_comparisons)
      # Generate the plotOutput names based on the current dataset and comparisons chosen (including the title output)
      #current_dataset_acc <- dataset_comparisons$dataset_acc
      dataset_plots_list <- list(paste0(comparisons[i], "_plot_title"))
      dataset_plots <- list()
      for (gene in genes) {
        dataset_plots <- append(dataset_plots, paste0(gene, "_", comparisons[i], "_expression_plot_", c(1:nrow(dataset_comparisons))))
      }
      dataset_plots_list[[2]] <- dataset_plots
      plot_index_list[[i]] <- dataset_plots_list
      print(dataset_plots_list)
      print(plot_index_list)
    }
  }
  
  return(plot_index_list)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate a boxplot of gene expression for every selected comparison
# Returns a list of plots given selected comparisons from the same dataset 
#==============================================================================================================
# plot_expression <- function(gene, current_dataset_acc, comparisons, db) {
#   
#   # Split the comparisons into dataset ID's and the comparison ID's
#   comparisons <- lapply(comparisons, function(x) str_split(x, "_")[[1]])
#   # Convert the list of dataset ID's and comparison ID's to a dataframe and filter for the given dataset ID
#   comparison_df <- as.data.frame(do.call(rbind, comparisons))
#   comparison_df <- comparison_df %>% 
#     setNames(c("dataset_acc", "comparison_id")) %>%
#     filter(dataset_acc == current_dataset_acc)
#   
#   # Pull the expression data and comparison results associated with the selected comparisons
#   expression_data <- get_selected_comparison_expression(gene, comparison_df, db)
#   comparison_results <- get_selected_comparison_results(gene, comparison_df, db)
#   
  # # If the given gene is not found in the given dataset, return plots that state the given gene was not found
  # if (nrow(expression_data) == 0) {
  #   plot_list <- lapply(1:nrow(comparison_df), function(x) {
  #     exp_plot <- ggplot(data.frame("x" = 1, "y" = 1), aes(x = x, y = y)) +
  #       geom_point() +
  #       xlab(NULL) +
  #       ylab("Expression") +
  #       ggtitle(gene) +
  #       geom_label(label = paste(gene, "not found in dataset"), 
  #                  x = 1,
  #                  y = 1,
  #                  label.padding = unit(0.55, "lines"), 
  #                  size = 7,
  #                  label.size = 1,
  #                  color = "red") +
  #       theme(panel.background = element_rect(fill = 'white'),
  #             axis.text.x = element_blank(),
  #             axis.text.y = element_blank(),
  #             panel.border = element_rect(colour = "black", fill=NA),
  #             aspect.ratio = 1)
  # 
  #     return(exp_plot)
  #   })
  # 
  #   return(plot_list)
  # }
#   
#   # Include a space after every ; in the annotation column - necessary for the tick labels to wrap in the plot
#   expression_data$annotation <- gsub(";", "; ", expression_data$annotation)
#   comparison_results$case_ann <- gsub(";", "; ", comparison_results$case_ann)
#   comparison_results$control_ann <- gsub(";", "; ", comparison_results$control_ann)
#   
#   # Extend the color palette if necessary
#   num_colors <- max(aggregate(annotation ~ dataset_acc, expression_data, function(x) length(unique(x)))$annotation)
#   plot_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
#   
#   # Create a plot for every dataset 
#   plot_list <- lapply(1:nrow(comparison_results), function(x) {
#     
#     # Filter the expression data and comparison results for the given comparison 
#     dataset_comparisons <- comparison_results[x,]
#     dataset_data <- expression_data %>%
#       filter(annotation == dataset_comparisons$case_ann | annotation == dataset_comparisons$control_ann)
#     
#     # Generate the y coordinates for each p-value
#     group_comparisons <- lapply(1:nrow(dataset_comparisons), function(x) {
#       c(dataset_comparisons$case_ann[x], dataset_comparisons$control_ann[x])
#     })
#     y_positions <- get_y_position(dataset_data, 
#                                   expr ~ annotation, 
#                                   comparisons = group_comparisons,
#                                   step.increase = 0.1)
#     y_positions <- select(y_positions, -c("groups"))
#     # Increase the y position of the p-value
#     y_positions$y.position <- y_positions$y.position + 0.25
#     
#     # Create a dataframe used to display the p-values in the expression boxplot
#     ### Must str_wrap group1 and group2 annotations in the dataset comparisons dataframe ###
#     ### Otherwise p-values will not display properly ###
#     dataset_comparisons <- dataset_comparisons %>% 
#       select(case_ann, control_ann, p_value, p_value_adj) %>%
#       rename(group1 = case_ann, group2 = control_ann) %>%
#       mutate(p_value = formatC(p_value, format = "e", digits = 2), p_value_adj = formatC(p_value_adj, format = "e", digits = 2)) %>%
#       left_join(y_positions, by = c("group1" = "group1", "group2" = "group2")) %>%
#       mutate(group1 = stringr::str_wrap(group1, 10), group2 = stringr::str_wrap(group2, 10))
#     
#     # Boxplot of expression for the given dataset and gene
#     exp_plot <- ggplot(dataset_data, aes(x = stringr::str_wrap(annotation, 10), y = expr, fill = annotation)) + 
#       geom_boxplot() +
#       xlab(NULL) +
#       ylab("Expression") +
#       ggtitle(gene) +
#       scale_fill_manual(values = plot_colors) +
#       stat_pvalue_manual(dataset_comparisons, label = "p_value", size = 4, inherit.aes = FALSE) +
#       scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
#       theme(panel.background = element_rect(fill = 'white'), 
#             plot.title = element_text(size = 12),
#             axis.title.y = element_text(size = 10),
#             axis.text.x = element_text(size = 10, angle = 30, hjust = 1), 
#             axis.text.y = element_text(size = 10),
#             legend.position = 'None',
#             panel.border = element_rect(colour = "black", fill=NA),
#             plot.margin = margin(1, 1, 1, 40))
#     
#     return(exp_plot)
#   })
#   
#   return(plot_list)
# }
#==============================================================================================================  


#==============================================================================================================
# Function used to generate a boxplot of gene expression for every selected comparison
# Returns a list of plots given selected comparisons from the same dataset 
#==============================================================================================================
plot_expression <- function(gene, gene_map, current_comparison, datasets, comparisons, sample_ann, pval_option, jitter, db) {
  print("Plot expression")
  # Split the comparisons into dataset ID's and the comparison ID's
  # print(current_comparison)
  comparison_ids <- lapply(current_comparison, function(x) str_split(x, "_")[[1]])
  # print(comparison_ids)
  # Convert the list of dataset ID's and comparison ID's to a dataframe and filter for the given dataset ID
  comparison_df <- as.data.frame(do.call(rbind, comparison_ids))
  #cat("DIM:",ncol(comparison_df),"\n")
  print(comparison_df)
  comparison_df <- comparison_df %>% 
    setNames(c("dataset_acc", "comparison_id")) 
  current_dataset_acc <- comparison_df$dataset_acc
  # Pull the expression data and comparison results associated with the selected comparisons
  expression_data <- get_selected_comparison_expression(gene, gene_map, datasets, comparisons, sample_ann, db)
  comparison_results <- get_selected_comparison_results(gene, gene_map, comparison_df, db)
  
  # If the given gene is not found in the given dataset, return plots that state the given gene was not found
  if (nrow(expression_data) == 0) {
    plot_list <- lapply(1:nrow(comparison_df), function(x) {
      exp_plot <- ggplot(data.frame("x" = 1, "y" = 1), aes(x = x, y = y)) +
        geom_point() +
        xlab(NULL) +
        ylab("Expression") +
        ggtitle(gene) +
        geom_label(label = paste(gene, "not found in dataset"), 
                   x = 1,
                   y = 1,
                   label.padding = unit(0.55, "lines"), 
                   size = 7,
                   label.size = 1,
                   color = "red") +
        theme(panel.background = element_rect(fill = 'white'), 
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              aspect.ratio = 1)
      
      return(exp_plot)
    })
    
    return(plot_list)
  }
  
  # Include a space after every ; in the annotation column - necessary for the tick labels to wrap in the plot
  expression_data$annotation <- gsub(";", "; ", expression_data$annotation)
  comparison_results$case_ann <- gsub(";", "; ", comparison_results$case_ann)
  comparison_results$control_ann <- gsub(";", "; ", comparison_results$control_ann)
  
  # Define the color palette
  plot_colors <- brewer.pal(8, "Dark2")[1:2]
  
  # Create a plot for every dataset 
  plot_list <- lapply(1:nrow(comparison_results), function(x) {
    print(x)
    # Filter the expression data and comparison results for the given comparison 
    dataset_comparisons <- comparison_results[x,]
    dataset_data <- expression_data %>%
      filter(annotation == dataset_comparisons$case_ann | annotation == dataset_comparisons$control_ann) %>%
      mutate(annotation = factor(annotation, levels = c(dataset_comparisons$case_ann, dataset_comparisons$control_ann)))
    
    # Generate the y coordinates for the p-value and p-value line
    comparison_exp <- pull(dataset_data, expr)
    p_val_line_position <- max(comparison_exp) + (max(comparison_exp) - min(comparison_exp))*0.1
    p_val_position <- max(comparison_exp) + (max(comparison_exp) - min(comparison_exp))*0.15
    
    # Create a dataframe used to display the p-values in the expression boxplot
    dataset_comparisons <- dataset_comparisons %>% 
      select(case_ann, control_ann, p_value, p_value_adj) %>%
      mutate(p_value = formatC(p_value, format = "e", digits = 2), p_value_adj = formatC(p_value_adj, format = "e", digits = 2))
    
    # Boxplot of expression for the given dataset and gene
    exp_plot <- ggplot(dataset_data, aes(x = stringr::str_wrap(annotation, 10), y = expr, fill = annotation)) + 
      geom_boxplot() +
      annotate("text", x = 1.5, y = p_val_position, label = dataset_comparisons[[pval_option]], size = 4) +
      geom_segment(aes(x = 1, y = p_val_line_position, xend = 2, yend = p_val_line_position)) +
      xlab(NULL) +
      ylab("Expression") +
      ggtitle(gene) +
      scale_fill_manual(values = plot_colors) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      theme(panel.background = element_rect(fill = 'white'), 
            plot.title = element_text(size = 12),
            axis.title.y = element_text(size = 10),
            axis.text.x = element_text(size = 10, angle = 30, hjust = 1), 
            axis.text.y = element_text(size = 10),
            legend.position = 'None',
            panel.border = element_rect(colour = "black", fill=NA),
            plot.margin = margin(1, 1, 1, 40))
    if (jitter) {
      exp_plot <- exp_plot + geom_jitter()
    }
    
    return(exp_plot)
  })
  
  return(plot_list)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate all boxplots of gene expression for every selected comparison
# This is used by the downloadHandler when users download plots
# Returns a list of lists containing plots from every comparison selected across datasets
#==============================================================================================================
# plots_for_download <- function(genes, comparisons, db) {
#   
#   # Split the comparisons into dataset ID's and the comparison ID's
#   comparisons <- lapply(comparisons, function(x) str_split(x, "_")[[1]])
#   # Convert the list of dataset ID's and comparison ID's to a dataframe and filter for the given dataset ID
#   comparison_df_total <- as.data.frame(do.call(rbind, comparisons))
#   comparison_df_total <- comparison_df_total %>% 
#     setNames(c("dataset_acc", "comparison_id"))
#   
#   # Define the unique dataset ID's
#   unique_dataset_acc <- unique(comparison_df_total$dataset_acc)
#     
#   # Initialize a list used to store the expression plots
#   plot_list <- list()
#   
#   ##########
#   for (gene in genes) {
#     
#     dataset_plot_list <- list()
#     
#     ##########
#     for (current_dataset_acc in unique_dataset_acc) {
#       
#       # Filter the comparison dataframe for the given dataset ID
#       comparison_df <- filter(comparison_df_total, dataset_acc == current_dataset_acc)
#       
#       # Pull the expression data and comparison results associated with the selected comparisons
#       expression_data <- get_selected_comparison_expression(gene, comparison_df, db)
#       comparison_results <- get_selected_comparison_results(gene, comparison_df, db)
#       
#       # If the given gene is not found in the given dataset, continue to next iteration
#       if (nrow(expression_data) == 0) {
#         next
#       }
#       
#       # Include a space after every ; in the annotation column - necessary for the tick labels to wrap in the plot
#       expression_data$annotation <- gsub(";", "; ", expression_data$annotation)
#       comparison_results$case_ann <- gsub(";", "; ", comparison_results$case_ann)
#       comparison_results$control_ann <- gsub(";", "; ", comparison_results$control_ann)
#       
#       # Extend the color palette if necessary
#       num_colors <- max(aggregate(annotation ~ dataset_acc, expression_data, function(x) length(unique(x)))$annotation)
#       plot_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(num_colors)
#       
#       ### Create a plot for every comparison for the given dataset ###
#       exp_plots <- lapply(1:nrow(comparison_results), function(x) {
#         
#         # Filter the expression data and comparison results for the given comparison 
#         dataset_comparisons <- comparison_results[x,]
#         dataset_data <- expression_data %>%
#           filter(annotation == dataset_comparisons$case_ann | annotation == dataset_comparisons$control_ann)
#         
#         # Generate the y coordinates for each p-value
#         group_comparisons <- lapply(1:nrow(dataset_comparisons), function(x) {
#           c(dataset_comparisons$case_ann[x], dataset_comparisons$control_ann[x])
#         })
#         y_positions <- get_y_position(dataset_data, 
#                                       expr ~ annotation, 
#                                       comparisons = group_comparisons,
#                                       step.increase = 0.1)
#         y_positions <- select(y_positions, -c("groups"))
#         # Increase the y position of the p-value
#         y_positions$y.position <- y_positions$y.position + 0.25
#         
#         # Create a dataframe used to display the p-values in the expression boxplot
#         ### Must str_wrap group1 and group2 annotations in the dataset comparisons dataframe ###
#         ### Otherwise p-values will not display properly ###
#         dataset_comparisons <- dataset_comparisons %>% 
#           select(case_ann, control_ann, p_value, p_value_adj) %>%
#           rename(group1 = case_ann, group2 = control_ann) %>%
#           mutate(p_value = formatC(p_value, format = "e", digits = 2), p_value_adj = formatC(p_value_adj, format = "e", digits = 2)) %>%
#           left_join(y_positions, by = c("group1" = "group1", "group2" = "group2")) %>%
#           mutate(group1 = stringr::str_wrap(group1, 10), group2 = stringr::str_wrap(group2, 10))
#         
#         # Boxplot of expression for the given dataset and gene
#         exp_plot <- ggplot(dataset_data, aes(x = stringr::str_wrap(annotation, 10), y = expr, fill = annotation)) + 
#           geom_boxplot() +
#           xlab(NULL) +
#           ylab("Expression") +
#           ggtitle(paste0(gene, " (", current_dataset_acc, ")")) +
#           scale_fill_manual(values = plot_colors) +
#           stat_pvalue_manual(dataset_comparisons, label = "p_value", size = 4, inherit.aes = FALSE) +
#           scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
#           theme(panel.background = element_rect(fill = 'white'), 
#                 plot.title = element_text(size = 12),
#                 axis.title.y = element_text(size = 10),
#                 axis.text.x = element_text(size = 10, angle = 30, hjust = 1), 
#                 axis.text.y = element_text(size = 10),
#                 legend.position = 'None',
#                 panel.border = element_rect(colour = "black", fill=NA),
#                 plot.margin = margin(1, 1, 1, 40))
#         
#         return(exp_plot)
#       })
#       ##########
#       
#       # Append the plot(s) to dataset_plot_list
#       dataset_plot_list[[current_dataset_acc]] <- exp_plots
#     }
#     ##########
#     
#     # Append the plot(s) for each dataset to plot_list for the given gene
#     plot_list[[gene]] <- dataset_plot_list
#     
#   }
#   ##########
#   
#   return(plot_list)
# }
#==============================================================================================================  


#==============================================================================================================
# Function used to generate all boxplots of gene expression for every selected comparison
# This is used by the downloadHandler when users download plots
# Returns a list of lists containing plots from every comparison selected across datasets
#==============================================================================================================
plots_for_download <- function(genes, gene_map, datasets, comparisons, sample_ann, pval_option, jitter, db) {
  # print(comparisons)
  # Split the comparisons into dataset ID's and the comparison ID's
  comparisons_ids <- lapply(comparisons, function(x) str_split(x, "_")[[1]])
  # print(comparisons_ids)
  # Convert the list of dataset ID's and comparison ID's to a dataframe and filter for the given dataset ID
  comparison_df_total <- as.data.frame(do.call(rbind, comparisons_ids))
  comparison_df_total <- comparison_df_total %>% 
    setNames(c("dataset_acc", "comparison_id"))
  print(comparison_df_total)
  # Define the unique dataset ID's
  comparison_df <- comparison_df_total$comparison_id
  
  # Initialize a list used to store the expression plots
  plot_list <- list()
  
  ##########
  for (gene in genes) {
    print(gene)
    dataset_plot_list <- list()
    expression_data <- get_selected_comparison_expression(gene, gene_map, datasets, comparisons, sample_ann, db)
    print("exp_data")
    print(head(expression_data))
    expression_data$annotation <- gsub(";", "; ", expression_data$annotation)
    ##########
    for (current_comparison in comparison_df) {
      print(current_comparison)
      # Filter the comparison dataframe for the given dataset ID
      #comparison_df <- filter(comparison_df_total, dataset_acc == current_dataset_acc)
      
      # Pull the expression data and comparison results associated with the selected comparisons
      cmf <- filter(comparison_df_total,comparison_id==current_comparison)
      comparison_results <- get_selected_comparison_results(gene, gene_map, cmf, db)
      
      # If the given gene is not found in the given dataset, continue to next iteration
      if (nrow(expression_data) == 0) {
        next
      }
      
      # Include a space after every ; in the annotation column - necessary for the tick labels to wrap in the plot
      print("1")
      comparison_results$case_ann <- gsub(";", "; ", comparison_results$case_ann)
      comparison_results$control_ann <- gsub(";", "; ", comparison_results$control_ann)
      print("2")
      
      # Define the color palette
      plot_colors <- brewer.pal(8, "Dark2")[1:2]
      
      ### Create a plot for every comparison for the given dataset ###
      exp_plots <- lapply(1:nrow(comparison_results), function(x) {
      print("3")
        
        # Filter the expression data and comparison results for the given comparison 
        dataset_comparisons <- comparison_results[x,]
        dataset_data <- expression_data %>%
          filter(annotation == dataset_comparisons$case_ann | annotation == dataset_comparisons$control_ann)
        print("4")  
        # Generate the y coordinates for the p-value and p-value line
        comparison_exp <- pull(dataset_data, expr)
        p_val_line_position <- max(comparison_exp) + (max(comparison_exp) - min(comparison_exp))*0.1
        p_val_position <- max(comparison_exp) + (max(comparison_exp) - min(comparison_exp))*0.15
        
        # Create a dataframe used to display the p-values in the expression boxplot
        dataset_comparisons <- dataset_comparisons %>% 
          select(case_ann, control_ann, p_value, p_value_adj) %>%
          mutate(p_value = formatC(p_value, format = "e", digits = 2), p_value_adj = formatC(p_value_adj, format = "e", digits = 2))
        print("5")
        # Boxplot of expression for the given dataset and gene
        exp_plot <- ggplot(dataset_data, aes(x = stringr::str_wrap(annotation, 10), y = expr, fill = annotation)) + 
          geom_boxplot() +
          annotate("text", x = 1.5, y = p_val_position, label = dataset_comparisons[[pval_option]], size = 4) +
          geom_segment(aes(x = 1, y = p_val_line_position, xend = 2, yend = p_val_line_position)) +
          xlab(NULL) +
          ylab("Expression") +
          ggtitle(paste0(gene, " (", cmf$dataset_acc, ")")) +
          scale_fill_manual(values = plot_colors) +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
          theme(panel.background = element_rect(fill = 'white'), 
                plot.title = element_text(size = 12),
                axis.title.y = element_text(size = 10),
                axis.text.x = element_text(size = 10, angle = 30, hjust = 1), 
                axis.text.y = element_text(size = 10),
                legend.position = 'None',
                panel.border = element_rect(colour = "black", fill=NA),
                plot.margin = margin(1, 1, 1, 40))
        if (jitter) {
          exp_plot <- exp_plot + geom_jitter()
        }
        
        return(exp_plot)
      })
      ##########
      
      # Append the plot(s) to dataset_plot_list
      dataset_plot_list[[current_comparison]] <- exp_plots
    }
    ##########
    
    # Append the plot(s) for each dataset to plot_list for the given gene
    plot_list[[gene]] <- dataset_plot_list
    
  }
  ##########
  
  return(plot_list)
}
#==============================================================================================================  
